"""How far PlasmidScope's ORFs and annotations can stand in for ours.

1. ORF concordance. Calls genes with our gene caller (darkorf.genecall, same settings as
   config orf.min_call_length_aa) on a seeded random sample of analysis-set plasmids, and
   compares each ORF with PlasmidScope's ORFs on the same plasmid:
     IDENTICAL  same protein sequence (the join our pipeline would use)
     SAME_STOP  same strand and stop codon, different start
     ABSENT     no PlasmidScope ORF on that strand and stop
2. Annotation agreement. Cross-tabulates PlasmidScope's ps_class against our cascade's
   functional_class for every protein of a finished pipeline run (by seq_id).

Usage:
  python tools/plasmidscope_concordance.py \
      data/PlasmidScope/annotation/analysis_set_orfs.tsv.gz \
      data/plasmidscope_primary/provenance/working_set.fna.gz \
      data/plasmidscope_primary/analysis_set.tsv \
      results_test/05_annotation_cascade/protein_annotation.tsv \
      <n_sample> <seed>
"""
import collections
import csv
import gzip
import multiprocessing
import pathlib
import random
import statistics
import sys

import yaml

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from darkorf import genecall  # noqa: E402
from plasmidann.dereplicate import _seq_id  # noqa: E402
from plasmidann.fasta import iter_fasta  # noqa: E402

ps_orfs, fasta, analysis_set, cascade, n_sample, seed = sys.argv[1:7]
MIN_AA = yaml.safe_load(open(ROOT / "config/config.yaml"))["orf"]["min_call_length_aa"]
csv.field_size_limit(sys.maxsize)


def stop_of(start, end, strand):
    return int(end) if strand in ("+", "1") else int(start)


# --- PlasmidScope ORFs, and each protein's class --------------------------------------
ps_by_plasmid = collections.defaultdict(list)
with gzip.open(ps_orfs, "rt", newline="") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        ps_by_plasmid[r["plasmid_id"]].append(
            {k: r[k] for k in ("orf_source", "start", "end", "strand", "seq_id")})
ps_class_by_seq = {}
with gzip.open(ps_orfs.replace("_orfs.", "_proteins."), "rt", newline="") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        ps_class_by_seq[r["seq_id"]] = f"{r['ps_class']}/named={r['product_named']}"

# --- 1. ORF concordance on a random sample --------------------------------------------
with open(analysis_set, newline="") as fh:
    topology = {r["plasmid_id"]: r["topology"] for r in csv.DictReader(fh, delimiter="\t")}
eligible = sorted(p for p in topology if p in ps_by_plasmid)
sample = set(random.Random(int(seed)).sample(eligible, int(n_sample)))

records = [((pid, seq), topology[pid]) for pid, seq in iter_fasta([fasta]) if pid in sample]
with multiprocessing.Pool(initializer=genecall.configure,
                          initargs=((MIN_AA + 1) * 3,)) as pool:
    called = pool.map(genecall.call_record, records, chunksize=16)

ours = collections.Counter()
ours_len = collections.defaultdict(list)
ps_matched = collections.Counter()
ps_total = collections.Counter()
for pid, genes, _ in called:
    theirs = ps_by_plasmid[pid]
    # PlasmidScope called ORFs itself only for some sources; the rest are deposited CDS.
    group = "prodigal" if theirs[0]["orf_source"].startswith("Prodigal") else "deposited"
    ps_total[group] += len(theirs)
    their_seqs = {r["seq_id"] for r in theirs}
    their_stops = {(r["strand"], stop_of(r["start"], r["end"], r["strand"])) for r in theirs}
    our_seqs = set()
    for g in genes:
        sid = _seq_id(g["seq"])
        our_seqs.add(sid)
        strand = "+" if g["strand"] == 1 else "-"
        if sid in their_seqs:
            kind = "IDENTICAL"
        elif (strand, stop_of(g["start"], g["end"], strand)) in their_stops:
            kind = "SAME_STOP"
        else:
            kind = "ABSENT"
        kind += "_origin" if g["origin_spanning"] else ""
        ours[(group, kind)] += 1
        ours_len[(group, kind)].append(len(g["seq"]))
    ps_matched[group] += sum(r["seq_id"] in our_seqs for r in theirs)

print(f"== ORF concordance: {len(called):,} random analysis-set plasmids (seed {seed})")
for group in ("prodigal", "deposited"):
    n_ours = sum(v for (g, _), v in ours.items() if g == group)
    print(f"-- PlasmidScope ORFs from {group}: our ORFs {n_ours:,}, "
          f"PlasmidScope ORFs {ps_total[group]:,}")
    for (g, kind) in sorted(k for k in ours if k[0] == group):
        n = ours[(g, kind)]
        print(f"  ours {kind:<18} {n:>8,}  {n / n_ours:6.1%}"
              f"   median {statistics.median(ours_len[(g, kind)]):.0f} aa")
    print(f"  PlasmidScope ORFs with an identical protein among ours "
          f"{ps_matched[group]:,}  {ps_matched[group] / ps_total[group]:.1%}")

# --- 2. PlasmidScope class against our cascade ----------------------------------------
table = collections.Counter()
with open(cascade, newline="") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        table[(ps_class_by_seq.get(r["seq_id"], "NOT_IN_PS"), r["functional_class"])] += 1

rows = [f"{c}/named={n}" for c in ("ANNOTATED", "UNKNOWN_ORTHOLOG", "NONE") for n in (1, 0)]
rows.append("NOT_IN_PS")
cols = sorted({c for _, c in table})
print(f"\n== ps_class (rows) x our functional_class (columns), unique proteins of {cascade}")
print(f"{'':<26}" + "".join(f"{c:>25}" for c in cols) + f"{'total':>8}")
for r in rows:
    n = sum(table[(r, c)] for c in cols)
    print(f"{r:<26}" + "".join(f"{table[(r, c)]:>25,}" for c in cols) + f"{n:>8,}")
