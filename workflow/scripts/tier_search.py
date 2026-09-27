"""S3: search one tier, then hand the next tier only what stayed unexplained.

Runs the tier's tool (hmmsearch, pharokka or DIAMOND) over the proteins the previous tier
left unexplained, rejects hits above the tier's max_evalue before they can contribute
anything, adds the spans of informative hits to each protein's explained fraction, and
writes hits.tsv, the cumulative spans.tsv and unresolved.faa: the proteins still below
narrow_at, plus the sweep cohort, which is never narrowed away.
"""
import csv
import os
import pathlib
import subprocess

import _ctx  # noqa: F401  - puts src/ on sys.path for the plasmidann package

from plasmidann import pharokka, scratch
from plasmidann.cascade import (
    as_float,
    explained_fraction,
    is_informative,
    narrow_by_explained,
    passes_significance,
)

spec = snakemake.params.spec              # one entry from cascade.yaml: tiers
faa = snakemake.input.faa                 # this tier's queries
max_evalue = spec["max_evalue"]           # None for bit-score tiers such as T1 (--cut_ga)
hmmer_z = snakemake.params.hmmer_z        # pinned -Z / --domZ, see docs/annotation_statistics.md
max_target_seqs = snakemake.params.max_target_seqs

ids = [l[1:].split()[0] for l in open(faa) if l[0] == ">"]

# Proteins an earlier tier already NAMED from a source this tier defers to are not searched
# here (skip_if_named_by in config/cascade.yaml): a protein Pfam or Swiss-Prot named keeps
# that curated name, and nr is spent only on what they could not name. They still pass
# through the narrowing below with the explained fraction they arrived with. The sweep
# cohort is not exempt: the rule is about which database may name a protein.
skip_sources = set(spec.get("skip_if_named_by") or [])
skipped = set()
if skip_sources:
    id_set0 = set(ids)
    for path in snakemake.input.named:
        with open(path, newline="") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                if r["informative"] == "True" and r["source"] in skip_sources \
                        and r["query"] in id_set0:
                    skipped.add(r["query"])
search_ids = [i for i in ids if i not in skipped]

# Proteins that bypass narrowing and are searched by every tier, so that the cost of
# narrowing can be measured on them. An empty file (fraction 0) disables the cohort.
sweep_ids = set()
if getattr(snakemake.input, "sweep", None):
    sweep_ids = {l.strip() for l in open(snakemake.input.sweep) if l.strip()}

# What earlier tiers already explained. The spans hold only informative alignments: a
# 'hypothetical protein' hit explains nothing, however well it aligns.
spans, qlen = {}, {}
prev = snakemake.input.spans
for f in ([prev] if isinstance(prev, str) else list(prev)):
    with open(f, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            qlen[r["seq_id"]] = int(r["qlen"])
            spans[r["seq_id"]] = [tuple(map(int, iv.split("-")))
                                  for iv in filter(None, r["intervals"].split(";"))]

named = {}     # every informative hit per query, this tier
best = {}      # which of them supplies the winning label
unnamed = {}   # every uninformative hit, kept as evidence someone else has seen this protein
n_rejected = 0 # hits dropped for insignificance - reported, so the gate is visible


def record(q, label, qcov, tcov, ev, start, end, tlen, accession="", source=None,
           category="", identity="", align_length="", bitscore="", target_length=""):
    """Admit one parsed hit, or reject it when it is above the tier's max_evalue.

    `start`/`end` are None for a tier that reports no alignment span (pharokka): the hit
    then adds nothing to the explained fraction and its columns are written empty.
    identity, align_length, bitscore and target_length are written empty where the tool
    reports none. Every informative hit is kept, not one per tier, so that
    n_informative_hits counts domains and explained_fraction can be recomputed from
    hits.tsv. The winning label is the most significant hit, not the widest.
    """
    global n_rejected
    if not passes_significance(ev, max_evalue):
        n_rejected += 1
        return

    qlen[q] = tlen
    spanned = start is not None and end is not None
    informative = is_informative(label)
    hit = {"query": q, "label": label, "target_accession": accession,
           "coverage": round(qcov, 4) if qcov is not None else "",
           "target_coverage": round(tcov, 4) if tcov is not None else "",
           "evalue": ev, "informative": informative,
           "start": start if spanned else "", "end": end if spanned else "",
           "source": source or spec["source"], "category": category,
           "identity": identity, "align_length": align_length, "bitscore": bitscore,
           "target_length": target_length}
    if not informative:
        # Kept with its coordinates: they become dark_covered_fraction at cascade_resolve.
        unnamed.setdefault(q, []).append(hit)
        return
    if spanned:
        spans.setdefault(q, []).append((start, end))
    named.setdefault(q, []).append(hit)
    if q not in best or as_float(ev) < as_float(best[q]["evalue"]):
        best[q] = hit


if search_ids:
    # A stable name: each tier has its own output directory, so no two jobs share it, and
    # a killed attempt's directory (a partial res.m8 of the nr tier is ~GBs) is cleared by
    # the rerun instead of accumulating beside the new one.
    tmp = scratch.scratch_dir(pathlib.Path(snakemake.output.hits).parent, "search_tmp")
    query = faa
    if skipped:
        query = f"{tmp}/query.faa"
        with open(query, "w") as out:
            emit = False
            for line in open(faa):
                if line[0] == ">":
                    emit = line[1:].split()[0] not in skipped
                if emit:
                    out.write(line)

    if spec["method"] == "hmmer":
        raw = f"{tmp}/dom.tbl"
        # -Z / --domZ pin the search space. hmmsearch reports E = (sequences searched) x
        # P(score | null), and each tier's input is the previous tier's residue, so without
        # them an E-value would mean something different on every tier. Pinned, it means
        # "expected false positives across the whole study" (docs/annotation_statistics.md).
        cmd = (f"hmmsearch {spec['args']} -Z {hmmer_z} --domZ {hmmer_z} "
               f"--noali --cpu {snakemake.threads} --domtblout {raw} {spec['db']} {query}")
        subprocess.run(cmd, shell=True, check=True, stdout=subprocess.DEVNULL)

        # --domtblout columns (0-indexed). hmmsearch searches profiles against sequences,
        # so the "target" is the protein and the "query" is the Pfam family:
        #   0  target name  = the protein       2  tlen = protein length
        #   3  query name   = the family        4  query accession (PF06970.19)
        #   5  qlen         = profile length    12 i-Evalue of this domain
        #  13  domain score (bits)              15,16 hmm from,to
        #  17,18 ali from,to = coordinates on the protein
        # The i-Evalue, not the full-sequence E-value (6), governs an individual domain, as
        # the HMMER user guide advises. The accession is kept because it is stable across
        # Pfam releases, while a family name can be changed by a curator.
        for line in open(raw):
            if line.startswith("#"):
                continue
            f = line.split()
            tlen, hlen = int(f[2]), int(f[5])
            a, b = int(f[17]), int(f[18])
            record(f[0], f[3],
                   (b - a + 1) / tlen if tlen else 0.0,
                   (int(f[16]) - int(f[15]) + 1) / hlen if hlen else 0.0,
                   f[12], a, b, tlen, accession=f[4], bitscore=f[13],
                   target_length=hlen)

    elif spec["method"] == "pharokka":
        # The phage tier: pharokka in protein mode, used as its authors ship it. It runs in
        # its own environment, named by `exe`; pharokka checks for its helpers on PATH, so
        # that environment's bin directory is put there for the call. See
        # plasmidann.pharokka for what the output does and does not carry.
        exe = pathlib.Path(spec["exe"]).resolve()
        outdir = f"{tmp}/pharokka"
        evalue_flag = "" if max_evalue is None else f"-e {max_evalue} "
        subprocess.run(
            f"{exe} proteins -i {query} -o {outdir} -d {spec['db']} -t {snakemake.threads} "
            f"-f -p tier {evalue_flag}{spec.get('args', '')}",
            shell=True, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
            env={**os.environ, "PATH": f"{exe.parent}:{os.environ.get('PATH', '')}"})

        merged = pathlib.Path(outdir) / "tier_full_merged_output.tsv"
        if not merged.exists():
            raise SystemExit(
                f"[{spec['id']}] pharokka exited 0 but wrote no {merged.name}. See "
                f"{outdir}/logs; the scratch directory is kept.")
        for hit in pharokka.parse_merged(merged.read_text()):
            # No span: start and end are None. The label goes through is_informative as at
            # every other tier, so a family named 'hypothetical protein' names nothing.
            record(hit["query"], hit["label"], None, None, hit["evalue"], None, None,
                   hit["query_length"], accession=hit["family_id"],
                   source=hit["source"], category=hit["category"])

    elif spec["method"] == "diamond":
        raw = f"{tmp}/res.m8"
        # The E-value threshold is declared per tier rather than left at DIAMOND's default
        # of 0.001. With max_evalue: null the flag is omitted: `--evalue None` is accepted
        # by DIAMOND, returns zero hits and exits 0.
        evalue_flag = "" if max_evalue is None else f"--evalue {max_evalue} "
        # --tmpdir puts DIAMOND's spill space on node-local disk; by default it is the
        # output directory on the shared filesystem, where the 3,910-query benchmark ran
        # 15% over the local-disk prediction. sseqid (the accession) is kept beside stitle
        # (free text) because only the accession joins to UniProt or RefSeq.
        cmd = (f"diamond blastp -q {query} -d {spec['db']} -o {raw} {spec['args']} "
               f"{evalue_flag}--threads {snakemake.threads} "
               f"--tmpdir {snakemake.resources.tmpdir} "
               f"--max-target-seqs {max_target_seqs} "
               f"--outfmt 6 qseqid sseqid stitle qcovhsp scovhsp evalue qstart qend qlen "
               f"pident length bitscore slen --quiet")
        subprocess.run(cmd, shell=True, check=True)
        for line in open(raw):
            (q, sid, title, qc, tc, ev, qs, qe, ql,
             pid, alen, bits, slen) = line.rstrip("\n").split("\t")
            record(q, title, float(qc) / 100, float(tc) / 100, ev,
                   int(qs), int(qe), int(ql), accession=sid, identity=float(pid) / 100,
                   align_length=int(alen), bitscore=bits, target_length=int(slen))

    else:
        raise SystemExit(
            f"tier {spec['id']} declares method {spec['method']!r}, which is not one of "
            "hmmer, pharokka or diamond.")

# ------------------------------------------------------------------------------------
# Outputs
# ------------------------------------------------------------------------------------
# start/end give cascade_resolve the coordinates of uninformative hits (dark coverage).
# `is_best` marks the hit that supplies the winning label for this protein at this tier;
# the other hits stay, because the merged coverage is built from all of them. `source` is
# the database the tier searched, declared per tier in config; `category` is the
# functional group a database assigns to a family where it has one (pharokka, CARD).
# Every row carries the tier arguments and max_evalue that produced it.
cols = ["query", "label", "target_accession", "coverage", "target_coverage", "evalue",
        "informative", "is_best", "start", "end", "tier", "source", "category",
        "threshold", "max_evalue", "identity", "align_length", "bitscore", "target_length"]
with open(snakemake.output.hits, "w", newline="") as out:
    w = csv.DictWriter(out, fieldnames=cols, delimiter="\t")
    w.writeheader()
    all_hits = ([x for v in named.values() for x in v]
                + [x for v in unnamed.values() for x in v])
    winners = {id(h) for h in best.values()}
    for h in all_hits:
        w.writerow({**h, "is_best": int(id(h) in winners), "tier": spec["id"],
                    "threshold": spec["args"], "max_evalue": max_evalue})

explained = {q: explained_fraction(qlen.get(q, 0), spans.get(q, [])) for q in qlen}
with open(snakemake.output.spans, "w", newline="") as out:
    w = csv.writer(out, delimiter="\t")
    w.writerow(["seq_id", "qlen", "intervals", "explained_fraction"])
    for q in qlen:
        w.writerow([q, qlen[q], ";".join(f"{a}-{b}" for a, b in spans.get(q, [])),
                    explained.get(q, 0.0)])

# narrow_at, NOT min_explained: narrowing decides what to keep searching; min_explained is
# applied afterwards, at cascade_resolve. The sweep cohort is never narrowed away.
keep = set(narrow_by_explained(ids, explained, snakemake.params.narrow_at))
cohort = sweep_ids & set(ids)
keep |= cohort

with open(snakemake.output.unresolved, "w") as out:
    emit = False
    for line in open(faa):
        if line[0] == ">":
            emit = line[1:].split()[0] in keep
        if emit:
            out.write(line)

print(f"[{spec['id']}] received={len(ids)} skipped_named_by_{'_'.join(sorted(skip_sources)) or 'none'}"
      f"={len(skipped)} queried={len(search_ids)} informative={len(best)} "
      f"informative_hits={sum(len(v) for v in named.values())} unnamed={len(unnamed)} "
      f"rejected_insignificant={n_rejected} carried_forward={len(keep)} "
      f"(sweep_cohort={len(cohort)})")

# The scratch directory is removed only here, on the ordinary path: a script that raised
# never reaches this line, and the raw domtbl or m8 is what a tier failure is diagnosed
# from. Guarded: an empty tier created no directory.
if search_ids:
    scratch.release(tmp)
