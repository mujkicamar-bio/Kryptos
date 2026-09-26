"""S2s: the proteins the cascade annotates, and the ones it actually searches.

See src/plasmidann/selection.py for the rule. Two steps:

1. SELECT by family (clustering.primary, intermediate): proteins neither annotated by
   Tier 0 nor flagged by AntiFam (S2b), in a family holding a small-plasmid protein that
   is neither.
2. SEARCH the representatives of a 90% clustering of the selected proteins (identity and
   coverage in config/cascade.yaml search_clustering). Coverage is required of BOTH
   sequences (cov-mode 0): a member then has the representative's length to within 20%,
   so the representative's result - which is copied to the member at cascade_resolve -
   cannot describe a domain the member does not have.

Output selection.tsv, one row per unique protein:
    seq_id, on_small, role, search_representative
role is one of
    artefact_antifam  flagged by AntiFam (S2b) as a probable non-protein; skips every
                      annotation tier, Tier 0 included. Other artefact flags (low
                      complexity) skip nothing
    plasmidscope      annotated by Tier 0; not searched
    representative    searched by the cascade
    member            not searched; takes its representative's result
    not_selected      in no family with an unexplained small-plasmid protein; not searched
"""
import collections
import csv
import pathlib
import subprocess

import _ctx  # noqa: F401

from plasmidann import scratch
from plasmidann.fasta import iter_fasta
from plasmidann.selection import select

cfg = snakemake.params.search

small_plasmids = {l.strip() for l in open(snakemake.input.small_ids) if l.strip()}
on_small = set()
with open(snakemake.input.map) as fh:
    for line in fh:
        sid, orf_ids = line.rstrip("\n").split("\t")
        if any(o.rsplit("|", 1)[0] in small_plasmids for o in orf_ids.split(",")):
            on_small.add(sid)

with open(snakemake.input.ps, newline="") as fh:
    ps_annotated = {r["seq_id"] for r in csv.DictReader(fh, delimiter="\t")
                    if r["ps_class"] == "ANNOTATED"}

clusters = collections.defaultdict(list)
with open(snakemake.input.families) as fh:
    for line in fh:
        rep, mem = line.rstrip("\n").split("\t")
        clusters[rep].append(mem)

# AntiFam flags only: a low-complexity flag marks composition, not a known artefact family,
# and such a protein is still searched.
with open(snakemake.input.artefact, newline="") as fh:
    antifam = {r["seq_id"] for r in csv.DictReader(fh, delimiter="\t") if r["antifam_family"]}

selected = select(clusters, on_small, ps_annotated | antifam)

# ---- the search clustering, over the selected proteins only -------------------------
work = scratch.scratch_dir(pathlib.Path(snakemake.output.tsv).parent, "search_clustering")
selected_faa = work / "selected.faa"
all_ids = []
with open(selected_faa, "w") as out:
    for sid, seq in iter_fasta([snakemake.input.faa]):
        all_ids.append(sid)
        if sid in selected:
            out.write(f">{sid}\n{seq}\n")

representative_of = {}
if selected:
    prefix = work / "search"
    subprocess.run(
        f"mmseqs easy-cluster {selected_faa} {prefix} {work / 'tmp'} "
        f"--min-seq-id {cfg['min_seq_id']} -c {cfg['coverage']} "
        f"--cov-mode {cfg['cov_mode']} --cluster-mode {cfg['cluster_mode']} "
        f"--threads {snakemake.threads} -v 1",
        shell=True, check=True)
    with open(f"{prefix}_cluster.tsv") as fh:
        for line in fh:
            rep, mem = line.rstrip("\n").split("\t")
            representative_of[mem] = rep
    # The FASTA of the representatives is the cascade's query (prepare_control).
    (work / "search_rep_seq.fasta").replace(snakemake.output.faa)
else:
    open(snakemake.output.faa, "w").close()
assert set(representative_of) == selected, "the search clustering lost selected proteins"

roles = collections.Counter()
with open(snakemake.output.tsv, "w", newline="") as out:
    w = csv.writer(out, delimiter="\t")
    w.writerow(["seq_id", "on_small", "role", "search_representative"])
    for sid in all_ids:
        rep = representative_of.get(sid, "")
        role = ("artefact_antifam" if sid in antifam else
                "plasmidscope" if sid in ps_annotated else
                "representative" if rep == sid else
                "member" if rep else "not_selected")
        roles[role] += 1
        w.writerow([sid, int(sid in on_small), role, rep])
scratch.release(work)

print(f"unique proteins {len(all_ids)}, on a small plasmid {len(on_small)}; "
      f"selected {len(selected)} in families with an unexplained small-plasmid protein; "
      + " ".join(f"{k}={v}" for k, v in sorted(roles.items())))
