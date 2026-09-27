"""The dark set: the proteins that are candidates at all.

Target-eligible proteins (rule target_eligibility: UNCHARACTERIZED_HOMOLOG or NONE, not
artefact-flagged), minus proteins that occur only as partial ORFs. A fragment aligns to only
part of a domain, drifts toward DOMAIN_ONLY or out of annotation, and looks like a novel
dark protein, so a protein without one complete ORF is set aside and counted.
"""
import csv

import _ctx  # noqa: F401

eligible = set()
with open(snakemake.input.flags, newline="") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        if r["target_eligible"] == "1":
            eligible.add(r["seq_id"])

# Any unique protein represented ONLY by partial ORFs is a fragment, whatever else is true
# of it. A protein with at least one complete ORF representative is kept.
complete = set()
partial_only = set()
with open(snakemake.input.map) as fh:
    seq_of_orf = {}
    for line in fh:
        sid, ids = line.rstrip("\n").split("\t")
        for oid in ids.split(","):
            seq_of_orf[oid] = sid

with open(snakemake.input.index, newline="") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        # Dereplication is lossless, so every ORF of the index is in the map.
        sid = seq_of_orf[r["orf_id"]]
        if r["partial"] == "0":
            complete.add(sid)
        else:
            partial_only.add(sid)
partial_only -= complete

dark = eligible - partial_only

n_written = 0
with open(snakemake.output.faa, "w") as out:
    emit = False
    for line in open(snakemake.input.faa):
        if line[0] == ">":
            emit = line[1:].split()[0] in dark
            n_written += emit
        if emit:
            out.write(line)

with open(snakemake.output.ids, "w") as out:
    for sid in sorted(dark):
        out.write(sid + "\n")

print(f"eligible={len(eligible)} excluded_partial_only={len(eligible & partial_only)} "
      f"dark_set={n_written}")
assert n_written > 0, ("the dark set is empty - check the target-eligibility flags and the "
                       "cascade output")
