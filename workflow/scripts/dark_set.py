"""S6a: assemble the dark set - the proteins that are candidates at all.

Dark = the cascade named nothing (UNCHARACTERIZED_HOMOLOG or NONE), minus everything S5
excluded: backbone, artefact, and ORFs still flagged partial after S1 origin repair.

Partial ORFs are excluded rather than repaired at this point because a fragment is a real
protein that is not whole: it aligns to only part of a domain, drifts toward DOMAIN_ONLY or
out of annotation entirely, and looks exactly like a novel dark protein. Synthesising half
a protein guarantees a dead well. S1 recovers the ones it can; whatever remains partial is
reported and set aside.
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
        sid = seq_of_orf.get(r["orf_id"])
        if sid is None:
            continue
        if r["partial"] == "0":
            complete.add(sid)
        else:
            partial_only.add(sid)
partial_only -= complete

# Spiked controls are instrumentation, not candidates. They are excluded by prefix rather
# than by any property of their annotation, so a control that failed the gate cannot leak
# into the target list.
dark = {s for s in (eligible - partial_only) if not s.startswith("CTRL_")}

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
assert n_written > 0, "the dark set is empty - check S5 flags and the cascade output"
