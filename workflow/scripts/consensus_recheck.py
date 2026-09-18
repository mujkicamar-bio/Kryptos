"""S7c: is the family collectively novel, or only individually unmatched?

THE 6.5% THIS CATCHES

A protein can miss every per-sequence threshold while its family is collectively
recognisable. The shared signal is spread across members and no single one of them carries
enough of it to clear a cut, so every member looks dark and the family looks like a
discovery. Pavlopoulos et al. removed 6.5% of their clusters exactly this way, by searching
the family CONSENSUS back against the reference databases.

For this pipeline that 6.5% is the fraction of the library that would otherwise reach the
bench described as novel when it is not. Every other test here asks about one sequence at a
time; this is the only one that asks about the family.

WHAT IT IS NOT

It is not a filter. A family whose consensus hits Pfam keeps its row, its members and its
place in the table - it gains a label, `collectively_novel = 0`, and the name of what its
consensus matched. Deciding what to do with that belongs to the report, not here.

The search uses Pfam's curated gathering thresholds, the same authority as T1. A consensus
is a synthetic sequence, so a global E-value would be even less meaningful for it than for
a real protein.
"""
import _ctx  # noqa: F401
import csv
import pathlib
import subprocess

from plasmidann import scratch

db = snakemake.params.db
ids = [l[1:].split()[0] for l in open(snakemake.input.consensus) if l[0] == ">"]

hits = {}
if ids:
    tmp = scratch.scratch_dir(pathlib.Path(snakemake.output[0]).parent)
    dom = f"{tmp}/consensus.domtbl"
    # -Z pinned to the same reference as the cascade, so an E-value reported here means the
    # same thing as one reported there even though this input is a few tens of thousands of
    # sequences rather than 3.5 million.
    subprocess.run(
        f"hmmsearch {snakemake.params.args} -Z {snakemake.params.hmmer_z} "
        f"--domZ {snakemake.params.hmmer_z} --noali --cpu {snakemake.threads} "
        f"--domtblout {dom} {db} {snakemake.input.consensus}",
        shell=True, check=True, stdout=subprocess.DEVNULL)

    for line in open(dom):
        if line.startswith("#"):
            continue
        f = line.split()
        family, label, ievalue = f[0], f[3], f[12]
        # Strongest match wins, so the label names the best explanation of the consensus.
        if family not in hits or float(ievalue) < float(hits[family][1]):
            hits[family] = (label, ievalue)

cols = ["family_id", "consensus_hit", "consensus_label", "consensus_ievalue",
        "collectively_novel"]
with open(snakemake.output[0], "w", newline="") as out:
    w = csv.DictWriter(out, fieldnames=cols, delimiter="\t")
    w.writeheader()
    for fid in ids:
        label, ievalue = hits.get(fid, ("", ""))
        w.writerow({"family_id": fid, "consensus_hit": int(bool(label)),
                    "consensus_label": label, "consensus_ievalue": ievalue,
                    "collectively_novel": int(not label)})

print(f"consensus re-check: {len(ids)} families, {len(hits)} recognised collectively "
      f"({100 * len(hits) / max(len(ids), 1):.1f}%) - labelled, not removed")

# The scratch directory is removed only here, on the ordinary path. A script that raised
# never reaches this line, and its intermediates are what the failure is diagnosed from.
# Guarded, because with no families to re-check no directory was created.
if ids:
    scratch.release(tmp)

