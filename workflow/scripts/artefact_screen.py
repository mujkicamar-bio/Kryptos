"""S2b: flag proteins that are probably not proteins.

WHY THIS STAGE EXISTS

Gene callers make systematic, repeatable mistakes:

  * Shadow ORFs. On the reverse-complement strand of a real gene there is often a long
    ORF purely by chance, and it inherits coding-like statistics because it IS the reverse
    complement of genuine coding sequence.
  * Translated structural RNA. rRNA and tRNA read in some frame yield a consistent
    "protein" across thousands of genomes.
  * Repeat-derived ORFs from transposon inverted repeats and CRISPR arrays.

WHY IT MATTERS DISPROPORTIONATELY HERE

Our selection criterion is "nothing named it". A shadow ORF is BY CONSTRUCTION something
nothing named, because it is not a protein and no database contains it - so it passes the
entire cascade cleanly and lands in the dark set with a perfect score.

Worse, these artefacts are CONSERVED, because the real feature underneath them (the gene
on the opposite strand, the rRNA) is conserved. They therefore also survive the
multi-lineage and purifying-selection tests at S7. They look like ideal candidates all the
way to the plate.

AntiFam is Pfam's companion database, curated specifically as a blocklist of the artefact
families researchers kept independently rediscovering and reporting as exciting novel
conserved hypothetical proteins. It is small - 278 profiles against Pfam-A's ~21,000 - so
this costs minutes.

It does NOT save cascade compute, and an earlier version of this docstring said it did. The
cascade queries unique_proteins.faa and never reads these flags, so a flagged sequence is
still searched through every tier. That follows from P5 below: the flag is reversible and
countable precisely because nothing was removed from the searched set.

Plasmids are high-yield for these artefacts: gene-dense, GC-skewed, saturated with mobile
elements.

DESIGN PRINCIPLE P5: FLAG, NEVER DISCARD

Nothing is deleted. A flagged protein stays in every table and every count; it is excluded
from target eligibility at S5 and the exclusion is reversible and countable.
"""
import csv
import pathlib
import subprocess

import _ctx  # noqa: F401

from plasmidann import scratch
from plasmidann.cascade import passes_significance

cfg = snakemake.params.artefact
faa = snakemake.input.faa

seq_ids = [l[1:].split()[0] for l in open(faa) if l[0] == ">"]
tmp = scratch.scratch_dir(pathlib.Path(snakemake.output[0]).parent)

# ------------------------------------------------------------------------------------
# AntiFam
# ------------------------------------------------------------------------------------
# THRESHOLD: AntiFam's own curated gathering thresholds, not a blanket E-value.
#
# All 278 AntiFam profiles ship a GA line, and 274 of them (98.6%) are LOOSER than
# E=1e-5 at Z=3,497,616 - the median curated cut corresponds to E=7.3e-4, roughly 73x
# looser. An earlier version imposed -E 1e-5 here, which silently overrode the curator on
# essentially the whole database, in the one screen whose job is to stop shadow ORFs
# entering the dark set with a perfect score. It is the same argument that already sets
# --cut_ga on T1, applied where the consequence is worse: a missed artefact is not a missed
# annotation, it is a non-protein sent to the bench.
#
# -Z is still pinned for the same reason as in the cascade: the E-values REPORTED
# alongside each flag must not depend on how many sequences happened to be in this input.
antifam_hits = {}
dom = f"{tmp}/antifam.domtbl"
max_evalue = cfg.get("antifam_max_evalue")
subprocess.run(
    f"hmmsearch {cfg['antifam_args']} -Z {snakemake.params.hmmer_z} "
    f"--domZ {snakemake.params.hmmer_z} --noali --cpu {snakemake.threads} "
    f"--domtblout {dom} {cfg['antifam_db']} {faa}",
    shell=True, check=True, stdout=subprocess.DEVNULL)

for line in open(dom):
    if line.startswith("#"):
        continue
    f = line.split()
    protein, family, i_evalue = f[0], f[3], f[12]
    # None means "the tool's own curated threshold decided this", so nothing further is
    # applied. Shared with the cascade so the two gates cannot drift apart.
    if not passes_significance(i_evalue, max_evalue):
        continue
    # Keep the strongest AntiFam match, so the reason for the flag is reportable.
    if protein not in antifam_hits or float(i_evalue) < float(antifam_hits[protein][1]):
        antifam_hits[protein] = (family, i_evalue)

# ------------------------------------------------------------------------------------
# Low complexity
# ------------------------------------------------------------------------------------
# tantan masks low-complexity regions in place, writing lowercase. A protein that is mostly
# masked is a compositional artefact - homopolymer runs, simple repeats - rather than a
# folded domain, and it will match many things weakly and nothing well.
masked_fraction = {}
masked = f"{tmp}/masked.faa"
with open(masked, "w") as out:
    subprocess.run(f"tantan -p {faa}", shell=True, check=True, stdout=out)

current, lower, total = None, 0, 0
def _flush():
    if current is not None and total:
        masked_fraction[current] = round(lower / total, 4)

for line in open(masked):
    if line[0] == ">":
        _flush()
        current, lower, total = line[1:].split()[0], 0, 0
    else:
        s = line.strip()
        total += len(s)
        lower += sum(1 for c in s if c.islower())
_flush()

# ------------------------------------------------------------------------------------
# Output: one row per protein, flag plus the evidence for it
# ------------------------------------------------------------------------------------
cols = ["seq_id", "artefact_flag", "antifam_family", "antifam_ievalue",
        "low_complexity_fraction", "artefact_reason"]
n_flagged = 0
with open(snakemake.output[0], "w", newline="") as out:
    w = csv.DictWriter(out, fieldnames=cols, delimiter="\t")
    w.writeheader()
    for sid in seq_ids:
        family, ievalue = antifam_hits.get(sid, ("", ""))
        lc = masked_fraction.get(sid, 0.0)
        reasons = []
        if family:
            reasons.append("antifam")
        if lc > cfg["max_low_complexity_fraction"]:
            reasons.append("low_complexity")
        n_flagged += bool(reasons)
        w.writerow({"seq_id": sid, "artefact_flag": int(bool(reasons)),
                    "antifam_family": family, "antifam_ievalue": ievalue,
                    "low_complexity_fraction": lc,
                    "artefact_reason": ",".join(reasons)})

print(f"screened={len(seq_ids)} antifam_hits={len(antifam_hits)} "
      f"flagged={n_flagged} ({100 * n_flagged / max(len(seq_ids), 1):.3f}%)")

# The scratch directory is removed only here, on the ordinary path. A script that raised
# never reaches this line, and its intermediates are what the failure is diagnosed from.
scratch.release(tmp)

