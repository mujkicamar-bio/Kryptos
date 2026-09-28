"""Rule artefact_screen: flag unique proteins that are probably not proteins.

Two screens run on unique_proteins.faa:
  * AntiFam (Eberhardt et al. 2012, Database bas003), Pfam's database of families known
    to be spurious ORFs (shadow ORFs, translated RNA, repeat-derived ORFs), searched with
    hmmsearch --cut_ga, the curators' own gathering thresholds;
  * tantan: the fraction of each protein masked as low complexity, flagged above
    artefact_screen.max_low_complexity_fraction.
No other artefact test is made; in particular an ORF on the opposite strand of a real gene
that matches no AntiFam family is not flagged.

Output artefact_flags.tsv, one row per protein with the flag and its evidence. Nothing is
removed from any table. cascade_selection reads the flags: an AntiFam-flagged protein skips
every annotation tier, Tier 0 included, and does not open its family for searching; a
protein flagged for low complexity only is still searched.
"""
import csv
import pathlib
import subprocess

import _ctx  # noqa: F401

from plasmidann import scratch
from plasmidann.cascade import passes_significance
from plasmidann.fasta import iter_fasta

cfg = snakemake.params.artefact
faa = snakemake.input.faa

seq_ids = [l[1:].split()[0] for l in open(faa) if l[0] == ">"]
tmp = scratch.scratch_dir(pathlib.Path(snakemake.output[0]).parent)

# ------------------------------------------------------------------------------------
# AntiFam
# ------------------------------------------------------------------------------------
# Threshold: AntiFam's curated gathering thresholds (--cut_ga). -Z is pinned so that the
# E-values reported with each flag do not depend on the size of this input.
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

for sid, s in iter_fasta(masked):
    if s:
        masked_fraction[sid] = round(sum(c.islower() for c in s) / len(s), 4)

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

