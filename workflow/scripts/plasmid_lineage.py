"""Stage 6: cluster plasmids by sequence similarity into lineages.

Input: the analysis-set FASTA. Mash sketches every plasmid, `mash dist` compares the sketch
database against itself, and the connected components of the linked pairs are the
lineages (plasmidann.lineage). Output plasmid_lineage.tsv: one row per plasmid with its
lineage id and the Mash distance, k-mer size and sketch size that produced it.
"""
import csv
import pathlib
import subprocess

import _ctx  # noqa: F401

from plasmidann.lineage import lineage_clusters, parse_mash_dist

cfg = snakemake.params.lineage
out = pathlib.Path(snakemake.output.tsv)
work = out.parent / "mash"
work.mkdir(parents=True, exist_ok=True)

names = []
for line in open(snakemake.input.fasta):
    if line[0] == ">":
        names.append(line[1:].split()[0])
if not names:
    raise SystemExit("plasmid_lineage: no sequences found in the input - every downstream "
                     "independence count would be empty.")

# -i sketches each sequence rather than the whole file.
sketch = work / "plasmids"
subprocess.run(
    f"mash sketch -i -k {cfg['kmer']} -s {cfg['sketch_size']} "
    f"-p {snakemake.threads} -o {sketch} {snakemake.input.fasta}",
    shell=True, check=True)

# -d applies the distance threshold at the tool, so only candidate pairs are written.
dist_out = work / "dist.tsv"
with open(dist_out, "w") as fh:
    subprocess.run(
        f"mash dist -p {snakemake.threads} -d {cfg['max_distance']} "
        f"{sketch}.msh {sketch}.msh",
        shell=True, check=True, stdout=fh)

# Streamed: with large clonal groups the all-against-all output can reach tens of GB.
with open(dist_out) as fh:
    clusters = lineage_clusters(names, parse_mash_dist(
        fh, max_distance=cfg["max_distance"], max_pvalue=cfg["max_pvalue"]))

with open(out, "w", newline="") as fh:
    writer = csv.writer(fh, delimiter="\t")
    writer.writerow(["plasmid_id", "plasmid_lineage_cluster", "lineage_max_distance",
                     "lineage_kmer", "lineage_sketch_size"])
    for plasmid_id in sorted(clusters):
        writer.writerow([plasmid_id, clusters[plasmid_id], cfg["max_distance"],
                         cfg["kmer"], cfg["sketch_size"]])

print(f"lineage: {len(names)} plasmid records -> {len(set(clusters.values()))} lineages "
      f"at mash distance <= {cfg['max_distance']}")
