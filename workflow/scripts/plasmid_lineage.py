"""Stage 6: cluster plasmids by sequence similarity into lineages (spec section 33).

Independence, not typing. See src/plasmidann/lineage.py for why this is separate from MOB
class and what single linkage costs.

Mash sketches every plasmid, compares all pairs, and the connected components below the
configured distance are the lineage clusters. `mash dist` on a sketch database is all
against all in one pass, so this is one job rather than a per-shard fan-out.
"""
import _ctx  # noqa: F401
import csv
import pathlib
import subprocess

from plasmidann.lineage import lineage_clusters, parse_mash_dist

cfg = snakemake.params.lineage
out = pathlib.Path(snakemake.output.tsv)
work = out.parent / "mash"
work.mkdir(parents=True, exist_ok=True)

# ------------------------------------------------------------------------------------
# Sketch. -i sketches each SEQUENCE separately rather than each file, which is what makes
# one shard of many plasmids into many sketches.
# ------------------------------------------------------------------------------------
sketch = work / "plasmids"
shards = " ".join(str(p) for p in snakemake.input.shards)
subprocess.run(
    f"mash sketch -i -k {cfg['kmer']} -s {cfg['sketch_size']} "
    f"-p {snakemake.threads} -o {sketch} {shards}",
    shell=True, check=True)

# ------------------------------------------------------------------------------------
# All-against-all. -d filters at the tool, so the distance threshold is applied once here
# and the same number is recorded on every row below.
# ------------------------------------------------------------------------------------
dist_out = work / "dist.tsv"
with open(dist_out, "w") as fh:
    subprocess.run(
        f"mash dist -p {snakemake.threads} -d {cfg['max_distance']} "
        f"{sketch}.msh {sketch}.msh",
        shell=True, check=True, stdout=fh)

names = []
for shard in snakemake.input.shards:
    for line in open(shard):
        if line[0] == ">":
            names.append(line[1:].split()[0])

edges = parse_mash_dist(open(dist_out).read(), max_distance=cfg["max_distance"],
                        max_pvalue=cfg["max_pvalue"])
clusters = lineage_clusters(names, edges)

with open(out, "w", newline="") as fh:
    writer = csv.writer(fh, delimiter="\t")
    # The threshold travels on every row: the cluster count is a direct function of it, so
    # any statement about independence that cites this table can be traced to the number
    # that produced it (design principle P4).
    writer.writerow(["plasmid_id", "plasmid_lineage_cluster", "lineage_max_distance",
                     "lineage_kmer", "lineage_sketch_size"])
    for plasmid_id in sorted(clusters):
        writer.writerow([plasmid_id, clusters[plasmid_id], cfg["max_distance"],
                         cfg["kmer"], cfg["sketch_size"]])

n_clusters = len(set(clusters.values()))
print(f"lineage: {len(names)} plasmids -> {n_clusters} independent lineage clusters "
      f"at mash distance <= {cfg['max_distance']}")
# Spec section 33.2: these are different biological properties, so both are reported and
# neither is presented as the other.
print(f"         plasmid records {len(names)}, independent lineages {n_clusters}")

if not names:
    raise SystemExit("plasmid_lineage: no sequences found in the shards - every downstream "
                     "independence count would be empty.")
