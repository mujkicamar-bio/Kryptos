#!/usr/bin/env python
"""Per-family stats for one re-clustering of the dark plasmid proteome.

Mirrors the columns of data/dark_orf_run/darkfam_stats.tsv (written by
scripts/summarize_dark_orf_clusters.py) so the sweep is directly comparable to the original
30% run, but takes the clustering as an argument instead of hard-coding it.

Lineage spread uses mob_cluster (MOB-suite primary cluster) as the lineage proxy, per
reports/small_cryptic_methodology.md; habitat spread uses hab_sub.

Usage: recluster_family_stats.py CLUSTER_TSV REP_SEQ_FASTA OUT_TSV
"""
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
MASTER = ROOT / "data" / "plasmidscope_primary" / "plasmid_metadata_master.tsv"


def rep_lengths(fasta):
    lens, cur = {}, None
    for line in open(fasta):
        if line.startswith(">"):
            cur = line[1:].strip().split()[0]
            lens[cur] = 0
        else:
            lens[cur] += len(line.strip())
    return lens


def main(cluster_tsv, rep_fasta, out_tsv):
    cl = pd.read_csv(cluster_tsv, sep="\t", header=None, names=["rep", "member"])
    cl["plasmid_id"] = cl["member"].str.rsplit("|", n=1).str[0]
    met = pd.read_csv(MASTER, sep="\t", low_memory=False,
                      usecols=["plasmid_id", "mob_cluster", "hab_sub"])
    cl = cl.merge(met, on="plasmid_id", how="left")

    g = cl.groupby("rep")
    fam = pd.DataFrame({
        "proteins": g.size(),
        "plasmids": g["plasmid_id"].nunique(),
        "lineages": g["mob_cluster"].nunique(),
        "habitats": g["hab_sub"].nunique(),
    })
    fam["rep_len_aa"] = pd.Series(rep_lengths(rep_fasta))
    n_subset = cl["plasmid_id"].nunique()
    fam["pct_of_subset"] = fam["plasmids"] / n_subset * 100
    fam = fam.sort_values("proteins", ascending=False)
    fam.index.name = "rep"
    fam.to_csv(out_tsv, sep="\t")
    print(f"{Path(cluster_tsv).name}: {len(fam):,} families, {int(fam['proteins'].sum()):,} "
          f"proteins, {n_subset:,} carrier plasmids -> {out_tsv}")


if __name__ == "__main__":
    sys.exit(main(*sys.argv[1:4]))
