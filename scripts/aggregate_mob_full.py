#!/usr/bin/env python3
"""Concatenate the full MOB-suite mob_typer output (595 shard TSVs) into one compact table the typing
notebook reads — keeping the rich fields not carried into the master (mash novelty distance, host-range
rank/breadth, secondary cluster, oriT, size/gc). One row per plasmid; sample_id = raw plasmid_id.

Output: data/plasmidscope_primary/mob_full.tsv.gz
Run: python3 scripts/aggregate_mob_full.py
"""
import csv
import glob
import gzip

KEEP = ["sample_id", "size", "gc", "rep_type(s)", "relaxase_type(s)", "mpf_type", "orit_type(s)",
        "predicted_mobility", "mash_neighbor_distance", "primary_cluster_id", "secondary_cluster_id",
        "predicted_host_range_overall_rank", "predicted_host_range_overall_name"]
OUT = "data/plasmidscope_primary/mob_full.tsv.gz"


def main():
    n = 0
    with gzip.open(OUT, "wt") as out:
        w = csv.writer(out, delimiter="\t")
        w.writerow(KEEP)
        for fp in sorted(glob.glob("data/typing/mob/shard_*.txt")):
            with open(fp, newline="") as fh:
                for r in csv.DictReader(fh, delimiter="\t"):
                    if not r.get("sample_id"):
                        continue
                    w.writerow([r.get(k, "") for k in KEEP])
                    n += 1
    print(f"wrote {OUT}: {n:,} plasmids x {len(KEEP)} cols")


if __name__ == "__main__":
    main()
