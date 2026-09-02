#!/usr/bin/env python3
"""Concatenate one shard's per-plasmid PlasAnn annotation CSVs into a single compact gzipped TSV,
prefixed with plasmid_id and with the bulky Translation column dropped. Runs inside each SLURM task on
node-local output; only the small result is written back to the project FS.

Usage: python3 plasann_harvest_shard.py <plasann_out_dir> <out.tsv.gz>
       <plasann_out_dir> holds one subfolder per plasmid, each with <pid>/<pid>_annotations.csv
"""
import csv
import glob
import gzip
import os
import sys

KEEP = ["Gene Name", "Product", "Start", "End", "Strand", "Category", "feature type"]


def main():
    out_dir, out_tsv = sys.argv[1], sys.argv[2]
    os.makedirs(os.path.dirname(out_tsv), exist_ok=True)
    n_plasmids = n_rows = 0
    with gzip.open(out_tsv, "wt", compresslevel=6) as out:
        w = csv.writer(out, delimiter="\t")
        w.writerow(["plasmid_id"] + KEEP)
        for csv_path in sorted(glob.glob(f"{out_dir}/*/*_annotations.csv")):
            pid = os.path.basename(os.path.dirname(csv_path))
            n_plasmids += 1
            with open(csv_path, newline="") as fh:
                for r in csv.DictReader(fh):
                    w.writerow([pid] + [r.get(c, "") for c in KEEP])
                    n_rows += 1
    print(f"  harvested {n_plasmids} plasmids, {n_rows} feature rows -> {out_tsv}")


if __name__ == "__main__":
    main()
