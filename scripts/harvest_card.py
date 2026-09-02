#!/usr/bin/env python3
"""Harvest RGI/CARD per-shard outputs into one flat hits table.

Reads every data/card_run/rgi/shard_*.txt (RGI 'main' tabular output, produced with --include_loose
so all cut-offs are present) and writes data/plasmidscope_primary/card_hits.tsv with one row per
ORF-level ARG hit, keeping the analysis-relevant columns. The RGI 'Contig' column == working-set
plasmid_id (the join key to the master). Streaming/stdlib only (frugal over 595 files).

Run:  python3 scripts/harvest_card.py
"""
import csv
import glob
import os
import sys

RGI_DIR = "data/card_run/rgi"
OUT = "data/plasmidscope_primary/card_hits.tsv"

# RGI column name -> output column name (RGI 'main' header, stable for a given RGI version).
COLS = [
    ("Contig", "plasmid_id"),
    ("ORF_ID", "orf_id"),
    ("Cut_Off", "cut_off"),
    ("Best_Hit_ARO", "best_hit_aro"),
    ("ARO", "aro"),
    ("Model_type", "model_type"),
    ("Drug Class", "drug_class"),
    ("Resistance Mechanism", "resistance_mechanism"),
    ("AMR Gene Family", "amr_gene_family"),
    ("Best_Identities", "best_identities"),
    ("Best_Hit_Bitscore", "best_hit_bitscore"),
    ("Percentage Length of Reference Sequence", "pct_length_ref"),
    ("Nudged", "nudged"),
]


def main():
    shards = sorted(glob.glob(f"{RGI_DIR}/shard_*.txt"))
    if not shards:
        sys.exit(f"no shard outputs found in {RGI_DIR}")
    src = [c[0] for c in COLS]
    n_rows = n_shards = 0
    cutoff_counts = {}
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w", newline="") as fo:
        w = csv.writer(fo, delimiter="\t")
        w.writerow([c[1] for c in COLS])
        for sh in shards:
            with open(sh, newline="") as fi:
                r = csv.DictReader(fi, delimiter="\t")
                missing = [c for c in src if c not in (r.fieldnames or [])]
                if missing:
                    sys.exit(f"{sh}: missing expected columns {missing}")
                for row in r:
                    w.writerow([row[c] for c in src])
                    n_rows += 1
                    co = row["Cut_Off"]
                    cutoff_counts[co] = cutoff_counts.get(co, 0) + 1
            n_shards += 1
    print(f"harvested {n_shards} shards -> {OUT}")
    print(f"total hit rows: {n_rows:,}")
    print("cut-off distribution:", dict(sorted(cutoff_counts.items())))


if __name__ == "__main__":
    main()
