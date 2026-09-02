#!/usr/bin/env python3
"""Aggregate the 595 per-shard PlasAnn feature TSVs into one per-plasmid annotation summary, joinable
to the master by plasmid_id.

Reads:  data/plasann_run/annot/shard_*.tsv.gz   (cols: plasmid_id, Gene Name, Product, Start, End,
        Strand, Category, feature type)
Writes: data/plasmidscope_primary/plasann_features.tsv   (one row per annotated plasmid)

PlasAnn's replicon calls are recorded here as an ANNOTATION layer (plasann_replicons); the authoritative
Inc typing comes separately from the MOB-suite/PlasmidFinder run. Coverage vs the 208,248 working set is
reported so any un-annotated plasmids (e.g. the 3 missing from the FASTA, or PlasAnn per-file failures)
are visible and can be re-run.

Run: python3 scripts/harvest_plasann.py
"""
import csv
import glob
import gzip
import os
import re
import sys
from collections import defaultdict

ANNOT = "data/plasann_run/annot"
WS = "data/plasmidscope_primary/working_set.tsv"
OUT = "data/plasmidscope_primary/plasann_features.tsv"

# Category -> summary counter column
CAT_COL = {
    "Antibiotic Resistance": "plasann_n_amr",
    "Metal and Biocide Resistance": "plasann_n_metal_biocide",
    "Conjugation": "plasann_n_conjugation",
    "Non-conjugative DNA mobility": "plasann_n_mob_dna",
    "Mobile Element": "plasann_n_mobile_element",
    "Toxin-Antitoxin System": "plasann_n_toxin_antitoxin",
    "Virulence and Defense Mechanism": "plasann_n_virulence",
    "Non coding RNA/Regulatory elements": "plasann_n_ncrna",
    "Plasmid Maintenance, Replication and Regulation": "plasann_n_maintenance",
}
COUNT_COLS = list(dict.fromkeys(CAT_COL.values()))
OUT_COLS = (["plasmid_id", "plasann_annotated", "plasann_n_features", "plasann_n_cds",
             "plasann_replicons", "plasann_n_replicons", "plasann_has_oriv", "plasann_has_orit"]
            + COUNT_COLS)


def main():
    csv.field_size_limit(sys.maxsize)
    # PlasAnn names outputs by the SANITIZED fasta filename (non [A-Za-z0-9._-] -> '_'), so mMGE ids
    # with '~' arrive sanitized. Rebuild a sanitized->raw map from the working set (injective) to
    # restore the exact plasmid_id used everywhere else.
    ws_ids = [r["plasmid_id"] for r in csv.DictReader(open(WS, newline=""), delimiter="\t")] \
        if os.path.exists(WS) else []
    san2raw = {re.sub(r"[^A-Za-z0-9._-]", "_", i): i for i in ws_ids}

    rec = defaultdict(lambda: dict(n=0, cds=0, reps=[], oriv=0, orit=0,
                                   **{c: 0 for c in COUNT_COLS}))
    files = sorted(glob.glob(f"{ANNOT}/*.tsv.gz"))   # shard_* (first run) + mshard_* (missing re-run)
    print(f"reading {len(files)} shard files ...", flush=True)
    for fp in files:
        with gzip.open(fp, "rt") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                pid = san2raw.get(r["plasmid_id"], r["plasmid_id"])
                d = rec[pid]
                d["n"] += 1
                cat = r.get("Category", "")
                ft = r.get("feature type", "")
                if ft == "CDS":
                    d["cds"] += 1
                if cat == "Replicon":
                    d["reps"].append(r.get("Gene Name", ""))
                elif cat == "Origin of Replication":
                    d["oriv"] = 1
                elif cat == "Origin of Transfer":
                    d["orit"] = 1
                col = CAT_COL.get(cat)
                if col:
                    d[col] += 1

    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w", newline="") as out:
        w = csv.DictWriter(out, fieldnames=OUT_COLS, delimiter="\t")
        w.writeheader()
        for pid, d in rec.items():
            reps = sorted(set(x for x in d["reps"] if x))
            row = {"plasmid_id": pid, "plasann_annotated": 1, "plasann_n_features": d["n"],
                   "plasann_n_cds": d["cds"], "plasann_replicons": ";".join(reps),
                   "plasann_n_replicons": len(reps), "plasann_has_oriv": d["oriv"],
                   "plasann_has_orit": d["orit"]}
            row.update({c: d[c] for c in COUNT_COLS})
            w.writerow(row)

    # coverage vs the working set
    ws = {r["plasmid_id"] for r in csv.DictReader(open(WS, newline=""), delimiter="\t")} \
        if os.path.exists(WS) else set()
    n_ann = len(rec)
    with_rep = sum(1 for d in rec.values() if d["reps"])
    print(f"wrote {OUT}: {n_ann:,} annotated plasmids")
    if ws:
        missing = ws - set(rec)
        print(f"  coverage: {n_ann:,}/{len(ws):,} ({100*n_ann/len(ws):.1f}%); "
              f"missing {len(missing):,}")
    print(f"  with >=1 replicon: {with_rep:,} ({100*with_rep/max(n_ann,1):.1f}%)")


if __name__ == "__main__":
    main()
