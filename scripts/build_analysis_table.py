#!/usr/bin/env python3
"""Build the single flat analysis table the data-investigation notebook reads.

Takes the Phase-A working set (`working_set.tsv`, 208,248 kept plasmids) and adds the two
per-plasmid fields that live only in the PlasmidScope master (`all_metadata.tsv`) and were not
carried into the working set: **GC content** and the **MOB relaxase family set**. The join key is
the representative id: every working-set `plasmid_id` is the first member of a comma-joined
provenance set in `complete_provenance.tsv`, whose full string is exactly the `Plasmid_ID` key in
`all_metadata.tsv` (verified: 208,360/208,360 match). Nothing is filtered here — this is a pure
left-join enrichment so the notebook loads one small, self-describing table.

Adds columns: gc_percent, mob_families (unique MOB* families, ';'-joined; '' = none/non-typeable),
n_mob (count of distinct families), n_source_dbs (size of the provenance set).

Inputs:  data/plasmidscope_primary/{working_set.tsv, complete_provenance.tsv, all_metadata.tsv}
Output:  data/plasmidscope_primary/analysis_table.tsv
Run:     python3 scripts/build_analysis_table.py
"""
import csv
import re
import sys

csv.field_size_limit(sys.maxsize)
PP = "data/plasmidscope_primary"
WS = f"{PP}/working_set.tsv"
PROV = f"{PP}/complete_provenance.tsv"
META = f"{PP}/all_metadata.tsv"
OUT = f"{PP}/analysis_table.tsv"

MOB_RE = re.compile(r"MOB[A-Z]")


def main():
    # rep_id -> full provenance string (the all_metadata key) + predicted_mobility (prov-only field)
    rep2full = {}
    rep2mob = {}
    with open(PROV, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            rep = r["plasmid_id"].split(",")[0]
            rep2full[rep] = r["plasmid_id"]
            rep2mob[rep] = r.get("predicted_mobility", "")

    # full string -> (gc_percent, mob_families set)
    meta = {}
    with open(META, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            try:
                gc = round(float(r["GC"]) * 100, 1)  # stored as fraction 0.1-0.9 -> percent
            except (TypeError, ValueError):
                gc = ""
            fams = sorted(set(MOB_RE.findall(r.get("MOB_type(s)", "") or "")))
            meta[r["Plasmid_ID"]] = (gc, fams)

    n = matched = 0
    with open(WS, newline="") as fh, open(OUT, "w", newline="") as out:
        r = csv.DictReader(fh, delimiter="\t")
        extra = ["predicted_mobility", "gc_percent", "mob_families", "n_mob", "n_source_dbs"]
        w = csv.DictWriter(out, fieldnames=r.fieldnames + extra, delimiter="\t")
        w.writeheader()
        for row in r:
            n += 1
            full = rep2full.get(row["plasmid_id"])
            gc, fams = meta.get(full, ("", []))
            if full in meta:
                matched += 1
            row["predicted_mobility"] = rep2mob.get(row["plasmid_id"], "")
            row["gc_percent"] = gc
            row["mob_families"] = ";".join(fams)
            row["n_mob"] = len(fams)
            row["n_source_dbs"] = len(row["sources"].split(";"))
            w.writerow(row)

    print(f"wrote {OUT}: {n:,} rows; GC/MOB matched {matched:,} ({100*matched/n:.1f}%)")


if __name__ == "__main__":
    main()
