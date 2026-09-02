#!/usr/bin/env python3
"""Aggregate RGI/CARD hits to one row per plasmid (Perfect+Strict only; Loose excluded).

Input:  data/plasmidscope_primary/card_hits.tsv          (harvest_card.py, all cut-offs)
Output: data/plasmidscope_primary/card_per_plasmid.tsv   (confident AMR summary per plasmid)

Only plasmids with >=1 Perfect/Strict hit appear here; the master fold fills every other plasmid
with 0 / empty. Drug Class, Resistance Mechanism and AMR Gene Family are ';'-delimited multi-value
fields within a single hit -> split, strip, and union across all of a plasmid's hits.

Run (any python with pandas, e.g. the genesis env):  python3 scripts/aggregate_card.py
"""
import pandas as pd

HITS = "data/plasmidscope_primary/card_hits.tsv"
OUT = "data/plasmidscope_primary/card_per_plasmid.tsv"


def split_union(series):
    vals = set()
    for cell in series.dropna():
        for part in str(cell).split(";"):
            p = part.strip()
            if p and p.lower() != "n/a":
                vals.add(p)
    return sorted(vals)


def main():
    df = pd.read_csv(HITS, sep="\t", dtype=str)
    ps = df[df["cut_off"].isin(["Perfect", "Strict"])].copy()
    rows = []
    for pid, g in ps.groupby("plasmid_id"):
        aros = sorted(set(g["best_hit_aro"].dropna()))
        drug = split_union(g["drug_class"])
        mech = split_union(g["resistance_mechanism"])
        fam = split_union(g["amr_gene_family"])
        rows.append({
            "plasmid_id": pid,
            "card_n_arg": len(g),                    # Perfect+Strict ORF-level hits
            "card_n_arg_unique": len(aros),          # distinct ARO gene names
            "card_aro_list": ";".join(aros),
            "card_drug_classes": ";".join(drug),
            "card_n_drug_classes": len(drug),
            "card_resistance_mechanisms": ";".join(mech),
            "card_amr_gene_families": ";".join(fam),
            "card_multidrug": int(len(drug) >= 2),
        })
    out = pd.DataFrame(rows).sort_values("plasmid_id")
    out.to_csv(OUT, sep="\t", index=False)
    print(f"plasmids with >=1 Perfect/Strict ARG: {len(out):,}")
    print(f"wrote {OUT}")
    if len(out):
        print(f"  multidrug (>=2 drug classes): {int(out.card_multidrug.sum()):,}")
        print(f"  median ARGs among carriers:   {out.card_n_arg.median()}")
        print(f"  max ARGs on one plasmid:      {out.card_n_arg.max()}")


if __name__ == "__main__":
    main()
