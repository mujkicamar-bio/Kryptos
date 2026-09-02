#!/usr/bin/env python
"""Pool the Pfam results and roll them up per family, for every threshold in the sweep.

Two sources of Pfam hits, both hmmsearch --cut_ga against Pfam-A 38.2, so directly poolable:
    data/pfam_run/domtbl/            the original run (representatives of the 30% clustering)
    data/pfam_run/recluster/domtbl/  the 75,544 representatives the sweep introduced

Best hit per sequence is the highest-scoring domain row, the same rule as
scripts/summarize_pfam.py line 76, so the pooled table matches the original convention.

Writes data/dark_orf_run/recluster/{tag}_fampfam.tsv for each threshold.
"""
import glob
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
R = ROOT / "data" / "dark_orf_run" / "recluster"


def load_domtbl(pattern):
    rows = []
    for f in sorted(glob.glob(pattern)):
        for line in open(f):
            if line.startswith("#"):
                continue
            p = line.split(None, 22)
            rows.append((p[0], int(p[2]), p[3], p[4], float(p[6]), float(p[7])))
    return pd.DataFrame(rows, columns=["target", "tlen", "pfam_name", "pfam_acc",
                                       "evalue", "score"])


def main():
    hits = pd.concat([
        load_domtbl(str(ROOT / "data" / "pfam_run" / "domtbl" / "*.domtbl")),
        load_domtbl(str(ROOT / "data" / "pfam_run" / "recluster" / "domtbl" / "*.domtbl")),
    ], ignore_index=True)
    hits["seq"] = hits["target"].str.split("::").str[1]
    hits = hits[hits["target"].str.startswith(("DARKREP::", "RECLUSTREP::"))]
    best = (hits.sort_values("score", ascending=False)
            .drop_duplicates("seq").set_index("seq"))
    print(f"pooled Pfam hits over dark representatives: {len(best):,} sequences with a hit")

    for tag in ("dark30", "dark50", "dark70", "dark90"):
        fam = pd.read_csv(R / f"{tag}_famstats.tsv", sep="\t").set_index("rep")
        fam["pfam_hit"] = fam.index.isin(best.index)
        fam["pfam_acc"] = best["pfam_acc"].reindex(fam.index)
        fam["pfam_name"] = best["pfam_name"].reindex(fam.index)
        fam.to_csv(R / f"{tag}_fampfam.tsv", sep="\t")
        print(f"  {tag}: {len(fam):,} families, {int(fam['pfam_hit'].sum()):,} Pfam-named "
              f"({fam['pfam_hit'].mean() * 100:.1f}%), covering "
              f"{fam.loc[fam['pfam_hit'], 'proteins'].sum() / fam['proteins'].sum() * 100:.1f}% "
              f"of dark ORFs")


if __name__ == "__main__":
    main()
