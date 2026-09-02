#!/usr/bin/env python
"""Aggregate the DefenseFinder / AntiDefenseFinder run over the payload-free small plasmidome.

Unlike the Pfam and dbAPIS searches, this one uses GENOMIC CONTEXT: MacSyFinder only calls a
system when its component genes co-localise, so a lone homolog does not become a "system".

Writes into data/defensefinder_run/:
    df_systems.tsv    every system called, with the real plasmid_id restored
    df_genes.tsv      every gene assigned to a system
    df_summary.txt    the printed summary
"""
import glob
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
R = ROOT / "data" / "defensefinder_run"
N_PLASMID_IN = 65_805          # plasmids with PlasAnn genes (of the 71,414 payload-free small)


def cat(kind):
    parts = []
    for f in sorted(glob.glob(str(R / "out" / "c*" / f"*_defense_finder_{kind}.tsv"))):
        try:
            d = pd.read_csv(f, sep="\t")
        except pd.errors.EmptyDataError:
            continue
        if len(d):
            parts.append(d)
    return pd.concat(parts, ignore_index=True) if parts else pd.DataFrame()


def main():
    out = []

    def say(*a):
        line = " ".join(str(x) for x in a)
        print(line)
        out.append(line)

    idmap = pd.read_csv(R / "plasmid_id_map.tsv", sep="\t")
    g2p = dict(zip(idmap["gembase_id"], idmap["plasmid_id"]))

    sysd = cat("systems")
    genes = cat("genes")
    say(f"chunks parsed: {len(glob.glob(str(R / 'out' / 'c*')))}/33")
    if sysd.empty:
        say("no systems called"); return

    # sys_id looks like PL000018_...; recover the replicon then the real plasmid id.
    sysd["gembase"] = sysd["sys_id"].str.split("_").str[0]
    sysd["plasmid_id"] = sysd["gembase"].map(g2p)
    genes["plasmid_id"] = genes["replicon"].map(g2p)

    n_pl = sysd["plasmid_id"].nunique()
    say(f"\nsystems called: {len(sysd):,} on {n_pl:,} plasmids "
        f"({n_pl / N_PLASMID_IN * 100:.2f}% of the {N_PLASMID_IN:,} searched)")

    # DefenseFinder marks anti-defense systems in the `activity` column.
    say(f"\nby activity:")
    say(sysd["activity"].value_counts(dropna=False).to_string())

    for act in sysd["activity"].dropna().unique():
        s = sysd[sysd["activity"].eq(act)]
        say(f"\n--- activity = {act}: {len(s):,} systems on "
            f"{s['plasmid_id'].nunique():,} plasmids ---")
        say(s.groupby("type").agg(systems=("sys_id", "size"),
                                  plasmids=("plasmid_id", "nunique"))
            .sort_values("systems", ascending=False).head(20).to_string())

    say(f"\nsubtypes (top 25 overall):")
    say(sysd.groupby(["activity", "subtype"]).agg(systems=("sys_id", "size"),
                                                  plasmids=("plasmid_id", "nunique"))
        .sort_values("systems", ascending=False).head(25).to_string())

    say(f"\ngenes assigned to systems: {len(genes):,} "
        f"({len(genes) / 456_949 * 100:.2f}% of the 456,949 proteins searched)")

    sysd.to_csv(R / "df_systems.tsv", sep="\t", index=False)
    genes.to_csv(R / "df_genes.tsv", sep="\t", index=False)
    (R / "df_summary.txt").write_text("\n".join(out) + "\n")
    say(f"\nwrote {R/'df_systems.tsv'}, {R/'df_genes.tsv'}, {R/'df_summary.txt'}")


if __name__ == "__main__":
    main()
