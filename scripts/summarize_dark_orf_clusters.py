#!/usr/bin/env python
"""Summarise the dark-ORF MMseqs2 clusterings and write the result tables.

Reads the two clusterings from scripts/cluster_dark_orfs.sh and the master table,
and writes into data/dark_orf_run/:

    darkfam_stats.tsv        one row per dark-only-clustering family
    darkfam_top200.tsv       the 200 most lineage-widespread families
    dark_orf_summary.txt     the printed summary (same text as stdout)

Lineage spread uses mob_cluster (MOB-suite primary cluster) as the lineage proxy,
per reports/small_cryptic_methodology.md; habitat spread uses hab_sub.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
D = ROOT / "data" / "dark_orf_run"
MASTER = ROOT / "data" / "plasmidscope_primary" / "plasmid_metadata_master.tsv"
N_TARGET = 71_414


def rep_lengths(fasta):
    lens, cur = {}, None
    for line in open(fasta):
        if line.startswith(">"):
            cur = line[1:].strip().split()[0]
            lens[cur] = 0
        else:
            lens[cur] += len(line.strip())
    return lens


def main():
    out = []

    def say(*a):
        line = " ".join(str(x) for x in a)
        print(line)
        out.append(line)

    master = pd.read_csv(MASTER, sep="\t", low_memory=False,
                         usecols=["plasmid_id", "mob_cluster", "hab_sub", "size_bp"])
    meta = master.set_index("plasmid_id")

    # ---- dark-only clustering -------------------------------------------------
    cl = pd.read_csv(D / "dark30_cluster.tsv", sep="\t", header=None, names=["rep", "member"])
    cl["plasmid"] = cl["member"].str.rsplit("|", n=1).str[0]
    cl = cl.join(meta, on="plasmid")
    n_prot = len(cl)
    size = cl.groupby("rep").size()

    say(f"dark ORF proteins: {n_prot:,}   families: {size.size:,} "
        f"({size.size / n_prot * 100:.1f}% of proteins are their own family)")
    say(f"singleton families: {(size == 1).sum():,} ({(size == 1).mean() * 100:.1f}% of families), "
        f"holding {(size == 1).sum() / n_prot * 100:.1f}% of dark ORFs")
    say("")
    say(f"{'members':>12}  {'families':>8}  {'dark ORFs':>10}  {'% of dark':>9}")
    for t in (2, 5, 10, 50, 100, 1000):
        sel = size >= t
        say(f"{'>= ' + str(t):>12}  {sel.sum():>8,}  {size[sel].sum():>10,}  "
            f"{size[sel].sum() / n_prot * 100:>8.1f}%")

    fam = cl.groupby("rep").agg(proteins=("member", "size"), plasmids=("plasmid", "nunique"),
                                lineages=("mob_cluster", "nunique"),
                                habitats=("hab_sub", "nunique"))
    fam["rep_len_aa"] = pd.Series(rep_lengths(D / "dark30_rep_seq.fasta"))
    fam["pct_of_subset"] = fam["plasmids"] / N_TARGET * 100

    multi = fam[fam.proteins >= 2]
    say("")
    say(f"families with >=2 members: {len(multi):,}")
    say(f"  spanning >=2 MOB lineages: {(multi.lineages >= 2).sum():,} "
        f"({(multi.lineages >= 2).mean() * 100:.1f}%)")
    say(f"  spanning >=5 MOB lineages: {(multi.lineages >= 5).sum():,} "
        f"({(multi.lineages >= 5).mean() * 100:.1f}%)")
    say(f"  confined to one lineage:   {(multi.lineages == 1).sum():,} "
        f"({(multi.lineages == 1).mean() * 100:.1f}%)")
    cross = multi.loc[multi.lineages >= 2, "proteins"].sum()
    say(f"  dark ORFs in cross-lineage families: {cross:,} ({cross / n_prot * 100:.1f}% of dark ORFs)")

    big = fam[fam.proteins >= 10]
    say("")
    say(f"substantial families (>=10 members): {len(big):,}")
    say(f"  spanning >=2 lineages: {(big.lineages >= 2).sum():,} "
        f"({(big.lineages >= 2).mean() * 100:.1f}%)")
    say(f"  median lineages spanned: {big.lineages.median():.0f}   "
        f"median habitats: {big.habitats.median():.0f}")

    # ---- dark + named co-clustering ------------------------------------------
    mx = pd.read_csv(D / "mix30_cluster.tsv", sep="\t", header=None, names=["rep", "member"])
    mx["kind"] = mx["member"].str.rsplit("|", n=1).str[1]
    g = mx.groupby("rep")["kind"].agg(n="size", nD=lambda s: (s == "D").sum(),
                                      nN=lambda s: (s == "N").sum())
    mixed = g[(g.nD > 0) & (g.nN > 0)]
    pureD = g[(g.nD > 0) & (g.nN == 0)]
    say("")
    say(f"co-clustering all {len(mx):,} CDS proteins -> {len(g):,} families")
    say(f"  families with both dark and named members: {len(mixed):,}, "
        f"absorbing {mixed.nD.sum():,} dark ORFs ({mixed.nD.sum() / n_prot * 100:.1f}%)")
    say(f"  dark-only families: {len(pureD):,}, holding {pureD.nD.sum():,} dark ORFs "
        f"({pureD.nD.sum() / n_prot * 100:.1f}%)")
    say(f"  dark-only families with >=10 members: {(pureD.n >= 10).sum():,}, holding "
        f"{pureD.loc[pureD.n >= 10, 'nD'].sum():,} dark ORFs "
        f"({pureD.loc[pureD.n >= 10, 'nD'].sum() / n_prot * 100:.1f}%)")

    # ---- protein length -------------------------------------------------------
    lens = np.array([len(l.strip()) for l in open(D / "dark_orfs.faa") if not l.startswith(">")])
    say("")
    say(f"dark protein length (aa): min {lens.min()}  q25 {np.percentile(lens, 25):.0f}  "
        f"median {np.median(lens):.0f}  q75 {np.percentile(lens, 75):.0f}  max {lens.max()}")
    say(f"  under 50 aa: {(lens < 50).mean() * 100:.1f}%   under 100 aa: {(lens < 100).mean() * 100:.1f}%")

    fam.sort_values("proteins", ascending=False).to_csv(D / "darkfam_stats.tsv", sep="\t")
    fam.nlargest(200, "lineages").to_csv(D / "darkfam_top200.tsv", sep="\t")
    (D / "dark_orf_summary.txt").write_text("\n".join(out) + "\n")
    say("")
    say(f"wrote {D/'darkfam_stats.tsv'}, {D/'darkfam_top200.tsv'}, {D/'dark_orf_summary.txt'}")


if __name__ == "__main__":
    sys.exit(main())
