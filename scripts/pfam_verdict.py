#!/usr/bin/env python
"""Decide the question: is the dark plasmidome a PlasAnn artifact, or genuinely uncharacterized?

Cross-tabs the Pfam-A result against the mix30 co-clustering (which defined "dark-only") and
against MOB-lineage spread, and classifies the *rescued* families by function -- the sharpest
artifact test, because a rescued family whose Pfam domain falls inside one of PlasAnn's own
categories is an annotation failure, not novelty.

Writes data/pfam_run/pfam_verdict.txt and data/pfam_run/dark_novel_families.tsv.
"""
import re
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
P = ROOT / "data" / "pfam_run"
D = ROOT / "data" / "dark_orf_run"
MASTER = ROOT / "data" / "plasmidscope_primary" / "plasmid_metadata_master.tsv"

# Pfam name patterns that fall inside a category PlasAnn already has. A hit here means
# PlasAnn had a slot for the protein and still labelled it "Open reading frame".
PLASANN_SCOPE = [
    ("replication", r"^(Rep_|Rep\d|Rep3|RepL|RepA|Replicase|Rep_OBD|Rep_trans|RPA|RepB|Prim)"),
    ("mobilization / conjugation",
     r"^(Mob_|MobA|MobC|Relaxase|NikA|TraD|TraM|TrwB|T4SS|VirD|CagX|Conjug|MobP|MobQ)"),
    ("toxin-antitoxin", r"(toxin|antitoxin|RelB|ParE|ParD|HigA|HigB|YoeB|YafQ|VapB|VapC|MazE|MazF|Gp49)"),
    ("transposition / recombination",
     r"^(Transpos|DDE|Resolvase|Phage_integrase|Integrase|Recombinase|HTH_Tnp|rve|IS\d)"),
    ("partition / maintenance", r"^(ParA|ParB|CbiA|ParM|SopA|SopB|CcdB|StbA|Soj)"),
]


def scope_of(name):
    for label, pat in PLASANN_SCOPE:
        if re.search(pat, name, re.I):
            return label
    return "outside PlasAnn's categories"


def main():
    out = []

    def say(*a):
        line = " ".join(str(x) for x in a)
        print(line)
        out.append(line)

    per = pd.read_csv(P / "pfam_per_seq.tsv", sep="\t")
    darkhit = per[per["set"] == "DARKREP"].copy()
    darkhit["rep"] = darkhit["target"].str.split("::").str[1]
    hitset = set(darkhit["rep"])

    dc = pd.read_csv(D / "dark30_cluster.tsv", sep="\t", header=None, names=["rep", "member"])
    size = dc.groupby("rep").size()
    N_DARK = len(dc)

    # --- 1. cross-tab against the mix30 "dark-only" definition -----------------
    mx = pd.read_csv(D / "mix30_cluster.tsv", sep="\t", header=None, names=["mrep", "member"])
    mx["kind"] = mx["member"].str.rsplit("|", n=1).str[1]
    k = mx.groupby("mrep")["kind"].agg(nD=lambda s: (s == "D").sum(), nN=lambda s: (s == "N").sum())
    dark_only_reps = set(k[(k.nD > 0) & (k.nN == 0)].index)
    d_members = mx[mx["kind"].eq("D")].copy()
    d_members["dark_only"] = d_members["mrep"].isin(dark_only_reps)
    # member id in mix30 carries a |D suffix; strip it to match dark30 member ids
    d_members["dmem"] = d_members["member"].str.rsplit("|", n=1).str[0]
    dmap = dict(zip(dc["member"], dc["rep"]))
    d_members["drep"] = d_members["dmem"].map(dmap)
    d_members["pfam"] = d_members["drep"].isin(hitset)

    ct = pd.crosstab(d_members["dark_only"], d_members["pfam"])
    ct.index = ["has a PlasAnn-named homolog", "dark-only (no named homolog)"]
    ct.columns = ["no Pfam match", "Pfam match"]
    say("=== dark ORFs: named-homolog status x Pfam status ===")
    say(ct.to_string())
    say("\nas % of all dark ORFs:")
    say((ct / ct.values.sum() * 100).round(1).to_string())

    do = d_members[d_members["dark_only"]]
    say(f"\nAmong the {len(do):,} dark ORFs with NO PlasAnn-named homolog, "
        f"{do['pfam'].sum():,} ({do['pfam'].mean() * 100:.1f}%) still match a Pfam domain.")
    say(f"Genuinely unmatched by BOTH: {(~do['pfam']).sum():,} "
        f"({(~do['pfam']).sum() / N_DARK * 100:.1f}% of the dark proteome).")

    # --- 2. what did Pfam rescue? ---------------------------------------------
    darkhit["members"] = darkhit["rep"].map(size)
    darkhit["scope"] = darkhit["pfam_name"].map(scope_of)
    sc = (darkhit.groupby("scope")
          .agg(families=("rep", "size"), dark_ORFs=("members", "sum"))
          .sort_values("dark_ORFs", ascending=False))
    sc["% of rescued ORFs"] = sc["dark_ORFs"] / sc["dark_ORFs"].sum() * 100
    say("\n=== what Pfam rescued, by whether PlasAnn has a category for it ===")
    say(sc.round(1).to_string())
    inside = sc.loc[sc.index != "outside PlasAnn's categories", "dark_ORFs"].sum()
    say(f"\n{inside:,} of {sc['dark_ORFs'].sum():,} rescued dark ORFs "
        f"({inside / sc['dark_ORFs'].sum() * 100:.1f}%) carry a domain that falls squarely inside a "
        f"category PlasAnn already annotates.")
    say(f"That is {inside / N_DARK * 100:.1f}% of the whole dark proteome — annotation failure, "
        f"not novelty.")

    # --- 3. the residue: Pfam-negative families and their reach ---------------
    meta = (pd.read_csv(MASTER, sep="\t", low_memory=False,
                        usecols=["plasmid_id", "mob_cluster", "hab_sub"]).set_index("plasmid_id"))
    dc["plasmid"] = dc["member"].str.rsplit("|", n=1).str[0]
    dc = dc.join(meta, on="plasmid")
    fam = dc.groupby("rep").agg(proteins=("member", "size"), plasmids=("plasmid", "nunique"),
                                lineages=("mob_cluster", "nunique"),
                                habitats=("hab_sub", "nunique"))
    fam["pfam_hit"] = fam.index.isin(hitset)
    novel = fam[(~fam.pfam_hit) & (fam.proteins >= 10)]
    say(f"\n=== the residue: substantial families with NO Pfam match ===")
    say(f"families >=10 members and Pfam-negative: {len(novel):,} "
        f"(of {int((fam.proteins >= 10).sum()):,} substantial families)")
    say(f"  they hold {novel.proteins.sum():,} dark ORFs "
        f"({novel.proteins.sum() / N_DARK * 100:.1f}% of the dark proteome)")
    say(f"  spanning >=2 MOB lineages: {(novel.lineages >= 2).sum():,} "
        f"({(novel.lineages >= 2).mean() * 100:.1f}%)")
    say(f"  median {novel.lineages.median():.0f} lineages / {novel.habitats.median():.0f} habitats")
    say("\ntop Pfam-negative families by lineage spread:")
    say(novel.nlargest(10, "lineages")[["proteins", "plasmids", "lineages", "habitats"]]
        .to_string())

    novel.sort_values("lineages", ascending=False).to_csv(P / "dark_novel_families.tsv", sep="\t")
    (P / "pfam_verdict.txt").write_text("\n".join(out) + "\n")
    say(f"\nwrote {P/'dark_novel_families.tsv'}, {P/'pfam_verdict.txt'}")


if __name__ == "__main__":
    main()
