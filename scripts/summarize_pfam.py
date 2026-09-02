#!/usr/bin/env python
"""Summarise the Pfam-A search and decide the dark-plasmidome novelty question.

Reads data/pfam_run/domtbl/*.domtbl (hmmsearch --cut_ga --domtblout) and writes into
data/pfam_run/:

    pfam_per_seq.tsv       one row per query sequence with a hit: best family, score, coverage
    darkfam_pfam.tsv       per dark family: rep Pfam status + size/lineage/habitat spread
    pfam_summary.txt       the printed summary

In hmmsearch domtblout the TARGET is our sequence and the QUERY is the Pfam profile.
"""
import glob
import gzip
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
P = ROOT / "data" / "pfam_run"
D = ROOT / "data" / "dark_orf_run"
PFAM = ROOT / "data" / "refs" / "pfam"
N_DARK = 378_552


def load_domtbl():
    rows = []
    for f in sorted(glob.glob(str(P / "domtbl" / "*.domtbl"))):
        for line in open(f):
            if line.startswith("#"):
                continue
            p = line.split(None, 22)
            rows.append((p[0], int(p[2]), p[3], p[4], float(p[6]), float(p[7]),
                         int(p[17]), int(p[18])))
    return pd.DataFrame(rows, columns=["target", "tlen", "pfam_name", "pfam_acc",
                                       "evalue", "score", "ali_from", "ali_to"])


def pfam_descriptions():
    desc = {}
    acc = name = None
    with gzip.open(PFAM / "Pfam-A.hmm.dat.gz", "rt") as fh:
        for line in fh:
            if line.startswith("#=GF AC"):
                acc = line.split()[-1]
            elif line.startswith("#=GF ID"):
                name = line.split(None, 2)[2].strip()
            elif line.startswith("#=GF DE") and acc:
                desc[acc.split(".")[0]] = (name, line.split(None, 2)[2].strip())
                acc = name = None
    return desc


def main():
    out = []

    def say(*a):
        line = " ".join(str(x) for x in a)
        print(line)
        out.append(line)

    man = dict(l.split("\t") for l in (P / "query_manifest.txt").read_text().splitlines())
    n_chunks_done = len(glob.glob(str(P / "domtbl" / "*.domtbl")))
    say(f"chunks parsed: {n_chunks_done}/{man['CHUNKS']}")
    if n_chunks_done < int(man["CHUNKS"]):
        say("WARNING: incomplete run — numbers below are partial")

    hits = load_domtbl()
    hits["set"] = hits["target"].str.split("::").str[0]
    hits["cov"] = (hits["ali_to"] - hits["ali_from"] + 1) / hits["tlen"]
    say(f"domain rows: {len(hits):,}   distinct Pfam families hit: {hits.pfam_acc.nunique():,}")

    # Best hit per sequence, by full-sequence score.
    best = hits.sort_values("score", ascending=False).drop_duplicates("target")

    # ---- per-set hit rates ---------------------------------------------------
    counts = {"DARKREP": int(man["DARKREP"]), "NAMED": int(man["NAMED"]),
              "VALID": int(man["VALID"])}
    say("\n=== Pfam-A hit rate by query set (--cut_ga, Pfam's curated thresholds) ===")
    tab = []
    for s, n in counts.items():
        h = int((best["set"] == s).sum())
        tab.append({"set": s, "queries": n, "with Pfam hit": h, "% hit": h / n * 100})
    say(pd.DataFrame(tab).set_index("set").round(1).to_string())

    dark_rate = tab[0]["% hit"]
    named_rate = tab[1]["% hit"]
    say(f"\nPOSITIVE CONTROL: PlasAnn-named proteins hit Pfam at {named_rate:.1f}%.")
    say(f"Dark family representatives hit at {dark_rate:.1f}%.")
    say(f"Ratio dark/named = {dark_rate / named_rate:.2f}" if named_rate else "")

    # ---- weight the dark result by family size -------------------------------
    size = pd.read_csv(D / "dark30_cluster.tsv", sep="\t", header=None,
                       names=["rep", "member"]).groupby("rep").size()
    darkhit = set(best.loc[best["set"] == "DARKREP", "target"].str.split("::").str[1])
    hit_mask = size.index.isin(darkhit)
    say(f"\n=== weighted by family size (all {N_DARK:,} dark ORFs) ===")
    say(f"dark ORFs in a family whose representative hits Pfam: "
        f"{size[hit_mask].sum():,} ({size[hit_mask].sum() / N_DARK * 100:.1f}%)")
    say(f"dark ORFs with no Pfam match at family level:         "
        f"{size[~hit_mask].sum():,} ({size[~hit_mask].sum() / N_DARK * 100:.1f}%)")
    big = size >= 10
    say(f"among the {int(big.sum()):,} families with >=10 members: "
        f"{int((big & hit_mask).sum()):,} hit Pfam ({(big & hit_mask).sum() / big.sum() * 100:.1f}%)")

    # ---- propagation validity ------------------------------------------------
    val = best[best["set"] == "VALID"].copy()
    allv = pd.Series(sorted({t for t in
                             _valid_targets(P)}))
    if len(allv):
        vrep = allv.str.split("::").str[1]
        vhit = allv.isin(val["target"])
        agree = pd.DataFrame({"rep": vrep.values, "member_hit": vhit.values})
        agree["rep_hit"] = agree["rep"].isin(darkhit)
        conc = (agree["member_hit"] == agree["rep_hit"]).mean() * 100
        say(f"\n=== propagation check ({man['VALID_FAMILIES']} families, {len(agree):,} members) ===")
        say(f"member agrees with its representative's Pfam status: {conc:.1f}%")
        byfam = agree.groupby("rep").agg(n=("member_hit", "size"),
                                         frac_hit=("member_hit", "mean"),
                                         rep_hit=("rep_hit", "first"))
        mixed = ((byfam.frac_hit > 0.05) & (byfam.frac_hit < 0.95)).sum()
        say(f"families where members disagree among themselves (5-95% hit): "
            f"{mixed} of {len(byfam)}")

    # ---- what the dark families that DO hit look like ------------------------
    desc = pfam_descriptions()
    dk = best[best["set"] == "DARKREP"].copy()
    dk["members"] = dk["target"].str.split("::").str[1].map(size)
    top = (dk.groupby(["pfam_acc", "pfam_name"])
           .agg(families=("target", "size"), dark_ORFs=("members", "sum"))
           .sort_values("dark_ORFs", ascending=False).head(25).reset_index())
    top["description"] = top["pfam_acc"].str.split(".").str[0].map(
        lambda a: desc.get(a, ("", ""))[1][:52])
    say("\n=== top Pfam families among the dark representatives that DO hit ===")
    say(top.to_string(index=False))

    say(f"\nmedian alignment coverage of the dark protein by its Pfam domain: "
        f"{dk['cov'].median() * 100:.0f}%")

    # ---- write tables --------------------------------------------------------
    best.to_csv(P / "pfam_per_seq.tsv", sep="\t", index=False)
    fam = pd.DataFrame({"members": size})
    fam["pfam_hit"] = fam.index.isin(darkhit)
    dkk = dk.set_index(dk["target"].str.split("::").str[1])
    fam["pfam_acc"] = dkk["pfam_acc"]
    fam["pfam_name"] = dkk["pfam_name"]
    fam.sort_values("members", ascending=False).to_csv(P / "darkfam_pfam.tsv", sep="\t")
    (P / "pfam_summary.txt").write_text("\n".join(out) + "\n")
    say(f"\nwrote {P/'pfam_per_seq.tsv'}, {P/'darkfam_pfam.tsv'}, {P/'pfam_summary.txt'}")


def _valid_targets(P):
    """Every VALID query we submitted (hit or not) — needed for the denominator."""
    for f in sorted(glob.glob(str(P / "chunks" / "*.faa"))):
        for line in open(f):
            if line.startswith(">VALID::"):
                yield line[1:].strip()


if __name__ == "__main__":
    main()
