#!/usr/bin/env python
"""Summarise the dbAPIS anti-defense search over the dark plasmidome.

dbAPIS models carry NO curated gathering thresholds (unlike Pfam), so a cutoff must be chosen.
Rather than pick one and present it as fact, results are reported across three stringency tiers
with a domain-coverage filter; the qualitative answer should not depend on the tier.

Writes data/dbapis_run/dbapis_hits.tsv and data/dbapis_run/dbapis_summary.txt.
"""
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
R = ROOT / "data" / "dbapis_run"
DB = ROOT / "data" / "refs" / "dbapis"
D = ROOT / "data" / "dark_orf_run"
N_DARK = 378_552

TIERS = [("permissive", 1e-3), ("standard", 1e-5), ("strict", 1e-10)]
MIN_COV = 0.5


def main():
    out = []

    def say(*a):
        line = " ".join(str(x) for x in a)
        print(line)
        out.append(line)

    rows = []
    for line in open(R / "dbapis.domtbl"):
        if line.startswith("#"):
            continue
        p = line.split(None, 22)
        rows.append((p[0], int(p[2]), p[3], int(p[5]), float(p[6]), float(p[11]),
                     int(p[15]), int(p[16])))
    h = pd.DataFrame(rows, columns=["target", "tlen", "model", "qlen", "evalue",
                                    "c_evalue", "hmm_from", "hmm_to"])
    h["set"] = h["target"].str.split("::").str[0]
    h["hmm_cov"] = (h["hmm_to"] - h["hmm_from"] + 1) / h["qlen"]
    h["kind"] = h["model"].str.startswith("Acr").map({True: "anti-CRISPR (Acr)",
                                                      False: "anti-defense (APIS)"})

    meta = pd.read_csv(DB / "seed_and_familyrep_all_infor.tsv", sep="\t", low_memory=False)
    sysmap = dict(zip(meta["APIS families"], meta["Defense systems"]))
    h["inhibits"] = h["model"].map(sysmap)

    size = pd.read_csv(D / "dark30_cluster.tsv", sep="\t", header=None,
                       names=["rep", "member"]).groupby("rep").size()

    say(f"dbAPIS: 290 HMMs (200 APIS + 90 Acr). {len(h):,} raw domain rows at E<=1e-3.")
    say(f"Coverage filter: HMM coverage >= {MIN_COV}.\n")

    say(f"{'tier':<12}{'E-value':>9} | {'dark families':>14}{'dark ORFs':>11}{'% of dark':>10} | "
        f"{'named prots':>12}")
    say("-" * 74)
    keep = {}
    for name, ev in TIERS:
        s = h[(h["evalue"] <= ev) & (h["hmm_cov"] >= MIN_COV)]
        dk = s[s["set"] == "DARKREP"].copy()
        dk["rep"] = dk["target"].str.split("::").str[1]
        dk = dk.drop_duplicates("rep")
        dk["members"] = dk["rep"].map(size)
        nm = s[s["set"] == "NAMED"]["target"].nunique()
        say(f"{name:<12}{ev:>9.0e} | {len(dk):>14,}{dk['members'].sum():>11,}"
            f"{dk['members'].sum() / N_DARK * 100:>9.3f}% | {nm:>12,}")
        keep[name] = (s, dk)

    s, dk = keep["standard"]
    say(f"\n=== at the standard tier (E<=1e-5, cov>={MIN_COV}) ===")
    say(f"\nby model class:")
    cls = dk.groupby("kind").agg(families=("rep", "size"), dark_ORFs=("members", "sum"))
    say(cls.to_string())

    say("\ndark families by dbAPIS model (top 20):")
    t = (dk.groupby(["model", "inhibits"])
         .agg(families=("rep", "size"), dark_ORFs=("members", "sum"))
         .sort_values("dark_ORFs", ascending=False).head(20))
    say(t.to_string())

    say("\nwhich defense system is being inhibited (dark ORF mass):")
    inh = dk.groupby(dk["inhibits"].fillna("(Acr model — CRISPR-Cas)"))["members"].sum() \
        .sort_values(ascending=False)
    say(inh.head(15).to_string())

    nmset = s[s["set"] == "NAMED"]
    say(f"\nNAMED positive control at the same tier: {nmset['target'].nunique():,} proteins, "
        f"top models:")
    say(nmset.drop_duplicates("target").groupby("model").size()
        .sort_values(ascending=False).head(10).to_string())

    keep["standard"][0].to_csv(R / "dbapis_hits.tsv", sep="\t", index=False)
    (R / "dbapis_summary.txt").write_text("\n".join(out) + "\n")
    say(f"\nwrote {R/'dbapis_hits.tsv'}, {R/'dbapis_summary.txt'}")


if __name__ == "__main__":
    main()
