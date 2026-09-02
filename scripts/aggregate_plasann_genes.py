#!/usr/bin/env python3
"""Aggregate every PlasAnn feature (across the 1,072 shard TSVs) into compact frequency tables the
annotation notebook reads, so the notebook never re-streams the multi-GB feature set.

Rows within each shard TSV are grouped by plasmid (harvest wrote them plasmid-by-plasmid) and a plasmid
lives in exactly one file, so consecutive-grouping gives exact per-gene plasmid counts with bounded memory.

Outputs (data/plasmidscope_primary/):
  plasann_category_totals.tsv   category, occurrences, n_plasmids
  plasann_feature_types.tsv     feature_type, occurrences
  plasann_gene_freq.tsv         gene_name, top_category, occurrences, n_plasmids   (occ >= 5)
  plasann_replicon_freq.tsv     replicon, family, n_plasmids
Run: python3 scripts/aggregate_plasann_genes.py
"""
import csv
import glob
import gzip
import os
import re
import sys
from collections import Counter, defaultdict

ANNOT = "data/plasann_run/annot"
PP = "data/plasmidscope_primary"


def main():
    csv.field_size_limit(sys.maxsize)
    # remap sanitized mMGE ids back to raw, and dedup plasmids that appear in BOTH the first-run
    # (shard_*) and retry (mshard_*) files -> count each plasmid exactly once.
    ws = [r["plasmid_id"] for r in csv.DictReader(open(f"{PP}/working_set.tsv"), delimiter="\t")]
    san2raw = {re.sub(r"[^A-Za-z0-9._-]", "_", i): i for i in ws}

    cat_occ, cat_np = Counter(), Counter()
    ft_occ = Counter()
    gene_occ, gene_np = Counter(), Counter()
    gene_cat = defaultdict(Counter)          # gene -> category votes
    rep_np = Counter()                        # replicon Gene Name -> n plasmids
    seen = set()
    n_rows = n_plasmids = 0

    def flush(cur_pid, l_cat, l_gene, l_ft, l_repcat, cats, reps):
        nonlocal n_plasmids
        if cur_pid is None:
            return
        raw = san2raw.get(cur_pid, cur_pid)
        if raw in seen:                       # already counted from another file -> skip duplicate
            return
        seen.add(raw)
        n_plasmids += 1
        cat_occ.update(l_cat); ft_occ.update(l_ft); gene_occ.update(l_gene)
        for g, cc in l_repcat.items():
            gene_cat[g].update(cc)
        for g in l_gene:
            gene_np[g] += 1
        for c in cats:
            cat_np[c] += 1
        for rp in reps:
            rep_np[rp] += 1

    files = sorted(glob.glob(f"{ANNOT}/*.tsv.gz"))
    print(f"streaming {len(files)} shard files ...", flush=True)
    for fp in files:
        cur = None
        l_cat, l_gene, l_ft, l_repcat = Counter(), Counter(), Counter(), defaultdict(Counter)
        cats, reps = set(), set()
        with gzip.open(fp, "rt") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                pid = r["plasmid_id"]
                if pid != cur:
                    flush(cur, l_cat, l_gene, l_ft, l_repcat, cats, reps)
                    cur = pid
                    l_cat, l_gene, l_ft, l_repcat = Counter(), Counter(), Counter(), defaultdict(Counter)
                    cats, reps = set(), set()
                n_rows += 1
                gn = (r.get("Gene Name") or "").strip()
                cat = (r.get("Category") or "").strip()
                ft = (r.get("feature type") or "").strip()
                if ft:
                    l_ft[ft] += 1
                if cat:
                    l_cat[cat] += 1
                    cats.add(cat)
                if gn:
                    l_gene[gn] += 1
                    if cat:
                        l_repcat[gn][cat] += 1
                if cat == "Replicon" and gn:
                    reps.add(gn)
            flush(cur, l_cat, l_gene, l_ft, l_repcat, cats, reps)
    print(f"  {n_rows:,} feature rows across {n_plasmids:,} unique plasmids", flush=True)

    with open(f"{PP}/plasann_category_totals.tsv", "w", newline="") as o:
        w = csv.writer(o, delimiter="\t"); w.writerow(["category", "occurrences", "n_plasmids"])
        for c, n in cat_occ.most_common():
            w.writerow([c, n, cat_np[c]])

    with open(f"{PP}/plasann_feature_types.tsv", "w", newline="") as o:
        w = csv.writer(o, delimiter="\t"); w.writerow(["feature_type", "occurrences"])
        for c, n in ft_occ.most_common():
            w.writerow([c, n])

    with open(f"{PP}/plasann_gene_freq.tsv", "w", newline="") as o:
        w = csv.writer(o, delimiter="\t"); w.writerow(["gene_name", "top_category", "occurrences", "n_plasmids"])
        for g, n in gene_occ.most_common():
            if n < 5:
                break
            tc = gene_cat[g].most_common(1)[0][0] if gene_cat[g] else ""
            w.writerow([g, tc, n, gene_np[g]])

    def fam(x):
        return re.sub(r"\(.*?\)", "", x).strip()

    fam_np = Counter()
    for rp, n in rep_np.items():
        fam_np[fam(rp)] += n
    with open(f"{PP}/plasann_replicon_freq.tsv", "w", newline="") as o:
        w = csv.writer(o, delimiter="\t"); w.writerow(["replicon", "family", "n_plasmids"])
        for rp, n in rep_np.most_common():
            w.writerow([rp, fam(rp), n])
    with open(f"{PP}/plasann_replicon_family_freq.tsv", "w", newline="") as o:
        w = csv.writer(o, delimiter="\t"); w.writerow(["family", "n_plasmids"])
        for f, n in fam_np.most_common():
            w.writerow([f, n])

    print("wrote category_totals, feature_types, gene_freq, replicon_freq, replicon_family_freq")


if __name__ == "__main__":
    main()
