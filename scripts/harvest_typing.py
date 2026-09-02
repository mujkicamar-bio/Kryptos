#!/usr/bin/env python3
"""Harvest the authoritative typing runs into one per-plasmid table joinable to the master.

  MOB-suite mob_typer:  data/typing/mob/shard_*.txt  (one row per plasmid; sample_id = raw rep id)
  PlasmidFinder blastn: data/typing/plasmidfinder_blast.tsv.gz  (qseqid sseqid pident length slen)

PlasmidFinder call rule (matches PLSDB / CGE): keep a hit if pident >= 80 AND aligned length >= 60% of
the replicon-reference length; report the distinct replicon names (and Inc *families*, variants collapsed).

Output: data/plasmidscope_primary/typing_authoritative.tsv
  plasmid_id, pf_inc_types, pf_inc_families, pf_n_inc,
  mob_rep_types, mob_relaxase, mob_mpf, mob_orit, mob_mobility, mob_cluster, mob_host_range
Run: python3 scripts/harvest_typing.py
"""
import csv
import glob
import gzip
import os
import re
import sys
from collections import defaultdict

PP = "data/plasmidscope_primary"
MOB = "data/typing/mob"
PF = "data/typing/plasmidfinder_blast.tsv.gz"
OUT = f"{PP}/typing_authoritative.tsv"


def fam(x):
    return re.sub(r"\(.*?\)", "", x).strip()


def main():
    csv.field_size_limit(sys.maxsize)

    # --- PlasmidFinder: blast hits -> per-plasmid replicon calls (>=80% id, >=60% subj coverage) ---
    pf = defaultdict(set)
    if os.path.exists(PF):
        with gzip.open(PF, "rt") as fh:
            for line in fh:
                p = line.rstrip("\n").split("\t")
                if len(p) < 5:
                    continue
                q, s, pid, length, slen = p[0], p[1], float(p[2]), int(p[3]), int(p[4])
                if pid >= 80 and slen and length / slen >= 0.60:
                    # sseqid like 'IncFII(pCoo)_1__CP...' -> replicon name before the first '_<digit>'
                    name = re.split(r"_\d", s, 1)[0]
                    pf[q].add(name)
    print(f"PlasmidFinder: {len(pf):,} plasmids with >=1 replicon call", flush=True)

    # --- mob_typer: one row per plasmid ---
    mob = {}
    for fp in sorted(glob.glob(f"{MOB}/shard_*.txt")):
        with open(fp, newline="") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                sid = r.get("sample_id", "")
                if not sid:
                    continue
                def g(k):
                    v = (r.get(k, "") or "").strip()
                    return "" if v == "-" else v
                mob[sid] = dict(
                    rep=g("rep_type(s)"), relaxase=g("relaxase_type(s)"), mpf=g("mpf_type"),
                    orit=g("orit_type(s)"), mobility=g("predicted_mobility"),
                    cluster=g("primary_cluster_id"),
                    host=g("predicted_host_range_overall_name"))
    print(f"mob_typer: {len(mob):,} plasmids typed", flush=True)

    cols = ["plasmid_id", "pf_inc_types", "pf_inc_families", "pf_n_inc", "mob_rep_types",
            "mob_relaxase", "mob_mpf", "mob_orit", "mob_mobility", "mob_cluster", "mob_host_range"]
    ids = set(pf) | set(mob)
    with open(OUT, "w", newline="") as out:
        w = csv.DictWriter(out, fieldnames=cols, delimiter="\t")
        w.writeheader()
        for pid in ids:
            incs = sorted(pf.get(pid, []))
            fams = sorted(set(fam(x) for x in incs))
            m = mob.get(pid, {})
            w.writerow({"plasmid_id": pid, "pf_inc_types": ";".join(incs),
                        "pf_inc_families": ";".join(fams), "pf_n_inc": len(incs),
                        "mob_rep_types": m.get("rep", ""), "mob_relaxase": m.get("relaxase", ""),
                        "mob_mpf": m.get("mpf", ""), "mob_orit": m.get("orit", ""),
                        "mob_mobility": m.get("mobility", ""), "mob_cluster": m.get("cluster", ""),
                        "mob_host_range": m.get("host", "")})
    print(f"wrote {OUT}: {len(ids):,} plasmids (PF and/or mob typed)")


if __name__ == "__main__":
    main()
