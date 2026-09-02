#!/usr/bin/env python3
"""Join PLSDB 2024_05_31_v2 curated metadata onto every PLSDB-sourced plasmid in our working set.

PlasmidScope carries only a thin PLSDB slice (sequence + host + mobility). PLSDB's own META archive
(downloaded to data/PLSDB/plsdb2025_meta/) adds the fields PlasmidScope lacks, all keyed by NCBI
accession (NUCCORE_ACC): **PlasmidFinder replicon/Inc typing, MOB-suite typing, AMR genes, parsed
geo-coordinates, hierarchical ecosystem tags, disease ontology, and taxonomy lineage**.

Join: for each working-set plasmid, pull the PLSDB accession(s) out of its provenance set
(`complete_provenance.tsv`, `PLSDB_`/`PLDSB_` prefix), version-strip, and look up each PLSDB table.
biosample.csv keys on a LIST-valued NUCCORE_UID (one BioSample -> several nuccore), so it is exploded
via nuccore.csv's UID->ACC map.

Output: data/plasmidscope_primary/plsdb_enrichment.tsv  (one row per enriched working-set plasmid)
Run:    python3 scripts/enrich_plsdb_metadata.py
"""
import ast
import csv
import os
import sys
from collections import defaultdict

csv.field_size_limit(sys.maxsize)
PL = "data/PLSDB/plsdb2025_meta"
PROV = "data/plasmidscope_primary/complete_provenance.tsv"
OUT = "data/plasmidscope_primary/plsdb_enrichment.tsv"


def strip_ver(a):
    return a.rsplit(".", 1)[0]


def uniq_join(xs, sep="|"):
    return sep.join(dict.fromkeys(x for x in xs if x and x not in ("-", "")))


def main():
    # --- our PLSDB accessions, and rep plasmid_id they belong to ---------------------------------
    rep_acc = {}          # rep_id -> version-stripped PLSDB accession (first one)
    for r in csv.DictReader(open(PROV, newline=""), delimiter="\t"):
        ids = r["plasmid_id"].split(",")
        rep = ids[0]
        for pid in ids:
            if pid.startswith(("PLSDB_", "PLDSB_")):
                rep_acc[rep] = strip_ver(pid.split("_", 1)[1])
                break
    ours = set(rep_acc.values())
    print(f"working-set plasmids with a PLSDB accession: {len(rep_acc):,}", flush=True)

    # --- nuccore: UID<->ACC, plus completeness/source/topology/taxUID ----------------------------
    uid2acc = {}
    nuc = {}
    for r in csv.DictReader(open(f"{PL}/nuccore.csv", newline="")):
        acc = strip_ver(r["NUCCORE_ACC"])
        uid2acc[str(r["NUCCORE_UID"])] = acc
        if acc in ours:
            nuc[acc] = {"tax_uid": r.get("TAXONOMY_UID", ""), "completeness": r.get("NUCCORE_Completeness", ""),
                        "source": r.get("NUCCORE_Source", "")}

    # --- taxonomy lineage by TAXONOMY_UID --------------------------------------------------------
    tax = {}
    for r in csv.DictReader(open(f"{PL}/taxonomy.csv", newline="")):
        tax[str(r.get("TAXONOMY_UID", ""))] = (r.get("TAXONOMY_species") or r.get("TAXONOMY_taxon_name", ""),
                                               r.get("TAXONOMY_taxon_lineage", ""))

    # --- biosample: explode list-valued NUCCORE_UID -> geo / ecosystem / disease -----------------
    bio = {}
    for r in csv.DictReader(open(f"{PL}/biosample.csv", newline="")):
        raw = r.get("NUCCORE_UID", "")
        try:
            uids = ast.literal_eval(raw) if raw.startswith("[") else [raw]
        except Exception:
            uids = [raw]
        rec = {"lat": r.get("LOCATION_lat", ""), "lng": r.get("LOCATION_lng", ""),
               "loc": r.get("LOCATION_query", ""), "eco": r.get("ECOSYSTEM_tags", ""),
               "disease": r.get("DISEASE_tags", "")}
        for u in uids:
            acc = uid2acc.get(str(u))
            if acc in ours:
                bio[acc] = rec

    # --- plasmidfinder: replicon/Inc typing (many rows/acc) --------------------------------------
    incs = defaultdict(list)
    for r in csv.DictReader(open(f"{PL}/plasmidfinder.csv", newline="")):
        acc = strip_ver(r["NUCCORE_ACC"])
        if acc in ours and r.get("typing"):
            incs[acc].append(r["typing"])

    # --- MOB-suite typing (one row/acc) ----------------------------------------------------------
    typ = {}
    for r in csv.DictReader(open(f"{PL}/typing.csv", newline="")):
        acc = strip_ver(r.get("NUCCORE_ACC", ""))
        if acc in ours:
            typ[acc] = {"rep": r.get("rep_type(s)", ""), "relaxase": r.get("relaxase_type(s)", ""),
                        "mpf": r.get("mpf_type", ""), "orit": r.get("orit_type(s)", ""),
                        "mobility": r.get("predicted_mobility", ""),
                        "mob_cluster": r.get("primary_cluster_id", ""),
                        "pmlst": r.get("PMLST_sequence_type", "")}

    # --- AMR genes (many rows/acc) ---------------------------------------------------------------
    amr_genes = defaultdict(list)
    amr_drugs = defaultdict(list)
    for r in csv.DictReader(open(f"{PL}/amr.tsv", newline=""), delimiter="\t"):
        acc = strip_ver(r.get("NUCCORE_ACC", ""))
        if acc in ours:
            if r.get("gene_symbol"):
                amr_genes[acc].append(r["gene_symbol"])
            if r.get("drug_class"):
                amr_drugs[acc].append(r["drug_class"])

    # --- write joined enrichment -----------------------------------------------------------------
    cols = ["plasmid_id", "plsdb_acc", "plsdb_completeness", "plsdb_source", "plsdb_species",
            "plsdb_inc_types", "plsdb_n_inc", "plsdb_rep_types", "plsdb_relaxase", "plsdb_mpf",
            "plsdb_orit", "plsdb_mobility", "plsdb_mob_cluster", "plsdb_pmlst",
            "plsdb_amr_genes", "plsdb_n_amr", "plsdb_drug_classes",
            "plsdb_ecosystem_tags", "plsdb_disease_tags", "plsdb_lat", "plsdb_lng", "plsdb_location"]
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    n = 0
    cov = defaultdict(int)
    with open(OUT, "w", newline="") as out:
        w = csv.DictWriter(out, fieldnames=cols, delimiter="\t")
        w.writeheader()
        for rep, acc in rep_acc.items():
            if acc not in nuc:                       # accession not in PLSDB 2024_05_31_v2
                continue
            t = typ.get(acc, {})
            b = bio.get(acc, {})
            inc = incs.get(acc, [])
            sp, _lin = tax.get(nuc[acc]["tax_uid"], ("", ""))
            row = {
                "plasmid_id": rep, "plsdb_acc": acc,
                "plsdb_completeness": nuc[acc]["completeness"], "plsdb_source": nuc[acc]["source"],
                "plsdb_species": sp,
                "plsdb_inc_types": uniq_join(inc), "plsdb_n_inc": len(set(inc)),
                "plsdb_rep_types": t.get("rep", ""), "plsdb_relaxase": t.get("relaxase", ""),
                "plsdb_mpf": t.get("mpf", ""), "plsdb_orit": t.get("orit", ""),
                "plsdb_mobility": t.get("mobility", ""), "plsdb_mob_cluster": t.get("mob_cluster", ""),
                "plsdb_pmlst": t.get("pmlst", ""),
                "plsdb_amr_genes": uniq_join(amr_genes.get(acc, [])), "plsdb_n_amr": len(set(amr_genes.get(acc, []))),
                "plsdb_drug_classes": uniq_join(amr_drugs.get(acc, [])),
                "plsdb_ecosystem_tags": b.get("eco", ""), "plsdb_disease_tags": b.get("disease", ""),
                "plsdb_lat": b.get("lat", ""), "plsdb_lng": b.get("lng", ""), "plsdb_location": b.get("loc", ""),
            }
            w.writerow(row)
            n += 1
            if inc: cov["inc_typed"] += 1
            if t.get("mobility"): cov["mob_typed"] += 1
            if amr_genes.get(acc): cov["amr"] += 1
            if b.get("eco"): cov["ecosystem"] += 1
            if b.get("lat"): cov["geo"] += 1
            if b.get("disease"): cov["disease"] += 1

    print(f"wrote {OUT}: {n:,} PLSDB-matched plasmids "
          f"({100*n/len(rep_acc):.1f}% of our PLSDB set)")
    for k in ("inc_typed", "mob_typed", "amr", "ecosystem", "geo", "disease"):
        print(f"  with {k:<10}: {cov[k]:,}")


if __name__ == "__main__":
    main()
