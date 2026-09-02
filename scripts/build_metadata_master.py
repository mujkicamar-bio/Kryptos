#!/usr/bin/env python3
"""Build the metadata MASTER table: one row per working-set plasmid, joining every metadata layer we
hold, with a maximally-supplemented geography block. No PIPdb (excluded by decision).

Layers merged (each field keeps its source in its column prefix / a *_source tag):
  base        analysis_table.tsv                    provenance, type, GC, MOB family, raw env
  reconciled  environment_reconciled_PROPOSED.tsv   hab_top/sub, sample_nature (notebook §6)
  PLSDB       plsdb_enrichment.tsv                  Inc typing, MOB, AMR, ecosystem, disease, species
  geography   env_raw + plsdb + sra (this script)   country, admin1, lat, lng, precision, source, date

Geography is coalesced per plasmid with priority:
  1. explicit point coordinates   (BioSample `lat_lon`, else PLSDB lat/lng)          -> precision=point
  2. country centroid from a place name (BioSample geo_loc_name, else PLSDB location,
     else mMGE SRA geo_loc_name); centroids are EMPIRICAL (median of the real points we
     already hold per country), with a hand fallback for tail countries               -> precision=country_centroid
  3. place name but no centroid available                                             -> precision=country_only
Country is recorded even when only a point exists (from the accompanying place name).

Output: data/plasmidscope_primary/plasmid_metadata_master.tsv
Run:    python3 scripts/build_metadata_master.py
"""
import csv
import os
import statistics
import sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from geo_utils import parse_latlon, split_geo_name, FALLBACK_CENTROIDS  # noqa: E402

csv.field_size_limit(sys.maxsize)
PP = "data/plasmidscope_primary"
ANALYSIS = f"{PP}/analysis_table.tsv"
ENV_LOCAL = f"{PP}/env_local.tsv"
ENV_RAW = "data/processed/plasmidscope/biosample/env_raw.tsv"
SRA = f"{PP}/sra_env_raw.tsv"
PLSDB = f"{PP}/plsdb_enrichment.tsv"
RECON = f"{PP}/environment_reconciled.tsv"
PLASANN = f"{PP}/plasann_features.tsv"   # PlasAnn per-plasmid annotation summary (scripts/harvest_plasann.py)
TYPING = f"{PP}/typing_authoritative.tsv"  # PlasmidFinder(>=80/60) + MOB-suite mob_typer (scripts/harvest_typing.py)
CARD = f"{PP}/card_per_plasmid.tsv"        # RGI/CARD AMR summary, Perfect+Strict (scripts/aggregate_card.py)
PROV = f"{PP}/complete_provenance.tsv"
IMGPR_DATA = "data/Cus_PR/IMG_VR_2023-08-08_1/IMGPR_plasmid_data.tsv"
IMG_GEO = "data/GOLD/img_taxon_geo.tsv"          # scraped IMG/M sample geography (taxon_oid -> lat/lng)
OUT = f"{PP}/plasmid_metadata_master.tsv"


def strip_ver(a):
    return a.rsplit(".", 1)[0]


def img_country(country_field, geo_location):
    """Best-effort country from IMG fields: explicit Country, else 'Country: ...' prefix, else last
    comma token of a 'City, State, Country' string."""
    if country_field:
        return country_field
    if not geo_location:
        return ""
    if ":" in geo_location:
        return geo_location.split(":")[0].strip()
    parts = [p.strip() for p in geo_location.split(",") if p.strip()]
    return parts[-1] if parts else ""


def main():
    # --- BioSample: accession -> geo fields --------------------------------------------------------
    bios = {}
    country_pts = defaultdict(list)
    for r in csv.DictReader(open(ENV_RAW, newline=""), delimiter="\t"):
        acc = strip_ver(r["accession"])
        bios[acc] = r
        # accumulate empirical centroid inputs
        p = parse_latlon(r.get("lat_lon", ""))
        c, _ = split_geo_name(r.get("geo_loc_name", ""))
        if p and c:
            country_pts[c].append(p)
    # empirical country centroids (median lat, median lng); require the real points we hold
    centroids = {c: (round(statistics.median(x[0] for x in ps), 3),
                     round(statistics.median(x[1] for x in ps), 3)) for c, ps in country_pts.items()}
    print(f"empirical centroids from data: {len(centroids)} countries; "
          f"hand fallback adds {len(set(FALLBACK_CENTROIDS) - set(centroids))} more", flush=True)

    def centroid(country):
        return centroids.get(country) or FALLBACK_CENTROIDS.get(country)

    # --- SRA run -> geo_loc_name -------------------------------------------------------------------
    sra = {r["run"]: r for r in csv.DictReader(open(SRA, newline=""), delimiter="\t")} if os.path.exists(SRA) else {}
    # --- per-plasmid source ids --------------------------------------------------------------------
    local = {r["rep_id"]: r for r in csv.DictReader(open(ENV_LOCAL, newline=""), delimiter="\t")}
    # --- PLSDB enrichment + reconciliation, keyed by plasmid_id ------------------------------------
    plsdb = {r["plasmid_id"]: r for r in csv.DictReader(open(PLSDB, newline=""), delimiter="\t")}
    recon = {r["plasmid_id"]: r for r in csv.DictReader(open(RECON, newline=""), delimiter="\t")}
    plasann = {r["plasmid_id"]: r for r in csv.DictReader(open(PLASANN, newline=""), delimiter="\t")} \
        if os.path.exists(PLASANN) else {}
    typing = {r["plasmid_id"]: r for r in csv.DictReader(open(TYPING, newline=""), delimiter="\t")} \
        if os.path.exists(TYPING) else {}
    card = {r["plasmid_id"]: r for r in csv.DictReader(open(CARD, newline=""), delimiter="\t")} \
        if os.path.exists(CARD) else {}

    # --- IMG/M sample geography (scraped) -> plasmid, via IMG/PR taxon_oid -------------------------
    img_geo = {r["taxon_oid"]: r for r in csv.DictReader(open(IMG_GEO, newline=""), delimiter="\t")} \
        if os.path.exists(IMG_GEO) else {}
    imgpr_taxon = {r["plasmid_id"]: r.get("taxon_oid", "")
                   for r in csv.DictReader(open(IMGPR_DATA, newline=""), delimiter="\t")}
    rep2taxon = {}
    for r in csv.DictReader(open(PROV, newline=""), delimiter="\t"):
        ids = r["plasmid_id"].split(",")
        img = next((i for i in ids if i.startswith("IMGPR_plasmid_")), "")
        if img:
            rep2taxon[ids[0]] = imgpr_taxon.get(img, "")
    print(f"IMG/M geography: {len(img_geo):,} samples; mapped to {len(rep2taxon):,} IMG/PR plasmids", flush=True)

    def geography(rep):
        """Coalesce best geography for a plasmid -> dict."""
        lc = local.get(rep, {})
        accs = [a for a in lc.get("insdc_accessions", "").split(";") if a]
        runs = [x for x in lc.get("sra_runs", "").split(";") if x]
        pl = plsdb.get(rep, {})
        # gather candidate place names (for country/admin) in priority order
        names = []
        for a in accs:
            b = bios.get(strip_ver(a))
            if b and b.get("geo_loc_name"):
                names.append(b["geo_loc_name"])
        if pl.get("plsdb_location"):
            names.append(pl["plsdb_location"])
        for run in runs:
            s = sra.get(run)
            if s and s.get("geo_loc_name"):
                names.append(s["geo_loc_name"])
        country, admin = (None, "")
        for nm in names:
            c, ad = split_geo_name(nm)
            if c:
                country, admin = c, ad
                break
        # collection date (BioSample)
        date = ""
        for a in accs:
            b = bios.get(strip_ver(a))
            if b and b.get("collection_date"):
                date = b["collection_date"]
                break
        # 1) explicit point: BioSample lat_lon first, else PLSDB lat/lng
        for a in accs:
            b = bios.get(strip_ver(a))
            p = parse_latlon(b.get("lat_lon", "")) if b else None
            if p:
                return dict(country=country or "", admin1=admin, lat=p[0], lng=p[1],
                            precision="point", source="biosample_latlon", date=date)
        if pl.get("plsdb_lat") and pl.get("plsdb_lng"):
            try:
                return dict(country=country or "", admin1=admin, lat=float(pl["plsdb_lat"]),
                            lng=float(pl["plsdb_lng"]), precision="point", source="plsdb", date=date)
            except ValueError:
                pass
        # 1b) IMG/M sample coordinates (for metagenomic IMG/PR plasmids with no INSDC geo)
        ig = img_geo.get(rep2taxon.get(rep, ""))
        if ig:
            ic = img_country(ig.get("country", ""), ig.get("geo_location", ""))
            idate = date or ig.get("collection_date", "")
            if ig.get("lat"):
                try:
                    return dict(country=ic, admin1=ig.get("geo_location", ""), lat=float(ig["lat"]),
                                lng=float(ig["lng"]), precision="point", source="img_gold", date=idate)
                except ValueError:
                    pass
            if ic:  # IMG place name but no coordinate
                ct = centroid(ic)
                return dict(country=ic, admin1=ig.get("geo_location", ""),
                            lat=ct[0] if ct else "", lng=ct[1] if ct else "",
                            precision="country_centroid" if ct else "country_only",
                            source="img_gold", date=idate)
        # 2) country centroid from a place name
        if country:
            ct = centroid(country)
            if ct:
                return dict(country=country, admin1=admin, lat=ct[0], lng=ct[1],
                            precision="country_centroid", source="geoname_centroid", date=date)
            return dict(country=country, admin1=admin, lat="", lng="", precision="country_only",
                        source="geoname", date=date)
        return dict(country="", admin1="", lat="", lng="", precision="", source="", date=date)

    base_cols = ["plasmid_id", "sources", "topology", "size_bp", "gc_percent",
                 "predicted_mobility", "mob_families", "n_source_dbs", "env_status"]
    recon_cols = ["hab_top", "hab_sub", "hab_channel", "is_clinical"]
    plsdb_cols = ["plsdb_acc", "plsdb_inc_types", "plsdb_n_inc", "plsdb_rep_types", "plsdb_relaxase",
                  "plsdb_mobility", "plsdb_mob_cluster", "plsdb_amr_genes", "plsdb_n_amr",
                  "plsdb_drug_classes", "plsdb_ecosystem_tags", "plsdb_disease_tags", "plsdb_species"]
    geo_cols = ["geo_country", "geo_admin1", "geo_lat", "geo_lng", "geo_precision", "geo_source",
                "collection_date"]
    plasann_cols = ["plasann_annotated", "plasann_n_features", "plasann_n_cds", "plasann_replicons",
                    "plasann_n_replicons", "plasann_has_oriv", "plasann_has_orit", "plasann_n_amr",
                    "plasann_n_metal_biocide", "plasann_n_conjugation", "plasann_n_mob_dna",
                    "plasann_n_mobile_element", "plasann_n_toxin_antitoxin", "plasann_n_virulence",
                    "plasann_n_ncrna", "plasann_n_maintenance"]
    typing_cols = ["pf_inc_types", "pf_inc_families", "pf_n_inc", "mob_rep_types", "mob_relaxase",
                   "mob_mpf", "mob_orit", "mob_mobility", "mob_cluster", "mob_host_range"]
    card_cols = ["card_n_arg", "card_n_arg_unique", "card_aro_list", "card_drug_classes",
                 "card_n_drug_classes", "card_resistance_mechanisms", "card_amr_gene_families",
                 "card_multidrug"]
    card_num = {"card_n_arg", "card_n_arg_unique", "card_n_drug_classes", "card_multidrug"}
    out_cols = base_cols + recon_cols + plsdb_cols + geo_cols + plasann_cols + typing_cols + card_cols

    n = 0
    cov = defaultdict(int)
    prec = defaultdict(int)
    with open(ANALYSIS, newline="") as fh, open(OUT, "w", newline="") as out:
        r = csv.DictReader(fh, delimiter="\t")
        w = csv.DictWriter(out, fieldnames=out_cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for row in r:
            rep = row["plasmid_id"]
            g = geography(rep)
            rc = recon.get(rep, {})
            pl = plsdb.get(rep, {})
            # Simulated-communities artifact samples carry the JGI compute-site coordinate (Berkeley),
            # which is NOT a real location -> flag so it is never mapped/counted as real geography.
            # excluded (non-habitat) buckets: simulated-community artifact + lab-artifact (lab
            # enrichment / content-free engineered). Both are dropped from the real-geography stats.
            artifact = rc.get("hab_top") in ("Simulated-artifact", "Lab-artifact")
            if rc.get("hab_top") == "Simulated-artifact" and g["source"] == "img_gold":
                g = dict(g, source="img_gold_simulated", precision="artifact")
            out_row = {k: row.get(k, "") for k in base_cols}
            out_row.update({k: rc.get(k, "") for k in recon_cols})
            out_row.update({k: pl.get(k, "") for k in plsdb_cols})
            out_row.update({"geo_country": g["country"], "geo_admin1": g["admin1"],
                            "geo_lat": g["lat"], "geo_lng": g["lng"], "geo_precision": g["precision"],
                            "geo_source": g["source"], "collection_date": g["date"]})
            pa = plasann.get(rep, {})
            out_row.update({k: pa.get(k, "0" if k == "plasann_annotated" else "") for k in plasann_cols})
            ty = typing.get(rep, {})
            out_row.update({k: ty.get(k, "") for k in typing_cols})
            ca = card.get(rep, {})
            out_row.update({k: ca.get(k, "0" if k in card_num else "") for k in card_cols})
            w.writerow(out_row)
            n += 1
            if artifact:
                cov["artifact"] += 1
            if g["country"] and not artifact:
                cov["country_real"] += 1
            if g["lat"] != "" and not artifact:
                cov["coords_real"] += 1
            if g["date"] and not artifact:
                cov["date_real"] += 1
            prec[g["precision"] or "none"] += 1

    non_art = n - cov["artifact"]
    print(f"wrote {OUT}: {n:,} plasmids ({cov['artifact']:,} Simulated-communities artifact, "
          f"geography flagged/excluded)")
    print(f"  REAL geography (of {non_art:,} non-artifact plasmids): "
          f"country {cov['country_real']:,} ({100*cov['country_real']/non_art:.1f}%); "
          f"coordinates {cov['coords_real']:,} ({100*cov['coords_real']/non_art:.1f}%); "
          f"collection_date {cov['date_real']:,} ({100*cov['date_real']/non_art:.1f}%)")
    print("  precision:", {k: f"{v:,}" for k, v in sorted(prec.items(), key=lambda x: -x[1])})


if __name__ == "__main__":
    main()
