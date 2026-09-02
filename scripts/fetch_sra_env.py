#!/usr/bin/env python3
"""Phase A.5 (mMGE channel). Recover environment for mMGE plasmids from their SRA run accession via
NCBI efetch db=sra (retmode=xml), parsing the SAMPLE organism + attributes. mMGE is human-microbiome
metagenomic, so the useful signal is organism (e.g. 'human gut metagenome') + free-text sample
attributes (body site, host, isolation source). Kept RAW; resumable (skips runs already fetched).

Input:  data/plasmidscope_primary/fetch_sra_runs.txt   (one SRR/DRR/ERR per line)
Output: data/plasmidscope_primary/sra_env_raw.tsv
        run, status, organism, host, isolation_source, geo_loc_name, body_site, raw_attributes
Run:  python3 scripts/fetch_sra_env.py
"""
import csv
import os
import sys
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET

RUNS = "data/plasmidscope_primary/fetch_sra_runs.txt"
OUT = "data/plasmidscope_primary/sra_env_raw.tsv"
EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
EMAIL = "mujkicamar321@gmail.com"
BATCH = int(os.environ.get("SRA_BATCH", "100"))  # smaller batches avoid silent run omissions
SLEEP = 0.4
COLS = ["run", "status", "organism", "host", "isolation_source", "geo_loc_name",
        "body_site", "raw_attributes"]
# attribute tags we surface into typed columns (everything is also kept in raw_attributes)
HOST_TAGS = {"host", "host scientific name"}
ISO_TAGS = {"isolation_source", "isolation source", "sample comment", "source", "env_material",
            "env_medium"}
GEO_TAGS = {"geo_loc_name", "geographic location", "geo loc name", "country"}
BODY_TAGS = {"body site", "body_site", "isolation source host-associated", "host body site",
             "host_body_site", "env_local_scale", "biome"}


def efetch(ids):
    data = urllib.parse.urlencode({"db": "sra", "id": ",".join(ids), "rettype": "full",
                                   "retmode": "xml", "email": EMAIL, "tool": "plasmidscope_primary"}).encode()
    for attempt in range(5):
        try:
            with urllib.request.urlopen(EFETCH, data=data, timeout=180) as resp:
                return resp.read().decode("utf-8", "replace")
        except Exception as e:  # noqa: BLE001
            wait = 2 ** attempt
            print(f"  [retry {attempt+1}] {e}; sleeping {wait}s", file=sys.stderr)
            time.sleep(wait)
    raise RuntimeError(f"sra efetch failed for {len(ids)} ids")


def parse(xml_text):
    """run_accession -> dict of env fields, from EXPERIMENT_PACKAGE_SET."""
    out = {}
    try:
        root = ET.fromstring(xml_text)
    except ET.ParseError:
        return out
    for pkg in root.findall(".//EXPERIMENT_PACKAGE"):
        sample = pkg.find(".//SAMPLE")
        organism = ""
        attrs = {}
        if sample is not None:
            sn = sample.find(".//SCIENTIFIC_NAME")
            if sn is not None and sn.text:
                organism = sn.text.strip()
            for a in sample.findall(".//SAMPLE_ATTRIBUTE"):
                t = a.findtext("TAG", "").strip()
                v = a.findtext("VALUE", "").strip()
                if t:
                    attrs[t.lower()] = v
        def pick(tags):
            for t, v in attrs.items():
                if t in tags and v:
                    return v
            return ""
        rec = {
            "organism": organism,
            "host": pick(HOST_TAGS),
            "isolation_source": pick(ISO_TAGS),
            "geo_loc_name": pick(GEO_TAGS),
            "body_site": pick(BODY_TAGS),
            "raw_attributes": "; ".join(f"{t}={v}" for t, v in attrs.items() if v)[:1000],
        }
        for run in pkg.findall(".//RUN"):
            acc = run.get("accession")
            if acc:
                out[acc] = rec
    return out


def main():
    runs = [l.strip() for l in open(RUNS) if l.strip()]
    done = set()
    if os.path.exists(OUT):
        with open(OUT, newline="") as fh:
            done = {r["run"] for r in csv.DictReader(fh, delimiter="\t")}
    todo = [r for r in runs if r not in done]
    print(f"runs={len(runs):,} done={len(done):,} todo={len(todo):,}", flush=True)

    new = not os.path.exists(OUT)
    with open(OUT, "a", newline="") as out:
        w = csv.DictWriter(out, fieldnames=COLS, delimiter="\t")
        if new:
            w.writeheader()
        for i in range(0, len(todo), BATCH):
            batch = todo[i:i + BATCH]
            recs = parse(efetch(batch))
            for run in batch:
                r = recs.get(run)
                if r is None:
                    w.writerow({"run": run, "status": "not_returned",
                                **{c: "" for c in COLS[2:]}})
                else:
                    hit = any(r[k] for k in ("host", "isolation_source", "body_site", "organism"))
                    w.writerow({"run": run, "status": "env" if hit else "no_env", **r})
            out.flush()
            if (i // BATCH) % 5 == 0:
                print(f"  {min(i+BATCH,len(todo)):,}/{len(todo):,}", flush=True)
            time.sleep(SLEEP)
    print(f"done -> {OUT}")


if __name__ == "__main__":
    main()
