#!/usr/bin/env python3
"""Recover sample geography for IMG/PR plasmids by scraping IMG/M taxon pages (authenticated).

IMG/PR gives only a `taxon_oid` (the IMG/M dataset a plasmid was found in) — no coordinates. IMG's
taxon pages DO expose Latitude / Longitude / Geographic Location / Country / Sample Collection Date
(sourced from GOLD). The 136,316 IMG/PR working-set plasmids collapse to 26,149 distinct taxa, so we
fetch each taxon once and later fan the result back out to its plasmids.

Auth: a JGI session token is required (GOLD/IMG are login-gated). Pass it via env var IMG_SESSION
(value like `/api/sessions/<hash>`); it is sent as the `jgi_session` cookie. Never hard-coded here.

Endpoint by taxon id: `3300*` -> IMG/M MetaDetail; everything else -> TaxonDetail. Resumable (skips
taxa already in the output). Polite (sleep between requests); aborts if the session looks expired.

Input:  data/GOLD/img_taxa_to_fetch.tsv   (taxon_oid, source_type, n_plasmids)
Output: data/GOLD/img_taxon_geo.tsv       (taxon_oid, page, status, lat, lng, geo_location, country,
                                            collection_date, habitat)
Run:    IMG_SESSION='/api/sessions/<hash>' python3 scripts/scrape_img_geo.py
"""
import csv
import os
import re
import sys
import time
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed

IN = os.environ.get("IMG_TAXA_IN", "data/GOLD/img_taxa_to_fetch.tsv")
OUT = "data/GOLD/img_taxon_geo.tsv"
BASE = "https://img.jgi.doe.gov/cgi-bin/mer/main.cgi"
UA = "Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 Chrome/120 Safari/537.36"
WORKERS = int(os.environ.get("IMG_WORKERS", "12"))  # IMG pages render slowly (~6s); concurrency paces
COLS = ["taxon_oid", "page", "status", "lat", "lng", "geo_location", "country",
        "collection_date", "habitat"]


def field(html, label):
    m = re.search(r">" + re.escape(label) + r"</th>\s*<td[^>]*>(.*?)</td>", html, re.S)
    if not m:
        return ""
    return re.sub(r"\s+", " ", re.sub(r"<[^>]+>", " ", m.group(1))).strip()


def num(s):
    try:
        v = float(s)
        return v if -180 <= v <= 180 else ""
    except (TypeError, ValueError):
        return ""


def fetch(session, taxon, page):
    section = "MetaDetail" if page == "metaDetail" else "TaxonDetail"
    url = f"{BASE}?section={section}&page={page}&taxon_oid={taxon}"
    req = urllib.request.Request(url, headers={"User-Agent": UA, "Cookie": f"jgi_session={session}"})
    for attempt in range(4):
        try:
            with urllib.request.urlopen(req, timeout=90) as resp:
                return resp.getcode(), resp.read().decode("utf-8", "replace")
        except urllib.error.HTTPError as e:
            return e.code, ""
        except Exception as e:  # noqa: BLE001
            time.sleep(2 ** attempt)
    return 0, ""


def main():
    session = os.environ.get("IMG_SESSION", "").strip()
    if not session:
        sys.exit("set IMG_SESSION env var to your JGI session token (/api/sessions/<hash>)")

    todo = list(csv.DictReader(open(IN, newline=""), delimiter="\t"))
    done = set()
    if os.path.exists(OUT):
        done = {r["taxon_oid"] for r in csv.DictReader(open(OUT, newline=""), delimiter="\t")}
    todo = [r for r in todo if r["taxon_oid"] not in done]
    print(f"taxa total={len(done)+len(todo):,} done={len(done):,} todo={len(todo):,} "
          f"workers={WORKERS}", flush=True)

    def work(r):
        tox = r["taxon_oid"]
        page = "metaDetail" if tox.startswith("3300") else "taxonDetail"
        code, html = fetch(session, tox, page)
        if code != 200 or ("taxon_oid" not in html and "Taxon" not in html):
            return {"taxon_oid": tox, "page": page, "status": f"http_{code}",
                    **{c: "" for c in COLS[3:]}}
        rec = {"taxon_oid": tox, "page": page, "status": "ok",
               "lat": num(field(html, "Latitude")), "lng": num(field(html, "Longitude")),
               "geo_location": field(html, "Geographic Location"),
               "country": field(html, "Country"),
               "collection_date": field(html, "Sample Collection Date"),
               "habitat": field(html, "Habitat")}
        return rec

    new = not os.path.exists(OUT)
    n_geo = n_done = recent_fail = 0
    t0 = time.time()
    with open(OUT, "a", newline="") as out, ThreadPoolExecutor(max_workers=WORKERS) as ex:
        w = csv.DictWriter(out, fieldnames=COLS, delimiter="\t")
        if new:
            w.writeheader()
        futs = {ex.submit(work, r): r for r in todo}
        for fut in as_completed(futs):
            rec = fut.result()
            w.writerow(rec)
            n_done += 1
            if rec["status"] != "ok":
                recent_fail += 1
            else:
                recent_fail = 0
                if rec["lat"] != "" or rec["geo_location"] or rec["country"]:
                    n_geo += 1
            if recent_fail >= 40 and n_geo == 0:
                out.flush()
                sys.exit("40 failures with no data — session likely expired; refresh IMG_SESSION and rerun (resumable).")
            if n_done % 100 == 0:
                out.flush()
                rate = n_done / (time.time() - t0)
                eta = (len(todo) - n_done) / rate / 60 if rate else 0
                print(f"  {n_done:,}/{len(todo):,}  with-geo={n_geo:,}  {rate:.1f}/s  eta {eta:.0f}m", flush=True)
    print(f"done: fetched {n_done:,}, with geography {n_geo:,} -> {OUT}")


if __name__ == "__main__":
    main()
