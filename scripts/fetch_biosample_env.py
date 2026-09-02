#!/usr/bin/env python3
"""Stage 3 step 2: fetch environmental metadata for the Stage-3 target accessions from NCBI.

Method: batched efetch db=nuccore rettype=gb with seq_start=1&seq_stop=2 (truncates the sequence to
2 bp but preserves the FEATURES `source` qualifiers -- isolation_source, host, geo_loc_name/country,
lat_lon, collection_date, organism). ~150 accessions/request, ~840 requests total. Runs on the
LOGIN NODE (compute nodes have no internet). Resumable: appends to env_raw.tsv and skips accessions
already present, so it can be re-run after an interruption.

Inputs:  data/processed/plasmidscope/stage3_targets.tsv
Output:  data/processed/plasmidscope/biosample/env_raw.tsv
         accession, status, organism, isolation_source, host, geo_loc_name, lat_lon, collection_date

Run:  python3 scripts/fetch_biosample_env.py
"""
import csv
import os
import re
import sys
import time
import urllib.parse
import urllib.request

# TARGETS may be overridden with a file whose only requirement is an `accession` column, so the same
# fetcher tops up env_raw.tsv from any accession list (e.g. the PlasmidScope-primary INSDC gap list).
TARGETS = sys.argv[1] if len(sys.argv) > 1 else "data/processed/plasmidscope/stage3_targets.tsv"
OUTDIR = "data/processed/plasmidscope/biosample"
OUT = os.path.join(OUTDIR, "env_raw.tsv")
EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
EMAIL = "mujkicamar321@gmail.com"
BATCH = 150
SLEEP = 0.4           # ~2.5 req/s, under the 3/s no-key cap
FIELDS = ["organism", "isolation_source", "host", "geo_loc_name", "lat_lon", "collection_date"]
COLS = ["accession", "status"] + FIELDS


def strip_ver(acc):
    return acc.rsplit(".", 1)[0]


def load_targets():
    accs = []
    with open(TARGETS, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            accs.append(row["accession"])
    return accs


def load_done():
    done = set()
    if os.path.exists(OUT):
        with open(OUT, newline="") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                done.add(strip_ver(row["accession"]))
    return done


def efetch(ids):
    data = urllib.parse.urlencode({
        "db": "nuccore", "id": ",".join(ids), "rettype": "gb", "retmode": "text",
        "seq_start": "1", "seq_stop": "2", "email": EMAIL, "tool": "plasmidscope_enrich",
    }).encode()
    for attempt in range(4):
        try:
            with urllib.request.urlopen(EFETCH, data=data, timeout=120) as resp:
                return resp.read().decode("utf-8", "replace")
        except Exception as e:  # noqa: BLE001 -- transient network/HTTP, retry with backoff
            wait = 2 ** attempt
            print(f"  [retry {attempt+1}] {e}; sleeping {wait}s", file=sys.stderr)
            time.sleep(wait)
    raise RuntimeError(f"efetch failed after retries for {len(ids)} ids")


def fetch_records(ids):
    """Return {acc_no_version: record}. On persistent failure (e.g. one record makes NCBI truncate
    the response), split the batch to isolate the culprit; a single failing id is dropped so the
    run continues (it will be recorded with status=error by the caller)."""
    try:
        text = efetch(ids)
        return {strip_ver(r["accession"]): r for r in parse_records(text) if r["accession"]}
    except RuntimeError:
        if len(ids) == 1:
            print(f"  [drop] unfetchable accession {ids[0]}", file=sys.stderr)
            return {}
        mid = len(ids) // 2
        time.sleep(SLEEP)
        left = fetch_records(ids[:mid])
        time.sleep(SLEEP)
        right = fetch_records(ids[mid:])
        left.update(right)
        return left


_SRC_START = re.compile(r"^ {5}source ")
_FEAT_KEY = re.compile(r"^ {5}\S")
_QUAL = re.compile(r'^ {21}/(\w+)=(.*)$')


def parse_records(text):
    """Yield dict(accession, organism, isolation_source, ...) per LOCUS record in a gb response."""
    records = re.split(r"(?m)^LOCUS ", text)
    for rec in records[1:]:
        acc = None
        m = re.search(r"(?m)^VERSION\s+(\S+)", rec)
        if m:
            acc = m.group(1)
        lines = rec.splitlines()
        # isolate the source feature block
        in_src = False
        quals = {}
        cur_key = None
        for ln in lines:
            if _SRC_START.match(ln):
                in_src = True
                continue
            if in_src and _FEAT_KEY.match(ln):  # next feature -> source block ended
                break
            if in_src:
                qm = _QUAL.match(ln)
                if qm:
                    cur_key = qm.group(1)
                    val = qm.group(2).strip().strip('"')
                    quals[cur_key] = quals.get(cur_key, "") + val
                elif cur_key and ln.startswith(" " * 21):
                    quals[cur_key] += " " + ln.strip().strip('"')
        out = {"accession": acc or "", "organism": quals.get("organism", "")}
        for f in FIELDS[1:]:
            v = quals.get(f, "")
            if f == "geo_loc_name" and not v:
                v = quals.get("country", "")  # legacy qualifier name
            out[f] = v
        yield out


def main():
    os.makedirs(OUTDIR, exist_ok=True)
    targets = load_targets()
    done = load_done()
    todo = [a for a in targets if strip_ver(a) not in done]
    print(f"targets={len(targets):,} already_done={len(done):,} todo={len(todo):,}")

    new_file = not os.path.exists(OUT)
    with open(OUT, "a", newline="") as out:
        w = csv.DictWriter(out, fieldnames=COLS, delimiter="\t")
        if new_file:
            w.writeheader()
        for i in range(0, len(todo), BATCH):
            batch = todo[i:i + BATCH]
            by_acc = fetch_records(batch)
            for acc in batch:
                r = by_acc.get(strip_ver(acc))
                if r is None:
                    # not in response: either NCBI didn't return it, or it was dropped as unfetchable
                    w.writerow({"accession": acc, "status": "not_returned",
                                **{f: "" for f in FIELDS}})
                else:
                    has_env = any(r.get(f) for f in ("isolation_source", "host", "geo_loc_name"))
                    w.writerow({"accession": acc, "status": "env" if has_env else "no_env",
                                **{f: r.get(f, "") for f in FIELDS}})
            out.flush()
            if (i // BATCH) % 20 == 0:
                print(f"  {min(i+BATCH,len(todo)):,}/{len(todo):,} fetched", flush=True)
            time.sleep(SLEEP)
    print(f"done -> {OUT}")


if __name__ == "__main__":
    main()
