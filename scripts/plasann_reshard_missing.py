#!/usr/bin/env python3
"""Reshard ONLY the plasmids not yet annotated (missing from data/plasann_run/annot/*.tsv.gz) into
smaller gzipped shards, so the re-run stays under the memory limit (PlasAnn's batch mode leaks
~25-30 MB/plasmid; ~120 plasmids/shard peaks ~5 GB, safe under 8 GB).

Reads:  data/plasmidscope_primary/working_set.fna.gz  +  data/plasann_run/annot/*.tsv.gz (done set)
Writes: data/plasann_run/shards2/mshard_XXXX.fna.gz
Run:    python3 scripts/plasann_reshard_missing.py [PER_SHARD]   (default 120)
"""
import csv
import glob
import gzip
import os
import sys

PER = int(sys.argv[1]) if len(sys.argv) > 1 else 120
SRC = "data/plasmidscope_primary/working_set.fna.gz"
OUTDIR = "data/plasann_run/shards2"


def main():
    csv.field_size_limit(sys.maxsize)
    done = set()
    for f in glob.glob("data/plasann_run/annot/*.tsv.gz"):
        with gzip.open(f, "rt") as fh:
            r = csv.reader(fh, delimiter="\t"); next(r, None)
            for row in r:
                if row:
                    done.add(row[0])
    print(f"already annotated: {len(done):,}", flush=True)

    os.makedirs(OUTDIR, exist_ok=True)
    # first pass: collect missing ids in file order so we can size the shard count
    missing = []
    with gzip.open(SRC, "rt") as src:
        for line in src:
            if line.startswith(">"):
                pid = line[1:].split()[0].strip()
                if pid not in done:
                    missing.append(pid)
    nsh = max(1, (len(missing) + PER - 1) // PER)
    print(f"missing: {len(missing):,} -> {nsh} shards of ~{PER}", flush=True)
    shard_of = {pid: (i // PER) for i, pid in enumerate(missing)}

    handles = {}
    cur = None
    keep = False
    with gzip.open(SRC, "rt") as src:
        for line in src:
            if line.startswith(">"):
                pid = line[1:].split()[0].strip()
                s = shard_of.get(pid)
                if s is None:
                    keep = False
                    continue
                keep = True
                if s not in handles:
                    handles[s] = gzip.open(f"{OUTDIR}/mshard_{s:04d}.fna.gz", "wt", compresslevel=1)
                cur = handles[s]
                cur.write(line if line.endswith("\n") else line + "\n")
            elif keep:
                cur.write(line)
    for h in handles.values():
        h.close()
    print(f"wrote {len(handles)} shards to {OUTDIR}")


if __name__ == "__main__":
    main()
