#!/usr/bin/env python
"""Build MacSyFinder 'gembase' input for DefenseFinder from the PlasAnn proteome.

DefenseFinder needs genomic CONTEXT (MacSyFinder clusters co-localised genes into systems), so
this must run on whole plasmids in gene order -- not on the clustered family representatives
used for the Pfam/dbAPIS searches.

Gembase ids are parsed by macsypy as: replicon = everything before the LAST underscore, then
`.rstrip('ib')` (a Gembase quirk). Real plasmid ids contain underscores and dots and can end in
'b'/'i', so they are unsafe. Plasmids are therefore renamed PL###### (ending in a digit) with a
mapping file written alongside.

Usage: build_defensefinder_input.py [--per-chunk 2000]
"""
import argparse
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SRC = ROOT / "data" / "dark_orf_run" / "all_cds.faa"
OUT = ROOT / "data" / "defensefinder_run"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--per-chunk", type=int, default=2000, help="plasmids per chunk")
    args = ap.parse_args()

    # Read once, grouping by plasmid. Building an index up front matters: scanning all keys
    # per plasmid is O(n_plasmids * n_proteins) and does not finish on this dataset.
    from collections import defaultdict
    by_plasmid = defaultdict(list)
    pid = idx = None
    buf = []
    def flush():
        if pid is not None:
            by_plasmid[pid].append((idx, "".join(buf)))
    for line in open(SRC):
        if line.startswith(">"):
            flush()
            hdr = line[1:].strip()
            pid, idx_s, _ = hdr.rsplit("|", 2)
            idx = int(idx_s)
            buf = []
        else:
            buf.append(line.strip())
    flush()

    plasmids = sorted(by_plasmid)
    pmap = {p: f"PL{i:06d}" for i, p in enumerate(plasmids, 1)}
    (OUT / "plasmid_id_map.tsv").write_text(
        "gembase_id\tplasmid_id\n" + "".join(f"{v}\t{k}\n" for k, v in pmap.items()))

    chunks = OUT / "chunks"
    chunks.mkdir(parents=True, exist_ok=True)
    for old in chunks.glob("*.faa"):
        old.unlink()

    n_written = 0
    for c in range(0, len(plasmids), args.per_chunk):
        part = plasmids[c:c + args.per_chunk]
        with open(chunks / f"c{c // args.per_chunk:04d}.faa", "w") as fh:
            for p in part:
                # gene order = ascending CDS index, exactly as PlasAnn emitted them
                for i, (_, seq) in enumerate(sorted(by_plasmid[p]), 1):
                    fh.write(f">{pmap[p]}_{i:05d}\n{seq}\n")
                    n_written += 1
    n_chunks = (len(plasmids) + args.per_chunk - 1) // args.per_chunk
    (OUT / "manifest.txt").write_text(
        f"PLASMIDS\t{len(plasmids)}\nPROTEINS\t{n_written}\nCHUNKS\t{n_chunks}\n"
        f"PER_CHUNK\t{args.per_chunk}\n")
    print(f"{len(plasmids):,} plasmids, {n_written:,} proteins -> {n_chunks} chunks")
    print(f"array range: 0-{n_chunks - 1}")


if __name__ == "__main__":
    main()
