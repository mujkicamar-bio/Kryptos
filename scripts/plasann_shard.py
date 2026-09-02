#!/usr/bin/env python3
"""Split the working-set FASTA into NSHARDS gzipped multi-FASTA shard files (inode-frugal: NSHARDS
files, not 208k). A SLURM array task later expands its own shard into per-plasmid files on node-local
$SNIC_TMP, so the quota'd project FS never holds the 208k transient files.

Round-robin assignment (plasmid index % NSHARDS) spreads the rare megaplasmids evenly across shards.

Input:  data/plasmidscope_primary/working_set.fna.gz   (header = rep id)
Output: data/plasann_run/shards/shard_XXXX.fna.gz       (XXXX = 0..NSHARDS-1)
Run:    python3 scripts/plasann_shard.py [NSHARDS]       (default 595)
"""
import gzip
import os
import sys

NSHARDS = int(sys.argv[1]) if len(sys.argv) > 1 else 595
SRC = "data/plasmidscope_primary/working_set.fna.gz"
OUT = "data/plasann_run/shards"


def main():
    os.makedirs(OUT, exist_ok=True)
    handles = [gzip.open(f"{OUT}/shard_{i:04d}.fna.gz", "wt", compresslevel=1) for i in range(NSHARDS)]
    idx = -1
    cur = None
    counts = [0] * NSHARDS
    with gzip.open(SRC, "rt") as src:
        for line in src:
            if line.startswith(">"):
                idx += 1
                cur = handles[idx % NSHARDS]
                counts[idx % NSHARDS] += 1
                cur.write(line if line.endswith("\n") else line + "\n")
            elif cur is not None:
                cur.write(line)
    for h in handles:
        h.close()
    print(f"done: {idx + 1:,} plasmids -> {NSHARDS} shards "
          f"(min {min(counts)}, max {max(counts)} plasmids/shard)")


if __name__ == "__main__":
    main()
