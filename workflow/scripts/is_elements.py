"""S8e: insertion sequence (IS) elements with ISEScan.

Transposases are already named per protein by the cascade (Pfam, pharokka, Swiss-Prot);
what nothing else provides is the ELEMENT - its boundaries, its IS family, and whether it
is complete. S8c records, as cons_is_element, whether a dark ORF overlaps an IS element
(often a degenerate transposase fragment or an IS-borne accessory gene). It is context
evidence and does not remove the ORF from the dark set.

ISEScan calls its own genes (FragGeneScan), so its ORFs are not ours. Elements are joined
to our ORFs by coordinates in S8c, exactly as integron arrays are.

Default settings, as published (Xie & Tang 2017): partial elements are reported, because
a partial IS on a plasmid is still an IS-derived region.
"""
import concurrent.futures
import csv
import pathlib
import subprocess

import _ctx  # noqa: F401

from plasmidann import scratch
from plasmidann.fasta import split_fasta
from plasmidann.isescan import COLUMNS, parse_isescan

out = pathlib.Path(snakemake.output[0])
work = scratch.scratch_dir(out.parent, "isescan_tmp")

# ISEScan's own threading kept 3.8 of 16 cores busy on the test run (1,693 CPU-s in 441 s),
# so the analysis set is dealt into one chunk per core and each chunk runs on one thread.
# ~0.8 GB per process (measured), ~75 GB at 96.
chunks = split_fasta(snakemake.input.fasta, snakemake.threads, work / "chunks")


# check=True: a tool failure must stop the run rather than produce an empty table.
def run(chunk):
    subprocess.run(
        f"isescan.py --seqfile {chunk} --output {work / chunk.stem} --nthread 1",
        shell=True, check=True, stdout=subprocess.DEVNULL)
    # ISEScan writes <output>/<input parent dir>/<input name>.tsv. A chunk with no
    # element writes no table at all, which is a legitimate result: the header-only
    # output below then says "searched, nothing found", distinct from a stage that never
    # ran.
    tables = list((work / chunk.stem).rglob(chunk.name + ".tsv"))
    if len(tables) > 1:
        raise SystemExit(f"ISEScan wrote {len(tables)} result tables for {chunk.name}, "
                         f"expected one: {tables}")
    return parse_isescan(tables[0]) if tables else []


with concurrent.futures.ThreadPoolExecutor(snakemake.threads) as pool:
    rows = [r for chunk_rows in pool.map(run, chunks) for r in chunk_rows]

with open(out, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=COLUMNS, delimiter="\t")
    w.writeheader()
    w.writerows(rows)

n_complete = sum(r["complete"] for r in rows)
print(f"is_elements: {len(rows)} IS elements ({n_complete} complete) on "
      f"{len({r['plasmid_id'] for r in rows})} plasmids")
scratch.release(work)
