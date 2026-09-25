"""S2p: PlasmidScope's eggNOG result for every one of our proteins it contains.

See src/plasmidann/plasmidscope.py for why the class uses the eggNOG fields only and why
the join is identical sequence. The whole `ALL` table is read, not only the analysis-set
plasmids: an identical sequence has the same eggNOG result wherever it occurs.

Output: one row per unique protein of ours that has an identical PlasmidScope protein.
Proteins with no row here are not in PlasmidScope and go through the whole cascade.
"""
import csv
import gzip
import sys

import _ctx  # noqa: F401

from plasmidann.fasta import iter_fasta
from plasmidann.plasmidscope import CLASSES, FIELDS, reduce_rows

csv.field_size_limit(sys.maxsize)       # the Sequence column exceeds csv's default limit

ours = {sid for sid, _ in iter_fasta([snakemake.input.faa])}

path = snakemake.input.ps
opener = gzip.open if path.endswith(".gz") else open
with opener(path, "rt", newline="") as fh:
    records = reduce_rows(csv.DictReader(fh, delimiter="\t"), ours)

cols = ["seq_id", "ps_class", *FIELDS.values(), "orf_source", "n_orfs"]
with open(snakemake.output[0], "w", newline="") as out:
    w = csv.DictWriter(out, fieldnames=cols, delimiter="\t")
    w.writeheader()
    for sid in sorted(records):
        w.writerow(records[sid])

counts = {c: sum(r["ps_class"] == c for r in records.values()) for c in CLASSES}
print(f"plasmidscope: {len(records)} of {len(ours)} unique proteins found "
      f"({100 * len(records) / max(len(ours), 1):.1f}%); "
      + " ".join(f"{c}={n}" for c, n in counts.items()))
