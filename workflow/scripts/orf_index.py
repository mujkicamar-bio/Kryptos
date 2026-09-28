import csv

import _ctx  # noqa: F401

from plasmidann.orfindex import assign_orf_ids

orfs = []
for f in snakemake.input:
    with open(f, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            r["start"], r["end"] = int(r["start"]), int(r["end"])
            orfs.append(r)

indexed = assign_orf_ids(orfs)
# spans_origin is carried through from orf_call: a gene reconstructed across the cut point
# of a circular plasmid runs begin..length then 1..end, so a naive end - begin is negative.
cols = ["orf_id", "plasmid_id", "start", "end", "strand", "partial", "partial_begin",
        "partial_end", "spans_origin", "translation_table", "seq"]
with open(snakemake.output[0], "w", newline="") as out:
    w = csv.DictWriter(out, fieldnames=cols, delimiter="\t")
    w.writeheader()
    w.writerows(indexed)
