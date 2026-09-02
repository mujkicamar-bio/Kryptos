import _ctx  # noqa: F401
import csv
from plasmidann.orfindex import assign_orf_ids

orfs = []
for f in snakemake.input:
    with open(f, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            r["start"], r["end"] = int(r["start"]), int(r["end"])
            orfs.append(r)

indexed = assign_orf_ids(orfs)
cols = ["orf_id", "plasmid_id", "start", "end", "strand", "partial", "seq"]
with open(snakemake.output[0], "w", newline="") as out:
    w = csv.DictWriter(out, fieldnames=cols, delimiter="\t")
    w.writeheader()
    w.writerows(indexed)
