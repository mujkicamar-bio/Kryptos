import _ctx  # noqa: F401
import csv, collections
from plasmidann.cascade import classify

by_query = collections.defaultdict(list)
for f in snakemake.input.hits:
    with open(f, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            r["coverage"] = float(r["coverage"])
            by_query[r["query"]].append(r)

seq_ids = [l[1:].split()[0] for l in open(snakemake.input.faa) if l[0] == ">"]
cols = ["seq_id", "annot_tier", "annot_label", "functional_class", "homology_depth"]
with open(snakemake.output[0], "w", newline="") as out:
    w = csv.DictWriter(out, fieldnames=cols, delimiter="\t")
    w.writeheader()
    for sid in seq_ids:
        w.writerow({"seq_id": sid, **classify(by_query.get(sid, []))})
