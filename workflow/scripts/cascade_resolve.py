import _ctx  # noqa: F401
import csv, collections
from plasmidann.cascade import classify, explained_fraction

by_query = collections.defaultdict(list)
spans = collections.defaultdict(list)
lengths = {}
for f in snakemake.input.hits:
    with open(f, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            r["coverage"] = float(r["coverage"])
            r["target_coverage"] = float(r["target_coverage"]) if r["target_coverage"] else None
            by_query[r["query"]].append(r)
            if r.get("qlen"):
                lengths[r["query"]] = int(r["qlen"])
            for iv in filter(None, (r.get("intervals") or "").split(";")):
                a, b = iv.split("-")
                spans[r["query"]].append((int(a), int(b)))

seq_ids = [l[1:].split()[0] for l in open(snakemake.input.faa) if l[0] == ">"]
cols = ["seq_id", "annot_tier", "annot_label", "functional_class", "homology_depth",
        "annot_qcov", "annot_tcov", "annot_evalue", "explained_fraction"]
with open(snakemake.output[0], "w", newline="") as out:
    w = csv.DictWriter(out, fieldnames=cols, delimiter="\t")
    w.writeheader()
    for sid in seq_ids:
        w.writerow({"seq_id": sid, **classify(by_query.get(sid, [])),
                    "explained_fraction": explained_fraction(lengths.get(sid, 0),
                                                             spans.get(sid, []))})
