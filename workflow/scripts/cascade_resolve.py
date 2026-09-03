import _ctx  # noqa: F401
import csv, collections
from plasmidann.cascade import classify, completeness, dark_evidence, is_informative

# hits from every tier
by_query = collections.defaultdict(list)
unnamed = collections.defaultdict(list)
for f in snakemake.input.hits:
    with open(f, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            r["coverage"] = float(r["coverage"])
            r["target_coverage"] = float(r["target_coverage"]) if r["target_coverage"] else None
            by_query[r["query"]].append(r)
            if not is_informative(r["label"]):
                unnamed[r["query"]].append(r)

# explained fraction comes from the final tier's cumulative spans, which hold
# informative alignments only
explained, qlen = {}, {}
with open(snakemake.input.spans, newline="") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        explained[r["seq_id"]] = float(r["explained_fraction"])
        qlen[r["seq_id"]] = int(r["qlen"])

seq_ids = [l[1:].split()[0] for l in open(snakemake.input.faa) if l[0] == ">"]
cols = ["seq_id", "annot_tier", "annot_label", "functional_class", "homology_depth",
        "annot_qcov", "annot_tcov", "annot_evalue", "explained_fraction",
        "annot_completeness", "dark_evidence", "uninformative_labels",
        "uninformative_tiers"]

with open(snakemake.output[0], "w", newline="") as out:
    w = csv.DictWriter(out, fieldnames=cols, delimiter="\t")
    w.writeheader()
    for sid in seq_ids:
        ef = explained.get(sid, 0.0)
        u = unnamed.get(sid, [])
        labels = [x["label"] for x in u]
        w.writerow({
            "seq_id": sid,
            **classify(by_query.get(sid, [])),
            "explained_fraction": ef,
            "annot_completeness": completeness(ef),
            "dark_evidence": dark_evidence(labels),
            # the labels themselves are kept: someone else has described this protein
            "uninformative_labels": " | ".join(dict.fromkeys(labels)),
            "uninformative_tiers": ",".join(dict.fromkeys(x["tier"] for x in u)),
        })
