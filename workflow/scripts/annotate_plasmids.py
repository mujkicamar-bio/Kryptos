import _ctx  # noqa: F401
import csv

with open(snakemake.input.prot, newline="") as fh:
    annot = {r["seq_id"]: r for r in csv.DictReader(fh, delimiter="\t")}

orf_to_seq = {}
for line in open(snakemake.input.map):
    sid, ids = line.rstrip("\n").split("\t")
    for oid in ids.split(","):
        orf_to_seq[oid] = sid

cols = ["orf_id", "plasmid_id", "start", "end", "strand", "partial",
        "annot_tier", "annot_label", "functional_class", "homology_depth",
        "annot_qcov", "annot_tcov", "annot_evalue", "explained_fraction", "annot_completeness"]
n = 0
with open(snakemake.input.index, newline="") as fh, open(snakemake.output[0], "w", newline="") as out:
    w = csv.DictWriter(out, fieldnames=cols, delimiter="\t")
    w.writeheader()
    for r in csv.DictReader(fh, delimiter="\t"):
        a = annot[orf_to_seq[r["orf_id"]]]
        w.writerow({k: r.get(k, a.get(k)) for k in cols})
        n += 1
print(f"annotated {n} ORFs", file=open(snakemake.log[0], "w"))
