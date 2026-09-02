import _ctx  # noqa: F401
import csv

exclude = set(snakemake.params.exclude)
with open(snakemake.input[0], newline="") as fh, open(snakemake.output[0], "w") as out:
    for row in csv.DictReader(fh, delimiter="\t"):
        if row["hab_top"] not in exclude:
            out.write(row["plasmid_id"] + "\n")
