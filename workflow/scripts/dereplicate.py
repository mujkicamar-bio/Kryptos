import csv

import _ctx  # noqa: F401

from plasmidann.dereplicate import dereplicate

with open(snakemake.input.index, newline="") as fh:
    orfs = list(csv.DictReader(fh, delimiter="\t"))

uniques, mapping = dereplicate(orfs)
assert sum(len(v) for v in mapping.values()) == len(orfs), "dereplication lost ORFs"

with open(snakemake.output.faa, "w") as faa:
    for u in uniques:
        faa.write(f">{u['seq_id']}\n{u['seq']}\n")
with open(snakemake.output.map, "w") as m:
    for sid, ids in mapping.items():
        m.write(f"{sid}\t{','.join(ids)}\n")
