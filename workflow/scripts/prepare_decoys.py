"""S2d: build the negative controls and make them part of the cascade's query set.

WHY THIS STAGE EXISTS

config/config.yaml has declared 250 shuffled and 250 reverse-complement decoys since the
controls block was written, and nothing built them. plasmidann.decoys held the
construction, tested, and no rule called it - so the pipeline ran with one of its two
controls missing and nothing in any output said so.

The missing direction is the one that matters most here. Positive controls ask "can the
cascade still recover known biology?". Negative controls ask "can it be induced to name
something that is not a protein?" - and the deliverable of this project is the DARK set,
which is defined as the complement of what the cascade could name. If that complement
admits sequences the cascade named for the wrong reason, every statement about dark
proteins inherits the error, and no positive control can detect it.

WHY THE SOURCE IS REAL CDS FROM THIS COLLECTION

A decoy drawn from a uniform random model is too easy: it fails to resemble anything for
reasons that have nothing to do with the cascade's thresholds. Both constructions start
from real plasmid CDS so that they carry this collection's composition and length
distribution, and differ from a real protein only in the one property being tested.

ORIGIN-SPANNING GENES ARE NOT USED AS SOURCES

A gene S1 reconstructed across the cut point runs start..length then 1..end, and rebuilding
it here would duplicate logic that extract_cds already owns and that is easy to get subtly
wrong. A decoy tests composition, not position, so drawing only from genes with ordinary
coordinates costs the construction nothing.

THE DECOYS TRAVEL THROUGH THE CASCADE

Like the positive controls, they are spiked into the query set and searched by every tier
under the identical thresholds. A control that took a different code path would be testing
a different pipeline from the one that produced the results.
"""
import csv
import random

import _ctx  # noqa: F401

from plasmidann import decoys
from plasmidann.fasta import iter_fasta

rng = random.Random(snakemake.params.seed)
min_length = int(snakemake.params.min_length)

# Which CDS sit on which plasmid. Ordinary coordinates only - see the module docstring.
wanted = {}
with open(snakemake.input.index, newline="") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        if r.get("spans_origin") == "1":
            continue
        wanted.setdefault(r["plasmid_id"], []).append(r)

cds_records = []
for pid, sequence in iter_fasta([snakemake.input.fasta]):
    for r in wanted.get(pid, ()):
        nt = sequence[int(r["start"]) - 1:int(r["end"])]
        if r["strand"] in ("-1", "-"):
            nt = decoys.reverse_complement(nt)
        cds_records.append((r["orf_id"], nt))

built = decoys.build_decoys(
    cds_records,
    n_shuffled=int(snakemake.params.n_shuffled),
    n_reverse_complement=int(snakemake.params.n_reverse_complement),
    rng=rng,
    min_length=min_length)

with open(snakemake.output.faa, "w") as out:
    for decoy_id, protein, _class in built:
        out.write(f">{decoy_id}\n{protein}\n")

n_shuf = sum(1 for _, _, c in built if c == "shuffled")
n_rc = len(built) - n_shuf
print(f"negative controls: {len(built)} decoys from {len(cds_records)} source CDS "
      f"({n_shuf} shuffled, {n_rc} reverse-complement)")

# Asking for decoys and getting none is the failure this stage exists to prevent, and it
# is silent: the cascade would run, the gate would report nothing, and the dark set would
# look exactly the same. Requested-but-absent is therefore an error, not a warning.
requested = int(snakemake.params.n_shuffled) + int(snakemake.params.n_reverse_complement)
if requested and not built:
    raise SystemExit(
        f"{requested} negative controls were requested and none could be built from "
        f"{len(cds_records)} source CDS. Every candidate was shorter than "
        f"min_length={min_length}, carried an ambiguity code, or had no open frame on the "
        "opposite strand. Without decoys the cascade has no false-positive measurement.")
