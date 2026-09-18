"""S2c: build the positive control set and spike it into the cascade input.

WHY SPIKE-IN RATHER THAN A SELF-CONTROL

An earlier version used proteins the cascade itself had labelled as backbone. That is
circular, and it fails in the one direction that matters: a protein the cascade MISSED
never enters the control set, so the control could not detect the failure it exists to
catch. It could only see a protein sliding from FUNCTIONAL to DOMAIN_ONLY.

These proteins come from Swiss-Prot - manually reviewed, of known function, and entirely
independent of anything this pipeline computed. They are added to the query set with a
CTRL_ prefix and travel through every tier exactly as a real protein does, which means the
control tests the real code path including narrowing, sharding and thresholds.

WHAT A FAILURE MEANS

If a reviewed Swiss-Prot protein of known function comes out of this cascade unnamed, then
"unnamed" carries no information, and a target list built on it would send unannotatable
noise to the bench while looking entirely reasonable.

Controls are excluded from the dark set at S6 by their prefix: they are instrumentation,
not candidates.

THE NEGATIVE CONTROLS ARE SPIKED HERE TOO

S2d builds the decoys; this stage is where both control sets enter the query. They have to
be spiked at the same point and into the same file, because the property that makes a
control a control is that it traverses the identical code path - the same narrowing, the
same shards, the same thresholds - as a real protein.
"""
import _ctx  # noqa: F401
import random

from plasmidann.cascade import is_informative

CONTROL_PREFIX = "CTRL_"

n_wanted = snakemake.params.n_controls
rng = random.Random(snakemake.params.seed)

records = []
name, buf = None, []
for line in open(snakemake.input.raw):
    if line[0] == ">":
        if name:
            records.append((name, "".join(buf)))
        name, buf = line[1:].rstrip("\n"), []
    else:
        buf.append(line.strip())
if name:
    records.append((name, "".join(buf)))

# A control must be a protein the pipeline is EXPECTED to annotate. Two filters:
#
# 1. Length and composition - no fragments, nothing absurd, no ambiguity codes.
#
# 2. The description must be informative. 14.3% of reviewed Swiss-Prot plasmid entries are
#    titled "Uncharacterized protein ..." or "UPF0102 protein ...", and cascade.is_informative
#    correctly rejects those - they are exactly the labels the pipeline is built to treat as
#    dark. Sampling them into the control set capped achievable recall at 0.862 against a
#    required 0.99, so the gate would halt EVERY run after the full four-tier cascade and
#    blame the cascade for a defect in its own instrumentation.
#
#    A protein whose own curators could not name it is not a test of our recall.
def _title(header):
    """Swiss-Prot defline: sp|ACC|ID Description OS=... - take the description."""
    tail = header.split(None, 1)[1] if " " in header else ""
    for marker in (" OS=", " OX=", " GN=", " PE=", " SV="):
        if marker in tail:
            tail = tail.split(marker)[0]
    return tail.strip()


records = [(h, s) for h, s in records
           if 50 <= len(s) <= 2000 and "X" not in s and is_informative(_title(h))]
rng.shuffle(records)
records = records[:n_wanted]

with open(snakemake.output.control, "w") as out:
    for i, (header, seq) in enumerate(records, start=1):
        # The accession is kept in the header so a failure can be looked up, but the id is
        # prefixed so every later stage can recognise a control without a lookup table.
        acc = header.split("|")[1] if "|" in header else f"u{i}"
        out.write(f">{CONTROL_PREFIX}{i:05d}_{acc}\n{seq}\n")

# The spiked file is what the cascade actually searches: the unique proteins, the positive
# controls, and the negative controls from S2d.
n_real = n_decoys = 0
with open(snakemake.output.spiked, "w") as out:
    for line in open(snakemake.input.faa):
        out.write(line)
        n_real += line[0] == ">"
    for i, (header, seq) in enumerate(records, start=1):
        acc = header.split("|")[1] if "|" in header else f"u{i}"
        out.write(f">{CONTROL_PREFIX}{i:05d}_{acc}\n{seq}\n")
    for line in open(snakemake.input.decoys):
        out.write(line)
        n_decoys += line[0] == ">"

print(f"controls={len(records)} decoys={n_decoys} spiked into {n_real} proteins "
      f"(total {n_real + len(records) + n_decoys})")
# Enough controls to measure recall at the resolution the gate demands: at
# min_control_recall = 0.99, fewer than 100 makes a single failure a 1% swing.
# From config/targets.yaml, not a literal here: this governs a run-halting gate.
MIN_CONTROLS = snakemake.params.min_controls
if len(records) < MIN_CONTROLS and n_wanted >= MIN_CONTROLS:
    raise SystemExit(
        f"only {len(records)} usable control sequences after filtering, need "
        f"{MIN_CONTROLS}. The gate cannot measure recall at 1% resolution below that. "
        "Check data/refs/control/raw.faa and the informative-title filter.")
