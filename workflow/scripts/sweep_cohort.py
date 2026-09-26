"""S3: choose the proteins that bypass narrowing and are searched by every tier.

WHY A COHORT IS NEEDED

The cascade is self-narrowing: a protein explained past `narrow_at` is physically removed
from later tiers' input. That is a compute optimisation, and it costs the counterfactual -
for a narrowed protein there is no T4 result, so "what would the answer have been at a
different threshold?" is unanswerable without a complete re-run, and each run is weeks.

Measured on v1: `min_explained` 0.3 -> 0.8 moved the deep-tier set by 76%, and the chosen
value sat exactly at the 25th percentile of the observed distribution - the densest
possible place to put a hard cut, and the one place it could not be checked.

For proteins in this cohort the counterfactual DOES exist, because they are searched by
every tier regardless of how well they were explained. That makes it possible to report,
with a confidence interval, what the threshold choice actually cost:

    "narrowing at narrow_at = 0.7 rather than searching every tier changes the deep-tier
     set by X% (95% CI ...)"

At 2% of 3.5M proteins this is ~70,000 sequences and roughly 2% of the compute.

The sample is deterministic given the seed, so the cohort is identical across re-runs and
the comparison is stable.
"""
import random

fraction = snakemake.params.fraction
seed = snakemake.params.seed

ids = [l[1:].split()[0] for l in open(snakemake.input.faa) if l[0] == ">"]

# A uniform random sample, not the first N: the FASTA is ordered by sequence hash, which
# correlates with nothing biological, but a systematic slice would still be a slice.
rng = random.Random(seed)
k = int(round(len(ids) * fraction))
cohort = sorted(rng.sample(ids, k)) if k else []

with open(snakemake.output[0], "w") as out:
    for i in cohort:
        out.write(i + "\n")

print(f"sweep cohort: {len(cohort)} of {len(ids)} proteins "
      f"({100 * len(cohort) / max(len(ids), 1):.2f}%), seed={seed}")
