"""S3: choose the proteins that bypass narrowing and are searched by every tier.

The cascade narrows: a protein explained up to `narrow_at` is removed from later tiers'
input, so for it there is no deeper result and the effect of the threshold cannot be seen.
For the proteins in this cohort the deeper results exist, because every tier searches them
whatever their explained fraction; comparing their results with and without the deeper
tiers measures what narrowing at `narrow_at` changes. The cohort is a uniform random sample
of `fraction` of the cascade query set (search_representatives.faa), drawn from the sorted ids, so it
is the same for the same seed and the same query set, whatever the FASTA order.
"""
import random

fraction = snakemake.params.fraction
seed = snakemake.params.seed

ids = sorted(l[1:].split()[0] for l in open(snakemake.input.faa) if l[0] == ">")

rng = random.Random(seed)
k = int(round(len(ids) * fraction))
cohort = sorted(rng.sample(ids, k)) if k else []

with open(snakemake.output[0], "w") as out:
    for i in cohort:
        out.write(i + "\n")

print(f"sweep cohort: {len(cohort)} of {len(ids)} proteins "
      f"({100 * len(cohort) / max(len(ids), 1):.2f}%), seed={seed}")
