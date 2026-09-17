"""S3: split the cascade query set into independently searchable, independently resumable
units.

WHY THE CASCADE IS SHARDED AT ALL

T4 searches whatever the earlier tiers could not explain against nr, 375 GB. Run as a
single job that is several days of work in one unit of failure: a node eviction, a full
filesystem or a wall-clock overrun loses all of it, and the retry starts from zero. Every
other long stage in this pipeline is already sharded for exactly this reason; the cascade
was the one that was not.

A protein stays in the SAME shard for the whole cascade. That is what keeps the accounting
correct: explained spans accumulate tier by tier within a shard, and no protein is ever
compared against a partial span set assembled somewhere else.

THE COST, STATED HONESTLY

Sharding a DIAMOND search multiplies the fixed cost of reading the database by the number
of shards, because DIAMOND streams the whole thing per invocation. That is the trade-off
n_cascade_shards controls, and it is why the number is in config with the reasoning next to
it rather than hardcoded here.

Assignment is round-robin over the input order, which is deterministic and gives shards
within one sequence of each other in size. Nothing biological depends on which shard a
protein lands in.
"""
import _ctx  # noqa: F401

handles = [open(f, "w") for f in snakemake.output]
n = 0
try:
    with open(snakemake.input[0]) as fh:
        for line in fh:
            if line[0] == ">":
                out = handles[n % len(handles)]
                n += 1
            out.write(line)
finally:
    for h in handles:
        h.close()

print(f"cascade input split: {n} proteins across {len(handles)} shards")
