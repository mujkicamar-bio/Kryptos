# The 100-plasmid smoke set

A small, seeded, stratified plasmid set for checking that the pipeline runs end to end and
produces every output table with its declared columns. It is not a benchmark and it measures no
biology: 100 plasmids cannot estimate a dark fraction, a background rate or a recall.

## Building it

```bash
envs/plasmidann/bin/python tools/make_test_set.py \
  --master data/plasmidscope_primary/analysis_set.tsv \
  --fasta  data/plasmidscope_primary/provenance/working_set.fna.gz \
  --out data/test_plasmids.fna --n 100
```

The selection is seeded, so the same command reproduces the same 100 plasmids. The file
is ignored by git with the rest of `/data/`; the generator is committed, which is what
makes the set reproducible rather than merely archived.

## Running it

```bash
envs/plasmidann/bin/snakemake -s workflow/Snakefile --configfile config/test/config.yaml -c 8
```

Outputs land in `results_test/`, so a smoke run can never overwrite a production run.

## What the sample is chosen to exercise

Stratified by topology and size band, round-robin rather than proportionally. A
proportional sample of this collection would be about 98 circular and 2 linear, and the
linear branch would barely be touched.

| topology | what it tests |
|---|---|
| `circular` | origin repair RUNS: genes broken by linearising the molecule are reconstructed |
| `direct terminal repeat` | also treated as closed, so repair runs |
| `inverted terminal repeat` | must be treated as LINEAR - an ITR marks a genuinely linear replicon with hairpin telomeres, and joining its ends would fabricate a gene across a junction that does not exist in the cell |
| `linear` | repair must NOT run |

Measured on the current set, which is what correct behaviour looks like:

| topology | ORFs | origin-spanning | partial |
|---|---|---|---|
| circular | 2,837 | 8 | 0 |
| direct terminal repeat | 665 | 18 | 0 |
| inverted terminal repeat | 724 | 0 | 12 |
| linear | 1,447 | 0 | 20 |

The two closed topologies produce origin-spanning genes and no partials, because the
fragments the cut created were reconstructed. The two linear topologies produce no
origin-spanning genes and keep their real edge partials.

Size bands matter for the same reason: on a small cryptic plasmid a plus or minus three
neighbourhood is the entire molecule, which is the statistical trap the context rules exist
to avoid, and a sample of only large plasmids would never reach it.

## How this configuration differs from production

Only in what the set size forces. Everything else is identical, deliberately: a smoke run
with different thresholds is not testing the pipeline that produces results.

| setting | test | production | why it must differ |
|---|---|---|---|
| `hmmer_z` | 5297 | 3401393 | hmmsearch reports `E = Z x P(score \| null)`, and Z counts every unique protein, including the proteins PlasmidScope annotates and the cascade does not search. A `-Z` pinned to the production size would rescale every E-value in the run. Rule `check_hmmer_z` stops the run after dereplication if the count differs by more than 2%. |
| tiers | T1..T4 | T1..T5 | T5 searches ClusteredNR (208 GB), and DIAMOND streams the whole database whatever the query size, so including it would make a smoke run cost a production run. T1 and T2 exercise hmmer, T3 pharokka and T4 DIAMOND, so every search method stays covered. |
| `outdir` | `results_test` | `results` | a smoke run must not overwrite a production run. |

## What a passing smoke run does and does not tell you

It tells you every stage executes, every file contract holds, the identifiers are stable
and the topology branches behave. It tells you nothing about the dark fraction, the
background rates, or whether any threshold is well chosen - those need the full collection.
