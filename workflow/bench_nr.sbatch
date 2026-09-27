#!/bin/bash -l
#SBATCH -A uppmax2025-2-42
#SBATCH -p pelle
#SBATCH -n 1
#SBATCH -c 96
#SBATCH --mem=600G
#SBATCH -t 48:00:00
#SBATCH -J bench_nr
#SBATCH -o results/logs/bench_nr.%j.out
#
# Benchmark of DIAMOND against full nr (375 GB: 707 M sequences, 273 G residues), kept as
# the provenance of the nr tier's resources in workflow/rules/annotation_cascade.smk
# (memory 220 GB; about 45 min per database pass and 1.7 s per query on 96 threads). It
# measures full nr at --max-target-seqs 5, not the ClusteredNR database and -k 25 of the
# production tier, so rerunning it does not check the production settings.
#
# It runs on a compute node: DIAMOND streams the whole database, and the constant cost of
# one pass over nr dominates until the query set is large, so the cost cannot be
# extrapolated from a small database. A first attempt (job 6823981: 2,000 queries, 16
# threads, DIAMOND's default -b 2 -c 4) did not finish inside 12 hours: with those defaults
# nr is ~137 blocks, each indexed four times per seed shape, whatever the query count.
# This is why each cascade tier runs as one job, and why this run uses the whole node,
# -c 1 (one index pass per block) and -b 16 (~17 blocks).
#
# Three measurements, in this order:
#
#   1. test5   cold  -> 152 proteins from 5 plasmids of the 100-plasmid smoke set: the
#                       constant cost of one pass
#   2. test5   warm  -> the same queries again, nr now in page cache (772 GB node): if this
#                       is not much faster than 1, the pass is CPU-bound indexing, not disk
#   3. test100       -> the 3,808 smoke-set proteins T4 left unresolved: slope = (t3 - t1) / 3,656
#
# maxRSS of measurement 3 is what set resources.mem_mb for the nr tier.
#
# Query sets: exactly what the production nr tier receives - the proteins the earlier
# tiers could not explain. results_test/05_annotation_cascade/T4/unresolved.faa of the
# smoke run, minus the 100 decoys and 2 controls (data/bench/test100_unresolved.faa), and
# its subset from five plasmids sampled with seed 20260922 (data/bench/test5_unresolved.faa,
# plasmid ids in test5_plasmids.ids).
set -u
ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
export PATH=$ROOT/envs/plasmidann/bin:$PATH
# The query sets must live on the shared filesystem: /scratch is node-local, so a query
# file there is not found on the compute node.
SP=$ROOT/data/bench
WORK=${TMPDIR:-/tmp}/bench_nr.$SLURM_JOB_ID
mkdir -p "$WORK"; cd "$WORK"

# The alignments are written to the shared filesystem, because $WORK is node-local scratch
# purged when the job ends, and a timing without the hits cannot say what matched. $WORK
# is DIAMOND's spill space.
OUT=$ROOT/results/bench/nr
mkdir -p "$OUT"

DB=/sw/data/diamond_databases/Blast/latest/nr.dmnd
THREADS=$SLURM_CPUS_PER_TASK
run () {   # run <tag> <query fasta>
  echo "=== nr, $1, $THREADS threads, -b 16 -c 1, $(date) ==="
  /usr/bin/time -f "$1 wall=%e s maxRSS=%M KB" \
    diamond blastp -q "$2" -d "$DB" -o "$OUT/$1.m8" \
      --very-sensitive --evalue 1e-5 --threads $THREADS --max-target-seqs 5 \
      -b 16 -c 1 --tmpdir "$WORK" \
      --outfmt 6 qseqid stitle qcovhsp scovhsp evalue qstart qend qlen --quiet
  echo "queries with a hit: $(cut -f1 "$OUT/$1.m8" | sort -u | wc -l) / $(grep -c '>' "$2")"
}

run test5_cold  "$SP/test5_unresolved.faa"
run test5_warm  "$SP/test5_unresolved.faa"
run test100     "$SP/test100_unresolved.faa"
echo "=== DONE $(date) ==="
