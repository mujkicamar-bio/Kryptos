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
# Benchmark the nr tier (DIAMOND against nr, 375 GB: 707 M sequences, 273 G residues) on
# the node and with the settings the production tier will use.
#
# This is deliberately NOT run on a login node: DIAMOND streams the whole database, and
# 375 GB of I/O from a login node is antisocial on a shared filesystem. It is also the one
# tier whose cost cannot be extrapolated from a small database, because the fixed cost of
# one pass over nr dominates until the query set is large.
#
# WHAT THE FIRST ATTEMPT SHOWED
#
# Job 6823981: 2,000 queries, 16 threads, DIAMOND's default block size and index chunking
# (-b 2 -c 4), did not finish inside 12 hours. That is the fixed cost of one pass: with the
# defaults nr is ~137 blocks, each indexed four times per seed shape, whatever the query
# count. It is why the cascade now runs each tier as ONE job (spec section 72), and why
# this run uses the whole node, -c 1 (one index pass per block) and -b 16 (~17 blocks).
#
# THREE MEASUREMENTS, IN THIS ORDER, AND WHY
#
#   1. test5   cold  -> 152 proteins from 5 plasmids of the 100-plasmid smoke set: the
#                       fixed cost of one pass with the production settings
#   2. test5   warm  -> the same queries again, nr now in page cache (772 GB node): if this
#                       is not much faster than 1, the pass is CPU-bound indexing, not disk
#   3. test100       -> the 3,808 smoke-set proteins T4 left unresolved: slope = (t3 - t1) / 3,656
#
# maxRSS from measurement 3 is what sets resources.mem_mb for the nr tier in
# workflow/rules/annotation_cascade.smk, where 400000 is currently a declared placeholder.
#
# Query sets: exactly what the production nr tier receives - the proteins the earlier
# tiers could not explain. results_test/05_annotation_cascade/T4/unresolved.faa of the
# smoke run, minus the 100 decoys and 2 controls (data/bench/test100_unresolved.faa), and
# its subset from five plasmids sampled with seed 20260922 (data/bench/test5_unresolved.faa,
# plasmid ids in test5_plasmids.ids).
set -u
ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
export PATH=$ROOT/envs/plasmidann/bin:$PATH
# The query sets must live on the SHARED filesystem. A session scratchpad under
# /scratch is node-local: the first attempt at this job (6823965) failed in 11 s
# with "No such file or directory" on every query file.
SP=$ROOT/data/bench
WORK=${TMPDIR:-/tmp}/bench_nr.$SLURM_JOB_ID
mkdir -p "$WORK"; cd "$WORK"

# The ALIGNMENTS must outlive the job, and $WORK does not: it is node-local scratch, which
# is purged when the job ends. Run 6928733 measured its timings correctly and then lost
# every hit it had found, so the run could say that 3,689 of 3,808 queries matched nr but
# not what any of them matched - and a hit to "hypothetical protein" is not an annotation.
# $WORK stays as DIAMOND's spill space, which is exactly what node-local disk is for.
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
