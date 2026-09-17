#!/bin/bash -l
#SBATCH -A uppmax2025-2-42
#SBATCH -p pelle
#SBATCH -n 1
#SBATCH -c 16
#SBATCH --mem=200G
#SBATCH -t 12:00:00
#SBATCH -J bench_nr
#SBATCH -o results/logs/bench_nr.%j.out
#
# Benchmark T4 (DIAMOND against nr, 375 GB) on a compute node.
#
# This is deliberately NOT run on a login node: DIAMOND streams the whole database, and
# 375 GB of I/O from a login node is antisocial on a shared filesystem. It is also the one
# tier whose cost cannot be extrapolated from a small database, because the fixed cost of
# reading nr dominates until the query set is large.
#
# THREE MEASUREMENTS, IN THIS ORDER, AND WHY
#
#   1. n=2,000  cold  -> the fixed cost of one DIAMOND invocation against nr
#   2. n=20,000       -> slope = (t2 - t1) / 18,000 residues-weighted queries
#   3. n=2,000  warm  -> the same query set again, with 375 GB now in page cache
#
# The third is the one that decides whether sharding T4 64 ways is survivable. The cascade
# runs tier_search as 64 independent (tier, shard) jobs at 4 threads each, so on a 96-core
# node roughly 24 DIAMOND processes stream nr concurrently. If the fixed cost is dominated
# by reading the database from disk, that sharding multiplies it; if the node's 772 GB of
# RAM holds nr in page cache after the first pass, it does not. Cold-versus-warm is the
# measurement that tells them apart.
#
# maxRSS from measurement 2 is what sets resources.mem_mb for the T4 tier in
# workflow/rules/s3_cascade.smk, where 64000 is currently a declared over-allocation.
#
# Query set: 2,000 and 20,000 unique PlasAnn-dark plasmid proteins sampled (seed 20260910)
# from 8,000 non-simulated plasmids; median length 205 aa, mean 263 aa. These are NOT
# Pfam-filtered, so their hit rate is an upper bound on the true T4 input's - T4 sees only
# what survived T1-T3. Timing is therefore an upper bound on the per-query slope.
set -u
ROOT=/gorilla/proj/h-mel-phylo/h-mel-phylo/private/amar/plasmid_analyis
export PATH=$ROOT/envs/plasmidann/bin:$PATH
# The query sets must live on the SHARED filesystem. A session scratchpad under
# /scratch is node-local: the first attempt at this job (6823965) failed in 11 s
# with "No such file or directory" on every query file.
SP=$ROOT/data/bench
WORK=${TMPDIR:-/tmp}/bench_nr.$SLURM_JOB_ID
mkdir -p "$WORK"; cd "$WORK"

DB=/sw/data/diamond_databases/Blast/latest/nr.dmnd
run () {   # run <tag> <query fasta>
  echo "=== nr, $1, 16 threads ==="
  /usr/bin/time -f "$1 wall=%e s maxRSS=%M KB" \
    diamond blastp -q "$2" -d "$DB" -o "$1.m8" \
      --very-sensitive --evalue 1e-5 --threads 16 --max-target-seqs 5 \
      --outfmt 6 qseqid stitle qcovhsp scovhsp evalue qstart qend qlen --quiet
  echo "queries with a hit: $(cut -f1 "$1.m8" | sort -u | wc -l) / $(grep -c '>' "$2")"
}

run n2000_cold  "$SP/bench2k.faa"
run n20000      "$SP/bench20k.faa"
run n2000_warm  "$SP/bench2k.faa"
echo "=== DONE ==="
