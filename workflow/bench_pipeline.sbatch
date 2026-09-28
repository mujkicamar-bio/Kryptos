#!/bin/bash -l
#SBATCH -A uppmax2025-2-42
#SBATCH -p pelle
#SBATCH -n 1
#SBATCH -c 96
#SBATCH --mem=600G
#SBATCH -t 24:00:00
#SBATCH -J bench_pipeline
#SBATCH -o slurm-%x-%j.out
#
# Cost measurement: the whole pipeline on the 100-plasmid set, every production tier,
# with every rule benchmarked.
#
# WHAT THIS IS FOR
#
# Sizing the production run. Every rule carries a benchmark: directive, so each one
# writes results_bench/benchmarks/<rule>.tsv with wall seconds, max RSS, max VMS, I/O and
# mean load: the basis for each rule's mem_mb and runtime.
#
# It is NOT the smoke run. config/test/ checks that the pipeline is correct; this measures
# what it costs, and the two need different tier sets. config/bench/cascade.yaml therefore
# carries five tiers including the production T5 on ClusteredNR, which
# config/test/cascade.yaml omits. Output goes to results_bench, so
# neither run can overwrite the other. The Slurm log is written to the submit directory.
#
# WHAT THE NUMBERS SUPPORT, AND WHAT THEY DO NOT
#
# Memory transfers. Peak RSS is set by the database and the block size, not by the query
# count, for every search stage here - measured on nr, where 25x the queries moved peak
# RSS by 36% (129.6 -> 176.7 GB, job 6928733). A per-rule mem_mb read off these files is
# therefore a real number for production, with headroom added.
#
# Wall time does NOT transfer by multiplication. Three reasons, in order of size:
#
#   1. The streamed-database tiers are mostly a per-pass cost. nr measured at ~41 min per
#      pass plus 1.71 s per query on a full node. At 100 plasmids the per-pass cost is
#      essentially the whole figure; at 143,503 it is a rounding error. Extrapolating the
#      test figure linearly overestimates production by orders of magnitude, and
#      extrapolating the slope alone underestimates it.
#   2. Several stages are superlinear in the protein count by construction - clustering,
#      orthology and synteny all compare proteins to each other, not to a database. A
#      per-protein cost measured at 5,304 proteins does not hold at 3.5 million.
#   3. The dark set drives S6 to S8, and its SIZE is a result of the cascade, not an input
#      to it. Adding nr changes how many proteins survive to be clustered.
#
# So: read mem_mb off these files directly; treat wall time as a lower bound and a ranking
# of which stages deserve attention, not as a figure to multiply by 1,435.
#
# HOW TO RUN IT
#
#   sbatch workflow/bench_pipeline.sbatch                                  # 1. CPU, S0-S8c
#   sbatch --export=CONFIGFILE=config/bench/config.yaml \
#          workflow/structure_gpu.sbatch                                   # 2. GPU, S8d
#   sbatch --export=PHASE=finish workflow/bench_pipeline.sbatch            # 3. CPU, report
#
# Step 2 is on the GPU partition because ProstT5 on CPU dominates every other stage and
# would make the measurement describe a run nobody intends to perform. Step 3 is the
# single rule that consumes the structural hits; it is minutes, not hours.
set -euo pipefail

# Slurm runs a copy of this script from its spool directory, so the repository root is the
# directory the job was submitted from: submit from the repository root.
ROOT=${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}
[ -f "$ROOT/workflow/Snakefile" ] || { echo "submit from the repository root (no $ROOT/workflow/Snakefile)" >&2; exit 1; }
cd "$ROOT"
# Snakemake runs a `script:` directive with the interpreter running Snakemake, not with
# the first python on PATH, so the environment holding the tools must also run Snakemake.
SNAKEMAKE="$ROOT/envs/plasmidann/bin/snakemake"
export PATH="$ROOT/envs/plasmidann/bin:$PATH"

# See run_pipeline.sbatch: in local execution Snakemake ignores every declared mem_mb
# unless a ceiling is supplied here, and schedules on the core budget alone. Without it
# orf_call, integrons and the nr tier can become schedulable at once inside one 600 GB
# allocation. The nr tier alone peaked at 177 GB.
MEM_MB=550000
CONFIG=config/bench/config.yaml

case "${PHASE:-cpu}" in
  cpu)
    # Pre-flight in its own invocation, so a missing tool or database fails in seconds
    # rather than after the cascade has spent a node for hours.
    # The target follows --nolock, a boolean flag: after --resources or --configfile,
    # which take one or more values, it would be read as one of those values.
    "$SNAKEMAKE" -s workflow/Snakefile --configfile "$CONFIG" \
      --resources mem_mb=$MEM_MB -j 1 --nolock \
      results_bench/05_annotation_cascade/preflight.tsv

    # --forceall because a benchmark of a stage Snakemake skipped is an empty file. Every
    # rule must actually execute for its cost to be measured.
    #
    # --omit-from structure_search holds back S8d and the one rule downstream of it
    # (annotation_report), which is what leaves them for the GPU job and for PHASE=finish.
    #
    # No --keep-going: unlike a production run, a failed stage here invalidates the
    # measurement of everything that would have followed it, so stopping is more useful.
    "$SNAKEMAKE" -s workflow/Snakefile --configfile "$CONFIG" \
      -j 96 \
      --resources mem_mb=$MEM_MB \
      --forceall \
      --omit-from structure_search \
      --nolock \
      --printshellcmds \
      --rerun-triggers mtime
    ;;
  finish)
    # After the GPU job. No --forceall: everything except annotation_report exists, and
    # rebuilding it would discard the measurement this script just spent a node on.
    "$SNAKEMAKE" -s workflow/Snakefile --configfile "$CONFIG" \
      -j 96 \
      --resources mem_mb=$MEM_MB \
      --nolock \
      --printshellcmds \
      --rerun-triggers mtime
    ;;
  *)
    echo "PHASE must be 'cpu' or 'finish', not '${PHASE}'" >&2
    exit 1
    ;;
esac

echo "=== benchmarks written ==="
head -n 1 -q results_bench/benchmarks/*.tsv 2>/dev/null | head -1
for f in results_bench/benchmarks/*.tsv; do
    printf '%-28s %s\n' "$(basename "$f" .tsv)" "$(sed -n '2p' "$f")"
done
echo "=== DONE $(date) ==="
