#!/bin/bash -l
#SBATCH -A uppmax2025-2-422
#SBATCH -p pelle
#SBATCH -n 1
#SBATCH -c 96
#SBATCH --mem=600G
#SBATCH -t 10-00:00:00
#SBATCH -J dark_plasmidome
#SBATCH --dependency=singleton
#SBATCH -o slurm-%x-%j.out
#
# The whole pipeline: 143,503 plasmids -> a complete, evidence-rich annotation of every
# plasmid ORF and its dark subset. It produces no ranking and no candidate list. The Slurm
# log is written to the submit directory (the repository root).
#
# ONE SUBMISSION ANNOTATES EVERYTHING. `sbatch workflow/run_pipeline.sbatch` runs every
# stage in this one job: the cascade, the plasmid label databases and AMRFinderPlus
# (label_databases), DefenseFinder, CONJScan (conjugation_systems), IntegronFinder,
# ISEScan, synteny, the context terms, the structure search and the report.
#
# ENVIRONMENT
#
# envs/plasmidann, built from workflow/envs/plasmidann.yaml, runs Snakemake and holds the
# tools; it is put on PATH, so the conda: directives match the running environment. Do not
# add --use-conda as well: pre-flight would then check one environment while the searches
# run in another.
#
# RESUMING
#
# Every stage is one job. Re-submitting this script resumes: --rerun-incomplete discards
# the partial output of the job that was killed and Snakemake rebuilds only what is
# missing, so a failure costs one stage rather than the run. DIAMOND cannot checkpoint, so
# a tier killed by the time limit restarts from zero: -t is the partition maximum (10 days)
# because the serial estimate is ~120-230 h.
#
# --dependency=singleton, with the same job name in structure_gpu.sbatch: both invocations
# below run --nolock, so a second submission would otherwise write the same outputs at the
# same time. SLURM starts it only when the first has ended.
#
# STRUCTURE SEARCH (S8d) runs HERE, on CPU, by default (~15-35 h on 16 CPU threads for the
# expected dark representatives), so that one submission produces the whole annotation.
# ProstT5 is far faster on a GPU, so a GPU split remains available as an option:
#
#   STRUCTURE_ON_GPU=1 sbatch workflow/run_pipeline.sbatch   # stops before structure_search
#   sbatch workflow/structure_gpu.sbatch                     # S8d on the GPU partition
#   sbatch workflow/run_pipeline.sbatch                      # the rest, resuming
#
# With STRUCTURE_ON_GPU=1 the job stops before structure_search (--omit-from) until
# structure_hits.tsv exists. The split trades queue waits for CPU hours; the default
# trades CPU hours for a single run, which is what was asked for.
#
# TOOLS IN THEIR OWN ENVIRONMENTS. AMRFinderPlus (envs/amrfinder) and the MacSyFinder 2.1.6
# that runs CONJScan (envs/conjscan) are not on PATH: config names them by path
# (amrfinder.executable, conjugation.exe), as it names pharokka, and pre-flight checks them
# there - including that the MacSyFinder is new enough for CONJScan's model grammar.
set -euo pipefail

# MEMORY. Every rule declares mem_mb, and in local execution Snakemake ignores all of it
# unless a ceiling is supplied here: it schedules on the core budget alone. Without this the
# 600 GB allocation is respected by luck rather than by accounting - orf_call, integrons and
# the nr tier can all become schedulable at once. 550 GB of the 600 leaves
# headroom for the OS and for tool temporary files.
MEM_MB=550000
CORES=${SLURM_CPUS_PER_TASK:-96}
OMIT=()
if [ -n "${STRUCTURE_ON_GPU:-}" ] && [ ! -e results/12_context_and_structure/structure_hits.tsv ]; then
    OMIT=(--omit-from structure_search)
fi

# Slurm runs a copy of this script from its spool directory, so the repository root is the
# directory the job was submitted from: submit from the repository root.
ROOT=${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}
[ -f "$ROOT/workflow/Snakefile" ] || { echo "submit from the repository root (no $ROOT/workflow/Snakefile)" >&2; exit 1; }
cd "$ROOT"

# Snakemake runs a `script:` directive with the interpreter running Snakemake, not with the
# first python on PATH, so the environment that holds the tools also runs Snakemake.
SNAKEMAKE="$ROOT/envs/plasmidann/bin/snakemake"
export PATH="$ROOT/envs/plasmidann/bin:$PATH"
# Tool caches off the home directory, whose file quota is 94% used.
export XDG_CACHE_HOME="${XDG_CACHE_HOME:-$ROOT/.cache}"

if [ ! -x "$SNAKEMAKE" ]; then
    echo "no snakemake at $SNAKEMAKE - see README.md, section Installation" >&2
    exit 1
fi

# Pre-flight runs first in its own invocation so a missing tool or database fails in
# seconds, before the whole node is spent on gene calling.
#
# The target follows --nolock, a boolean flag: after --resources or --configfile, which
# take one or more values, it would be read as one of those values.
"$SNAKEMAKE" -s workflow/Snakefile \
  --resources mem_mb=$MEM_MB -j 1 --nolock \
  results/05_annotation_cascade/preflight.tsv

"$SNAKEMAKE" -s workflow/Snakefile \
  -j "$CORES" \
  --resources mem_mb=$MEM_MB \
  --rerun-incomplete \
  --nolock \
  --keep-going \
  --printshellcmds \
  --rerun-triggers mtime \
  "${OMIT[@]}"
