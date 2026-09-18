#!/bin/bash -l
#SBATCH -A uppmax2025-2-42
#SBATCH -p pelle
#SBATCH -n 1
#SBATCH -c 96
#SBATCH --mem=600G
#SBATCH -t 7-00:00:00
#SBATCH -J dark_plasmidome
#SBATCH -o results/logs/pipeline_%j.out
#
# The whole pipeline: 143,503 plasmids -> a complete, evidence-rich annotation of every
# plasmid ORF and its dark subset. It produces no ranking and no candidate list; spec
# section 76 puts experimental prioritisation in a separate downstream workflow.
#
# S8d (structural screening) is GPU-bound and has its own submission script,
# workflow/structure_gpu.sbatch. Running it here works but spends the ProstT5 cost on
# CPU, which dominates the run.
#
# ENVIRONMENT
#
# envs/plasmidann is a materialised copy of workflow/envs/plasmidann.yaml and holds all
# eleven executables the workflow invokes. Putting it on PATH is what makes the conda:
# directives meaningful: an earlier version of this script activated nothing, so all 26
# directives were inert and the run silently depended on whatever the login shell had.
#
# To let Snakemake build the environment itself instead, drop the PATH line and add
# --use-conda. Do not do both: two environments means the pre-flight rule checks one and
# the searches run in the other, which is how a missing tool reaches hour 45 unnoticed.
#
# RESUMING
#
# Every stage is sharded, so a failure costs one shard rather than the stage. Re-submitting
# this script resumes: --rerun-incomplete discards partial outputs from the killed job and
# Snakemake rebuilds only what is missing.
set -euo pipefail

# MEMORY. Every rule declares mem_mb, and in local execution Snakemake ignores all of it
# unless a ceiling is supplied here: it schedules on the core budget alone. Without this the
# 600 GB allocation is respected by luck rather than by accounting - orf_call, integrons and
# 64 concurrent tier_search jobs can all become schedulable at once. 550 GB of the 600 leaves
# headroom for the OS and for tool temporary files.
MEM_MB=550000

ROOT=/gorilla/proj/h-mel-phylo/h-mel-phylo/private/amar/plasmid_analyis
cd "$ROOT"
export PATH="$ROOT/envs/plasmidann/bin:$PATH"

# Pre-flight runs first in its own invocation so a missing tool or database fails in
# seconds, before 96 cores are spent on gene calling.
.venv/bin/snakemake -s workflow/Snakefile -j 1 --nolock \
  --resources mem_mb=$MEM_MB results/05_annotation_cascade/preflight.tsv

.venv/bin/snakemake -s workflow/Snakefile \
  -j 96 \
  --resources mem_mb=$MEM_MB \
  --rerun-incomplete \
  --nolock \
  --keep-going \
  --printshellcmds \
  --rerun-triggers mtime
