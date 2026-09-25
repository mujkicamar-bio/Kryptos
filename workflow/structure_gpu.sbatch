#!/bin/bash -l
#SBATCH -A uppmax2025-2-42
#SBATCH -p gpu
#SBATCH --gres=gpu:l40s:1
#SBATCH -n 1
#SBATCH -c 16
#SBATCH --mem=120G
#SBATCH -t 12:00:00
#SBATCH -J dark_plasmidome
#SBATCH --dependency=singleton
#SBATCH -o results/logs/structure_gpu_%j.out
#
# S8d structural screening on a GPU.
#
# WHY THIS STAGE GETS ITS OWN JOB
#
# ProstT5 is a transformer that predicts the 3Di structural alphabet for every dark
# protein, and it is the cost of this stage - the Foldseek search that follows it is not.
# On CPU it dominates the whole pipeline; on an L40S it is a small fraction of it. Nothing
# else in the pipeline uses a GPU, so allocating one for the whole run would idle it for
# the entire cascade.
#
# Run the rest of the pipeline on the CPU partition and this stage here:
#
#   sbatch workflow/run_pipeline.sbatch          # everything up to S8d (--omit-from)
#   sbatch workflow/structure_gpu.sbatch         # S8d
#   sbatch workflow/run_pipeline.sbatch          # the rest, resuming
#
# Snakemake resumes from what exists, so the split costs nothing but the queue wait. All
# three share the job name dark_plasmidome and --dependency=singleton, so they can be
# submitted together and SLURM runs them one after another.
#
# --rerun-triggers mtime, as in run_pipeline.sbatch: targets.gpu.yaml changes the params
# of rules upstream of structure_search (preflight among them), and with the default
# triggers that would re-run preflight and everything after it - the cascade included.
#
# WHY --gpu IS SET IN CONFIG RATHER THAN DETECTED
#
# A stage that silently behaves differently depending on which node it landed on cannot
# have its results traced to the settings that produced them, and this pipeline records the
# thresholds behind every row precisely so that it can. A run that asked for a GPU and did
# not get one should fail visibly rather than take ten times longer without saying so.
set -euo pipefail

ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
cd "$ROOT"
export PATH="$ROOT/envs/plasmidann/bin:$PATH"
export XDG_CACHE_HOME="${XDG_CACHE_HOME:-$ROOT/.cache}"

nvidia-smi --query-gpu=name,memory.total --format=csv || {
    echo "no GPU visible in this allocation - refusing to fall back to CPU silently" >&2
    exit 1
}

snakemake \
    --configfile "${CONFIGFILE:-config/config.yaml}" \
    --config targets_config=config/targets.gpu.yaml \
    -j 16 --rerun-incomplete \
    --rerun-triggers mtime \
    --until structure_search
