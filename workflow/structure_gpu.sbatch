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
# Rule structure_search on a GPU. ProstT5 inference is the cost of this stage; nothing
# else in the pipeline uses a GPU, so the stage gets its own job rather than a GPU idle for
# the whole run.
#
# Optional: by default workflow/run_pipeline.sbatch runs structure_search on CPU. To run it here:
#
#   STRUCTURE_ON_GPU=1 sbatch workflow/run_pipeline.sbatch   # up to structure_search (--omit-from)
#   sbatch workflow/structure_gpu.sbatch                     # structure_search
#   sbatch workflow/run_pipeline.sbatch                      # the rest, resuming
#
# All three share the job name dark_plasmidome and --dependency=singleton, so SLURM runs
# them one after another. --rerun-triggers mtime: targets.gpu.yaml changes the params of
# rules upstream of structure_search, which would otherwise re-run them. The GPU is set in
# config rather than detected, and a job without one fails rather than running on CPU.
set -euo pipefail

# Slurm runs a copy of this script from its spool directory, so the repository root is the
# directory the job was submitted from: submit from the repository root.
ROOT=${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}
[ -f "$ROOT/workflow/Snakefile" ] || { echo "submit from the repository root (no $ROOT/workflow/Snakefile)" >&2; exit 1; }
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
