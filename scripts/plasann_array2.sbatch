#!/bin/bash
# RETRY: annotate only the plasmids that OOM'd in the first pass, using smaller re-shards (shards2/)
# so memory stays under the limit. Array size is set at submit time to the number of mshard_*.fna.gz.
# Submit:  sbatch --array=0-<N-1>%300 scripts/plasann_array2.sbatch
#SBATCH -A uppmax2025-2-42
#SBATCH -p pelle
#SBATCH -n 1
#SBATCH -t 08:00:00
#SBATCH --mem=6G
#SBATCH -J plasann2
#SBATCH -o data/plasann_run/logs/plasann2_%A_%a.out
set -uo pipefail
cd /gorilla/proj/h-mel-phylo/h-mel-phylo/private/amar/plasmid_analyis
export SHARDS_SUBDIR=shards2 SHARD_PREFIX=mshard
bash scripts/plasann_task.sh "$SLURM_ARRAY_TASK_ID"
