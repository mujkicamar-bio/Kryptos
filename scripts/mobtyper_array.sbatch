#!/bin/bash
# MOB-suite mob_typer over the full working set (595 shards).
# Submit:  sbatch scripts/mobtyper_array.sbatch
#SBATCH -A uppmax2025-2-42
#SBATCH -p pelle
#SBATCH -n 1
#SBATCH -t 02:00:00
#SBATCH --mem=8G
#SBATCH -J mobtyper
#SBATCH --array=0-594%300
#SBATCH -o data/plasann_run/logs/mobtyper_%A_%a.out
set -uo pipefail
cd /gorilla/proj/h-mel-phylo/h-mel-phylo/private/amar/plasmid_analyis
bash scripts/mobtyper_task.sh "$SLURM_ARRAY_TASK_ID"
