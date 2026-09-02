#!/bin/bash
# RGI/CARD antibiotic-resistance annotation over the full working set (595 shards, ~350 plasmids each).
# Submit:  sbatch scripts/card_array.sbatch
#SBATCH -A uppmax2025-2-42
#SBATCH -p pelle
#SBATCH -n 1
#SBATCH -c 2
#SBATCH -t 01:00:00
#SBATCH --mem=8G
#SBATCH -J card_rgi
#SBATCH --array=0-594%300
#SBATCH -o data/card_run/logs/card_%A_%a.out
set -uo pipefail
cd /gorilla/proj/h-mel-phylo/h-mel-phylo/private/amar/plasmid_analyis
bash scripts/card_task.sh "$SLURM_ARRAY_TASK_ID"
