#!/bin/bash
# DefenseFinder 3.0.0 + AntiDefenseFinder over the payload-free small plasmidome.
# Build input first:  python scripts/build_defensefinder_input.py --per-chunk 2000
# Submit:             sbatch scripts/defensefinder_array.sbatch
#SBATCH -A uppmax2025-2-42
#SBATCH -p pelle
#SBATCH -n 1
#SBATCH -c 8
#SBATCH -t 00:30:00
#SBATCH --mem=8G
#SBATCH -J dfinder
#SBATCH --array=0-32%33
#SBATCH -o data/defensefinder_run/logs/df_%A_%a.out
set -uo pipefail
cd /gorilla/proj/h-mel-phylo/h-mel-phylo/private/amar/plasmid_analyis
bash scripts/defensefinder_task.sh "$SLURM_ARRAY_TASK_ID"
