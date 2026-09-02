#!/bin/bash
# Pfam-A 38.2 search over the dark-plasmidome query set (DARKREP + NAMED control + VALID).
# Build the chunks first:  python scripts/build_pfam_query.py
# Submit:                  sbatch scripts/pfam_array.sbatch
#SBATCH -A uppmax2025-2-42
#SBATCH -p pelle
#SBATCH -n 1
#SBATCH -c 4
#SBATCH -t 02:00:00
#SBATCH --mem=12G
#SBATCH -J pfam_dark
#SBATCH --array=0-87%44
#SBATCH -o data/pfam_run/logs/pfam_%A_%a.out
set -uo pipefail
cd /gorilla/proj/h-mel-phylo/h-mel-phylo/private/amar/plasmid_analyis
bash scripts/pfam_task.sh "$SLURM_ARRAY_TASK_ID"
