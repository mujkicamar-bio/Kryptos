#!/bin/bash
#SBATCH -A uppmax2025-2-42
#SBATCH -p pelle
#SBATCH -n 1
#SBATCH -c 96
#SBATCH --mem=600G
#SBATCH -t 12:00:00
#SBATCH -J plasmid_s1s2
#SBATCH -o results/logs/s1s2_%j.out
set -euo pipefail
cd /gorilla/proj/h-mel-phylo/h-mel-phylo/private/amar/plasmid_analyis
.venv/bin/snakemake -s workflow/Snakefile -j 96 --rerun-incomplete --nolock
