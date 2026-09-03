#!/bin/bash
#SBATCH -A uppmax2025-2-42
#SBATCH -p pelle
#SBATCH -n 1
#SBATCH -c 96
#SBATCH --mem=600G
#SBATCH -t 3-00:00:00
#SBATCH -J plasmid_annot
#SBATCH -o results/logs/s1s2_%j.out
set -euo pipefail
source /usr/share/lmod/lmod/init/bash
module load HMMER/3.4-gompi-2024a
module load MMseqs2/18-8cc5c-gompi-2024a
cd /gorilla/proj/h-mel-phylo/h-mel-phylo/private/amar/plasmid_analyis
.venv/bin/snakemake -s workflow/Snakefile -j 96 --rerun-incomplete --nolock
