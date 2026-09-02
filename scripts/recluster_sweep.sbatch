#!/bin/bash
# Threshold sweep for the dark plasmid proteome. Precedent for sweeping rather than picking one
# number: Yu et al., Nat Microbiol 2024 (PlasX/MobMess) ran 12 thresholds on plasmid gene
# families. The four here are chosen to be citable landmarks:
#   0.3  the existing run, rebuilt WITH --cluster-reassign so it is comparable
#   0.5  UniRef50 (50% identity, 80% overlap)
#   0.7  NMPFamsDB (NAR 2024)
#   0.9  UniRef90
# Submit: sbatch scripts/recluster_sweep.sbatch
#SBATCH -A uppmax2025-2-42
#SBATCH -p pelle
#SBATCH -n 1
#SBATCH -c 24
#SBATCH -t 08:00:00
#SBATCH --mem=64G
#SBATCH -J recluster
#SBATCH -o data/dark_orf_run/recluster/logs/sweep_%j.out
set -euo pipefail
cd /gorilla/proj/h-mel-phylo/h-mel-phylo/private/amar/plasmid_analyis
for ID in 0.3 0.5 0.7 0.9; do
    THREADS=24 bash scripts/recluster_dark_orfs.sh "$ID"
done
echo "SWEEP DONE"
