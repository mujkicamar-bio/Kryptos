#!/bin/bash
# dbAPIS (anti-prokaryotic immune system) HMM search over the dark-plasmidome query set.
# Only 290 models vs Pfam's 30,134, so this runs as a single job rather than an array.
#
# NOTE: dbAPIS HMMs carry NO curated GA/TC/NC thresholds, unlike Pfam. --cut_ga is therefore
# impossible; we report a permissive E-value here and apply stringency tiers downstream, so
# the threshold choice is visible and testable rather than baked in.
#SBATCH -A uppmax2025-2-42
#SBATCH -p pelle
#SBATCH -n 1
#SBATCH -c 16
#SBATCH -t 02:00:00
#SBATCH --mem=16G
#SBATCH -J dbapis
#SBATCH -o data/dbapis_run/logs/dbapis_%j.out
set -uo pipefail
cd /gorilla/proj/h-mel-phylo/h-mel-phylo/private/amar/plasmid_analyis
export PATH=/gorilla/home/amujkic/.conda/envs/panaroo/bin:$PATH
hmmsearch -E 1e-3 --domE 1e-3 --noali --cpu 16 \
    --domtblout data/dbapis_run/dbapis.domtbl.partial \
    data/refs/dbapis/dbAPIS.hmm data/dbapis_run/query_all.faa > /dev/null || exit 1
mv data/dbapis_run/dbapis.domtbl.partial data/dbapis_run/dbapis.domtbl
echo "done: $(grep -vc '^#' data/dbapis_run/dbapis.domtbl) domain rows"
