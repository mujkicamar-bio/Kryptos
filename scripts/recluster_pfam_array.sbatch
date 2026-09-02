#!/bin/bash
# Pfam-A 38.2 over the cluster representatives introduced by the threshold sweep.
#
# The original Pfam run (scripts/pfam_array.sbatch) searched the representatives of the 30%
# clustering. A different threshold produces different representatives, so 75,544 of the reps
# across dark30/50/70/90 had never been searched. Same tool, same --cut_ga criterion, so the
# results are directly poolable with data/pfam_run/pfam_per_seq.tsv.
#
# Build chunks first (done inline; 38 chunks of 2,000), then: sbatch scripts/recluster_pfam_array.sbatch
#SBATCH -A uppmax2025-2-42
#SBATCH -p pelle
#SBATCH -n 1
#SBATCH -c 4
#SBATCH -t 04:00:00
#SBATCH --mem=12G
#SBATCH -J pfam_recl
#SBATCH --array=0-37%38
#SBATCH -o data/pfam_run/recluster/logs/pfam_%A_%a.out
set -uo pipefail
cd /gorilla/proj/h-mel-phylo/h-mel-phylo/private/amar/plasmid_analyis
export PATH=/gorilla/home/amujkic/.conda/envs/panaroo/bin:$PATH
i=$(printf "%04d" "$SLURM_ARRAY_TASK_ID")
Q="data/pfam_run/recluster/chunks/q$i.faa"
O="data/pfam_run/recluster/domtbl/q$i.domtbl"
[ -s "$Q" ] || { echo "missing $Q"; exit 1; }
hmmsearch --cut_ga --noali --cpu 4 --domtblout "$O.partial" \
    data/refs/pfam/Pfam-A.hmm "$Q" > /dev/null || exit 1
mv "$O.partial" "$O"
echo "chunk $i done: $(grep -vc '^#' "$O") domain rows"
