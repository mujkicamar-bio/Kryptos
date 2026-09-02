#!/bin/bash
# Task T7 of plans/conservative_reclustering.md — cluster the WHOLE dark plasmid proteome, not
# just the cryptic small subset: 378,552 small-plasmid dark ORFs + 5,135,708 from plasmids
# >= 20 kb (extracted by scripts/fam_bigplasmid_search.sbatch) = 5,514,260, zero ID collisions.
#
# Two thresholds: 30% for comparability with everything published so far, 50% (UniRef50) as the
# conservative standard. Both with --cluster-reassign and -c 0.8 --cov-mode 0.
# Submit: sbatch scripts/recluster_union.sbatch
#SBATCH -A uppmax2025-2-42
#SBATCH -p pelle
#SBATCH -n 1
#SBATCH -c 24
#SBATCH -t 24:00:00
#SBATCH --mem=192G
#SBATCH -J recl_union
#SBATCH -o data/dark_orf_run/union/union_%j.out
set -euo pipefail
cd /gorilla/proj/h-mel-phylo/h-mel-phylo/private/amar/plasmid_analyis
export PATH=/gorilla/home/amujkic/.conda/envs/panaroo/bin:$PATH
IN=data/dark_orf_run/union/all_dark_orfs.faa
for ID in 0.5 0.3; do
    TAG="union$(printf '%.0f' "$(echo "$ID * 100" | bc -l)")"
    mkdir -p "data/dark_orf_run/union/tmp_$TAG"
    echo "[$(date +%H:%M:%S)] $TAG  min-seq-id=$ID -c 0.8 --cov-mode 0 --cluster-reassign"
    mmseqs easy-cluster "$IN" "data/dark_orf_run/union/$TAG" \
        "data/dark_orf_run/union/tmp_$TAG" \
        --min-seq-id "$ID" -c 0.8 --cov-mode 0 --cluster-reassign 1 --threads 24 -v 1
    rm -rf "data/dark_orf_run/union/tmp_$TAG" "data/dark_orf_run/union/${TAG}_all_seqs.fasta"
    echo "[$(date +%H:%M:%S)] $TAG done: $(cut -f1 "data/dark_orf_run/union/${TAG}_cluster.tsv" | sort -u | wc -l) families"
done
echo "UNION DONE"
