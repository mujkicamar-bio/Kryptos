#!/bin/bash
# PlasAnn annotation of the full working set: one array task per shard, node-local $SNIC_TMP for the
# transient per-plasmid files, compact outputs (annot/*.tsv.gz + gbk/*.tar.gz) back to project FS.
# Submit:  sbatch scripts/plasann_array.sbatch
# Resume:  sbatch --array=<comma-list>%300 scripts/plasann_array.sbatch   (re-run failed/missing shards)
#SBATCH -A uppmax2025-2-42
#SBATCH -p pelle
#SBATCH -n 1
#SBATCH -t 12:00:00
#SBATCH --mem=8G
#SBATCH -J plasann
#SBATCH --array=0-594%300
#SBATCH -o data/plasann_run/logs/plasann_%A_%a.out
set -uo pipefail
cd /gorilla/proj/h-mel-phylo/h-mel-phylo/private/amar/plasmid_analyis
bash scripts/plasann_task.sh "$SLURM_ARRAY_TASK_ID"
