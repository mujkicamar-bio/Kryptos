#!/bin/bash
# PlasmidFinder replicon typing of the whole working set: blastn vs the PlasmidFinder DB at the standard
# >=80% identity threshold (coverage filtered to >=60% of the replicon reference at parse time), matching
# PLSDB's plsdb_inc_types method. One multithreaded job (small DB -> fast).
# Submit:  sbatch scripts/plasmidfinder_run.sbatch
#SBATCH -A uppmax2025-2-42
#SBATCH -p pelle
#SBATCH -n 6
#SBATCH -t 12:00:00
#SBATCH --mem=24G
#SBATCH -J plasmidfinder
#SBATCH -o data/plasann_run/logs/plasmidfinder_%j.out
set -uo pipefail
cd /gorilla/proj/h-mel-phylo/h-mel-phylo/private/amar/plasmid_analyis
export PATH=/gorilla/home/amujkic/.conda/envs/plasann_env/bin:$PATH

TMP="${SNIC_TMP:-$(mktemp -d)}"
echo "host=$(hostname) start=$(date +%F_%T)"
zcat data/plasmidscope_primary/working_set.fna.gz > "$TMP/ws.fna"
echo "queries: $(grep -c '^>' "$TMP/ws.fna")  decompressed=$(date +%T)"

# raw blast: qseqid sseqid pident length slen (slen = replicon-reference length, for coverage)
blastn -query "$TMP/ws.fna" -db data/typing/plasmidfinder -task blastn \
  -perc_identity 80 -evalue 1e-5 -max_target_seqs 200 -num_threads 6 -mt_mode 1 \
  -outfmt "6 qseqid sseqid pident length slen" 2>"$TMP/err" | gzip > data/typing/plasmidfinder_blast.tsv.gz
echo "blast done=$(date +%F_%T); hits=$(zcat data/typing/plasmidfinder_blast.tsv.gz | wc -l)"
tail -3 "$TMP/err" 2>/dev/null || true
