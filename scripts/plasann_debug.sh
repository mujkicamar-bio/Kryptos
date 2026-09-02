#!/bin/bash
set -uo pipefail
ROOT=/gorilla/proj/h-mel-phylo/h-mel-phylo/private/amar/plasmid_analyis
export PATH=/gorilla/home/amujkic/.conda/envs/plasann_env/bin:$PATH
echo "host=$(hostname) SNIC_TMP=${SNIC_TMP:-unset}"
echo "DB: $(ls -d ~/.plasann/Database 2>&1 | head -1)"
T="${SNIC_TMP:-/tmp}/dbg"; rm -rf "$T"; mkdir -p "$T/in" "$T/out"
zcat "$ROOT/data/plasann_run/shards/shard_0000.fna.gz" | awk -v D="$T/in" '
  /^>/ { n++; if(n>2) exit; id=substr($1,2); gsub(/[^A-Za-z0-9._-]/,"_",id); f=D"/"id".fasta"; print > f; next }
  n<=2 { print >> f }'
echo "inputs=$(ls "$T/in" | wc -l)"
cd "$T"                       # run from node-local scratch (hypothesis: PlasAnn writes temp in CWD)
PlasAnn -i "$T/in" -o "$T/out" -t fasta 2>&1 | tail -25
echo "RESULT csvs=$(find "$T/out" -name '*_annotations.csv' | wc -l)"
