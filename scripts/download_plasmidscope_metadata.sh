#!/usr/bin/env bash
# Stage 1 of PLASMIDSCOPE_INTEGRATION_PLAN.md: download per-source PlasmidScope metadata tables
# (plasmid_list.download.tsv) only. Metadata only -- no sequences. Complete-only filtering happens
# downstream in scripts/plasmidscope_composition.py using the Completeness column.
set -euo pipefail

BASE="https://plasmidapi.deepomics.org/api/database/files"
OUT="data/PlasmidScope/metadata"
mkdir -p "$OUT"

DATASETS=(IMG-PR PLSDB COMPASS GenBank RefSeq ENA DDBJ mMGE Kraken2)

for ds in "${DATASETS[@]}"; do
  dest="$OUT/${ds}.plasmid_list.tsv"
  if [[ -s "$dest" ]]; then
    echo "[skip] $ds already present ($(wc -l < "$dest") lines)"
    continue
  fi
  url="$BASE/${ds}/data/${ds}.plasmid_list.download.tsv"
  echo "[get ] $ds <- $url"
  curl -fsSL -m 900 "$url" -o "$dest"
  echo "[ok  ] $ds -> $dest ($(wc -l < "$dest") lines)"
done

echo "All metadata tables downloaded to $OUT/"
