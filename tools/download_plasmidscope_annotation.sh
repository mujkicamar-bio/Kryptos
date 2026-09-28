#!/usr/bin/env bash
# Download PlasmidScope's per-protein annotation tables for its non-redundant ALL set
# (Li et al., Nucleic Acids Res. 2025, 53:D179; https://plasmid.deepomics.org/download).
#
# protein_list is the pipeline input (config input.plasmidscope_proteins).
#
# The API refuses HEAD requests and needs the trailing slash. Tables are streamed through
# gzip. Record the download date in the comment above input.plasmidscope_proteins in
# config/config.yaml.
set -euo pipefail

BASE="https://plasmidapi.deepomics.org/api/database/files/ALL/data"
OUT="data/PlasmidScope/annotation"
mkdir -p "$OUT"

for table in protein_list; do
  dest="$OUT/ALL.$table.tsv.gz"
  if [[ -s "$dest" ]]; then
    echo "[skip] $dest"
    continue
  fi
  curl -fsSL --retry 3 "$BASE/ALL.$table.tsv/" | gzip -6 > "$dest.part"
  mv "$dest.part" "$dest"
  echo "[ok  ] $dest ($(stat -c %s "$dest") bytes)"
done
