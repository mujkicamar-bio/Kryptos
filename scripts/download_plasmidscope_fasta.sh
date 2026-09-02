#!/usr/bin/env bash
# Stage 2 download: per-source PlasmidScope nucleotide FASTA archives for the NON-IMG/PR datasets
# (IMG-PR is skipped -- its complete set is identical to our existing data/processed/complete_plasmids.fna).
# Archives contain the whole dataset (all completeness); complete-only subsetting happens after
# extraction in the Stage-2 dereplication step. See PLASMIDSCOPE_INTEGRATION_PLAN.md.
set -euo pipefail

BASE="https://plasmidapi.deepomics.org/api/database/files"
OUT="data/PlasmidScope/fasta"
mkdir -p "$OUT"

# smallest-first so we fail fast and gauge per-plasmid size before the big pulls
DATASETS=(Kraken2 DDBJ ENA COMPASS PLSDB RefSeq GenBank mMGE)

for ds in "${DATASETS[@]}"; do
  dest="$OUT/${ds}.fasta.tar.gz"
  if [[ -s "$dest" ]]; then
    echo "[skip] $ds present ($(du -h "$dest" | cut -f1))"
    continue
  fi
  url="$BASE/${ds}/${ds}.fasta.tar.gz"
  echo "[get ] $(date +%H:%M:%S) $ds <- $url"
  curl -fsSL -m 7200 "$url" -o "$dest.part"
  mv "$dest.part" "$dest"
  echo "[ok  ] $(date +%H:%M:%S) $ds -> $dest ($(du -h "$dest" | cut -f1))"
done
echo "[done] all non-IMG/PR FASTA archives in $OUT/"
