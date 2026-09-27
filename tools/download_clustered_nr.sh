#!/usr/bin/env bash
# NCBI ClusteredNR as the T5 DIAMOND database (config/cascade.yaml, tier T5).
#
# ClusteredNR is NCBI's nr clustered at 90% identity and 90% length, one representative per
# cluster, published as a BLAST database only
# (https://ftp.ncbi.nlm.nih.gov/blast/db/clustered_nr-prot-metadata.json). DIAMOND needs
# FASTA, so the representatives are exported with blastdbcmd and indexed with makedb.
#
#   tools/download_clustered_nr.sh download   login node (needs the network); resumable,
#                                             every volume checked against NCBI's md5
#   tools/download_clustered_nr.sh build      compute node; export + diamond makedb
#
# Each step skips what is already complete. Record the metadata's last-updated date as
# references.nr_version in config/config.yaml.
set -euo pipefail

OUT="data/refs/clustered_nr"
DB="$OUT/blastdb"
BASE="https://ftp.ncbi.nlm.nih.gov/blast/db"
THREADS="${SLURM_CPUS_PER_TASK:-8}"
mkdir -p "$DB"

download() {
  curl -fsSL --retry 3 "$BASE/clustered_nr-prot-metadata.json" > "$OUT/metadata.json"
  grep -E '"(last-updated|number-of-sequences|number-of-letters)"' "$OUT/metadata.json"
  # The volume list is the metadata's own; ftp:// links are fetched over https.
  grep -o 'clustered_nr\.[0-9]*\.tar\.gz' "$OUT/metadata.json" | sort -u > "$OUT/volumes.txt"
  echo "$(wc -l < "$OUT/volumes.txt") volumes"
  while read -r vol; do
    if [[ -e "$DB/.$vol.done" ]]; then continue; fi
    echo "$BASE/$vol"; echo "  out=$vol"
    echo "$BASE/$vol.md5"; echo "  out=$vol.md5"
  done < "$OUT/volumes.txt" > "$OUT/urls.txt"
  if [[ -s "$OUT/urls.txt" ]]; then
    # NCBI answers 503 to more than a few connections at once, so four, one each, and
    # every failure retried rather than abandoned.
    aria2c -d "$DB" -i "$OUT/urls.txt" -j 4 -x 1 -c --auto-file-renaming=false \
      --allow-overwrite=true --max-tries=0 --retry-wait=60 \
      --console-log-level=warn --summary-interval=600
  fi
  # Verify, unpack, and drop the archive; a marker makes the step resumable.
  while read -r vol; do
    if [[ -e "$DB/.$vol.done" ]]; then continue; fi
    (cd "$DB" && md5sum -c --quiet "$vol.md5")
    tar -xzf "$DB/$vol" -C "$DB"
    rm "$DB/$vol" "$DB/$vol.md5"
    touch "$DB/.$vol.done"
  done < "$OUT/volumes.txt"
  echo "download complete: $(ls "$DB"/.clustered_nr.*.done | wc -l) volumes"
}

build() {
  n_vol=$(wc -l < "$OUT/volumes.txt")
  n_done=$(ls "$DB"/.clustered_nr.*.done 2>/dev/null | wc -l)
  if [[ "$n_done" -ne "$n_vol" ]]; then
    echo "only $n_done of $n_vol volumes downloaded; run 'download' first" >&2
    exit 1
  fi
  if [[ ! -s "$OUT/clustered_nr.dmnd" ]]; then
    # Plain FASTA, one record per representative: accession and title in the header,
    # which is what the cascade parses (plasmidann.labels.parse_ncbi_title).
    # One blastdbcmd per volume, in parallel: over the whole alias database it stalled at
    # ~0.15 MB/s after 27 GB (job 6971996); one volume exports in ~41 s at 1.6 GB RSS.
    mkdir -p "$OUT/faa_parts"
    sed 's/\.tar\.gz$//' "$OUT/volumes.txt" | xargs -P "$THREADS" -I{} sh -c \
      'blastdbcmd -db "$1/$2" -entry all -outfmt %f > "$3/$2.faa"' _ "$DB" {} "$OUT/faa_parts"
    cat "$OUT"/faa_parts/clustered_nr.*.faa > "$OUT/clustered_nr.faa"
    rm -r "$OUT/faa_parts"
    # makedb writes to a temporary name, so a killed build cannot leave a truncated
    # clustered_nr.dmnd that the -s test above would accept as finished.
    diamond makedb --in "$OUT/clustered_nr.faa" -d "$OUT/clustered_nr.tmp" \
      --threads "$THREADS"
    mv "$OUT/clustered_nr.tmp.dmnd" "$OUT/clustered_nr.dmnd"
    rm "$OUT/clustered_nr.faa"
  fi
  diamond dbinfo -d "$OUT/clustered_nr.dmnd"
}

case "${1:-}" in
  download) download ;;
  build) build ;;
  *) echo "usage: $0 download|build" >&2; exit 2 ;;
esac
