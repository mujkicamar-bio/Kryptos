#!/usr/bin/env bash
# Re-cluster the dark plasmid proteome at a given identity threshold.
#
# Differs from scripts/cluster_dark_orfs.sh in three ways, all deliberate:
#   1. --min-seq-id is an argument, so the same code produces the whole threshold sweep.
#   2. --cluster-reassign is set. MMseqs2's own help: cascaded clustering "can cluster
#      sequence that do not fulfill the clustering criteria"; reassignment corrects that.
#      Audited on the 30% run, 13.2% of one family's members violated the stated criterion.
#   3. Output goes to data/dark_orf_run/recluster/, never over dark30_*/mix30_*, which back
#      three notebooks and two reports.
#
# Coverage is held at -c 0.8 --cov-mode 0 throughout: 80% of BOTH sequences. That matches
# UniRef's 80% overlap rule (added 2013 to stop partial-sequence merging) and is already
# stricter than Rodriguez del Rio et al. Nature 2024 (-c 0.5 --cov-mode 1).
#
# Usage:  bash scripts/recluster_dark_orfs.sh <min_seq_id> [input.faa]
#   e.g.  bash scripts/recluster_dark_orfs.sh 0.5
set -euo pipefail

ID="${1:?usage: recluster_dark_orfs.sh <min_seq_id> [input.faa]}"
IN="${2:-data/dark_orf_run/dark_orfs.faa}"
THREADS="${THREADS:-24}"
COV="${COV:-0.8}"
export PATH=/gorilla/home/amujkic/.conda/envs/panaroo/bin:$PATH

TAG="dark$(printf '%.0f' "$(echo "$ID * 100" | bc -l)")"
OUT="data/dark_orf_run/recluster"
mkdir -p "$OUT/tmp_$TAG"

echo "[$(date +%H:%M:%S)] $TAG  min-seq-id=$ID  -c $COV --cov-mode 0 --cluster-reassign"
echo "  input: $IN  ($(grep -c '^>' "$IN") sequences)"

mmseqs easy-cluster "$IN" "$OUT/$TAG" "$OUT/tmp_$TAG" \
    --min-seq-id "$ID" -c "$COV" --cov-mode 0 --cluster-reassign 1 \
    --threads "$THREADS" -v 1

rm -rf "$OUT/tmp_$TAG"
# _all_seqs.fasta is a redundant re-emission of the input; drop it.
rm -f "$OUT/${TAG}_all_seqs.fasta"

echo "[$(date +%H:%M:%S)] $TAG done: $(cut -f1 "$OUT/${TAG}_cluster.tsv" | sort -u | wc -l) families"
