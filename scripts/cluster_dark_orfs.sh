#!/usr/bin/env bash
# Cluster the dark (unannotatable) plasmid proteome with MMseqs2.
#
# Two clusterings, both at 30% identity / 80% bidirectional coverage:
#   dark30 -- dark ORF proteins only        -> does the dark proteome recur?
#   mix30  -- dark + named proteins together -> do dark families have named homologs?
#
# Requires the `panaroo` conda env (mmseqs 18.8cc5c). Run from the project root
# AFTER scripts/extract_plasann_proteins.py has populated faa_dark/ and faa_all/.
set -euo pipefail

D=data/dark_orf_run
THREADS=${THREADS:-24}
export PATH=/gorilla/home/amujkic/.conda/envs/panaroo/bin:$PATH

# Deduplicate on the way in. 144 plasmids are present in BOTH PlasAnn shard series
# (shard_* and mshard_*), so globbing both -- which is required, or ~27% of the data is
# lost -- emits their CDS twice. Records are exactly two lines, so first-occurrence-wins
# on the header is safe and deterministic given the sorted glob.
dedup() { awk '/^>/{keep = !seen[$0]++} keep' "$@"; }
dedup "$D"/faa_dark/*.faa > "$D/dark_orfs.faa"
dedup "$D"/faa_all/*.faa  > "$D/all_cds.faa"

echo "after dedup: $(grep -c '^>' "$D/dark_orfs.faa") dark, $(grep -c '^>' "$D/all_cds.faa") all-CDS proteins"

rm -rf "$D/tmp_dark" "$D/tmp_mix" "$D"/dark30_* "$D"/mix30_*

mmseqs easy-cluster "$D/dark_orfs.faa" "$D/dark30" "$D/tmp_dark" \
    --min-seq-id 0.3 -c 0.8 --cov-mode 0 --threads "$THREADS" -v 1

mmseqs easy-cluster "$D/all_cds.faa" "$D/mix30" "$D/tmp_mix" \
    --min-seq-id 0.3 -c 0.8 --cov-mode 0 --threads "$THREADS" -v 1

rm -rf "$D/tmp_dark" "$D/tmp_mix"
# _all_seqs.fasta is a redundant re-emission of the input; drop it (~190 MB).
rm -f "$D"/dark30_all_seqs.fasta "$D"/mix30_all_seqs.fasta

echo "clustering done:"
ls -la "$D"/dark30_* "$D"/mix30_*
