#!/bin/bash
# One chunk of the Pfam-A search. Called by scripts/pfam_array.sbatch.
set -uo pipefail
i=$(printf "%04d" "$1")
ROOT=/gorilla/proj/h-mel-phylo/h-mel-phylo/private/amar/plasmid_analyis
export PATH=/gorilla/home/amujkic/.conda/envs/panaroo/bin:$PATH

Q="$ROOT/data/pfam_run/chunks/q$i.faa"
O="$ROOT/data/pfam_run/domtbl/q$i.domtbl"
[ -s "$Q" ] || { echo "missing $Q"; exit 1; }

# --cut_ga : Pfam's curated per-family gathering thresholds. This is the standard
#            "is this a real Pfam match" criterion -- no arbitrary e-value picked by us.
# --noali  : we only need the tabular domain output.
hmmsearch --cut_ga --noali --cpu 4 \
    --domtblout "$O.partial" \
    "$ROOT/data/refs/pfam/Pfam-A.hmm" "$Q" > /dev/null || exit 1
mv "$O.partial" "$O"          # atomic: a truncated file never looks complete
echo "chunk $i done: $(grep -vc '^#' "$O") domain rows"
