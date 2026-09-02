#!/bin/bash
# One chunk of the DefenseFinder run. Called by scripts/defensefinder_array.sbatch.
set -uo pipefail
i=$(printf "%04d" "$1")
ROOT=/gorilla/proj/h-mel-phylo/h-mel-phylo/private/amar/plasmid_analyis
export PATH=/gorilla/home/amujkic/.conda/envs/panaroo/bin:$PATH   # hmmsearch
cd "$ROOT"

Q="data/defensefinder_run/chunks/c$i.faa"
O="data/defensefinder_run/out/c$i"
[ -s "$Q" ] || { echo "missing $Q"; exit 1; }
rm -rf "$O.partial" && mkdir -p "$O.partial"

# -a  : also run AntiDefenseFinder (anti-CRISPR / anti-restriction / anti-SOS models).
#       This is an INDEPENDENT second opinion on the dbAPIS result in antidefense.ipynb.
# --db-type gembase : whole plasmids in gene order, so MacSyFinder can use genomic context.
envs/defensefinder/bin/defense-finder run \
    -o "$O.partial" -w 8 --db-type gembase -a "$Q" || exit 1
mv "$O.partial" "$O"
echo "chunk $i done: $(( $(wc -l < "$O/c${i}_defense_finder_systems.tsv") - 1 )) systems"
