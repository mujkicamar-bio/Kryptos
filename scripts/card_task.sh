#!/bin/bash
# RGI (CARD) antibiotic-resistance annotation on ONE shard of the working set.
#   rgi main --input_type contig -a DIAMOND   against CARD v4.0.1 (local DB in data/card/localDB).
# Design mirrors mobtyper_task.sh: copy the CARD localDB to NODE-LOCAL disk once per node (flock),
# run RGI with output + all intermediates on node-local $SNIC_TMP so the quota'd project FS never
# holds RGI's transient files, then copy back only the final .txt. --include_loose keeps every
# cut-off (Perfect/Strict/Loose) in the raw table; the master later uses Perfect+Strict only.
# Usage: card_task.sh <shard_id>
set -uo pipefail
ROOT=/gorilla/proj/h-mel-phylo/h-mel-phylo/private/amar/plasmid_analyis
export MAMBA_ROOT_PREFIX="$ROOT/envs/.micromamba_root"
ENV="$MAMBA_ROOT_PREFIX/envs/rgi_env"
export PATH="$ENV/bin:$PATH"
export PYTHONNOUSERSITE=1                 # do not let ~/.local shadow the env's packages
SRCDB="$ROOT/data/card/localDB"

SID=$1
SHARD=$(printf "shard_%04d" "$SID")
NLOCAL="${SNIC_TMP:-$(mktemp -d)}"
TMP="$NLOCAL/$SHARD"
mkdir -p "$TMP" "$ROOT/data/card_run/rgi"

# Copy CARD localDB to node-local disk, once per node (flock). RGI --local reads ./localDB from CWD.
# The prebuilt DIAMOND index (protein.db.dmnd) makes localDB read-only at run time, so one shared
# node-local copy serves every concurrent task on the node.
LOCALDB="$NLOCAL/localDB"
(
  flock 9
  if [ ! -f "$LOCALDB/protein.db.dmnd" ]; then
    rm -rf "$LOCALDB.tmp"; cp -r "$SRCDB" "$LOCALDB.tmp" && mv "$LOCALDB.tmp" "$LOCALDB"
  fi
) 9>"$NLOCAL/.carddb.lock"

zcat "$ROOT/data/plasann_run/shards/$SHARD.fna.gz" > "$TMP/in.fna"
NIN=$(grep -c '^>' "$TMP/in.fna")
cd "$NLOCAL"                               # CWD must contain ./localDB for rgi --local
echo "[$SHARD] host=$(hostname) start=$(date +%T) inputs=$NIN"

# Output + intermediates (<out>.temp) all node-local; only the final .txt is copied back.
rgi main -i "$TMP/in.fna" -o "$TMP/out" --input_type contig \
  -a DIAMOND -n "${SLURM_CPUS_PER_TASK:-2}" --clean --local -d plasmid --include_loose \
  > "$TMP/log" 2>&1
RC=$?
[ -f "$TMP/out.txt" ] && cp "$TMP/out.txt" "$ROOT/data/card_run/rgi/$SHARD.txt"

NROW=$(tail -n +2 "$ROOT/data/card_run/rgi/$SHARD.txt" 2>/dev/null | wc -l)
NPS=$(tail -n +2 "$ROOT/data/card_run/rgi/$SHARD.txt" 2>/dev/null \
       | awk -F'\t' '$6=="Perfect"||$6=="Strict"{print $2}' | sort -u | wc -l)
echo "[$SHARD] done=$(date +%T) rc=$RC inputs=$NIN rows=$NROW plasmids_with_PS_hit=$NPS"
if [ $RC -ne 0 ] || [ ! -f "$ROOT/data/card_run/rgi/$SHARD.txt" ]; then
  echo "[$SHARD] --- log tail ---"; tail -25 "$TMP/log"
fi
exit 0
