#!/bin/bash
# MOB-suite mob_typer on ONE shard: replicon(rep_type), relaxase, MPF, oriT, mobility, MOB cluster,
# host range. Uses the self-contained mobtyper_env (its own blast+mash+DB). Node-local temp.
# Usage: mobtyper_task.sh <shard_id>
set -uo pipefail
ROOT=/gorilla/proj/h-mel-phylo/h-mel-phylo/private/amar/plasmid_analyis
NE=$ROOT/envs/mobtyper_env
DB=$(ls -d "$NE"/lib/python*/site-packages/mob_suite/databases | head -1)
export PATH="$NE/bin:$PATH"

SID=$1
SHARD=$(printf "shard_%04d" "$SID")
NLOCAL="${SNIC_TMP:-$(mktemp -d)}"
TMP="$NLOCAL/$SHARD"
export TMPDIR="$TMP"                 # mob_typer writes blast/mash intermediates to $TMPDIR
mkdir -p "$TMP" "$ROOT/data/typing/mob"

# Copy the MOB database to NODE-LOCAL disk. The shared-FS taxa.sqlite (host-range) fails under 300
# concurrent readers (sqlite disk I/O error) and network-FS sqlite is pathologically slow (~81s/plasmid).
# A private local copy fixes both. flock => one copy per node (shared $SNIC_TMP), others reuse it.
LOCALDB="$NLOCAL/mobdb"
(
  flock 9
  if [ ! -f "$LOCALDB/status.txt" ]; then
    rm -rf "$LOCALDB.tmp"; cp -r "$DB" "$LOCALDB.tmp" && mv "$LOCALDB.tmp" "$LOCALDB"
  fi
) 9>"$NLOCAL/.mobdb.lock"
DB="$LOCALDB"

zcat "$ROOT/data/plasann_run/shards/$SHARD.fna.gz" > "$TMP/in.fna"
NIN=$(grep -c '^>' "$TMP/in.fna")
cd "$TMP"
echo "[$SHARD] host=$(hostname) start=$(date +%T) inputs=$NIN"

mob_typer --multi --num_threads 1 -d "$DB" \
  --infile "$TMP/in.fna" --out_file "$ROOT/data/typing/mob/$SHARD.txt" > "$TMP/log" 2>&1
RC=$?
NOUT=$(tail -n +2 "$ROOT/data/typing/mob/$SHARD.txt" 2>/dev/null | wc -l)
echo "[$SHARD] done=$(date +%T) rc=$RC inputs=$NIN typed=$NOUT"
if [ "$NOUT" -lt "$NIN" ]; then echo "[$SHARD] --- log tail ---"; tail -15 "$TMP/log"; fi
exit 0
