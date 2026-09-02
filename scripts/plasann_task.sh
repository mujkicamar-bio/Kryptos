#!/bin/bash
# Annotate ONE shard with PlasAnn, using node-local scratch for the transient per-plasmid files so the
# quota'd project FS only ever receives two small files per shard (compact TSV + GenBank tarball).
# Usage: plasann_task.sh <shard_id> [tmp_dir]   (tmp_dir defaults to $SNIC_TMP, else a mktemp dir)
set -uo pipefail

ROOT=/gorilla/proj/h-mel-phylo/h-mel-phylo/private/amar/plasmid_analyis
ENVBIN=/gorilla/home/amujkic/.conda/envs/plasann_env/bin
GPY=/gorilla/home/amujkic/.conda/envs/genesis/bin/python
export PATH="$ENVBIN:$PATH"

SID=$1
# shard set is parametrizable so a retry can point at the smaller re-shards:
#   SHARDS_SUBDIR (default shards) + SHARD_PREFIX (default shard)
SHARDS_SUBDIR=${SHARDS_SUBDIR:-shards}
SHARD_PREFIX=${SHARD_PREFIX:-shard}
SHARD=$(printf "%s_%04d" "$SHARD_PREFIX" "$SID")
TMP=${2:-${SNIC_TMP:-$(mktemp -d)}}
IN="$TMP/$SHARD/in"; OUT="$TMP/$SHARD/out"
mkdir -p "$IN" "$OUT" "$ROOT/data/plasann_run/annot" "$ROOT/data/plasann_run/gbk"

SHARD_GZ="$ROOT/data/plasann_run/$SHARDS_SUBDIR/$SHARD.fna.gz"
[ -s "$SHARD_GZ" ] || { echo "[$SHARD] missing shard gz"; exit 1; }

# expand the shard's multi-FASTA into per-plasmid files ON NODE-LOCAL DISK
zcat "$SHARD_GZ" | awk -v D="$IN" '
  /^>/ { id=substr($1,2); gsub(/[^A-Za-z0-9._-]/,"_",id); f=D"/"id".fasta"; print > f; next }
  f    { print >> f }'
NIN=$(ls "$IN" | wc -l)
echo "[$SHARD] host=$(hostname) start=$(date +%F_%T) inputs=$NIN tmp=$TMP"

# annotate FROM node-local scratch as CWD: PlasAnn writes its temp_dir_blast_* into the working
# directory, so 300 concurrent tasks must NOT share the project dir as CWD (collisions + inode quota).
cd "$TMP/$SHARD"
PlasAnn -i "$IN" -o "$OUT" -t fasta > "$TMP/$SHARD/plasann.log" 2>&1
cd "$ROOT"
NANN=$(find "$OUT" -name '*_annotations.csv' | wc -l)
# surface PlasAnn errors into the SLURM log if anything failed
[ "$NANN" -lt "$NIN" ] && { echo "[$SHARD] --- plasann.log tail ---"; tail -20 "$TMP/$SHARD/plasann.log"; }

# compact per-feature TSV (drops Translation) + GenBank tarball, written back to project FS
"$GPY" "$ROOT/scripts/plasann_harvest_shard.py" "$OUT" "$ROOT/data/plasann_run/annot/$SHARD.tsv.gz"
( cd "$OUT" && tar -czf "$ROOT/data/plasann_run/gbk/$SHARD.gbk.tar.gz" */*_genbank.gbk 2>/dev/null ) || true

echo "[$SHARD] done=$(date +%F_%T) inputs=$NIN annotated=$NANN"
[ "$NANN" -eq "$NIN" ] || echo "[$SHARD] WARN: $((NIN-NANN)) plasmid(s) produced no annotation"
