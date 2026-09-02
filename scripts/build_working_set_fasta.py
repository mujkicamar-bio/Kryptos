#!/usr/bin/env python3
"""Consolidate the 208,248 working-set plasmid sequences into one FASTA (gzipped).

PlasmidScope deduplicated at 100% identity/coverage, so every source member of a plasmid is the same
sequence — we pull ONE per working-set plasmid, headed by its rep id. Sequences come from the per-source
PlasmidScope FASTA archives (header first-token = the exact provenance id, e.g. `>GenBank_CP028544.1`)
and, for IMG/PR plasmids (no archive), from the IMG/PR complete FASTA (`>IMGPR_plasmid_..|taxon|scaf`).

This is the sequence file the replicon-typing run needs; it lets us delete the 9.5 GB of raw archives.

Inputs:  data/plasmidscope_primary/{working_set.tsv, complete_provenance.tsv}
         data/PlasmidScope/fasta/*.fasta.tar.gz   (per-source)
         data/processed/complete_plasmids.fna     (IMG/PR)
Output:  data/plasmidscope_primary/working_set.fna.gz   + a .missing.txt if any rep unresolved
Run:     python3 scripts/build_working_set_fasta.py
"""
import csv
import gzip
import os
import subprocess
import sys

csv.field_size_limit(sys.maxsize)
PP = "data/plasmidscope_primary"
WS = f"{PP}/working_set.tsv"
PROV = f"{PP}/complete_provenance.tsv"
IMGPR_FNA = "data/processed/complete_plasmids.fna"
ARCH_DIR = "data/PlasmidScope/fasta"
OUT = f"{PP}/working_set.fna.gz"
# process order: IMG/PR first (covers ~136k), then archives for the remaining isolate/mixed reps
ARCHIVES = ["RefSeq", "GenBank", "PLSDB", "COMPASS", "DDBJ", "ENA", "Kraken2", "mMGE"]


def main():
    reps = {r["plasmid_id"] for r in csv.DictReader(open(WS, newline=""), delimiter="\t")}
    print(f"working-set reps: {len(reps):,}", flush=True)
    # source id -> rep (only for working-set reps)
    id2rep = {}
    for r in csv.DictReader(open(PROV, newline=""), delimiter="\t"):
        rep = r["plasmid_id"].split(",")[0]
        if rep in reps:
            for sid in r["plasmid_id"].split(","):
                id2rep[sid] = rep
    print(f"source ids mapped: {len(id2rep):,}", flush=True)

    written = set()
    # pipe through system `gzip -1` (fast) instead of Python's slow default-level gzip
    outfh = open(OUT, "wb")
    gz = subprocess.Popen(["gzip", "-1"], stdin=subprocess.PIPE, stdout=outfh, text=True)
    out = gz.stdin

    def consume(stream, id_of_header):
        keep = False
        for line in stream:
            if line.startswith(">"):
                hid = id_of_header(line[1:])
                rep = id2rep.get(hid)
                if rep is not None and rep not in written:
                    written.add(rep)
                    out.write(">" + rep + "\n")
                    keep = True
                else:
                    keep = False
            elif keep:
                out.write(line)

    # 1) IMG/PR complete FASTA (header: IMGPR_plasmid_..|taxon|scaffold -> id before first '|')
    if os.path.exists(IMGPR_FNA):
        with open(IMGPR_FNA) as fh:
            consume(fh, lambda h: h.split("|", 1)[0].strip())
        print(f"after IMG/PR: {len(written):,} written", flush=True)

    # 2) per-source archives. Each archive stores ONE .fasta file per sequence, so stream the WHOLE
    # archive to stdout (`tar -xzO` with no member = all members concatenated) — a valid multi-fasta.
    # header: '<id> optional description' -> first whitespace token.
    for src in ARCHIVES:
        arch = f"{ARCH_DIR}/{src}.fasta.tar.gz"
        if not os.path.exists(arch):
            continue
        p = subprocess.Popen(["tar", "-xzOf", arch], stdout=subprocess.PIPE, text=True, bufsize=1 << 20)
        consume(p.stdout, lambda h: h.split()[0].strip() if h.split() else "")
        p.wait()
        print(f"after {src}: {len(written):,} written", flush=True)

    out.close()
    gz.wait()
    outfh.close()
    missing = reps - written
    print(f"\nDONE: {len(written):,}/{len(reps):,} sequences written -> {OUT}")
    if missing:
        with open(f"{PP}/working_set.fna.missing.txt", "w") as fh:
            fh.write("\n".join(sorted(missing)) + "\n")
        print(f"  ⚠ {len(missing):,} reps had no sequence (listed in working_set.fna.missing.txt)")
    else:
        print("  all reps resolved.")


if __name__ == "__main__":
    main()
