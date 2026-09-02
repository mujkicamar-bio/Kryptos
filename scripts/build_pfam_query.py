#!/usr/bin/env python
"""Build the Pfam query set for the dark-plasmidome novelty test, and chunk it for SLURM.

Three sets go into one FASTA, tagged in the header as `{SET}::{original_id}`:

  DARKREP  the 92,750 dark-family representatives (dark30_rep_seq.fasta)
           -> the test set: do the dark families match any Pfam-A domain?

  NAMED    all 78,545 proteins PlasAnn *did* name
           -> POSITIVE CONTROL. Without this a low DARKREP hit rate is uninterpretable:
              it could mean the dark families are novel, or merely that plasmid proteins
              in general are poorly covered by Pfam. The control separates the two.

  VALID    every member of a sample of dark families with 10-50 members
           -> validates propagating a representative's Pfam status to its whole family.
              If members disagree with their representative, family-level claims are unsafe.

Usage: build_pfam_query.py [--chunk-size 2000] [--valid-families 200] [--seed 20260828]
"""
import argparse
import random
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
D = ROOT / "data" / "dark_orf_run"
OUT = ROOT / "data" / "pfam_run"


def read_fasta(path):
    name, buf = None, []
    for line in open(path):
        if line.startswith(">"):
            if name:
                yield name, "".join(buf)
            name, buf = line[1:].strip(), []
        else:
            buf.append(line.strip())
    if name:
        yield name, "".join(buf)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--chunk-size", type=int, default=2000)
    ap.add_argument("--valid-families", type=int, default=200)
    ap.add_argument("--seed", type=int, default=20260828)
    args = ap.parse_args()

    OUT.mkdir(parents=True, exist_ok=True)
    chunks = OUT / "chunks"
    chunks.mkdir(exist_ok=True)
    for old in chunks.glob("*.faa"):
        old.unlink()

    records = []

    # --- DARKREP -------------------------------------------------------------
    dark_reps = {n: s for n, s in read_fasta(D / "dark30_rep_seq.fasta")}
    records += [(f"DARKREP::{n}", s) for n, s in dark_reps.items()]
    print(f"DARKREP {len(dark_reps):,}")

    # --- NAMED (positive control) -------------------------------------------
    n_named = 0
    for n, s in read_fasta(D / "all_cds.faa"):
        if n.endswith("|N"):
            records.append((f"NAMED::{n}", s))
            n_named += 1
    print(f"NAMED   {n_named:,}")

    # --- VALID (propagation check) ------------------------------------------
    members = defaultdict(list)
    for line in open(D / "dark30_cluster.tsv"):
        rep, mem = line.rstrip("\n").split("\t")
        members[rep].append(mem)
    eligible = [r for r, m in members.items() if 10 <= len(m) <= 50]
    rng = random.Random(args.seed)
    picked = rng.sample(eligible, min(args.valid_families, len(eligible)))
    dark_seqs = {n: s for n, s in read_fasta(D / "dark_orfs.faa")}
    n_valid = 0
    for rep in picked:
        for mem in members[rep]:
            records.append((f"VALID::{rep}::{mem}", dark_seqs[mem]))
            n_valid += 1
    print(f"VALID   {n_valid:,} members from {len(picked)} families "
          f"(of {len(eligible):,} eligible, seed {args.seed})")

    # --- write chunks --------------------------------------------------------
    n = 0
    for i in range(0, len(records), args.chunk_size):
        part = records[i:i + args.chunk_size]
        with open(chunks / f"q{i // args.chunk_size:04d}.faa", "w") as fh:
            for name, seq in part:
                fh.write(f">{name}\n{seq}\n")
        n += 1
    (OUT / "query_manifest.txt").write_text(
        f"DARKREP\t{len(dark_reps)}\nNAMED\t{n_named}\nVALID\t{n_valid}\n"
        f"TOTAL\t{len(records)}\nCHUNKS\t{n}\nCHUNK_SIZE\t{args.chunk_size}\n"
        f"VALID_FAMILIES\t{len(picked)}\nSEED\t{args.seed}\n")
    print(f"\nTOTAL   {len(records):,} sequences -> {n} chunks of {args.chunk_size}")
    print(f"array range: 0-{n - 1}")


if __name__ == "__main__":
    main()
