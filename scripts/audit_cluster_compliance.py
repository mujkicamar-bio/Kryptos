#!/usr/bin/env python
"""Do the members of an MMseqs2 clustering actually meet the criterion it claims?

MMseqs2's own help states that cascaded clustering "can cluster sequence that do not fulfill
the clustering criteria", and offers --cluster-reassign to correct it. This script measures how
often that happens, by aligning sampled members back to their OWN cluster representative and
testing identity and coverage directly -- no transitive links allowed.

Usage:
    audit_cluster_compliance.py CLUSTER_TSV FASTA MIN_SEQ_ID [--cov 0.8] [--n 200] [--seed 0]

Prints the compliant fraction and a breakdown of how non-compliant members fail.
"""
import argparse
import random
import subprocess
import sys
import tempfile
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
MMSEQS = "/gorilla/home/amujkic/.conda/envs/panaroo/bin/mmseqs"


def read_faa(path):
    out, name, cur = {}, None, []
    for line in open(path):
        if line[0] == ">":
            if name is not None:
                out[name] = "".join(cur)
            name, cur = line[1:].strip().split()[0], []
        else:
            cur.append(line.strip())
    out[name] = "".join(cur)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cluster_tsv")
    ap.add_argument("fasta")
    ap.add_argument("min_seq_id", type=float)
    ap.add_argument("--cov", type=float, default=0.8)
    ap.add_argument("--n", type=int, default=200)
    ap.add_argument("--seed", type=int, default=0)
    a = ap.parse_args()

    cl = pd.read_csv(a.cluster_tsv, sep="\t", header=None, names=["rep", "member"])
    sizes = cl.groupby("rep").size()
    multi = sizes[sizes >= 2].index.tolist()
    random.seed(a.seed)
    picked = set(random.sample(multi, min(a.n, len(multi))))
    sub = cl[cl["rep"].isin(picked)]
    # Members only -- a representative trivially matches itself at 100%.
    sub = sub[sub["rep"] != sub["member"]]

    seqs = read_faa(a.fasta)
    with tempfile.TemporaryDirectory(dir="/scratch") as td:
        td = Path(td)
        with open(td / "q.faa", "w") as fh:
            for m in sub["member"]:
                fh.write(f">{m}\n{seqs[m]}\n")
        with open(td / "t.faa", "w") as fh:
            for r in sorted(picked):
                fh.write(f">{r}\n{seqs[r]}\n")
        # Permissive search: we want the alignment even when it fails the criterion, so the
        # thresholds are applied afterwards in pandas rather than by mmseqs.
        subprocess.run(
            [MMSEQS, "easy-search", str(td / "q.faa"), str(td / "t.faa"), str(td / "hits.tsv"),
             str(td / "tmp"), "--min-seq-id", "0.0", "-c", "0.0", "-s", "7.5", "-e", "1000",
             "--max-seqs", "5000", "--threads", "8", "-v", "1", "--format-output",
             "query,target,fident,qcov,tcov,evalue"],
            check=True, stdout=subprocess.DEVNULL)
        h = pd.read_csv(td / "hits.tsv", sep="\t", header=None,
                        names=["member", "rep", "fident", "qcov", "tcov", "evalue"])

    own = sub.merge(h, on=["member", "rep"], how="left")
    aligned = own["fident"].notna()
    ok_id = own["fident"] >= a.min_seq_id
    ok_cov = (own["qcov"] >= a.cov) & (own["tcov"] >= a.cov)
    ok = aligned & ok_id & ok_cov

    print(f"{a.cluster_tsv}")
    print(f"  criterion: identity >= {a.min_seq_id}, coverage >= {a.cov} both ways")
    print(f"  sampled {len(picked):,} families with >=2 members -> {len(own):,} member-rep pairs")
    print(f"  COMPLIANT: {ok.sum():,} / {len(own):,}  ({ok.mean() * 100:.2f}%)")
    print(f"    no alignment to own rep : {(~aligned).sum():,}")
    print(f"    fails identity          : {(aligned & ~ok_id).sum():,}")
    print(f"    fails coverage          : {(aligned & ok_id & ~ok_cov).sum():,}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
