"""Build a small, representative plasmid set for exercising the pipeline end to end.

WHY A FIXED SAMPLE RATHER THAN THE FIRST N RECORDS

The first N records of the working set are whatever order the FASTA happens to be in,
which is neither representative nor stable. This selects a STRATIFIED sample so that the
properties the pipeline branches on are all present:

  * topology       circular and linear both appear, because origin repair runs on one and
                   must not run on the other (spec section 8.4);
  * size           small cryptic plasmids and large ones both appear, because the plus or
                   minus three neighbourhood on a six-gene plasmid is the whole molecule
                   and that is the statistical trap S8 exists to avoid;
  * the locked exclusion is applied, so the sample cannot contain a simulated record.

The selection is seeded, so the same command reproduces the same sample.

OUTPUT

One FASTA, ready for config.input.fasta.
"""
import argparse
import collections
import gzip
import pathlib
import random
import sys


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--master", required=True, help="the analysis-set TSV")
    ap.add_argument("--fasta", required=True, help="the working-set FASTA, gzipped")
    ap.add_argument("--out", required=True, help="the FASTA to write")
    ap.add_argument("--n", type=int, default=100, help="plasmids to select")
    ap.add_argument("--seed", type=int, default=20260917)
    ap.add_argument("--exclude-hab-top", default="Simulated-artifact")
    args = ap.parse_args()

    # ---- read the master table, keeping only what the stratification needs -----------
    rows = []
    with open(args.master) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        col = {name: i for i, name in enumerate(header)}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) <= max(col["topology"], col["size_bp"]):
                continue
            if f[col.get("hab_top", 0)] == args.exclude_hab_top:
                continue
            try:
                size = int(f[col["size_bp"]])
            except ValueError:
                continue
            rows.append((f[col["plasmid_id"]], f[col["topology"]].strip().lower(), size))

    if not rows:
        sys.exit(f"no usable rows in {args.master}")

    # ---- stratify: topology x size band ----------------------------------------------
    # Small is below 10 kb, which is where the cryptic plasmids the project cares about
    # sit; the band boundary is descriptive here and decides nothing downstream.
    def band(size):
        return "small" if size < 10_000 else "large"

    strata = collections.defaultdict(list)
    for pid, topology, size in rows:
        strata[(topology or "unknown", band(size))].append(pid)

    rng = random.Random(args.seed)
    # Take round-robin across strata so every stratum present contributes before any one is
    # taken twice. A proportional sample would give 98 circular and 2 linear, and the linear
    # branch - the one that must NOT get origin repair - would barely be exercised.
    order = sorted(strata)
    for key in order:
        rng.shuffle(strata[key])
    chosen, i = [], 0
    while len(chosen) < args.n and any(strata[k] for k in order):
        key = order[i % len(order)]
        if strata[key]:
            chosen.append(strata[key].pop())
        i += 1

    wanted = set(chosen)
    print(f"selected {len(wanted)} plasmids from {len(rows)} eligible")
    for key in order:
        n = sum(1 for pid in chosen if pid in set(strata[key]) | {pid} and True)
    counts = collections.Counter()
    for pid, topology, size in rows:
        if pid in wanted:
            counts[(topology or "unknown", band(size))] += 1
    for key in sorted(counts):
        print(f"  {key[0]:<24} {key[1]:<6} {counts[key]}")

    # ---- stream the FASTA once, writing each wanted record ---------------------------
    out = pathlib.Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    written, keep = 0, False
    with gzip.open(args.fasta, "rt") as fh, open(out, "w") as fout:
        for line in fh:
            if line[0] == ">":
                pid = line[1:].split()[0]
                keep = pid in wanted
                if keep:
                    written += 1
                    wanted.discard(pid)
            if keep:
                fout.write(line)

    print(f"wrote {written} records to {out}")
    if wanted:
        # A selected plasmid absent from the FASTA means the master table and the sequence
        # file disagree, which every downstream count would inherit silently.
        sys.exit(f"{len(wanted)} selected plasmids were not found in {args.fasta}, "
                 f"e.g. {sorted(wanted)[:3]}")


if __name__ == "__main__":
    main()
