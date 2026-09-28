"""Seeded, stratified sample (topology x size band) of the analysis set, excluding
simulated records; writes one FASTA for config input.fasta."""
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
            if len(f) <= max(col["topology"], col["size_bp"], col["hab_top"]):
                continue
            if f[col["hab_top"]] == args.exclude_hab_top:
                continue
            try:
                size = int(f[col["size_bp"]])
            except ValueError:
                continue
            rows.append((f[col["plasmid_id"]], f[col["topology"]].strip().lower(), size))

    if not rows:
        sys.exit(f"no usable rows in {args.master}")

    # ---- stratify: topology x size band ----------------------------------------------
    # The band only spreads the sample over sizes and decides nothing downstream. It is not
    # the pipeline's small-plasmid cut-off (input.max_plasmid_size_bp, 20 kb); it stays at
    # 10 kb so that the seeded command reproduces the existing test set.
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
