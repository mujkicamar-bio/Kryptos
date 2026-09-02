#!/usr/bin/env python
"""Extract CDS protein translations from the PlasAnn GenBank output.

PlasAnn writes one GenBank file per plasmid, tarred per shard, in
data/plasann_run/gbk/{shard,mshard}_*.gbk.tar.gz. Every CDS carries /category=
and /translation=, so the unannotatable ORFs ("Open reading frame") have protein
sequence just like the named genes -- this is the only place in the pipeline
where the dark proteome is recoverable.

Usage:
    extract_plasann_proteins.py --ids IDS.txt --out OUTDIR [--mode dark|all]
                                [--jobs N] [--limit N]

--mode dark : only CDS with category "Open reading frame"
              FASTA header  >{plasmid_id}|{cds_index}
--mode all  : every CDS that has a category and a translation
              FASTA header  >{plasmid_id}|{cds_index}|{D|N}   (D = dark ORF, N = named)

cds_index counts EVERY translated CDS on the plasmid, in file order, regardless of mode.
This is load-bearing: it makes the two outputs joinable, so a dark ORF has the same id in
both files. (Indexing only the emitted subset silently breaks that join for any plasmid
carrying a named CDS -- the failure is invisible because ids still look well-formed.)

Writes one .faa per shard into OUTDIR (shards contributing nothing are skipped).
"""
import argparse
import glob
import os
import tarfile
from multiprocessing import Pool
from pathlib import Path

GBK_DIR = Path(__file__).resolve().parent.parent / "data" / "plasann_run" / "gbk"
DARK_CATEGORY = "Open reading frame"

_TARGET: set[str] = set()
_MODE = "dark"
_OUT: Path = Path(".")


def parse_stream(fh, target, mode):
    """Yield (plasmid_id, idx, category, protein) from a concatenated GenBank byte stream."""
    pid = None
    keep = False
    cat = trans = None
    in_trans = False
    idx = 0

    for raw in fh:
        line = raw.decode("utf-8", "replace").rstrip("\n")
        if line.startswith("LOCUS"):
            pid = line.split()[1]
            keep = pid in target
            idx = 0
            cat = trans = None
            in_trans = False
            continue
        if not keep:
            continue
        if in_trans:
            trans += line.strip().rstrip('"')
            if line.rstrip().endswith('"'):
                in_trans = False
                if cat is not None and trans:
                    idx += 1
                    yield pid, idx, cat, trans.rstrip("*")
                cat = trans = None
            continue
        stripped = line.strip()
        if line.startswith("     CDS "):
            cat = trans = None
        elif stripped.startswith('/category="'):
            cat = stripped[11:].rstrip('"')
        elif stripped.startswith('/translation="'):
            trans = stripped[14:].rstrip('"')
            if not stripped.endswith('"') or len(stripped) == 14:
                in_trans = True
            else:
                if cat is not None and trans:
                    idx += 1
                    yield pid, idx, cat, trans.rstrip("*")
                cat = trans = None


def do_shard(tarpath):
    name = os.path.basename(tarpath).replace(".gbk.tar.gz", "")
    outp = _OUT / f"{name}.faa"
    n = 0
    with tarfile.open(tarpath, "r:gz") as tf, open(outp, "w") as out:
        for member in tf:
            if not member.isfile():
                continue
            fh = tf.extractfile(member)
            if fh is None:
                continue
            for pid, idx, cat, prot in parse_stream(fh, _TARGET, _MODE):
                if _MODE == "dark":
                    if cat != DARK_CATEGORY:
                        continue
                    out.write(f">{pid}|{idx}\n{prot}\n")
                else:
                    tag = "D" if cat == DARK_CATEGORY else "N"
                    out.write(f">{pid}|{idx}|{tag}\n{prot}\n")
                n += 1
    if n == 0:
        outp.unlink()
    return name, n


def _init(target, mode, out):
    global _TARGET, _MODE, _OUT
    _TARGET, _MODE, _OUT = target, mode, out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ids", required=True, help="file with one plasmid_id per line")
    ap.add_argument("--out", required=True, help="output directory for per-shard .faa")
    ap.add_argument("--mode", choices=["dark", "all"], default="dark")
    ap.add_argument("--jobs", type=int, default=24)
    ap.add_argument("--limit", type=int, help="only process the first N shards (smoke test)")
    args = ap.parse_args()

    target = {line.strip() for line in open(args.ids) if line.strip()}
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    tars = sorted(glob.glob(str(GBK_DIR / "*.gbk.tar.gz")))
    if args.limit:
        tars = tars[: args.limit]
    print(f"{len(tars)} shard tarballs, {len(target):,} target plasmids, mode={args.mode}",
          flush=True)

    total = 0
    with Pool(args.jobs, initializer=_init, initargs=(target, args.mode, out)) as pool:
        for i, (_, n) in enumerate(pool.imap_unordered(do_shard, tars), 1):
            total += n
            if i % 200 == 0 or i == len(tars):
                print(f"  {i}/{len(tars)} shards, {total:,} proteins", flush=True)
    print(f"TOTAL {total:,}")


if __name__ == "__main__":
    main()
