"""S8b: locate integron cassette arrays with IntegronFinder.

THE STRONGEST PLASMID-SPECIFIC SIGNAL AVAILABLE, and one neither Nature study could use,
because both worked on genomes and metagenomes.

An integron captures gene cassettes, and a large fraction of cassettes are of unknown
function. A dark ORF in a cassette array is a real gene BY CONSTRUCTION: it carries an attC
recombination site, and it has been physically excised, mobilised and re-integrated - and
then retained. That is direct physical evidence of both existence and selection, obtained
without any homology at all.
"""
import concurrent.futures
import csv
import pathlib
import shutil
import subprocess

import _ctx  # noqa: F401

from plasmidann.fasta import split_fasta

fasta = pathlib.Path(snakemake.input.fasta)
outdir = pathlib.Path(snakemake.output[0]).parent / "integron_finder"
# A rerun starts clean: result files from an interrupted run would be read below.
shutil.rmtree(outdir, ignore_errors=True)
outdir.mkdir(parents=True)

# IntegronFinder walks replicons one at a time and threads only its HMM searches: the test
# run used 11.6 CPU-s in 57 s of wall time on 14 threads, which scales to ~25 h holding the
# whole node. So the analysis set is dealt into one chunk per core and each chunk runs on
# one thread; ~0.3 GB per process (measured), ~30 GB at 96.
chunks = split_fasta(fasta, snakemake.threads, outdir / "chunks")

# --local-max is the sensitive mode: it searches for attC sites beyond those adjacent to a
# detected integrase, which is what finds CALIN elements - cassette arrays whose integrase
# has been lost. Those are common on plasmids and are precisely the arrays whose cassettes
# nobody has characterised.
# `--pdf-off` is not an IntegronFinder option - the real flag is `--pdf`, and it is
# opt-IN, so it is simply omitted. Passing the invented flag made argparse exit 2 before
# any work was done, and check=False swallowed it: the stage wrote a header-only TSV and
# the rule reported success.
#
# `--circ` sets circular topology, which matters because 94% of these plasmids are closed
# and an integron spanning the origin is invisible under linear topology.
#
# `--keep-tmp` is dropped: it would create one directory per replicon, 143,503 of them.
#
# check=True: a tool failure must stop the run rather than produce an empty table.
def run(chunk):
    subprocess.run(
        f"integron_finder --local-max --circ --cpu 1 "
        f"--outdir {outdir / chunk.stem} {chunk}",
        shell=True, check=True, stdout=subprocess.DEVNULL)


with concurrent.futures.ThreadPoolExecutor(snakemake.threads) as pool:
    list(pool.map(run, chunks))

# Column order verified against integron_finder/results.py, not guessed. The earlier
# version read f[8] as the integron type; f[8] is `annotation` and the type is f[10].
#
#   0 ID_integron   1 ID_replicon   2 element        3 pos_beg    4 pos_end
#   5 strand        6 evalue        7 type_elt       8 annotation 9 model
#  10 type         11 default      12 distance_2attC 13 considered_topology
COLUMNS = ["ID_integron", "ID_replicon", "element", "pos_beg", "pos_end", "strand",
           "evalue", "type_elt", "annotation", "model", "type", "default",
           "distance_2attC", "considered_topology"]

rows = []
for path in outdir.rglob("*.integrons"):
    with open(path) as fh:
        for line in fh:
            # A replicon with no integron yields a file of comments only, which is a
            # legitimate result rather than a failure.
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if f[0] == "ID_integron":
                continue
            if len(f) != len(COLUMNS):
                raise SystemExit(
                    f"{path}: expected {len(COLUMNS)} columns, got {len(f)}. The "
                    "IntegronFinder output contract has changed; re-verify against "
                    "integron_finder/results.py before trusting any of these rows.")
            rec = dict(zip(COLUMNS, f))
            rows.append({"plasmid_id": rec["ID_replicon"],
                         "integron_id": rec["ID_integron"],
                         "element": rec["element"],
                         "start": rec["pos_beg"], "end": rec["pos_end"],
                         "integron_type": rec["type"],
                         "annotation": rec["annotation"],
                         "type_elt": rec["type_elt"]})

with open(snakemake.output[0], "w", newline="") as out:
    w = csv.DictWriter(out, fieldnames=["plasmid_id", "integron_id", "element",
                                        "start", "end", "integron_type", "annotation",
                                        "type_elt"],
                       delimiter="\t")
    w.writeheader()
    w.writerows(rows)

print(f"integron elements={len(rows)} from {len(chunks)} chunks")
shutil.rmtree(outdir / "chunks")
