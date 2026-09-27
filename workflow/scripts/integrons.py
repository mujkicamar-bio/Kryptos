"""S8b: integron elements per plasmid, with IntegronFinder (--local-max).

Input: the analysis-set FASTA, split into one chunk per core, and each plasmid's registry
topology (darkorf.circular.is_circular), passed as --topology-file. Output: integrons.tsv, one
row per IntegronFinder element (integrase, attC site, promoter, attI or cassette protein)
with its integron's type (complete, In0 or CALIN).
"""
import concurrent.futures
import csv
import pathlib
import shutil
import subprocess

import _ctx  # noqa: F401

from darkorf.circular import is_circular
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

# One "<replicon> <topology>" per line (integron_finder.topology). An integron spanning the
# origin of a circular plasmid is found only under circular topology.
topology_file = outdir / "topology.txt"
with open(snakemake.input.master, newline="") as fh, open(topology_file, "w") as out:
    for r in csv.DictReader(fh, delimiter="\t"):
        out.write(f"{r['plasmid_id']} {'circ' if is_circular(r.get('topology')) else 'lin'}\n")


# --local-max is the sensitive mode: it searches for attC sites beyond those adjacent to a
# detected integrase, which finds CALIN elements - cassette arrays whose integrase has been
# lost, common on plasmids. PDF output (--pdf) and --keep-tmp (one directory per replicon)
# are not requested. check=True: a tool failure must stop the run rather than produce an
# empty table.
def run(chunk):
    subprocess.run(
        f"integron_finder --local-max --topology-file {topology_file} --cpu 1 "
        f"--outdir {outdir / chunk.stem} {chunk}",
        shell=True, check=True, stdout=subprocess.DEVNULL)


with concurrent.futures.ThreadPoolExecutor(snakemake.threads) as pool:
    list(pool.map(run, chunks))

# Column order from integron_finder/results.py:
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
