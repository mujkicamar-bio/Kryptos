"""S8a phase 2: call defence systems from gene adjacency.

MacSyFinder is driven directly rather than through `defense-finder run`, for one reason:
the wrapper does not pass `--replicon-topology` through, and 94% of these plasmids are
circular. Under linear topology a system spanning the origin is invisible, and origin-
spanning genes are exactly what S1 worked to reconstruct.

`--db-type gembase` lets a single run hold every candidate replicon and still treat each
separately, which turns tens of thousands of per-plasmid invocations into one job. That one
job is split again into one chunk of whole replicons per core, each run as its own
MacSyFinder process with one worker: MacSyFinder's --worker parallelises only the profile
searches, and the test run kept 0.6 of 16 cores busy (2.2 CPU-s in 99 s), which scales to
~35-40 h for the production candidates. A system never spans two replicons, so splitting
between replicons changes no call.

Thresholds are MacSyFinder's and DefenseFinder's own published defaults (Tesson et al.
2022; Abby et al. 2014 for MacSyFinder). The quorum and co-localisation rules come from the
711 shipped model definitions - referenced by construction, since they ARE the published
models rather than our reinterpretation of them.
"""
import concurrent.futures
import csv
import pathlib
import shutil
import subprocess
import sys

import _ctx  # noqa: F401

from darkorf import status

outdir = pathlib.Path(snakemake.output.tsv).parent / "phase2"

# The same NOT_RUN contract as phase 1. Phase 2 cannot call a system from models that are
# not installed, and spec section 7.2 separates "not run" from "ran and found nothing": an
# empty systems table with no status would read as "this collection carries no defence
# systems", which is a biological claim this run has not earned.
models_dir = pathlib.Path(snakemake.params.models_dir)
if not (models_dir.is_dir() and any(models_dir.iterdir())):
    if snakemake.params.get("required", False):
        sys.exit(f"S8a phase 2: DefenseFinder models are not installed at {models_dir}.")
    print(f"S8a phase 2: models absent at {models_dir}; recording NOT_RUN. Defence system "
          "calls are absent-because-not-searched, not absent-because-searched.")
    with open(snakemake.output.tsv, "w", newline="") as out:
        csv.DictWriter(out, fieldnames=["orf_id", "plasmid_id", "gembase_id", "system",
                                        "system_id", "component", "hit_evalue",
                                        "hit_status", "sys_wholeness", "hit_gene_ref",
                                        "hit_profile_cov", "status"],
                       delimiter="\t").writeheader()
    sys.exit(0)

# skip_run exists for the parser test, which pre-populates the output tree. It is never
# set by the workflow: a missing MacSyFinder run in production must fail, not be skipped.
if not snakemake.params.get("skip_run", False):
    # A rerun starts clean: results from an interrupted run would be read below.
    shutil.rmtree(outdir, ignore_errors=True)
    chunk_dir = outdir / "chunks"
    chunk_dir.mkdir(parents=True)
    # Whole replicons to the chunk with the fewest genes so far. defence_gembase writes each
    # replicon's genes together, and the replicon is the gembase id up to the last '_'.
    n = snakemake.threads
    handles = [open(chunk_dir / f"chunk_{i:03d}.faa", "w") for i in range(n)]
    genes = [0] * n
    current, target = None, 0
    with open(snakemake.input.faa) as fh:
        for line in fh:
            if line.startswith(">"):
                replicon = line[1:].split()[0].rsplit("_", 1)[0]
                if replicon != current:
                    current, target = replicon, genes.index(min(genes))
                genes[target] += 1
            handles[target].write(line)
    for h in handles:
        h.close()
    chunks = [chunk_dir / f"chunk_{i:03d}.faa" for i in range(n) if genes[i]]

    def run(chunk):
        subprocess.run(
            f"macsyfinder --models-dir {snakemake.params.models_dir} "
            f"--models defense-finder-models all "
            f"--sequence-db {chunk} "
            f"--db-type gembase --replicon-topology circular "
            f"--worker 1 --out-dir {outdir / chunk.stem} --mute",
            shell=True, check=True)

    with concurrent.futures.ThreadPoolExecutor(n) as pool:
        list(pool.map(run, chunks))
    shutil.rmtree(chunk_dir)

# Map gembase ids back to our orf_ids.
back = {}
with open(snakemake.input.map, newline="") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        back[r["gembase_id"]] = r

rows = []
for path in outdir.rglob("best_solution.tsv"):
    with open(path) as fh:
        header = None
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if header is None:
                header = fields
                continue
            rec = dict(zip(header, fields))
            gid = rec.get("hit_id", "")
            src = back.get(gid)
            if not src:
                continue
            rows.append({"orf_id": src["orf_id"], "plasmid_id": src["plasmid_id"],
                         "gembase_id": gid, "system": rec.get("model_fqn", ""),
                         "system_id": rec.get("sys_id", ""),
                         "component": rec.get("gene_name", ""),
                         "hit_evalue": rec.get("hit_i_eval", ""),
                         # A mandatory component of a complete system and a neutral
                         # component of a fragment are different evidence and arrived as
                         # identical rows. hit_status is the model's own declaration of
                         # which it is, and sys_wholeness is how much of the model was
                         # found - both are MacSyFinder's numbers, not our
                         # reinterpretation of them.
                         "hit_status": rec.get("hit_status", ""),
                         "sys_wholeness": rec.get("sys_wholeness", ""),
                         "hit_gene_ref": rec.get("hit_gene_ref", ""),
                         "hit_profile_cov": rec.get("hit_profile_cov", ""),
                         "status": status.SUCCESS})

with open(snakemake.output.tsv, "w", newline="") as out:
    w = csv.DictWriter(out, fieldnames=["orf_id", "plasmid_id", "gembase_id", "system",
                                        "system_id", "component", "hit_evalue",
                                        "hit_status", "sys_wholeness", "hit_gene_ref",
                                        "hit_profile_cov", "status"],
                       delimiter="\t")
    w.writeheader()
    w.writerows(rows)

print(f"phase 2: {len(rows)} genes in {len({r['system_id'] for r in rows})} systems "
      f"across {len({r['plasmid_id'] for r in rows})} plasmids")
