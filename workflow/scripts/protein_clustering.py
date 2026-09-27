"""S2f: cluster every unique protein into families, before any annotation.

A family is a sequence cluster and nothing else, so it needs only the sequences and can be
made first. Made first, it serves two stages:

  * cascade_selection (S2s) reads the primary resolution to decide which proteins the
    cascade searches: the families that hold a small-plasmid protein;
  * protein_families (Stage 5) reads all three resolutions after the cascade and adds
    what the annotation says about each family.

Every unique protein is clustered - small and large plasmids, annotated and not - so no
member is lost to a filter applied before the clustering.

Three resolutions (close, intermediate, broad), thresholds in config/targets.yaml. MMseqs2
easy-cluster writes <prefix>_cluster.tsv (representative<TAB>member) and
<prefix>_rep_seq.fasta per resolution; those files are the outputs.
"""
import pathlib
import subprocess

import _ctx  # noqa: F401

from plasmidann import scratch

cfg = snakemake.params.clustering
out_dir = pathlib.Path(snakemake.output.clusters[0]).parent

if not cfg["resolutions"]:
    raise SystemExit("protein_clustering: no resolutions configured - "
                     "clustering.resolutions is empty.")

for resolution, thresholds in sorted(cfg["resolutions"].items()):
    tmp = scratch.scratch_dir(out_dir, f"mmseqs_tmp_{resolution}")
    prefix = out_dir / f"families_{resolution}"
    subprocess.run(
        f"mmseqs easy-cluster {snakemake.input.faa} {prefix} {tmp} "
        f"--min-seq-id {thresholds['min_seq_id']} -c {thresholds['coverage']} "
        f"--cov-mode {cfg['cov_mode']} --cluster-mode {cfg['cluster_mode']} "
        f"--cluster-reassign --threads {snakemake.threads} -v 1",
        shell=True, check=True)
    # MMseqs2 also writes every member's sequence again (<prefix>_all_seqs.fasta); nothing
    # reads it, and at 3.5M proteins it is a gigabyte per resolution.
    pathlib.Path(f"{prefix}_all_seqs.fasta").unlink(missing_ok=True)
    # Released per resolution: three resolutions hold their databases at once otherwise.
    scratch.release(tmp)
    n = sum(1 for _ in open(f"{prefix}_rep_seq.fasta") if _.startswith(">"))
    print(f"{resolution}: {n} clusters")
