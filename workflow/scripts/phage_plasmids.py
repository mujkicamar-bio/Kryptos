"""Rule phage_plasmids: phage-plasmids, from geNomad (Camargo et al. 2024, Nat Biotechnol 42:1303).

geNomad end-to-end runs as is, with its default filters, on the analysis-set FASTA. Every
plasmid gets one row, labelled only from geNomad's own output:

  genomad_virus         the sequence names geNomad's virus summary lists for the plasmid:
                        the whole plasmid, or '<plasmid>|provirus_<start>_<end>' for a
                        provirus inside it; empty when none
  n_virus_hallmarks     genes geNomad marks virus_hallmark = 1 on the plasmid (genes.tsv)
  phage_plasmid         1 when genomad_virus is not empty or n_virus_hallmarks >= 1, else 0
"""
import collections
import csv
import os
import pathlib
import shutil
import subprocess

import _ctx  # noqa: F401

from plasmidann.fasta import iter_fasta

fasta = pathlib.Path(snakemake.input.fasta)
outdir = pathlib.Path(snakemake.output.tsv).parent / "genomad"
# A rerun starts clean: result files from an interrupted run would be read below.
shutil.rmtree(outdir, ignore_errors=True)
outdir.mkdir(parents=True)
# geNomad calls MMseqs2 and ARAGORN, which its own environment provides.
exe = pathlib.Path(snakemake.params.exe)
env = dict(os.environ, PATH=f"{exe.parent.resolve()}{os.pathsep}{os.environ['PATH']}")
subprocess.run([str(exe), "end-to-end", "--cleanup", "--threads", str(snakemake.threads),
                str(fasta), str(outdir), snakemake.params.db], env=env, check=True)

prefix = fasta.name.split(".")[0]
virus = collections.defaultdict(list)
with open(outdir / f"{prefix}_summary" / f"{prefix}_virus_summary.tsv", newline="") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        virus[r["seq_name"].split("|")[0]].append(r["seq_name"])

# A gene is named '<plasmid>_<n>'.
hallmarks = collections.Counter()
with open(outdir / f"{prefix}_annotate" / f"{prefix}_genes.tsv", newline="") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        if r["virus_hallmark"] == "1":
            hallmarks[r["gene"].rsplit("_", 1)[0]] += 1

n = n_phage = 0
with open(snakemake.output.tsv, "w", newline="") as out:
    w = csv.writer(out, delimiter="\t")
    w.writerow(["plasmid_id", "genomad_virus", "n_virus_hallmarks", "phage_plasmid"])
    for pid, _ in iter_fasta(fasta):
        phage = int(bool(virus[pid]) or hallmarks[pid] >= 1)
        w.writerow([pid, ",".join(virus[pid]), hallmarks[pid], phage])
        n += 1
        n_phage += phage

print(f"phage-plasmids: {n_phage} of {n} plasmids ({sum(map(bool, virus.values()))} with a "
      f"geNomad virus call, {len(hallmarks)} with a virus hallmark gene)")
