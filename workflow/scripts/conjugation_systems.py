"""S8f: conjugation and mobilisation systems (CONJScan) on every plasmid.

CONJScan 2.1.0 uses model grammar 2.1, which needs MacSyFinder >= 2.1.6, while
DefenseFinder pins MacSyFinder 2.1.4. CONJScan therefore runs from its own environment,
named by path in config (`conjugation.exe`).

Input: the ORF index. Every ORF of every plasmid is written in genomic order as one gembase
database (plasmidann.defence.gembase_records), not only the plasmids with a defence
component: on the 100-plasmid test set 31 of the 55 plasmids with a CONJScan system carry no
defence component. Models `CONJScan/Plasmids all`: the package recommends running all models
of one set together, and the Plasmids set is the one built for plasmids (Coluzzi et al.
2022). Topology is circular; a per-replicon topology file gave identical calls on the test
set (212 of 212 ORFs), because the T4SS and dCONJ plasmid models allow up to 500
intervening genes.

One MacSyFinder process searches the whole database, with the cores given to --worker.
HMMER's independent e-value scales with the number of sequences in the database and
MacSyFinder selects hits by it (--i-evalue-sel 0.001, its default): on the test set a
16.9-bit T4SS_MOBP1 hit had i-evalue 0.0014 against all 5,673 ORFs and 0.00014 against a
553-ORF chunk, so per-core chunks would make the calls depend on the core count.

Outputs: conjugation_systems.tsv, one row per system component (status SUCCESS), and
conjugation_plasmid_class.tsv, one mobility class per plasmid. When the models or the
executable are absent and conjugation.required is false, the systems table holds one row
with status NOT_RUN and no orf_id, and the class table is empty: pMOBless for every
plasmid would be a claim the run has not made.
"""
import collections
import csv
import os
import pathlib
import shutil
import subprocess
import sys

import _ctx  # noqa: F401

from darkorf import status
from plasmidann.conjscan import (
    CLASS_COLUMNS,
    COLUMNS,
    MODEL_SET,
    installed_version,
    plasmid_class,
    read_best_solution,
)
from plasmidann.defence import gembase_records

models_dir = pathlib.Path(snakemake.params.models_dir)
exe = snakemake.params.exe
outdir = pathlib.Path(snakemake.output.systems).parent / "conjscan"


def write_tables(rows, classes):
    with open(snakemake.output.systems, "w", newline="") as out:
        w = csv.DictWriter(out, fieldnames=COLUMNS, delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    with open(snakemake.output.classes, "w", newline="") as out:
        w = csv.writer(out, delimiter="\t")
        w.writerow(CLASS_COLUMNS)
        w.writerows(classes)


version = installed_version(models_dir)
missing = []
if version is None:
    missing.append(f"CONJScan models are not installed at {models_dir}")
if shutil.which(exe) is None:
    missing.append(f"the MacSyFinder executable {exe} is absent or not executable")
if missing:
    if snakemake.params.required:
        sys.exit("S8f: " + "; ".join(missing) + ".")
    print("S8f: " + "; ".join(missing) + ". Recording NOT_RUN.")
    write_tables([{"status": status.NOT_RUN}], [])
    sys.exit(0)

# The version is what the methods cite. 2.1.0 added the MOBM relaxase profile and raised the
# MOBH threshold, so a different release gives different calls.
if version != str(snakemake.params.version):
    sys.exit(f"S8f: CONJScan {version} is installed at {models_dir}, but config "
             f"conjugation.version is {snakemake.params.version}.")

# A rerun starts clean: results from an interrupted run would be read below.
shutil.rmtree(outdir, ignore_errors=True)
outdir.mkdir(parents=True)
gembase = outdir / "gembase.faa"
# gembase_id -> orf_id, for mapping the calls back. The index is streamed one plasmid at a
# time, so the sequences are never all held in memory.
orf_of = {}
with open(snakemake.input.index, newline="") as fh, open(gembase, "w") as out:
    for gid, _, orf in gembase_records(csv.DictReader(fh, delimiter="\t")):
        out.write(f">{gid}\n{orf['seq']}\n")
        orf_of[gid] = orf["orf_id"]
plasmids = sorted({oid.rsplit("|", 1)[0] for oid in orf_of.values()})

# The executable's own directory goes first on PATH, so MacSyFinder finds the hmmsearch of
# its own environment (envs/conjscan pins hmmer) rather than whichever the rule's
# environment puts first; a venv without hmmer still falls back to the rule's PATH.
exe_dir = str(pathlib.Path(exe).absolute().parent)
env = dict(os.environ, PATH=os.pathsep.join([exe_dir, os.environ["PATH"]]))
subprocess.run(
    [exe, "--models-dir", str(models_dir), "--models", MODEL_SET, "all",
     "--sequence-db", str(gembase), "--db-type", "gembase",
     "--replicon-topology", "circular", "--worker", str(snakemake.threads),
     "--out-dir", str(outdir / "run"), "--mute"],
    check=True, env=env)
# The gembase is rebuilt from the ORF index on every run; only MacSyFinder's results stay.
for path in outdir.glob("gembase.faa*"):
    path.unlink()

rows = []
types_of = collections.defaultdict(set)
for rec in read_best_solution(outdir / "run" / "best_solution.tsv"):
    orf_id = orf_of[rec["hit_id"]]
    plasmid = orf_id.rsplit("|", 1)[0]
    system = rec["model_fqn"].rsplit("/", 1)[-1]
    types_of[plasmid].add(system)
    rows.append({"orf_id": orf_id, "plasmid_id": plasmid, "system": system,
                 "system_id": rec["sys_id"], "component": rec["gene_name"],
                 "hit_status": rec["hit_status"], "sys_wholeness": rec["sys_wholeness"],
                 "conjscan_version": version, "status": status.SUCCESS})

# One class per plasmid in the ORF index; a plasmid with no ORF has no row.
classes = [(p, plasmid_class(types_of.get(p, ()))) for p in plasmids]
write_tables(rows, classes)

print(f"S8f: CONJScan {version}: {len({r['orf_id'] for r in rows})} ORFs in "
      f"{len({r['system_id'] for r in rows})} systems across {len(types_of)} of "
      f"{len(plasmids)} plasmids; classes "
      f"{dict(collections.Counter(c for _, c in classes))}")
