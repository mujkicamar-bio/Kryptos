"""S8f: conjugation and mobilisation systems (CONJScan) on every plasmid.

CONJScan 2.1.0 uses model grammar 2.1, which needs MacSyFinder >= 2.1.6; DefenseFinder pins
MacSyFinder 2.1.4 and imports a package 2.1.6 no longer ships. CONJScan therefore runs from
its own environment, named by path in config (`conjugation.exe`), as pharokka does.

Unlike defence phase 2, the input is NOT pruned: the gembase holds every ORF of every
plasmid in genomic order (the order_orfs / gembase_id helpers of the defence stage). On the
100-plasmid test set, 31 of the 55 plasmids with a CONJScan system carry no defence
component, so the defence gembase would miss them; and MacSyFinder counts intervening genes,
so every ORF of a plasmid must be present.

Models `CONJScan/Plasmids all`: the package recommends running all models of one set
together, and the Plasmids set is the one built for plasmids (Coluzzi et al. 2022). Topology
is circular, as defence_systems.py passes; a per-replicon topology file gave identical calls
on the test set (212 of 212 ORFs), because the T4SS and dCONJ plasmid models allow up to
500 intervening genes.

ONE DATABASE, NOT ONE CHUNK PER CORE. Defence phase 2 splits its input into one chunk of
whole replicons per core. That split is not neutral here: HMMER's independent e-value scales
with the number of target sequences in the database, and MacSyFinder selects hits by it
(--i-evalue-sel 0.001, its default). On the test set a T4SS_MOBP1 hit with the same score
(16.9 bits) had i-evalue 0.0014 against all 5,673 ORFs and 0.00014 against a 553-ORF chunk,
so 8 chunks called 214 ORFs in 56 systems where one database called 212 in 55: the calls
depended on the core count. One MacSyFinder process over the whole gembase, with --worker
parallelising the 125 profile searches, gives calls that do not depend on -c; on 113,460
ORFs it took 16.4 s with 16 workers against 87.7 s with one. The e-values still depend on
the size of the collection searched, as with any single MacSyFinder database: the same
113,460-ORF set called 211 of each copy's 212 ORFs.
"""
import collections
import csv
import os
import pathlib
import shutil
import subprocess
import sys

import _ctx  # noqa: F401

from plasmidann.conjscan import (CLASS_COLUMNS, COLUMNS, MODEL_SET, installed_version,
                                 plasmid_class, read_best_solution)
from plasmidann.defence import gembase_id, order_orfs

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


# The NOT_RUN contract of the defence stage (spec section 7.2): a missing optional database
# must not be fatal, and an empty table must not read as "no conjugation systems here". The
# class table stays empty too - pMOBless for every plasmid would be exactly that claim.
version = installed_version(models_dir)
missing = []
if version is None:
    missing.append(f"CONJScan models are not installed at {models_dir}")
if shutil.which(exe) is None:
    missing.append(f"the MacSyFinder executable {exe} is absent or not executable")
if missing:
    if snakemake.params.required:
        sys.exit("S8f: " + "; ".join(missing) + ".")
    print("S8f: " + "; ".join(missing) + ". Recording NOT_RUN: conjugation system calls "
          "are absent-because-not-searched, not absent-because-searched.")
    write_tables([], [])
    sys.exit(0)

# The version is what the methods cite. 2.1.0 added the MOBM relaxase profile and raised the
# MOBH threshold, so a different release changes the calls under the old citation.
if version != str(snakemake.params.version):
    sys.exit(f"S8f: CONJScan {version} is installed at {models_dir}, but config "
             f"conjugation.version is {snakemake.params.version}.")

by_plasmid = collections.defaultdict(list)
with open(snakemake.input.index, newline="") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        by_plasmid[r["plasmid_id"]].append(r)
replicons = {p.replace("_", "-") for p in by_plasmid}
assert len(replicons) == len(by_plasmid), \
    "two plasmid ids map to the same gembase replicon after '_' -> '-'"

# A rerun starts clean: results from an interrupted run would be read below.
shutil.rmtree(outdir, ignore_errors=True)
outdir.mkdir(parents=True)
gembase = outdir / "gembase.faa"
back = {}
with open(gembase, "w") as out:
    for plasmid in sorted(by_plasmid):
        for position, orf in enumerate(order_orfs(by_plasmid[plasmid]), start=1):
            gid = gembase_id(plasmid, position)
            out.write(f">{gid}\n{orf['seq']}\n")
            back[gid] = (orf["orf_id"], plasmid)

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
    orf_id, plasmid = back[rec["hit_id"]]
    system = rec["model_fqn"].rsplit("/", 1)[-1]
    types_of[plasmid].add(system)
    # hit_status (mandatory / accessory / neutral) is the model's own declaration of the
    # component's role, and sys_wholeness how much of the model was found.
    rows.append({"orf_id": orf_id, "plasmid_id": plasmid, "system": system,
                 "system_id": rec["sys_id"], "component": rec["gene_name"],
                 "hit_status": rec["hit_status"], "sys_wholeness": rec["sys_wholeness"],
                 "conjscan_version": version})

# One class per plasmid in the ORF index; a plasmid with no ORF has no row.
classes = [(p, plasmid_class(types_of.get(p, ()))) for p in sorted(by_plasmid)]
write_tables(rows, classes)

print(f"S8f: CONJScan {version}: {len({r['orf_id'] for r in rows})} ORFs in "
      f"{len({r['system_id'] for r in rows})} systems across {len(types_of)} of "
      f"{len(by_plasmid)} plasmids; classes "
      f"{dict(collections.Counter(c for _, c in classes))}")
