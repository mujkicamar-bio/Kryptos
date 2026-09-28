"""Rule structure_search: structural homology for dark-family representatives (or every
dark protein), via Foldseek with ProstT5.

Sequence search has already failed on these proteins by definition. Structure is
conserved longer than sequence, so a family with no sequence homolog can still match a
characterised fold, which suggests a mechanism to test. Foldseek derives the 3Di
structural alphabet directly from sequence with its ProstT5 model, which avoids predicting
a structure for every family first. This is a screening step, not a structure
determination.

SCOPE. `structure.scope: representatives` (the default) searches one sequence per dark
family. ProstT5 is a transformer and its cost grows with the number of queries, which makes
this the difference between the largest job in the pipeline and a modest one. `scope: all`
searches every dark protein: a representative is chosen by MMseqs2 on sequence criteria,
and at 30% identity only some members of a family may reach a recognisable fold.

Output: structure_hits.tsv, the best match per query (lowest e-value) with the target's
description line (`theader`), since the accession alone (12as-assembly1_A) names nothing;
status SUCCESS. When structure.required is false and the Foldseek database or the ProstT5
model is absent, one row with status NOT_RUN and no seq_id.
"""
import csv
import pathlib
import subprocess
import sys

import _ctx  # noqa: F401

from darkorf import status
from plasmidann import scratch

COLUMNS = ["seq_id", "target", "target_description", "fident", "alnlen", "evalue", "bits",
           "status"]
cfg = snakemake.params.structure
target_db = snakemake.params.target_db


def write_table(rows):
    with open(snakemake.output[0], "w", newline="") as out:
        w = csv.DictWriter(out, fieldnames=COLUMNS, delimiter="\t")
        w.writeheader()
        w.writerows(rows)


# Pre-flight fails the run when structure.required is true and a reference is missing, so
# a missing reference here means structure is optional: the stage records NOT_RUN.
absent = [p for p in (target_db, snakemake.params.prostt5) if not pathlib.Path(p).exists()]
if absent and not cfg["required"]:
    print(f"structure_search: {', '.join(absent)} absent and structure.required is false; "
          "recording NOT_RUN.")
    write_table([{"status": status.NOT_RUN}])
    sys.exit(0)

# Foldseek's working files, the query set and Foldseek's raw table are all intermediates.
tmp = scratch.scratch_dir(pathlib.Path(snakemake.output[0]).parent, "foldseek_tmp")
raw = tmp / "structure_hits.raw.tsv"

# ---- the query set: one sequence per dark family, or every dark protein ----------------
scope = cfg.get("scope", "representatives")
if scope not in ("representatives", "all"):
    sys.exit(f"structure.scope must be 'representatives' or 'all', not {scope!r}")

query_faa = snakemake.input.faa
if scope == "representatives":
    wanted = set()
    with open(snakemake.input.families, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row.get("representative"):
                wanted.add(row["representative"])

    query_faa = str(tmp / "structure_query.faa")
    n_written = 0
    with open(query_faa, "w") as out:
        emit = False
        for line in open(snakemake.input.faa):
            if line[0] == ">":
                emit = line[1:].split()[0] in wanted
                n_written += emit
            if emit:
                out.write(line)

    # A representative named in the families table but absent from the dark FASTA means the
    # two disagree about what the dark set is, and every structural count would inherit it.
    if n_written != len(wanted):
        sys.exit(f"structure_search: {len(wanted)} family representatives declared but {n_written} found "
                 f"in {snakemake.input.faa} - the families table and the dark set disagree.")
    print(f"structure_search: scope=representatives, searching {n_written} of "
          f"{sum(1 for l in open(snakemake.input.faa) if l[0] == '>')} dark proteins")
else:
    print("structure_search: scope=all, searching every dark protein")

# GPU use is configured (structure.gpu), not detected, so that the settings recorded for a
# run say how it ran; ProstT5 inference is the cost of this stage, and foldseek takes --gpu
# for both the ProstT5 step and the search.
gpu_flag = " --gpu 1" if cfg.get("gpu", False) else ""

subprocess.run(
    f"foldseek easy-search {query_faa} {target_db} "
    f"{raw} {tmp} --prostt5-model {snakemake.params.prostt5} "
    f"-e {cfg['max_evalue']} --threads {snakemake.threads}{gpu_flag} "
    # `prob` is not requested: it is derived from C-alpha coordinates, which a ProstT5
    # query database does not carry, and requesting it makes foldseek exit 1. `theader`
    # is the target's description line.
    f"--format-output query,target,theader,fident,alnlen,evalue,bits "
    f"--max-seqs 10 -v 1",
    shell=True, check=True)

rows = []
if raw.exists():
    best = {}
    for line in open(raw):
        q, t, theader, fident, alnlen, ev, bits = line.rstrip("\n").split("\t")
        if q not in best or float(ev) < float(best[q]["evalue"]):
            best[q] = {"seq_id": q, "target": t, "target_description": theader,
                       "fident": fident, "alnlen": alnlen, "evalue": ev, "bits": bits,
                       "status": status.SUCCESS}
    rows = list(best.values())
else:
    # foldseek exited 0 without writing its output: a tool failure, not "no hits".
    raise SystemExit(f"foldseek reported success but wrote no output to {raw}.")

write_table(rows)

print(f"structural matches={len(rows)}")

# The scratch directory is removed only here, on the ordinary path. A script that raised
# never reaches this line, and its intermediates are what the failure is diagnosed from.
scratch.release(tmp)

