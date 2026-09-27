"""Search every unique protein against the plasmid-specific label databases.

TADB, BacMet, oriTDB, mobileOG-db, dbAPIS and Anti-CRISPRdb are searched with DIAMOND and
labelled at PlasAnn's identity and coverage tiers; CARD's protein homolog models with RGI's
Perfect / Strict rules; AMRFinderPlus runs as itself, with its own curated rules. The rules
and the reasons for them are in plasmidann.labeldb; this script runs the tools and writes
one long table, 08_protein_labels/protein_labels_plasmid.tsv, which protein_labels merges.

A database is searched over every unique protein, not only the dark ones: these labels are
the vocabulary of a dark protein's neighbours, and the neighbours are the named proteins.

NOT RUN IS NOT "NOTHING FOUND"

A required database that is absent stops the run at pre-flight. One that is not required
and absent is recorded as NOT_RUN in label_databases_status.tsv, as the defence stage does.
Without that file an absent Anti-CRISPRdb and an Anti-CRISPRdb that matched nothing would
both read as "no anti-CRISPR on any plasmid".
"""
import csv
import os
import pathlib
import shutil
import subprocess
import sys

import _ctx  # noqa: F401

from darkorf import status
from plasmidann import labeldb

cfg = snakemake.params.labels
amr = snakemake.params.amrfinder
faa = snakemake.input.faa
threads = snakemake.threads
refs = pathlib.Path(cfg["dir"])
workdir = pathlib.Path(snakemake.output.tsv).parent / "label_databases"

# ------------------------------------------------------------------------------------
# What is installed. A required database or AMRFinderPlus that is absent has already
# stopped the run at pre-flight (rule input preflight.tsv); the rest are NOT_RUN.
# ------------------------------------------------------------------------------------
absent = [db for db in labeldb.DATABASES if not (refs / db).is_dir()]
amr_db = pathlib.Path(amr["database"])
amr_absent = []
if shutil.which(amr["executable"]) is None:
    amr_absent.append(f"the amrfinder executable {amr['executable']} is absent or not "
                      "executable")
if not (amr_db / "version.txt").is_file():
    amr_absent.append(f"the AMRFinderPlus database {amr_db} (version.txt) is absent")


def version_of(db):
    """The one-line VERSION file. A database without one cannot be cited, so it halts."""
    path = refs / db / "VERSION"
    text = path.read_text().strip() if path.is_file() else ""
    if not text:
        sys.exit(f"label_databases: {path} is missing or empty; every label must carry "
                 "its release.")
    return text.splitlines()[0]


def diamond(db, entries, args):
    """DIAMOND blastp of all unique proteins against `entries`; yields hit dicts."""
    ref = workdir / f"{db}.faa"
    with open(ref, "w") as out:
        for e in entries:
            out.write(f">{e['key']}\n{e['sequence']}\n")
    subprocess.run(["diamond", "makedb", "--in", str(ref), "-d", str(workdir / db),
                    "--threads", str(threads), "--quiet"], check=True)
    raw = workdir / f"{db}.hits.tsv"
    # --tmpdir to node-local disk, as in tier_search: DIAMOND otherwise spills next to
    # the output, on the shared filesystem.
    subprocess.run(
        f"diamond blastp -q {faa} -d {workdir / db} -o {raw} {args} "
        f"--threads {threads} --tmpdir {snakemake.resources.tmpdir} --quiet "
        f"--outfmt 6 qseqid sseqid pident qcovhsp scovhsp bitscore sstart send slen",
        shell=True, check=True)
    # Streamed, not collected: at full scale the tiered searches return tens of millions
    # of rows, and only the best hit per protein is kept.
    with open(raw) as fh:
        for line in fh:
            q, s, pid, qc, sc, bits, ss, se, sl = line.rstrip("\n").split("\t")
            yield {"query": q, "key": s, "pident": float(pid), "qcov": float(qc),
                   "scov": float(sc), "bitscore": float(bits), "sstart": int(ss),
                   "send": int(se), "slen": int(sl)}


# A rerun starts clean: hits from an interrupted run would be read below.
shutil.rmtree(workdir, ignore_errors=True)
workdir.mkdir(parents=True)

rows, report = [], []
for db in labeldb.DATABASES:
    if db in absent:
        print(f"label_databases: {db} not installed under {refs}; recording NOT_RUN. Its "
              "labels are absent-because-not-searched, not absent-because-searched.")
        report.append([db, status.NOT_RUN, "", ""])
        continue
    version = version_of(db)
    if db == "card":
        models = labeldb.card_models(refs / db / "card.json")
        best = labeldb.best_card_hits(
            diamond(db, models.values(), labeldb.DIAMOND_CARD_ARGS), models)
        entries = models
    else:
        listed = labeldb.load_reference(db, refs / db)
        best = labeldb.best_tiered_hits(diamond(db, listed, labeldb.DIAMOND_TIERED_ARGS))
        entries = {e["key"]: e for e in listed}
    found = labeldb.label_rows(db, best, entries, version)
    rows += found
    n = len({r["seq_id"] for r in found})
    report.append([db, status.SUCCESS if n else status.NO_HIT, version, n])
    print(f"label_databases: {db} {version}: {len(entries)} reference entries, "
          f"{n} proteins labelled")

# ------------------------------------------------------------------------------------
# AMRFinderPlus, protein mode with --plus: AMR, stress (metal, biocide, acid, heat) and
# virulence elements, by NCBI's curated rules (Feldgarden et al. 2021, Sci. Rep.).
# ------------------------------------------------------------------------------------
if amr_absent:
    print("label_databases: " + "; ".join(amr_absent) + ". Recording AMRFinderPlus as NOT_RUN.")
    report.append(["amrfinder", status.NOT_RUN, "", ""])
else:
    version = (amr_db / "version.txt").read_text().strip().splitlines()[0]
    table = workdir / "amrfinder.tsv"
    # The executable's own directory first on PATH, so amrfinder finds the blastp and
    # hmmsearch of its own environment rather than those of the rule's environment.
    exe_dir = str(pathlib.Path(amr["executable"]).absolute().parent)
    env = dict(os.environ, PATH=os.pathsep.join([exe_dir, os.environ["PATH"]]))
    subprocess.run([amr["executable"], "-p", faa, "--plus", "--database", str(amr_db),
                    "--threads", str(threads), "-o", str(table)], check=True, env=env)
    found = labeldb.parse_amrfinder(table, version)
    rows += found
    n = len({r["seq_id"] for r in found})
    report.append(["amrfinder", status.SUCCESS if n else status.NO_HIT, version, n])
    print(f"label_databases: AMRFinderPlus database {version}: {n} proteins labelled")

with open(snakemake.output.tsv, "w", newline="") as out:
    w = csv.DictWriter(out, fieldnames=labeldb.COLUMNS, delimiter="\t")
    w.writeheader()
    w.writerows(sorted(rows, key=lambda r: (r["seq_id"], r["source"], r["label"],
                                            r["sub_label"])))
with open(snakemake.output.status, "w", newline="") as out:
    w = csv.writer(out, delimiter="\t")
    w.writerow(["database", "status", "version", "n_proteins"])
    w.writerows(report)
# The raw DIAMOND tables, the reference FASTA files and their databases are parsed into the
# tables above; at full scale the raw tables hold tens of millions of rows.
shutil.rmtree(workdir)

print(f"label_databases: {len(rows)} label rows on {len({r['seq_id'] for r in rows})} proteins")
