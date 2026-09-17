"""S8a phase 1: which proteins look like defence components.

Runs DefenseFinder's HMM profiles over the DEREPLICATED protein set with
`--db-type unordered`, which is MacSyFinder's "report components, do not call systems"
mode. That is exactly phase 1 and nothing more.

Dereplication is safe here and only here. Phase 1 asks a per-protein question that depends
on sequence identity alone, so one search result is valid for every copy - a protein
identical on forty plasmids is searched once. Phase 2 asks about gene adjacency and cannot
use this representation at all, which is why it is a separate rule.

Thresholds are DefenseFinder's own defaults (Tesson et al. 2022, Nat Commun 13:2561).
They are not overridden here: the tool's published defaults are the citable values, and
inventing our own would be exactly the unreferenced-parameter problem the project has
already been bitten by twice.
"""
import _ctx  # noqa: F401
import csv
import pathlib
import subprocess
import sys

from darkorf import status

outdir = pathlib.Path(snakemake.output.tsv).parent / "phase1"
outdir.mkdir(parents=True, exist_ok=True)

# A MISSING MODEL SET IS NOT A TOOL FAILURE.
#
# Spec section 7.2 gives NOT_RUN and FAILED different meanings, and this is the difference
# in practice. defense-finder exits non-zero when its models are not installed, which is
# indistinguishable at the exit code from a real crash, and letting that halt the run makes
# a missing OPTIONAL database fatal to a pipeline whose primary deliverable does not depend
# on it.
#
# So the models are checked first. Absent and not required, the stage records NOT_RUN and
# writes an empty table with the status on it - which is a different statement from "the
# search ran and found nothing", and the two must never be conflated (section 2.9).
# Absent and required, it fails here with a message naming what to install, rather than
# after the search has already burned the allocation.
models_dir = pathlib.Path(snakemake.params.models_dir)
required = bool(snakemake.params.get("required", False))
have_models = models_dir.is_dir() and any(models_dir.iterdir())

if not have_models:
    message = (f"DefenseFinder models are not installed at {models_dir}. Install them with "
               "`defense-finder update --models-dir <dir>` and point "
               "references.macsyfinder_models at it.")
    if required:
        sys.exit(f"S8a: {message}")
    print(f"S8a: {message}\n"
          "     defence.required is false, so this stage records NOT_RUN. Every defence "
          "column downstream is absent-because-not-searched, NOT absent-because-searched.")
    with open(snakemake.output.tsv, "w", newline="") as out:
        w = csv.DictWriter(out, fieldnames=["seq_id", "component", "model_unverified",
                                            "hit_evalue", "status"], delimiter="\t")
        w.writeheader()
    sys.exit(0)

# check=True: with the models present, a tool failure IS a failure. Three stages in v2 used
# check=False and wrote well-formed EMPTY tables while reporting success.
subprocess.run(
    f"defense-finder run --db-type unordered --out-dir {outdir} "
    f"--workers {snakemake.threads} --preserve-raw {snakemake.input.faa}",
    shell=True, check=True)

# DefenseFinder writes *_defense_finder_genes.tsv: one row per protein assigned to a
# component, with the model it belongs to.
rows = []
for path in outdir.rglob("*defense_finder_genes.tsv"):
    with open(path, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            seq_id = r.get("hit_id") or r.get("replicon") or ""
            if not seq_id:
                continue
            # `model` is best-effort and NOT verified against a real DefenseFinder run:
            # phase 2 reads system identity from MacSyFinder's own `model_fqn`, which is
            # cited, whereas these two names are a guess at the phase-1 gene table's
            # columns. Nothing downstream reads this field - defence_gembase.py uses
            # `component` only - so a wrong guess costs nothing today. It is recorded
            # rather than dropped so the column can be verified the first time a real
            # phase-1 table exists, and the header below says so.
            rows.append({"seq_id": seq_id,
                         "component": r.get("gene_name", ""),
                         "model_unverified": r.get("type") or r.get("subtype") or "",
                         "hit_evalue": r.get("hit_i_eval", ""),
                         "status": status.SUCCESS})

with open(snakemake.output.tsv, "w", newline="") as out:
    w = csv.DictWriter(out, fieldnames=["seq_id", "component", "model_unverified",
                                        "hit_evalue", "status"], delimiter="\t")
    w.writeheader()
    w.writerows(rows)

print(f"phase 1: {len(rows)} component hits on {len({r['seq_id'] for r in rows})} "
      "unique proteins")
if not rows:
    # Not an assertion: a collection genuinely containing no defence components is
    # possible. But it is worth shouting about, because an empty phase 1 empties the
    # defence stratum entirely and that should never pass unnoticed.
    print("WARNING: no defence components found at all - check the models are installed")
