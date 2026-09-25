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

WHY THE TOOL'S NON-ZERO EXIT IS NOT AUTOMATICALLY A FAILURE

`--db-type unordered` is the mode this phase needs, and MacSyFinder does not write
best_solution.tsv in it. defense-finder's post-treatment step opens that file
unconditionally and raises FileNotFoundError, so a search that has just found systems in
every model family exits 1.

The search itself is complete at that point: MacSyFinder has already written all_systems.tsv
per model family, and that is what this stage reads. So the exit code alone cannot decide
whether the stage succeeded - the OUTPUT decides. No all_systems.tsv anywhere means the
search really did fail, and that is still fatal.
"""
import csv
import pathlib
import subprocess
import sys

import _ctx  # noqa: F401

from darkorf import status
from plasmidann import defence

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
        w = csv.DictWriter(out, fieldnames=["seq_id", "component", "model", "hit_evalue",
                                            "status"], delimiter="\t")
        w.writeheader()
    sys.exit(0)

# check=False, then judged on the OUTPUT - see the note above. Three stages in v2 used
# check=False and wrote well-formed EMPTY tables while reporting success, so the exit code
# is still recorded and an absent all_systems.tsv is still fatal.
completed = subprocess.run(
    f"defense-finder run --db-type unordered --out-dir {outdir} "
    f"--workers {snakemake.threads} --preserve-raw {snakemake.input.faa}",
    shell=True)

# MacSyFinder writes one all_systems.tsv per model family under the preserved raw output.
system_tables = sorted(outdir.rglob("all_systems.tsv"))
if not system_tables:
    sys.exit(
        f"S8a: defense-finder exited {completed.returncode} and wrote no all_systems.tsv "
        f"under {outdir}. The search itself failed, so there is nothing to read - this is "
        "not the known post-treatment crash, which leaves those tables in place.")

if completed.returncode != 0:
    print(f"S8a: defense-finder exited {completed.returncode}. MacSyFinder wrote "
          f"{len(system_tables)} all_systems.tsv tables, so the SEARCH completed; the "
          "failure is defense-finder's post-treatment step, which opens a best_solution.tsv "
          "that --db-type unordered never produces. The component hits below come from "
          "MacSyFinder's own output and are unaffected.")

rows = [{**r, "status": status.SUCCESS}
        for r in defence.parse_all_systems(system_tables)]

with open(snakemake.output.tsv, "w", newline="") as out:
    w = csv.DictWriter(out, fieldnames=["seq_id", "component", "model", "hit_evalue",
                                        "status"], delimiter="\t")
    w.writeheader()
    w.writerows(rows)

print(f"phase 1: {len(rows)} component hits on {len({r['seq_id'] for r in rows})} "
      "unique proteins")
if not rows:
    # Not an assertion: a collection genuinely containing no defence components is
    # possible. But it is worth shouting about, because an empty phase 1 empties the
    # defence stratum entirely and that should never pass unnoticed.
    print("WARNING: no defence components found at all - check the models are installed")
