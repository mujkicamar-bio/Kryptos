"""S8a phase 1: which proteins look like defence components.

Runs DefenseFinder's HMM profiles over the dereplicated protein set with
`--db-type unordered`, MacSyFinder's mode that reports components without calling systems.
Dereplication is valid here because the question is per protein and depends on sequence
identity alone; phase 2 asks about gene adjacency and runs on ordered ORFs instead.

Thresholds are DefenseFinder's own defaults (Tesson et al. 2022, Nat Commun 13:2561),
not overridden.

MacSyFinder writes no best_solution.tsv in unordered mode, and defense-finder's
post-treatment step opens that file unconditionally and exits 1 after a complete search.
The exit code therefore cannot decide success: the search succeeded when MacSyFinder wrote
the all_systems.tsv table of each model family (DefenseFinder, RM, CasFinder), which are
what this stage reads, and failed when any is missing.

Output: defence_components.tsv, one row per component hit, status SUCCESS; or one row with
status NOT_RUN and no seq_id when the models are absent and defence.required is false.
"""
import csv
import pathlib
import shutil
import subprocess
import sys

import _ctx  # noqa: F401

from darkorf import status
from plasmidann import defence

COLUMNS = ["seq_id", "component", "model", "hit_evalue", "status"]
outdir = pathlib.Path(snakemake.output.tsv).parent / "phase1"
outdir.mkdir(parents=True, exist_ok=True)


def write_table(rows):
    with open(snakemake.output.tsv, "w", newline="") as out:
        w = csv.DictWriter(out, fieldnames=COLUMNS, delimiter="\t")
        w.writeheader()
        w.writerows(rows)


# defense-finder also exits non-zero when its models are absent, so they are checked
# before the search. Absent and optional, the stage records NOT_RUN, which downstream
# stages read as "not searched" rather than "searched, nothing found".
models_dir = pathlib.Path(snakemake.params.models_dir)
if not (models_dir.is_dir() and any(models_dir.iterdir())):
    message = (f"DefenseFinder models are not installed at {models_dir}. Install them with "
               "`defense-finder update --models-dir <dir>` and point "
               "references.macsyfinder_models at it.")
    if snakemake.params.get("required", False):
        sys.exit(f"S8a: {message}")
    print(f"S8a: {message}\n     defence.required is false, so this stage records NOT_RUN.")
    write_table([{"status": status.NOT_RUN}])
    sys.exit(0)

# defense-finder runs MacSyFinder once per model family, in turn, and stops at the first
# failure; each run writes its all_systems.tsv under the preserved raw output. Tables left
# by an earlier run are removed first, so that only this run's tables are read.
raw = outdir / "defense-finder-tmp"
shutil.rmtree(raw, ignore_errors=True)
completed = subprocess.run(
    f"defense-finder run --db-type unordered --out-dir {outdir} "
    f"--workers {snakemake.threads} --preserve-raw {snakemake.input.faa}",
    shell=True)

system_tables = [raw / family / "all_systems.tsv" for family in ("DefenseFinder", "RM", "Cas")]
missing = [t.parent.name for t in system_tables if not t.is_file()]
if missing:
    sys.exit(
        f"S8a: defense-finder exited {completed.returncode} and wrote no all_systems.tsv "
        f"for the {', '.join(missing)} models under {raw}: the search itself failed.")

if completed.returncode != 0:
    print(f"S8a: defense-finder exited {completed.returncode} after MacSyFinder wrote "
          "every all_systems.tsv table: the search completed, and the "
          "failure is the post-treatment step that opens a best_solution.tsv "
          "--db-type unordered never writes.")

rows = [{**r, "status": status.SUCCESS}
        for r in defence.parse_all_systems(system_tables)]
write_table(rows)

print(f"phase 1: {len(rows)} component hits on {len({r['seq_id'] for r in rows})} "
      "unique proteins")
if not rows:
    # A collection without defence components is possible, but it empties every defence
    # column downstream, so it is reported.
    print("WARNING: no defence components found at all - check the models are installed")
