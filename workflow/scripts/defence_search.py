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

outdir = pathlib.Path(snakemake.output.tsv).parent / "phase1"
outdir.mkdir(parents=True, exist_ok=True)

# check=True: a tool failure must stop the run. Three stages in v2 used check=False and
# wrote well-formed EMPTY tables while reporting success, which is how a 175-construct
# stratum can silently vanish.
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
                         "hit_evalue": r.get("hit_i_eval", "")})

with open(snakemake.output.tsv, "w", newline="") as out:
    w = csv.DictWriter(out, fieldnames=["seq_id", "component", "model_unverified",
                                        "hit_evalue"], delimiter="\t")
    w.writeheader()
    w.writerows(rows)

print(f"phase 1: {len(rows)} component hits on {len({r['seq_id'] for r in rows})} "
      "unique proteins")
if not rows:
    # Not an assertion: a collection genuinely containing no defence components is
    # possible. But it is worth shouting about, because an empty phase 1 empties the
    # defence stratum entirely and that should never pass unnoticed.
    print("WARNING: no defence components found at all - check the models are installed")
