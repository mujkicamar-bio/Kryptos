"""Rule defence_systems, DefenseFinder phase 2: call defence systems from gene adjacency.

MacSyFinder is run directly rather than through `defense-finder run`, because the wrapper
passes no replicon topology: under linear topology a system spanning the origin of a
circular plasmid is not found, and under circular topology the genes at the two ends of a
linear plasmid can be joined into one system. Each replicon gets its registry topology
(darkorf.circular.is_circular) through `--topology-file`. `--db-type gembase` holds every
candidate replicon in one database and still treats each separately.

The model families and their options are those of `defense-finder run` (DefenseFinder
3.0.0, Tesson et al. 2022): the 486 DefenseFinder and 24 RM models with --coverage-profile
0.4 and --exchangeable-weight 1, and the 44 CasFinder models without options, so that the
CasFinder package's own configuration applies. The AntiDefenseFinder models, which
DefenseFinder searches only on request, are not searched. Other thresholds and the quorum
and co-localisation rules are those of MacSyFinder (Abby et al. 2014) and of each model.

Each family is one MacSyFinder process over all candidate replicons, with the cores given
to --worker. HMMER's independent e-value is the p-value times the number of sequences in the
database, and MacSyFinder keeps a hit only below --i-evalue-sel (0.001, its default), so
splitting the database into per-core chunks would make the calls depend on the core count.
The e-values still depend on the size of the candidate set, as for any single database.

Output: defence_systems.tsv, one row per system component with MacSyFinder's hit_status
(mandatory, accessory, neutral) and sys_wholeness, status SUCCESS; or one row with status
NOT_RUN and no orf_id when the models are absent and defence.required is false.
"""
import csv
import pathlib
import shutil
import subprocess
import sys

import _ctx  # noqa: F401

from darkorf import status
from darkorf.circular import is_circular

COLUMNS = ["orf_id", "plasmid_id", "gembase_id", "system", "system_id", "component",
           "hit_evalue", "hit_status", "sys_wholeness", "hit_gene_ref", "hit_profile_cov",
           "status"]
outdir = pathlib.Path(snakemake.output.tsv).parent / "phase2"
# Output directory -> MacSyFinder options, as `defense-finder run` passes them.
FAMILIES = {
    "DefenseFinder": ["--models", "defense-finder-models/DefenseFinder", "all",
                      "--coverage-profile", "0.4", "--exchangeable-weight", "1"],
    "RM": ["--models", "defense-finder-models/RM", "all",
           "--coverage-profile", "0.4", "--exchangeable-weight", "1"],
    "Cas": ["--models", "CasFinder", "all"],
}


def write_table(rows):
    with open(snakemake.output.tsv, "w", newline="") as out:
        w = csv.DictWriter(out, fieldnames=COLUMNS, delimiter="\t")
        w.writeheader()
        w.writerows(rows)


models_dir = pathlib.Path(snakemake.params.models_dir)
if not (models_dir.is_dir() and any(models_dir.iterdir())):
    if snakemake.params.get("required", False):
        sys.exit(f"defence_systems: DefenseFinder models are not installed at {models_dir}.")
    print(f"defence_systems: models absent at {models_dir}; recording NOT_RUN.")
    write_table([{"status": status.NOT_RUN}])
    sys.exit(0)

# skip_run is set only by the parser test, which pre-populates the output tree.
if not snakemake.params.get("skip_run", False):
    # A rerun starts clean: results from an interrupted run would be read below.
    shutil.rmtree(outdir, ignore_errors=True)
    outdir.mkdir(parents=True)
    with open(snakemake.input.master, newline="") as fh:
        circular = {r["plasmid_id"] for r in csv.DictReader(fh, delimiter="\t")
                    if is_circular(r.get("topology"))}
    with open(snakemake.input.map, newline="") as fh:
        plasmids = {r["plasmid_id"] for r in csv.DictReader(fh, delimiter="\t")}
    # MacSyFinder 2.1.4 reads one "<replicon>: <topology>" per line (macsypy.database);
    # the gembase replicon name is the plasmid id with each '_' written as '-'.
    topology_file = outdir / "topology.txt"
    topology_file.write_text("".join(
        f"{p.replace('_', '-')}: {'circular' if p in circular else 'linear'}\n"
        for p in sorted(plasmids)))
    for family, options in FAMILIES.items():
        subprocess.run(
            ["macsyfinder", "--models-dir", str(models_dir), *options,
             "--sequence-db", snakemake.input.faa,
             "--db-type", "gembase", "--topology-file", str(topology_file),
             "--worker", str(snakemake.threads), "--out-dir", str(outdir / family),
             "--mute"],
            check=True)

back = {}
with open(snakemake.input.map, newline="") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        back[r["gembase_id"]] = r

rows = []
for path in outdir.rglob("best_solution.tsv"):
    with open(path) as fh:
        header = None
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if header is None:
                header = fields
                continue
            rec = dict(zip(header, fields))
            gid = rec.get("hit_id", "")
            src = back.get(gid)
            if not src:
                continue
            rows.append({"orf_id": src["orf_id"], "plasmid_id": src["plasmid_id"],
                         "gembase_id": gid, "system": rec.get("model_fqn", ""),
                         "system_id": rec.get("sys_id", ""),
                         "component": rec.get("gene_name", ""),
                         "hit_evalue": rec.get("hit_i_eval", ""),
                         "hit_status": rec.get("hit_status", ""),
                         "sys_wholeness": rec.get("sys_wholeness", ""),
                         "hit_gene_ref": rec.get("hit_gene_ref", ""),
                         "hit_profile_cov": rec.get("hit_profile_cov", ""),
                         "status": status.SUCCESS})
write_table(rows)

print(f"phase 2: {len(rows)} genes in {len({r['system_id'] for r in rows})} systems "
      f"across {len({r['plasmid_id'] for r in rows})} plasmids")
