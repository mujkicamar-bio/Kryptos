"""S8a phase 2: call defence systems from gene adjacency.

MacSyFinder is run directly rather than through `defense-finder run`, because the wrapper
does not pass `--replicon-topology` and 94% of these plasmids are circular; under linear
topology a system spanning the origin is not found. `--db-type gembase` holds every
candidate replicon in one database and still treats each separately.

One MacSyFinder process searches all candidate replicons, with the cores given to
--worker. HMMER's independent e-value is the p-value times the number of sequences in the
database, and MacSyFinder keeps a hit only below --i-evalue-sel (0.001, its default), so
splitting the database into per-core chunks would make the calls depend on the core count.
The e-values still depend on the size of the candidate set, as for any single database.

Thresholds and the quorum and co-localisation rules are those of MacSyFinder (Abby et al.
2014) and the 711 shipped DefenseFinder models (Tesson et al. 2022), not overridden.

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

COLUMNS = ["orf_id", "plasmid_id", "gembase_id", "system", "system_id", "component",
           "hit_evalue", "hit_status", "sys_wholeness", "hit_gene_ref", "hit_profile_cov",
           "status"]
outdir = pathlib.Path(snakemake.output.tsv).parent / "phase2"


def write_table(rows):
    with open(snakemake.output.tsv, "w", newline="") as out:
        w = csv.DictWriter(out, fieldnames=COLUMNS, delimiter="\t")
        w.writeheader()
        w.writerows(rows)


models_dir = pathlib.Path(snakemake.params.models_dir)
if not (models_dir.is_dir() and any(models_dir.iterdir())):
    if snakemake.params.get("required", False):
        sys.exit(f"S8a phase 2: DefenseFinder models are not installed at {models_dir}.")
    print(f"S8a phase 2: models absent at {models_dir}; recording NOT_RUN.")
    write_table([{"status": status.NOT_RUN}])
    sys.exit(0)

# skip_run is set only by the parser test, which pre-populates the output tree.
if not snakemake.params.get("skip_run", False):
    # A rerun starts clean: results from an interrupted run would be read below.
    shutil.rmtree(outdir, ignore_errors=True)
    outdir.mkdir(parents=True)
    subprocess.run(
        f"macsyfinder --models-dir {snakemake.params.models_dir} "
        f"--models defense-finder-models all "
        f"--sequence-db {snakemake.input.faa} "
        f"--db-type gembase --replicon-topology circular "
        f"--worker {snakemake.threads} --out-dir {outdir / 'run'} --mute",
        shell=True, check=True)

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
