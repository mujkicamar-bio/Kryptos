"""S3 pre-flight: fail in seconds if a tool or database is missing, not after 45 hours.

v1's most expensive failure was a job that ran for 45 hours and then died because DIAMOND
was not on PATH - the module load line had been written for HMMER and MMseqs2 and never
updated when the tier list changed. The identical failure occurred twice.

This rule is a dependency of every tier, so the DAG cannot start a search until every tool
and every database the run needs has been confirmed present. It costs about a second.

IT CHECKS THE WHOLE PIPELINE, NOT JUST THE CASCADE.

An earlier version walked only the configured tier list, which is four tools of the eleven
the workflow runs. That left the original failure class wide open one stage further down: a
missing mafft still killed S7b after the cascade had run for days, a missing prodigal made
all 600 IntegronFinder shards fail, and a missing foldseek silently emptied the novel_fold
stratum. The registry now lives in plasmidann.tools and a test scans workflow/scripts/ for
subprocess calls, so it cannot fall behind the code again.
"""
import _ctx  # noqa: F401
import os
import shutil

from plasmidann.tools import required_tools

tiers = snakemake.params.tiers
artefact = snakemake.params.artefact
structure = snakemake.params.get("structure") or {}
structure_required = bool(structure.get("required"))
orthology = snakemake.params.get("orthology") or {}
orthology_required = bool(orthology.get("required"))

# method -> the executable that method actually invokes
EXECUTABLE = {"hmmer": "hmmsearch", "hhsearch": "hhblits_omp", "diamond": "diamond"}

# An hhsearch tier's `db` is a PREFIX, not a file: HH-suite splits a database across
# <prefix>_hhm.ff{data,index} and <prefix>_cs219.ff{data,index}. Checking the prefix with
# os.path.exists reports "database not found" for a database that is present and complete.
HHSUITE_PARTS = ("_hhm.ffdata", "_hhm.ffindex", "_cs219.ffdata", "_cs219.ffindex")

problems = []
resolved = {}

# ------------------------------------------------------------------------------------
# Every executable the whole workflow runs, not only the cascade's.
# ------------------------------------------------------------------------------------
for entry in required_tools(structure_required=structure_required,
                            orthology_required=orthology_required):
    path = shutil.which(entry["name"])
    resolved[entry["name"]] = path
    if path is None:
        problems.append(
            f"{entry['stage']}: executable {entry['name']!r} not found on PATH. "
            f"{entry['why']}. Activate envs/plasmidann or load the module.")

# ------------------------------------------------------------------------------------
# Cascade databases. A tier naming an unknown method is a config error, not a missing tool.
# ------------------------------------------------------------------------------------
checked_dbs = set()
for tier in tiers:
    if tier["method"] not in EXECUTABLE:
        problems.append(f"{tier['id']}: unknown method {tier['method']!r}")
        continue
    db = tier["db"]
    if db in checked_dbs:
        continue
    checked_dbs.add(db)
    if tier["method"] == "hhsearch":
        # The cs219 pair is the prefilter. Without it hhblits scores every query against
        # all 38,880 profiles, which is not a failure - it is the same search, orders of
        # magnitude slower - so its absence is reported as the missing piece it is.
        missing = [part for part in HHSUITE_PARTS if not os.path.exists(db + part)]
        if missing:
            problems.append(
                f"{tier['id']}: HH-suite database {db} is incomplete, missing "
                f"{', '.join(missing)} - build it with tools/build_phrogs_hhsuite.sh")
        continue
    if not os.path.exists(db):
        problems.append(f"{tier['id']}: database not found: {db}")
    elif tier["method"] == "hmmer" and not os.path.exists(db + ".h3i"):
        # An unpressed HMM library makes hmmsearch re-parse a multi-gigabyte flat file on
        # every shard. Pfam-A is 2.2 GB; the pressed index is what makes it tractable.
        problems.append(f"{tier['id']}: {db} is not pressed - run `hmmpress {db}`")

# The artefact screen has its own database and runs before the cascade.
if not os.path.exists(artefact["antifam_db"]):
    problems.append(f"artefact screen: AntiFam not found: {artefact['antifam_db']}")
elif not os.path.exists(artefact["antifam_db"] + ".h3i"):
    problems.append(f"artefact screen: {artefact['antifam_db']} is not pressed")

# Structure evidence. When declared required, a missing database fails HERE, in seconds,
# rather than producing an empty table after the cascade has already run for days. When not
# required the stage is skipped entirely and structure is recorded as NOT_RUN - which is
# honest - rather than "searched and found nothing", which is not.
if orthology_required and not os.path.isdir(orthology.get("data_dir", "")):
    problems.append(
        f"orthology.required is true but the eggNOG data directory is missing: "
        f"{orthology.get('data_dir')}. Run `download_eggnog_data.py --data_dir "
        f"{orthology.get('data_dir')}` (~50 GB), or set orthology.required to false in "
        "config/targets.yaml to run without KEGG and COG terms.")

if structure_required:
    for label, path in (("foldseek target database", snakemake.params.foldseek_db),
                        ("ProstT5 model", snakemake.params.prostt5)):
        if not os.path.exists(path):
            problems.append(
                f"structure.required is true but the {label} is missing: {path}. "
                "Download it, or set structure.required to false in config/targets.yaml "
                "to run without structural evidence.")

if problems:
    raise SystemExit(
        "pre-flight failed - the run would have died mid-search:\n  "
        + "\n  ".join(problems))

with open(snakemake.output[0], "w") as out:
    out.write("tool\tresolved_path\n")
    for tool in sorted(resolved):
        out.write(f"{tool}\t{resolved[tool]}\n")
    out.write("\ndatabase\tsize_bytes\n")
    for db in sorted(checked_dbs | {artefact["antifam_db"]}):
        # An HH-suite database is a prefix; its size is the profiles it holds.
        path = db if os.path.exists(db) else db + "_hhm.ffdata"
        out.write(f"{db}\t{os.path.getsize(path) if os.path.exists(path) else 0}\n")
    out.write(f"\nstructure_required\t{structure_required}\n")

print(f"pre-flight OK: {len(resolved)} tool(s), {len(checked_dbs) + 1} database(s)")
