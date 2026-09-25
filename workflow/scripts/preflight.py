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
import os
import pathlib
import shutil
import subprocess

import _ctx  # noqa: F401

from plasmidann import labeldb
from plasmidann.conjscan import installed_version
from plasmidann.tools import (grammar_problem, macsyfinder_version, model_grammars,
                              required_tools)

tiers = snakemake.params.tiers
artefact = snakemake.params.artefact
structure = snakemake.params.get("structure") or {}
structure_required = bool(structure.get("required"))
orthology = snakemake.params.get("orthology") or {}
orthology_required = bool(orthology.get("required"))

# method -> the executable that method actually invokes
EXECUTABLE = {"hmmer": "hmmsearch", "diamond": "diamond", "pharokka": None}

# What a pharokka database directory must hold for protein mode: the phage families as an
# MMseqs2 profile database and as HMMER3 profiles, the annotation table that names and
# categorises them, and the CARD and VFDB databases searched in the same run.
PHAROKKA_DB_FILES = ("phrogs_profile_db", "all_phrogs.h3m", "phrog_annot_v4.tsv",
                     "CARD", "vfdb")

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
    if tier["method"] == "pharokka":
        # pharokka lives in its own environment and is named by path, not found on PATH.
        exe = tier.get("exe", "")
        if not (exe and os.access(exe, os.X_OK)):
            problems.append(f"{tier['id']}: pharokka executable not found or not "
                            f"executable: {exe!r} - build workflow/envs/pharokka.yaml")
        missing = [f for f in PHAROKKA_DB_FILES if not os.path.exists(os.path.join(db, f))]
        if missing:
            problems.append(f"{tier['id']}: pharokka database {db} is missing "
                            f"{', '.join(missing)} - run `pharokka install -o {db}`")
        continue
    if not os.path.exists(db):
        problems.append(f"{tier['id']}: database not found: {db}")
    elif tier["method"] == "diamond" and tier.get("expected_sequences"):
        # Existence is not completeness: a makedb killed part-way leaves a non-empty
        # .dmnd, which would first fail - or silently search a fraction of the database -
        # when the tier starts, a day or more into the run.
        info = subprocess.run(["diamond", "dbinfo", "-d", db], capture_output=True,
                              text=True).stdout
        n = next((int(l.split()[-1]) for l in info.splitlines()
                  if l.strip().startswith("Sequences")), None)
        if n != tier["expected_sequences"]:
            problems.append(f"{tier['id']}: {db} holds {n} sequences, expected "
                            f"{tier['expected_sequences']} - the index is incomplete; "
                            "rebuild it")
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

# ------------------------------------------------------------------------------------
# S4d plasmid label databases, and AMRFinderPlus in its own environment. The same contract
# as the stage (label_databases.py): an absent database directory fails only when the
# databases are required, and otherwise becomes NOT_RUN; a directory that is present but
# incomplete always fails, because the stage would stop on it hours into the run.
# ------------------------------------------------------------------------------------
labels_cfg = snakemake.params.labels
labels_dir = labels_cfg["dir"]
for db in labeldb.DATABASES:
    d = os.path.join(labels_dir, db)
    if not os.path.isdir(d):
        if labels_cfg["required"]:
            problems.append(f"label database {db}: {d} is not installed - run "
                            "tools/download_label_dbs.py, or set labels.required to false")
        continue
    version = os.path.join(d, "VERSION")
    if not (os.path.isfile(version) and open(version).read().strip()):
        problems.append(f"label database {db}: {version} is missing or empty; every label "
                        "carries its release")
    # CARD is read from card.json (sequences and curated cut-offs); the others from the
    # FASTA the installer writes; BacMet also needs its compound mapping.
    needed = ["card.json"] if db == "card" else [f"{db}.faa"]
    for name in needed:
        path = os.path.join(d, name)
        if not (os.path.isfile(path) and os.path.getsize(path) > 0):
            problems.append(f"label database {db}: {path} is missing or empty")
    if db == "bacmet" and not list(pathlib.Path(d).rglob("*mapping*.txt")):
        problems.append(f"label database bacmet: no BacMet mapping file (*mapping*.txt) "
                        f"under {d}")

amr = snakemake.params.amrfinder
amr_missing = []
if not os.access(amr["executable"], os.X_OK):
    amr_missing.append(f"AMRFinderPlus executable not found or not executable: "
                       f"{amr['executable']!r} - run tools/install_tool_envs.py")
if not os.path.isfile(os.path.join(amr["database"], "version.txt")):
    amr_missing.append(f"AMRFinderPlus database {amr['database']} has no version.txt - "
                       "run tools/install_tool_envs.py")
if amr["required"]:
    problems.extend(amr_missing)

# ------------------------------------------------------------------------------------
# S8f CONJScan, run by the MacSyFinder named in conjugation.exe. The models' grammar
# decides which MacSyFinder can read them, and the wrong one stops with a parse error
# (plasmidann.tools.MIN_MACSYFINDER). Checked whenever both are installed, required or
# not: the stage runs then, and would fail after the cascade rather than here.
# ------------------------------------------------------------------------------------
conj = snakemake.params.conjugation
conj_models = snakemake.params.conjscan_models
conj_exe = conj["exe"]
conj_version = installed_version(conj_models)
conj_exe_ok = os.access(conj_exe, os.X_OK)
if conj["required"]:
    if not conj_exe_ok:
        problems.append(f"conjugation.required is true but the MacSyFinder executable is "
                        f"absent or not executable: {conj_exe!r} - run "
                        "tools/install_tool_envs.py")
    if conj_version is None:
        problems.append(f"conjugation.required is true but CONJScan is not installed "
                        f"under {conj_models}")
macsyfinder_found = ""
if conj_exe_ok and conj_version is not None:
    if conj_version != str(conj["version"]):
        problems.append(f"CONJScan {conj_version} is installed at {conj_models}, but "
                        f"conjugation.version is {conj['version']}")
    grammars = model_grammars(conj_models)
    # MacSyFinder's start-up imports took up to 2 min 53 s on /gorilla, so a slow answer
    # is waited for; only an hour without one is a failure.
    try:
        answer = subprocess.run([conj_exe, "--version"], capture_output=True, text=True,
                                timeout=3600)
        found = macsyfinder_version(answer.stdout + answer.stderr)
    except subprocess.TimeoutExpired:
        found = None
    if not grammars:
        problems.append(f"CONJScan at {conj_models} holds no model definition")
    elif found is None:
        problems.append(f"{conj_exe} --version reported no MacSyFinder version")
    else:
        macsyfinder_found = ".".join(map(str, found))
        problem = grammar_problem(grammars, found)
        if problem:
            problems.append(f"CONJScan at {conj_models}: {problem} ({conj_exe}). Use the "
                            "MacSyFinder of envs/conjscan (tools/install_tool_envs.py).")

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
        # A pharokka database is a directory; report the size of its profile database.
        path = os.path.join(db, "phrogs_profile_db") if os.path.isdir(db) else db
        out.write(f"{db}\t{os.path.getsize(path) if os.path.exists(path) else 0}\n")
    out.write(f"\nstructure_required\t{structure_required}\n")
    # What the S4d and S8f stages will find, so a NOT_RUN in their outputs can be traced.
    out.write(f"labels_dir\t{labels_dir}\n")
    out.write(f"amrfinder\t{'; '.join(amr_missing) or amr['executable']}\n")
    out.write(f"conjscan\t{conj_version or 'not installed'}\n")
    out.write(f"conjscan_macsyfinder\t{macsyfinder_found or 'not checked'}\n")

print(f"pre-flight OK: {len(resolved)} tool(s), {len(checked_dbs) + 1} database(s)")
