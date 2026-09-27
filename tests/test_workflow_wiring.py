"""The workflow wiring: the rule graph Snakemake resolves on the test configuration, the
conda environments the rules declare, and the batch scripts."""
import os
import pathlib
import re
import subprocess

import pytest
import yaml
from conftest import requires

from plasmidann.tools import REQUIRED_TOOLS

ROOT = pathlib.Path(__file__).resolve().parents[1]
WORKFLOW = ROOT / "workflow"
SMK = [WORKFLOW / "Snakefile"] + sorted((WORKFLOW / "rules").glob("*.smk"))


def smk_text():
    return "\n".join(p.read_text() for p in SMK)


def declared_envs():
    """Every path named by a `conda:` directive, resolved relative to workflow/."""
    envs = set()
    for path in SMK:
        text = path.read_text()
        for match in re.finditer(r"conda:\s*\n\s*\"([^\"]+)\"", text):
            envs.add((path.parent / match.group(1)).resolve())
    return envs


def test_every_conda_directive_points_at_a_file_that_exists():
    missing = [str(e) for e in declared_envs() if not e.exists()]
    assert not missing, f"conda: directives naming absent environment files: {missing}"


def test_the_declared_environment_provides_every_required_tool():
    """Every declared environment provides every tool in plasmidann.tools.REQUIRED_TOOLS,
    including those a tool shells out to (IntegronFinder runs prodigal and cmsearch)."""
    # conda package name -> the executables it provides, where they differ from the
    # package name. A package may provide several.
    PROVIDES = {"hmmer": ("hmmsearch",), "mmseqs2": ("mmseqs",),
                "infernal": ("cmsearch",),
                "mdmparis-defense-finder": ("defense-finder",),
                "macsyfinder": ("macsyfinder",),
                "integron_finder": ("integron_finder",),
                "rnacode": ("RNAcode",),
                "eggnog-mapper": ("emapper.py",),
                "isescan": ("isescan.py",)}

    for env_path in declared_envs():
        spec = yaml.safe_load(env_path.read_text())
        packages = []
        for dep in spec["dependencies"]:
            if isinstance(dep, dict):
                packages.extend(dep.get("pip", []))
            else:
                packages.append(dep)
        provided = set()
        for pkg in packages:
            base = re.split(r"[=<>]", pkg)[0].strip()
            provided.update(PROVIDES.get(base, (base,)))

        for entry in REQUIRED_TOOLS:
            assert entry["name"] in provided, (
                f"{env_path.name} does not provide {entry['name']}, which "
                f"{entry['stage']} invokes: {entry['why']}")


def test_the_submission_script_activates_the_environment_it_declares():
    """The conda: directives take effect only if the submission script puts the built
    environment on PATH or passes --use-conda."""
    sbatch = (WORKFLOW / "run_pipeline.sbatch").read_text()
    assert ("envs/plasmidann/bin" in sbatch) or ("--use-conda" in sbatch), (
        "the submission script neither puts the built environment on PATH nor passes "
        "--use-conda, so the conda: directives do nothing")


def test_each_cascade_tier_is_one_job_with_every_core():
    """A search against a streamed database has a constant cost per invocation (DIAMOND reads
    the whole database each time), so a tier is one job and takes every core the run has."""
    assert "threads: workflow.cores" in _rule("tier_search"), (
        "tier_search does not take every core, so a single-job tier runs on a fraction "
        "of the allocation")


def test_only_s0_reads_the_configured_fasta():
    """The configured FASTA may hold the whole working set, simulated plasmids included, so
    only analysis_set reads it; every other rule takes the analysis-set FASTA it writes."""
    offenders = []
    for path in SMK:
        text = path.read_text()
        for match in re.finditer(
                r"config\[[\"']input[\"']\]\[[\"']fasta[\"']\]", text):
            rule = text[:match.start()].rsplit("\nrule ", 1)[-1].split(":", 1)[0]
            if rule != "analysis_set":
                offenders.append(f"{path.name}: rule {rule}")
    assert not offenders, (
        "these rules read the configured FASTA rather than the analysis-set FASTA S0 "
        "wrote, so their output covers plasmids outside the analysis scope: "
        + "; ".join(offenders))


def test_every_script_that_takes_a_scratch_directory_releases_it():
    """A scratch directory holds the intermediate databases of a clustering or a Foldseek
    run, so every script that takes one releases it (plasmidann.scratch). tempfile.mkdtemp
    is refused: plasmidann.scratch also clears the directory a failed attempt left behind,
    which mmseqs needs."""
    offenders = []
    for path in sorted((WORKFLOW / "scripts").glob("*.py")):
        text = path.read_text()
        if "tempfile.mkdtemp" in text:
            offenders.append(f"{path.name} calls tempfile.mkdtemp directly")
            continue
        if "scratch.scratch_dir" in text and "scratch.release" not in text:
            offenders.append(f"{path.name} takes a scratch directory and never releases it")
    assert not offenders, "; ".join(offenders)


def _rule(name):
    """The text of one rule, from its `rule` line to the next."""
    text = smk_text()
    body = text[text.index(f"rule {name}:"):]
    nxt = body.find("\nrule ", 1)
    return body if nxt < 0 else body[:nxt]


def test_one_submission_runs_every_stage_including_structure_search():
    """The user asked for one run that annotates everything. The structure search is
    omitted only when the optional GPU split is asked for explicitly."""
    sbatch = (WORKFLOW / "run_pipeline.sbatch").read_text()
    omit = sbatch[sbatch.index("OMIT=()"):sbatch.index("fi\n", sbatch.index("OMIT=()"))]
    assert "STRUCTURE_ON_GPU" in omit, (
        "structure_search is omitted from the default submission")


def _operative_text(path):
    """What a file DOES, without its prose: comments and docstrings removed.

    Python through ast (every identifier and every string that is not a docstring), rule
    files by stripping triple-quoted blocks and comments, YAML by loading it. Prose may say
    why PlasAnn is not used; an operative reference is what must not exist."""
    import ast
    text = path.read_text()
    if path.suffix == ".py":
        tree = ast.parse(text)
        docstrings = {id(n.value) for n in ast.walk(tree) if isinstance(n, ast.Expr)
                      and isinstance(getattr(n, "value", None), ast.Constant)}
        parts = []
        for node in ast.walk(tree):
            if isinstance(node, ast.Constant) and isinstance(node.value, str) \
                    and id(node) not in docstrings:
                parts.append(node.value)
            elif isinstance(node, ast.Name):
                parts.append(node.id)
            elif isinstance(node, ast.Attribute):
                parts.append(node.attr)
            elif isinstance(node, (ast.FunctionDef, ast.ClassDef, ast.alias)):
                parts.append(getattr(node, "name", ""))
        return "\n".join(parts)
    if path.suffix == ".yaml":
        return yaml.safe_dump(yaml.safe_load(text))
    text = re.sub(r'"""(?:.|\n)*?"""', "", text)
    return "\n".join(line.split("#", 1)[0] for line in text.splitlines())


def test_the_pipeline_never_uses_plasann_labels_or_kegg_context_terms():
    """User decision 2026-09-25: PlasAnn's database and labels are not used (only its
    published tier thresholds), and KEGG gives no context term. The KEGG KO list maps Tier 0
    KOs to gene symbols for the disagreement table, and nothing else."""
    from plasmidann import context_terms
    code = [WORKFLOW / "Snakefile", *sorted((WORKFLOW / "rules").glob("*.smk")),
            *sorted((WORKFLOW / "scripts").glob("*.py")),
            *sorted((ROOT / "src").rglob("*.py")),
            *sorted((ROOT / "config").rglob("*.yaml"))]
    offenders = [str(p.relative_to(ROOT)) for p in code
                 if "plasann" in _operative_text(p).lower()]
    assert not offenders, f"PlasAnn referenced operatively in: {offenders}"

    # KEGG: no context term type, and the KO list reaches only protein_labels.
    assert not any("kegg" in t.lower() for t in context_terms.TERM_PREFIX.values())
    assert not any("kegg" in k.lower() for k in context_terms.TERM_PREFIX)
    readers = [p.name for p in code if "kegg_ko_list" in _operative_text(p)
               and p.suffix != ".yaml"]
    assert readers == ["annotation_cascade.smk"], readers
    assert "kegg_ko_list" in _rule("protein_labels")
    context_code = _operative_text(WORKFLOW / "scripts" / "context_features.py") \
        + _operative_text(ROOT / "src" / "plasmidann" / "context_terms.py")
    assert "kegg" not in context_code.lower(), "KEGG reaches the context terms"


# Inputs that no rule produces, on the test configuration. The rule graph is built over
# empty stand-ins for them, so the test needs no reference data.
EXTERNAL_INPUTS = [("input", "master_table"), ("input", "fasta"),
                   ("input", "plasmidscope_proteins"), ("input", "host_provenance"),
                   ("input", "working_set"), ("references", "pfam_dat"),
                   ("references", "kegg_ko_list")]


@pytest.fixture(scope="module")
def rule_graph(tmp_path_factory):
    """(producer, consumer) rule pairs of the graph Snakemake resolves for `all` on the
    test configuration, in a scratch copy of the repository (symlinks), so nothing is
    written into the checkout. Snakemake's source cache goes to the scratch directory too:
    the home directory on this cluster is at its file quota."""
    from snakemake.utils import update_config
    work = tmp_path_factory.mktemp("dag")
    for name in ("config", "src", "workflow"):
        (work / name).symlink_to(ROOT / name)
    config = yaml.safe_load((ROOT / "config" / "config.yaml").read_text())
    update_config(config, yaml.safe_load((ROOT / "config" / "test" / "config.yaml").read_text()))
    for path in [config[a][b] for a, b in EXTERNAL_INPUTS]:
        (work / path).parent.mkdir(parents=True, exist_ok=True)
        (work / path).touch()
    snakemake = ROOT / "envs" / "plasmidann" / "bin" / "snakemake"
    proc = subprocess.run(
        [str(snakemake) if snakemake.exists() else "snakemake", "-s", "workflow/Snakefile",
         "--configfile", "config/test/config.yaml", "-n", "-c", "1", "--rulegraph"],
        cwd=work, env={**os.environ, "XDG_CACHE_HOME": str(work / ".cache")},
        capture_output=True, text=True, timeout=600)
    assert proc.returncode == 0, f"the DAG does not resolve:\n{proc.stderr[-3000:]}"
    names = dict(re.findall(r'^\s*(\d+)\[label = "(\w+)"', proc.stdout, re.M))
    return {(names[a], names[b])
            for a, b in re.findall(r"^\s*(\d+) -> (\d+)", proc.stdout, re.M)}


@requires("snakemake")
def test_every_rule_is_reached_by_all(rule_graph):
    defined = set(re.findall(r"^rule (\w+):", smk_text(), re.M)) - {"annotate_only"}
    in_graph = {rule for edge in rule_graph for rule in edge}
    assert in_graph == defined, (
        f"not reached by `all`: {sorted(defined - in_graph)}; "
        f"unknown: {sorted(in_graph - defined)}")


@requires("snakemake")
@pytest.mark.parametrize("producer, consumer", [
    # No search starts before pre-flight has confirmed every tool and database.
    ("preflight", "tier_search"), ("preflight", "artefact_screen"),
    ("preflight", "label_databases"), ("preflight", "conjugation_systems"),
    # No search uses an unconfirmed -Z.
    ("check_hmmer_z", "artefact_screen"),
    # AntiFam-flagged proteins are left out of the cascade selection, and the DIAMOND
    # tiers skip every artefact-flagged protein.
    ("artefact_screen", "cascade_selection"), ("artefact_screen", "tier_search"),
    # The search representatives are the first tier's query.
    ("cascade_selection", "tier_search"),
    # The dark set is the target-eligible proteins.
    ("target_eligibility", "dark_set"),
    # The label table merges the label databases, defence and CONJScan components.
    ("label_databases", "protein_labels"), ("defence_systems", "protein_labels"),
    ("conjugation_systems", "protein_labels"),
    # Context terms: the labels, the systems, every family and the lineages.
    ("protein_labels", "context_features"), ("conjugation_systems", "context_features"),
    ("protein_families", "context_features"), ("plasmid_lineage", "context_features"),
    # Synteny: one cluster table per level, counted over lineages, dark proteins only.
    ("protein_clustering", "synteny"), ("plasmid_lineage", "synteny"),
    ("dark_set", "synteny"),
    ("protein_families", "dark_cooccurrence"), ("plasmid_lineage", "dark_cooccurrence"),
    # The report joins every per-ORF and per-family evidence table.
    ("label_databases", "annotation_report"), ("conjugation_systems", "annotation_report"),
    ("protein_clustering", "annotation_report"), ("dark_cooccurrence", "annotation_report"),
    ("recurrence", "annotation_report"), ("synteny", "annotation_report"),
    ("rarity", "annotation_report"), ("structure_search", "annotation_report"),
])
def test_the_rule_graph_has_the_dependency(rule_graph, producer, consumer):
    assert (producer, consumer) in rule_graph


@pytest.mark.parametrize("script", ["run_pipeline.sbatch", "bench_pipeline.sbatch"])
def test_the_slurm_log_goes_to_the_submit_directory(script):
    """Slurm opens the -o file before the script runs and fails without output when its
    directory does not exist, as results/logs does not on a fresh checkout."""
    (path,) = re.findall(r"^#SBATCH -o (\S+)", (WORKFLOW / script).read_text(), re.M)
    assert "/" not in path, f"{script} writes its Slurm log into {path}"
