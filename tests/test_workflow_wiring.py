"""The workflow files themselves must be internally consistent.

These are cheap text checks over Snakefile, the rule modules and the batch scripts. They
exist because three of the defects that would have killed a real run were not in any
function - they were in the wiring: conda environments that no submission script ever
activated, a rule that declared an environment missing the tool it invokes, and a
cascade tier that streamed a 375 GB database once per shard.
"""
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
    """A per-rule environment that omits a tool the rule invokes is worse than no
    environment at all: the rule fails at run time with the tool visibly installed on the
    machine. IntegronFinder failed exactly this way - its environment omitted prodigal and
    infernal, which it shells out to."""
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
    """26 rules carried conda directives and no submission script passed --use-conda, so
    every one of them was inert and the run silently depended on whatever happened to be on
    PATH. Either the script activates the environment on PATH or it asks Snakemake to
    build it; doing neither is the failure this checks for."""
    sbatch = (WORKFLOW / "run_pipeline.sbatch").read_text()
    assert ("envs/plasmidann/bin" in sbatch) or ("--use-conda" in sbatch), (
        "the submission script neither puts the built environment on PATH nor passes "
        "--use-conda, so the conda: directives do nothing")


def test_each_cascade_tier_is_one_job_with_every_core():
    """A search against a streamed database has a fixed cost per invocation - DIAMOND
    reads the whole of nr each time it runs, and 2,000 queries took more than 12 hours
    on 16 threads. Sharding the query set multiplies that cost by the shard count, so a
    tier is one job, and it must be given every core the run has or it runs on one."""
    text = (WORKFLOW / "rules" / "annotation_cascade.smk").read_text()
    rule = text[text.index("rule tier_search"):]
    rule = rule[:rule.index("\nrule ")]
    assert "{cshard}" not in rule, (
        "tier_search is sharded again: every shard streams the database")
    assert "threads: workflow.cores" in rule, (
        "tier_search does not take every core, so a single-job tier runs on a fraction "
        "of the allocation")


def test_the_codon_model_stage_is_gone():
    """S7d - FastTree plus HyPhy BUSTED - was removed from the analysis.

    Removal has to be complete across four files or the run breaks somewhere different
    each time: a script invoking a tool the environment no longer pins, a rule writing an
    output nothing consumes, a script reading an input no rule declares, or a config key
    nothing reads. The fifth failure - a pre-flight that halts on a tool no stage needs -
    is covered generically by test_tools_registry's stale-registry guard.

    Prose is not checked. A docstring may say why the stage was removed; what must not
    survive is an operative reference.
    """
    offenders = []

    for path in sorted((WORKFLOW / "scripts").glob("*.py")):
        for tool in re.findall(r"subprocess\.run\(\s*\n?\s*f?\"([A-Za-z_][\w.-]*)",
                               path.read_text()):
            if tool in ("FastTree", "FastTreeMP", "hyphy", "HYPHYMPI"):
                offenders.append(f"{path.name} invokes {tool}")

    spec = yaml.safe_load((WORKFLOW / "envs" / "plasmidann.yaml").read_text())
    for dep in spec["dependencies"]:
        if isinstance(dep, str) and re.split(r"[=<>]", dep)[0] in ("fasttree", "hyphy"):
            offenders.append(f"the environment still pins {dep}")

    if "rule busted_confirm" in smk_text():
        offenders.append("rule busted_confirm is still defined")

    targets = yaml.safe_load((ROOT / "config" / "targets.yaml").read_text())
    for key in targets.get("evolution", {}):
        if "busted" in key:
            offenders.append(f"config/targets.yaml still declares evolution.{key}")

    assert not offenders, "S7d was not fully removed: " + "; ".join(offenders)


@requires("snakemake")
@pytest.mark.slow
def test_the_resolved_dag_matches_the_rules_that_exist(tmp_path):
    """The text checks above cannot see a dangling input.

    Leave `busted=f"{OUT}/11_distribution_and_evolution/busted.tsv"` in rule annotation_report's input block after
    deleting rule busted_confirm and every other test in this file still passes: the
    substring "rule busted_confirm" is gone, no script invokes the tool, and the smoke test
    drives annotation_report.py through a hand-built input dict that never had that key.
    Only Snakemake's own DAG resolution fails, with a missing-input error, and until now
    that check lived in a shell history rather than the suite.

    A dry-run is the cheapest thing that exercises the resolved DAG. Snakemake keeps a
    source cache under $XDG_CACHE_HOME, which on this cluster sits in a home directory at
    its file quota, so the cache is pointed at tmp_path.
    """
    env = {**os.environ, "XDG_CACHE_HOME": str(tmp_path)}
    proc = subprocess.run(["snakemake", "-n", "-c", "1", "--quiet", "rules"], cwd=ROOT, env=env,
                          capture_output=True, text=True, timeout=600)
    assert proc.returncode == 0, f"dry-run failed:\n{proc.stderr[-3000:]}"

    # `--quiet rules` prints the job table; rule names are the first token of each row.
    jobs = {line.split()[0] for line in proc.stdout.splitlines()
            if line and not line.startswith(("Job stats", "job", "-", "total"))}
    assert "busted_confirm" not in jobs, "the DAG still schedules the removed stage"
    for rule in ("family_evolution", "consensus_recheck", "annotation_report"):
        assert rule in jobs, f"{rule} is missing from the resolved DAG"


def test_only_s0_reads_the_configured_fasta():
    """The analysis set defines the scope; the configured FASTA does not.

    rule feature_files once streamed the configured corpus and wrote one record per
    sequence it found there. On the 100-plasmid test configuration that produced 208,245
    GenBank records - 11.9 GB - for a run that had been asked to look at 100 plasmids.
    The TSV written by the same stage was correctly scoped, which is what kept the defect
    out of sight: every count a reader checks came from the TSV.

    The configured FASTA may hold the whole working set, simulated plasmids included, so
    exactly one rule may read it: analysis_set, which applies the exclusion and writes the
    FASTA every other stage takes.
    """
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
    """Six scripts created a working directory beside their output and none removed it.

    On the 100-plasmid test set that left seven directories behind. At production scale
    the same directories hold the intermediate databases of a 3.5M-protein clustering and
    a foldseek run, so the leak is hundreds of gigabytes - and it accumulates across
    reruns, because each attempt makes its own.

    mkdtemp is banned outright rather than merely paired with a cleanup: plasmidann.scratch
    also clears the stale directory a previous FAILURE deliberately left behind, which a
    bare mkdtemp does not, and mmseqs refuses to start against a mismatched tmp directory.
    """
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


def test_the_label_and_conjugation_stages_are_rules_with_declared_resources():
    """S4d and S8f run inside the one submission, so the scheduler must know what each
    takes, and neither may start a search before pre-flight has passed."""
    for name, script in (("label_databases", "label_databases.py"),
                         ("conjugation_systems", "conjugation_systems.py")):
        rule = _rule(name)
        assert f'"../scripts/{script}"' in rule
        assert "threads:" in rule and "mem_mb=" in rule and "runtime=" in rule, name
        assert "preflight.tsv" in rule, f"{name} can start before pre-flight passes"


def test_protein_labels_merges_the_label_databases_and_writes_the_disagreements():
    rule = _rule("protein_labels")
    for needed in ("protein_labels_plasmid.tsv", "label_disagreements.tsv",
                   "kegg_ko_list", "defence_systems.tsv", "conjugation_systems.tsv",
                   "protein_map.tsv"):
        assert needed in rule, f"protein_labels does not declare {needed}"


def test_synteny_reads_one_cluster_table_per_level_and_the_lineages():
    rule = _rule("synteny")
    assert 'res=targets["synteny"]["levels"]' in rule
    assert "plasmid_lineage.tsv" in rule and "dark_ids.txt" in rule
    assert 'synteny=targets["synteny"]' in rule


def test_context_features_reads_the_labels_and_conjugation_and_writes_the_terms():
    rule = _rule("context_features")
    for needed in ("conjugation_systems.tsv", "protein_labels.tsv", "protein_families.tsv",
                   "plasmid_lineage.tsv", "family_context_terms.tsv",
                   'primary=targets["clustering"]["primary"]'):
        assert needed in rule, f"context_features does not declare {needed}"


def test_the_report_reads_the_new_evidence():
    rule = _rule("annotation_report")
    for needed in ("conjugation_systems.tsv", "conjugation_plasmid_class.tsv",
                   "protein_labels_plasmid.tsv", "families_close_cluster.tsv"):
        assert needed in rule, f"annotation_report does not declare {needed}"


def test_one_submission_runs_every_stage_including_structure_search():
    """The user asked for one run that annotates everything. The structure search is
    omitted only when the optional GPU split is asked for explicitly."""
    sbatch = (WORKFLOW / "run_pipeline.sbatch").read_text()
    omit = sbatch[sbatch.index("OMIT=()"):sbatch.index("fi\n", sbatch.index("OMIT=()"))]
    assert "STRUCTURE_ON_GPU" in omit, (
        "structure_search is omitted from the default submission")
    assert "--until" not in sbatch.split("\nset -euo pipefail", 1)[1]


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
