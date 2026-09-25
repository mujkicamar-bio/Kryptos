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
