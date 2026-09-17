"""Test harness for the Snakemake script layer.

WHY THIS EXISTS

172 unit tests covered src/plasmidann/ and not one line of workflow/scripts/. Five
independent reviews then found that essentially every serious defect lived in that
untested layer - including four stages that produced no output while reporting success,
and a NameError that had never fired because a missing tool masked it.

The scripts are not importable modules: Snakemake injects a global named `snakemake` and
executes them. This harness supplies that global, so a script can be run against a ten-row
fixture in milliseconds. That is the whole trick, and it is what was missing.
"""
import os
import pathlib
import shutil
import sys
from types import SimpleNamespace

import pytest

ROOT = pathlib.Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "workflow" / "scripts"
TOOLBIN = ROOT / "envs" / "plasmidann" / "bin"


class NamedList(list):
    """Stand-in for snakemake.io.Namedlist.

    Supports positional indexing, attribute access by name, and .get() - all three of
    which the scripts use. Named values are also appended to the list, matching
    Snakemake's own behaviour, so `snakemake.output[0]` works alongside
    `snakemake.output.hits`.
    """

    def __init__(self, *positional, **named):
        items = list(positional)
        for value in named.values():
            if isinstance(value, (list, tuple)):
                items.extend(value)
            else:
                items.append(value)
        super().__init__(items)
        self._named = dict(named)
        for key, value in named.items():
            setattr(self, key, value)

    def get(self, key, default=None):
        return self._named.get(key, default)


def requires(*tools):
    """Skip a test when an executable it drives is not installed.

    Some smoke tests run the real tool - a three-sequence DIAMOND database, mafft over three
    peptides, mmseqs over two. On a machine with the environment built they are the most
    valuable tests here, because they exercise the code path the cluster will take. On a
    fresh clone with no environment yet they should SKIP, not fail: a red suite on checkout
    tells a reader nothing about the code.

    The check looks in the project environment as well as on PATH, because the session
    fixture that prepends it has not run yet at collection time.
    """
    missing = [t for t in tools
               if shutil.which(t) is None and not (TOOLBIN / t).exists()]
    return pytest.mark.skipif(
        bool(missing),
        reason=f"not installed: {', '.join(missing)} - build workflow/envs/plasmidann.yaml")


def _as_namedlist(spec):
    """Build a Namedlist from either a dict (named) or a list/tuple (positional)."""
    if spec is None:
        return NamedList()
    if isinstance(spec, dict):
        return NamedList(**spec)
    return NamedList(*spec)


class FakeSnakemake:
    """The object Snakemake injects into a script's global namespace."""

    def __init__(self, input=None, output=None, params=None, log=None,
                 threads=1, wildcards=None):
        # Rules declare inputs and outputs either positionally or by name, and scripts
        # access them the same way, so the harness has to accept both forms.
        self.input = _as_namedlist(input)
        self.output = _as_namedlist(output)
        self.params = _as_namedlist(params)
        self.log = NamedList(*(log or []))
        self.threads = threads
        self.wildcards = SimpleNamespace(**(wildcards or {}))


def run_script(name, snake):
    """Execute a workflow script with `snakemake` bound, exactly as Snakemake does.

    Any exception propagates: a script that raises here would have raised on the cluster,
    which is the entire point. A script that writes nothing is caught by the caller's
    assertions rather than by silence.
    """
    path = SCRIPTS / name
    for p in (str(SCRIPTS), str(ROOT / "src")):
        if p not in sys.path:
            sys.path.insert(0, p)
    # Snakemake runs each script in its own process, so _ctx's one-time setup - the src
    # path and the redirection of stdout into the rule's log - happens once per script.
    # Inside one pytest process the module would be cached after the first test and the
    # setup would silently not run for any later one, so the import is made fresh here.
    sys.modules.pop("_ctx", None)
    saved_out, saved_err = sys.stdout, sys.stderr
    globals_dict = {"snakemake": snake, "__name__": "__main__", "__file__": str(path)}
    try:
        exec(compile(path.read_text(), str(path), "exec"), globals_dict)
    finally:
        sys.stdout, sys.stderr = saved_out, saved_err
    return globals_dict


def write_tsv(path, header, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as fh:
        fh.write("\t".join(header) + "\n")
        for r in rows:
            fh.write("\t".join(str(x) for x in r) + "\n")


def write_fasta(path, records):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as fh:
        for name, seq in records:
            fh.write(f">{name}\n{seq}\n")


def read_tsv(path):
    import csv
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


@pytest.fixture(scope="session", autouse=True)
def project_tools_on_path():
    """Tools live in the project environment, not on the default PATH.

    Session-scoped and autouse so it also covers fixtures that build test data with a tool
    (a three-sequence DIAMOND database, for instance) before any script runs. Putting the
    real binaries on PATH means a missing tool fails the test rather than silently routing
    a script down an error branch that exits 0 - which is precisely how the collections
    NameError in family_evolution.py hid through five reviews."""
    os.environ["PATH"] = f"{TOOLBIN}{os.pathsep}{os.environ['PATH']}"


@pytest.fixture
def fixture_dir(tmp_path):
    return tmp_path
