"""The pre-flight registry must name every tool the workflow invokes, and no other.

The scan reads the scripts and the library themselves rather than a hand-maintained list,
so a tool added to a script without being declared fails here.
"""
import pathlib
import re

from plasmidann.tools import (
    REQUIRED_TOOLS,
    grammar_problem,
    macsyfinder_version,
    model_grammars,
)

ROOT = pathlib.Path(__file__).resolve().parents[1]
SOURCES = [ROOT / "workflow" / "scripts", ROOT / "src" / "plasmidann"]
TOOL_NAMES = {t["name"] for t in REQUIRED_TOOLS}

# The executable is the first token of a command given to subprocess directly - a string
# (f"mmseqs ...") or a list (["diamond", ...]) - or assigned first to `cmd = (f"...`.
INVOCATION = re.compile(r"""(?:subprocess\.\w+\(|\bcmd\s*=\s*\(?)\s*\[?f?["']([A-Za-z_][\w.-]*)""")

# Tools no script names because something else shells out to them: IntegronFinder runs
# prodigal and cmsearch.
INDIRECT = {"prodigal", "cmsearch"}


def invoked_executables():
    found = {}
    for source in SOURCES:
        for path in sorted(source.glob("*.py")):
            for tool in INVOCATION.findall(path.read_text()):
                found.setdefault(tool, set()).add(path.name)
    return found


def test_the_scanner_sees_every_form_of_invocation():
    """Positive control: a string command, a `cmd` variable and a list-form call."""
    invoked = invoked_executables()
    assert "tier_search.py" in invoked.get("hmmsearch", ()), "cmd = (f\"hmmsearch ..."
    assert "preflight.py" in invoked.get("diamond", ()), 'subprocess.run(["diamond", ...'
    assert "protein_clustering.py" in invoked.get("mmseqs", ()), 'subprocess.run(f"mmseqs'


def test_every_invoked_executable_is_in_the_registry():
    undeclared = {t: sorted(f) for t, f in invoked_executables().items()
                  if t not in TOOL_NAMES}
    assert not undeclared, (
        f"tools invoked but absent from plasmidann.tools.REQUIRED_TOOLS, so pre-flight "
        f"cannot check them: {undeclared}")


def test_the_registry_names_no_tool_the_workflow_does_not_invoke():
    """A stale entry halts the run on a missing tool that no stage needs."""
    stale = sorted(TOOL_NAMES - set(invoked_executables()) - INDIRECT)
    assert not stale, (
        f"pre-flighted but invoked nowhere, and not declared as an indirect dependency: "
        f"{stale}")


def test_the_registry_names_the_indirect_dependencies():
    """IntegronFinder shells out to prodigal and cmsearch, which no script names."""
    assert INDIRECT <= TOOL_NAMES


def test_every_registry_entry_states_the_rules_that_need_it():
    for entry in REQUIRED_TOOLS:
        assert entry["stage"], f"{entry['name']} does not say which rule needs it"
        assert entry["why"], f"{entry['name']} does not say what breaks without it"


# --- MacSyFinder model grammar (pre-flight's CONJScan check) --------------------------

def test_the_macsyfinder_version_is_read_whatever_its_capitalisation():
    assert macsyfinder_version("MacSyFinder 2.1.6 \nusing:\n- Python 3.12") == (2, 1, 6)
    assert macsyfinder_version("Macsyfinder 2.1.4 \nusing:") == (2, 1, 4)
    assert macsyfinder_version("command not found") is None


def test_grammar_2_1_needs_macsyfinder_2_1_6():
    """The exact failure measured: CONJScan 2.1.0 under MacSyFinder 2.1.4."""
    assert "MacSyFinder >= 2.1.6" in grammar_problem({"2.1"}, (2, 1, 4))
    assert grammar_problem({"2.1"}, (2, 1, 6)) == ""
    # The DefenseFinder models (grammar 2.0) run under both.
    assert grammar_problem({"2.0"}, (2, 1, 4)) == ""
    assert "not known" in grammar_problem({"3.0"}, (2, 1, 6))


def test_model_grammars_are_read_from_the_definitions(tmp_path):
    d = tmp_path / "CONJScan" / "definitions" / "Plasmids"
    d.mkdir(parents=True)
    (d / "MOB.xml").write_text('<model min_genes_required="1" vers="2.1">\n</model>\n')
    (tmp_path / "CONJScan" / "profiles").mkdir()
    (tmp_path / "CONJScan" / "profiles" / "x.xml").write_text('<model vers="9.9">')
    assert model_grammars(tmp_path) == {"2.1"}
