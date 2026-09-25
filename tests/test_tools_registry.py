"""The pre-flight registry must name every tool the workflow actually invokes.

v1 lost a 45-hour job twice to a tool that was missing from PATH. The pre-flight rule was
written to prevent exactly that, and then covered only the four tools the cascade uses -
so a missing mafft, mmseqs, macsyfinder, integron_finder, prodigal, cmsearch or foldseek
still killed the run days in, which is the failure class pre-flight exists to remove.

This test reads the scripts themselves rather than a hand-maintained list, so adding a new
tool to a script without declaring it fails here in milliseconds.
"""
import pathlib
import re

from plasmidann.tools import (REQUIRED_TOOLS, grammar_problem, macsyfinder_version,
                              model_grammars, tool_names)

SCRIPTS = pathlib.Path(__file__).resolve().parents[1] / "workflow" / "scripts"

# subprocess.run(f"<tool> ...) - the executable is the first token of the command string.
INVOCATION = re.compile(r"subprocess\.run\(\s*\n?\s*f?\"([A-Za-z_][\w.-]*)", re.M)


def invoked_executables():
    found = {}
    for path in sorted(SCRIPTS.glob("*.py")):
        for tool in INVOCATION.findall(path.read_text()):
            found.setdefault(tool, set()).add(path.name)
    return found


def test_every_invoked_executable_is_in_the_registry():
    invoked = invoked_executables()
    assert invoked, "the scanner found no subprocess calls at all - the regex is stale"
    undeclared = {t: sorted(f) for t, f in invoked.items() if t not in tool_names()}
    assert not undeclared, (
        f"tools invoked by scripts but absent from plasmidann.tools.REQUIRED_TOOLS, so "
        f"pre-flight cannot check them: {undeclared}")


# Tools no script names because something else shells out to them.
INDIRECT = {"prodigal", "cmsearch"}


def test_the_registry_names_no_tool_the_workflow_stopped_using():
    """A stale registry entry is not harmless: pre-flight halts the run on a missing tool
    that no stage needs, which is the same lost day as the failure pre-flight exists to
    prevent, arriving from the opposite direction.

    This fired when S7d was removed - FastTree and hyphy stayed in the registry after the
    script that invoked them was deleted."""
    # Not the subprocess scanner: tier_search builds its command with the tool name
    # inside an f-string that the scanner's leading-token regex cannot see. Any mention
    # of the name in a script counts as use, which is loose enough to keep diamond and
    # strict enough to catch a tool whose only caller was deleted.
    mentioned = set()
    for path in SCRIPTS.glob("*.py"):
        text = path.read_text()
        mentioned.update(n for n in tool_names() if re.search(rf"\b{re.escape(n)}\b", text))
    stale = sorted(n for n in tool_names() if n not in mentioned and n not in INDIRECT)
    assert not stale, (
        f"pre-flighted but invoked by no script, and not declared as an indirect "
        f"dependency: {stale}")


def test_the_registry_names_the_indirect_dependencies():
    """IntegronFinder shells out to prodigal and cmsearch. Neither appears in any script,
    and their absence is invisible until IntegronFinder itself fails - which it did."""
    names = tool_names()
    for tool in ("prodigal", "cmsearch"):
        assert tool in names, f"{tool} is an indirect dependency and must be pre-flighted"


def test_every_registry_entry_states_the_stage_that_needs_it():
    for entry in REQUIRED_TOOLS:
        assert entry["stage"], f"{entry['name']} does not say which stage needs it"
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
