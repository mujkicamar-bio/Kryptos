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

from plasmidann.tools import REQUIRED_TOOLS, tool_names

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


# ------------------------------------------------------------------------------------
# The PHROGs profile database must report a PHROG identifier.
# ------------------------------------------------------------------------------------

def test_the_phrogs_name_line_carries_the_phrog_id():
    """Spec section 18 lists phrog_id as a required field of a PHROGs hit.

    The first build did not produce one. `mmseqs convertprofiledb` takes each profile's
    header from its NAME line, and the PHROGs NAME line names the seed protein
    ('p428256 VI_04636'), not the PHROG. It also does not assign database keys in ffindex
    order - key 0's header is phrog_2402's NAME, not phrog_1's - so the identifier could
    not be recovered from the key either. A hit would have carried a phage protein name
    with no traceable PHROG, which is not the field the spec asks for and cannot be joined
    to the PHROGs annotation table when that becomes available.

    Prefixing the id onto the NAME line puts it where mmseqs will report it, as the first
    whitespace-delimited token of the target header.
    """
    import importlib.util
    import pathlib

    spec_ = importlib.util.spec_from_file_location(
        "build_phrogs_mmseqs",
        pathlib.Path(__file__).resolve().parents[1] / "tools" / "build_phrogs_mmseqs.py")
    module = importlib.util.module_from_spec(spec_)
    spec_.loader.exec_module(module)

    profile = ("HH-suite HHM file\nNAME  p428256 VI_04636\nFAM   \n"
               "LENG  120 match states\n")

    rewritten = module.label_with_phrog_id(profile, "phrog_1")

    assert "NAME  phrog_1 p428256 VI_04636\n" in rewritten, (
        f"the PHROG id is not the first token of the NAME line:\n{rewritten}")
    # Everything else must survive: the profile itself is what is being searched.
    assert "LENG  120 match states" in rewritten
    assert rewritten.count("NAME") == 1, "the original NAME line was duplicated"


def test_a_profile_already_carrying_its_id_is_not_labelled_twice():
    """The build is re-run whenever PHROGs is updated, and a double prefix would change
    the identifier a hit reports without anything failing."""
    import importlib.util
    import pathlib

    spec_ = importlib.util.spec_from_file_location(
        "build_phrogs_mmseqs",
        pathlib.Path(__file__).resolve().parents[1] / "tools" / "build_phrogs_mmseqs.py")
    module = importlib.util.module_from_spec(spec_)
    spec_.loader.exec_module(module)

    once = module.label_with_phrog_id("NAME  p428256 VI_04636\n", "phrog_1")
    twice = module.label_with_phrog_id(once, "phrog_1")

    assert once == twice
