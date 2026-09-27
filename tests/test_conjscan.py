"""CONJScan conjugation systems: parser, plasmid mobility class, and the stage contract.

The stage runs MacSyFinder from a configured path, so the script tests drive it with a
stand-in executable that writes the same best_solution.tsv MacSyFinder does. What is under
test is everything around the tool: the gembase of every plasmid as one database, the
mapping back to orf ids, the version stamp, the class per plasmid, and the NOT_RUN / halt
contract shared with the defence stage.
"""
import subprocess
import sys

import pytest
from conftest import FakeSnakemake, read_tsv, run_script, write_tsv

from plasmidann.conjscan import (
    CLASS_COLUMNS,
    COLUMNS,
    installed_version,
    plasmid_class,
    read_best_solution,
)

BEST_HEADER = ["replicon", "hit_id", "gene_name", "hit_pos", "model_fqn", "sys_id",
               "sys_loci", "locus_num", "sys_wholeness", "sys_score", "sys_occ",
               "hit_gene_ref", "hit_status", "hit_seq_len", "hit_i_eval", "hit_score",
               "hit_profile_cov", "hit_seq_cov", "hit_begin_match", "hit_end_match",
               "counterpart", "used_in"]


# --- parser ---------------------------------------------------------------------------

def _best_solution(path, rows, comments=("# macsyfinder 2.1.6", "# Systems found:")):
    with open(path, "w") as fh:
        for c in comments:
            fh.write(c + "\n")
        fh.write("\n")
        if rows:
            fh.write("\t".join(BEST_HEADER) + "\n")
            for r in rows:
                fh.write("\t".join(r.get(k, "") for k in BEST_HEADER) + "\n")


def test_best_solution_rows_are_read_below_the_comment_header(tmp_path):
    path = tmp_path / "best_solution.tsv"
    _best_solution(path, [{"replicon": "p-1", "hit_id": "p-1_00003",
                           "gene_name": "T4SS_MOBP1",
                           "model_fqn": "CONJScan/Plasmids/MOB", "sys_id": "p-1_MOB_1",
                           "sys_wholeness": "1.000", "hit_status": "mandatory"}])
    rows = read_best_solution(path)
    assert len(rows) == 1
    assert rows[0]["hit_id"] == "p-1_00003"
    assert rows[0]["model_fqn"] == "CONJScan/Plasmids/MOB"


def test_a_run_that_found_no_system_contributes_no_rows(tmp_path):
    """MacSyFinder writes a file of comments alone when no replicon carries a system."""
    path = tmp_path / "best_solution.tsv"
    _best_solution(path, [], comments=("# macsyfinder 2.1.6", "# No Systems found"))
    assert read_best_solution(path) == []


# --- models version ---------------------------------------------------------------------

def _models(root, version="2.1.0"):
    meta = root / "CONJScan" / "metadata.yml"
    meta.parent.mkdir(parents=True)
    meta.write_text("short_desc: CONJScan models\n"
                    f"vers: {version}\n"
                    "license: CC BY-NC-SA 4.0\n")
    return root


def test_the_installed_version_is_read_from_the_package_metadata(tmp_path):
    assert installed_version(_models(tmp_path / "m", "2.1.0")) == "2.1.0"


def test_no_installed_version_when_the_models_are_absent(tmp_path):
    assert installed_version(tmp_path / "nothing") is None


# --- plasmid mobility class (Coluzzi et al. 2022) ---------------------------------------

@pytest.mark.parametrize("systems,expected", [
    ({"T4SS_typeF", "MOB"}, "pCONJ"),
    ({"T4SS_typeT"}, "pCONJ"),
    # A complete machinery outranks a decayed one on the same plasmid.
    ({"T4SS_typeF", "dCONJ_typeF"}, "pCONJ"),
    ({"dCONJ_typeG"}, "pdCONJ"),
    ({"dCONJ_typeB", "MOB"}, "pdCONJ"),
    ({"MOB"}, "pMOB"),
    (set(), "pMOBless"),
])
def test_plasmid_class_follows_the_most_complete_system(systems, expected):
    assert plasmid_class(systems) == expected


def test_an_unknown_system_type_is_an_error_not_a_silent_class():
    """A model renamed in a later CONJScan release must stop the stage, not fall through
    to pMOBless, which would be a biological claim the run has not earned."""
    with pytest.raises(ValueError):
        plasmid_class({"T6SS_typeX"})


# --- the stage script -------------------------------------------------------------------

FAKE_MSF = r'''#!{python}
"""Stand-in for macsyfinder: writes best_solution.tsv from a table of planned calls."""
import csv, pathlib, sys
args = sys.argv[1:]
opt = lambda k: args[args.index(k) + 1]
models = pathlib.Path(opt("--models-dir"))
out = pathlib.Path(opt("--out-dir"))
out.mkdir(parents=True)
(out / "argv.txt").write_text("\n".join(args))
ids = [l[1:].split()[0] for l in open(opt("--sequence-db")) if l.startswith(">")]
(out / "ids.txt").write_text("\n".join(ids))
plan = {r["hit_id"]: r for r in csv.DictReader(open(models / "planned_calls.tsv"),
                                                 delimiter="\t")}
header = {header}
with open(out / "best_solution.tsv", "w") as fh:
    fh.write("# macsyfinder 2.1.6\n# Systems found:\n\n" + "\t".join(header) + "\n")
    for gid in ids:
        if gid in plan:
            p = plan[gid]
            rep = gid.rsplit("_", 1)[0]
            model = "CONJScan/Plasmids/" + p["system"]
            row = dict(replicon=rep, hit_id=gid, gene_name=p["component"],
                       model_fqn=model, sys_id=rep + "_" + p["system"] + "_1",
                       sys_wholeness=p["wholeness"], hit_status=p["hit_status"])
            fh.write("\t".join(row.get(k, "") for k in header) + "\n")
'''


def _fake_exe(tmp_path, fail=False):
    exe = tmp_path / "bin" / "macsyfinder"
    exe.parent.mkdir(parents=True)
    body = (FAKE_MSF.replace("{python}", sys.executable)
            .replace("{header}", repr(BEST_HEADER)))
    if fail:
        body = f"#!{sys.executable}\nimport sys\nsys.exit(3)\n"
    exe.write_text(body)
    exe.chmod(0o755)
    return exe


def _orf_index(tmp_path):
    """Three plasmids: pA conjugative, p_B mobilisable (underscore in its id), pC none."""
    rows = []
    for plasmid, n in (("pA", 4), ("p_B", 3), ("pC", 2)):
        for i in range(1, n + 1):
            # Listed out of genomic order on purpose; the gembase must reorder them.
            start = (n + 1 - i) * 1000
            rows.append([f"{plasmid}|{i}", plasmid, start, start + 300, 1, 0, 0, 11,
                         "MKV"])
    path = tmp_path / "orf_index.tsv"
    write_tsv(path, ["orf_id", "plasmid_id", "start", "end", "strand", "partial",
                     "spans_origin", "translation_table", "seq"], rows)
    return path


def _planned(models):
    # Gembase positions are genomic: pA|4 has the smallest start, so it is pA_00001.
    write_tsv(models / "planned_calls.tsv",
              ["hit_id", "system", "component", "hit_status", "wholeness"],
              [["pA_00001", "T4SS_typeF", "T4SS_MOBF", "mandatory", "0.846"],
               ["pA_00002", "T4SS_typeF", "T4SS_virb4", "mandatory", "0.846"],
               ["pA_00003", "T4SS_typeF", "T4SS_F_traL", "accessory", "0.846"],
               ["p-B_00002", "MOB", "T4SS_MOBP1", "mandatory", "1.000"]])


def _run(tmp_path, *, models=None, exe=None, required=False, version="2.1.0", threads=1):
    out = tmp_path / "out" / "conjugation_systems.tsv"
    classes = tmp_path / "out" / "conjugation_plasmid_class.tsv"
    out.parent.mkdir(parents=True)  # Snakemake creates output directories before a job
    run_script("conjugation_systems.py", FakeSnakemake(
        input={"index": str(_orf_index(tmp_path))},
        output={"systems": str(out), "classes": str(classes)},
        params={"models_dir": str(models or tmp_path / "no_models"),
                "exe": str(exe or tmp_path / "no_exe" / "macsyfinder"),
                "required": required, "version": version},
        threads=threads))
    return out, classes


def test_calls_are_mapped_back_to_orfs_and_stamped_with_the_models_version(tmp_path):
    models = _models(tmp_path / "models")
    _planned(models)
    out, classes = _run(tmp_path, models=models, exe=_fake_exe(tmp_path))

    with open(out) as fh:
        assert fh.readline().rstrip("\n").split("\t") == COLUMNS
    rows = {r["orf_id"]: r for r in read_tsv(out)}
    assert set(rows) == {"pA|4", "pA|3", "pA|2", "p_B|2"}
    assert rows["pA|4"]["plasmid_id"] == "pA"
    assert rows["pA|4"]["system"] == "T4SS_typeF"
    assert rows["pA|4"]["component"] == "T4SS_MOBF"
    assert rows["pA|4"]["system_id"] == "pA_T4SS_typeF_1"
    assert rows["pA|2"]["hit_status"] == "accessory"
    assert rows["pA|4"]["sys_wholeness"] == "0.846"
    assert rows["p_B|2"]["plasmid_id"] == "p_B"
    assert {r["conjscan_version"] for r in rows.values()} == {"2.1.0"}
    assert {r["status"] for r in rows.values()} == {"SUCCESS"}

    with open(classes) as fh:
        assert fh.readline().rstrip("\n").split("\t") == CLASS_COLUMNS
    cls = {r["plasmid_id"]: r["class"] for r in read_tsv(classes)}
    # Every plasmid gets a class, including the one on which nothing was found.
    assert cls == {"pA": "pCONJ", "p_B": "pMOB", "pC": "pMOBless"}


def test_every_plasmid_is_searched_with_the_plasmid_models_and_circular_topology(tmp_path):
    """The defence gembase is pruned to plasmids with a defence component; CONJScan must
    see all of them, and every ORF of each, because MacSyFinder counts intervening genes."""
    models = _models(tmp_path / "models")
    _planned(models)
    out, _ = _run(tmp_path, models=models, exe=_fake_exe(tmp_path))

    run = out.parent / "conjscan" / "run"
    ids = (run / "ids.txt").read_text().split()
    assert len(ids) == 9
    argv = (run / "argv.txt").read_text().split("\n")
    assert argv[argv.index("--models") + 1:argv.index("--models") + 3] == \
        ["CONJScan/Plasmids", "all"]
    assert argv[argv.index("--db-type") + 1] == "gembase"
    assert argv[argv.index("--replicon-topology") + 1] == "circular"
    # The rebuilt gembase is not kept; MacSyFinder's own results are.
    assert not list((out.parent / "conjscan").glob("gembase.faa*"))


def test_one_database_is_searched_whatever_the_core_count(tmp_path):
    """HMMER's i-evalue scales with the number of sequences searched, so per-core chunks
    would make the calls depend on -c. The cores go to MacSyFinder's --worker instead."""
    models = _models(tmp_path / "models")
    _planned(models)
    out, _ = _run(tmp_path, models=models, exe=_fake_exe(tmp_path), threads=3)

    runs = [p for p in (out.parent / "conjscan").iterdir() if p.is_dir()]
    assert [p.name for p in runs] == ["run"]
    assert len((runs[0] / "ids.txt").read_text().split()) == 9
    argv = (runs[0] / "argv.txt").read_text().split("\n")
    assert argv[argv.index("--worker") + 1] == "3"


def _assert_not_run(out, classes):
    """One NOT_RUN row without an ORF, and no class: an empty systems table would read as
    'searched, no conjugation system', and pMOBless for every plasmid likewise."""
    rows = read_tsv(out)
    assert [(r["orf_id"], r["status"]) for r in rows] == [("", "NOT_RUN")]
    assert classes.read_text().splitlines() == ["\t".join(CLASS_COLUMNS)]


def test_absent_models_record_not_run_when_not_required(tmp_path):
    with pytest.raises(SystemExit) as done:
        _run(tmp_path, exe=_fake_exe(tmp_path), required=False)
    assert done.value.code in (0, None)
    _assert_not_run(tmp_path / "out" / "conjugation_systems.tsv",
                    tmp_path / "out" / "conjugation_plasmid_class.tsv")


def test_an_absent_executable_records_not_run_when_not_required(tmp_path):
    models = _models(tmp_path / "models")
    with pytest.raises(SystemExit) as done:
        _run(tmp_path, models=models, required=False)
    assert done.value.code in (0, None)
    _assert_not_run(tmp_path / "out" / "conjugation_systems.tsv",
                    tmp_path / "out" / "conjugation_plasmid_class.tsv")


def test_absent_models_halt_the_stage_when_required(tmp_path):
    with pytest.raises(SystemExit) as done:
        _run(tmp_path, exe=_fake_exe(tmp_path), required=True)
    assert "CONJScan models" in str(done.value.code)


def test_an_absent_executable_halts_the_stage_when_required(tmp_path):
    models = _models(tmp_path / "models")
    with pytest.raises(SystemExit) as done:
        _run(tmp_path, models=models, required=True)
    assert "executable" in str(done.value.code)


def test_a_models_version_other_than_the_configured_one_fails(tmp_path):
    """The version is what the methods cite; a silent upgrade would change the calls
    (2.1.0 added the MOBM profile) under the old citation."""
    models = _models(tmp_path / "models", version="2.0.1")
    _planned(models)
    with pytest.raises(SystemExit) as done:
        _run(tmp_path, models=models, exe=_fake_exe(tmp_path), required=False)
    assert "2.0.1" in str(done.value.code) and "2.1.0" in str(done.value.code)


def test_a_failed_macsyfinder_run_fails_the_stage(tmp_path):
    models = _models(tmp_path / "models")
    with pytest.raises(subprocess.CalledProcessError):
        _run(tmp_path, models=models, exe=_fake_exe(tmp_path, fail=True))
