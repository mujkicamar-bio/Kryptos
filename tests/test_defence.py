"""DefenseFinder's two phases on the representation each needs.

Phase 1 (per-protein HMM search) runs on dereplicated proteins; phase 2 (MacSyFinder system
calling, whose models count intervening genes) runs on every ORF of a candidate plasmid in
genomic order, written as one gembase database.
"""
import os
import sys

import pytest
from conftest import FakeSnakemake, read_tsv, run_script, write_fasta, write_tsv

from plasmidann.defence import (
    candidate_plasmids,
    gembase_id,
    gembase_records,
    order_orfs,
    parse_all_systems,
    propagate_components,
)

ORF_TO_SEQ = {
    "p1|1": "aaa", "p1|2": "bbb", "p1|3": "ccc",
    "p2|1": "aaa", "p2|2": "zzz",
    "p3|1": "yyy",
}


def test_a_component_hit_propagates_to_every_orf_sharing_the_sequence():
    """A protein identical on several plasmids is searched once and applies to every copy."""
    assert propagate_components({"aaa"}, ORF_TO_SEQ) == {"p1|1", "p2|1"}


def test_propagation_does_not_invent_labels_for_unhit_sequences():
    assert propagate_components(set(), ORF_TO_SEQ) == set()


# --- pruning: only plasmids that could possibly carry a system ------------------------

def test_only_plasmids_carrying_a_component_become_candidates():
    """System calling needs ordered proteins, which is the expensive representation. A
    plasmid with no component hit cannot produce a system, so it never needs to be
    written out in genomic order."""
    assert candidate_plasmids({"p1|1", "p2|1"}) == {"p1", "p2"}


def test_a_plasmid_with_no_components_is_pruned():
    assert "p3" not in candidate_plasmids({"p1|1"})


# --- ordering ------------------------------------------------------------------------

def test_orfs_are_returned_in_genomic_order():
    """MacSyFinder counts intervening genes. If the order is wrong, inter_gene_max_space
    is measuring nothing."""
    orfs = [
        {"orf_id": "p1|3", "start": 900, "end": 1200, "spans_origin": "0"},
        {"orf_id": "p1|1", "start": 10, "end": 300, "spans_origin": "0"},
        {"orf_id": "p1|2", "start": 400, "end": 800, "spans_origin": "0"},
    ]

    assert [o["orf_id"] for o in order_orfs(orfs)] == ["p1|1", "p1|2", "p1|3"]


def test_an_origin_spanning_gene_sorts_to_the_start_of_the_replicon():
    """A gene reconstructed across the cut runs start..length then 1..end, so start > end;
    it is placed first, so that positions follow the molecule from its origin."""
    orfs = [
        {"orf_id": "p1|2", "start": 400, "end": 800, "spans_origin": "0"},
        {"orf_id": "p1|1", "start": 4800, "end": 200, "spans_origin": "1"},
        {"orf_id": "p1|3", "start": 900, "end": 1200, "spans_origin": "0"},
    ]

    assert [o["orf_id"] for o in order_orfs(orfs)] == ["p1|1", "p1|2", "p1|3"]


# --- gembase naming ------------------------------------------------------------------

def test_gembase_ids_encode_replicon_and_position():
    """MacSyFinder's gembase mode splits an id on the last underscore to recover the
    replicon, so one run can hold many plasmids and still treat each separately."""
    assert gembase_id("COMPASS_AB007909.1", 7) == "COMPASS-AB007909.1_00007"


def test_gembase_ids_are_sortable_in_genomic_order():
    """Zero padding matters: position 2 must sort before position 10."""
    ids = [gembase_id("p1", i) for i in (10, 2, 1)]
    assert sorted(ids) == [gembase_id("p1", 1), gembase_id("p1", 2), gembase_id("p1", 10)]


def test_underscores_in_a_plasmid_id_do_not_break_the_replicon_split():
    """Plasmid ids here contain underscores (COMPASS_AB007909.1). MacSyFinder splits on
    the LAST underscore, so an unescaped id would make every plasmid its own malformed
    replicon name."""
    gid = gembase_id("COMPASS_AB007909.1", 1)
    replicon, position = gid.rsplit("_", 1)
    assert replicon == "COMPASS-AB007909.1"
    assert position == "00001"


# --- phase 1 reads MacSyFinder's own all_systems.tsv ------------------------------------
# MacSyFinder writes no best_solution.tsv under --db-type unordered.

ALL_SYSTEMS = """\
# macsyfinder 2.1.4
# models : defense-finder-models-3.1.0
# defense-finder run --db-type unordered
# Likely Systems found:

replicon\thit_id\tgene_name\thit_pos\tmodel_fqn\tsys_id\tsys_wholeness\thit_i_eval
unique_proteins\t2171fc6c\tClover__CloA\t55\tdefense-finder-models/DefenseFinder/Clover/Clover\tunique_proteins_Clover_41\t1.000\t3.3e-188
unique_proteins\t284ca90c\tClover__CloB\t54\tdefense-finder-models/DefenseFinder/Clover/Clover\tunique_proteins_Clover_41\t1.000\t3.6e-82
"""


def test_component_hits_are_read_from_all_systems(tmp_path):
    """hit_id names the protein, gene_name the component, model_fqn the model and
    hit_i_eval the significance."""
    path = tmp_path / "all_systems.tsv"
    path.write_text(ALL_SYSTEMS)

    rows = parse_all_systems([path])

    assert len(rows) == 2
    assert rows[0]["seq_id"] == "2171fc6c"
    assert rows[0]["component"] == "Clover__CloA"
    assert rows[0]["model"] == "defense-finder-models/DefenseFinder/Clover/Clover"
    assert rows[0]["hit_evalue"] == "3.3e-188"


def test_the_macsyfinder_comment_header_is_not_parsed_as_data(tmp_path):
    """all_systems.tsv opens with four '#' lines and a blank line before its header. A
    parser that took the first line as the header would read every row as one field and
    silently produce nothing."""

    path = tmp_path / "all_systems.tsv"
    path.write_text(ALL_SYSTEMS)

    rows = parse_all_systems([path])

    assert all(r["seq_id"] and not r["seq_id"].startswith("#") for r in rows)


def test_a_model_family_that_found_nothing_contributes_no_rows(tmp_path):
    """MacSyFinder writes an all_systems.tsv with only its comment header for a model
    family that matched nothing; it reads as zero rows, not as an error."""

    empty = tmp_path / "empty_all_systems.tsv"
    empty.write_text("# macsyfinder 2.1.4 \n# models : CasFinder-3.1.0\n"
                     "# No Systems found\n")
    full = tmp_path / "all_systems.tsv"
    full.write_text(ALL_SYSTEMS)

    assert parse_all_systems([empty]) == []
    assert len(parse_all_systems([empty, full])) == 2


# --- phase 2 is one MacSyFinder database, whatever the core count ----------------------
# HMMER's independent e-value is the p-value times the number of sequences searched, and
# MacSyFinder keeps a hit only below --i-evalue-sel (0.001), so per-core chunks would make
# the calls depend on the core count.

# Stand-in for macsyfinder that reproduces the dependence under test: a planned hit is
# kept when p-value x (sequences in --sequence-db) < 0.001, as HMMER and MacSyFinder do.
# The planned hits are DefenseFinder systems, so only that family's run reports them.
FAKE_MACSYFINDER = r'''#!{python}
import csv, pathlib, sys
args = sys.argv[1:]
opt = lambda k: args[args.index(k) + 1]
out = pathlib.Path(opt("--out-dir"))
out.mkdir(parents=True)
(out / "argv.txt").write_text("\n".join(args))
ids = [l[1:].split()[0] for l in open(opt("--sequence-db")) if l.startswith(">")]
plan = {r["hit_id"]: r for r in csv.DictReader(open({plan!r}), delimiter="\t")}
if opt("--models") != "defense-finder-models/DefenseFinder":
    plan = {}
header = ["replicon", "hit_id", "gene_name", "hit_pos", "model_fqn", "sys_id",
          "sys_wholeness", "sys_score", "hit_gene_ref", "hit_status", "hit_i_eval",
          "hit_profile_cov"]
with open(out / "best_solution.tsv", "w") as fh:
    fh.write("# macsyfinder 2.1.4\n# Systems found:\n\n" + "\t".join(header) + "\n")
    for gid in ids:
        p = plan.get(gid)
        if p is None or float(p["pvalue"]) * len(ids) >= 0.001:
            continue
        rep = gid.rsplit("_", 1)[0]
        row = dict(replicon=rep, hit_id=gid, gene_name=p["component"],
                   model_fqn="defense-finder-models/DefenseFinder/" + p["system"],
                   sys_id=rep + "_" + p["system"] + "_1", sys_wholeness="1.000",
                   hit_status="mandatory", hit_i_eval=str(float(p["pvalue"]) * len(ids)))
        fh.write("\t".join(row.get(k, "") for k in header) + "\n")
'''


def _phase2(tmp_path, monkeypatch, threads):
    """Run defence_systems.py over eight replicons of 25 genes with the stand-in."""
    base = tmp_path / f"t{threads}"
    records, mapping = [], []
    for r in range(8):
        for pos in range(1, 26):
            gid = gembase_id(f"p{r}", pos)
            records.append((gid, "MKV"))
            mapping.append([gid, f"p{r}|{pos}", f"p{r}"])
    # 200 sequences in one database: 4e-6 x 200 = 8e-4 is kept, 1e-5 x 200 = 2e-3 is not.
    # A chunk of 50 sequences would keep both, which is the core-count dependence.
    plan = [[gembase_id("p0", 3), "Gabija", "GajA", 4e-6],
            [gembase_id("p5", 7), "Zorya", "ZorA", 1e-5]]
    faa, gmap, master = base / "cand.faa", base / "map.tsv", base / "master.tsv"
    write_fasta(faa, records)
    write_tsv(gmap, ["gembase_id", "orf_id", "plasmid_id"], mapping)
    # p1 is linear; 'direct terminal repeat' is a circular molecule (darkorf.circular).
    write_tsv(master, ["plasmid_id", "topology"],
              [["p0", "circular"], ["p1", "linear"], ["p2", "direct terminal repeat"]]
              + [[f"p{r}", "circular"] for r in range(3, 8)])
    write_tsv(base / "plan.tsv", ["hit_id", "system", "component", "pvalue"], plan)
    exe = base / "bin" / "macsyfinder"
    exe.parent.mkdir(parents=True)
    exe.write_text(FAKE_MACSYFINDER.replace("{python}", sys.executable)
                   .replace("{plan!r}", repr(str(base / "plan.tsv"))))
    exe.chmod(0o755)
    monkeypatch.setenv("PATH", f"{exe.parent}:{os.environ['PATH']}")
    models = base / "models"
    models.mkdir()
    (models / "defense-finder-models").mkdir()
    out = base / "out" / "defence_systems.tsv"
    out.parent.mkdir(parents=True)
    run_script("defence_systems.py", FakeSnakemake(
        input={"faa": str(faa), "map": str(gmap), "master": str(master)},
        output={"tsv": str(out)},
        params={"models_dir": str(models), "required": True}, threads=threads))
    return out


def test_phase2_calls_do_not_depend_on_the_core_count(tmp_path, monkeypatch):
    one = read_tsv(_phase2(tmp_path, monkeypatch, threads=1))
    many_path = _phase2(tmp_path, monkeypatch, threads=4)
    many = read_tsv(many_path)

    assert one == many, "the defence calls changed with the number of threads"
    assert [r["orf_id"] for r in one] == ["p0|3"]



def test_phase2_searches_the_model_families_of_defense_finder(tmp_path, monkeypatch):
    """`defense-finder run` (DefenseFinder 3.0.0) runs one MacSyFinder process per family:
    DefenseFinder and RM with --coverage-profile 0.4 and --exchangeable-weight 1, CasFinder
    without options (its package configuration applies), and the AntiDefenseFinder models
    only on request. Each process covers every candidate replicon, each with its registry
    topology, with all cores."""
    phase2 = _phase2(tmp_path, monkeypatch, threads=4).parent / "phase2"
    topology = str(phase2 / "topology.txt")
    # The format MacSyFinder 2.1.4 parses: "<replicon>: <topology>" (macsypy.database).
    assert (phase2 / "topology.txt").read_text().splitlines() == [
        "p0: circular", "p1: linear", "p2: circular", "p3: circular", "p4: circular",
        "p5: circular", "p6: circular", "p7: circular"]
    calls = {}
    for run in (d for d in phase2.iterdir() if d.is_dir()):
        argv = (run / "argv.txt").read_text().split("\n")
        opt = {k: argv[argv.index(k) + 1] if k in argv else None
               for k in ("--models", "--coverage-profile", "--exchangeable-weight",
                         "--worker", "--topology-file")}
        assert argv[argv.index("--models") + 2] == "all"
        calls[run.name] = opt
    assert calls == {
        "DefenseFinder": {"--models": "defense-finder-models/DefenseFinder",
                          "--coverage-profile": "0.4", "--exchangeable-weight": "1",
                          "--worker": "4", "--topology-file": topology},
        "RM": {"--models": "defense-finder-models/RM",
               "--coverage-profile": "0.4", "--exchangeable-weight": "1",
               "--worker": "4", "--topology-file": topology},
        "Cas": {"--models": "CasFinder",
                "--coverage-profile": None, "--exchangeable-weight": None,
                "--worker": "4", "--topology-file": topology},
    }


# --- gembase records ------------------------------------------------------------------

def _index_row(orf_id, start, spans="0"):
    return {"orf_id": orf_id, "plasmid_id": orf_id.rsplit("|", 1)[0], "start": start,
            "spans_origin": spans, "seq": "MKV"}


def test_gembase_records_number_each_plasmid_in_genomic_order():
    rows = [_index_row("p_1|1", 500), _index_row("p_1|2", 100), _index_row("p2|1", 10)]
    assert [(g, o["orf_id"]) for g, _, o in gembase_records(rows)] == [
        ("p-1_00001", "p_1|2"), ("p-1_00002", "p_1|1"), ("p2_00001", "p2|1")]
    assert [o["orf_id"] for _, _, o in gembase_records(rows, keep={"p2"})] == ["p2|1"]


def test_two_plasmids_with_one_replicon_name_are_refused():
    """p_1 and p-1 are both replicon p-1; merging them would join two molecules' ORFs."""
    with pytest.raises(ValueError, match="occurs twice"):
        list(gembase_records([_index_row("p_1|1", 1), _index_row("p-1|1", 1)]))


def test_an_index_not_grouped_by_plasmid_is_refused():
    with pytest.raises(ValueError, match="occurs twice"):
        list(gembase_records([_index_row("a|1", 1), _index_row("b|1", 1),
                              _index_row("a|2", 9)]))


# --- phase 1 needs the tables of every model family ------------------------------------

# Stand-in for defense-finder: writes all_systems.tsv for the families in FAMILIES, where
# MacSyFinder's raw output lies under --preserve-raw, and exits 1 as the real wrapper does
# after an unordered search.
FAKE_DEFENSE_FINDER = r'''#!{python}
import pathlib, sys
args = sys.argv[1:]
raw = pathlib.Path(args[args.index("--out-dir") + 1]) / "defense-finder-tmp"
for family in {families!r}:
    (raw / family).mkdir(parents=True)
    (raw / family / "all_systems.tsv").write_text({table!r})
sys.exit(1)
'''


def _phase1(tmp_path, monkeypatch, families):
    exe = tmp_path / "bin" / "defense-finder"
    exe.parent.mkdir()
    exe.write_text(FAKE_DEFENSE_FINDER.replace("{python}", sys.executable)
                   .replace("{families!r}", repr(families))
                   .replace("{table!r}", repr(ALL_SYSTEMS)))
    exe.chmod(0o755)
    monkeypatch.setenv("PATH", f"{exe.parent}:{os.environ['PATH']}")
    models = tmp_path / "models"
    (models / "defense-finder-models").mkdir(parents=True)
    faa, out = tmp_path / "unique.faa", tmp_path / "defence_components.tsv"
    write_fasta(faa, [("2171fc6c", "MKV")])
    run_script("defence_search.py", FakeSnakemake(
        input={"faa": str(faa)}, output={"tsv": str(out)},
        params={"models_dir": str(models), "required": True}, threads=1))
    return out


def test_phase1_reads_every_family_despite_the_wrappers_exit_code(tmp_path, monkeypatch):
    out = _phase1(tmp_path, monkeypatch, ["DefenseFinder", "RM", "Cas"])
    rows = read_tsv(out)
    assert len(rows) == 6
    assert {r["status"] for r in rows} == {"SUCCESS"}


def test_phase1_fails_when_a_model_family_wrote_no_table(tmp_path, monkeypatch):
    """defense-finder runs the families in turn and stops at the first failure, so a crash
    in the CasFinder run leaves the DefenseFinder and RM tables behind."""
    with pytest.raises(SystemExit, match="Cas"):
        _phase1(tmp_path, monkeypatch, ["DefenseFinder", "RM"])
