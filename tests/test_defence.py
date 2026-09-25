"""S8a: splitting DefenseFinder's two phases across the representation each one needs.

DefenseFinder runs in two phases. Phase 1 is an HMM search over 1,887 profiles: per
protein, order-independent. Phase 2 is MacSyFinder system calling over 711 model
definitions, each carrying a co-localisation constraint - a real one reads

    <model inter_gene_max_space="3" min_mandatory_genes_required="2" ...>

so it is entirely about gene adjacency.

v2 ran the whole thing on unique_proteins.faa: dereplicated and ordered by SHA-256 hash,
with the default --db-type ordered_replicon. Phase 1 was fine. Phase 2 was told that a
cryptographic hash ordering was genomic order. Demonstrated in review: the same two
proteins adjacent gave 1 system, separated by 5 decoys gave 0, with identical HMM hits.

The fix searches the dereplicated set once, propagates component labels to every ORF that
shares the sequence, and calls systems on ordered per-plasmid gene lists.
"""
from plasmidann.defence import candidate_plasmids, gembase_id, order_orfs, propagate_components

ORF_TO_SEQ = {
    "p1|1": "aaa", "p1|2": "bbb", "p1|3": "ccc",
    "p2|1": "aaa", "p2|2": "zzz",
    "p3|1": "yyy",
}


def test_a_component_hit_propagates_to_every_orf_sharing_the_sequence():
    """This is the whole point of dereplicating first: a protein identical on forty
    plasmids is searched once and labelled forty times."""
    hits = {"aaa": "RM_type_II"}

    labels = propagate_components(hits, ORF_TO_SEQ)

    assert labels["p1|1"] == "RM_type_II"
    assert labels["p2|1"] == "RM_type_II"
    assert "p1|2" not in labels


def test_propagation_does_not_invent_labels_for_unhit_sequences():
    assert propagate_components({}, ORF_TO_SEQ) == {}


# --- pruning: only plasmids that could possibly carry a system ------------------------

def test_only_plasmids_carrying_a_component_become_candidates():
    """System calling needs ordered proteins, which is the expensive representation. A
    plasmid with no component hit cannot produce a system, so it never needs to be
    written out in genomic order."""
    labels = {"p1|1": "RM_type_II", "p2|1": "RM_type_II"}

    assert candidate_plasmids(labels) == {"p1", "p2"}


def test_a_plasmid_with_no_components_is_pruned():
    assert "p3" not in candidate_plasmids({"p1|1": "RM_type_II"})


# --- ordering: the thing the hash ordering destroyed ---------------------------------

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
    """A gene reconstructed across the cut runs start..length then 1..end, so start > end.
    Sorting naively on start would place it last, when on the circle it is adjacent to the
    first gene. 94% of these plasmids are circular, so this is the common case, not an
    edge case."""
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


# ------------------------------------------------------------------------------------
# Phase 1 reads MacSyFinder's own all_systems.tsv.
#
# `--db-type unordered` is the mode phase 1 needs - report components, do not call systems
# - and MacSyFinder does not write best_solution.tsv in it. Two consequences followed, and
# both were invisible in the output:
#
#   1. defense-finder's post-treatment step opens best_solution.tsv unconditionally and
#      raises FileNotFoundError, so a search that had just found systems in all three
#      model families exited non-zero.
#   2. phase 1 then looked for *defense_finder_genes.tsv, which is a POST-TREATMENT
#      output and therefore never existed. Even without the crash it would have written
#      an empty table, and an empty defence table reads as "this collection has no
#      defence systems" rather than as "the parser found no file".
# ------------------------------------------------------------------------------------

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
    """The columns phase 1 needs are all in MacSyFinder's own output: hit_id names the
    protein, gene_name the component, model_fqn the model it belongs to, hit_i_eval the
    significance. model_fqn in particular is the tool's cited model identity rather than
    the guess the previous parser recorded as `model_unverified`."""
    from plasmidann.defence import parse_all_systems

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
    from plasmidann.defence import parse_all_systems

    path = tmp_path / "all_systems.tsv"
    path.write_text(ALL_SYSTEMS)

    rows = parse_all_systems([path])

    assert all(r["seq_id"] and not r["seq_id"].startswith("#") for r in rows)


def test_a_model_family_that_found_nothing_contributes_no_rows(tmp_path):
    """MacSyFinder writes an all_systems.tsv with only its comment header for a family
    that matched nothing. That is a real and common outcome - it is what crashed
    defense-finder's post-treatment - and it must read as zero rows, not as an error."""
    from plasmidann.defence import parse_all_systems

    empty = tmp_path / "empty_all_systems.tsv"
    empty.write_text("# macsyfinder 2.1.4 \n# models : CasFinder-3.1.0\n"
                     "# No Systems found\n")
    full = tmp_path / "all_systems.tsv"
    full.write_text(ALL_SYSTEMS)

    assert parse_all_systems([empty]) == []
    assert len(parse_all_systems([empty, full])) == 2


# ------------------------------------------------------------------------------------
# Phase 2 is ONE MacSyFinder database, whatever the core count.
#
# HMMER's independent e-value is the score's p-value times the number of sequences in the
# database searched, and MacSyFinder keeps a hit only below --i-evalue-sel (0.001 by
# default). Phase 2 once split its input into one chunk of replicons per core: a smaller
# database gives a smaller i-evalue, so the calls depended on -c. Measured with CONJScan on
# the test set (leaf 1.3): 8 chunks called 214 ORFs in 56 systems, one database 212 in 55.
# ------------------------------------------------------------------------------------

import os  # noqa: E402
import sys  # noqa: E402

from conftest import FakeSnakemake, read_tsv, run_script, write_fasta, write_tsv  # noqa: E402

# Stand-in for macsyfinder that reproduces the dependence under test: a planned hit is
# kept when p-value x (sequences in --sequence-db) < 0.001, as HMMER and MacSyFinder do.
FAKE_MACSYFINDER = r'''#!{python}
import csv, pathlib, sys
args = sys.argv[1:]
opt = lambda k: args[args.index(k) + 1]
out = pathlib.Path(opt("--out-dir"))
out.mkdir(parents=True)
(out / "argv.txt").write_text("\n".join(args))
ids = [l[1:].split()[0] for l in open(opt("--sequence-db")) if l.startswith(">")]
plan = {r["hit_id"]: r for r in csv.DictReader(open({plan!r}), delimiter="\t")}
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
    faa, gmap = base / "cand.faa", base / "map.tsv"
    write_fasta(faa, records)
    write_tsv(gmap, ["gembase_id", "orf_id", "plasmid_id"], mapping)
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
        input={"faa": str(faa), "map": str(gmap)}, output={"tsv": str(out)},
        params={"models_dir": str(models), "required": True}, threads=threads))
    return out


def test_phase2_calls_do_not_depend_on_the_core_count(tmp_path, monkeypatch):
    one = read_tsv(_phase2(tmp_path, monkeypatch, threads=1))
    many_path = _phase2(tmp_path, monkeypatch, threads=4)
    many = read_tsv(many_path)

    assert one == many, "the defence calls changed with the number of threads"
    assert [r["orf_id"] for r in one] == ["p0|3"]

    # One MacSyFinder process over every candidate replicon, the cores to --worker.
    runs = [p for p in (many_path.parent / "phase2").iterdir() if p.is_dir()]
    assert [p.name for p in runs] == ["run"]
    argv = (runs[0] / "argv.txt").read_text().split("\n")
    assert argv[argv.index("--worker") + 1] == "4"
    assert argv[argv.index("--replicon-topology") + 1] == "circular"
