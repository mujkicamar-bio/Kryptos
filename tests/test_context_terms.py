"""S8c context terms: what a family's neighbours are, counted per independent lineage.

Two neighbour rules, deliberately different. A gene label (amr:, ta:, metal:, conj_role:,
mge:, antidefence:) counts only from a neighbour in the same directon - same strand, gaps
at most context.max_operon_gap - within +-3 genes, the FESNov rule (Rodriguez del Rio et
al. 2024 Nature 626:377). A system (defence:, conj:) counts when the ORF is itself a
component or a component lies within +-3 on EITHER strand, because systems mix strands.
"""
import pytest
from conftest import FakeSnakemake, read_tsv, run_script, write_tsv

from plasmidann.context import directons, flanks
from plasmidann.context_terms import (
    COLUMNS,
    TOO_FEW_LINEAGES,
    family_term_rows,
    label_term,
    orf_term_sources,
    system_term,
    window_covers_plasmid,
)


def _g(oid, start, end, strand):
    return {"orf_id": oid, "start": start, "end": end, "strand": strand}


def _unit(genes, oid, circular=False, length=None):
    for unit in directons(genes, max_gap=100, circular=circular, length=length):
        if oid in unit:
            return unit
    raise AssertionError(oid)


# --- label kinds to terms ------------------------------------------------------------

@pytest.mark.parametrize("kind,prefix", [
    ("card_amr_family", "amr"), ("tadb_ta", "ta"), ("bacmet_compound", "metal"),
    ("oritdb_role", "conj_role"), ("mobileog_category", "mge"),
    ("dbapis_family", "antidefence"), ("acrdb_family", "antidefence")])
def test_every_contracted_label_kind_maps_to_its_prefix(kind, prefix):
    assert label_term(kind, "X1", "") == f"{prefix}:X1"


def test_amrfinder_terms_follow_the_element_type():
    """AMR elements are amr:, STRESS metal and biocide are metal:, and the other stress
    subtypes and virulence factors have no context term in the contract."""
    assert label_term("amrfinder_gene", "blaTEM-1", "AMR/AMR") == "amr:blaTEM-1"
    assert label_term("amrfinder_gene", "merA", "STRESS/METAL") == "metal:merA"
    assert label_term("amrfinder_gene", "qacE", "STRESS;BIOCIDE") == "metal:qacE"
    assert label_term("amrfinder_gene", "asr", "STRESS/ACID") is None
    assert label_term("amrfinder_gene", "iutA", "VIRULENCE/VIRULENCE") is None


def test_an_amrfinder_label_without_an_element_type_is_refused():
    """Without the element type the prefix is a guess; a silent guess would put a metal
    resistance gene under amr:."""
    with pytest.raises(ValueError):
        label_term("amrfinder_gene", "merA", "")


def test_kinds_outside_the_contract_give_no_term():
    """KEGG is not a context source, and amr: comes from card_amr_family, not pharokka."""
    for kind in ("kegg_ko", "gene_symbol", "card_gene_family", "pfam_family"):
        assert label_term(kind, "x", "") is None


def test_a_system_term_uses_the_model_name_not_its_path():
    assert system_term("defence", "defense-finder-models/DefenseFinder/Clover/Clover") == \
        "defence:Clover"
    assert system_term("conj", "T4SS_typeF") == "conj:T4SS_typeF"


# --- the two neighbour rules -----------------------------------------------------------

GENES = [_g("p|1", 100, 400, 1), _g("p|2", 430, 700, 1),       # directon with p|1
         _g("p|3", 1500, 1800, 1),                             # same strand, gap 800
         _g("p|4", 1850, 2100, -1),                            # opposite strand
         _g("p|5", 2200, 2500, 1)]


def _near(genes, focal, circular=False):
    left, right = flanks(genes, window=3, circular=circular)[focal]
    return left + right


def _sources(focal, label_terms=None, system_terms=None, genes=GENES):
    return orf_term_sources(focal, _near(genes, focal), directon=_unit(genes, focal),
                            label_terms=label_terms or {}, system_terms=system_terms or {})


def test_a_gene_label_counts_from_a_same_directon_neighbour():
    assert _sources("p|1", label_terms={"p|2": {"amr:blaTEM-1"}}) == [
        ("amr:blaTEM-1", "p|2")]


def test_a_gene_label_on_the_same_strand_but_another_directon_does_not_count():
    """p|3 is within +-3 and on the same strand, but 800 nt away: not co-transcribed."""
    assert _sources("p|1", label_terms={"p|3": {"amr:blaTEM-1"}}) == []


def test_a_gene_label_on_the_opposite_strand_does_not_count():
    assert _sources("p|3", label_terms={"p|4": {"ta:relE"}}) == []


def test_a_gene_label_beyond_the_window_does_not_count_even_in_the_directon():
    genes = [_g(f"p|{i}", 100 * i + 1, 100 * i + 90, 1) for i in range(1, 7)]
    assert _sources("p|1", label_terms={"p|5": {"ta:relE"}}, genes=genes) == []
    assert _sources("p|1", label_terms={"p|4": {"ta:relE"}}, genes=genes) == [
        ("ta:relE", "p|4")]


def test_a_system_component_counts_on_either_strand():
    """p|4 is on the opposite strand and in no directon with p|3; the system rule still
    applies, because a defence or conjugation system spans both strands."""
    assert _sources("p|3", system_terms={"p|4": {"defence:Clover"}}) == [
        ("defence:Clover", "p|4")]


def test_an_orf_inside_a_system_carries_the_term_itself():
    assert _sources("p|3", system_terms={"p|3": {"conj:T4SS_typeF"}}) == [
        ("conj:T4SS_typeF", None)]


def test_an_orfs_own_gene_label_is_not_its_context():
    assert _sources("p|1", label_terms={"p|1": {"amr:blaTEM-1"}}) == []


def test_the_neighbour_window_wraps_on_a_circular_plasmid():
    genes = [_g("p|1", 11, 300, 1), _g("p|2", 1500, 1800, -1), _g("p|3", 4900, 4990, 1)]
    sources = orf_term_sources("p|1", _near(genes, "p|1", circular=True),
                               directon=_unit(genes, "p|1", circular=True, length=5000),
                               label_terms={"p|3": {"mge:transfer"}}, system_terms={})
    assert sources == [("mge:transfer", "p|3")]


# --- window covers the plasmid ---------------------------------------------------------

def test_the_window_covers_a_small_circular_plasmid_but_not_a_large_one():
    small = [_g(f"p|{i}", 100 * i + 1, 100 * i + 90, 1) for i in range(1, 8)]
    large = [_g(f"p|{i}", 100 * i + 1, 100 * i + 90, 1) for i in range(1, 9)]
    assert window_covers_plasmid(_near(small, "p|1", circular=True), len(small))
    assert not window_covers_plasmid(_near(large, "p|1", circular=True), len(large))


def test_on_a_linear_record_only_the_central_gene_sees_everything():
    genes = [_g(f"p|{i}", 100 * i + 1, 100 * i + 90, 1) for i in range(1, 8)]
    assert window_covers_plasmid(_near(genes, "p|4"), len(genes))
    assert not window_covers_plasmid(_near(genes, "p|1"), len(genes))


# --- per family: the lineage is the unit ------------------------------------------------

LINEAGE = {"A1": "A", "A2": "A", "B1": "B", "C1": "C"}


def _row(rows, term):
    (row,) = [r for r in rows if r["term"] == term]
    return row


def test_clonal_copies_in_one_lineage_count_once():
    """Two plasmids of lineage A carry the term, the one plasmid of lineage B does not:
    one lineage of two, not two occurrences of three."""
    occurrences = [("A1|1", [("amr:x", "A1|2")], False),
                   ("A2|1", [("amr:x", "A2|2")], False),
                   ("B1|1", [], False)]
    rows = family_term_rows("F", "dark", occurrences, family_of_orf={}, lineage_of=LINEAGE)
    row = _row(rows, "amr:x")
    assert (row["n_lineages"], row["n_lineages_with_term"], row["conservation"]) == (2, 1, 0.5)
    assert row["status"] == "SUCCESS" and row["term_type"] == "amr"
    assert row["family_set"] == "dark"
    assert list(row) == COLUMNS


def test_a_neighbour_in_the_focal_family_is_a_tandem_paralogue_and_is_excluded():
    """B1|2 is another copy of family F beside the focal ORF; its label says what F is,
    not what F's context is."""
    occurrences = [("A1|1", [("amr:x", "A1|2")], False),
                   ("B1|1", [("amr:x", "B1|2")], False)]
    rows = family_term_rows("F", "known", occurrences, family_of_orf={"B1|2": "F",
                                                                      "A1|2": "G"},
                            lineage_of=LINEAGE)
    assert _row(rows, "amr:x")["n_lineages_with_term"] == 1


def test_an_orf_inside_a_system_is_never_excluded_as_a_paralogue():
    occurrences = [("A1|1", [("defence:Clover", None)], False),
                   ("B1|1", [("defence:Clover", None)], False)]
    rows = family_term_rows("F", "dark", occurrences, family_of_orf={"A1|1": "F",
                                                                     "B1|1": "F"},
                            lineage_of=LINEAGE)
    assert _row(rows, "defence:Clover")["conservation"] == 1.0


def test_one_lineage_is_too_few():
    occurrences = [("A1|1", [("ta:relE", "A1|2")], True), ("A2|1", [], False)]
    rows = family_term_rows("F", "dark", occurrences, family_of_orf={}, lineage_of=LINEAGE)
    row = _row(rows, "ta:relE")
    assert row["status"] == TOO_FEW_LINEAGES
    assert row["n_lineages"] == 1 and row["conservation"] == ""
    assert row["window_covers_plasmid_fraction"] == 0.5


def test_a_family_with_no_terms_writes_no_rows():
    assert family_term_rows("F", "dark", [("A1|1", [], False)], family_of_orf={},
                            lineage_of=LINEAGE) == []


def test_a_plasmid_without_a_lineage_is_an_error():
    with pytest.raises(KeyError):
        family_term_rows("F", "dark", [("Z9|1", [("amr:x", "Z9|2")], False)],
                         family_of_orf={}, lineage_of=LINEAGE)


# --- the stage script ------------------------------------------------------------------

def _run_stage(d):
    """Two plasmids in two lineages, each carrying one dark ORF of family D and one ORF of
    known family K, each in a directon with a blaTEM-labelled ORF (family T). On pl1 the
    known ORF K1 is a CONJScan component; on pl2 the dark ORF sits beside a DefenseFinder
    component on the opposite strand."""
    ann = d / "plasmid_annotation.tsv"
    write_tsv(ann, ["plasmid_id", "orf_id", "start", "end", "strand", "annot_label",
                    "functional_class"],
              [["pl1", "pl1|1", 100, 400, "+", "", "NONE"],
               ["pl1", "pl1|2", 430, 700, "+", "bla", "FUNCTIONAL"],
               ["pl1", "pl1|3", 730, 1000, "+", "k", "FUNCTIONAL"],
               ["pl1", "pl1|4", 5000, 5300, "-", "", "NONE"],
               ["pl2", "pl2|1", 100, 400, "+", "", "NONE"],
               ["pl2", "pl2|2", 430, 700, "+", "bla", "FUNCTIONAL"],
               ["pl2", "pl2|3", 730, 1000, "+", "k", "FUNCTIONAL"],
               ["pl2", "pl2|4", 1100, 1300, "-", "def", "FUNCTIONAL"]])
    dark = d / "dark_families.tsv"
    write_tsv(dark, ["family_id", "representative", "n_members", "family_class",
                     "members"], [["intermediate:D", "D", 1, "ORPHAN", "D"]])
    fams = d / "protein_families.tsv"
    write_tsv(fams, ["family_id", "family_resolution", "representative", "members"],
              [["intermediate:D", "intermediate", "D", "D"],
               ["intermediate:K", "intermediate", "K", "K"],
               ["intermediate:T", "intermediate", "T", "T"],
               ["close:K", "close", "K", "K"]])
    pmap = d / "protein_map.tsv"
    pmap.write_text("D\tpl1|1,pl2|1\nT\tpl1|2,pl2|2\nK\tpl1|3,pl2|3\n"
                    "N\tpl1|4\nE\tpl2|4\n")
    labels = d / "protein_labels.tsv"
    write_tsv(labels, ["protein_id", "source", "tier", "kind", "label", "sub_label"],
              [["T", "card", "Strict", "card_amr_family", "TEM beta-lactamase", ""],
               ["T", "eggnog", "S4b", "kegg_ko", "ko:K18698", ""],
               ["K", "amrfinder", "", "amrfinder_gene", "merA", "STRESS/METAL"]])
    defence = d / "defence_systems.tsv"
    write_tsv(defence, ["orf_id", "plasmid_id", "system", "status"],
              [["pl2|4", "pl2", "defense-finder-models/DefenseFinder/Clover/Clover",
                "SUCCESS"]])
    conj = d / "conjugation_systems.tsv"
    write_tsv(conj, ["orf_id", "plasmid_id", "system", "system_id", "component", "status"],
              [["pl1|3", "pl1", "T4SS_typeF", "s1", "T4SS_F_traL", "SUCCESS"]])
    integrons = d / "integrons.tsv"
    write_tsv(integrons, ["plasmid_id", "start", "end"], [])
    is_el = d / "is_elements.tsv"
    write_tsv(is_el, ["plasmid_id", "start", "end"], [])
    master = d / "master.tsv"
    write_tsv(master, ["plasmid_id", "topology"], [["pl1", "linear"], ["pl2", "linear"]])
    lengths = d / "plasmid_lengths.tsv"
    write_tsv(lengths, ["plasmid_id", "length_bp"], [["pl1", 20000], ["pl2", 20000]])
    lineage = d / "plasmid_lineage.tsv"
    write_tsv(lineage, ["plasmid_id", "plasmid_lineage_cluster"], [["pl1", "L1"],
                                                                   ["pl2", "L2"]])
    out_fam, out_terms = d / "family_context.tsv", d / "family_context_terms.tsv"
    run_script("context_features.py", FakeSnakemake(
        input={"annotation": str(ann), "families": str(dark), "map": str(pmap),
               "defence": str(defence), "conjugation": str(conj),
               "integrons": str(integrons), "is_elements": str(is_el),
               "master": str(master), "lengths": str(lengths), "labels": str(labels),
               "all_families": str(fams), "lineage": str(lineage)},
        output={"families": str(out_fam), "terms": str(out_terms)},
        params={"context": {"max_operon_gap": 100, "neighbourhood_window": 3},
                "primary": "intermediate"}))
    return read_tsv(out_fam), read_tsv(out_terms)


def test_the_stage_keeps_the_rates_and_adds_cons_conj(fixture_dir):
    families, _ = _run_stage(fixture_dir)
    (row,) = families
    assert row["family_id"] == "intermediate:D" and row["n_units"] == "2"
    assert row["cons_conj"] == "0.0" and row["cons_defence"] == "0.0"
    assert row["cons_operon_with_annotated"] == "1.0"


def test_the_stage_writes_terms_for_dark_and_known_families(fixture_dir):
    _, terms = _run_stage(fixture_dir)
    assert list(terms[0]) == COLUMNS
    by = {(r["family_id"], r["term"]): r for r in terms}

    dark_amr = by[("intermediate:D", "amr:TEM beta-lactamase")]
    assert dark_amr["family_set"] == "dark"
    assert (dark_amr["n_lineages"], dark_amr["n_lineages_with_term"]) == ("2", "2")
    assert dark_amr["conservation"] == "1.0" and dark_amr["status"] == "SUCCESS"
    # K is two genes away in the same directon; its AMRFinderPlus STRESS/METAL label is metal:.
    assert by[("intermediate:D", "metal:merA")]["n_lineages_with_term"] == "2"
    # conj: from K1 on pl1 (same strand, a component within +-3); defence: from the
    # opposite-strand component on pl2.
    assert by[("intermediate:D", "conj:T4SS_typeF")]["conservation"] == "0.5"
    assert by[("intermediate:D", "defence:Clover")]["conservation"] == "0.5"

    known = by[("intermediate:K", "conj:T4SS_typeF")]
    assert known["family_set"] == "known" and known["conservation"] == "0.5"
    assert by[("intermediate:K", "amr:TEM beta-lactamase")]["family_set"] == "known"
    # The dark family is not also written as a known one.
    assert sum(r["family_id"] == "intermediate:D" and r["term"] == "amr:TEM beta-lactamase"
               for r in terms) == 1
    # Only the primary resolution, and no KEGG term.
    assert not any(fid.startswith("close:") for fid, _ in by)
    assert not any(t.startswith("kegg") or t.startswith("ko:") for _, t in by)
