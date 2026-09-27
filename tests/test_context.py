"""Genomic context (plasmidann.context) - the stage that decides whether the screen finds anything.

Both experimental successes in the FESNov study were selected by genomic context, not by
any novelty measure: one sat in the canonical che operon, the other beside antibiotic
resistance genes. Plasmids suit this better than metagenomes - small, gene-dense, modular,
with cargo organised into recognisable islands.
"""
import pytest

from plasmidann.context import directons, flanks, overlapping_islands


def _g(oid, start, end, strand):
    return {"orf_id": oid, "start": start, "end": end, "strand": strand}


# --- directons: putative transcriptional units --------------------------------------

def test_consecutive_genes_on_one_strand_with_short_gaps_form_one_directon():
    """A dark ORF INSIDE a transcriptional unit whose other members are annotated is a
    much stronger claim than one merely nearby: it is predicted to be co-transcribed."""
    genes = [_g("a", 1, 300, 1), _g("b", 320, 600, 1), _g("c", 640, 900, 1)]

    assert directons(genes, max_gap=100) == [["a", "b", "c"]]


def test_a_strand_change_ends_a_directon():
    genes = [_g("a", 1, 300, 1), _g("b", 320, 600, -1), _g("c", 640, 900, -1)]

    assert directons(genes, max_gap=100) == [["a"], ["b", "c"]]


def test_a_long_intergenic_gap_ends_a_directon():
    """Beyond ~100 nt the genes are unlikely to share a promoter."""
    genes = [_g("a", 1, 300, 1), _g("b", 900, 1200, 1)]

    assert directons(genes, max_gap=100) == [["a"], ["b"]]


def test_overlapping_genes_are_in_the_same_directon():
    """Overlapping start codons are common in tightly packed operons; the gap is negative."""
    genes = [_g("a", 1, 300, 1), _g("b", 295, 600, 1)]

    assert directons(genes, max_gap=100) == [["a", "b"]]


def test_a_gap_of_exactly_max_gap_stays_in_the_directon():
    """FESNov breaks a unit at intergenic regions ABOVE 100 nt, so 100 itself joins."""
    genes = [_g("a", 1, 300, 1), _g("b", 401, 700, 1)]

    assert directons(genes, max_gap=100) == [["a", "b"]]


# --- directons across the origin of a circular plasmid ---------------------------------

def test_a_directon_crossing_the_origin_of_a_circular_plasmid_is_one_unit():
    """The record's last run and first run are one operon on the molecule when the strand
    matches and the gap across the origin is small: here 50 nt (1000-950) + 10 nt."""
    genes = [_g("a", 11, 300, 1), _g("b", 320, 600, 1),
             _g("m", 700, 800, -1),
             _g("y", 820, 900, 1), _g("z", 910, 950, 1)]

    assert directons(genes, max_gap=100, circular=True, length=1000) == [
        ["y", "z", "a", "b"], ["m"]]


def test_a_linear_record_is_not_merged_across_its_ends():
    genes = [_g("a", 11, 300, 1), _g("m", 400, 500, -1), _g("z", 910, 950, 1)]

    assert directons(genes, max_gap=100) == [["a"], ["m"], ["z"]]


def test_a_strand_change_at_the_origin_keeps_the_runs_apart():
    genes = [_g("a", 11, 300, -1), _g("m", 800, 900, 1), _g("z", 910, 950, 1)]

    assert directons(genes, max_gap=100, circular=True, length=1000) == [
        ["a"], ["m", "z"]]


def test_a_long_gap_across_the_origin_keeps_the_runs_apart():
    """150 nt to the end of the molecule plus 10 nt past the origin is 160 > 100."""
    genes = [_g("a", 11, 300, 1), _g("m", 400, 500, -1), _g("z", 700, 850, 1)]

    assert directons(genes, max_gap=100, circular=True, length=1000) == [
        ["a"], ["m"], ["z"]]


def test_a_gene_spanning_the_origin_joins_the_run_after_it():
    """z runs 950..1000 then 1..40 (start > end); a starts 30 nt after it ends."""
    genes = [_g("a", 71, 300, 1), _g("m", 400, 500, -1), _g("z", 950, 40, 1)]

    assert directons(genes, max_gap=100, circular=True, length=1000) == [
        ["z", "a"], ["m"]]


def test_a_circle_that_is_one_directon_is_not_merged_with_itself():
    genes = [_g("a", 11, 300, 1), _g("b", 320, 600, 1), _g("c", 640, 990, 1)]

    assert directons(genes, max_gap=100, circular=True, length=1000) == [["a", "b", "c"]]


def test_a_circular_directon_needs_the_molecule_length():
    with pytest.raises(ValueError):
        directons([_g("a", 11, 300, 1), _g("m", 400, 500, -1)], circular=True)


# --- flanks ----------------------------------------------------------------------------

def test_the_flanks_are_the_genes_on_either_side_nearest_first():
    genes = [_g(x, i * 100 + 1, i * 100 + 90, 1) for i, x in enumerate("abcdefg")]

    assert flanks(genes, window=2)["d"] == (["c", "b"], ["e", "f"])


def test_the_flanks_at_the_end_of_a_linear_record_are_truncated_not_wrapped():
    genes = [_g(x, i * 100 + 1, i * 100 + 90, 1) for i, x in enumerate("abcde")]

    assert flanks(genes, window=2)["a"] == ([], ["b", "c"])


def test_genes_are_ordered_by_coordinate_not_by_input_order():
    genes = [_g("c", 201, 290, 1), _g("a", 1, 90, 1), _g("b", 101, 190, 1)]

    assert flanks(genes, window=1)["b"] == (["a"], ["c"])


# --- islands ---------------------------------------------------------------------------

def test_an_orf_inside_an_island_is_detected():
    islands = [{"name": "defence", "start": 1000, "end": 5000}]

    assert overlapping_islands(_g("x", 2000, 2500, 1), islands) == islands


def test_an_orf_partly_overlapping_an_island_counts_as_inside():
    islands = [{"name": "integron", "start": 1000, "end": 5000}]

    assert overlapping_islands(_g("x", 4800, 5400, 1), islands) == islands


# --- islands: a gene can sit in more than one, and can cross the origin ---------------

def test_a_gene_inside_two_islands_reports_both():
    """A gene in a defence system inside a cassette array is in both, whatever the order in
    which the islands are listed."""
    gene = {"orf_id": "o1", "start": 500, "end": 800, "strand": 1}
    islands = [{"name": "integron", "start": 100, "end": 1200},
               {"name": "defence", "start": 450, "end": 900}]
    assert {i["name"] for i in overlapping_islands(gene, islands)} == {"integron",
                                                                      "defence"}


def test_a_gene_in_no_island_reports_none():
    gene = {"orf_id": "o1", "start": 5000, "end": 5200, "strand": 1}
    assert overlapping_islands(gene, [{"name": "integron", "start": 1, "end": 400}]) == []


def test_a_gene_crossing_the_origin_still_matches_its_island():
    """An origin-spanning gene is written start > end (start..L then 1..end); it overlaps
    islands on either side of the origin."""
    wrapped = {"orf_id": "o1", "start": 4900, "end": 120, "strand": 1}
    assert [i["name"] for i in overlapping_islands(
        wrapped, [{"name": "defence", "start": 4800, "end": 4950}])] == ["defence"]
    assert [i["name"] for i in overlapping_islands(
        wrapped, [{"name": "integron", "start": 50, "end": 300}])] == ["integron"]
    assert overlapping_islands(
        wrapped, [{"name": "integron", "start": 1000, "end": 2000}]) == []


def test_a_circular_neighbourhood_wraps_across_the_origin():
    """On a circle the first gene's left neighbours are the last genes of the record; a
    gene is never its own neighbour, however small the circle."""
    genes = [{"orf_id": x, "start": 100 * i + 1, "end": 100 * i + 90, "strand": 1}
             for i, x in enumerate("abcde")]
    assert flanks(genes, window=2, circular=True)["a"] == (["e", "d"], ["b", "c"])
    assert flanks(genes, window=2)["a"] == ([], ["b", "c"])
    two = genes[:2]
    assert flanks(two, window=3, circular=True)["a"] == (["b"], ["b"])


# --- the neighbourhood diagram tool draws the pipeline's neighbourhood --------------------

def test_the_diagram_row_wraps_across_the_origin_of_a_circular_plasmid():
    """The focal gene a is first on a 1,000 bp circle; its left neighbours z and y are the
    record's last genes, drawn to its left, as context.flanks counts them."""
    import importlib.util
    import pathlib
    pytest.importorskip("matplotlib")
    tool = pathlib.Path(__file__).resolve().parents[1] / "tools" / "draw_neighbourhoods.py"
    spec = importlib.util.spec_from_file_location("draw_neighbourhoods", tool)
    draw = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(draw)
    genes = [{"orf_id": o, "start": str(s), "end": str(e), "strand": st}
             for o, s, e, st in (("a", 11, 100, "1"), ("b", 201, 300, "1"),
                                 ("y", 701, 800, "1"), ("z", 951, 5, "-1"))]

    circle = draw.row_layout(genes, "a", 2, True, 1000)
    assert [(r["orf_id"], x0, x1, fwd) for r, x0, x1, fwd in circle] == [
        ("y", -310, -211, True), ("z", -60, -6, False), ("a", 0, 89, True),
        ("b", 190, 289, True), ("y", 690, 789, True)]
    linear = draw.row_layout(genes, "a", 2, False, 1000)
    assert [r["orf_id"] for r, *_ in linear] == ["a", "b", "y"]
