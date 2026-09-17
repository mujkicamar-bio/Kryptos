"""S8: genomic context - the stage that decides whether the screen finds anything.

Both experimental successes in the FESNov study were selected by genomic context, not by
any novelty measure: one sat in the canonical che operon, the other beside antibiotic
resistance genes. Plasmids suit this better than metagenomes - small, gene-dense, modular,
with cargo organised into recognisable islands.
"""
import pytest
from plasmidann.context import (overlapping_islands, directons, neighbourhood, overlapping_island,
                                context_conservation, background_rate, enrichment)


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


# --- neighbourhood -------------------------------------------------------------------

def test_the_neighbourhood_is_the_genes_on_either_side():
    genes = [_g(x, i * 100 + 1, i * 100 + 90, 1) for i, x in enumerate("abcdefg")]

    assert neighbourhood(genes, "d", window=2) == ["b", "c", "e", "f"]


def test_a_neighbourhood_at_the_end_of_a_plasmid_is_truncated_not_wrapped():
    """Wrapping would be right for a circular molecule but must be an explicit choice,
    not an accident of indexing."""
    genes = [_g(x, i * 100 + 1, i * 100 + 90, 1) for i, x in enumerate("abcde")]

    assert neighbourhood(genes, "a", window=2) == ["b", "c"]


# --- islands -------------------------------------------------------------------------

def test_an_orf_inside_a_defence_island_is_detected():
    """A dark ORF inside a defence island is a defence-system candidate."""
    islands = [{"name": "defence", "start": 1000, "end": 5000}]

    assert overlapping_island(_g("x", 2000, 2500, 1), islands)["name"] == "defence"


def test_an_orf_outside_every_island_returns_nothing():
    islands = [{"name": "defence", "start": 1000, "end": 5000}]

    assert overlapping_island(_g("x", 6000, 6500, 1), islands) is None


def test_an_orf_partly_overlapping_an_island_counts_as_inside():
    islands = [{"name": "integron", "start": 1000, "end": 5000}]

    assert overlapping_island(_g("x", 4800, 5400, 1), islands)["name"] == "integron"


# --- family-level aggregation --------------------------------------------------------

def test_context_conservation_is_the_fraction_of_members_sharing_an_association():
    """The family-level statistic is CONSERVATION of context, not a single instance.
    FESNov thresholded this at 90% for its high-confidence set."""
    members = [
        {"orf_id": "a", "context": {"defence"}},
        {"orf_id": "b", "context": {"defence"}},
        {"orf_id": "c", "context": {"amr"}},
        {"orf_id": "d", "context": {"defence"}},
    ]

    assert context_conservation(members, "defence") == 0.75
    assert context_conservation(members, "amr") == 0.25


def test_an_empty_family_has_no_conservation():
    assert context_conservation([], "defence") == 0.0


# --- the small-plasmid trap ----------------------------------------------------------

def test_co_occurrence_is_measured_against_a_background():
    """On a 5 kb cryptic plasmid with six genes, +/-3 neighbours IS the whole plasmid, so
    everything co-occurs with everything. Raw frequency would rank the smallest plasmids
    as the most informative when they are the least."""
    # 'defence' appears near 60% of family members, but near 60% of ALL proteins too.
    assert enrichment(observed=0.60, background=0.60) == pytest.approx(1.0)
    # The same 60% is meaningful when the background is 5%.
    assert enrichment(observed=0.60, background=0.05) == pytest.approx(12.0)


def test_background_rate_is_computed_over_the_whole_corpus():
    all_orfs = [{"context": {"defence"}}, {"context": set()},
                {"context": {"amr"}}, {"context": {"defence"}}]

    assert background_rate(all_orfs, "defence") == 0.5


def test_enrichment_against_a_zero_background_does_not_divide_by_zero():
    assert enrichment(observed=0.4, background=0.0) == float("inf")
    assert enrichment(observed=0.0, background=0.0) == 1.0


# --- islands: a gene can sit in more than one, and can cross the origin ---------------

def test_a_gene_inside_two_islands_reports_both():
    """`overlapping_island` returned the FIRST match and stopped. Defence intervals are
    appended after integron intervals, so a dark ORF inside a defence system that also sits
    in a cassette array was only ever labelled 'integron' - and the defence_island stratum,
    175 of the 1,000 constructs, was unreachable behind it."""
    gene = {"orf_id": "o1", "start": 500, "end": 800, "strand": 1}
    islands = [{"name": "integron", "start": 100, "end": 1200},
               {"name": "defence", "start": 450, "end": 900}]
    assert {i["name"] for i in overlapping_islands(gene, islands)} == {"integron",
                                                                      "defence"}


def test_a_gene_in_no_island_reports_none():
    gene = {"orf_id": "o1", "start": 5000, "end": 5200, "strand": 1}
    assert overlapping_islands(gene, [{"name": "integron", "start": 1, "end": 400}]) == []


def test_a_gene_crossing_the_origin_still_matches_its_island():
    """S1 reconstructs genes broken by linearising a circular plasmid, and writes them in
    the GenBank join() convention: start > end, read as start..L then 1..end. 160,375 ORFs
    are in that state. A plain interval test asks `start <= island.end and end >=
    island.start`, which for such a gene is simply false, so every reconstructed
    origin-spanning gene was invisible to every island - on a molecule where 94% of records
    are circular."""
    wrapped = {"orf_id": "o1", "start": 4900, "end": 120, "strand": 1}
    assert [i["name"] for i in overlapping_islands(
        wrapped, [{"name": "defence", "start": 4800, "end": 4950}])] == ["defence"]
    assert [i["name"] for i in overlapping_islands(
        wrapped, [{"name": "integron", "start": 50, "end": 300}])] == ["integron"]
    assert overlapping_islands(
        wrapped, [{"name": "integron", "start": 1000, "end": 2000}]) == []
