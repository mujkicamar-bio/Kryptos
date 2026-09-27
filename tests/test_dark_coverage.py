"""dark_covered_fraction: the coverage of a dark protein by the hits that named nothing.

annot_completeness is NONE for every dark protein by construction, because the explained
fraction is built from informative spans only (measured: 71/71 and 703/703). The
uninformative spans give the dark set its own coverage measure.
"""
from plasmidann.cascade import explained_fraction, n_dark_databases, uninformative_spans


def test_a_full_length_unnamed_protein_is_separated_from_a_fragment():
    """95% covered by 'hypothetical protein' is a real, conserved, full-length protein
    nobody has named. A 20-aa fragment hit is much weaker evidence."""
    full = [{"label": "hypothetical protein", "start": 5, "end": 195}]
    fragment = [{"label": "hypothetical protein", "start": 5, "end": 25}]

    assert explained_fraction(200, uninformative_spans(full)) == 0.955
    assert explained_fraction(200, uninformative_spans(fragment)) == 0.105


def test_named_hits_do_not_count_toward_dark_coverage():
    """The two fractions are complements over disjoint evidence, not the same number."""
    hits = [
        {"label": "relaxase MobA", "start": 1, "end": 100},
        {"label": "hypothetical protein", "start": 101, "end": 200},
    ]

    assert uninformative_spans(hits) == [(101, 200)]


def test_overlapping_unnamed_hits_are_merged_not_summed():
    """Three databases calling the same region hypothetical cover that region once."""
    hits = [
        {"label": "hypothetical protein", "start": 1, "end": 120},
        {"label": "MULTISPECIES: hypothetical protein", "start": 100, "end": 200},
    ]

    assert explained_fraction(200, uninformative_spans(hits)) == 1.0


def test_independent_databases_calling_a_protein_unnamed_are_counted():
    """Three databases agreeing a protein is real and unnamed is stronger evidence than
    one. Counted per tier, since two hits from one tier are not independent."""
    hits = [
        {"label": "hypothetical protein", "tier": "T3"},
        {"label": "hypothetical protein", "tier": "T4"},
        {"label": "DUF1234 domain-containing protein", "tier": "T4"},
        {"label": "relaxase MobA", "tier": "T1"},
    ]

    assert n_dark_databases(hits) == 2


def test_a_protein_nothing_has_seen_has_no_dark_databases():
    assert n_dark_databases([]) == 0
    assert n_dark_databases([{"label": "relaxase MobA", "tier": "T1"}]) == 0
