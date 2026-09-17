"""Significance for the genomic-context association, over independent units.

Three defects in the statistic this replaces:

  1. No test. Enrichment was observed rate over background rate, so a family with three
     members all beside a relaxase scored identically to one with three hundred.
  2. The wrong unit. The rate was over family MEMBERS, which are homologs, often on
     near-identical plasmids. A protein family present on forty copies of one sequenced
     plasmid is one observation, not forty.
  3. No correction. Six hand-picked features needed none; the open label vocabulary means
     thousands of tests per family, where an uncorrected p-value is meaningless.
"""
import pytest

from plasmidann import enrich


def test_a_strong_association_is_significant():
    """Eighteen of twenty independent plasmids carrying this family have the category in
    the neighbourhood, against 5% of the corpus."""
    result = enrich.fisher_enrichment(k=18, n=20, K=500, N=10_000)

    assert float(result["p_value"]) < 1e-10
    assert result["enrichment"] > 3
    assert result["status"] == "SUCCESS"


def test_an_association_at_the_background_rate_is_not_significant():
    """The small-plasmid trap: on a six-gene plasmid a plus or minus three neighbourhood is
    the whole molecule, so everything co-occurs with everything. Observed equal to
    background must come out unremarkable, whatever the absolute rate."""
    result = enrich.fisher_enrichment(k=5, n=10, K=5_000, N=10_000)

    assert float(result["p_value"]) > 0.5
    assert result["enrichment"] == pytest.approx(1.0, abs=0.01)


def test_a_single_unit_cannot_be_tested():
    """One plasmid is an anecdote. Returning a p-value for n=1 would let a family found
    once rank alongside one found on two hundred independent plasmids."""
    result = enrich.fisher_enrichment(k=1, n=1, K=100, N=10_000)

    assert result["status"] == "TOO_FEW_MEMBERS"
    assert result["p_value"] == ""


def test_a_family_carrying_more_of_a_category_than_the_corpus_is_an_error():
    """The family's plasmids are a subset of the corpus, so the corpus count can never be
    below the family count. Violating that means the caller built its two counts over
    different sets, which would produce a plausible-looking p-value from nonsense - and it
    is the one inconsistency that cannot be seen by reading the output."""
    with pytest.raises(ValueError, match="not a subset"):
        enrich.fisher_enrichment(k=3, n=10, K=0, N=10_000)


def test_no_observation_and_no_background_is_no_signal():
    """Zero over zero is 1.0, not an error and not infinity. Because K >= k always, a zero
    background forces a zero observation, so the enrichment is never infinite."""
    result = enrich.fisher_enrichment(k=0, n=10, K=0, N=10_000)

    assert result["enrichment"] == 1.0


def test_an_impossible_count_is_an_error_not_a_silent_result():
    """k > n means the caller's units disagree with its counts, which would produce a
    plausible-looking p-value from nonsense."""
    with pytest.raises(ValueError):
        enrich.fisher_enrichment(k=11, n=10, K=100, N=10_000)
    with pytest.raises(ValueError):
        enrich.fisher_enrichment(k=1, n=10, K=100, N=10)


def test_the_family_is_tested_against_the_rest_of_the_corpus():
    """The two-by-two table compares the family's plasmids against the plasmids that are
    NOT the family's. Including them in both margins would test the family against a
    background it contributes to, which shrinks any real effect."""
    result = enrich.fisher_enrichment(k=10, n=10, K=10, N=1_000)

    assert float(result["p_value"]) < 1e-10, (
        "a category unique to this family came out unremarkable, so the family is being "
        "tested against a background that includes it")


def test_benjamini_hochberg_preserves_input_order():
    """The caller joins q-values back to rows by position, so a sorted return would
    silently attach every q-value to the wrong category."""
    q = enrich.benjamini_hochberg([0.5, 0.001, 0.04, 0.2])

    assert len(q) == 4
    assert q[1] == min(q), "the smallest p-value did not get the smallest q-value"


def test_benjamini_hochberg_matches_the_published_procedure():
    """Benjamini and Hochberg 1995, J R Stat Soc B 57:289. With m=4, the q-value of the
    i-th smallest p is min over j >= i of (m/j) * p_j, which is monotone by construction."""
    q = enrich.benjamini_hochberg([0.01, 0.02, 0.03, 0.04])

    assert q == pytest.approx([0.04, 0.04, 0.04, 0.04])


def test_benjamini_hochberg_ignores_untested_entries():
    """A category that could not be tested has no p-value, and treating an empty string as
    zero would give it the strongest q-value in the table."""
    q = enrich.benjamini_hochberg([0.01, "", 0.02])

    assert q[1] == ""
    assert q[0] != ""


def test_no_q_value_exceeds_one():
    q = enrich.benjamini_hochberg([0.9, 0.95, 0.99])

    assert all(v <= 1.0 for v in q)
