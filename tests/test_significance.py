"""The E-value gates whether a hit exists at all; coverage then decides the class.

The dangerous quadrant is high coverage with a weak E-value: a spurious long alignment
silently removes a genuine dark protein from the screening pool, and no downstream stage
can recover it. That is the worst error this project can make, so significance is applied
before a hit is allowed to contribute a span, a label or a row.
"""
import pytest
from plasmidann.cascade import passes_significance


def test_a_hit_worse_than_the_tier_threshold_is_rejected():
    """T2 reports domains at i-Evalues up to 2000; 7.2% of its resolutions depended on them."""
    assert passes_significance("2000", 1e-5) is False


def test_a_hit_at_the_tier_threshold_is_kept():
    """The threshold is inclusive: a hit exactly at the declared cutoff counts."""
    assert passes_significance("1e-5", 1e-5) is True


def test_a_strong_hit_is_kept():
    assert passes_significance("6.5e-40", 1e-5) is True


def test_a_tier_with_curated_thresholds_applies_no_evalue_floor():
    """T1 uses Pfam gathering thresholds, set per family by its curator. For some short
    families GA is looser than any global cut, so a blanket floor would override curation -
    degrading the highest-quality signal in the cascade."""
    assert passes_significance("0.5", None) is True


def test_an_unparsable_evalue_is_rejected_when_a_floor_is_set():
    """A missing E-value where one was required means the parser and the tool disagree."""
    assert passes_significance(None, 1e-5) is False
    assert passes_significance("", 1e-5) is False
    assert passes_significance("not-a-number", 1e-5) is False


def test_an_unparsable_evalue_is_kept_when_no_floor_is_set():
    """Bit-score tiers legitimately have no E-value to check."""
    assert passes_significance(None, None) is True


@pytest.mark.parametrize("ev", ["1e-300", "0", "0.0"])
def test_numerically_extreme_values_do_not_raise(ev):
    assert passes_significance(ev, 1e-5) is True
