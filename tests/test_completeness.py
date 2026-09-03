import pytest
from plasmidann.cascade import completeness, narrow_by_explained


def test_a_nearly_complete_explanation_is_full():
    assert completeness(0.95) == "FULL"


def test_a_half_explained_protein_is_partial_not_full():
    """0.51 and 0.99 must not collapse into one label."""
    assert completeness(0.60) == "PARTIAL"
    assert completeness(0.51) == "PARTIAL"


def test_a_lone_small_domain_is_a_fragment():
    assert completeness(0.15) == "FRAGMENT"


def test_nothing_explained_is_none():
    assert completeness(0.0) == "NONE"


def test_narrowing_keeps_searching_until_half_the_protein_is_explained():
    """A 15%-covered protein must keep descending the cascade, not stop at tier 1."""
    explained = {"a": 0.95, "b": 0.15, "c": 0.0}

    still_open = narrow_by_explained(["a", "b", "c"], explained, 0.5)

    assert still_open == ["b", "c"]


def test_narrowing_refuses_an_id_it_never_queried():
    with pytest.raises(ValueError):
        narrow_by_explained(["a"], {"zzz": 0.9}, 0.5)
