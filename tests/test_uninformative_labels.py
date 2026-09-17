"""Labelled test set for the uninformative-label detector.

The design requires this list to ship with measured recall, because the previous
project's equivalent regex had good precision and poor recall (RelE missing while
RelB was present; TrfA, RepC, TrwC, Rop filed as novel), which moved a published
percentage by more than 7 points.
"""
import pytest
from plasmidann.cascade import UNINFORMATIVE

UNINFORMATIVE_LABELS = [
    "hypothetical protein",
    "MULTISPECIES: hypothetical protein",
    "Uncharacterized protein",
    "putative uncharacterized protein",
    "conserved hypothetical protein",
    "DUF1234 domain-containing protein",
    "DUF domain-containing protein",
    "UPF0102 protein",
    "protein of unknown function",
    "predicted protein",
    "conserved protein",
    "unnamed protein product",
    "ORF",
]

FUNCTIONAL_LABELS = [
    "relaxase MobA",
    "Rep_1",
    "type I restriction enzyme HsdR",
    "putative relaxase",
    "TrfA",
    "RepC",
    "TrwC",
    "Rop",
    "RelE",
    "RelB",
    "conjugal transfer protein TraD",
    "ParE toxin",
]


@pytest.mark.parametrize("label", UNINFORMATIVE_LABELS)
def test_uninformative_labels_are_detected(label):
    assert UNINFORMATIVE.search(label), f"missed: {label}"


@pytest.mark.parametrize("label", FUNCTIONAL_LABELS)
def test_real_function_names_are_not_flagged(label):
    assert not UNINFORMATIVE.search(label), f"false positive: {label}"


def test_a_missing_label_is_not_informative():
    """v1 returned True here, because `not UNINFORMATIVE.search(None or "")` is True.

    A label that failed to parse would therefore be treated as naming a function, silently
    promoting the protein out of the dark set. Failing toward "we do not know" is the safe
    direction for a discovery pipeline.
    """
    from plasmidann.cascade import is_informative

    assert is_informative(None) is False
    assert is_informative("") is False
    assert is_informative("   ") is True   # whitespace is a label we cannot judge, not an absence
