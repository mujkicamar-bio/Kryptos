"""Labelled cases for the uninformative-label detector, each asserted individually.

The functional names include the ones a comparable regex of an earlier project filed as
novel (RelE missing while RelB was present; TrfA, RepC, TrwC, Rop), which moved a
published percentage by more than 7 points.
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
    "putative protein",
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


def test_a_missing_or_blank_label_is_not_informative():
    """A label that failed to parse must leave the protein unknown, not name it."""
    from plasmidann.cascade import is_informative

    assert is_informative(None) is False
    assert is_informative("") is False
    assert is_informative("   ") is False
