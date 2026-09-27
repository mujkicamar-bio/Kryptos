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
    # The whole name inside a DIAMOND title: accession, [organism], Swiss-Prot 'Full='.
    "WP_1.1 putative protein [Escherichia coli]",
    "WP_1.1 MULTISPECIES: conserved protein [Enterobacteriaceae]",
    "P1.1 RecName: Full=Putative protein; AltName: Full=X [Escherichia coli]",
    "P2.1 predicted protein",
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
    # putative / predicted / conserved protein as part of a longer name
    "WP_1.1 putative protein kinase [Escherichia coli]",
    "Putative protein-export protein",
    "P1 predicted protein kinase",
    "WP_1.1 conserved protein of the ABC transporter family",
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
