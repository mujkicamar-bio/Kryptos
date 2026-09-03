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
