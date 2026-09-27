import pytest

from plasmidann.cascade import dark_evidence

CASES = [
    # strongest: a curator built and named a family for it
    (["DUF4054 domain-containing protein"], "CURATED_FAMILY"),
    (["UPF0102 protein"], "CURATED_FAMILY"),
    (["DUF2201_N"], "CURATED_FAMILY"),
    # homologs exist, no function
    (["conserved hypothetical protein"], "CONSERVED"),
    (["conserved protein"], "CONSERVED"),
    # one algorithm's output and nothing more
    (["hypothetical protein"], "PREDICTED_ONLY"),
    (["MULTISPECIES: hypothetical protein"], "PREDICTED_ONLY"),
    (["predicted protein"], "PREDICTED_ONLY"),
    (["unnamed protein product"], "PREDICTED_ONLY"),
    # nothing named it at any tier
    ([], "NONE"),
]


@pytest.mark.parametrize("labels,expected", CASES)
def test_each_label_maps_to_its_rung(labels, expected):
    assert dark_evidence(labels) == expected


def test_the_strongest_rung_across_all_tiers_wins():
    """A protein called 'hypothetical' at one tier and DUF at another is a DUF protein."""
    assert dark_evidence(["hypothetical protein", "DUF1234 domain-containing protein"]) \
        == "CURATED_FAMILY"


def test_conserved_outranks_plain_hypothetical_in_one_label():
    assert dark_evidence(["conserved hypothetical protein"]) == "CONSERVED"
