"""S2s: which proteins the cascade annotates (plasmidann.selection).

A broad family is annotated when a small plasmid carries a protein in it that Tier 0 has
not explained; then every member Tier 0 has not explained is annotated, whatever plasmid
it is on. Nothing else is.
"""
from plasmidann.selection import select


def test_a_family_with_an_unexplained_small_protein_is_annotated_on_both_sides():
    clusters = {"s": ["s", "large_dark", "large_ps"]}
    assert select(clusters, on_small={"s"}, ps_annotated={"large_ps"}) == {
        "s", "large_dark"}


def test_a_family_whose_small_proteins_tier0_explains_is_not_annotated():
    """The small plasmid's protein is already known, so the family says nothing new about
    small plasmids; its unexplained large-plasmid members are not searched."""
    clusters = {"s": ["s", "large_dark"]}
    assert select(clusters, on_small={"s"}, ps_annotated={"s"}) == set()


def test_a_family_of_large_plasmid_proteins_alone_is_not_annotated():
    clusters = {"a": ["a", "b"], "s": ["s"]}
    assert select(clusters, on_small={"s"}, ps_annotated=set()) == {"s"}
