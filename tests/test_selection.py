"""S2s: which proteins the cascade annotates (plasmidann.selection).

A broad family is annotated when a small plasmid carries a protein in it that Tier 0 has
not explained; then every member Tier 0 has not explained is annotated, whatever plasmid
it is on. Nothing else is. An AntiFam-flagged protein is never searched and, like a Tier 0
protein, does not open its family for searching.
"""
from plasmidann.selection import select


def test_a_family_with_an_unexplained_small_protein_is_annotated_on_both_sides():
    clusters = {"s": ["s", "large_dark", "large_ps"]}
    assert select(clusters, on_small={"s"}, skipped={"large_ps"}) == {
        "s", "large_dark"}


def test_a_family_whose_small_proteins_tier0_explains_is_not_annotated():
    """The small plasmid's protein is already known, so the family says nothing new about
    small plasmids; its unexplained large-plasmid members are not searched."""
    clusters = {"s": ["s", "large_dark"]}
    assert select(clusters, on_small={"s"}, skipped={"s"}) == set()


def test_a_family_of_large_plasmid_proteins_alone_is_not_annotated():
    clusters = {"a": ["a", "b"], "s": ["s"]}
    assert select(clusters, on_small={"s"}, skipped=set()) == {"s"}


def test_an_antifam_flagged_small_protein_does_not_open_its_family():
    """The skipped set holds Tier 0 proteins AND AntiFam-flagged ones: an artefact on a
    small plasmid is not a protein still to be explained, so its family's large-plasmid
    members are not searched on its account, and the artefact itself never is."""
    clusters = {"art": ["art", "large_dark"], "s": ["s", "art2"]}
    assert select(clusters, on_small={"art", "s", "art2"},
                  skipped={"art", "art2"}) == {"s"}
