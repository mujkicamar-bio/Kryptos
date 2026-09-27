"""Producers that write a status must use the vocabulary's spelling."""
from darkorf import status


def test_the_lineage_count_statuses_use_the_central_spelling():
    """Synteny (Stage 9) and the context terms (S8c) both report a measurement over too few
    independent lineages. They must write the one string the vocabulary declares, or a
    reader filtering on the status would see two different statements."""
    from plasmidann import context_terms, synteny

    assert context_terms.TOO_FEW_LINEAGES == status.TOO_FEW_LINEAGES
    one_lineage = [{"left": ["a"], "right": ["b"], "operon": False, "lineage": "L1"}] * 3
    assert synteny.conservation(one_lineage)["status"] == status.TOO_FEW_LINEAGES


def test_synteny_without_named_neighbours_writes_the_central_no_context():
    """Synteny (Stage 9) reports NO_CONTEXT when no occurrence has a named neighbour. The
    string must be one the vocabulary declares."""
    from plasmidann import synteny

    no_neighbours = [{"left": [], "right": [], "operon": False, "lineage": "L1"}]
    assert synteny.conservation(no_neighbours)["status"] == status.NO_CONTEXT
