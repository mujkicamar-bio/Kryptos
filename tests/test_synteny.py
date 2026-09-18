"""Stage 9: synteny and context conservation (spec section 42).

Stage 8 asks what a dark ORF sits next to, once. This asks whether the arrangement RECURS.
The spec's example is the distinction:

    A - B - DARK - C - D
    A - B - DARK - C - D
    A - B - DARK - C
    A - B - DARK - C - D

Not four observations of "near B" but one conserved gene order seen four times, which is a
far stronger claim: order survives because the arrangement matters.
"""
from plasmidann import synteny


def spec_example():
    """The four occurrences from section 42, verbatim."""
    return [
        {"left": ["B", "A"], "right": ["C", "D"], "operon": True},
        {"left": ["B", "A"], "right": ["C", "D"], "operon": True},
        {"left": ["B", "A"], "right": ["C"], "operon": True},
        {"left": ["B", "A"], "right": ["C", "D"], "operon": True},
    ]


def test_the_spec_example_is_fully_conserved():
    result = synteny.conservation(spec_example())

    assert result["left_neighbor_conservation"] == 1.0
    assert result["right_neighbor_conservation"] == 1.0
    assert result["synteny_conservation"] == 1.0
    assert result["modal_synteny"] == "B|C"
    assert result["status"] == "SUCCESS"


def test_a_conserved_left_and_a_variable_right_are_reported_separately():
    """A real and common arrangement - the left gene may be the promoter-sharing partner.
    A single averaged context score would hide it, which is why section 42 names six
    measurements rather than one."""
    occurrences = [
        {"left": ["B"], "right": ["C"]},
        {"left": ["B"], "right": ["X"]},
        {"left": ["B"], "right": ["Y"]},
        {"left": ["B"], "right": ["Z"]},
    ]

    result = synteny.conservation(occurrences)

    assert result["left_neighbor_conservation"] == 1.0
    assert result["right_neighbor_conservation"] == 0.25
    assert result["synteny_conservation"] == 0.25, (
        "ordered synteny cannot exceed its weakest side")


def test_ordered_synteny_is_stricter_than_the_unordered_neighbourhood():
    """The same two genes in swapped positions are one neighbourhood but not one synteny."""
    occurrences = [
        {"left": ["B"], "right": ["C"]},
        {"left": ["C"], "right": ["B"]},
    ]

    result = synteny.conservation(occurrences)

    assert result["neighborhood_conservation"] == 1.0
    assert result["synteny_conservation"] == 0.5


def test_a_missing_neighbour_is_not_counted_as_a_failure_to_conserve():
    """An ORF at the end of a contig has no left neighbour. Counting that against the
    family would penalise it for where the assembler cut, not for its biology."""
    occurrences = [
        {"left": [], "right": ["C"]},
        {"left": ["B"], "right": ["C"]},
        {"left": ["B"], "right": ["C"]},
    ]

    result = synteny.conservation(occurrences)

    assert result["left_neighbor_conservation"] == 1.0, (
        "the truncated occurrence was counted as a conservation failure")
    assert result["right_neighbor_conservation"] == 1.0


def test_a_single_occurrence_is_not_measured():
    """Conservation over one occurrence is 1.0 by construction and means nothing. Emitting
    it would put every singleton at the top of a ranking of conserved context."""
    result = synteny.conservation([{"left": ["B"], "right": ["C"]}])

    assert result["status"] == "TOO_FEW_MEMBERS"
    assert result["synteny_conservation"] == ""


def test_a_family_with_no_context_at_all_is_reported_as_such():
    """Distinct from 'measured and found variable'."""
    result = synteny.conservation([{"left": [], "right": []},
                                   {"left": [], "right": []}])

    assert result["context_recurrence"] == 0
    assert result["status"] == "TOO_FEW_MEMBERS"


def test_the_occurrence_count_travels_with_the_measurement():
    """Conservation is computed over occurrences, so a family on forty redepositions of one
    plasmid shows perfect synteny from a single biological event. The count is reported so
    it can be read against the independent-lineage count from Stage 7."""
    result = synteny.conservation(spec_example())

    assert result["n_occurrences"] == 4
    assert result["context_recurrence"] == 4


def test_operon_conservation_is_its_own_measurement():
    occurrences = [
        {"left": ["B"], "right": ["C"], "operon": True},
        {"left": ["B"], "right": ["C"], "operon": False},
    ]

    result = synteny.conservation(occurrences)

    assert result["operon_like_conservation"] == 0.5
    assert result["synteny_conservation"] == 1.0, (
        "operon membership and gene order are different measurements")
