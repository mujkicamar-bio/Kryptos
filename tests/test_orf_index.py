import pytest

from plasmidann.orfindex import assign_orf_ids


def test_ids_follow_coordinate_order_not_input_order():
    """ORF ids number every ORF on a plasmid by coordinate, whatever order they arrive in."""
    orfs = [
        {"plasmid_id": "p1", "start": 900, "end": 1200},
        {"plasmid_id": "p1", "start": 100, "end": 400},
        {"plasmid_id": "p1", "start": 500, "end": 800},
    ]

    result = assign_orf_ids(orfs)

    assert [o["orf_id"] for o in result] == ["p1|1", "p1|2", "p1|3"]
    assert [o["start"] for o in result] == [100, 500, 900]


def test_reindexing_an_already_indexed_set_is_refused():
    """A subset must not be renumbered from 1."""
    orfs = [
        {"plasmid_id": "p1", "start": 100, "end": 400},
        {"plasmid_id": "p1", "start": 500, "end": 800},
        {"plasmid_id": "p1", "start": 900, "end": 1200},
    ]
    indexed = assign_orf_ids(orfs)
    subset = [o for o in indexed if o["start"] >= 500]

    with pytest.raises(ValueError):
        assign_orf_ids(subset)
