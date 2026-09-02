import pytest
from plasmidann.cascade import narrow


def test_resolved_and_unresolved_partition_the_input_exactly():
    """Each tier hands the next only what it could not name - and loses nothing."""
    all_ids = ["a", "b", "c", "d"]
    hit_ids = {"b", "d"}

    unresolved = narrow(all_ids, hit_ids)

    assert unresolved == ["a", "c"]
    assert len(unresolved) + len(hit_ids) == len(all_ids)


def test_a_hit_for_an_unknown_id_is_an_error():
    """A hit id absent from the query set means the search and the index disagree."""
    with pytest.raises(ValueError):
        narrow(["a", "b"], {"b", "zzz"})
