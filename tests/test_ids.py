from darkorf.ids import family_id


def test_family_id_is_derived_from_its_representative():
    assert family_id("close", "abc123") == "close:abc123"
