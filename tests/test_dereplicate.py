from plasmidann.dereplicate import dereplicate


def test_identical_sequences_collapse_and_every_orf_is_accounted_for():
    """Dereplication is a compute optimisation only: no ORF may be lost or double-counted."""
    orfs = [
        {"orf_id": "p1|1", "seq": "MKV"},
        {"orf_id": "p1|2", "seq": "MAAA"},
        {"orf_id": "p2|1", "seq": "MKV"},
        {"orf_id": "p3|7", "seq": "MKV"},
    ]

    uniques, mapping = dereplicate(orfs)

    assert len(uniques) == 2
    assert sum(len(v) for v in mapping.values()) == len(orfs)
    assert sorted(o for v in mapping.values() for o in v) == ["p1|1", "p1|2", "p2|1", "p3|7"]
    assert sorted(u["seq"] for u in uniques) == ["MAAA", "MKV"]
    assert all(len(u["seq_id"]) == 32 for u in uniques)

