from plasmidann.dereplicate import dereplicate


def test_identical_sequences_collapse_and_every_orf_is_accounted_for():
    """S2 is a compute optimisation only: no ORF may be lost or double-counted."""
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


def test_two_different_sequences_never_share_an_identity():
    """Protein identity is a truncated SHA-256. At 16 hex characters that is 64 bits, and
    across 3.5M unique proteins the birthday probability of a collision is about 3e-7 per
    run. If it fires, `setdefault` keeps the first sequence and the second, genuinely
    different protein is silently replaced by it everywhere downstream - and the existing
    losslessness assertion cannot catch it, because it checks that ORF counts are preserved,
    not that sequences are.

    A wider digest costs nothing. The check below is what makes the guarantee explicit."""
    orfs = [{"orf_id": "o1", "seq": "MKVL"}, {"orf_id": "o2", "seq": "MKVA"},
            {"orf_id": "o3", "seq": "MKVL"}]
    uniques, mapping = dereplicate(orfs)

    by_id = {u["seq_id"]: u["seq"] for u in uniques}
    assert len(by_id) == 2
    assert sorted(by_id.values()) == ["MKVA", "MKVL"]
    # 128 bits: at 3.5M sequences the collision probability is below 1e-26.
    assert all(len(sid) >= 32 for sid in by_id), (
        f"identity is only {min(len(s) for s in by_id)} hex characters wide")
