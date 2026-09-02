import hashlib


def _seq_id(seq):
    return hashlib.sha256(seq.encode()).hexdigest()[:16]


def dereplicate(orfs):
    """Collapse identical sequences. Lossless: every orf_id appears in exactly one group.

    Returns (uniques, mapping) where uniques is [{"seq_id", "seq"}] and mapping is
    seq_id -> [orf_id, ...].
    """
    mapping = {}
    uniques = {}
    for o in orfs:
        sid = _seq_id(o["seq"])
        mapping.setdefault(sid, []).append(o["orf_id"])
        uniques.setdefault(sid, {"seq_id": sid, "seq": o["seq"]})
    return list(uniques.values()), mapping
