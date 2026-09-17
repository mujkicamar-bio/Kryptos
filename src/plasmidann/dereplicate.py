import hashlib

# Width of the protein identity, in hex characters.
#
# 16 characters is 64 bits, and across 3.5M unique proteins the birthday probability of a
# collision is about 3e-7 per run. Small, but the failure is silent and unrecoverable: the
# first sequence wins and the second, genuinely different protein is replaced by it
# everywhere downstream. The losslessness assertion in workflow/scripts/dereplicate.py
# cannot catch it either - it checks that ORF counts are preserved, not that sequences are.
#
# 32 characters is 128 bits, which puts the same probability below 1e-26, and costs 16 bytes
# per protein.
ID_HEX_CHARS = 32


def _seq_id(seq):
    return hashlib.sha256(seq.encode()).hexdigest()[:ID_HEX_CHARS]


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
