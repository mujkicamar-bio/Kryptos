"""Protein identity and dereplication: identical sequences share one seq_id."""
import hashlib

# 128-bit truncated SHA-256: collision probability about 2e-26 for 3.5 M sequences.
ID_HEX_CHARS = 32


def sequence_id(seq):
    """The pipeline's identity of a protein sequence: its truncated SHA-256."""
    return hashlib.sha256(seq.encode()).hexdigest()[:ID_HEX_CHARS]


def dereplicate(orfs):
    """Collapse identical sequences. Lossless: every orf_id appears in exactly one group.

    Returns (uniques, mapping) where uniques is [{"seq_id", "seq"}] and mapping is
    seq_id -> [orf_id, ...].
    """
    mapping = {}
    uniques = {}
    for o in orfs:
        sid = sequence_id(o["seq"])
        mapping.setdefault(sid, []).append(o["orf_id"])
        uniques.setdefault(sid, {"seq_id": sid, "seq": o["seq"]})
    return list(uniques.values()), mapping
