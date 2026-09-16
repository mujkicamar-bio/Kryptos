"""Deterministic identifier construction (spec §5).

Every identifier in this pipeline is a pure function of biological content. Nothing is
numbered by position, enumeration order, or shard layout, because all three change between
runs for reasons that have nothing to do with the biology.
"""
import hashlib

# 32 hex characters is 128 bits. At 3.5 million proteins the collision probability is
# around 1e-26, so the truncation is safe and the identifiers stay readable in a terminal.
_ID_LENGTH = 32


def _normalize_protein_sequence(sequence):
    """Canonical form of a protein sequence for hashing.

    Upper-cased, whitespace removed, and a single trailing stop codon dropped - Prodigal
    emits '*' and most other tools do not, and the same protein must not hash differently
    because of which tool produced the string.
    """
    cleaned = "".join(sequence.split()).upper()
    return cleaned[:-1] if cleaned.endswith("*") else cleaned


def protein_id(sequence):
    """Stable identifier for a unique protein sequence (spec §5.3)."""
    canonical = _normalize_protein_sequence(sequence)
    return hashlib.sha256(canonical.encode()).hexdigest()[:_ID_LENGTH]


def occurrence_id(plasmid_id, start, end, strand):
    """Stable identifier for one ORF occurrence (spec §5.2).

    Coordinate-based rather than ordinal, so the identifier does not depend on how many ORFs
    were called before it or on which shard called it. Origin-spanning ORFs keep their
    unrotated coordinates, where end < start is legitimate (spec §8.4).
    """
    sign = "+" if strand > 0 else "-"
    return f"{plasmid_id}:{start}-{end}:{sign}"


def family_id(resolution, representative_protein_id):
    """Stable identifier for a sequence family (spec §5.4).

    Derived from the cluster representative, so adding an unrelated protein elsewhere in the
    dataset cannot renumber this family.
    """
    return f"{resolution}:{representative_protein_id}"


def observation_id(protein_id_value, tool, database, database_version, target, qstart, qend):
    """Stable identifier for one annotation/context/structure observation (spec §5.5).

    Tool version is deliberately excluded: it is recorded as provenance, but including it
    would give every observation a new identity after a patch release.
    """
    key = "\t".join(
        [protein_id_value, tool, database, database_version, target, str(qstart), str(qend)]
    )
    return hashlib.sha256(key.encode()).hexdigest()[:_ID_LENGTH]
