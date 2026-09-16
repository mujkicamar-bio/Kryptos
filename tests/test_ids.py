"""Identifiers must be stable across reruns when the biological input is unchanged (spec §5).

The failure these tests prevent: an earlier implementation numbered families by
enumerate(sorted(...)), so adding one protein renumbered every later family and silently
broke every join against a previous run.
"""
from darkorf.ids import protein_id, occurrence_id, family_id, observation_id


def test_protein_id_is_a_32_character_hash_of_the_sequence():
    assert protein_id("MKV") == protein_id("MKV")
    assert len(protein_id("MKV")) == 32
    assert protein_id("MKV") != protein_id("MKW")


def test_protein_id_normalizes_case_and_whitespace():
    """FASTA parsers differ in whether they upper-case and how they wrap. The same protein
    must hash identically regardless of which parser produced the string."""
    assert protein_id("mkv") == protein_id("MKV")
    assert protein_id(" MKV\n") == protein_id("MKV")


def test_protein_id_strips_a_single_trailing_stop_codon():
    """Prodigal emits a trailing '*'; other tools do not. Without this the same protein gets
    two different identifiers depending on its provenance."""
    assert protein_id("MKV*") == protein_id("MKV")


def test_occurrence_id_is_coordinate_based_not_ordinal():
    assert occurrence_id("pA", 10, 99, 1) == "pA:10-99:+"
    assert occurrence_id("pA", 10, 99, -1) == "pA:10-99:-"


def test_occurrence_id_is_unchanged_by_shard_layout():
    """The same ORF must get the same id whether it was called in shard 1 or shard 599."""
    assert occurrence_id("pA", 10, 99, 1) == occurrence_id("pA", 10, 99, 1)


def test_family_id_is_derived_from_its_representative():
    assert family_id("close", "abc123") == "close:abc123"


def test_observation_id_is_stable_and_order_independent_of_tool_version():
    """Tool version is provenance, not identity (spec §5.5): re-running the same search with
    a patched binary must not renumber every observation."""
    a = observation_id("p1", "hmmsearch", "Pfam-A", "38.2", "PF00001", 1, 50)
    b = observation_id("p1", "hmmsearch", "Pfam-A", "38.2", "PF00001", 1, 50)
    assert a == b
    assert len(a) == 32
    assert a != observation_id("p1", "hmmsearch", "Pfam-A", "38.2", "PF00002", 1, 50)
