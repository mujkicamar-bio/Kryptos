"""Every table has a declared column contract (spec §4.1). A stage that writes a column
nobody declared is a stage whose output nothing downstream can rely on."""
import pytest
import pyarrow as pa
from darkorf import schemas


def test_every_spec_table_is_declared():
    for table in ["plasmids", "orf_occurrences", "proteins", "controls", "run_manifest"]:
        assert table in schemas.TABLES
        assert len(schemas.columns(table)) > 0


def test_orf_occurrences_carries_the_cds_nucleotide_sequence():
    """Spec §8.3. Without it RNAcode and dN/dS become one-shot: re-running them later means
    re-extracting CDS from every plasmid in the collection."""
    assert "cds_nucleotide_sequence" in schemas.columns("orf_occurrences")


def test_validate_frame_rejects_a_missing_column():
    frame = pa.table({"plasmid_id": ["p1"]})
    with pytest.raises(ValueError, match="missing"):
        schemas.validate_frame("orf_occurrences", frame)


def test_validate_frame_rejects_an_undeclared_column():
    columns = {name: ["x"] for name in schemas.columns("plasmids")}
    columns["surprise"] = ["x"]
    with pytest.raises(ValueError, match="undeclared"):
        schemas.validate_frame("plasmids", pa.table(columns))


def test_the_protein_labels_table_is_declared():
    """The substrate for the functional grouping is a declared table, not a loose file.
    An undeclared column here is a label kind nothing downstream can group."""
    cols = schemas.columns("protein_labels")

    for required in ("protein_id", "source", "kind", "label", "database",
                     "database_version"):
        assert required in cols, f"protein_labels has no {required} column"
