"""The relational store (spec §4.2): Parquet on disk, DuckDB for queries, and never the
whole occurrence table in pandas."""
import pyarrow as pa
import pytest
from darkorf import store, schemas


def _plasmids_frame():
    return pa.table({name: ["x"] for name in schemas.columns("plasmids")})


def test_write_then_read_round_trips(tmp_path):
    path = store.table_path(tmp_path, "plasmids")
    store.write_table(path, "plasmids", _plasmids_frame())
    assert store.read_table(path).column_names == schemas.columns("plasmids")


def test_write_rejects_a_frame_that_breaks_the_contract(tmp_path):
    """Schema validation happens at write time, so a bad stage fails where it is, not three
    stages later in something that tried to join against it."""
    bad = pa.table({"plasmid_id": ["p1"]})
    with pytest.raises(ValueError):
        store.write_table(store.table_path(tmp_path, "plasmids"), "plasmids", bad)


def test_connect_exposes_each_written_table_as_a_queryable_view(tmp_path):
    store.write_table(store.table_path(tmp_path, "plasmids"), "plasmids", _plasmids_frame())
    with store.connect(tmp_path) as connection:
        assert connection.execute("SELECT count(*) FROM plasmids").fetchone()[0] == 1
