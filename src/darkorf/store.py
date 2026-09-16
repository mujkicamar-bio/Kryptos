"""Parquet store with DuckDB query access (spec §4.2).

Why not one big TSV: the integrated occurrence view is roughly 9.3 million rows by ~120
columns, of which about 90 are constant within protein_id or family_id. Materializing that
as text is 9-14 GB and needs 50-110 GB of RAM to open. Parquet keeps each table at its own
level and DuckDB reconstructs the wide view on demand, without loading it (spec §63).
"""
import contextlib
import pathlib

import duckdb
import pyarrow.parquet as parquet

from darkorf import schemas


def table_path(outdir, table_name):
    """Canonical on-disk location of one table."""
    return pathlib.Path(outdir) / "tables" / f"{table_name}.parquet"


def write_table(path, table_name, frame):
    """Validate against the declared contract, then write Parquet.

    Validation is deliberately at write time: a stage that emits the wrong columns must
    fail in its own rule, where the log says which stage it was.
    """
    schemas.validate_frame(table_name, frame)
    path = pathlib.Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    # ZSTD because these tables are read far more often than written, and the sequence
    # columns compress by roughly an order of magnitude.
    parquet.write_table(frame, path, compression="zstd")


def read_table(path):
    return parquet.read_table(path)


@contextlib.contextmanager
def connect(outdir):
    """Open DuckDB with every Parquet table registered as a view of the same name.

    Views rather than imports: the Parquet files stay the single copy of the data, so a
    stage that rewrites a table is immediately visible to every later query.
    """
    connection = duckdb.connect()
    try:
        for path in sorted((pathlib.Path(outdir) / "tables").glob("*.parquet")):
            connection.execute(
                f"CREATE VIEW {path.stem} AS SELECT * FROM read_parquet('{path}')"
            )
        yield connection
    finally:
        connection.close()
