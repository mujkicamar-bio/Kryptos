"""Discovering input shards from a directory rather than generating them.

Shards are prepared outside the pipeline and added one at a time. The pipeline reads
whatever is in the shard directory, so that a new batch is one more file rather than a
re-partition of the whole collection.
"""
import gzip

import pytest

from plasmidann import shards


def test_shards_are_discovered_from_the_directory(tmp_path):
    (tmp_path / "batch01.fna").write_text(">p1\nATG\n")
    (tmp_path / "batch02.fna").write_text(">p2\nATG\n")

    found = shards.discover_shards(tmp_path)

    assert sorted(found) == ["batch01", "batch02"]
    assert found["batch01"].endswith("batch01.fna")


def test_a_gzipped_shard_is_found_and_named_without_both_extensions(tmp_path):
    """'.fna.gz' must be stripped whole. Leaving a stranded '.fna' in the name would make
    the shard's outputs land under a different path from the plain-text equivalent."""
    with gzip.open(tmp_path / "batch03.fna.gz", "wt") as fh:
        fh.write(">p3\nATG\n")

    assert list(shards.discover_shards(tmp_path)) == ["batch03"]


def test_the_order_is_sorted_not_filesystem_order(tmp_path):
    """The DAG must be identical between runs on the same inputs. Filesystem order would
    change the job order, and every log that records it, for no biological reason."""
    for name in ("c", "a", "b"):
        (tmp_path / f"{name}.fna").write_text(">p\nATG\n")

    assert list(shards.discover_shards(tmp_path)) == ["a", "b", "c"]


def test_a_non_fasta_file_is_ignored(tmp_path):
    """Shard directories collect README files, checksums and job scripts. Treating one as
    a shard would send it to the gene caller."""
    (tmp_path / "batch01.fna").write_text(">p1\nATG\n")
    (tmp_path / "README.md").write_text("notes")
    (tmp_path / "batch01.fna.md5").write_text("abc")

    assert list(shards.discover_shards(tmp_path)) == ["batch01"]


def test_an_empty_directory_is_an_error_not_an_empty_run(tmp_path):
    """With no shards every aggregating rule has an empty input, every stage succeeds, and
    the run writes well-formed empty tables. Four stages of this pipeline once did exactly
    that and reported success."""
    with pytest.raises(SystemExit, match="no FASTA shards"):
        shards.discover_shards(tmp_path)


def test_a_missing_directory_is_an_error(tmp_path):
    with pytest.raises(SystemExit, match="does not exist"):
        shards.discover_shards(tmp_path / "absent")


def test_two_files_giving_one_shard_name_is_an_error(tmp_path):
    """batch01.fna and batch01.fna.gz both name shard 'batch01', and their outputs would
    overwrite each other - silently, because the second run would look like a resume."""
    (tmp_path / "batch01.fna").write_text(">p1\nATG\n")
    with gzip.open(tmp_path / "batch01.fna.gz", "wt") as fh:
        fh.write(">p1\nATG\n")

    with pytest.raises(SystemExit, match="both give the shard name"):
        shards.discover_shards(tmp_path)


def test_a_name_that_cannot_be_a_wildcard_is_an_error(tmp_path):
    """A shard name becomes a Snakemake wildcard and a path component. A space in it either
    breaks the DAG or silently changes where an output lands."""
    (tmp_path / "batch 01.fna").write_text(">p1\nATG\n")

    with pytest.raises(SystemExit, match="cannot"):
        shards.discover_shards(tmp_path)
