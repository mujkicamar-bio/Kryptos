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


# ------------------------------------------------------------------------------------
# iter_fasta: the shards are the analysis scope, so every stage that needs sequence
# reads them rather than the corpus. Two stages do, which is why this lives here and not
# in either script.
# ------------------------------------------------------------------------------------

def test_iter_fasta_reads_plain_and_gzipped_shards_alike(tmp_path):
    """A shard directory may hold either, and SHARD_SUFFIXES already accepts both.

    If the reader handled only one, adding a compressed batch to an uncompressed
    collection would make the pipeline skip it silently - the file is discovered, so the
    shard exists, and it simply yields no sequence.
    """
    plain = tmp_path / "a.fna"
    plain.write_text(">p1 some description\nACGT\nACGT\n>p2\nTTTT\n")
    zipped = tmp_path / "b.fna.gz"
    with gzip.open(zipped, "wt") as fh:
        fh.write(">p3\nGGGG\n")

    records = list(shards.iter_fasta([plain, zipped]))
    assert records == [("p1", "ACGTACGT"), ("p2", "TTTT"), ("p3", "GGGG")]


def test_iter_fasta_takes_the_identifier_only_from_the_header(tmp_path):
    """Everything after the first whitespace is description, not identity.

    plasmid_id is the join key for every table in the run, so a description carried into
    it makes every lookup miss - and that failure presents as an unannotated plasmid
    rather than as an error.
    """
    path = tmp_path / "c.fna"
    path.write_text(">NZ_CP012345.1 Escherichia coli plasmid pX, complete sequence\nACGT\n")

    assert list(shards.iter_fasta([path])) == [("NZ_CP012345.1", "ACGT")]


def test_iter_fasta_yields_nothing_for_an_empty_shard(tmp_path):
    """An empty shard is not an error here. discover_shards already refuses a directory
    with no shards in it; a single empty file is a batch that happened to be filtered to
    nothing upstream, and the stages that read it should produce no rows rather than
    stop the run."""
    path = tmp_path / "d.fna"
    path.write_text("")

    assert list(shards.iter_fasta([path])) == []
