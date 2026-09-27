"""Streaming the plasmid FASTA.

The analysis-set FASTA written at S0 is what every stage that needs sequence reads, so
one reader serves them all and is tested once here.
"""
import gzip

import pytest

from plasmidann import fasta


def test_iter_fasta_reads_plain_and_gzipped_alike(tmp_path):
    """The working set is delivered gzipped and the smoke set is written plain. If the
    reader handled only one, the other would yield no sequence - silently."""
    plain = tmp_path / "a.fna"
    plain.write_text(">p1 some description\nACGT\nACGT\n>p2\nTTTT\n")
    zipped = tmp_path / "b.fna.gz"
    with gzip.open(zipped, "wt") as fh:
        fh.write(">p3\nGGGG\n")

    records = list(fasta.iter_fasta([plain, zipped]))
    assert records == [("p1", "ACGTACGT"), ("p2", "TTTT"), ("p3", "GGGG")]


def test_iter_fasta_takes_the_identifier_only_from_the_header(tmp_path):
    """Everything after the first whitespace is description, not identity.

    plasmid_id is the join key for every table in the run, so a description carried into
    it makes every lookup miss - and that failure presents as an unannotated plasmid
    rather than as an error.
    """
    path = tmp_path / "c.fna"
    path.write_text(">NZ_CP012345.1 Escherichia coli plasmid pX, complete sequence\nACGT\n")

    assert list(fasta.iter_fasta([path])) == [("NZ_CP012345.1", "ACGT")]


def test_iter_fasta_yields_nothing_for_an_empty_file(tmp_path):
    path = tmp_path / "d.fna"
    path.write_text("")

    assert list(fasta.iter_fasta([path])) == []


def test_split_fasta_keeps_every_record_once_and_balances_length(tmp_path):
    """Chunks for the one-replicon-at-a-time tools: every record exactly once, and no
    chunk far longer than the rest, or one process finishes hours after the others."""
    src = tmp_path / "in.fna"
    lengths = [900, 800, 100, 100, 100, 100, 100, 100]
    src.write_text("".join(f">p{i}\n{'A' * n}\n" for i, n in enumerate(lengths)))

    chunks = fasta.split_fasta(src, 3, tmp_path / "chunks")

    records = [r for c in chunks for r in fasta.iter_fasta([c])]
    assert sorted(name for name, _ in records) == sorted(f"p{i}" for i in range(8))
    totals = [sum(len(s) for _, s in fasta.iter_fasta([c])) for c in chunks]
    assert max(totals) <= 900
    # Fewer records than chunks: only the chunks that received one are returned.
    assert len(fasta.split_fasta(src, 20, tmp_path / "many")) == 8


def test_an_empty_header_names_the_file_and_line(tmp_path):
    path = tmp_path / "e.fna"
    path.write_text(">p1\nACGT\n>\nACGT\n")

    with pytest.raises(ValueError, match="e.fna line 3"):
        list(fasta.iter_fasta([path]))
