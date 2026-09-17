"""S4: the feature files.

The design names GFF3 and GenBank as S4 deliverables and records that v1 never wrote them.
The hard part is not the format, it is the 160,375 ORFs that S1 reconstructed across the
origin of a circular plasmid. Those carry start > end, and the two formats disagree about
how to say that - GFF3 forbids it outright, GenBank has a dedicated syntax for it. Getting
either wrong produces a file that loads without complaint and places genes in the wrong
part of the molecule.
"""
import pytest

from plasmidann.features import gff3_attributes, gff3_features, genbank_location


# --- GenBank: the join() convention S1 already writes ---------------------------------

def test_an_ordinary_gene_is_a_plain_range():
    assert genbank_location(100, 400, strand=1, length=5000) == "100..400"


def test_a_reverse_strand_gene_is_complemented():
    assert genbank_location(100, 400, strand=-1, length=5000) == "complement(100..400)"


def test_an_origin_spanning_gene_uses_join():
    """start > end means the gene runs start..L then 1..end - the convention S1 writes and
    the one GenBank defines. A naive `start..end` here would be rejected by every parser,
    and `end..start` would silently place the gene across the wrong 4,700 bases."""
    assert genbank_location(4900, 120, strand=1, length=5000) == "join(4900..5000,1..120)"


def test_an_origin_spanning_reverse_gene_complements_the_whole_join():
    """complement(join(a,b)), not join(complement(a),complement(b)) - the second reverses
    the order of the segments and yields a different protein."""
    assert genbank_location(4900, 120, strand=-1, length=5000) == (
        "complement(join(4900..5000,1..120))")


# --- GFF3: no join(), so an origin-spanning gene is two lines sharing one ID -----------

def test_an_ordinary_gene_is_one_gff3_feature():
    rows = gff3_features({"plasmid_id": "p1", "orf_id": "p1|3", "start": 100,
                          "end": 400, "strand": 1}, length=5000)
    assert len(rows) == 1
    assert rows[0][:5] == ("p1", "plasmidann", "CDS", 100, 400)
    assert rows[0][6] == "+"


def test_an_origin_spanning_gene_becomes_two_features_with_one_id():
    """GFF3 requires start <= end, so the join() cannot be written as one line. The spec's
    answer is a discontinuous feature: several lines carrying the SAME ID, which parsers
    reassemble into one gene. Emitting a single line with start > end produces a file that
    some tools reject and others silently reinterpret."""
    rows = gff3_features({"plasmid_id": "p1", "orf_id": "p1|9", "start": 4900,
                          "end": 120, "strand": 1}, length=5000)
    assert len(rows) == 2
    assert (rows[0][3], rows[0][4]) == (4900, 5000)
    assert (rows[1][3], rows[1][4]) == (1, 120)
    ids = {dict(kv.split("=", 1) for kv in r[8].split(";"))["ID"] for r in rows}
    assert len(ids) == 1, "the two segments must share one ID or they are two genes"
    assert all(r[3] <= r[4] for r in rows), "GFF3 requires start <= end"


# --- attributes: free text from nr goes into a structured field -----------------------

def test_attribute_values_are_escaped():
    """A DIAMOND stitle is free text and routinely contains ; = , and %, every one of which
    is a GFF3 attribute delimiter. Unescaped, one label silently becomes several attributes
    and the rest of the column is misparsed."""
    got = gff3_attributes({"ID": "p1|3",
                           "product": "protein A; subunit=2, 50% identity"})
    assert got == "ID=p1%7C3;product=protein A%3B subunit%3D2%2C 50%25 identity"


def test_an_empty_attribute_is_omitted_not_written_blank():
    """`product=` is a syntactically valid but meaningless attribute; a dark ORF has no
    product and the file should say nothing rather than say nothing loudly."""
    assert gff3_attributes({"ID": "x", "product": "", "note": None}) == "ID=x"
