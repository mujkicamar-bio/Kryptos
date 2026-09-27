"""GenBank and GFF3 encoding of CDS features, including origin-spanning ones."""

from plasmidann.features import genbank_location, gff3_attributes, gff3_features

# --- GenBank: join() for a gene across the origin ------------------------------------

def test_an_ordinary_gene_is_a_plain_range():
    assert genbank_location(100, 400, strand=1, length=5000) == "100..400"


def test_a_reverse_strand_gene_is_complemented():
    assert genbank_location(100, 400, strand=-1, length=5000) == "complement(100..400)"


def test_an_origin_spanning_gene_uses_join():
    """start > end means the gene runs start..L then 1..end. `end..start` would place the
    gene across the other 4,700 bases."""
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
    """GFF3 requires start <= end, so the join() is written as a discontinuous feature:
    several rows carrying the same ID, which parsers reassemble into one gene."""
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
    got = gff3_attributes({"ID": "p1_3",
                           "product": "protein A; subunit=2, 50% identity"})
    assert got == "ID=p1_3;product=protein A%3B subunit%3D2%2C 50%25 identity"


def test_an_empty_attribute_is_omitted_not_written_blank():
    """`product=` is a syntactically valid but meaningless attribute; a dark ORF has no
    product and the file should say nothing rather than say nothing loudly."""
    assert gff3_attributes({"ID": "x", "product": "", "note": None}) == "ID=x"


def test_the_second_segment_of_an_origin_spanning_gene_carries_the_reading_frame():
    """GFF3 phase is the number of bases to skip before the next codon. The segment
    translated first has phase 0; the other continues its frame."""
    def phases(strand):
        rows = gff3_features({"plasmid_id": "p", "orf_id": "p|1", "start": 9990, "end": 7,
                              "strand": strand}, length=10000)
        return [(r[3], r[4], r[7]) for r in rows]

    # + strand: 9990..10000 (11 nt) is translated first, so 1..7 starts 1 base into a codon.
    assert phases("+") == [(9990, 10000, 0), (1, 7, 1)]
    # - strand: 7..1 (7 nt) is translated first, so 10000..9990 skips 2 bases.
    assert phases("-") == [(9990, 10000, 2), (1, 7, 0)]


def test_a_pipe_in_an_identifier_is_not_escaped():
    """GFF3 reserves only tab, newline, carriage return, %, control characters and ; = & ,
    in column 9. The ORF id keeps the same spelling as in every other table."""
    assert gff3_attributes({"ID": "p1|3"}) == "ID=p1|3"
