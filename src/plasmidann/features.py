"""S4: turning the annotation table into the two feature formats.

The design lists GFF3 and GenBank as S4 deliverables and records that v1 never wrote them.
Formatting is the easy half. The hard half is the 160,375 ORFs that S1 reconstructed across
the origin of a circular plasmid: those carry start > end, and the two formats disagree
about how to express that.

  GenBank  has dedicated syntax - join(start..L,1..end) - and a strand wrapper that must
           enclose the WHOLE join, because complement(join(a,b)) and join(complement(a),
           complement(b)) describe different proteins.
  GFF3     forbids start > end outright. The spec's answer is a discontinuous feature:
           several lines carrying the same ID, which parsers reassemble into one gene.

Getting either wrong yields a file that loads without complaint and puts genes in the wrong
part of the molecule, which is the class of error nobody notices until an experiment fails.
"""

# Characters that carry structural meaning in a GFF3 attribute column. Escaped per the
# GFF3 specification, which requires percent-encoding for exactly these. A DIAMOND stitle
# is free text and routinely contains several of them.
_GFF3_ESCAPE = {"%": "%25", ";": "%3B", "=": "%3D", "&": "%26", ",": "%2C",
                "\t": "%09", "\n": "%0A", "\r": "%0D", "|": "%7C"}


def _escape(value):
    return "".join(_GFF3_ESCAPE.get(c, c) for c in str(value))


def gff3_attributes(attributes):
    """The attribute column, escaped, with empty values omitted.

    Omitting rather than writing `product=` keeps the file honest: a dark ORF has no
    product, and an empty attribute asserts that it has one which happens to be blank.

    Order is the caller's, because GFF3 conventionally leads with ID and readers scan for
    it; Python dictionaries preserve insertion order, so the caller controls this.
    """
    return ";".join(f"{_escape(k)}={_escape(v)}"
                    for k, v in attributes.items() if v not in (None, ""))


def genbank_location(start, end, strand, length):
    """A GenBank location string, honouring the origin-spanning convention.

    `start > end` means the feature runs start..length then 1..end, which is what S1 writes
    for a gene it reconstructed across the cut point of a circular plasmid.

    The complement wrapper encloses the entire join. complement(join(a,b)) reads the
    segments in the order b then a, reverse-complemented, which is the actual gene;
    join(complement(a),complement(b)) reads them a then b and yields a different protein.
    """
    if start <= end:
        span = f"{start}..{end}"
    else:
        span = f"join({start}..{length},1..{end})"
    return f"complement({span})" if strand in (-1, "-", "-1") else span


def gff3_features(gene, length, source="plasmidann", feature_type="CDS",
                  attributes=None):
    """One gene as one or more GFF3 rows, as 9-tuples ready to be tab-joined.

    An origin-spanning gene becomes TWO rows sharing one ID - a discontinuous feature. A
    single row with start > end is invalid GFF3: some parsers reject the file and others
    silently reinterpret the coordinates, which is worse.

    `phase` is 0 on every row. These are complete CDS features called by Pyrodigal from
    their own start codon, so the first base of the feature is the first base of a codon.
    For the second segment of an origin-spanning gene that is not strictly true, but GFF3
    has no way to express a phase carried across segments of a discontinuous feature and
    every parser recomputes it from the reassembled feature anyway.
    """
    strand_char = "-" if gene["strand"] in (-1, "-", "-1") else "+"
    attrs = dict(attributes or {})
    attrs.setdefault("ID", gene["orf_id"])

    if gene["start"] <= gene["end"]:
        segments = [(gene["start"], gene["end"])]
    else:
        segments = [(gene["start"], length), (1, gene["end"])]

    column9 = gff3_attributes(attrs)
    return [(gene["plasmid_id"], source, feature_type, s, e, ".", strand_char, "0",
             column9)
            for s, e in segments]
