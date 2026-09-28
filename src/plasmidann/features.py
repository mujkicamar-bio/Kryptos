"""GenBank locations and GFF3 rows for CDS features, including origin-spanning ones.

An ORF reconstructed across the origin of a circular plasmid has start > end: it runs
start..length, then 1..end. GenBank writes it as join(start..length,1..end), wrapped whole
in complement() on the minus strand. GFF3 does not allow start > end, so the ORF becomes a
discontinuous feature: one row per segment, all rows sharing one ID.
"""

# Characters with a reserved meaning in GFF3 column 9, percent-encoded as the Sequence
# Ontology's GFF3 format definition (version 1.26) requires. A DIAMOND title is free text
# and often contains several of them.
_GFF3_ESCAPE = {"%": "%25", ";": "%3B", "=": "%3D", "&": "%26", ",": "%2C",
                "\t": "%09", "\n": "%0A", "\r": "%0D"}


def _escape(value):
    return "".join(_GFF3_ESCAPE.get(c, c) for c in str(value))


def gff3_attributes(attributes):
    """The attribute column, escaped, in the caller's order, with empty values omitted.

    An empty value is left out rather than written as `product=`: a dark ORF has no product.
    """
    return ";".join(f"{_escape(k)}={_escape(v)}"
                    for k, v in attributes.items() if v not in (None, ""))


def genbank_location(start, end, strand, length, partial_begin=False, partial_end=False):
    """A GenBank location string; `start > end` means the feature runs start..length, 1..end.

    The complement wrapper encloses the entire join: complement(join(a,b)) reads b then a,
    reverse-complemented, which is the gene; join(complement(a),complement(b)) is not.
    A feature that runs off the sequence at its start or end coordinate is marked with '<'
    or '>' there, whatever its strand (INSDC Feature Table Definition, location syntax).
    """
    first = f"<{start}" if partial_begin else str(start)
    last = f">{end}" if partial_end else str(end)
    if start <= end:
        span = f"{first}..{last}"
    else:
        span = f"join({first}..{length},1..{last})"
    return f"complement({span})" if strand in (-1, "-", "-1") else span


def gff3_features(gene, length, attributes=None):
    """One gene as GFF3 rows (9-tuples): one row, or two sharing one ID if it spans the origin.

    Phase is the number of bases to skip at a segment's 5' end before the next codon. The
    genes are complete CDS called from their own start codon, so the segment translated
    first has phase 0 and the other segment continues its reading frame. On the minus strand
    translation starts in the 1..end segment.
    """
    strand_char = "-" if gene["strand"] in (-1, "-", "-1") else "+"
    attrs = dict(attributes or {})
    attrs.setdefault("ID", gene["orf_id"])

    if gene["start"] <= gene["end"]:
        segments = [(gene["start"], gene["end"])]
    else:
        segments = [(gene["start"], length), (1, gene["end"])]
    first = segments[-1] if strand_char == "-" else segments[0]
    carried = (first[1] - first[0] + 1) % 3

    column9 = gff3_attributes(attrs)
    return [(gene["plasmid_id"], "plasmidann", "CDS", s, e, ".", strand_char,
             0 if (s, e) == first else (3 - carried) % 3, column9)
            for s, e in segments]
