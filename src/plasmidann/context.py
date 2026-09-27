"""S8: the genomic neighbourhood of a gene - transcriptional units, flanking genes and the
islands (defence systems, integron arrays, IS elements) a gene lies in.

Genes are dicts with orf_id, start, end and strand, one plasmid at a time. A gene
reconstructed across the origin of a circular plasmid is written start > end (start..L
then 1..end, the GenBank join() convention). The same definitions serve S8 context and
S9 synteny.
"""

# Maximum intergenic distance, in nucleotides, for two consecutive same-strand genes to be
# treated as one transcriptional unit: the value FESNov used (Rodriguez del Rio et al. 2024,
# Nature 626:377).
DEFAULT_MAX_GAP = 100


def directons(genes, max_gap=DEFAULT_MAX_GAP, circular=False, length=None):
    """Putative transcriptional units: maximal runs of consecutive same-strand genes.

    Consecutive genes on the same strand separated by at most `max_gap` nucleotides share a
    unit; overlapping genes (a negative gap) do too, since overlapping start and stop codons
    are common in operons. Genes are sorted by start here. On a circular molecule the
    record's last and first units merge when they share a strand and the gap across the
    origin is at most `max_gap`; `length`, the molecule length, gives that gap, and the
    merged unit is listed first.

    Returns a list of lists of orf_id, in coordinate order.
    """
    if circular and length is None:
        raise ValueError("a circular directon needs the molecule length")
    ordered = sorted(genes, key=lambda g: (g["start"], g["end"]))
    units, current = [], []
    for gene in ordered:
        if not current:
            current = [gene]
            continue
        previous = current[-1]
        gap = gene["start"] - previous["end"] - 1
        if gene["strand"] == previous["strand"] and gap <= max_gap:
            current.append(gene)
        else:
            units.append(current)
            current = [gene]
    if current:
        units.append(current)
    if circular and len(units) > 1:
        last, first = units[-1][-1], units[0][0]
        if last["start"] > last["end"]:
            # The last gene already crosses the origin.
            gap = first["start"] - last["end"] - 1
        else:
            gap = (length - last["end"]) + (first["start"] - 1)
        if last["strand"] == first["strand"] and gap <= max_gap:
            units = [units[-1] + units[0]] + units[1:-1]
    return [[g["orf_id"] for g in unit] for unit in units]


def flanks(genes, window=3, circular=False):
    """{orf_id: (left, right)}: up to `window` flanking orf_ids on each side, nearest first.

    Left and right are in record coordinates. The window stops at the ends of a linear
    record and wraps across the origin of a circular one (on the test set 66% of
    small-plasmid ORFs had a window cut short at a record end). A gene is never its own
    neighbour, so on a circle of n genes each side holds at most n - 1, and on a small
    circle the two sides can share genes.
    """
    ids = [g["orf_id"] for g in sorted(genes, key=lambda g: (g["start"], g["end"]))]
    n = len(ids)
    out = {}
    for i, orf_id in enumerate(ids):
        if circular:
            k = min(window, n - 1)
            out[orf_id] = ([ids[(i - d) % n] for d in range(1, k + 1)],
                           [ids[(i + d) % n] for d in range(1, k + 1)])
        else:
            out[orf_id] = (ids[max(0, i - window):i][::-1], ids[i + 1:i + 1 + window])
    return out


def _segments(feature):
    """A feature as one or two plain intervals; an origin-spanning one is start..inf, 1..end.

    The open upper bound stands for the molecule end: no island lies beyond it, so the
    overlap answer is the same as with the true length.
    """
    start, end = feature["start"], feature["end"]
    if start <= end:
        return [(start, end)]
    return [(start, float("inf")), (1, end)]


def overlapping_islands(gene, islands):
    """Every island (dict with start and end) the gene overlaps, in the given order.

    Partial overlap counts, and all overlapping islands are returned, so a gene in a defence
    system inside an integron array is reported in both. Origin-spanning genes and islands
    are both handled.
    """
    hits = []
    for island in islands:
        if any(gs <= ie and ge >= is_
               for gs, ge in _segments(gene)
               for is_, ie in _segments(island)):
            hits.append(island)
    return hits
