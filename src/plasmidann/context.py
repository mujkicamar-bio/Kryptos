"""S8: genomic context - turning a dark ORF into a testable hypothesis.

WHY THIS IS THE STAGE THAT MATTERS

Both experimentally validated discoveries in the FESNov study (Nature 626:377) were
selected by genomic context, not by any novelty measure: one sat in the canonical *che*
operon and was confirmed by a swimming assay, the other sat beside antibiotic resistance
genes and turned out to be a 36-residue antimicrobial peptide.

You cannot screen a protein for "function". You screen it for an activity. Context is what
supplies the specific prediction to test, which is why a ranked list of dark proteins with
no hypothesis attached is not a screening design.

WHY PLASMIDS ARE ESPECIALLY WELL SUITED

They are small, gene-dense and modular, and their cargo is organised into recognisable
islands. A dark ORF inside a defence island is a defence-system candidate; one adjacent to
a toxin-antitoxin pair is a candidate antitoxin; one inside an integron cassette array is
a real gene by construction, because it carries an attC site and has been physically
mobilised and retained under selection.

THE SMALL-PLASMID CAVEAT

On a 5 kb cryptic plasmid carrying six genes, a +/-3 neighbourhood IS the entire plasmid,
so everything co-occurs with everything. The pipeline reports context as descriptive rates
with no background correction (the enrichment test was removed on 2026-09-25), so a high
neighbour rate on small plasmids is expected and is not in itself a signal.
"""

# Maximum intergenic distance, in nucleotides, for two consecutive same-strand genes to be
# treated as one transcriptional unit. 100 nt is the value FESNov used; bacterial operon
# members are typically separated by less, and beyond it a shared promoter is unlikely.
DEFAULT_MAX_GAP = 100


def directons(genes, max_gap=DEFAULT_MAX_GAP, circular=False, length=None):
    """Group genes into putative transcriptional units.

    A directon is a maximal run of consecutive genes on the same strand separated by at
    most `max_gap` nucleotides (<= 100 nt by default). Returns a list of lists of orf_id,
    in coordinate order.

    Membership of a directon is a much stronger contextual claim than mere adjacency: the
    genes are predicted to be co-transcribed, and therefore functionally coupled. A dark
    ORF sitting inside an otherwise fully annotated operon inherits that operon's
    hypothesis directly.

    Overlapping genes give a negative gap and stay in the same unit - overlapping start and
    stop codons are common in tightly packed operons and are evidence of coupling, not
    against it.

    Genes are sorted by start here rather than assumed sorted, because callers assemble
    them from a table that may be grouped by plasmid but not ordered within it.

    On a circular molecule (`circular`, decided by the caller who knows the topology) the
    record's last run and first run are one unit when they share a strand and the gap
    across the origin is at most `max_gap`; the merged unit is listed first, in molecule
    order. `length` is the length of the sequence the genes were called on, which the gap
    across the origin needs. Without this, an operon that the linearisation of the deposit
    happened to cut was reported as two.
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
            # The last gene already crosses the origin (start > end, see _segments).
            gap = first["start"] - last["end"] - 1
        else:
            gap = (length - last["end"]) + (first["start"] - 1)
        if last["strand"] == first["strand"] and gap <= max_gap:
            units = [units[-1] + units[0]] + units[1:-1]
    return [[g["orf_id"] for g in unit] for unit in units]


def flanks(genes, orf_id, window=3, circular=False):
    """(left, right): the orf_ids of up to `window` genes on each side, nearest first.

    Truncated at the ends of a linear record. On a circular one (`circular`, decided by the
    caller who knows the topology) the window wraps across the origin, so a gene near the
    record's start has neighbours where the molecule has them rather than where the deposit
    happened to be linearised: on the test set 66% of small-plasmid ORFs had a window cut
    short at a record end. A gene is never its own neighbour, so on a circle of n genes
    each side holds at most n - 1, and on a small circle the two sides can share genes -
    as they do on the molecule.
    """
    ordered = sorted(genes, key=lambda g: (g["start"], g["end"]))
    ids = [g["orf_id"] for g in ordered]
    try:
        i = ids.index(orf_id)
    except ValueError:
        return [], []
    if circular:
        k = min(window, len(ids) - 1)
        return ([ids[(i - d) % len(ids)] for d in range(1, k + 1)],
                [ids[(i + d) % len(ids)] for d in range(1, k + 1)])
    return ids[max(0, i - window):i][::-1], ids[i + 1:i + 1 + window]


def neighbourhood(genes, orf_id, window=3, circular=False):
    """The orf_ids of the `window` genes on either side of `orf_id` (see flanks)."""
    left, right = flanks(genes, orf_id, window, circular)
    return left[::-1] + right


def _segments(feature):
    """A feature as one or two plain intervals, honouring the origin-spanning convention.

    S1 reconstructs genes that linearising a circular plasmid had broken in two, and writes
    them the way GenBank does: start > end, read as start..L followed by 1..end. 160,375
    ORFs in this collection are in that state, on a set where 94% of records are circular.

    The first segment runs to the end of the molecule, which is written here as an open
    upper bound rather than passing the length in. That is exact for this use: every island
    is an annotation ON the molecule, so no island coordinate can lie beyond its end, and an
    open bound and the true length therefore give the same overlap answer.
    """
    start, end = feature["start"], feature["end"]
    if start <= end:
        return [(start, end)]
    return [(start, float("inf")), (1, end)]


def overlapping_islands(gene, islands):
    """Every island a gene overlaps, in the order the islands were given.

    Islands come from DefenseFinder, IntegronFinder, CARD and TnCentral as intervals on the
    plasmid. Partial overlap counts: a gene straddling an island boundary is part of the
    module, and cassette and system boundaries are themselves predictions with error bars.

    ALL of them, not the first. A predecessor returned the first match and stopped, and
    since defence intervals are appended after integron intervals, a dark ORF inside a
    defence system that also sat in a cassette array could only ever be labelled 'integron'.
    A dark ORF in a defence system could never be reported as such behind it.

    Origin-spanning genes and origin-spanning islands are both handled, because either can
    be written start > end.
    """
    hits = []
    for island in islands:
        if any(gs <= ie and ge >= is_
               for gs, ge in _segments(gene)
               for is_, ie in _segments(island)):
            hits.append(island)
    return hits


def overlapping_island(gene, islands):
    """The first island a gene overlaps, or None.

    Kept for callers that genuinely want one answer. Prefer overlapping_islands: a gene can
    sit in a defence system inside a cassette array, and which of those a caller sees should
    not depend on the order the intervals happened to be appended in.
    """
    hits = overlapping_islands(gene, islands)
    return hits[0] if hits else None


def context_conservation(members, feature):
    """Fraction of a family's members whose context includes `feature`.

    This is the family-level statistic, and it is the one that matters. A single member
    sitting beside a defence system is a coincidence; eighty per cent of members doing so
    across independent plasmids is a hypothesis. FESNov required 90% for its
    high-confidence set, which yielded 4,349 families out of 404,085.
    """
    if not members:
        return 0.0
    n = sum(1 for m in members if feature in m.get("context", ()))
    return round(n / len(members), 4)

