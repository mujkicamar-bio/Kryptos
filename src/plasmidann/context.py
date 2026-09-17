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

THE STATISTICAL TRAP THIS MODULE EXISTS TO AVOID

On a 5 kb cryptic plasmid carrying six genes, a +/-3 neighbourhood IS the entire plasmid.
Everything co-occurs with everything, and a raw co-occurrence frequency would rank the
smallest plasmids as the most informative when they are the least - which would be
catastrophic here, since small cryptic plasmids are a stratum of interest.

So co-occurrence is always reported as ENRICHMENT over a corpus-wide background, and the
family-level statistic is CONSERVATION across members rather than a single instance.
"""

# Maximum intergenic distance, in nucleotides, for two consecutive same-strand genes to be
# treated as one transcriptional unit. 100 nt is the value FESNov used; bacterial operon
# members are typically separated by less, and beyond it a shared promoter is unlikely.
DEFAULT_MAX_GAP = 100


def directons(genes, max_gap=DEFAULT_MAX_GAP):
    """Group genes into putative transcriptional units.

    A directon is a maximal run of consecutive genes on the same strand separated by less
    than `max_gap` nucleotides. Returns a list of lists of orf_id, in coordinate order.

    Membership of a directon is a much stronger contextual claim than mere adjacency: the
    genes are predicted to be co-transcribed, and therefore functionally coupled. A dark
    ORF sitting inside an otherwise fully annotated operon inherits that operon's
    hypothesis directly.

    Overlapping genes give a negative gap and stay in the same unit - overlapping start and
    stop codons are common in tightly packed operons and are evidence of coupling, not
    against it.

    Genes are sorted by start here rather than assumed sorted, because callers assemble
    them from a table that may be grouped by plasmid but not ordered within it.
    """
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
    return [[g["orf_id"] for g in unit] for unit in units]


def neighbourhood(genes, orf_id, window=3):
    """The orf_ids of the `window` genes on either side of `orf_id`.

    Truncated at the ends of the record rather than wrapped. Wrapping would be correct for
    a circular molecule, but it must be an explicit decision made by the caller who knows
    the topology, not an accident of negative list indexing - which is how it would happen
    silently in Python.
    """
    ordered = sorted(genes, key=lambda g: (g["start"], g["end"]))
    ids = [g["orf_id"] for g in ordered]
    try:
        i = ids.index(orf_id)
    except ValueError:
        return []
    left = ids[max(0, i - window):i]
    right = ids[i + 1:i + 1 + window]
    return left + right


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
    The defence_island stratum - 175 of the 1,000 constructs - was unreachable behind it.

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


def background_rate(all_orfs, feature):
    """How often `feature` appears in the context of ANY protein in the corpus.

    The denominator for enrichment. Without it, a feature that is simply common - as
    transposases are on plasmids - looks like a discovery every time.
    """
    if not all_orfs:
        return 0.0
    n = sum(1 for o in all_orfs if feature in o.get("context", ()))
    return round(n / len(all_orfs), 6)


def enrichment(observed, background):
    """Observed co-occurrence relative to the corpus background.

    1.0 means the association is exactly what chance predicts. This is the number that
    protects against the small-plasmid trap: on a six-gene plasmid the observed rate is
    high for everything, and so is the background, so the ratio stays near 1.

    A zero background with a non-zero observation is infinitely enriched, which is the
    honest answer - the caller decides how to rank it. Zero over zero is 1.0: no signal,
    no surprise.
    """
    if background == 0:
        return 1.0 if observed == 0 else float("inf")
    return observed / background
