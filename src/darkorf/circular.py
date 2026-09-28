"""Circular-origin repair.

A circular molecule written as a FASTA record is cut at an arbitrary point, and a gene
crossing that point is split into a fragment at each end. For a closed molecule the first
min(MAX_OVERLAP_BP, L) bases are appended to the record before gene calling, so a gene
crossing the cut is contiguous and is called whole. Every gene in the appended copy is
then called twice; resolve_origin_genes keeps one copy per gene and maps it back onto the
record, writing a gene that crosses position L with start > end. Topologies other than
CIRCULAR_TOPOLOGIES, including unknown ones, are treated as linear.
"""

# Longer than the coding sequence of almost every plasmid gene. A gene crossing the cut
# and longer than this is kept as a partial (see resolve_origin_genes). The cap bounds the
# extra calling cost on large plasmids.
MAX_OVERLAP_BP = 5000

# Topologies that denote a closed molecule. 'direct terminal repeat' is a circular molecule
# an assembler reported linearly with the overlap between its two ends left in the record:
# measured on 400 such records, every one begins with an exact copy of its own last >= 20
# bp, and in 76% the copy is not a multiple of 3 long, so one copy is removed before gene
# calling (terminal_repeat_length, rule analysis_set).
#
# 'inverted terminal repeat' is not listed: it marks a genuinely linear replicon with
# hairpin or protein-capped telomeres (the Borrelia and Streptomyces linear plasmids), and
# joining its ends would invent a gene across a junction that does not exist in the cell.
# Treating a circular molecule as linear loses the genes that cross the origin; treating a
# linear molecule as circular invents genes, so anything not listed is linear.
CIRCULAR_TOPOLOGIES = frozenset({
    "circular",
    "direct terminal repeat",
})


def terminal_repeat_length(seq, min_repeat):
    """Length of the longest exact repeat between the start and the end of a record, or 0.

    The overlap an assembler leaves when it reports a circular molecule linearly: the
    record's first k bases equal its last k. Only repeats of at least `min_repeat` bp count
    (20 bp, CheckV's direct-terminal-repeat criterion; Nayfach et al. 2021, Nat.
    Biotechnol. 39:578), and at most half the record, so the two copies cannot overlap.

    Linear in the record length: only positions in the second half where the first
    `min_repeat` bases recur are candidates, and the earliest one that matches through to
    the end is the longest repeat.
    """
    n = len(seq)
    if n < 2 * min_repeat:
        return 0
    seed = seq[:min_repeat]
    pos = seq.find(seed, n - n // 2)
    while pos != -1:
        if seq[pos:] == seq[:n - pos]:
            return n - pos
        pos = seq.find(seed, pos + 1)
    return 0


def overlap_for(length):
    """How many bases of the record's start to append to its end: the whole molecule, so
    that any gene crossing the cut is contiguous, up to MAX_OVERLAP_BP."""
    return min(MAX_OVERLAP_BP, length)


# The calls still depend somewhat on where the circle was cut: pyrodigal's meta mode picks
# one of its models per call, and the appended copy can change which model wins, and with
# it start sites or the translation table. Measured on the first 3,000 circular plasmids
# < 20 kb of the analysis set, each called at rotations 0, L/3 and 2L/3 after the terminal
# repeat is trimmed: 595 plasmids (19.8%) give a different protein set at some rotation,
# and 1,779 of 16,232 distinct proteins (11.0%) depend on the rotation.
def resolve_origin_genes(genes, original_length, extended_length):
    """Map genes called on the extended sequence back onto the record.

    `genes` are dicts with start, end (1-based inclusive on the extended sequence of
    `extended_length` bases), strand and partial. Returns new dicts, each with
    `origin_spanning`; a gene crossing the cut has start > end and is read as
    start..original_length followed by 1..end (the GenBank join convention).

    Every base of the appended copy lies in the record too, so a gene there can be called
    twice: once in the record and once in the copy, the two calls lying original_length
    apart in the same reading frame on the same strand. Of the two, the call farther from
    the ends of the extended sequence is kept, because near an end the caller sees only
    part of the gene's context and truncates the gene there (a partial call). This drops
    the head fragment of a gene called whole across the cut, keeps the whole copy of a gene
    whose record copy was truncated at position 1, and leaves a gene longer than the
    appended copy partial. A partial call with no twin is dropped: a circle has no end to
    run off, and such a call is one the caller made only because edge genes need no start
    codon. A complete call with no twin that starts in the appended copy is kept when no
    kept call overlaps it: the caller then found a gene there that its record calls miss,
    such as a gene near position 1 that the record calls only as a partial on the other
    strand. Otherwise the record's call stands. On 3,000 small circular plasmids of the
    analysis set this keeps 14 genes that would otherwise be lost, and no gene that
    overlaps another.
    """
    length = original_length

    def margin(g):
        return min(g["start"] - 1, extended_length - g["end"])

    def twins(a, b):
        # b shifted back by one molecule length overlaps a in the same frame and strand.
        start, end = b["start"] - length, b["end"] - length
        return (a is not b and a["strand"] == b["strand"] and (a["end"] - end) % 3 == 0
                and a["start"] <= end and start <= a["end"])

    twinned, dropped = set(), set()
    for j, b in enumerate(genes):
        if b["end"] <= length:
            continue
        for i, a in enumerate(genes):
            if twins(a, b):
                twinned |= {i, j}
                dropped.add(i if margin(b) > margin(a) else j)
    dropped |= {i for i, g in enumerate(genes) if i not in twinned and g["partial"]}
    in_copy = {i for i, g in enumerate(genes)
               if i not in twinned and i not in dropped and g["start"] > length}

    def overlaps_kept(c):
        # Against each other kept call, at its own coordinates and one molecule length back.
        start, end = c["start"] - length, c["end"] - length
        return any(g["start"] - shift <= end
                   and start <= g["end"] - shift
                   for i, g in enumerate(genes) if i not in dropped and i not in in_copy
                   for shift in (0, length))

    dropped |= {i for i in in_copy if overlaps_kept(genes[i])}

    kept = []
    for i, g in enumerate(genes):
        if i in dropped:
            continue
        out = dict(g, origin_spanning=g["start"] <= length < g["end"])
        if g["start"] > length:
            out["start"] -= length
        if g["end"] > length:
            out["end"] -= length
        kept.append(out)
    return kept


def is_circular(topology):
    """Whether a topology string denotes a closed molecule needing origin repair."""
    return (topology or "").strip().lower() in CIRCULAR_TOPOLOGIES
