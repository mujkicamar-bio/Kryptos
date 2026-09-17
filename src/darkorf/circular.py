"""Circular-origin repair: making gene calling independent of where the circle was cut.

THE PROBLEM

A plasmid is physically a circle. A FASTA record is a line. To write the circle as a line
you cut it at an arbitrary point, and any gene spanning that cut is split into two
fragments - one at the start of the record, one at the end.

Measured on the current analysis set:

    160,375 partial ORFs        1.72% of 9,317,050, i.e. 1.12 per plasmid
     88,600 of them start at coordinate <= 3
    134,748 of 143,503 plasmids (93.9%) are effectively circular
                                  75,231 'circular' + 59,517 'direct terminal repeat'

1.12 per plasmid is what you would predict when a random cut point in a gene-dense
circular molecule usually lands inside a gene.

WHY IT MATTERS HERE SPECIFICALLY

These are REAL genes, artificially truncated - which is the opposite of the artefact class
AntiFam catches, and needs the opposite treatment. A truncated protein aligns to only part
of any domain, so its coverage drops, pushing it toward DOMAIN_ONLY; if the cut removes
the domain entirely it gets no hit at all and is classed dark. Either way, half of a
perfectly well-understood RepA can enter the screening pool as a novel dark protein - and
synthesising half a protein guarantees a dead well.

THE APPROACH

Append the first `overlap` bases to the end, call genes on the extended sequence, then map
coordinates back. A gene that straddled the cut is now wholly inside the extended sequence
and is called intact, with the correct translation.

    original    [1 .................................... L]
    extended    [1 .................................... L][1 ... overlap]
                                              ^gene now contiguous^

An analogy with text makes the mechanism concrete. Consider a sentence written around a
ring, so that it has no first or last word. Transcribing the ring onto a line requires a
cut, and the cut usually falls inside a word:

    transcribed     ELLO WORLD HOW ARE YOU H
                                           ^ the cut

A reader of the line sees two fragments, "ELLO" and "H", and has no way to tell that they
are two halves of one word. Continue transcribing past the cut, repeating the first few
words, and the word reappears whole:

    extended        ELLO WORLD HOW ARE YOU H ELLO WORLD
                                           ^ HELLO is now contiguous

Read the words off the extended line and then discard what the repetition introduced:
"ELLO" at the start, because it is the head of a word already read whole, and the second
"WORLD", because it is a copy of the first. The ring is the plasmid, the words are genes,
the fragments are the two partial ORFs that a linear caller reports, and the pipeline
never joins the fragments; it discards them and reads the gene whole from the extended
sequence.

Four cases when mapping back, handled by resolve_origin_genes:

    end <= L            an ordinary gene            keep unchanged
    start <= 3, partial truncated at the record   the head fragment of a gene now called
                        start                      intact across the cut; drop
    start <= L < end    the gene crossed the cut    keep, wrap end, flag origin_spanning
    start > L           wholly in the appended tail duplicate of one already kept; drop

THE INVARIANT (success criterion SC6)

Calling genes on any rotation of a circular sequence must yield the same protein set. If it
does not, the coordinate system is contributing biology, which it must never do. Tested in
tests/test_circular.py.
"""

# Long enough to contain any plausible plasmid gene - the largest known plasmid proteins
# are well under 5 kb of coding sequence - and short enough that the extra calling cost is
# negligible against a 143,503-plasmid set.
MAX_OVERLAP_BP = 5000

# Topologies that denote a closed molecule. 'direct terminal repeat' is the signature of a
# circular molecule that an assembler resolved and reported linearly; measured, those
# records have their repeats already trimmed (0.9% intra-plasmid duplicate rate against
# 32.2% for 'circular', which is genuine multi-copy IS biology).
#
# 'inverted terminal repeat' is deliberately NOT here. An ITR is the signature of a
# genuinely linear replicon with hairpin or protein-capped telomeres - the Borrelia and
# Streptomyces linear plasmids, phi29, adenovirus - so joining its ends would fabricate a
# gene across a junction that does not exist in the cell.
#
# Anything not listed, including a missing or unrecognised value, is treated as linear. The
# two errors are not symmetric: calling a circular molecule linear loses the ~1.7% of genes
# that cross the origin, while calling a linear molecule circular invents genes outright.
# Losing real data is recoverable; inventing it is not.
CIRCULAR_TOPOLOGIES = frozenset({
    "circular",
    "direct terminal repeat",
})


def overlap_for(length):
    """How much sequence to append to the end, for a molecule of this length.

    Capped at MAX_OVERLAP_BP, and never more than half the molecule. The half limit is not
    cosmetic: appending more than half would let a single short gene be called three times
    - once in place, twice in the tail - which resolve_origin_genes deduplicates only for
    the tail region as a whole, not for repeated copies within it.
    """
    return min(MAX_OVERLAP_BP, length // 2)


def rotate(seq, offset):
    """Rotate a circular sequence so that `offset` becomes the new origin.

    Used only for testing the invariant; the pipeline never rotates real data. Kept here
    rather than in the test file so that the definition of "rotation" the invariant is
    stated against lives beside the code it constrains.
    """
    if not seq:
        return seq
    offset %= len(seq)
    return seq[offset:] + seq[:offset]


def resolve_origin_genes(genes, original_length):
    """Map genes called on the extended sequence back onto the original coordinates.

    `genes` are dicts with at least start, end (1-based inclusive, on the EXTENDED
    sequence), strand, partial and the reconstructed CDS. Returns a new list; inputs are
    not mutated.

    Every returned gene carries `origin_spanning`: True if it crosses the cut point, else
    False. For such a gene, start > end, and the feature is read as start..original_length
    followed by 1..end - the GenBank join() convention. Coordinates are reported unrotated,
    because ids.occurrence_id is built from them; normalising them into a linear-looking
    range would change the identifier. Downstream writers (GFF3, GenBank) must honour the
    convention; a naive end - start length calculation would be negative and is a bug.

    Genes lying wholly beyond original_length are duplicates of genes already called near
    the start of the record and are dropped. Without this, every circular plasmid would
    gain phantom ORFs equal to whatever fits in the appended tail.

    Genes truncated at the start of the record are dropped for the mirror reason. The
    caller only ever flags a gene partial at coordinate 1-3 when it ran off the left edge,
    and on a closed molecule that edge is not real: the gene is the one the extension has
    just called intact across the cut, and keeping both would emit every origin-spanning
    protein twice, once whole and once as its stub.
    """
    kept = []
    for g in genes:
        start, end = g["start"], g["end"]

        if start > original_length:
            # Wholly inside the appended tail: the same gene was already called at the
            # start of the record.
            continue
        if start <= 3 and g.get("partial"):
            # Truncated at the record start: the head fragment of a gene that the
            # extension calls whole across the origin.
            continue

        out = dict(g)
        if end > original_length:
            # Crossed the cut. The translation from the extended sequence is the real,
            # intact protein - which is the entire point of the exercise.
            out["end"] = ((end - 1) % original_length) + 1
            # orf_occurrences.origin_spanning is a declared column (spec 8.4, 9.4), so the
            # flag is carried explicitly on every gene rather than inferred downstream from
            # end < start.
            out["origin_spanning"] = True
            # A gene reconstructed across the origin is by definition no longer truncated
            # by the coordinate system.
            out["partial"] = 0
        else:
            out["origin_spanning"] = False
        kept.append(out)
    return kept


def is_circular(topology):
    """Whether a topology string denotes a closed molecule needing origin repair.

    Unknown or missing topology is treated as NOT circular. Extending a genuinely linear
    molecule would fabricate a junction that does not exist and could invent a chimeric
    gene across the two ends - a worse failure than leaving a real gene truncated.
    """
    return (topology or "").strip().lower() in CIRCULAR_TOPOLOGIES
