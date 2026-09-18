"""Stage 9: synteny and context conservation (spec section 42).

WHAT THIS MEASURES, AND WHY IT IS NOT THE SAME AS STAGE 8

Stage 8 asks what a dark ORF sits NEXT TO, once, on one plasmid. Stage 9 asks whether the
same arrangement recurs across the family's occurrences. The spec's own example:

    A - B - DARK - C - D
    A - B - DARK - C - D
    A - B - DARK - C
    A - B - DARK - C - D

That is not four observations of "near B" - it is one conserved gene order seen four times.
Conserved gene order is a far stronger claim than adjacency: it survives because the
arrangement matters, so a dark ORF locked between B and C across independent plasmids is
predicted to be functionally coupled to them.

SIX MEASUREMENTS, KEPT APART (section 42)

    left_neighbor_conservation    the most common immediate left neighbour's frequency
    right_neighbor_conservation   the same on the right
    neighborhood_conservation     the most common UNORDERED neighbour set
    operon_like_conservation      same strand AND close, across occurrences
    synteny_conservation          the most common ORDERED left-right pair
    context_recurrence            how many occurrences had any usable context at all

They are separate because they fail separately. A dark ORF with a conserved left neighbour
and a variable right one is a real and common arrangement - the left gene may be its
promoter-sharing partner - and a single averaged "context score" would hide it.

Section 42's closing line governs the whole stage: "these measurements remain separate from
annotation." A conserved neighbour is evidence about a protein; it is not a name for it.

THE UNIT IS THE OCCURRENCE, AND THAT IS A KNOWN LIMITATION

Conservation here is computed over the family's occurrences, so a family on forty copies of
one redeposited plasmid will show perfect synteny from a single biological event. The
independent-lineage count from Stage 7 is what a reader must check alongside it, and this
module reports `n_occurrences` so the two can never be confused. Weighting by lineage is a
better statistic and is not done here; doing it silently would change what the number means
without saying so.
"""
import collections

# What is reported when a family has too few usable occurrences to compare. Conservation
# over one occurrence is 1.0 by construction and means nothing, so it must not be emitted
# as a measurement - spec section 2.9, and the same reasoning as TOO_FEW_MEMBERS elsewhere.
MIN_OCCURRENCES = 2


def neighbour_pairs(context):
    """(left, right) for one occurrence, from its ordered neighbour list.

    `context` is (left_neighbours, right_neighbours), each already in genomic order and
    nearest-first. Missing neighbours - the ORF sits at the end of a record - are the empty
    string rather than being dropped, so "no left neighbour" stays distinguishable from
    "left neighbour we could not name".
    """
    left, right = context
    return (left[0] if left else "", right[0] if right else "")


def _modal_fraction(values):
    """Frequency of the most common non-empty value, and what it was.

    Empty values are excluded from the NUMERATOR and the DENOMINATOR both. An ORF at the
    end of a contig has no left neighbour, and counting that as a failure to conserve one
    would penalise a family for where the assembler cut, not for its biology.
    """
    present = [v for v in values if v]
    if not present:
        return 0.0, ""
    value, count = collections.Counter(present).most_common(1)[0]
    return round(count / len(present), 4), value


def conservation(occurrences):
    """Every section 42 measurement for one family.

    `occurrences` is a list of dicts, one per occurrence of the family, each with:

        left        ordered left neighbour labels, nearest first
        right       ordered right neighbour labels, nearest first
        operon      True when the ORF shares a transcriptional unit with a named neighbour

    Returns a dict with the six measurements, the modal values behind them, the occurrence
    count they were computed over, and a status.
    """
    usable = [o for o in occurrences if o.get("left") or o.get("right")]
    result = {
        "n_occurrences": len(occurrences),
        "context_recurrence": len(usable),
        "left_neighbor_conservation": "",
        "right_neighbor_conservation": "",
        "neighborhood_conservation": "",
        "operon_like_conservation": "",
        "synteny_conservation": "",
        "modal_left": "",
        "modal_right": "",
        "modal_synteny": "",
    }
    if len(usable) < MIN_OCCURRENCES:
        # One occurrence is perfectly conserved with itself. Emitting 1.0 would put a
        # singleton at the top of any ranking of conserved context.
        result["status"] = "TOO_FEW_MEMBERS"
        return result

    lefts = [(o.get("left") or [""])[0] for o in usable]
    rights = [(o.get("right") or [""])[0] for o in usable]

    result["left_neighbor_conservation"], result["modal_left"] = _modal_fraction(lefts)
    result["right_neighbor_conservation"], result["modal_right"] = _modal_fraction(rights)

    # UNORDERED: the same two genes flanking the ORF count as the same neighbourhood even
    # when the plasmid is written in the opposite orientation, which is a presentation
    # artefact rather than a biological difference.
    neighbourhoods = [
        ",".join(sorted({n for n in list(o.get("left") or []) + list(o.get("right") or [])
                         if n}))
        for o in usable
    ]
    result["neighborhood_conservation"], _ = _modal_fraction(neighbourhoods)

    # ORDERED: left-then-right as one string. This is the spec's A-B-DARK-C-D pattern and
    # is strictly stronger than the unordered version above.
    syntenies = [f"{left}|{right}" for left, right in zip(lefts, rights)
                 if left or right]
    fraction, modal = _modal_fraction(syntenies)
    result["synteny_conservation"], result["modal_synteny"] = fraction, modal

    operon = [bool(o.get("operon")) for o in usable]
    result["operon_like_conservation"] = round(sum(operon) / len(operon), 4)

    result["status"] = "SUCCESS"
    return result
