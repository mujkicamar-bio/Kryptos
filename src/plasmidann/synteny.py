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

    lineage_left_conservation          the most common immediate left neighbour's share
    lineage_right_conservation         the same on the right
    lineage_neighborhood_conservation  the most common UNORDERED neighbour set
    lineage_operon_like_conservation   same strand AND close, across lineages
    lineage_synteny_conservation       the most common ORDERED left-right pair
    context_recurrence                 how many occurrences had any usable context at all

They are separate because they fail separately. A dark ORF with a conserved left neighbour
and a variable right one is a real and common arrangement - the left gene may be its
promoter-sharing partner - and a single averaged "context score" would hide it.

Section 42's closing line governs the whole stage: "these measurements remain separate from
annotation." A conserved neighbour is evidence about a protein; it is not a name for it.

THE UNIT IS THE PLASMID LINEAGE, NOT THE OCCURRENCE

Counted over occurrences, a family on forty copies of one redeposited plasmid showed perfect
synteny from a single biological event. Each occurrence therefore carries the Stage 6
lineage of its plasmid (Mash distance <= 0.05, single linkage), and a measurement is a
FRACTIONAL VOTE over lineages:

    1. Empty values (no neighbour on that side) are dropped; a lineage whose values are all
       empty does not vote on that measurement.
    2. Each voting lineage has weight 1, split equally over its k non-empty values (1/k).
    3. A value's score is the sum of the weights it receives. The conservation is the
       largest score divided by the number of voting lineages; the modal value is the
       value with that score, ties broken by the lexicographically smallest value.

Copies of one lineage share one vote. A lineage whose copies disagree (one A|B, one C|D)
gives 1/2 to each and is counted in n_lineages_discordant, so a within-lineage
rearrangement is visible rather than averaged away. When every lineage holds one copy the
values are exactly the occurrence statistic this stage reported before, so the reading of
the numbers carries over for the common case. Scores are summed as exact fractions: ten
copies worth 1/10 each are one vote, not 0.9999999999999999.

`n_occurrences` and `context_recurrence` still count occurrences, beside `n_lineages`, the
lineages with usable context. That can be fewer than Stage 7's lineage count for the whole
family, which counts every member and not only the dark occurrences with context.
"""
import collections
from fractions import Fraction

# Fewer voting lineages than this gives TOO_FEW_LINEAGES and no value. Two is the arithmetic
# minimum for a comparison, not a tuned threshold: one lineage is conserved with itself by
# construction (spec section 2.9). The pipeline passes `synteny.min_lineages` from config.
MIN_LINEAGES = 2


def neighbour_pairs(context):
    """(left, right) for one occurrence, from its ordered neighbour list.

    `context` is (left_neighbours, right_neighbours), each already in genomic order and
    nearest-first. Missing neighbours - the ORF sits at the end of a record - are the empty
    string rather than being dropped, so "no left neighbour" stays distinguishable from
    "left neighbour we could not name".
    """
    left, right = context
    return (left[0] if left else "", right[0] if right else "")


def _lineage_modal(pairs):
    """Fractional vote over lineages: (conservation, modal value), or ("", "") with no vote.

    `pairs` is (lineage, value) per occurrence. Empty values are excluded from the vote
    entirely - an ORF at the end of a contig has no left neighbour, and counting that as a
    failure to conserve one would penalise a family for where the assembler cut, not for
    its biology.
    """
    values_of = collections.defaultdict(list)
    for lineage, value in pairs:
        if value:
            values_of[lineage].append(value)
    if not values_of:
        return "", ""
    score = collections.defaultdict(Fraction)
    for values in values_of.values():
        for value in values:
            score[value] += Fraction(1, len(values))
    modal = min(score, key=lambda v: (-score[v], v))
    return round(float(score[modal] / len(values_of)), 4), modal


def conservation(occurrences, min_lineages=MIN_LINEAGES):
    """Every section 42 measurement for one family, counted over lineages.

    `occurrences` is a list of dicts, one per occurrence of the family, each with:

        left        ordered left neighbour labels, nearest first
        right       ordered right neighbour labels, nearest first
        operon      True when the ORF shares a transcriptional unit with a named neighbour
        lineage     the Stage 6 lineage of the occurrence's plasmid

    Returns a dict with the measurements, the modal values behind them, the occurrence and
    lineage counts they were computed over, and a status.
    """
    usable = [o for o in occurrences if o.get("left") or o.get("right")]
    lefts = [(o["lineage"], (o.get("left") or [""])[0]) for o in usable]
    rights = [(o["lineage"], (o.get("right") or [""])[0]) for o in usable]
    # ORDERED: left-then-right as one string. This is the spec's A-B-DARK-C-D pattern and
    # is strictly stronger than the unordered neighbourhood below.
    syntenies = [(lin, f"{left}|{right}" if left or right else "")
                 for (lin, left), (_, right) in zip(lefts, rights)]

    arrangements = collections.defaultdict(set)
    for lineage, value in syntenies:
        if value:
            arrangements[lineage].add(value)
    n_lineages = len({o["lineage"] for o in usable})
    result = {
        "n_occurrences": len(occurrences),
        "context_recurrence": len(usable),
        "n_lineages": n_lineages,
        "n_lineages_discordant": sum(len(v) > 1 for v in arrangements.values()),
        "lineage_left_conservation": "",
        "lineage_right_conservation": "",
        "lineage_neighborhood_conservation": "",
        "lineage_operon_like_conservation": "",
        "lineage_synteny_conservation": "",
        "modal_left": "",
        "modal_right": "",
        "modal_synteny": "",
    }
    if not usable:
        result["status"] = "NO_CONTEXT"
        return result
    if n_lineages < min_lineages:
        # One lineage is perfectly conserved with itself, however many copies it has.
        # Emitting 1.0 would put a redeposited singleton at the top of any ranking.
        result["status"] = "TOO_FEW_LINEAGES"
        return result

    result["lineage_left_conservation"], result["modal_left"] = _lineage_modal(lefts)
    result["lineage_right_conservation"], result["modal_right"] = _lineage_modal(rights)

    # UNORDERED: the same two genes flanking the ORF count as the same neighbourhood even
    # when the plasmid is written in the opposite orientation, which is a presentation
    # artefact rather than a biological difference.
    neighbourhoods = [
        (o["lineage"], ",".join(sorted({n for n in list(o.get("left") or [])
                                        + list(o.get("right") or []) if n})))
        for o in usable
    ]
    result["lineage_neighborhood_conservation"], _ = _lineage_modal(neighbourhoods)
    result["lineage_synteny_conservation"], result["modal_synteny"] = _lineage_modal(
        syntenies)

    # The mean over lineages of the within-lineage operon-like fraction.
    operon_of = collections.defaultdict(list)
    for o in usable:
        operon_of[o["lineage"]].append(bool(o.get("operon")))
    result["lineage_operon_like_conservation"] = round(
        sum(sum(v) / len(v) for v in operon_of.values()) / len(operon_of), 4)

    result["status"] = "SUCCESS"
    return result
