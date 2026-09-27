"""Stage 9: whether the arrangement around a dark ORF recurs across its family's lineages.

Six measurements, kept apart because they vary independently:

    lineage_left_conservation          share of the most common immediate left neighbour
    lineage_right_conservation         the same on the right
    lineage_neighborhood_conservation  share of the most common unordered neighbour set
    lineage_synteny_conservation       share of the most common ordered left|right pair
    lineage_operon_like_conservation   mean over lineages of the operon-like fraction
    context_recurrence                 occurrences with any neighbour at all

Each share is a fractional vote over Stage 6 plasmid lineages (see _lineage_modal), so that
many copies of one redeposited plasmid are one observation; fewer than min_lineages voting
lineages gives no value.
"""
import collections
from fractions import Fraction

# Fewer lineages than this gives no value. Two is the arithmetic minimum for a comparison,
# not a tuned threshold: one lineage is conserved with itself by construction. The pipeline
# passes `synteny.min_lineages` from config.
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


def _lineage_modal(pairs, min_lineages):
    """Fractional vote over lineages: (conservation, modal value), or ("", "") when fewer
    than `min_lineages` lineages vote.

    `pairs` is (lineage, value) per occurrence. Empty values are excluded: an ORF at the end
    of a linear record has no left neighbour, which says nothing about conservation, and a
    lineage with only empty values does not vote. Each voting lineage has weight 1, split
    equally over its k values; a lineage whose copies disagree gives 1/k to each. The
    conservation is the largest summed score divided by the number of voting lineages, and
    the modal value is the value with that score, ties broken by the smallest value. Scores
    are exact fractions, so ten weights of 1/10 are exactly one vote.
    """
    values_of = collections.defaultdict(list)
    for lineage, value in pairs:
        if value:
            values_of[lineage].append(value)
    if len(values_of) < min_lineages:
        return "", ""
    score = collections.defaultdict(Fraction)
    for values in values_of.values():
        for value in values:
            score[value] += Fraction(1, len(values))
    modal = min(score, key=lambda v: (-score[v], v))
    return round(float(score[modal] / len(values_of)), 4), modal


def conservation(occurrences, min_lineages=MIN_LINEAGES):
    """Every synteny measurement for one family, counted over lineages.

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
    # ORDERED: left-then-right as one string, strictly stronger than the unordered
    # neighbourhood below.
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
        # Lineages whose copies hold more than one ordered arrangement.
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
        # One lineage is conserved with itself, however many copies it has.
        result["status"] = "TOO_FEW_LINEAGES"
        return result

    result["lineage_left_conservation"], result["modal_left"] = _lineage_modal(
        lefts, min_lineages)
    result["lineage_right_conservation"], result["modal_right"] = _lineage_modal(
        rights, min_lineages)

    # UNORDERED: the same two genes flanking the ORF count as the same neighbourhood even
    # when the plasmid is written in the opposite orientation, which is a presentation
    # artefact rather than a biological difference.
    neighbourhoods = [
        (o["lineage"], ",".join(sorted({n for n in list(o.get("left") or [])
                                        + list(o.get("right") or []) if n})))
        for o in usable
    ]
    result["lineage_neighborhood_conservation"], _ = _lineage_modal(
        neighbourhoods, min_lineages)
    result["lineage_synteny_conservation"], result["modal_synteny"] = _lineage_modal(
        syntenies, min_lineages)

    # The mean over lineages of the within-lineage operon-like fraction.
    operon_of = collections.defaultdict(list)
    for o in usable:
        operon_of[o["lineage"]].append(bool(o.get("operon")))
    result["lineage_operon_like_conservation"] = round(
        sum(sum(v) / len(v) for v in operon_of.values()) / len(operon_of), 4)

    result["status"] = "SUCCESS"
    return result
