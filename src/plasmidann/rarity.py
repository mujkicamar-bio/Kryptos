"""Stage 14: rarity, conservation and rarefaction (spec sections 54 and 55).

RARITY AND CONSERVATION ARE SEPARATE DESCRIPTORS

Section 54 opens with it, and the distinction is easy to lose. A family can be:

    rare and conserved        seen on three plasmids, identical on all three
    common and variable       seen on four hundred, diverging freely
    rare and variable         three plasmids, three different versions
    common and conserved      four hundred plasmids, one sequence

Collapsing these into one axis would merge the first and the last, which are opposite
biological situations: one is a lineage-restricted system under strong constraint, the other
is a housekeeping-like gene. So rarity is measured from BREADTH and conservation from
SEQUENCE IDENTITY, and the two are reported side by side.

THE LABELS ARE NOT A RANKING

Section 54: "These labels are not experimental rankings." RARE is not better than
WIDELY_CONSERVED, and neither is a score. They are named states that a downstream
prioritisation workflow may select on (section 77), and a reader can see exactly which
threshold produced each one because the thresholds are configuration and travel with the
row.

BREADTH IS COUNTED IN LINEAGES, NOT RECORDS

Every threshold here reads independent_plasmid_cluster_count from Stage 7, never the raw
plasmid count. A family on four hundred redepositions of one plasmid is not WIDELY_CONSERVED
- it is one observation - and section 34.2 forbids treating record counts as independent
observations. Using the record count would relabel exactly the families where the
distinction matters most.

RAREFACTION (section 55)

"This estimates whether additional plasmids continue to reveal new dark families." The curve
is dark families discovered against plasmids sampled. If it is still climbing at the full
sample size, the collection has not saturated and the dark set is a lower bound; if it has
flattened, more plasmids of this kind will not add much.

Section 55 is explicit that the ~1,000-candidate figure is "an experimental-budget
objective, not a biological assumption", so nothing here treats the curve as a target.
"""
import random

# Bumped when a label's definition changes, so two runs' labels cannot be silently compared.
RARITY_VERSION = "1"

LABELS = (
    "RARE",
    "LINEAGE_SPECIFIC",
    "PLASMID_FAMILY_SPECIFIC",
    "WIDELY_CONSERVED",
    "CROSS_MOB",
    "CROSS_HOST",
    "CROSS_TAXON",
)


def rarity_labels(family, thresholds):
    """Every label that applies to one family. A family may carry several.

    `family` needs the Stage 7 distribution counts: independent_plasmid_cluster_count,
    MOB_count, host_count, genus_count, unique_plasmid_count.

    Several labels rather than one, because they describe different axes and a family can
    be genuinely CROSS_MOB and CROSS_HOST at once. Forcing one label would make the answer
    depend on evaluation order rather than on the biology.
    """
    def count(name):
        try:
            return int(family.get(name) or 0)
        except (TypeError, ValueError):
            return 0

    lineages = count("independent_plasmid_cluster_count")
    labels = []

    # Breadth, always in LINEAGES. See the module docstring: the raw plasmid count would
    # relabel exactly the families where redeposition is doing the work.
    if lineages and lineages <= thresholds["rare_max_lineages"]:
        labels.append("RARE")
    if lineages == 1:
        labels.append("LINEAGE_SPECIFIC")
    if lineages >= thresholds["widely_conserved_min_lineages"]:
        labels.append("WIDELY_CONSERVED")

    # Confined to one MOB cluster while occurring on several lineages: the family travels
    # with a plasmid type rather than with a host or an environment.
    if count("MOB_count") == 1 and lineages > 1:
        labels.append("PLASMID_FAMILY_SPECIFIC")

    if count("MOB_count") >= thresholds["cross_min_mob"]:
        labels.append("CROSS_MOB")
    if count("host_count") >= thresholds["cross_min_hosts"]:
        labels.append("CROSS_HOST")
    if count("genus_count") >= thresholds["cross_min_genera"]:
        labels.append("CROSS_TAXON")

    return labels


def rarefaction(plasmid_families, sample_sizes=None, n_replicates=20, seed=0):
    """Dark families discovered against plasmids sampled (section 55).

    `plasmid_families` maps plasmid_id to the set of dark family ids on it. Returns a list
    of dicts: n_plasmids, mean_families, min_families, max_families, n_replicates.

    Averaged over `n_replicates` random orderings rather than computed once. A single
    ordering is one arbitrary curve - starting with the most gene-rich plasmid makes
    discovery look fast, starting with cryptic ones makes it look slow - and the shape of
    the curve is the entire output, so it must not be an artefact of one shuffle.

    The replicate spread is reported, not just the mean, because a wide spread means the
    curve is unstable and its flattening cannot be trusted.
    """
    plasmids = sorted(plasmid_families)
    if not plasmids:
        return []

    if sample_sizes is None:
        # Ten roughly even steps, always including the full set, so the last point is the
        # observed total rather than an extrapolation.
        step = max(1, len(plasmids) // 10)
        sample_sizes = sorted(set(list(range(step, len(plasmids), step))
                                  + [len(plasmids)]))

    rng = random.Random(seed)
    curve = []
    for size in sample_sizes:
        counts = []
        for _ in range(n_replicates):
            sample = rng.sample(plasmids, size)
            seen = set()
            for plasmid in sample:
                seen |= plasmid_families[plasmid]
            counts.append(len(seen))
        curve.append({
            "n_plasmids": size,
            "mean_families": round(sum(counts) / len(counts), 2),
            "min_families": min(counts),
            "max_families": max(counts),
            "n_replicates": n_replicates,
        })
    return curve


def saturation(curve):
    """How much the last decile of sampling added, as a fraction of the total.

    A single descriptive number for "is this still climbing". Near 0 means the collection
    has saturated and more plasmids of this kind will add little; a large value means the
    dark set is a lower bound.

    It is NOT a target and nothing selects on it. Section 55: the ~1,000-candidate figure
    is an experimental-budget objective, not a biological assumption.
    """
    if len(curve) < 2:
        return ""
    last, previous = curve[-1], curve[-2]
    if not last["mean_families"]:
        return ""
    gained = last["mean_families"] - previous["mean_families"]
    return round(gained / last["mean_families"], 4)
