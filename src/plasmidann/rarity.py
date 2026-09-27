"""Stage 14: rarity labels per family, and the dark-family rarefaction curve.

Every label describes breadth; sequence conservation is not measured here. RARE,
LINEAGE_SPECIFIC and WIDESPREAD count independent Stage 6 lineages
(independent_plasmid_cluster_count), never plasmid records, so a family on four hundred
redeposits of one plasmid is one observation. The MOB labels count MOB-suite clusters and
the host labels observed host species and genera. The labels are descriptive and nothing
selects on them. The thresholds are configuration (config/targets.yaml, `rarity`); the two
lineage thresholds are recorded on every output row.

The rarefaction curve counts dark families discovered against Stage 6 lineages sampled,
averaged over random orderings: a lineage is one observation, so redeposited copies of one
plasmid do not flatten the curve. saturation() compares the final slope with the initial
slope, in families gained per lineage added.
"""
import random


def rarity_labels(family, thresholds):
    """Every label that applies to one family, from RARE, LINEAGE_SPECIFIC,
    WIDESPREAD, SINGLE_MOB, CROSS_MOB, SINGLE_HOST, CROSS_HOST and CROSS_TAXON.

    `family` holds the Stage 7 distribution counts: independent_plasmid_cluster_count,
    MOB_count, host_count, genus_count, n_plasmids_with_species, unique_plasmid_count.
    The labels describe different axes, so a family may carry several.
    """
    def count(name):
        try:
            return int(family.get(name) or 0)
        except (TypeError, ValueError):
            return 0

    lineages = count("independent_plasmid_cluster_count")
    labels = []

    if lineages and lineages <= thresholds["rare_max_lineages"]:
        labels.append("RARE")
    if lineages == 1:
        labels.append("LINEAGE_SPECIFIC")
    if lineages >= thresholds["widespread_min_lineages"]:
        labels.append("WIDESPREAD")

    # MOB-suite clusters are reported as a label only, never as a lineage count: MOB-suite
    # assigns the nearest reference's cluster however distant, so it is no measure of
    # independence. A family with no MOB-suite cluster carries neither label.
    if count("MOB_count") == 1:
        labels.append("SINGLE_MOB")
    if count("MOB_count") >= 2:
        labels.append("CROSS_MOB")
    # Observed hosts only (plasmidann.hosts); MOB-suite's predicted range is never a host.
    # One species is SINGLE_HOST only when EVERY plasmid is named to the species: one E.
    # coli plasmid beside four unhosted ones says nothing about the other four.
    if (count("host_count") == 1
            and count("n_plasmids_with_species") == count("unique_plasmid_count")):
        labels.append("SINGLE_HOST")
    if count("host_count") >= thresholds["cross_min_hosts"]:
        labels.append("CROSS_HOST")
    if count("genus_count") >= thresholds["cross_min_genera"]:
        labels.append("CROSS_TAXON")

    return labels


def rarefaction(lineage_families, sample_sizes=None, n_replicates=20, seed=0):
    """Dark families discovered against lineages sampled.

    `lineage_families` maps lineage id to the set of dark family ids on it. Returns a list
    of dicts: n_lineages, mean_families, min_families, max_families, n_replicates. Each
    point is the mean over `n_replicates` random samples, because one ordering gives one
    arbitrary curve; the replicate range shows how stable the curve is.
    """
    lineages = sorted(lineage_families)
    if not lineages:
        return []

    if sample_sizes is None:
        # Ten steps ending on the full set, equal to within one lineage, so the final slope
        # is measured over as many lineages as every other.
        n = len(lineages)
        sample_sizes = sorted({n * i // 10 for i in range(1, 11)} - {0})

    rng = random.Random(seed)
    curve = []
    for size in sample_sizes:
        counts = []
        for _ in range(n_replicates):
            sample = rng.sample(lineages, size)
            seen = set()
            for lineage in sample:
                seen |= lineage_families[lineage]
            counts.append(len(seen))
        curve.append({
            "n_lineages": size,
            "mean_families": round(sum(counts) / len(counts), 2),
            "min_families": min(counts),
            "max_families": max(counts),
            "n_replicates": n_replicates,
        })
    return curve


def saturation(curve):
    """Final slope / initial slope, each in families gained per lineage added.

    Near 0 the curve has flattened; near 1 discovery is as fast at the end as at the start,
    so the dark family count is a lower bound. "" when undefined: fewer than three points,
    or no gain over the first step.
    """
    if len(curve) < 3:
        return ""

    def slope(a, b):
        span = b["n_lineages"] - a["n_lineages"]
        return (b["mean_families"] - a["mean_families"]) / span if span else 0.0

    initial = slope(curve[0], curve[1])
    final = slope(curve[-2], curve[-1])
    if initial <= 0:
        return ""
    return round(final / initial, 4)
