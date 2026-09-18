"""Stage 14: rarity, conservation and rarefaction (spec sections 54 and 55).

Section 54 opens with the rule: rarity and conservation are SEPARATE descriptors. A family
seen on three plasmids and identical on all three, and one seen on four hundred with one
sequence, are opposite biological situations - a lineage-restricted system under strong
constraint, and a housekeeping-like gene. One axis would merge them.
"""
from plasmidann import rarity

THRESHOLDS = {
    "rare_max_lineages": 3,
    "widely_conserved_min_lineages": 50,
    "cross_min_mob": 3,
    "cross_min_hosts": 5,
    "cross_min_genera": 3,
}


def test_breadth_is_counted_in_lineages_not_plasmid_records():
    """Section 34.2: record counts are not independent observations. A family on four
    hundred redepositions of ONE plasmid is one observation, and calling it
    WIDELY_CONSERVED would be exactly the error the independence counting exists to
    prevent."""
    redeposited = {"unique_plasmid_count": 400,
                   "independent_plasmid_cluster_count": 1,
                   "MOB_count": 1, "host_count": 1, "genus_count": 1}

    labels = rarity.rarity_labels(redeposited, THRESHOLDS)

    assert "WIDELY_CONSERVED" not in labels, (
        "400 redepositions of one plasmid were called widely conserved")
    assert "LINEAGE_SPECIFIC" in labels
    assert "RARE" in labels


def test_a_genuinely_widespread_family_is_labelled_so():
    """The other direction: counting lineages must not flatten real breadth."""
    widespread = {"unique_plasmid_count": 400,
                  "independent_plasmid_cluster_count": 120,
                  "MOB_count": 6, "host_count": 40, "genus_count": 12}

    labels = rarity.rarity_labels(widespread, THRESHOLDS)

    assert "WIDELY_CONSERVED" in labels
    assert "RARE" not in labels
    assert {"CROSS_MOB", "CROSS_HOST", "CROSS_TAXON"} <= set(labels)


def test_a_family_can_carry_several_labels_at_once():
    """They describe different axes. Forcing one would make the answer depend on
    evaluation order rather than on the biology."""
    family = {"independent_plasmid_cluster_count": 2,
              "MOB_count": 4, "host_count": 6, "genus_count": 4}

    labels = rarity.rarity_labels(family, THRESHOLDS)

    assert "RARE" in labels
    assert "CROSS_MOB" in labels, "rare and cross-MOB are not mutually exclusive"


def test_a_family_confined_to_one_mob_type_across_lineages_is_plasmid_family_specific():
    """It travels with a plasmid type rather than with a host or an environment."""
    family = {"independent_plasmid_cluster_count": 20, "MOB_count": 1,
              "host_count": 8, "genus_count": 4}

    assert "PLASMID_FAMILY_SPECIFIC" in rarity.rarity_labels(family, THRESHOLDS)


def test_every_emitted_label_is_declared():
    """An undeclared label is a state nothing downstream knows how to read."""
    for family in ({"independent_plasmid_cluster_count": 1, "MOB_count": 1},
                   {"independent_plasmid_cluster_count": 200, "MOB_count": 9,
                    "host_count": 40, "genus_count": 9}):
        assert set(rarity.rarity_labels(family, THRESHOLDS)) <= set(rarity.LABELS)


def test_the_rarefaction_curve_rises_and_ends_at_the_observed_total():
    """Section 55: dark families discovered against plasmids sampled. The last point must
    be the observed total, not an extrapolation."""
    plasmid_families = {f"p{i}": {f"F{i}"} for i in range(20)}

    curve = rarity.rarefaction(plasmid_families, n_replicates=5, seed=1)

    assert curve[0]["mean_families"] < curve[-1]["mean_families"]
    assert curve[-1]["n_plasmids"] == 20
    assert curve[-1]["mean_families"] == 20


def test_a_saturated_collection_shows_a_flat_curve():
    """Every plasmid carries the same family, so sampling more adds nothing. This is the
    shape that says more plasmids of this kind will not reveal new dark families."""
    plasmid_families = {f"p{i}": {"F1"} for i in range(50)}

    curve = rarity.rarefaction(plasmid_families, n_replicates=3, seed=1)

    assert {point["mean_families"] for point in curve} == {1.0}
    # A curve that never climbed has no initial slope to compare against, so the ratio is
    # undefined rather than 0. Reporting 0 would claim "discovery has stopped", which
    # implies it started - and here nothing was ever discovered beyond the first plasmid.
    assert rarity.saturation(curve) == ""


def test_the_curve_is_averaged_over_replicates_not_one_ordering():
    """A single ordering is one arbitrary curve - starting with the most gene-rich plasmid
    makes discovery look fast. The shape is the entire output, so it must not be an
    artefact of one shuffle."""
    plasmid_families = {f"p{i}": {f"F{i % 5}"} for i in range(30)}

    curve = rarity.rarefaction(plasmid_families, n_replicates=10, seed=3)

    assert all(point["n_replicates"] == 10 for point in curve)
    assert any(point["min_families"] != point["max_families"] for point in curve), (
        "no replicate spread reported, so curve instability would be invisible")


def test_rarefaction_of_nothing_is_empty_not_an_error():
    assert rarity.rarefaction({}) == []
    assert rarity.saturation([]) == ""


def test_saturation_is_not_fooled_by_an_uneven_final_step():
    """The first implementation compared the raw gain between the last two points. On a
    curve whose final step was 2 plasmids wide where the others were 9, that reported a
    steeply climbing collection as saturated - the exact wrong answer, since it would say a
    dark set is complete when it is a lower bound.

    Slope is families per plasmid added, so step width cannot change the reading."""
    steep = [
        {"n_plasmids": 10, "mean_families": 100},
        {"n_plasmids": 20, "mean_families": 200},
        {"n_plasmids": 90, "mean_families": 900},
        # A deliberately narrow final step, as the real curve had.
        {"n_plasmids": 92, "mean_families": 920},
    ]

    assert rarity.saturation(steep) == 1.0, (
        "a curve climbing at a constant rate was not reported as still climbing")


def test_a_flattening_curve_gives_a_low_saturation_value():
    flattening = [
        {"n_plasmids": 10, "mean_families": 100},
        {"n_plasmids": 20, "mean_families": 200},
        {"n_plasmids": 90, "mean_families": 300},
        {"n_plasmids": 100, "mean_families": 301},
    ]

    value = rarity.saturation(flattening)

    assert value != "" and value < 0.05


def test_the_sample_sizes_end_on_an_even_step():
    """An uneven last interval is what produced the false reading, so the sizes themselves
    are built to avoid it."""
    plasmid_families = {f"p{i}": {f"F{i}"} for i in range(92)}

    curve = rarity.rarefaction(plasmid_families, n_replicates=2, seed=1)
    sizes = [point["n_plasmids"] for point in curve]

    assert sizes[-1] == 92, "the curve must end on the observed total"
    gaps = {sizes[i + 1] - sizes[i] for i in range(len(sizes) - 2)}
    assert len(gaps) == 1, f"uneven sampling steps before the final point: {sorted(gaps)}"
