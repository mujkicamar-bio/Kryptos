"""Rarity labels, the rarefaction curve and its saturation value (plasmidann.rarity)."""
from plasmidann import rarity

THRESHOLDS = {
    "rare_max_lineages": 3,
    "widespread_min_lineages": 50,
    "cross_min_hosts": 2,
    "cross_min_genera": 2,
}


def test_breadth_is_counted_in_lineages_not_plasmid_records():
    """A family on four hundred redeposits of one plasmid is one observation."""
    redeposited = {"unique_plasmid_count": 400,
                   "independent_plasmid_cluster_count": 1,
                   "MOB_count": 1, "host_count": 1, "genus_count": 1}

    labels = rarity.rarity_labels(redeposited, THRESHOLDS)

    assert "WIDESPREAD" not in labels, (
        "400 redepositions of one plasmid were called widespread")
    assert "LINEAGE_SPECIFIC" in labels
    assert "RARE" in labels


def test_a_genuinely_widespread_family_is_labelled_so():
    """Counting lineages keeps real breadth."""
    widespread = {"unique_plasmid_count": 400,
                  "independent_plasmid_cluster_count": 120,
                  "MOB_count": 6, "host_count": 40, "genus_count": 12}

    labels = rarity.rarity_labels(widespread, THRESHOLDS)

    assert "WIDESPREAD" in labels
    assert "RARE" not in labels
    assert {"CROSS_MOB", "CROSS_HOST", "CROSS_TAXON"} <= set(labels)


def test_the_widespread_threshold_is_the_nearest_rank_percentile_of_measured_counts():
    assert rarity.widespread_threshold(range(1, 101), 99) == 99
    assert rarity.widespread_threshold([1] * 95 + [2, 3, 4, 50, 60], 95) == 1
    # A count of 0 is not measured and does not dilute the distribution.
    assert rarity.widespread_threshold([0] * 1000 + [1, 2, 3, 4], 50) == 2
    assert rarity.widespread_threshold([0, 0], 99) is None
    # 7 * 100 / 100 is exactly rank 7; 0.07 * 100 in floating point would round up to 8.
    assert rarity.widespread_threshold(range(1, 101), 7) == 7


def test_a_family_without_a_measured_lineage_count_is_not_widespread():
    labels = rarity.rarity_labels({"independent_plasmid_cluster_count": 0},
                                  {**THRESHOLDS, "widespread_min_lineages": None})
    assert "WIDESPREAD" not in labels


def test_a_family_can_carry_several_labels_at_once():
    family = {"independent_plasmid_cluster_count": 2,
              "MOB_count": 4, "host_count": 6, "genus_count": 4}

    labels = rarity.rarity_labels(family, THRESHOLDS)

    assert "RARE" in labels
    assert "CROSS_MOB" in labels, "rare and cross-MOB are not mutually exclusive"


def test_mobsuite_is_reported_as_single_or_cross_mob_only():
    """MOB-suite clusters are a reported label, never a lineage count: exactly one cluster
    is SINGLE_MOB, two or more CROSS_MOB, and a family with none carries neither."""
    def labels(mob, lineages=20):
        return set(rarity.rarity_labels(
            {"independent_plasmid_cluster_count": lineages, "MOB_count": mob}, THRESHOLDS))

    assert "SINGLE_MOB" in labels(1) and "CROSS_MOB" not in labels(1)
    assert "SINGLE_MOB" in labels(1, lineages=1)
    assert "CROSS_MOB" in labels(2) and "SINGLE_MOB" not in labels(2)
    assert not labels(0) & {"SINGLE_MOB", "CROSS_MOB"}


def test_single_host_needs_every_plasmid_to_have_a_species():
    """One species on 1 of 5 plasmids says nothing about the other 4, so it is not
    SINGLE_HOST; one species on all 5 is."""
    def labels(with_species):
        return set(rarity.rarity_labels(
            {"independent_plasmid_cluster_count": 5, "unique_plasmid_count": 5,
             "host_count": 1, "genus_count": 1,
             "n_plasmids_with_species": with_species}, THRESHOLDS))

    assert "SINGLE_HOST" not in labels(1)
    assert "SINGLE_HOST" in labels(5)
    assert "CROSS_HOST" not in labels(5)


def test_two_hosts_or_two_genera_are_cross_host_and_cross_taxon():
    family = {"independent_plasmid_cluster_count": 2, "unique_plasmid_count": 2,
              "host_count": 2, "genus_count": 2, "n_plasmids_with_species": 2}

    labels = set(rarity.rarity_labels(family, THRESHOLDS))

    assert {"CROSS_HOST", "CROSS_TAXON"} <= labels
    assert "SINGLE_HOST" not in labels


def test_the_rarefaction_curve_rises_and_ends_at_the_observed_total():
    """The last point is the observed total, not an extrapolation."""
    lineage_families = {f"p{i}": {f"F{i}"} for i in range(20)}

    curve = rarity.rarefaction(lineage_families, n_replicates=5, seed=1)

    assert curve[0]["mean_families"] < curve[-1]["mean_families"]
    assert curve[-1]["n_lineages"] == 20
    assert curve[-1]["mean_families"] == 20


def test_a_saturated_collection_shows_a_flat_curve():
    """Every lineage carries the same family, so sampling more adds nothing."""
    lineage_families = {f"p{i}": {"F1"} for i in range(50)}

    curve = rarity.rarefaction(lineage_families, n_replicates=3, seed=1)

    assert {point["mean_families"] for point in curve} == {1.0}
    # A curve that never climbed has no initial slope, so the ratio is undefined, not 0.
    assert rarity.saturation(curve) == ""


def test_the_curve_is_averaged_over_replicates_not_one_ordering():
    lineage_families = {f"p{i}": {f"F{i % 5}"} for i in range(30)}

    curve = rarity.rarefaction(lineage_families, n_replicates=10, seed=3)

    assert all(point["n_replicates"] == 10 for point in curve)
    assert any(point["min_families"] != point["max_families"] for point in curve), (
        "no replicate spread reported, so curve instability would be invisible")


def test_rarefaction_of_nothing_is_empty_not_an_error():
    assert rarity.rarefaction({}) == []
    assert rarity.saturation([]) == ""


def test_saturation_does_not_depend_on_step_width():
    """Slope is families per lineage added, so a narrow final step reads the same."""
    steep = [
        {"n_lineages": 10, "mean_families": 100},
        {"n_lineages": 20, "mean_families": 200},
        {"n_lineages": 90, "mean_families": 900},
        {"n_lineages": 92, "mean_families": 920},
    ]

    assert rarity.saturation(steep) == 1.0, (
        "a curve climbing at a constant rate was not reported as still climbing")


def test_a_flattening_curve_gives_a_low_saturation_value():
    flattening = [
        {"n_lineages": 10, "mean_families": 100},
        {"n_lineages": 20, "mean_families": 200},
        {"n_lineages": 90, "mean_families": 300},
        {"n_lineages": 100, "mean_families": 301},
    ]

    value = rarity.saturation(flattening)

    assert value != "" and value < 0.05


def test_the_sample_sizes_are_ten_equal_steps_ending_on_the_total():
    """Every interval, the last included, is equal to within one lineage, so the final
    slope is measured over as many lineages as the others."""
    for n in (92, 15, 101):
        curve = rarity.rarefaction({f"p{i}": {f"F{i}"} for i in range(n)}, n_replicates=1)
        sizes = [0] + [point["n_lineages"] for point in curve]
        gaps = [b - a for a, b in zip(sizes, sizes[1:])]

        assert sizes[-1] == n and len(gaps) == 10
        assert max(gaps) - min(gaps) <= 1, f"uneven steps for {n} lineages: {gaps}"


def test_fewer_than_ten_lineages_give_one_point_per_distinct_size():
    curve = rarity.rarefaction({"p0": {"F0"}, "p1": {"F1"}}, n_replicates=1)

    assert [point["n_lineages"] for point in curve] == [1, 2]
