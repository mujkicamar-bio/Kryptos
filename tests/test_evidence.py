"""Reality lines and darkness state (plasmidann.evidence). Lineage breadth is counted in
Stage 6 lineages, not MOB-suite clusters.

MOB-suite gives every plasmid the cluster of its nearest reference however distant, so
unrelated novel plasmids share a cluster (AA379 holds 21,817 small plasmids of the analysis
set). Counting MOB-suite clusters would under-count independence exactly on small plasmids.
"""
from plasmidann import evidence

THRESHOLDS = evidence.reality_thresholds(
    {"dnds_purifying_max": 0.5, "min_members_for_dnds": 3})


def test_multi_lineage_fires_on_two_lineages():
    _, fired, _ = evidence.reality_lines(
        {"independent_plasmid_cluster_count": "2"}, THRESHOLDS)
    assert "multi_lineage" in fired


def test_multi_lineage_ignores_mobsuite_clusters():
    _, fired, _ = evidence.reality_lines(
        {"independent_plasmid_cluster_count": "1", "n_mob_clusters": "7"}, THRESHOLDS)
    assert "multi_lineage" not in fired


def test_an_unmeasured_lineage_count_does_not_fire():
    _, fired, _ = evidence.reality_lines({"independent_plasmid_cluster_count": ""},
                                         THRESHOLDS)
    assert "multi_lineage" not in fired


def test_purifying_selection_fires_only_on_a_measured_dnds_below_the_threshold():
    def fires(record):
        return "purifying_selection" in evidence.reality_lines(record, THRESHOLDS)[1]

    assert fires({"dnds_status": "MEASURED", "dnds_median": "0.2"})
    assert not fires({"dnds_status": "MEASURED", "dnds_median": "0.5"})
    assert not fires({"dnds_status": "NO_DIVERGENCE", "dnds_median": ""})
    assert not fires({"dnds_status": "TOO_SHORT", "dnds_median": "0.1"})


def test_is_family_is_implied_when_purifying_selection_fires():
    n, fired, implied = evidence.reality_lines(
        {"dnds_status": "MEASURED", "dnds_median": "0.2", "n_members": "5",
         "independent_plasmid_cluster_count": "3", "structural_match": "1abc_A"},
        THRESHOLDS)
    assert (n, fired, implied) == (3, ["purifying_selection", "multi_lineage", "folds"],
                                   ["is_family"])
    n, fired, implied = evidence.reality_lines({"n_members": "5"}, THRESHOLDS)
    assert (n, fired, implied) == (1, ["is_family"], [])


def test_the_darkness_state_follows_the_structural_match():
    assert evidence.darkness_state({"structural_match": "1abc_A"}) == "DARK_FOLD_KNOWN"
    assert evidence.darkness_state({"structural_match": ""}) == "DARK_NO_STRUCTURE"
