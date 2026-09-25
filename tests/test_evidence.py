"""Reality tests (plasmidann.evidence): lineage breadth is counted in Stage 6 lineages.

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
