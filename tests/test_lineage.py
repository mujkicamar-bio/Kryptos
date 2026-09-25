"""Stage 6: plasmid lineage clusters (spec section 33).

MOB class and sequence similarity are different concepts. Section 33.2's worked example is
the whole argument: 2,143 plasmid occurrences, 8 MOB clusters, 47 independent plasmid
clusters. A protein on 2,143 records is not 2,143 observations, and a MOB breadth of 8 does
not establish that it is, because one widespread lineage can carry one relaxase type.
"""
from plasmidann import lineage

MASH = "\n".join([
    # reference query distance p-value shared-hashes
    "pA\tpA\t0\t0\t1000/1000",
    "pA\tpB\t0.01\t0\t900/1000",
    "pB\tpA\t0.01\t0\t900/1000",
    "pB\tpC\t0.02\t0\t850/1000",
    "pC\tpB\t0.02\t0\t850/1000",
    "pA\tpD\t0.30\t0\t10/1000",
    "pD\tpA\t0.30\t0\t10/1000",
])


def test_pairs_below_the_threshold_become_edges():
    edges = set(lineage.parse_mash_dist(MASH, max_distance=0.05))

    assert ("pA", "pB") in edges
    assert ("pB", "pC") in edges
    assert ("pA", "pD") not in edges, "a distant pair was linked"


def test_self_comparisons_are_not_edges():
    """Mash reports every sequence against itself at distance 0. Admitted as an edge it
    would be harmless here, but it hides a real bug: a run where nothing matched would
    still produce one edge per plasmid and look like it had worked."""
    assert ("pA", "pA") not in set(lineage.parse_mash_dist(MASH, max_distance=0.05))


def test_a_low_distance_from_almost_no_shared_hashes_is_rejected():
    """Mash can report a small distance computed from two shared hashes out of a thousand,
    and says so with a large p-value. Filtering on distance alone would link plasmids that
    share almost nothing."""
    noise = "pX\tpY\t0.01\t0.9\t2/1000"

    assert list(lineage.parse_mash_dist(noise, max_distance=0.05)) == []
    assert list(lineage.parse_mash_dist(noise, max_distance=0.05, max_pvalue=1.0)) == \
        [("pX", "pY")]


def test_a_chain_of_similar_plasmids_is_one_lineage():
    """Single linkage is deliberate: A and C join through B because the question is whether
    they could descend from one recent ancestor, and a chain of near-identical
    intermediates is evidence that they could."""
    edges = lineage.parse_mash_dist(MASH, max_distance=0.05)

    clusters = lineage.lineage_clusters(["pA", "pB", "pC", "pD"], edges)

    assert clusters["pA"] == clusters["pB"] == clusters["pC"]
    assert clusters["pD"] != clusters["pA"], "a distant plasmid joined the lineage"


def test_a_plasmid_resembling_nothing_is_a_lineage_of_one():
    """Not absent from the table. A novel plasmid matching nothing is the common case for
    the biology this project is looking for, and dropping it would remove it from every
    independence count that follows."""
    clusters = lineage.lineage_clusters(["solo"], [])

    assert clusters == {"solo": "solo"}


def test_cluster_ids_do_not_depend_on_iteration_order():
    """Numbering that depended on order would change between runs on identical input, and
    every count derived from it would be uncomparable - the defect the family ids had."""
    names = ["pC", "pA", "pB"]
    edges = [("pA", "pB"), ("pB", "pC")]

    forward = lineage.lineage_clusters(names, edges)
    reverse = lineage.lineage_clusters(list(reversed(names)), list(reversed(edges)))

    assert forward == reverse
    assert set(forward.values()) == {"pA"}, "the id is not the smallest member"


def test_independence_is_not_the_same_as_occurrence_count():
    """Section 33.2, as a test. Five plasmid records that are all one lineage give an
    independence count of 1, which is the number any recurrence claim must use."""
    names = [f"p{i}" for i in range(5)]
    edges = [(f"p{i}", f"p{i + 1}") for i in range(4)]

    clusters = lineage.lineage_clusters(names, edges)

    assert len(names) == 5
    assert len(set(clusters.values())) == 1
