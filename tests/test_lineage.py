"""Plasmid lineages from Mash distances (plasmidann.lineage)."""
import io

import pytest

from plasmidann import lineage

MASH = "".join(f"{line}\n" for line in [
    # reference query distance p-value shared-hashes
    "pA\tpA\t0\t0\t1000/1000",
    "pA\tpB\t0.01\t0\t900/1000",
    "pB\tpA\t0.01\t0\t900/1000",
    "pB\tpC\t0.02\t0\t850/1000",
    "pC\tpB\t0.02\t0\t850/1000",
    "pA\tpD\t0.30\t0\t10/1000",
    "pD\tpA\t0.30\t0\t10/1000",
])


def _edges(text, **kw):
    return set(lineage.parse_mash_dist(io.StringIO(text), **kw))


def test_pairs_below_the_threshold_become_edges_and_self_hits_do_not():
    edges = _edges(MASH, max_distance=0.05)

    assert edges == {("pA", "pB"), ("pB", "pA"), ("pB", "pC"), ("pC", "pB")}


def test_a_low_distance_from_almost_no_shared_hashes_is_rejected():
    """Mash reports a small distance from two shared hashes with a large p-value."""
    noise = "pX\tpY\t0.01\t0.9\t2/1000\n"

    assert _edges(noise, max_distance=0.05) == set()
    assert _edges(noise, max_distance=0.05, max_pvalue=1.0) == {("pX", "pY")}


def test_a_malformed_line_raises():
    with pytest.raises(ValueError):
        _edges("pA\tpB\n", max_distance=0.05)


def test_a_chain_of_similar_plasmids_is_one_lineage():
    """Single linkage: pA and pC join through pB."""
    clusters = lineage.lineage_clusters(["pA", "pB", "pC", "pD"],
                                        _edges(MASH, max_distance=0.05))

    assert clusters == {"pA": "pA", "pB": "pA", "pC": "pA", "pD": "pD"}


def test_a_plasmid_resembling_nothing_is_a_lineage_of_one():
    assert lineage.lineage_clusters(["solo"], []) == {"solo": "solo"}


def test_cluster_ids_do_not_depend_on_iteration_order():
    names = ["pC", "pA", "pB"]
    edges = [("pA", "pB"), ("pB", "pC")]

    forward = lineage.lineage_clusters(names, edges)
    reverse = lineage.lineage_clusters(list(reversed(names)), list(reversed(edges)))

    assert forward == reverse == {"pA": "pA", "pB": "pA", "pC": "pA"}
