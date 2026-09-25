"""Family network rules, taken from Durairaj et al. (Nature 2023, Methods)."""
from plasmidann.network import (communities, edges_from_hits, member_brightness,
                                node_brightness)

RULE = {"max_out": 4, "min_cov": 0.5, "max_evalue": 1e-4}


def _hit(q, t, ev, qcov, tcov):
    return {"query": q, "target": t, "evalue": ev, "qcov": qcov, "tcov": tcov}


def test_coverage_of_either_protein_is_enough():
    """'covers at least 50% of one of the proteins': a short protein wholly inside a long
    one is linked, although the alignment covers little of the long one."""
    assert ("a", "b") in edges_from_hits([_hit("a", "b", 1e-10, 0.95, 0.2)], **RULE)
    assert edges_from_hits([_hit("a", "b", 1e-10, 0.4, 0.3)], **RULE) == {}


def test_evalue_must_be_below_the_cut_and_self_hits_are_dropped():
    assert edges_from_hits([_hit("a", "b", 1e-4, 0.9, 0.9)], **RULE) == {}
    assert edges_from_hits([_hit("a", "a", 1e-50, 1, 1)], **RULE) == {}


def test_each_node_keeps_its_four_best_outbound_edges():
    hits = [_hit("a", t, ev, 0.9, 0.9) for t, ev in
            [("b", 1e-9), ("c", 1e-30), ("d", 1e-20), ("e", 1e-10), ("f", 1e-5)]]
    edges = edges_from_hits(hits, **RULE)
    assert set(edges) == {("a", "c"), ("a", "d"), ("a", "e"), ("a", "b")}


def test_an_edge_found_both_ways_is_one_edge_with_the_better_evalue():
    edges = edges_from_hits([_hit("a", "b", 1e-8, 0.9, 0.9),
                             _hit("b", "a", 1e-12, 0.9, 0.9)], **RULE)
    assert edges == {("a", "b"): 1e-12}


def test_brightness_is_what_the_best_member_reaches():
    rows = [{"functional_class": "NONE", "explained_fraction": "0.0"},
            {"functional_class": "DOMAIN_ONLY", "explained_fraction": "0.3"}]
    assert node_brightness(rows) == 0.3


def test_an_unmeasured_functional_protein_is_bright_not_dark():
    """PlasmidScope Tier-0 and family-level pharokka rows have no span. Reading their
    empty explained_fraction as 0 would make known proteins look dark."""
    assert member_brightness({"functional_class": "FUNCTIONAL", "explained_fraction": "",
                              "annot_completeness": "NOT_MEASURED"}) == 1.0


def test_communities_repeat_with_the_same_seed():
    edges = {("a", "b"): 1e-20, ("b", "c"): 1e-20, ("d", "e"): 1e-20}
    nodes = ["a", "b", "c", "d", "e", "f"]
    first = communities(edges, nodes, seed=1)
    assert first == communities(edges, nodes, seed=1)
    assert first["a"] == first["b"] == first["c"] != first["d"]
    assert first["f"] not in (first["a"], first["d"]), "an isolated node is its own group"
