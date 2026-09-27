"""Dark family co-occurrence (plasmidann.cooccurrence).

Two dark families travel together in a lineage when some plasmid of that lineage carries a
member ORF of each. The test is the hypergeometric upper tail over lineages, with
Benjamini-Hochberg across every pair that shares a plasmid.
"""
import itertools
import math
from fractions import Fraction

import pytest

from plasmidann.cooccurrence import (
    benjamini_hochberg,
    cooccurrence,
    family_partners,
    hypergeom_sf,
)


def _exact_tail(k, N, K, n):
    """P(X >= k) from exact integer counts."""
    return Fraction(sum(math.comb(K, i) * math.comb(N - K, n - i)
                        for i in range(k, min(K, n) + 1)), math.comb(N, n))


def test_the_tail_matches_an_enumeration_of_every_draw():
    """Literal enumeration: every n-subset of N lineages, K of them marked."""
    N = 9
    for K in range(N + 1):
        marked = set(range(K))
        for n in range(N + 1):
            draws = list(itertools.combinations(range(N), n))
            for k in range(n + 2):
                hits = sum(1 for d in draws if len(marked.intersection(d)) >= k)
                assert hypergeom_sf(k, N, K, n) == pytest.approx(hits / len(draws),
                                                                 rel=1e-12, abs=1e-15)


@pytest.mark.parametrize("k,N,K,n", [
    (50, 1000, 50, 50),        # always together: p = 1 / C(1000, 50), about 1e-85
    (10, 5000, 10, 10),
    (30, 143503, 200, 300),    # a tail deep in the far end at full scale
    (3, 143503, 200, 300),     # near the mean (0.42): p is not small
    (0, 100, 10, 10),
    (1, 100, 10, 10),
    (7, 60, 30, 20),           # below the mode (10): the lower-tail complement
    (40, 20000, 500, 2000),    # the lower tail again (mode 50), at a larger N
])
def test_the_tail_matches_exact_counts_in_the_far_tail_and_near_the_mean(k, N, K, n):
    exact = _exact_tail(k, N, K, n)
    got = hypergeom_sf(k, N, K, n)
    if exact == 0 or float(exact) == 0.0:
        assert got == 0.0
    else:
        assert got == pytest.approx(float(exact), rel=1e-9)


def test_the_tail_is_one_at_or_below_the_support_and_zero_above_it():
    assert hypergeom_sf(0, 10, 3, 4) == 1.0
    # n + K - N = 2: at least two overlaps are certain.
    assert hypergeom_sf(2, 10, 6, 6) == 1.0
    assert hypergeom_sf(5, 10, 3, 4) == 0.0


def test_benjamini_hochberg_matches_a_hand_computation():
    # Sorted p: 0.01, 0.02, 0.03, 0.04, 0.2; m = 5.
    # raw p*m/rank: 0.05, 0.05, 0.05, 0.05, 0.2 - then the running minimum from the top.
    p = [0.04, 0.01, 0.2, 0.03, 0.02]
    assert benjamini_hochberg(p) == pytest.approx([0.05, 0.05, 0.2, 0.05, 0.05])
    # Monotone step-up: a larger p never gets a smaller q, and q is capped at 1.
    assert benjamini_hochberg([0.5, 0.9, 0.95]) == pytest.approx([0.95, 0.95, 0.95])
    assert benjamini_hochberg([]) == []


def _world(n_lineages=40):
    """Plasmids P<i> in lineage L<i>, one plasmid per lineage unless stated."""
    lineage_of = {f"P{i}": f"L{i}" for i in range(n_lineages)}
    plasmid_families = {f"P{i}": set() for i in range(n_lineages)}
    return lineage_of, plasmid_families


def test_always_together_is_significant_and_independent_is_not():
    lineage_of, pf = _world()
    for i in range(10):                 # A and B share a plasmid in ten lineages
        pf[f"P{i}"] |= {"A", "B"}
    for i in range(20):                 # C in lineages 0-19, D on every even lineage:
        pf[f"P{i}"].add("C")            # overlap 10, exactly the expectation 20*20/40
    for i in range(0, 40, 2):
        pf[f"P{i}"].add("D")
    rows = {(r["family_a"], r["family_b"]): r
            for r in cooccurrence(pf, lineage_of, n_lineages=40, min_lineages_together=2)}

    ab = rows[("A", "B")]
    assert (ab["n_lineages_a"], ab["n_lineages_b"], ab["n_lineages_together"],
            ab["n_lineages_total"]) == (10, 10, 10, 40)
    assert (ab["fraction_of_a"], ab["fraction_of_b"]) == (1.0, 1.0)
    assert ab["expected_together"] == pytest.approx(2.5)
    assert ab["p_value"] == pytest.approx(1 / math.comb(40, 10), rel=1e-9)
    assert ab["q_value"] < 0.05

    cd = rows[("C", "D")]
    assert cd["n_lineages_together"] == 10 and cd["expected_together"] == pytest.approx(10)
    assert cd["p_value"] > 0.4 and cd["q_value"] > 0.05
    # Every written row is a tested pair, so there is no status column.
    assert not any("status" in r for r in rows.values())


def test_clonal_copies_in_one_lineage_count_once():
    """Twenty redeposited copies of one plasmid are one observation, not twenty."""
    lineage_of, pf = _world()
    for j in range(20):
        lineage_of[f"clone{j}"] = "L0"
        pf[f"clone{j}"] = {"A", "B"}
    pf["P1"] |= {"A", "B"}
    rows = cooccurrence(pf, lineage_of, n_lineages=40, min_lineages_together=2)
    (ab,) = rows
    assert (ab["n_lineages_a"], ab["n_lineages_b"], ab["n_lineages_together"]) == (2, 2, 2)


def test_together_means_the_same_plasmid_not_the_same_lineage():
    """A on one plasmid of a lineage and B on another plasmid of it do not travel together."""
    lineage_of, pf = _world()
    for i in range(5):
        lineage_of[f"Q{i}"] = f"L{i}"
        pf[f"P{i}"].add("A")
        pf[f"Q{i}"] = {"B"}
    assert cooccurrence(pf, lineage_of, n_lineages=40, min_lineages_together=2) == []
    pf["P0"].add("B")
    pf["P1"].add("B")
    (ab,) = cooccurrence(pf, lineage_of, n_lineages=40, min_lineages_together=2)
    # Both families are in all five lineages; they share a plasmid in two.
    assert (ab["n_lineages_a"], ab["n_lineages_b"], ab["n_lineages_together"]) == (5, 5, 2)
    assert (ab["fraction_of_a"], ab["fraction_of_b"]) == (0.4, 0.4)


def test_a_pair_together_in_one_lineage_is_not_reported():
    """Two singletons on one plasmid would get p = 1/N from one observation."""
    lineage_of, pf = _world()
    pf["P0"] |= {"A", "B"}
    assert cooccurrence(pf, lineage_of, n_lineages=40, min_lineages_together=2) == []


def test_pairs_together_in_one_lineage_still_count_in_the_correction():
    """The k >= 2 cut filters on the test statistic, so it is applied after
    Benjamini-Hochberg, not before (Bourgon et al. 2010, PNAS 107:9546). C and D, each in
    two lineages, share a plasmid in one: not reported, but one of the m = 2 tests."""
    lineage_of, pf = _world()
    pf["P0"] |= {"A", "B"}
    pf["P1"] |= {"A", "B"}
    pf["P2"] |= {"C", "D"}
    pf["P3"].add("C")
    pf["P4"].add("D")
    (ab,) = cooccurrence(pf, lineage_of, n_lineages=40, min_lineages_together=2)
    p_ab = 1 / math.comb(40, 2)
    p_cd = 1 - math.comb(38, 2) / math.comb(40, 2)
    assert ab["p_value"] == pytest.approx(p_ab, rel=1e-9)
    assert ab["q_value"] == pytest.approx(benjamini_hochberg([p_ab, p_cd])[0], rel=1e-9)
    assert ab["q_value"] == pytest.approx(2 * p_ab, rel=1e-9)


def test_a_plasmid_without_a_lineage_is_an_error():
    lineage_of, pf = _world()
    pf["orphan"] = {"A", "B"}
    with pytest.raises(ValueError, match="orphan"):
        cooccurrence(pf, lineage_of, n_lineages=40, min_lineages_together=2)


def test_family_partners_counts_significant_partners_and_names_the_best():
    rows = [
        {"family_a": "A", "family_b": "B", "n_lineages_together": 4,
         "fraction_of_a": 0.8, "fraction_of_b": 1.0, "p_value": 1e-6, "q_value": 1e-5},
        {"family_a": "A", "family_b": "C", "n_lineages_together": 2,
         "fraction_of_a": 0.4, "fraction_of_b": 0.1, "p_value": 0.01, "q_value": 0.03},
        {"family_a": "C", "family_b": "D", "n_lineages_together": 2,
         "fraction_of_a": 0.1, "fraction_of_b": 0.5, "p_value": 0.3, "q_value": 0.6},
    ]
    got = family_partners(rows, fdr=0.05)
    assert got["A"] == {"n_cooccurring_partners": 2, "top_cooccurring_partner": "B",
                        "top_cooccurring_partner_q": 1e-5,
                        "top_cooccurring_partner_fraction": 0.8}
    # The fraction is of the focal family's lineages: B shares a plasmid with A in all.
    assert got["B"]["top_cooccurring_partner_fraction"] == 1.0
    assert got["C"]["n_cooccurring_partners"] == 1
    assert got["C"]["top_cooccurring_partner"] == "A"
    # Tested but nothing significant: the best partner is still named, with its q.
    assert got["D"] == {"n_cooccurring_partners": 0, "top_cooccurring_partner": "C",
                        "top_cooccurring_partner_q": 0.6,
                        "top_cooccurring_partner_fraction": 0.5}
