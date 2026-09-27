"""Dark sequence co-occurrence (plasmidann.cooccurrence).

Two unique dark protein sequences travel together in a lineage when some plasmid of that
lineage carries an ORF of each. The test is the hypergeometric upper tail over lineages,
with Benjamini-Hochberg across every enumerated pair that shares a plasmid.
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
    plasmid_seqs = {f"P{i}": set() for i in range(n_lineages)}
    return lineage_of, plasmid_seqs


def _pairs(pf, lineage_of, n_lineages=40):
    rows, _ = cooccurrence(pf, lineage_of, n_lineages=n_lineages, min_lineages_together=2)
    return rows


def test_always_together_is_significant_and_independent_is_not():
    lineage_of, pf = _world()
    for i in range(10):                 # A and B share a plasmid in ten lineages
        pf[f"P{i}"] |= {"A", "B"}
    for i in range(20):                 # C in lineages 0-19, D on every even lineage:
        pf[f"P{i}"].add("C")            # overlap 10, exactly the expectation 20*20/40
    for i in range(0, 40, 2):
        pf[f"P{i}"].add("D")
    rows = {(r["seq_a"], r["seq_b"]): r for r in _pairs(pf, lineage_of)}

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
    (ab,) = _pairs(pf, lineage_of)
    assert (ab["n_lineages_a"], ab["n_lineages_b"], ab["n_lineages_together"]) == (2, 2, 2)


def test_together_means_the_same_plasmid_not_the_same_lineage():
    """A on one plasmid of a lineage and B on another plasmid of it do not travel together."""
    lineage_of, pf = _world()
    for i in range(5):
        lineage_of[f"Q{i}"] = f"L{i}"
        pf[f"P{i}"].add("A")
        pf[f"Q{i}"] = {"B"}
    assert _pairs(pf, lineage_of) == []
    pf["P0"].add("B")
    pf["P1"].add("B")
    (ab,) = _pairs(pf, lineage_of)
    # Both sequences are in all five lineages; they share a plasmid in two.
    assert (ab["n_lineages_a"], ab["n_lineages_b"], ab["n_lineages_together"]) == (5, 5, 2)
    assert (ab["fraction_of_a"], ab["fraction_of_b"]) == (0.4, 0.4)


def test_a_pair_together_in_one_lineage_is_not_reported():
    """Two singletons on one plasmid would get p = 1/N from one observation."""
    lineage_of, pf = _world()
    pf["P0"] |= {"A", "B"}
    assert _pairs(pf, lineage_of) == []


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
    (ab,) = _pairs(pf, lineage_of)
    p_ab = 1 / math.comb(40, 2)
    p_cd = 1 - math.comb(38, 2) / math.comb(40, 2)
    assert ab["p_value"] == pytest.approx(p_ab, rel=1e-9)
    assert ab["q_value"] == pytest.approx(benjamini_hochberg([p_ab, p_cd])[0], rel=1e-9)
    assert ab["q_value"] == pytest.approx(2 * p_ab, rel=1e-9)


def test_a_plasmid_without_a_lineage_is_an_error():
    lineage_of, pf = _world()
    pf["orphan"] = {"A", "B"}
    with pytest.raises(ValueError, match="orphan"):
        _pairs(pf, lineage_of)


def test_a_sequence_in_fewer_lineages_than_the_minimum_is_not_enumerated():
    """C is in one lineage, so no pair with it can reach k >= 2: it is left out of the
    enumeration and the correction, and the pair count stays bounded."""
    lineage_of, pf = _world()
    pf["P0"] |= {"A", "B", "C"}
    pf["P1"] |= {"A", "B"}
    rows, stats = cooccurrence(pf, lineage_of, n_lineages=40, min_lineages_together=2)
    assert [(r["seq_a"], r["seq_b"]) for r in rows] == [("A", "B")]
    assert stats == {"n_sequences": 3, "n_enumerated": 2, "n_pair_occurrences": 2,
                     "n_tested": 1}


def test_family_partners_counts_significant_partner_sequences_and_names_the_best():
    """a1 and a2 are members of family FA, c of FC; b and d belong to no dark family."""
    family_of_seq = {"a1": "FA", "a2": "FA", "c": "FC"}
    rows = [
        {"seq_a": "a1", "seq_b": "b", "n_lineages_together": 4,
         "fraction_of_a": 0.8, "fraction_of_b": 1.0, "p_value": 1e-6, "q_value": 1e-5},
        {"seq_a": "a2", "seq_b": "b", "n_lineages_together": 2,
         "fraction_of_a": 0.5, "fraction_of_b": 0.5, "p_value": 1e-3, "q_value": 1e-2},
        {"seq_a": "a2", "seq_b": "c", "n_lineages_together": 2,
         "fraction_of_a": 0.4, "fraction_of_b": 0.1, "p_value": 0.01, "q_value": 0.03},
        {"seq_a": "c", "seq_b": "d", "n_lineages_together": 2,
         "fraction_of_a": 0.1, "fraction_of_b": 0.5, "p_value": 0.3, "q_value": 0.6},
    ]
    got = family_partners(rows, 0.05, family_of_seq)
    # b is a partner of both members of FA but one partner sequence; c is the other.
    assert got["FA"] == {"n_cooccurring_partners": 2, "top_cooccurring_partner": "b",
                         "top_cooccurring_partner_q": 1e-5,
                         "top_cooccurring_partner_fraction": 0.8}
    # Tested but significant only with a2: the best partner is named with its q, and the
    # fraction is of the member's lineages.
    assert got["FC"] == {"n_cooccurring_partners": 1, "top_cooccurring_partner": "a2",
                         "top_cooccurring_partner_q": 0.03,
                         "top_cooccurring_partner_fraction": 0.1}
    assert set(got) == {"FA", "FC"}, "a sequence outside the dark families has no row"
