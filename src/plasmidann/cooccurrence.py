"""Dark families that travel together (S8g, dark_cooccurrence).

THE QUESTION

Do two dark families occur together more often than chance predicts? A pair that does is a
candidate for a shared function or a shared mobile unit, and neither family has a name to
say so.

THE DEFINITION

  together  both families have a member ORF on the SAME plasmid. Being in the same
            lineage on different plasmids is not being together.
  unit      the Stage 6 lineage. A lineage counts as together when any of its plasmids
            carries both, so twenty redeposited copies of one plasmid are one observation.
  family    a dark family at the primary resolution, through its dark members.

THE TEST

The hypergeometric upper tail P(X >= k) over lineages: N lineages in the analysis set, K
carry family A, n carry family B, and k are lineages where A and B share a plasmid. Under
independence the lineages carrying B are a random n of the N, and the number of them that
also carry A is hypergeometric. k counts only lineages where the two share a PLASMID, which
is at most the number of lineages carrying both, so the test is conservative with respect
to lineage-level overlap.

Only pairs together in at least `min_lineages_together` lineages (2) are tested. A pair
together in one lineage is one observation of two genes on one plasmid, and every plasmid
carries many families; for two families seen once each, P(X >= 1) = 1/N, which would
declare every pair of singletons on a plasmid significant. The same minimum governs
synteny (targets.yaml synteny.min_lineages): conservation across one lineage is one
observation. Benjamini-Hochberg (Benjamini & Hochberg 1995, J R Stat Soc B 57:289) is
applied across the tested pairs.

Computed in log space with math.lgamma: scipy is not a declared dependency of the
pipeline environment, and the tails reach far below the smallest double (1/C(1000, 50) is
about 1e-85; at full scale much smaller).
"""
import collections
import itertools
import math

from darkorf import status


def _log_comb(a, b):
    return math.lgamma(a + 1) - math.lgamma(b + 1) - math.lgamma(a - b + 1)


def _log_sum_decreasing(terms):
    """log(sum(exp(t))) for log terms that decrease from the first, stopping at the first
    term too small to change the sum at double precision."""
    first = total = None
    for t in terms:
        if first is None:
            first, total = t, 1.0
            continue
        x = math.exp(t - first)
        total += x
        if x < total * 1e-17:
            break
    return first + math.log(total)


def hypergeom_sf(k, N, K, n):
    """P(X >= k) for X hypergeometric: n draws from N items of which K are marked.

    The tail away from the mode is summed, because its terms decrease monotonically: the
    upper tail directly when k is above the mode, otherwise 1 minus the lower tail, whose
    value is then not small, so the subtraction loses nothing that matters.
    """
    lo, hi = max(0, n + K - N), min(K, n)
    if k <= lo:
        return 1.0
    if k > hi:
        return 0.0
    log_total = _log_comb(N, n)

    def log_pmf(i):
        return _log_comb(K, i) + _log_comb(N - K, n - i) - log_total

    mode = (n + 1) * (K + 1) // (N + 2)
    if k > mode:
        return min(1.0, math.exp(_log_sum_decreasing(log_pmf(i) for i in range(k, hi + 1))))
    lower = math.exp(_log_sum_decreasing(log_pmf(i) for i in range(k - 1, lo - 1, -1)))
    return min(1.0, max(0.0, 1.0 - lower))


def benjamini_hochberg(pvalues):
    """Benjamini-Hochberg q-values, in the input order (step-up, monotone, capped at 1)."""
    m = len(pvalues)
    order = sorted(range(m), key=lambda i: pvalues[i])
    q = [0.0] * m
    running = 1.0
    for rank in range(m, 0, -1):
        i = order[rank - 1]
        running = min(running, pvalues[i] * m / rank)
        q[i] = running
    return q


def cooccurrence(plasmid_families, lineage_of, n_lineages, min_lineages_together):
    """Every pair of families together in >= min_lineages_together lineages, tested.

    plasmid_families  {plasmid_id: {family_id, ...}} - the dark families with a member ORF
                      on each plasmid
    lineage_of        {plasmid_id: lineage}; every plasmid above must have one
    n_lineages        N, the lineages in the analysis set

    Returns one dict per tested pair (family_a < family_b), sorted by p-value.
    """
    missing = sorted(p for p in plasmid_families if p not in lineage_of)
    if missing:
        raise ValueError(f"{len(missing)} plasmid(s) carrying a dark family have no Stage 6 "
                         f"lineage, e.g. {missing[:3]}")

    by_lineage = collections.defaultdict(list)
    for plasmid, fams in plasmid_families.items():
        if fams:
            by_lineage[lineage_of[plasmid]].append(fams)
    lineages_of_family = collections.Counter()
    for plasmids in by_lineage.values():
        lineages_of_family.update(set().union(*plasmids))

    # A family in fewer lineages than the minimum cannot be in a tested pair, so it is
    # left out of the enumeration: this is exact, and it is what keeps the pair count
    # bounded at full scale, where most dark families are in one lineage.
    eligible = {f for f, c in lineages_of_family.items() if c >= min_lineages_together}
    together = collections.Counter()
    for plasmids in by_lineage.values():
        pairs = set()
        for fams in plasmids:
            present = sorted(fams & eligible)
            pairs.update(itertools.combinations(present, 2))
        together.update(pairs)

    rows = []
    for (a, b), k in together.items():
        if k < min_lineages_together:
            continue
        K, n = lineages_of_family[a], lineages_of_family[b]
        rows.append({
            "family_a": a, "family_b": b,
            "n_lineages_a": K, "n_lineages_b": n,
            "n_lineages_together": k, "n_lineages_total": n_lineages,
            "fraction_of_a": round(k / K, 4), "fraction_of_b": round(k / n, 4),
            "expected_together": round(K * n / n_lineages, 4),
            "p_value": hypergeom_sf(k, n_lineages, K, n),
        })
    for row, q in zip(rows, benjamini_hochberg([r["p_value"] for r in rows])):
        row["q_value"] = q
        row["status"] = status.SUCCESS
    rows.sort(key=lambda r: (r["p_value"], r["family_a"], r["family_b"]))
    return rows


def family_partners(rows, fdr):
    """Per family: how many partners co-occur at q <= fdr, and the best partner.

    The best partner is the one with the lowest q (then lowest p, then most lineages
    together), named whether or not it is significant, with the fraction of THIS family's
    lineages in which the partner shares a plasmid.
    """
    best, n_sig = {}, collections.Counter()
    for r in rows:
        for focal, partner, fraction in ((r["family_a"], r["family_b"], r["fraction_of_a"]),
                                         (r["family_b"], r["family_a"], r["fraction_of_b"])):
            n_sig[focal] += r["q_value"] <= fdr
            key = (r["q_value"], r["p_value"], -r["n_lineages_together"], partner)
            if focal not in best or key < best[focal][0]:
                best[focal] = (key, partner, r["q_value"], fraction)
    return {f: {"n_cooccurring_partners": n_sig[f], "top_cooccurring_partner": partner,
                "top_cooccurring_partner_q": q, "top_cooccurring_partner_fraction": fraction}
            for f, (_, partner, q, fraction) in best.items()}
