"""Dark protein sequences that travel together (S8g), counted over Mash lineages.

The unit is the unique dark protein sequence (seq_id), not the family. Two sequences are
together in a lineage (plasmid_lineage.tsv) when one of its plasmids carries an ORF of
each, so redeposited copies of one plasmid are one observation.

Null model: the K lineages carrying sequence A and the n lineages carrying B are
independent, uniformly random subsets of the N lineages. Each pair is tested by the
hypergeometric upper tail P(X >= k), where k is the number of lineages in which A and B
share a plasmid. Every lineage is taken as equally likely to carry a sequence, although
lineages differ in gene content, so two sequences confined to large plasmids share
lineages more often than this null predicts.

Only sequences in at least `min_lineages_together` lineages enter the enumeration: a pair
needs k >= min_lineages_together to be reported, and k cannot exceed either sequence's
lineage count. This keeps the pair count bounded at full scale, where most of the unique
proteins are in one lineage. Every pair of such sequences that shares a plasmid enters the
Benjamini-Hochberg correction (Benjamini & Hochberg 1995, J R Stat Soc B 57:289); only pairs
with k >= min_lineages_together are reported. The cut on k is applied after the correction
because k is the test statistic, and filtering on it first would shrink the q-values
(Bourgon et al. 2010, PNAS 107:9546).

Computed in log space with math.lgamma: the tails fall below the smallest double, and scipy
is not a dependency of the pipeline environment.
"""
import collections
import itertools
import math


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


def cooccurrence(plasmid_seqs, lineage_of, n_lineages, min_lineages_together):
    """The pairs of sequences together in >= min_lineages_together lineages, tested.

    plasmid_seqs  {plasmid_id: {seq_id, ...}} - the dark sequences with an ORF on each
                  plasmid
    lineage_of    {plasmid_id: lineage}; every plasmid above must have one
    n_lineages    N, the lineages in the analysis set

    Returns (rows, stats): one dict per reported pair (seq_a < seq_b), sorted by p-value;
    stats counts the sequences enumerated, the pair occurrences on plasmids (an upper
    bound on the pairs held in memory) and the pairs tested.
    """
    missing = sorted(p for p in plasmid_seqs if p not in lineage_of)
    if missing:
        raise ValueError(f"{len(missing)} plasmid(s) carrying a dark sequence have no "
                         f"lineage, e.g. {missing[:3]}")

    by_lineage = collections.defaultdict(list)
    for plasmid, seqs in plasmid_seqs.items():
        if seqs:
            by_lineage[lineage_of[plasmid]].append(seqs)
    lineages_of_seq = collections.Counter()
    for plasmids in by_lineage.values():
        lineages_of_seq.update(set().union(*plasmids))

    eligible = {s for s, c in lineages_of_seq.items() if c >= min_lineages_together}
    together = collections.Counter()
    n_occurrences = 0
    for plasmids in by_lineage.values():
        pairs = set()
        for seqs in plasmids:
            present = sorted(seqs & eligible)
            n_occurrences += len(present) * (len(present) - 1) // 2
            pairs.update(itertools.combinations(present, 2))
        together.update(pairs)

    # Every counted pair enters the correction, but a row is built only for a reported
    # pair, so an unreported pair costs one p-value and one q-value in memory.
    counted = list(together.items())
    pvalues = [hypergeom_sf(k, n_lineages, lineages_of_seq[a], lineages_of_seq[b])
               for (a, b), k in counted]
    rows = []
    for ((a, b), k), p, q in zip(counted, pvalues, benjamini_hochberg(pvalues)):
        if k < min_lineages_together:
            continue
        K, n = lineages_of_seq[a], lineages_of_seq[b]
        rows.append({
            "seq_a": a, "seq_b": b,
            "n_lineages_a": K, "n_lineages_b": n,
            "n_lineages_together": k, "n_lineages_total": n_lineages,
            "fraction_of_a": round(k / K, 4), "fraction_of_b": round(k / n, 4),
            "expected_together": round(K * n / n_lineages, 4),
            "p_value": p, "q_value": q,
        })
    rows.sort(key=lambda r: (r["p_value"], r["seq_a"], r["seq_b"]))
    stats = {"n_sequences": len(lineages_of_seq), "n_enumerated": len(eligible),
             "n_pair_occurrences": n_occurrences, "n_tested": len(counted)}
    return rows, stats


def family_partners(rows, fdr, family_of_seq):
    """Per family: the partner sequences its members co-occur with at q <= fdr, and the
    best partner.

    n_cooccurring_partners counts distinct partner sequences over all the family's members.
    The best partner is the one with the lowest q (then lowest p, then most lineages
    together), named whether or not it is significant, with the fraction of the MEMBER's
    lineages in which the partner shares a plasmid.
    """
    best, significant = {}, collections.defaultdict(set)
    for r in rows:
        for focal, partner, fraction in ((r["seq_a"], r["seq_b"], r["fraction_of_a"]),
                                         (r["seq_b"], r["seq_a"], r["fraction_of_b"])):
            family = family_of_seq.get(focal)
            if family is None:
                continue
            if r["q_value"] <= fdr:
                significant[family].add(partner)
            key = (r["q_value"], r["p_value"], -r["n_lineages_together"], partner)
            if family not in best or key < best[family][0]:
                best[family] = (key, partner, r["q_value"], fraction)
    return {f: {"n_cooccurring_partners": len(significant[f]),
                "top_cooccurring_partner": partner,
                "top_cooccurring_partner_q": q, "top_cooccurring_partner_fraction": fraction}
            for f, (_, partner, q, fraction) in best.items()}
