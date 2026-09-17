"""Significance for a genomic-context association.

WHAT WAS WRONG WITH THE PREVIOUS STATISTIC

Context enrichment was observed rate divided by background rate, and nothing else. Three
consequences:

  * No test. A family with three members all sitting beside a relaxase produced the same
    number as a family with three hundred. The output could not distinguish a coincidence
    from a result, and the ranking that selects 1,000 constructs for synthesis read that
    number.

  * The wrong unit. The rate was over family MEMBERS. Members of a protein family are
    homologs, and they sit on plasmids that are frequently near-identical - the same
    clinical plasmid sequenced forty times is forty members and one observation. Counting
    members treats sequencing effort as evidence.

  * No correction. Six hand-picked features needed none. The label vocabulary is open, so
    a family is now tested against every category present in its neighbourhoods, which is
    thousands of tests, where an uncorrected p-value means nothing.

THE UNIT

The unit of observation is the PLASMID, not the family member. A family is counted as
carrying a category if at least one of its members on that plasmid has the category in its
plus or minus three neighbourhood, and the background is the fraction of all plasmids where
any ORF has it. This removes copy-number inflation within a plasmid.

It does NOT remove clonal redundancy between plasmids: forty independent depositions of the
same clinical plasmid remain forty units. That is a separate correction, and the registry
for it already exists in workflow/scripts/clonal_registry.py; applying it here is recorded
as outstanding in docs/PARAMETER_PROVENANCE.md rather than silently half-done.

THE TEST

Fisher's exact test on the two-by-two table of (family plasmids, rest of the corpus) against
(category present, absent). Exact rather than chi-squared because many categories are rare,
and the expected counts in a rare cell are far below the five that the chi-squared
approximation needs.

One-sided, testing for over-representation: a category a dark family AVOIDS is not a
screening hypothesis, and a two-sided test would spend half its power looking for one.

Multiple testing is controlled by the Benjamini-Hochberg procedure (Benjamini and Hochberg
1995, J R Stat Soc B 57:289), which controls the false discovery rate. FDR rather than
family-wise error because the purpose is to rank candidates for a 1,000-construct screen,
where a false positive costs a well and a false negative costs a discovery.
"""
from scipy.stats import fisher_exact

from darkorf import status


def fisher_enrichment(k, n, K, N, min_units=2):
    """Test whether a category is over-represented in one family's neighbourhoods.

    `k` of `n` plasmids carrying the family have the category nearby; `K` of `N` plasmids in
    the corpus have it. Returns the rates, the enrichment ratio, the odds ratio, a p-value
    and a status from the pipeline's single status vocabulary.

    A family on fewer than `min_units` plasmids is not tested. One plasmid is an anecdote,
    and returning a p-value for it would let a family seen once rank beside one seen on two
    hundred independent plasmids - which is the error the family-level statistic exists to
    prevent.
    """
    if not 0 <= k <= n:
        raise ValueError(f"family counts are impossible: k={k}, n={n}")
    if not 0 <= K <= N:
        raise ValueError(f"corpus counts are impossible: K={K}, N={N}")

    observed = k / n if n else 0.0
    background = K / N if N else 0.0
    # The enrichment is always finite. The family's plasmids are a subset of the corpus, so
    # K >= k: a zero background forces k = 0, and there is no observation to divide by it.
    # Zero over zero is no signal and no surprise.
    ratio = observed / background if background else 1.0

    result = {
        "observed_rate": round(observed, 4),
        "background_rate": round(background, 6),
        "enrichment": ratio,
        "odds_ratio": "",
        "p_value": "",
        "status": status.SUCCESS,
    }

    if n < min_units:
        result["status"] = status.TOO_FEW_MEMBERS
        return result
    if N - n <= 0:
        # The family covers the whole corpus, so there is no background to test against.
        result["status"] = status.NOT_APPLICABLE
        return result

    # The family's plasmids against the REST of the corpus. Testing against a background
    # the family contributes to shrinks any real effect, and for a family covering most of
    # the corpus it shrinks it to nothing.
    rest_with = K - k
    rest_without = (N - n) - rest_with
    if rest_with < 0 or rest_without < 0:
        raise ValueError(
            f"corpus counts are inconsistent with the family counts: k={k}, n={n}, "
            f"K={K}, N={N} - the family's plasmids are not a subset of the corpus")

    odds, p_value = fisher_exact([[k, n - k], [rest_with, rest_without]],
                                 alternative="greater")
    result["odds_ratio"] = round(float(odds), 4) if odds == odds else ""
    result["p_value"] = f"{p_value:.6g}"
    return result


def benjamini_hochberg(p_values):
    """False-discovery-rate q-values, in the input order.

    Benjamini and Hochberg 1995, J R Stat Soc B 57:289. FDR rather than a family-wise
    correction because these values rank candidates for a 1,000-construct screen: a false
    positive costs one well, a false negative costs a discovery.

    Entries that are not numbers - a category that could not be tested - pass through
    unchanged. Treating an empty value as zero would give an untested category the
    strongest q-value in the table.

    The order is preserved because the caller joins the result back to its rows by
    position; returning a sorted list would attach every q-value to the wrong category.
    """
    indexed = []
    for i, value in enumerate(p_values):
        try:
            indexed.append((float(value), i))
        except (TypeError, ValueError):
            continue

    out = list(p_values)
    m = len(indexed)
    if not m:
        return out

    indexed.sort()
    # Walk from the largest p-value down, keeping the running minimum. This enforces the
    # monotonicity the procedure requires: the q-value of the i-th smallest p is the
    # minimum over j >= i of (m / j) * p_j.
    running = 1.0
    for rank in range(m, 0, -1):
        p, i = indexed[rank - 1]
        running = min(running, p * m / rank)
        out[i] = round(min(running, 1.0), 6)
    return out
