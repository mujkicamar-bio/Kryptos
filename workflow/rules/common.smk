TIERS = cascade["tiers"]
TIER_IDS = [t["id"] for t in TIERS]
TIER_BY_ID = {t["id"]: t for t in TIERS}


def tier_query(wc):
    """A tier's query set is the previous tier's unresolved output, WITHIN THE SAME SHARD.

    The FIRST tier queries results/05_annotation_cascade/input/{cshard}.faa, which is the unique proteins plus
    the spiked positive controls, split into independently resumable units. Controls must
    traverse the identical code path - the same narrowing, the same shards, the same
    thresholds - or the gate at S5 would be testing a different pipeline from the one that
    produced the results.

    A protein never crosses shards. That is what keeps the explained fraction correct: the
    spans of a protein accumulate tier by tier inside one shard, so no tier ever decides
    whether to keep searching a protein against a span set assembled somewhere else.
    """
    i = TIER_IDS.index(wc.tier)
    if i == 0:
        return f"{OUT}/05_annotation_cascade/input/{wc.cshard}.faa"
    return f"{OUT}/05_annotation_cascade/{TIER_IDS[i - 1]}/{wc.cshard}/unresolved.faa"


def tier_spans(wc):
    """Cumulative explained spans from the previous tier of the same shard; empty first."""
    i = TIER_IDS.index(wc.tier)
    if i == 0:
        return []
    return f"{OUT}/05_annotation_cascade/{TIER_IDS[i - 1]}/{wc.cshard}/spans.tsv"
