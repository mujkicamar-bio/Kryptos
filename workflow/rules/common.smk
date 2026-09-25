TIERS = cascade["tiers"]
TIER_IDS = [t["id"] for t in TIERS]
TIER_BY_ID = {t["id"]: t for t in TIERS}


def tier_query(wc):
    """A tier's query set is the previous tier's unresolved output.

    The FIRST tier queries results/03_dereplication/cascade_input.faa, which is the unique
    proteins plus the spiked controls. Controls must traverse the identical code path - the
    same narrowing, the same thresholds - or the gate at S5 would be testing a different
    pipeline from the one that produced the results.
    """
    i = TIER_IDS.index(wc.tier)
    if i == 0:
        return f"{OUT}/03_dereplication/cascade_input.faa"
    return f"{OUT}/05_annotation_cascade/{TIER_IDS[i - 1]}/unresolved.faa"


def tier_named(wc):
    """Every earlier tier's hits, for a tier that skips proteins they already named
    (skip_if_named_by); empty otherwise."""
    if not TIER_BY_ID[wc.tier].get("skip_if_named_by"):
        return []
    i = TIER_IDS.index(wc.tier)
    return [f"{OUT}/05_annotation_cascade/{t}/hits.tsv" for t in TIER_IDS[:i]]


def tier_spans(wc):
    """Cumulative explained spans from the previous tier; empty for the first."""
    i = TIER_IDS.index(wc.tier)
    if i == 0:
        return []
    return f"{OUT}/05_annotation_cascade/{TIER_IDS[i - 1]}/spans.tsv"
