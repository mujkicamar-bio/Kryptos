TIERS = cascade["tiers"]
TIER_IDS = [t["id"] for t in TIERS]
TIER_BY_ID = {t["id"]: t for t in TIERS}


def tier_query(wc):
    """A tier's query set is the previous tier's unresolved output.

    The FIRST tier queries the search representatives (cascade_selection).
    """
    i = TIER_IDS.index(wc.tier)
    if i == 0:
        return f"{OUT}/03_dereplication/search_representatives.faa"
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
