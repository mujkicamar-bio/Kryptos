TIERS = cascade["tiers"]
TIER_IDS = [t["id"] for t in TIERS]
TIER_BY_ID = {t["id"]: t for t in TIERS}


def tier_query(wc):
    """A tier's query set is the previous tier's unresolved output."""
    i = TIER_IDS.index(wc.tier)
    if i == 0:
        return f"{OUT}/s2/unique_proteins.faa"
    return f"{OUT}/s3/{TIER_IDS[i - 1]}/unresolved.faa"
