import re

UNINFORMATIVE = re.compile(
    r"hypothetical|uncharacteri[sz]ed|\bDUF\d+|unknown function|^ORF$|putative protein",
    re.I,
)
MIN_COVERAGE = 0.5
TIERS = ["T1", "T2", "T3", "T4", "T5", "T6", "T7", "T8"]

_RANK = {"FUNCTIONAL": 0, "DOMAIN_ONLY": 1, "UNCHARACTERIZED_HOMOLOG": 2}


def _class_of(hit):
    if UNINFORMATIVE.search(hit["label"]):
        return "UNCHARACTERIZED_HOMOLOG"
    if hit["coverage"] < MIN_COVERAGE:
        return "DOMAIN_ONLY"
    return "FUNCTIONAL"


def classify(hits):
    """Resolve a protein against its cascade hits.

    functional_class is deliberately not a boolean: a homolog that is itself unnamed
    (UNCHARACTERIZED_HOMOLOG) is the class worth screening - certainly real, certainly unknown.
    Best class wins; ties break to the shallowest tier, since cascade order is authority order.
    """
    if not hits:
        return {"annot_tier": None, "annot_label": None, "functional_class": "NONE",
                "homology_depth": None, "annot_qcov": None, "annot_tcov": None,
                "annot_evalue": None}
    best = min(hits, key=lambda h: (_RANK[_class_of(h)], TIERS.index(h["tier"])))
    return {
        "annot_tier": best["tier"],
        "annot_label": best["label"],
        "functional_class": _class_of(best),
        "homology_depth": TIERS.index(best["tier"]) + 1,
        "annot_qcov": best["coverage"],
        "annot_tcov": best.get("target_coverage"),
        "annot_evalue": best.get("evalue"),
    }


def narrow(all_ids, hit_ids):
    """Ids a tier could not resolve. Refuses hits for ids that were never queried."""
    unknown = set(hit_ids) - set(all_ids)
    if unknown:
        raise ValueError(f"hits for {len(unknown)} id(s) not in the query set, e.g. {sorted(unknown)[:3]}")
    return [i for i in all_ids if i not in hit_ids]


def explained_fraction(length, intervals):
    """Fraction of a protein covered by the union of its domain hits.

    Overlapping domains are merged, not summed: a protein hit twice by the same family
    at 8-337 and 263-387 is 97% explained, not 168%. This is the metric that finds dark
    regions inside otherwise-annotated proteins.
    """
    if not intervals or not length:
        return 0.0
    covered, end = 0, 0
    for s, e in sorted(intervals):
        s, e = max(s, end + 1), max(e, end)
        if e >= s:
            covered += e - s + 1
            end = e
    return round(covered / length, 4)
