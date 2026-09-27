"""Annotation cascade: what a protein is, and how much of it is explained.

Every function answers one of three questions about a single protein: did a hit name it
(is_informative), how much of it do the naming hits cover (explained_fraction), and which
class does it therefore fall in (classify). The E-value decides whether a hit exists at
all; coverage merged over every naming hit of the protein then decides the class. A hit
with an uninformative label ("hypothetical protein") explains nothing and never stops the
cascade, but it is kept as evidence that the protein is real. The class thresholds arrive
as arguments from config/cascade.yaml and are written into the output rows they govern.
"""
import re

# ---------------------------------------------------------------------------------
# Labels: telling an annotation apart from a record that someone else saw this protein
# ---------------------------------------------------------------------------------

# A hit to any of these is a HOMOLOG, not an ANNOTATION.
#
# Recall matters more than precision: a missed pattern promotes an unknown protein to
# FUNCTIONAL and removes it from the dark set, while a false positive keeps a named protein
# in it. Checked case by case in tests/test_uninformative_labels.py.
UNINFORMATIVE = re.compile(
    r"""
      hypothetical                # 'hypothetical protein', the canonical case
    | uncharacteri[sz]ed          # both spellings; UniProt and NCBI differ
    | \bDUF\d*\b                  # DUF1234, and bare 'DUF domain-containing protein'
    | \bUPF\d+                    # UniProt uncharacterized protein family
    | unknown\ function
    | unnamed\ protein            # 'unnamed protein product', common in older GenBank
    # 'putative', 'predicted' or 'conserved protein' only as the whole name: at the start
    # of the label, after the accession (the first token holding a digit), after ':' or
    # 'Full=', and followed by the end, ';' or ' [organism]'. 'putative protein kinase'
    # names a function.
    | (?:^|^\S*\d\S*\s+|:\s*|Full=)
      (?:putative|predicted|conserved)\ protein
      \s*(?:$|;|\[)
    """,
    re.I | re.X,
)


# An informative name that describes a DOMAIN rather than the protein. NCBI's PGAP names a
# protein "X domain-containing protein" or "X family protein" when the evidence is a
# domain-level or family-level HMM, not a full-length functional assignment (Li W. et al.
# 2021, Nucleic Acids Res. 49:D1020). Such a name still says something, so it is
# informative and explains its span; but on its own it makes a protein DOMAIN_ONLY, never
# FUNCTIONAL. Measured on the nr benchmark: 286 of 1,887 FUNCTIONAL calls at T5 (15.2%)
# rested on "domain-containing" alone and 262 more on "family protein".
DOMAIN_NAMED = re.compile(r"domain[- ]containing\ protein|\bfamily\ protein\b", re.I | re.X)


def is_informative(label):
    """True when a label names a function, rather than recording that someone saw it.

    A missing or blank label is uninformative: a parser that produced no label must leave
    the protein unknown, not promote it out of the dark set.
    """
    if not (label or "").strip():
        return False
    return not UNINFORMATIVE.search(label)


# ---------------------------------------------------------------------------------
# Significance: whether a hit is allowed to exist at all
# ---------------------------------------------------------------------------------

def passes_significance(evalue, max_evalue):
    """Whether a hit passes the tier's E-value floor.

    None means the tier has no floor (T1 uses Pfam's curated GA thresholds). An unparsable
    E-value fails when a floor is set.
    """
    if max_evalue is None:
        return True
    try:
        return float(evalue) <= max_evalue
    except (TypeError, ValueError):
        return False


def as_float(value):
    """An E-value as a number, or infinity when it cannot be parsed.

    Infinity rather than an exception keeps bit-score-only tiers (T1) sortable: when no hit
    in a set has an E-value they all tie, and the caller's secondary key - cascade order,
    which is authority order - decides. Shared with the stages that rank hits, so that
    "unparsable sorts last" is one rule rather than a copy per script.
    """
    try:
        return float(value)
    except (TypeError, ValueError):
        return float("inf")


def _evalue_of(hit):
    """Sort key for hit significance. A hit without a parsable E-value sorts last."""
    return as_float(hit.get("evalue"))


# ---------------------------------------------------------------------------------
# Classification
# ---------------------------------------------------------------------------------

def classify(hits, explained, min_coverage, tier_order):
    """Resolve one protein against all of its cascade hits.

    hits         every significant hit for this protein, from every tier
    explained    merged coverage of the protein by its informative hits (explained_fraction)
    min_coverage explained fraction from which the protein counts as FUNCTIONAL
    tier_order   the configured tier ids in cascade order, which is authority order

    Class: FUNCTIONAL when `explained` reaches min_coverage or an informative hit has no
    span; DOMAIN_ONLY otherwise, and whenever every informative name is domain-level
    (DOMAIN_NAMED); UNCHARACTERIZED_HOMOLOG when no hit is informative; NONE without hits.
    Label: from the most authoritative tier that named the protein, and within it the hit
    with the best E-value, so that automated transfer does not outrank curated evidence;
    the label with the best E-value over all tiers is reported beside it.
    """
    if not hits:
        return {"annot_tier": None, "annot_label": None, "functional_class": "NONE",
                "homology_depth": None, "annot_qcov": None, "annot_tcov": None,
                "annot_evalue": None, "n_informative_hits": 0,
                "best_evalue_label": None, "best_evalue_tier": None,
                "named_by_domain_only": 0, "span_measured": 1}

    informative = [h for h in hits if is_informative(h.get("label"))]
    # A hit without coordinates comes from the pharokka tier, which reports no alignment
    # span. It is taken as a family assignment of the whole protein and makes the protein
    # FUNCTIONAL on its own; span_measured then records that completeness was not measured.
    family_level = [h for h in informative if not _has_span(h)]
    domain_only = bool(informative) and all(DOMAIN_NAMED.search(h["label"])
                                            for h in informative)
    if informative:
        best = min(informative, key=lambda h: (tier_order.index(h["tier"]), _evalue_of(h)))
        strongest = min(informative,
                        key=lambda h: (_evalue_of(h), tier_order.index(h["tier"])))
        if (explained >= min_coverage or family_level) and not domain_only:
            cls = "FUNCTIONAL"
        else:
            cls = "DOMAIN_ONLY"
        depth_from = informative
    else:
        # Nothing named it anywhere. Show the most authoritative record of having seen it.
        best = strongest = min(hits, key=lambda h: (tier_order.index(h["tier"]),
                                                    _evalue_of(h)))
        cls = "UNCHARACTERIZED_HOMOLOG"
        depth_from = hits

    return {
        # annot_tier: the most authoritative tier that named the protein.
        "annot_tier": best["tier"],
        "annot_label": best["label"],
        "functional_class": cls,
        # homology_depth: the shallowest tier with an informative hit, or with any hit when
        # none is informative - how far the cascade had to search to recognise the protein.
        "homology_depth": min(tier_order.index(h["tier"]) for h in depth_from) + 1,
        "annot_qcov": best["coverage"],
        "annot_tcov": best.get("target_coverage"),
        "annot_evalue": best.get("evalue"),
        # Informative hits over all tiers, as the tools reported them: a domain that T1
        # and T2 both report counts twice.
        "n_informative_hits": len(informative),
        "best_evalue_label": strongest["label"],
        "best_evalue_tier": strongest["tier"],
        "named_by_domain_only": int(domain_only),
        # 0 when the class rests on span-less hits alone: explained is then 0 because
        # nothing measured it, and cascade_resolve reports completeness as NOT_MEASURED.
        "span_measured": int(not (family_level and explained == 0)),
    }


def _has_span(hit):
    """Whether a hit carries alignment coordinates.

    hits.tsv writes start and end empty for a tier that reports no span, and csv.DictReader
    returns them as ''. A hit without the keys at all counts as spanned.
    """
    return hit.get("start", 0) != "" and hit.get("end", 0) != ""


# ---------------------------------------------------------------------------------
# Narrowing: what each tier hands the next
# ---------------------------------------------------------------------------------

def narrow_by_explained(all_ids, explained, threshold):
    """Ids still worth searching: those explained less than `threshold`.

    The caller passes narrow_at, how explained a protein must be for the search on it to
    stop. Every protein explained below narrow_at reaches every tier, except where a tier skips it:
    T5 skips proteins Pfam or Swiss-Prot named (skip_if_named_by) and the DIAMOND tiers
    skip artefact-flagged proteins. Narrowing on explained fraction rather than on "any
    hit" keeps a protein with a 15% domain match in the search for the other 85%.
    """
    return [i for i in all_ids if explained.get(i, 0.0) < threshold]


# ---------------------------------------------------------------------------------
# Coverage: how much of a protein anything can account for
# ---------------------------------------------------------------------------------

def explained_fraction(length, intervals):
    """Fraction of a protein covered by the union of a set of alignment intervals.

    Overlapping intervals are MERGED, not summed: a protein hit twice by the same family
    at 8-337 and 263-387 is 97% explained, not 168%. Summing would let a protein exceed
    full coverage and would make the narrowing threshold meaningless.

    Given the spans of informative hits it is the explained fraction of an annotated
    protein; given uninformative_spans(hits) it is the dark_covered_fraction of a dark one.

    The algorithm is a single sweep over sorted intervals, tracking the furthest right
    edge seen so far; `end` starts at 0 because coordinates are 1-based and inclusive.
    """
    if not intervals or not length:
        return 0.0
    covered, end = 0, 0
    for s, e in sorted(intervals):
        # Clip this interval to the part not already covered.
        s, e = max(s, end + 1), max(e, end)
        if e >= s:
            covered += e - s + 1
            end = e
    return round(covered / length, 4)


def completeness(fraction, full_at, partial_at):
    """Band an explained fraction into FULL, PARTIAL, FRAGMENT or NONE.

    For the dark set the explained fraction is 0 by construction (a dark protein has no
    informative span), so this band is NONE for every dark protein; dark_covered_fraction
    is the discriminating measure there.
    """
    if fraction >= full_at:
        return "FULL"
    if fraction >= partial_at:
        return "PARTIAL"
    if fraction > 0:
        return "FRAGMENT"
    return "NONE"


def uninformative_spans(hits):
    """Alignment intervals of the hits that named nothing - the dark protein's coverage.

    95% of a protein covered by 'hypothetical protein' hits is a conserved, full-length
    protein nobody has named; one 20-aa fragment hit is weak evidence. Fed through
    explained_fraction these spans give dark_covered_fraction, which separates the two.
    """
    return [(h["start"], h["end"]) for h in hits if not is_informative(h.get("label"))]


def n_dark_databases(hits):
    """How many independent tiers called this protein something uninformative.

    Counted per TIER, not per hit: two 'hypothetical protein' hits from the same DIAMOND
    search against the same database are one observation, not two. Three independent
    databases agreeing that a protein is real and unnamed is substantially stronger
    evidence than one, and this is the number that expresses it.
    """
    return len({h["tier"] for h in hits if not is_informative(h.get("label"))})


# ---------------------------------------------------------------------------------
# The dark evidence ladder
# ---------------------------------------------------------------------------------

# How strongly the databases support this being a real protein, kept separate from the
# fact that its function is unknown. EVERY rung means "function unknown"; they differ only
# in the risk that there is no protein there at all. The column is descriptive, never a
# filter: a CURATED_FAMILY protein is arguably less novel, since Pfam already
# recognised the family. How widely a protein occurs is measured by the pipeline itself,
# from the plasmids its family occurs on, not from database title prefixes.
_EVIDENCE_RUNGS = [
    # A curator built and named a family for it. Strongest evidence of reality.
    ("CURATED_FAMILY", re.compile(r"\bDUF\d*\b|\bUPF\d+", re.I)),
    # Homologs exist; no function.
    ("CONSERVED", re.compile(r"\bconserved\b", re.I)),
    # One algorithm's output and nothing more.
    ("PREDICTED_ONLY", re.compile(r".")),
]


def dark_evidence(labels):
    """Strongest evidence rung across every uninformative label a protein collected.

    A protein called 'hypothetical protein' at one tier and 'DUF1234' at another is a DUF
    protein: the strongest rung any tier reached wins, because the rungs describe evidence
    of existence and evidence accumulates.
    """
    if not labels:
        return "NONE"
    for rung, pattern in _EVIDENCE_RUNGS:
        if any(pattern.search(l or "") for l in labels):
            return rung
    return "NONE"


# ---------------------------------------------------------------------------------
# Threshold coherence
# ---------------------------------------------------------------------------------

REQUIRED_THRESHOLDS = ("narrow_at", "min_coverage", "full_at", "partial_at")


def check_thresholds(cfg):
    """Validate relationships between thresholds that a JSON schema cannot express.

    The schema checks that each threshold is a number in [0, 1]. It cannot check that the
    numbers make sense TOGETHER, which is where the interesting failures live. Raises
    ValueError with a message naming the offending pair.

    narrow_at >= min_coverage
        narrow_at withholds a protein from deeper tiers; min_coverage decides whether it is
        FUNCTIONAL. If narrowing were the stricter, the pipeline would stop searching
        proteins it then reports as DOMAIN_ONLY, without the deeper tiers' results.

    full_at >= partial_at
        Otherwise the completeness bands overlap and FULL becomes unreachable.
    """
    missing = [k for k in REQUIRED_THRESHOLDS if k not in cfg]
    if missing:
        raise ValueError(f"cascade config is missing threshold(s): {', '.join(missing)}")
    if cfg["narrow_at"] < cfg["min_coverage"]:
        raise ValueError(
            f"narrow_at ({cfg['narrow_at']}) < min_coverage ({cfg['min_coverage']}): the "
            "cascade would stop searching proteins it then reports as DOMAIN_ONLY")
    if cfg["full_at"] < cfg["partial_at"]:
        raise ValueError(
            f"full_at ({cfg['full_at']}) < partial_at ({cfg['partial_at']}): completeness "
            "bands overlap and FULL is unreachable")


# Allowed relative difference between hmmer_z and the unique-protein count. 2% is small
# against the orders of magnitude an E-value spans.
HMMER_Z_TOLERANCE = 0.02


def check_hmmer_z(declared, actual):
    """Raise when hmmer_z differs from the unique-protein count by more than
    HMMER_Z_TOLERANCE; the message names the correct value."""
    if not actual:
        raise ValueError("cannot check hmmer_z against an empty protein set")
    drift = abs(declared - actual) / actual
    if drift > HMMER_Z_TOLERANCE:
        raise ValueError(
            f"hmmer_z in config/cascade.yaml is {declared}, but the analysis set holds "
            f"{actual} unique proteins ({drift:.1%} drift, tolerance "
            f"{HMMER_Z_TOLERANCE:.0%}). Set hmmer_z: {actual} and re-run.")
