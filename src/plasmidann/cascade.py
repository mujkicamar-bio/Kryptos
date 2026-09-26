"""Annotation cascade: deciding what a protein is, and whether it is worth screening.

This module holds the intellectual core of the pipeline. Everything here answers one of
three questions about a single protein:

  1. Did anything *name* it?            -> is_informative, UNINFORMATIVE
  2. How much of it is *explained*?     -> explained_fraction, informative_spans
  3. What class does it therefore fall in, and how deep did we have to dig?
                                        -> classify

Three design commitments run through all of it, each of which was a defect once:

  * A class is a property of the PROTEIN, not of one alignment. A replication initiator
    carrying RepA_N and Bac_RepA_C is completely annotated even though neither domain
    alone covers half of it. Classifying per hit called 17% of Pfam-hit proteins
    DOMAIN_ONLY at >=50% explained, and DOMAIN_ONLY is target-eligible - so the plasmid
    backbone was walking into the screening pool.

  * The E-value gates whether a hit EXISTS; coverage then decides the class among the
    survivors. The dangerous case is a high-coverage hit with a weak E-value: a spurious
    long alignment silently removes a genuine dark protein from the pool, and no
    downstream stage can recover it. That is the worst error this project can make.

  * A "hypothetical protein" hit is a HOMOLOG, not an ANNOTATION. It never explains
    anything, however well it aligns, and it never stops the cascade - but it is recorded
    rather than discarded, because a protein several independent databases call
    hypothetical is real, widespread and genuinely uncharacterised, which is evidence FOR
    screening it, not against.

No function here reads a module-level threshold. Every threshold arrives as an argument
from config/cascade.yaml, is schema-validated, and is stamped into the output row it
governs (design principle P4). Module constants were how three of the six thresholds
escaped the config in v1.
"""
import re

# ---------------------------------------------------------------------------------
# Labels: telling an annotation apart from a record that someone else saw this protein
# ---------------------------------------------------------------------------------

# A hit to any of these is a HOMOLOG, not an ANNOTATION.
#
# Recall matters far more than precision here, and the asymmetry is worth stating: a
# missed pattern silently promotes an unknown protein to FUNCTIONAL and removes it from
# the screening set forever. A false positive merely keeps a named protein in the pool a
# little longer, where later evidence will demote it. Recall is measured against a
# labelled 25-case set in tests/test_uninformative_labels.py and CI fails if it drops.
#
# Patterns are deliberately verbose (re.X) so each can be read and argued with.
UNINFORMATIVE = re.compile(
    r"""
      hypothetical                # 'hypothetical protein', the canonical case
    | uncharacteri[sz]ed          # both spellings; UniProt and NCBI differ
    | \bDUF\d*\b                  # DUF1234, and bare 'DUF domain-containing protein'
    | \bUPF\d+                    # UniProt uncharacterized protein family
    | unknown\ function
    | unnamed\ protein            # 'unnamed protein product', common in older GenBank
    | predicted\ protein          # a gene caller's opinion, not an observation
    | conserved\ protein          # conserved, but conserved as what?
    | ^ORF$
    | putative\ protein
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

    A missing label counts as UNINFORMATIVE, not informative. The inverse - which is what
    v1 did, because `not UNINFORMATIVE.search(None or "")` evaluates to True - silently
    promoted proteins out of the dark set whenever a parser produced no label. Failing
    toward "we do not know" is the safe direction for a discovery pipeline.
    """
    if not label:
        return False
    return not UNINFORMATIVE.search(label)


# ---------------------------------------------------------------------------------
# Significance: whether a hit is allowed to exist at all
# ---------------------------------------------------------------------------------

def passes_significance(evalue, max_evalue):
    """Whether a hit is significant enough to contribute a span, a label or a row.

    This is applied at parse time, before anything else looks at the hit, because an
    insignificant hit must not be able to explain part of a protein and thereby withhold
    it from deeper tiers. Measured on T2: 7.2% of its resolutions depended on domains
    that were not individually significant, because `-E` sets only the SEQUENCE threshold
    and every domain of a passing sequence then lands in the output.

    `max_evalue` of None means "this tier has no E-value criterion" and everything passes.
    That is not laziness: T1 uses Pfam's per-family gathering thresholds, chosen by hand
    by each family's curator. For some short, low-information families GA corresponds to
    an E-value looser than any global cut, so imposing a blanket floor would override
    curation and degrade the highest-quality signal in the whole cascade.

    An unparsable E-value where a floor WAS set is rejected: it means the parser and the
    search tool disagree about the output format, and silently keeping such hits would
    hide that.
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

    Arguments
    ---------
    hits         every hit for this protein, from every tier, already significance-filtered
    explained    merged INFORMATIVE coverage of the protein, from explained_fraction
    min_coverage threshold above which the protein counts as functionally annotated
    tier_order   the configured tier ids in cascade order, e.g. ["T1","T2","T3","T4","T5"]

    The class comes from `explained`, not from any single hit's coverage. This is the fix
    for the backbone-contamination defect described in the module docstring: coverage is
    merged across informative hits before any decision is taken.

    The label comes from the most AUTHORITATIVE tier that named the protein - tier order,
    which config/cascade.yaml declares to be authority order - and within that tier from
    the most significant hit, not the widest one. Spec section 20: automated transfer must
    not outrank curated evidence. Ranking on E-value alone across tiers let an nr free-text
    title take the label from an informative Swiss-Prot hit on 39% of the proteins that had
    one in the nr benchmark. Within a tier, E-value: a longer alignment is not a better
    identification (ABC_membrane at E=1e-23 was once chosen over Peptidase_C39 at
    E=6.5e-40 purely because it aligned further). The label with the best E-value across
    all tiers is still reported, as best_evalue_label and best_evalue_tier.

    A protein whose only informative names are domain-level (DOMAIN_NAMED) is DOMAIN_ONLY
    whatever their coverage; named_by_domain_only records it, so the FUNCTIONAL count can
    be reported with and without the rule.

    functional_class is deliberately not a boolean. UNCHARACTERIZED_HOMOLOG - a protein
    whose only homologs are themselves unnamed - is the class worth screening: certainly
    real, certainly unknown.

    `tier_order` is passed in rather than read from a constant so that a four-tier cascade
    reports depths on a four-tier scale. v1 indexed into a hardcoded eight-tier list.
    """
    if not hits:
        return {"annot_tier": None, "annot_label": None, "functional_class": "NONE",
                "homology_depth": None, "annot_qcov": None, "annot_tcov": None,
                "annot_evalue": None, "n_informative_hits": 0,
                "best_evalue_label": None, "best_evalue_tier": None,
                "named_by_domain_only": 0, "span_measured": 1}

    informative = [h for h in hits if is_informative(h.get("label"))]
    # A hit with no coordinates is a FAMILY-LEVEL assignment from a tool that reports no
    # alignment span - the pharokka tier, whose families are whole-protein clusters and
    # whose raw alignments are deleted on exit. It says the whole protein belongs to a
    # named family. Classing it DOMAIN_ONLY because `explained` is 0 would say "a fragment
    # matched", the opposite of what was reported. So such a hit is FUNCTIONAL on its own,
    # and the row records that its completeness was NOT measured rather than measured as 0.
    family_level = [h for h in informative if not _has_span(h)]
    domain_only = bool(informative) and all(DOMAIN_NAMED.search(h["label"])
                                            for h in informative)
    if informative:
        # The most authoritative tier wins; within it, the strongest hit.
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
        "annot_tier": best["tier"],
        "annot_label": best["label"],
        "functional_class": cls,
        # HOW DEEP THE CASCADE HAD TO DIG, which is not the same question as which hit
        # named the protein best, and answering both with one field was a regression
        # introduced by ranking labels on E-value.
        #
        #   annot_tier      where the winning LABEL came from - the strongest alignment
        #   homology_depth  the SHALLOWEST tier that recognised this protein at all
        #
        # A protein Pfam named at T1 is a shallow, well-characterised protein even when nr
        # later produces a better E-value for the same assignment. Depth is used downstream
        # to describe how obscure a protein is, so it must be a property of the cascade
        # rather than of whichever database happened to align best.
        "homology_depth": min(tier_order.index(h["tier"]) for h in depth_from) + 1,
        "annot_qcov": best["coverage"],
        "annot_tcov": best.get("target_coverage"),
        "annot_evalue": best.get("evalue"),
        # 0.9 explained by one domain and 0.9 explained by six fragments are different
        # claims. Reported so the reader of the table can tell them apart.
        "n_informative_hits": len(informative),
        "best_evalue_label": strongest["label"],
        "best_evalue_tier": strongest["tier"],
        "named_by_domain_only": int(domain_only),
        # 0 when the class rests on a family-level assignment alone: explained_fraction
        # is then 0 because nothing MEASURED it, not because nothing matched, and
        # cascade_resolve reports completeness as NOT_MEASURED on that signal.
        "span_measured": int(not (family_level and explained == 0)),
    }


def _has_span(hit):
    """Whether a hit carries alignment coordinates.

    hits.tsv always has the start and end columns; a tier that reports no span writes them
    EMPTY, and csv.DictReader hands that back as ''. So the signal is an explicitly empty
    coordinate. A hit dict with no such key at all is treated as spanned - that is the
    shape of every hit before the column existed, and of the fixtures written for it.
    """
    return hit.get("start", 0) != "" and hit.get("end", 0) != ""


# ---------------------------------------------------------------------------------
# Narrowing: what each tier hands the next
# ---------------------------------------------------------------------------------

def narrow_by_explained(all_ids, explained, threshold):
    """Ids still worth searching: those not yet explained past `threshold`.

    Narrowing on explained fraction rather than on "got any hit" is what stops a 15%
    domain match from terminating the search over the other 85% of a protein.

    IMPORTANT - which threshold belongs here. The caller must pass `narrow_at`, NOT
    `min_explained`. They are two different numbers doing two different jobs:

      narrow_at (0.7)      how finished a protein must be before we stop searching it
      min_explained (0.5)  how explained a protein must be before we REPORT it as such

    v1 used one number for both. Because a protein removed at T2 has no T5 result, that
    made the threshold unsweepable: the counterfactual does not exist, and answering
    "what if it had been 0.7?" would need a full re-run per value. Measured, 0.3 -> 0.8
    moved the deep-tier set by 76%, and 0.5 sat exactly at the 25th percentile of the
    observed distribution - the densest possible place to put a hard cut.

    Keeping narrow_at above min_explained means every protein explained below narrow_at
    was searched by every tier, so min_explained can be swept up to narrow_at without a
    re-run; the sweep cohort measures what stopping at narrow_at costs. narrow_at 0.7 is
    a user decision (2026-09-25). See also check_thresholds, which enforces narrow_at >=
    min_explained.
    """
    unknown = set(explained) - set(all_ids)
    if unknown:
        raise ValueError(f"explained fractions for {len(unknown)} id(s) never queried, "
                         f"e.g. {sorted(unknown)[:3]}")
    return [i for i in all_ids if explained.get(i, 0.0) < threshold]


# ---------------------------------------------------------------------------------
# Coverage: how much of a protein anything can account for
# ---------------------------------------------------------------------------------

def explained_fraction(length, intervals):
    """Fraction of a protein covered by the union of a set of alignment intervals.

    Overlapping intervals are MERGED, not summed: a protein hit twice by the same family
    at 8-337 and 263-387 is 97% explained, not 168%. Summing would let a protein exceed
    full coverage and would make the narrowing threshold meaningless.

    This is the metric that finds dark regions inside otherwise-annotated proteins, and it
    is used for two different populations depending on which spans are passed in:

      informative_spans(hits)   -> explained_fraction, for annotated proteins
      uninformative_spans(hits) -> dark_covered_fraction, for dark proteins

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
    """Band an explained fraction. A 0.51 explanation is not a 0.99 one.

    Thresholds are arguments rather than module constants because in v1 they were
    FULL_AT/PARTIAL_AT literals in this file - two of the three thresholds that escaped
    config and schema entirely, and so could never be swept or recorded (P4).

    Note this field is informative about ANNOTATED proteins only. For the dark set it is
    constant NONE by construction, since explained_fraction there is built from
    informative spans and a dark protein has none. Measured 71/71 and 703/703. The
    discriminating axis for dark proteins is dark_covered_fraction.
    """
    if fraction >= full_at:
        return "FULL"
    if fraction >= partial_at:
        return "PARTIAL"
    if fraction > 0:
        return "FRAGMENT"
    return "NONE"


def informative_spans(hits):
    """Alignment intervals from hits that actually name a function.

    A 'hypothetical protein' hit explains nothing, however well it aligns, so it must not
    count toward explained_fraction and must not stop the cascade. Another database may
    still name the protein, and the whole point of the cascade is to give it that chance.
    """
    return [(h["start"], h["end"]) for h in hits if is_informative(h.get("label"))]


def uninformative_spans(hits):
    """Alignment intervals from hits that named nothing - the dark protein's own coverage.

    This is the evidence v1 threw away. Uninformative hits were kept as LABELS but their
    coordinates were discarded, which collapsed a real distinction:

      95% of the protein covered by 'hypothetical protein' across three databases
        -> a real, conserved, full-length protein that nobody has named. Strong target.

      one 20-aa 'hypothetical' fragment hit
        -> weak evidence, possibly a spurious call. Weak target.

    Fed through explained_fraction this yields dark_covered_fraction, the analogue of
    annot_completeness for the population the pipeline actually exists to characterise.
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

# How strongly the databases support this being a real protein, kept deliberately separate
# from the fact that its function is unknown. EVERY rung means "function unknown"; they
# differ only in the risk that there is no protein there at all.
#
# WEAK SIGNAL ONLY. This is a low-weight ranking feature at S9, never a gate: a
# PREDICTED_ONLY protein is not excluded from screening, and a CURATED_FAMILY protein is
# arguably LESS novel, since Pfam already recognised the family. Using it as a filter
# would invert its meaning.
_EVIDENCE_RUNGS = [
    # A curator built and named a family for it. Strongest evidence of reality.
    ("CURATED_FAMILY", re.compile(r"\bDUF\d*\b|\bUPF\d+", re.I)),
    # There used to be a MULTISPECIES rung here, read from NCBI's "MULTISPECIES:" title
    # prefix. It is removed: ClusteredNR titles are one representative's, and carry the
    # prefix in 1.8% of titles against 19.0% in full nr, so the rung would have measured
    # the database rather than the protein. Occurrence across hosts is measured by the
    # pipeline itself, from the plasmids a family occurs on (recurrence, rarity).
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

REQUIRED_THRESHOLDS = ("narrow_at", "min_explained", "min_coverage", "full_at", "partial_at")


def check_thresholds(cfg):
    """Validate relationships between thresholds that a JSON schema cannot express.

    The schema checks that each threshold is a number in [0, 1]. It cannot check that the
    numbers make sense TOGETHER, which is where the interesting failures live. Raises
    ValueError with a message naming the offending pair.

    narrow_at >= min_explained
        narrow_at withholds a protein from deeper tiers; min_explained decides what counts
        as explained in the output. If narrowing were the stricter of the two, the pipeline
        would stop searching proteins it then reports as unexplained - and, because the
        deeper tiers were never run on them, the contradiction could not be investigated.

    full_at >= partial_at
        Otherwise the completeness bands overlap and FULL becomes unreachable.
    """
    missing = [k for k in REQUIRED_THRESHOLDS if k not in cfg]
    if missing:
        raise ValueError(f"cascade config is missing threshold(s): {', '.join(missing)}")
    if cfg["narrow_at"] < cfg["min_explained"]:
        raise ValueError(
            f"narrow_at ({cfg['narrow_at']}) < min_explained ({cfg['min_explained']}): the "
            "cascade would stop searching proteins it then reports as unexplained, and the "
            "counterfactual could not be recovered")
    if cfg["full_at"] < cfg["partial_at"]:
        raise ValueError(
            f"full_at ({cfg['full_at']}) < partial_at ({cfg['partial_at']}): completeness "
            "bands overlap and FULL is unreachable")


# Tolerance on the declared -Z relative to the actual analysis-set size. -Z exists to fix
# the reference so an E-value means the same thing on every tier; it does not have to
# equal the input to the last sequence. 2% is well inside the noise of an E-value while
# being far tighter than any change a re-run of S1 would produce.
HMMER_Z_TOLERANCE = 0.02


def check_hmmer_z(declared, actual, tolerance=HMMER_Z_TOLERANCE):
    """Refuse a declared -Z that no longer describes the data it was computed from.

    hmmer_z is pinned in config so that E-values are comparable across tiers and across
    runs. But it is DERIVED from the analysis set: it is the number of unique protein
    sequences. Anything that changes the ORF set changes it - and S1 circular-origin repair
    changes the ORF set substantially, by reconstructing ~160,000 genes that the
    linearisation had broken in two.

    A stale value would silently rescale every E-value in the run, in exactly the way -Z
    was introduced to prevent. Raising here, with the correct number in the message, makes
    the config edit a deliberate one-line act rather than a thing to remember.
    """
    if not actual:
        raise ValueError("cannot check hmmer_z against an empty protein set")
    drift = abs(declared - actual) / actual
    if drift > tolerance:
        raise ValueError(
            f"hmmer_z in config/cascade.yaml is {declared}, but the analysis set holds "
            f"{actual} unique proteins ({drift:.1%} drift, tolerance {tolerance:.0%}). "
            f"Set hmmer_z: {actual} and re-run - a stale -Z rescales every E-value in the "
            "run, which is the failure -Z exists to prevent.")
