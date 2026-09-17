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


def _evalue_of(hit):
    """Sort key for hit significance. A hit without a parsable E-value sorts last.

    Returning infinity rather than raising keeps bit-score-only tiers (T1) sortable: when
    no hit in a set has an E-value they all tie, and the caller's secondary key - cascade
    order, which is authority order - decides.
    """
    try:
        return float(hit.get("evalue"))
    except (TypeError, ValueError):
        return float("inf")


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

    The label comes from the most significant hit, not the widest one, because a longer
    alignment is not a better identification. Measured: 15.4% of labels change, and the
    example that settled it was ABC_membrane at E=1e-23 being chosen over Peptidase_C39
    at E=6.5e-40 purely because it aligned further.

    functional_class is deliberately not a boolean. UNCHARACTERIZED_HOMOLOG - a protein
    whose only homologs are themselves unnamed - is the class worth screening: certainly
    real, certainly unknown.

    `tier_order` is passed in rather than read from a constant so that a four-tier cascade
    reports depths on a four-tier scale. v1 indexed into a hardcoded eight-tier list.
    """
    if not hits:
        return {"annot_tier": None, "annot_label": None, "functional_class": "NONE",
                "homology_depth": None, "annot_qcov": None, "annot_tcov": None,
                "annot_evalue": None, "n_informative_hits": 0}

    informative = [h for h in hits if is_informative(h.get("label"))]
    if informative:
        # Strongest evidence wins; cascade order breaks ties, since it is authority order.
        best = min(informative, key=lambda h: (_evalue_of(h), tier_order.index(h["tier"])))
        cls = "FUNCTIONAL" if explained >= min_coverage else "DOMAIN_ONLY"
        depth_from = informative
    else:
        # Nothing named it anywhere. Show the most authoritative record of having seen it.
        best = min(hits, key=lambda h: (tier_order.index(h["tier"]), _evalue_of(h)))
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
    }


# ---------------------------------------------------------------------------------
# Narrowing: what each tier hands the next
# ---------------------------------------------------------------------------------

def narrow(all_ids, hit_ids):
    """Ids a tier could not resolve. Refuses hits for ids that were never queried.

    The refusal is deliberate. A hit id absent from the query set means the search output
    and the ORF index disagree, which is a class of bug that would otherwise corrupt
    counts silently for the rest of the run.
    """
    unknown = set(hit_ids) - set(all_ids)
    if unknown:
        raise ValueError(f"hits for {len(unknown)} id(s) not in the query set, "
                         f"e.g. {sorted(unknown)[:3]}")
    return [i for i in all_ids if i not in hit_ids]


def narrow_by_explained(all_ids, explained, threshold):
    """Ids still worth searching: those not yet explained past `threshold`.

    Narrowing on explained fraction rather than on "got any hit" is what stops a 15%
    domain match from terminating the search over the other 85% of a protein.

    IMPORTANT - which threshold belongs here. The caller must pass `narrow_at`, NOT
    `min_explained`. They are two different numbers doing two different jobs:

      narrow_at (0.9)      how finished a protein must be before we stop searching it
      min_explained (0.5)  how explained a protein must be before we REPORT it as such

    v1 used one number for both. Because a protein removed at T2 has no T5 result, that
    made the threshold unsweepable: the counterfactual does not exist, and answering
    "what if it had been 0.7?" would need a full re-run per value. Measured, 0.3 -> 0.8
    moved the deep-tier set by 76%, and 0.5 sat exactly at the 25th percentile of the
    observed distribution - the densest possible place to put a hard cut.

    Keeping narrow_at permissive costs deep-tier compute and buys back the ability to
    check our own threshold. See also check_thresholds, which enforces narrow_at >=
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


def uninformative_hits(hits):
    """Hits that record only that someone else has seen this protein.

    Kept rather than discarded: a protein called 'hypothetical' by several independent
    databases is real, widespread and genuinely uncharacterised - which is evidence for it
    being a screening target, not against.
    """
    return [h for h in hits if not is_informative(h.get("label"))]


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
    # The identical sequence occurs in more than one species.
    #
    # NOT anchored with ^. DIAMOND emits `stitle`, which looks like
    #   WP_000123.1 MULTISPECIES: hypothetical protein [Enterobacteriaceae]
    # so the accession comes first and an anchored pattern never fires on real output.
    # In v1 this rung was unreachable for the entire run.
    ("MULTISPECIES", re.compile(r"\bMULTISPECIES\s*:", re.I)),
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
# the reference so an E-value means the same thing on every shard; it does not have to
# equal the input to the last sequence. 2% is well inside the noise of an E-value while
# being far tighter than any change a re-run of S1 would produce.
HMMER_Z_TOLERANCE = 0.02


def check_hmmer_z(declared, actual, tolerance=HMMER_Z_TOLERANCE):
    """Refuse a declared -Z that no longer describes the data it was computed from.

    hmmer_z is pinned in config so that E-values are comparable across shards and across
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
