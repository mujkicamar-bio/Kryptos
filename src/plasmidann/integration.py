"""Stage 15: evidence integration without a ranking (spec section 56).

WHAT "INTEGRATION" MEANS HERE

Section 56, first line: "Evidence integration combines results without producing an
experimental ranking." So this module counts and names what is known about a protein. It
produces no score, and section 2.3 says why: a composite destroys the information a
screening decision needs. Two proteins with the same total can be entirely different bets -
one with overwhelming evidence that it is a real protein and no idea what it does, the other
with a sharp hypothesis resting on almost nothing - and those demand different experiments.

THE EIGHT DIMENSIONS (section 56.1)

    ORF_QC                     is it a protein at all
    SEQUENCE_HOMOLOGY          Pfam, Swiss-Prot, nr, PHROGs
    ORTHOLOGY                  eggNOG
    GENOMIC_CONTEXT            neighbours, defence systems, integrons, synteny
    EVOLUTIONARY_CONSERVATION  dN/dS, RNAcode
    DISTRIBUTION               breadth over independent lineages
    STRUCTURAL_RELATIONSHIP    Foldseek
    PROTEIN_PROPERTIES         length, charge, transmembrane helices

DEPENDENCE IS THE POINT OF THE DIMENSIONS (section 56.2)

"The objective is not to manufacture 'independent evidence' but to preserve distinct
measurements and avoid double-counting them."

Pfam, Swiss-Prot, nr and PHROGs are four databases and ONE dimension. They share
evolutionary information: a protein found in Swiss-Prot is usually in nr, and a Pfam domain
is built from the same alignments. Counting four hits as four lines of evidence would make
a well-studied protein look four times better supported than an equally well-supported one
that happened to be in fewer databases. So the count is of DIMENSIONS PRESENT, never of
observations across databases - and `supporting_observations_count` is reported beside it as
a separate, explicitly non-independent number.

Foldseek is deliberately its own dimension despite also being homology. Structure reaches
further back than sequence, so a structural match where sequence found nothing is a genuinely
different measurement - but section 2.8 warns they are still related, which is why they are
two named dimensions and not two points on a scale.

THE HYPOTHESIS LAYER (section 57)

A hypothesis is a named association with traceable support, and it never changes what the
protein IS. Section 57 is explicit: "The protein remains DARK unless direct sequence/domain
evidence supports a specific function." A dark ORF beside a CBASS system is
defence_associated; it is not a defence protein, and section 2.7 makes that distinction a
design principle.
"""

# Section 56.1, verbatim. A dimension not in this set is one nothing downstream can read.
DIMENSIONS = (
    "ORF_QC",
    "SEQUENCE_HOMOLOGY",
    "ORTHOLOGY",
    "GENOMIC_CONTEXT",
    "EVOLUTIONARY_CONSERVATION",
    "DISTRIBUTION",
    "STRUCTURAL_RELATIONSHIP",
    "PROTEIN_PROPERTIES",
)

# Section 57. `unknown` is a real value, not a failure to assign: a dark protein with no
# contextual association is the most common case and the one the project is built around.
HYPOTHESES = (
    "defence_associated",
    "mobilization_associated",
    "toxin_antitoxin_associated",
    "replication_associated",
    "partition_associated",
    "regulatory_associated",
    "membrane_associated",
    "metabolic_associated",
    "structural_associated",
    "unknown",
)

# Which databases feed which dimension. Several databases to ONE dimension is the whole
# point (section 56.2): they are not independent, so they must not be counted separately.
_DIMENSION_SOURCES = {
    "SEQUENCE_HOMOLOGY": ("pfam", "swissprot", "nr", "phrogs"),
    "ORTHOLOGY": ("eggnog",),
    "STRUCTURAL_RELATIONSHIP": ("foldseek",),
}


def dimensions_present(record):
    """Which of the eight dimensions have a measurement for this protein.

    `record` is a flat dict of the columns the pipeline produced. A dimension counts as
    present when it was MEASURED, not when it was positive: a search that ran and found
    nothing is a measurement, and is exactly the measurement a dark protein is made of.
    """
    present = []

    # QC always runs, so this dimension is present whenever a flag was recorded at all.
    if record.get("artefact_flag") not in (None, ""):
        present.append("ORF_QC")

    if any(record.get(f"{source}_searched") for source in
           _DIMENSION_SOURCES["SEQUENCE_HOMOLOGY"]) or record.get("annot_tier"):
        present.append("SEQUENCE_HOMOLOGY")

    if record.get("eggnog_searched") or record.get("cog_category"):
        present.append("ORTHOLOGY")

    if record.get("context_status") == "SUCCESS" or record.get("top_hypothesis"):
        present.append("GENOMIC_CONTEXT")

    # Status rather than value: NO_DIVERGENCE and SATURATED are measurements, and a family
    # of identical sequences is what strong conservation looks like.
    if record.get("dnds_status") or record.get("rnacode_status"):
        present.append("EVOLUTIONARY_CONSERVATION")

    if record.get("independent_cluster_status") == "SUCCESS":
        present.append("DISTRIBUTION")

    if record.get("structure_status") or record.get("structural_match"):
        present.append("STRUCTURAL_RELATIONSHIP")

    if record.get("protein_length"):
        present.append("PROTEIN_PROPERTIES")

    return present


def supporting_observations(record):
    """How many individual observations support this protein.

    Reported because section 56.3 asks for it, and labelled everywhere it appears as NOT
    independent. It counts database hits, so a well-studied protein scores higher than an
    equally well-supported one that happens to be in fewer databases. It is a description
    of how much was found, never a measure of how strong the evidence is - the dimension
    count is the number that avoids double-counting.
    """
    try:
        return int(record.get("n_informative_hits") or 0)
    except (TypeError, ValueError):
        return 0


def evidence_summary(record):
    """The section 56.3 descriptive summary: dimensions present and observation count.

    No score, no rank, no total. Section 56: integration "combines results without producing
    an experimental ranking".
    """
    present = dimensions_present(record)
    return {
        "evidence_dimensions_present": ",".join(present),
        "evidence_dimension_count": len(present),
        "supporting_observations_count": supporting_observations(record),
        # Named so that nobody reading the table can mistake the count for a score, and so
        # that a column selected into a downstream ranking carries the warning with it.
        "supporting_observations_are_not_independent": 1,
    }


def functional_hypothesis(record):
    """A named association with traceable support (section 57).

    Returns (hypothesis, support) where `support` names the evidence behind it. The protein
    remains DARK: section 57 states that a hypothesis does not change what the protein IS
    unless direct sequence or domain evidence supports a specific function, and section 2.7
    makes it a design principle - a dark ORF repeatedly beside a defence system is
    defence_associated, never defence_protein.

    Only ONE hypothesis is returned, the best-supported, because a hypothesis is what
    someone takes to the bench and a list of ten is not a testable prediction. Every
    association that fired is still in the context table, so nothing is lost.
    """
    context = (record.get("top_hypothesis") or "").lower()
    if not context:
        return "unknown", ""

    support = []
    for field in ("top_conservation", "top_q_value", "synteny_conservation",
                  "context_recurrence"):
        if record.get(field) not in (None, ""):
            support.append(f"{field}={record[field]}")
    detail = "; ".join(support)

    # The mapping is intentionally narrow and reads the category names the label vocabulary
    # produces. An association this does not recognise stays `unknown` WITH its context
    # recorded, rather than being forced into the nearest hypothesis - a wrong hypothesis
    # sends someone to the bench to test the wrong thing.
    for token, hypothesis in (
            ("defence", "defence_associated"),
            ("defense", "defence_associated"),
            ("mobilis", "mobilization_associated"),
            ("mobiliz", "mobilization_associated"),
            ("conjug", "mobilization_associated"),
            ("toxin", "toxin_antitoxin_associated"),
            ("antitoxin", "toxin_antitoxin_associated"),
            ("replicat", "replication_associated"),
            ("partition", "partition_associated"),
    ):
        if token in context:
            return hypothesis, detail

    return "unknown", f"context={record['top_hypothesis']}" + (f"; {detail}" if detail else "")
