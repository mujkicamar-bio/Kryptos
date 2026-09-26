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
    SEQUENCE_HOMOLOGY          Pfam, Swiss-Prot, nr, pharokka
    ORTHOLOGY                  eggNOG
    GENOMIC_CONTEXT            neighbours, defence systems, integrons, synteny,
                               dark family co-occurrence
    EVOLUTIONARY_CONSERVATION  dN/dS, RNAcode
    DISTRIBUTION               breadth over independent lineages
    STRUCTURAL_RELATIONSHIP    Foldseek
    PROTEIN_PROPERTIES         length, charge, transmembrane helices

DEPENDENCE IS THE POINT OF THE DIMENSIONS (section 56.2)

"The objective is not to manufacture 'independent evidence' but to preserve distinct
measurements and avoid double-counting them."

Pfam, Swiss-Prot, nr and pharokka are four databases and ONE dimension. They share
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

# Which databases feed which dimension. Several databases to ONE dimension is the whole
# point (section 56.2): they are not independent, so they must not be counted separately.
_DIMENSION_SOURCES = {
    "SEQUENCE_HOMOLOGY": ("pfam", "swissprot", "nr", "pharokka"),
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

    # The S8c context rates (0 is a measurement), a synteny status or an S8g co-occurrence
    # status.
    if (record.get("cons_annotated_neighbour") not in (None, "")
            or record.get("synteny_status") or record.get("cooccurrence_status")):
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

