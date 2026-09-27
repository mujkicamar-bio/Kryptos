"""Which evidence dimensions were measured for a dark family. No score, no rank.

    ORF_QC                     artefact screen (AntiFam, low complexity)
    SEQUENCE_HOMOLOGY          the annotation cascade (Pfam, Swiss-Prot, ClusteredNR, pharokka)
    ORTHOLOGY                  eggNOG
    GENOMIC_CONTEXT            context rates, synteny, dark family co-occurrence
    EVOLUTIONARY_CONSERVATION  dN/dS, RNAcode
    DISTRIBUTION               breadth over independent lineages
    STRUCTURAL_RELATIONSHIP    Foldseek
    PROTEIN_PROPERTIES         representative length

Several databases feed one dimension because they share evolutionary information: counting
their hits separately would make a well-studied protein look better supported. A dimension
counts when it was measured, not when the measurement was positive. The count of
informative database hits over the family's members (supporting_observations_count) is
reported beside it, labelled as not independent.

The record is the family row of annotation_report, which supplies the measured fields
artefact_screened, cascade_searched, eggnog_searched, structure_searched,
representative_length_aa and n_informative_hits beside the stage columns.
"""
from darkorf import status

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


def _context_status(value):
    # TOO_FEW_LINEAGES says no test could be made, so it is not a context measurement.
    return value not in (None, "", status.TOO_FEW_LINEAGES)


def dimensions_present(record):
    """Which of the eight dimensions have a measurement for this family.

    A search that ran and found nothing is a measurement, and is exactly the measurement a
    dark protein is made of.
    """
    present = []
    if record.get("artefact_screened"):
        present.append("ORF_QC")
    if record.get("cascade_searched"):
        present.append("SEQUENCE_HOMOLOGY")
    if record.get("eggnog_searched"):
        present.append("ORTHOLOGY")
    # A context rate of 0 is a measurement.
    if (record.get("cons_annotated_neighbour") not in (None, "")
            or _context_status(record.get("synteny_status"))
            or _context_status(record.get("cooccurrence_status"))):
        present.append("GENOMIC_CONTEXT")
    # Status rather than value: a status without a value, such as TOO_SHORT, is a
    # measurement.
    if record.get("dnds_status") or record.get("rnacode_status"):
        present.append("EVOLUTIONARY_CONSERVATION")
    if record.get("independent_cluster_status") == status.SUCCESS:
        present.append("DISTRIBUTION")
    if record.get("structure_searched"):
        present.append("STRUCTURAL_RELATIONSHIP")
    if record.get("representative_length_aa"):
        present.append("PROTEIN_PROPERTIES")
    return present


def evidence_summary(record):
    """The dimensions present, their count and the observation count. No score, no rank."""
    present = dimensions_present(record)
    return {
        "evidence_dimensions_present": ",".join(present),
        "evidence_dimension_count": len(present),
        # How much was found, not how strong the evidence is: the hits are not independent.
        "supporting_observations_count": int(record.get("n_informative_hits") or 0),
        # Named so that nobody reading the table can mistake the count for a score, and so
        # that a column selected into a downstream ranking carries the warning with it.
        "supporting_observations_are_not_independent": 1,
    }
