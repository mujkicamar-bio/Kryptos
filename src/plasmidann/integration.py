"""Stage 15: which evidence dimensions were measured for a dark family. No score, no rank.

    ORF_QC                     artefact screen
    SEQUENCE_HOMOLOGY          Pfam, Swiss-Prot, nr, pharokka
    ORTHOLOGY                  eggNOG
    GENOMIC_CONTEXT            context rates, synteny, dark family co-occurrence
    EVOLUTIONARY_CONSERVATION  dN/dS, RNAcode
    DISTRIBUTION               breadth over independent lineages
    STRUCTURAL_RELATIONSHIP    Foldseek
    PROTEIN_PROPERTIES         protein length

Several databases feed one dimension because they share evolutionary information: counting
their hits separately would make a well-studied protein look better supported. A dimension
counts when it was measured, not when the measurement was positive. The count of database
hits (supporting_observations_count) is reported beside it, labelled as not independent.
"""

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

# Which databases feed which dimension.
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

    # The context rates (0 is a measurement), a synteny status or a co-occurrence status.
    if (record.get("cons_annotated_neighbour") not in (None, "")
            or record.get("synteny_status") or record.get("cooccurrence_status")):
        present.append("GENOMIC_CONTEXT")

    # Status rather than value: a status without a value, such as TOO_SHORT, is a
    # measurement.
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
    """The number of informative database hits: how much was found, not how strong the
    evidence is, because the hits are not independent."""
    try:
        return int(record.get("n_informative_hits") or 0)
    except (TypeError, ValueError):
        return 0


def evidence_summary(record):
    """The dimensions present, their count and the observation count. No score, no rank."""
    present = dimensions_present(record)
    return {
        "evidence_dimensions_present": ",".join(present),
        "evidence_dimension_count": len(present),
        "supporting_observations_count": supporting_observations(record),
        # Named so that nobody reading the table can mistake the count for a score, and so
        # that a column selected into a downstream ranking carries the warning with it.
        "supporting_observations_are_not_independent": 1,
    }

