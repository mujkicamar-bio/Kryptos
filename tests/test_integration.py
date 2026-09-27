"""Evidence dimensions counted without a score (plasmidann.integration), on the keys the
family row of annotation_report carries."""
from plasmidann import integration

# A searched dark family with every stage measured, as annotation_report assembles it.
MEASURED = {"artefact_screened": 1, "cascade_searched": 1, "eggnog_searched": 1,
            "cons_annotated_neighbour": "0.0", "synteny_status": "TOO_FEW_LINEAGES",
            "cooccurrence_status": "TOO_FEW_LINEAGES", "dnds_status": "TOO_SHORT",
            "rnacode_status": "TOO_FEW_MEMBERS", "independent_cluster_status": "SUCCESS",
            "structure_searched": 1, "representative_length_aa": 120,
            "n_informative_hits": 0}


def test_every_dimension_is_present_when_every_stage_measured_it():
    """A search that found nothing, a rate of 0 and a status without a value are all
    measurements, which is what a dark family is made of."""
    assert integration.dimensions_present(MEASURED) == list(integration.DIMENSIONS)


def test_a_stage_that_did_not_measure_leaves_its_dimension_out():
    for key, dim in (("artefact_screened", "ORF_QC"),
                     ("cascade_searched", "SEQUENCE_HOMOLOGY"),
                     ("eggnog_searched", "ORTHOLOGY"),
                     ("structure_searched", "STRUCTURAL_RELATIONSHIP"),
                     ("representative_length_aa", "PROTEIN_PROPERTIES")):
        record = {**MEASURED, key: 0}
        assert dim not in integration.dimensions_present(record), key
    assert "DISTRIBUTION" not in integration.dimensions_present(
        {**MEASURED, "independent_cluster_status": "TOO_FEW_MEMBERS"})


def test_too_few_lineages_is_not_a_context_measurement():
    """TOO_FEW_LINEAGES says no synteny or co-occurrence test was made."""
    no_context = {**MEASURED, "cons_annotated_neighbour": ""}
    assert "GENOMIC_CONTEXT" not in integration.dimensions_present(no_context)
    for key in ("synteny_status", "cooccurrence_status"):
        assert "GENOMIC_CONTEXT" in integration.dimensions_present(
            {**no_context, key: "SUCCESS"}), key


def test_the_summary_counts_and_labels_the_observations_and_holds_no_score():
    summary = integration.evidence_summary({**MEASURED, "n_informative_hits": 7})

    assert summary["evidence_dimension_count"] == len(integration.DIMENSIONS)
    assert summary["supporting_observations_count"] == 7
    assert summary["supporting_observations_are_not_independent"] == 1
    for forbidden in ("score", "rank", "novelty", "candidate", "priority"):
        assert not any(forbidden in key.lower() for key in summary), forbidden
