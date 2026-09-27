"""Stage 15: evidence dimensions counted without a score (plasmidann.integration)."""
from plasmidann import integration


def test_four_sequence_databases_are_one_dimension_not_four():
    """Pfam, Swiss-Prot, nr and pharokka share evolutionary information, so their hits
    fill one dimension, not four."""
    record = {"pfam_searched": 1, "swissprot_searched": 1, "nr_searched": 1,
              "pharokka_searched": 1, "annot_tier": "T1"}

    present = integration.dimensions_present(record)

    assert present.count("SEQUENCE_HOMOLOGY") == 1
    assert len(present) == 1, f"four databases produced {len(present)} dimensions"


def test_structure_is_a_separate_dimension_from_sequence():
    """Structure reaches further back than sequence, so a structural match is its own
    dimension."""
    record = {"annot_tier": "T1", "structural_match": "1abc_A"}

    present = integration.dimensions_present(record)

    assert "SEQUENCE_HOMOLOGY" in present
    assert "STRUCTURAL_RELATIONSHIP" in present


def test_a_dimension_is_present_when_measured_not_when_positive():
    """A search that ran and found nothing is a measurement, and it is exactly the
    measurement a dark protein is made of. Requiring a positive result would make the dark
    set look evidence-free by construction."""
    record = {"dnds_status": "NO_DIVERGENCE", "rnacode_status": "TOO_FEW_MEMBERS"}

    assert "EVOLUTIONARY_CONSERVATION" in integration.dimensions_present(record)


def test_the_summary_contains_no_score_or_rank():
    summary = integration.evidence_summary({"annot_tier": "T1", "protein_length": 120})

    for forbidden in ("score", "rank", "novelty", "candidate", "priority"):
        assert not any(forbidden in key.lower() for key in summary), (
            f"the evidence summary contains a {forbidden} field")


def test_the_observation_count_is_labelled_non_independent():
    """It counts database hits, so it rewards being well studied; the flag travels with
    the column."""
    summary = integration.evidence_summary({"n_informative_hits": 7})

    assert summary["supporting_observations_count"] == 7
    assert summary["supporting_observations_are_not_independent"] == 1


def test_every_emitted_dimension_is_declared():
    record = {"artefact_flag": 0, "annot_tier": "T1", "cog_category": "L",
              "cons_annotated_neighbour": 0.5, "dnds_status": "SUCCESS",
              "independent_cluster_status": "SUCCESS", "structural_match": "x",
              "protein_length": 100}

    assert set(integration.dimensions_present(record)) <= set(integration.DIMENSIONS)


def test_cooccurrence_is_genomic_context_and_not_a_new_dimension():
    """S8g asks which dark families share a plasmid: that is genomic context, so it
    fills GENOMIC_CONTEXT once, alone or beside synteny, and adds no ninth dimension."""
    alone = integration.dimensions_present({"cooccurrence_status": "TOO_FEW_LINEAGES"})
    both = integration.dimensions_present({"cooccurrence_status": "SUCCESS",
                                           "synteny_status": "SUCCESS"})

    assert alone == ["GENOMIC_CONTEXT"]
    assert both == ["GENOMIC_CONTEXT"]
