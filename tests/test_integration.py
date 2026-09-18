"""Stage 15: evidence integration without a ranking (spec sections 56 and 57).

Section 56, first line: "Evidence integration combines results without producing an
experimental ranking." Section 56.2 gives the harder half: Pfam, pharokka, Swiss-Prot, eggNOG,
nr, MMseqs2 and Foldseek "share evolutionary information to varying degrees", and the
objective "is not to manufacture 'independent evidence' but to preserve distinct
measurements and avoid double-counting them".
"""
from plasmidann import integration


def test_four_sequence_databases_are_one_dimension_not_four():
    """The central rule of section 56.2. Counting four database hits as four lines of
    evidence would make a well-studied protein look four times better supported than an
    equally well-supported one that happens to be in fewer databases."""
    record = {"pfam_searched": 1, "swissprot_searched": 1, "nr_searched": 1,
              "pharokka_searched": 1, "annot_tier": "T1"}

    present = integration.dimensions_present(record)

    assert present.count("SEQUENCE_HOMOLOGY") == 1
    assert len(present) == 1, f"four databases produced {len(present)} dimensions"


def test_structure_is_a_separate_dimension_from_sequence():
    """Section 2.8: structure reaches further back than sequence, so a structural match
    where sequence found nothing is a genuinely different measurement - while still being
    related, which is why they are two named dimensions rather than points on a scale."""
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
    """Section 56: integration combines results WITHOUT producing an experimental ranking,
    and section 2.3 forbids collapsing evidence into a single number."""
    summary = integration.evidence_summary({"annot_tier": "T1", "protein_length": 120})

    for forbidden in ("score", "rank", "novelty", "candidate", "priority"):
        assert not any(forbidden in key.lower() for key in summary), (
            f"the evidence summary contains a {forbidden} field")


def test_the_observation_count_is_labelled_non_independent():
    """It counts database hits, so it rewards being well studied. Reported because section
    56.3 asks for it, and flagged so a column selected into a downstream ranking carries
    the warning with it."""
    summary = integration.evidence_summary({"n_informative_hits": 7})

    assert summary["supporting_observations_count"] == 7
    assert summary["supporting_observations_are_not_independent"] == 1


def test_every_emitted_dimension_is_declared():
    record = {"artefact_flag": 0, "annot_tier": "T1", "cog_category": "L",
              "top_hypothesis": "defence", "dnds_status": "SUCCESS",
              "independent_cluster_status": "SUCCESS", "structural_match": "x",
              "protein_length": 100}

    assert set(integration.dimensions_present(record)) <= set(integration.DIMENSIONS)


def test_a_context_association_becomes_a_hypothesis_with_its_support():
    """Section 57: each hypothesis must have traceable support."""
    record = {"top_hypothesis": "defence", "top_conservation": 0.9,
              "synteny_conservation": 0.85}

    hypothesis, support = integration.functional_hypothesis(record)

    assert hypothesis == "defence_associated"
    assert "top_conservation=0.9" in support
    assert "synteny_conservation=0.85" in support


def test_a_hypothesis_is_an_association_never_an_identity():
    """Section 2.7 and section 57: a dark ORF repeatedly next to a defence system may be
    defence_associated; it must not be called a defence protein, and the protein remains
    DARK unless direct sequence or domain evidence says otherwise."""
    hypothesis, _ = integration.functional_hypothesis({"top_hypothesis": "defence"})

    assert hypothesis.endswith("_associated")
    assert hypothesis in integration.HYPOTHESES
    assert "protein" not in hypothesis


def test_an_unrecognised_association_stays_unknown_with_its_context_recorded():
    """Forcing it into the nearest hypothesis would send someone to the bench to test the
    wrong thing. The context is kept so the association is not lost."""
    hypothesis, support = integration.functional_hypothesis(
        {"top_hypothesis": "pfam_family:DUF1234"})

    assert hypothesis == "unknown"
    assert "DUF1234" in support


def test_no_context_gives_the_unknown_hypothesis():
    """The most common case, and the one the project is built around."""
    assert integration.functional_hypothesis({}) == ("unknown", "")
