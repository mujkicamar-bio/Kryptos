"""Threshold coherence: relationships a JSON schema cannot express.

The schema checks each threshold is a number in [0, 1]. It cannot check that the numbers
make sense together, which is where the interesting failures live.
"""
import pytest

from plasmidann.cascade import check_hmmer_z, check_thresholds, completeness

OK = {"narrow_at": 0.7, "min_coverage": 0.5, "full_at": 0.8, "partial_at": 0.5}


def test_a_coherent_threshold_block_is_accepted():
    check_thresholds(OK)


def test_narrowing_may_not_stop_below_the_functional_threshold():
    """If narrow_at < min_coverage a protein explained between the two stops being searched
    and is then reported DOMAIN_ONLY, although the deeper tiers might have explained it."""
    with pytest.raises(ValueError, match="min_coverage"):
        check_thresholds({**OK, "narrow_at": 0.6, "min_coverage": 0.8})


def test_narrowing_equal_to_the_functional_threshold_is_allowed():
    check_thresholds({**OK, "narrow_at": 0.5, "min_coverage": 0.5})


def test_completeness_bands_may_not_overlap():
    with pytest.raises(ValueError, match="full_at"):
        check_thresholds({**OK, "full_at": 0.4, "partial_at": 0.5})


def test_a_missing_threshold_is_named():
    """Silence here would mean a stage running on a default nobody declared."""
    cfg = {k: v for k, v in OK.items() if k != "min_coverage"}
    with pytest.raises(ValueError, match="min_coverage"):
        check_thresholds(cfg)


# --- completeness bands are arguments -------------------------------------------------

def test_moving_the_bands_moves_the_answer():
    """The same fraction bands differently under different declared bands."""
    assert completeness(0.75, full_at=0.8, partial_at=0.5) == "PARTIAL"
    assert completeness(0.75, full_at=0.7, partial_at=0.5) == "FULL"


# --- hmmer_z must not drift away from the data it describes --------------------------

def test_a_declared_hmmer_z_matching_the_data_is_accepted():
    check_hmmer_z(declared=3_497_616, actual=3_497_616)


def test_small_drift_is_tolerated():
    """-Z holds the reference constant so an E-value means the same thing on every tier. It
    need not equal the input size to the last sequence; it must not be materially wrong."""
    check_hmmer_z(declared=3_497_616, actual=3_500_000)


def test_material_drift_is_refused_and_names_the_correct_value():
    """A changed ORF set changes the unique-protein count; a stale hmmer_z would rescale
    every E-value in the run."""
    with pytest.raises(ValueError, match="3400000"):
        check_hmmer_z(declared=3_497_616, actual=3_400_000)
