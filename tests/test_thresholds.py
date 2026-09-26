"""Threshold coherence: relationships a JSON schema cannot express.

The schema checks each threshold is a number in [0, 1]. It cannot check that the numbers
make sense together, which is where the interesting failures live.
"""
import pytest

from plasmidann.cascade import check_thresholds, completeness

OK = {"narrow_at": 0.7, "min_explained": 0.5, "min_coverage": 0.5,
      "full_at": 0.8, "partial_at": 0.5}


def test_a_coherent_threshold_block_is_accepted():
    check_thresholds(OK)


def test_narrowing_may_not_be_stricter_than_reporting():
    """If narrow_at < min_explained the cascade stops searching proteins it then reports
    as unexplained - and since the deeper tiers never ran, the contradiction cannot be
    investigated afterwards."""
    with pytest.raises(ValueError, match="narrow_at"):
        check_thresholds({**OK, "narrow_at": 0.4, "min_explained": 0.5})


def test_narrowing_equal_to_reporting_is_allowed():
    """The boundary case is coherent, if unadventurous: it is what v1 did."""
    check_thresholds({**OK, "narrow_at": 0.5, "min_explained": 0.5})


def test_completeness_bands_may_not_overlap():
    with pytest.raises(ValueError, match="full_at"):
        check_thresholds({**OK, "full_at": 0.4, "partial_at": 0.5})


def test_a_missing_threshold_is_named():
    """Silence here would mean a stage running on a default nobody declared."""
    cfg = {k: v for k, v in OK.items() if k != "min_coverage"}
    with pytest.raises(ValueError, match="min_coverage"):
        check_thresholds(cfg)


# --- completeness now takes its bands as arguments, not module constants -------------

def test_completeness_bands_come_from_config_not_from_constants():
    """FULL_AT and PARTIAL_AT were two of the three thresholds that escaped config in v1,
    so they could never be swept or recorded."""
    assert completeness(0.95, full_at=0.8, partial_at=0.5) == "FULL"
    assert completeness(0.60, full_at=0.8, partial_at=0.5) == "PARTIAL"
    assert completeness(0.15, full_at=0.8, partial_at=0.5) == "FRAGMENT"
    assert completeness(0.0, full_at=0.8, partial_at=0.5) == "NONE"


def test_moving_the_bands_moves_the_answer():
    """The point of making them arguments: the same fraction bands differently under a
    different declared policy, and that policy is now visible and sweepable."""
    assert completeness(0.75, full_at=0.8, partial_at=0.5) == "PARTIAL"
    assert completeness(0.75, full_at=0.7, partial_at=0.5) == "FULL"


# --- hmmer_z must not drift away from the data it describes --------------------------

def test_a_declared_hmmer_z_matching_the_data_is_accepted():
    from plasmidann.cascade import check_hmmer_z

    check_hmmer_z(declared=3_497_616, actual=3_497_616)


def test_small_drift_is_tolerated():
    """-Z fixes the reference so an E-value means the same thing on every shard. It need
    not equal the input size to the last sequence; it must not be materially wrong."""
    from plasmidann.cascade import check_hmmer_z

    check_hmmer_z(declared=3_497_616, actual=3_500_000)


def test_material_drift_is_refused_and_names_the_correct_value():
    """S1 origin repair changes the ORF set, and therefore the unique-protein count. A
    stale hmmer_z would silently rescale every E-value in the run."""
    from plasmidann.cascade import check_hmmer_z

    with pytest.raises(ValueError, match="3400000"):
        check_hmmer_z(declared=3_497_616, actual=3_400_000)
