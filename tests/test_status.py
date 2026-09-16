"""Spec §2.9 and §7.2: four different reasons for a missing number must not collapse into
one value. A null dN/dS because the family had three members and a null because the
alignment saturated are different scientific statements."""
from darkorf import status


def test_the_vocabulary_is_centralized_and_complete():
    for expected in [
        "NOT_RUN", "NO_HIT", "TOO_FEW_MEMBERS", "NO_DIVERGENCE",
        "SATURATED", "NO_OUTPUT", "FAILED", "NOT_APPLICABLE", "SUCCESS",
    ]:
        assert expected in status.ALL
        assert getattr(status, expected) == expected


def test_na_is_not_a_status():
    """NA is for genuinely unavailable text only. Using it as a status is exactly the
    concealment spec §7.3 forbids."""
    assert "NA" not in status.ALL
