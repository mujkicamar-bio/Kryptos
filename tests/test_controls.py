"""The positive control for the annotation cascade (SC2).

Moved out of tests/test_backbone.py when the curated backbone family list was deleted. The
control is unrelated to that list and always was: it draws on reviewed Swiss-Prot proteins
spiked into the query set, which is what makes it able to detect a protein the cascade
missed.
"""
import pytest

from plasmidann.controls import control_recall


def test_control_recall_is_the_fraction_correctly_called_functional():
    """ECLIPSE recovered 99.2-100% of 246 virulence, 42 AMR and 75 essential genes as
    annotated. Anything known that comes out dark is a recall failure in the cascade."""
    rows = [{"functional_class": "FUNCTIONAL"}] * 99 + [{"functional_class": "NONE"}]

    assert control_recall(rows) == 0.99


def test_a_control_set_where_everything_is_dark_scores_zero():
    rows = [{"functional_class": "UNCHARACTERIZED_HOMOLOG"}] * 10

    assert control_recall(rows) == 0.0


def test_an_empty_control_set_is_an_error_not_a_pass():
    """Silently passing on an empty control is how a gate stops being a gate."""
    with pytest.raises(ValueError):
        control_recall([])
