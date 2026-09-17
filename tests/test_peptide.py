"""Peptide physicochemistry, used to assign screening strata.

Membrane and secretion signals come from hydrophobicity windows and net charge.
Licence-restricted topology predictors are deliberately not part of this pipeline: a stage
that cannot be installed from envs/ is a stage that cannot be reproduced.

These are coarse, and are used to assign strata and to rank, never to exclude. Every output
row records `topology_method` so a reader knows what produced the call.
"""
import pytest
from plasmidann.peptide import net_charge, gravy, max_hydrophobic_window, is_cationic_amphipathic


def test_a_polylysine_peptide_is_strongly_cationic():
    assert net_charge("KKKKKKKKKK") == pytest.approx(10.0)


def test_a_polyglutamate_peptide_is_strongly_anionic():
    assert net_charge("EEEEEEEEEE") == pytest.approx(-10.0)


def test_histidine_counts_as_partially_charged_at_physiological_ph():
    """His has a pKa near 6, so at pH 7.4 it is only fractionally protonated. Counting it
    as a full positive charge would misclassify His-rich proteins as antimicrobial."""
    assert 0 < net_charge("HHHHHHHHHH") < 10


def test_a_hydrophobic_stretch_is_detected():
    """A run of 19 hydrophobic residues is the signature of a transmembrane helix."""
    assert max_hydrophobic_window("KKKK" + "LLLLLLLLLLLLLLLLLLL" + "KKKK", window=19) > 2.0


def test_a_soluble_protein_has_no_hydrophobic_stretch():
    assert max_hydrophobic_window("KEKEKEKEKEKEKEKEKEKEKEKEKE", window=19) < 0.0


def test_gravy_is_the_mean_hydropathy():
    assert gravy("IIII") == pytest.approx(4.5)
    assert gravy("RRRR") == pytest.approx(-4.5)


def test_an_empty_sequence_does_not_divide_by_zero():
    assert gravy("") == 0.0
    assert net_charge("") == 0.0
    assert max_hydrophobic_window("", window=19) == 0.0


def test_a_short_cationic_amphipathic_peptide_is_flagged():
    """The NOVOQR9B class: a 36-residue peptide beside resistance genes that inhibited
    every Gram-positive species tested. Short, cationic, amphipathic."""
    assert is_cationic_amphipathic("KWKLFKKIGAVLKVLTTGLPALIS") is True


def test_a_long_acidic_protein_is_not_flagged():
    assert is_cationic_amphipathic("E" * 200) is False


def test_a_short_but_anionic_peptide_is_not_flagged():
    assert is_cationic_amphipathic("EEEEDDDDEEEEDDDD") is False
