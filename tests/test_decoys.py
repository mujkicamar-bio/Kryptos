"""Negative controls: decoys that must not be annotated (spec section 58.2).

Positive controls ask whether the pipeline can still recover known biology. These ask the
opposite, and it is the question this project cannot skip: the deliverable is the DARK set,
which is defined by what nothing could name. If the cascade can be induced to name a
sequence that is not a protein, the dark set is not what it claims to be.
"""
import random

from plasmidann import decoys


def test_a_shuffled_decoy_preserves_length_and_composition_exactly():
    """The decoy must differ from the real protein in ORDER alone, so that any hit to it is
    a hit to composition - which is not evidence of function."""
    protein = "MKVLATTLLGAAFAASSALAQKKWLVR"

    decoy = decoys.shuffled_decoy(protein, random.Random(1))

    assert len(decoy) == len(protein)
    assert sorted(decoy) == sorted(protein)


def test_a_shuffled_decoy_is_reproducible_from_the_seed():
    """A control set that changed between runs could not be compared between runs."""
    protein = "MKVLATTLLGAAFAASSALAQKKWLVR"

    assert (decoys.shuffled_decoy(protein, random.Random(7))
            == decoys.shuffled_decoy(protein, random.Random(7)))


def test_reverse_complement_is_its_own_inverse():
    dna = "ATGGCATTTCGATAG"

    assert decoys.reverse_complement(decoys.reverse_complement(dna)) == dna


def test_translation_stops_at_the_first_stop_codon():
    """A sequence carrying an internal '*' is not a protein the searches can handle: every
    tool would reject it or read the stop as an unknown residue. Truncating gives the
    longest real open frame, which is what a shadow ORF actually looks like."""
    assert decoys.translate("ATGGCATAAGGGCCC") == "MA"


def test_a_reverse_complement_decoy_is_not_the_original_protein():
    """The whole construction: read a real CDS on the opposite strand and you get a
    sequence with realistic codon statistics that is not a protein. This is the shadow-ORF
    artifact manufactured on purpose."""
    cds = "ATGGCATTTCGAGGCATTAAACGT" * 8

    forward = decoys.translate(cds)
    decoy = decoys.reverse_complement_decoy(cds)

    assert decoy
    assert decoy != forward


def test_build_decoys_produces_both_classes_in_the_requested_numbers():
    rng = random.Random(3)
    cds = [(f"c{i}", "ATGGCATTTCGAGGCATTAAACGTTTGCAAGAC" * 8) for i in range(50)]

    built = decoys.build_decoys(cds, n_shuffled=5, n_reverse_complement=5, rng=rng)

    classes = {}
    for _, _, cls in built:
        classes[cls] = classes.get(cls, 0) + 1
    assert classes.get("shuffled") == 5
    assert classes.get("reverse_complement") == 5


def test_every_decoy_id_carries_the_prefix_and_is_unique():
    """Every later stage recognises a decoy by its id alone, without a lookup table - the
    same contract the positive controls use. A duplicate id would silently merge two
    decoys in the dereplicated set."""
    rng = random.Random(4)
    cds = [(f"c{i}", "ATGGCATTTCGAGGCATTAAACGTTTGCAAGAC" * 8) for i in range(50)]

    built = decoys.build_decoys(cds, n_shuffled=5, n_reverse_complement=5, rng=rng)
    ids = [i for i, _, _ in built]

    assert all(i.startswith(decoys.DECOY_PREFIX) for i in ids)
    assert len(set(ids)) == len(ids)


def test_a_decoy_too_short_to_test_anything_is_not_built():
    """A 12-residue decoy passes the negative control for a reason that has nothing to do
    with the pipeline behaving correctly: nothing would name it either way."""
    rng = random.Random(5)
    cds = [("short", "ATGGCATTTCGA")] * 20

    assert decoys.build_decoys(cds, 5, 5, rng, min_length=50) == []


def test_a_cds_with_ambiguity_codes_is_not_used():
    """An N in the CDS translates to X, and a run of X is neither a real protein nor a
    decoy with realistic composition - it is a third thing that tests nothing."""
    rng = random.Random(6)
    cds = [("amb", "ATGNNNTTTCGAGGCATTAAACGTTTGCAAGAC" * 8)] * 20

    assert decoys.build_decoys(cds, 5, 5, rng) == []
