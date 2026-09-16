"""Origin handling for circular plasmids (spec §8.4).

The invariant: a circular sequence has no privileged starting point, so calling genes on
any rotation of it must give the same protein set. The measured stake is 160,375
origin-spanning ORFs (1.72% of the collection) that a naive linear caller either truncates
or misses entirely.
"""
import random
from darkorf import circular


def test_topologies_that_count_as_circular():
    assert circular.is_circular("circular")
    assert circular.is_circular("Circular")
    assert not circular.is_circular("linear")


def test_rotation_preserves_length_and_content():
    sequence = "ATGCGTACGT"
    rotated = circular.rotate(sequence, 4)
    assert len(rotated) == len(sequence)
    assert sorted(rotated) == sorted(sequence)
    assert circular.rotate(rotated, len(sequence) - 4) == sequence


def test_overlap_is_capped_for_long_plasmids():
    """The duplicated head must be long enough to contain any real gene crossing the origin,
    but duplicating a 400 kb megaplasmid wholesale would double the gene-calling cost."""
    assert circular.overlap_for(1_000) <= 1_000
    assert circular.overlap_for(1_000_000) == circular.MAX_OVERLAP_BP


def test_genes_duplicated_by_the_overlap_are_resolved_to_one_copy():
    """A gene lying inside the duplicated head is called twice - once at its true
    coordinates and once shifted by the sequence length. Exactly one copy must survive, or
    every downstream occurrence count is inflated."""
    length = 1_000
    genes = [
        {"start": 10, "end": 100, "strand": 1},
        {"start": 10 + length, "end": 100 + length, "strand": 1},
    ]
    resolved = circular.resolve_origin_genes(genes, length)
    assert len(resolved) == 1


def test_an_origin_spanning_gene_is_marked_and_keeps_unrotated_coordinates():
    """Spec §8.4: do not assume end > start. The occurrence id is built from these
    coordinates, so rewriting them to look linear would break the identifier."""
    length = 1_000
    genes = [{"start": 960, "end": 1_040, "strand": 1}]
    resolved = circular.resolve_origin_genes(genes, length)
    assert len(resolved) == 1
    assert resolved[0]["origin_spanning"] is True
    assert resolved[0]["start"] == 960
    assert resolved[0]["end"] == 40


def test_rotation_invariance_on_a_random_sequence():
    """The property the whole module exists to satisfy."""
    random.seed(7)
    sequence = "".join(random.choice("ACGT") for _ in range(3_000))
    baseline = circular.rotate(sequence, 0)
    for offset in (1, 137, 1_500, 2_999):
        assert sorted(circular.rotate(sequence, offset)) == sorted(baseline)
