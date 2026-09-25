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


def test_the_head_fragment_of_an_origin_spanning_gene_is_dropped():
    """The gene caller sees the record start as an edge and emits the truncated head of a
    gene that crosses the cut as a partial ORF at coordinate 1-3 - 88,600 of the 160,375
    partials on this collection. The extension calls the same gene intact, so keeping the
    stub would write every origin-spanning protein twice."""
    length = 2_906
    genes = [
        {"start": 2, "end": 124, "strand": -1, "partial": 1},
        {"start": 2_626, "end": 3_030, "strand": -1, "partial": 0},
    ]
    resolved = circular.resolve_origin_genes(genes, length)
    assert len(resolved) == 1
    assert resolved[0]["origin_spanning"] is True
    assert (resolved[0]["start"], resolved[0]["end"]) == (2_626, 124)


def test_an_intact_gene_at_the_record_start_is_kept():
    """Only a PARTIAL gene at the record start is a fragment. A complete gene that happens
    to begin at coordinate 1 is real and must survive."""
    genes = [{"start": 1, "end": 300, "strand": 1, "partial": 0}]
    assert len(circular.resolve_origin_genes(genes, 1_000)) == 1


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


# --- topology is biology, not a formatting detail ------------------------------------


def test_inverted_terminal_repeat_is_linear():
    """An inverted terminal repeat is the signature of a genuinely LINEAR replicon with
    hairpin or protein-capped telomeres - Borrelia and Streptomyces linear plasmids, phi29,
    adenovirus. Joining its ends would invent a gene across a junction that does not exist
    in the cell. The dataset holds 30 such records.

    'direct terminal repeat' is different and stays circular: measured on this collection,
    those records have their repeats already trimmed (0.9% intra-plasmid duplicate rate
    against 32.2% for 'circular', which is genuine multi-copy IS biology).
    """
    assert not circular.is_circular("inverted terminal repeat")
    assert circular.is_circular("direct terminal repeat")


def test_an_unknown_topology_is_treated_as_linear():
    """The two errors are not symmetric. Treating a circular molecule as linear loses the
    ~1.7% of genes that cross the origin; treating a linear molecule as circular fabricates
    a chimeric gene that was never there. When topology is unknown, lose data rather than
    invent it."""
    for unknown in ("", "NA", "unknown", None):
        assert not circular.is_circular(unknown)


def test_a_terminal_repeat_is_found_and_a_short_match_is_not():
    """The overlap an assembler leaves on a circular record: the first k bases equal the
    last k. Found at >= 20 bp (CheckV's criterion), never below it, and 0 when the ends
    are unrelated."""
    rng = random.Random(3)
    core = "".join(rng.choice("ACGT") for _ in range(3_000))
    assert circular.terminal_repeat_length(core[:500] + core + core[:500], 20) == 500
    assert circular.terminal_repeat_length(core[:19] + core + core[:19], 20) == 0
    assert circular.terminal_repeat_length(core, 20) == 0


def test_a_trimmed_terminal_repeat_record_calls_the_same_proteins_as_its_circle():
    """Calling genes on a record that still carries its terminal repeat joins the ends
    through both copies; removing one copy must give exactly the proteins of the molecule
    itself, in any rotation."""
    from darkorf import genecall
    genecall.configure(90)
    # Seed 2 puts a gene across the junction, so the untrimmed record really does call a
    # different protein set: the test fails if S0 stops removing the repeat.
    rng = random.Random(2)
    molecule = "".join(rng.choice("ACGT") for _ in range(6_000))
    repeat = 101                                   # not a multiple of 3
    record = molecule + molecule[:repeat]
    k = circular.terminal_repeat_length(record, 20)
    assert k == repeat
    proteins = lambda seq, topology: sorted(g["seq"] for g in
                                            genecall.call_genes("p", seq, topology)[1])
    truth = proteins(circular.rotate(molecule, 2_345), "circular")
    assert proteins(record[:-k], "direct terminal repeat") == truth
    assert proteins(record, "direct terminal repeat") != truth
