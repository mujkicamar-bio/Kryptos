"""Origin handling for circular plasmids.

A circular sequence has no privileged starting point, so a gene that crosses the point
where the circle was cut must be called whole, once, whatever the rotation.
"""
import random

from darkorf import circular, genecall

SENSE = [a + b + c for a in "ACGT" for b in "ACGT" for c in "ACGT"
         if a + b + c not in ("TAA", "TAG", "TGA")]


def rotate(seq, offset):
    """The same circle with `offset` as its first base."""
    return seq[offset:] + seq[:offset]


def proteins(seq, topology="circular"):
    return sorted(g["seq"] for g in genecall.call_genes("p", seq, topology)[1])


def test_topologies_that_count_as_circular():
    assert circular.is_circular("circular")
    assert circular.is_circular("Circular")
    assert not circular.is_circular("linear")


def test_the_whole_molecule_is_appended_up_to_the_cap():
    """A gene crossing the origin can be almost as long as the molecule, so a small
    molecule is appended whole; a large one gets MAX_OVERLAP_BP."""
    assert circular.overlap_for(1_906) == 1_906
    assert circular.overlap_for(1_000_000) == circular.MAX_OVERLAP_BP


def test_genes_duplicated_by_the_overlap_are_resolved_to_one_copy():
    """A gene inside the appended head is called twice, once shifted by the record length.
    Exactly one copy survives, at the record's own coordinates."""
    length = 1_000
    genes = [
        {"start": 10, "end": 100, "strand": 1, "partial": 0},
        {"start": 10 + length, "end": 100 + length, "strand": 1, "partial": 0},
    ]
    resolved = circular.resolve_origin_genes(genes, length, 2 * length)
    assert [(g["start"], g["end"]) for g in resolved] == [(10, 100)]


def test_the_head_fragment_of_an_origin_spanning_gene_is_dropped():
    """The caller sees the record start as an edge and emits the head of a gene crossing
    the cut as a partial ORF at coordinate 1-3. The extension calls the same gene whole;
    only the whole copy is kept."""
    length = 2_906
    genes = [
        {"start": 2, "end": 124, "strand": -1, "partial": 1},
        {"start": 2_626, "end": 3_030, "strand": -1, "partial": 0},
    ]
    resolved = circular.resolve_origin_genes(genes, length, 2 * length)
    assert len(resolved) == 1
    assert resolved[0]["origin_spanning"] is True
    assert (resolved[0]["start"], resolved[0]["end"]) == (2_626, 124)


def test_an_intact_gene_at_the_record_start_is_kept():
    genes = [{"start": 1, "end": 300, "strand": 1, "partial": 0}]
    assert len(circular.resolve_origin_genes(genes, 1_000, 2_000)) == 1


def test_an_origin_spanning_gene_is_marked_and_keeps_unrotated_coordinates():
    """start > end marks a gene read as start..length followed by 1..end."""
    genes = [{"start": 960, "end": 1_040, "strand": 1, "partial": 0}]
    resolved = circular.resolve_origin_genes(genes, 1_000, 2_000)
    assert [(g["start"], g["end"], g["origin_spanning"]) for g in resolved] == [
        (960, 40, True)]


def test_a_gene_near_the_record_start_called_as_an_edge_partial_is_not_lost():
    """Raw calls of COMPASS_AY362554.1 (5,674 bp, circular) on the record plus 2,837
    appended bases: the gene at 46..351 is called as a left-edge partial (1..351), and whole
    in the appended copy (5720..6025). No gene crosses the cut. The whole copy is kept,
    at the record's own coordinates."""
    length = 5_674
    genes = [
        {"start": 1, "end": 351, "strand": 1, "partial": 1},
        {"start": 348, "end": 1_757, "strand": 1, "partial": 0},
        {"start": 5_720, "end": 6_025, "strand": 1, "partial": 0},
        {"start": 6_022, "end": 7_431, "strand": 1, "partial": 0},
    ]
    resolved = circular.resolve_origin_genes(genes, length, length + 2_837)
    assert sorted((g["start"], g["end"], g["partial"]) for g in resolved) == [
        (46, 351, 0), (348, 1_757, 0)]


def test_a_whole_gene_called_only_in_the_appended_copy_is_kept():
    """Raw calls near position 1 of COMPASS_NC_007960.1 (188,318 bp, circular, 5,000 bases
    appended): the record gives only a left-edge partial on the + strand (2..628); the
    appended copy calls a complete 184-residue gene on the - strand (188488..189042) that
    has no twin in the record. The complete gene is kept at 170..724; a complete copy call
    that overlaps a kept record call is not."""
    length = 188_318
    genes = [
        {"start": 2, "end": 628, "strand": 1, "partial": 1},
        {"start": 1_000, "end": 1_600, "strand": 1, "partial": 0},
        {"start": 188_488, "end": 189_042, "strand": -1, "partial": 0},
        {"start": 189_200, "end": 189_500, "strand": -1, "partial": 0},
    ]
    resolved = circular.resolve_origin_genes(genes, length, length + 5_000)
    assert sorted((g["start"], g["end"], g["strand"]) for g in resolved) == [
        (170, 724, -1), (1_000, 1_600, 1)]


def test_a_partial_call_with_no_twin_is_dropped():
    """A left-edge partial with no call in the same frame in the appended copy exists only
    because edge genes need no start codon (COMPASS_EU999782.1, 1,549 bp)."""
    genes = [{"start": 1, "end": 153, "strand": 1, "partial": 1},
             {"start": 914, "end": 1_390, "strand": 1, "partial": 0}]
    resolved = circular.resolve_origin_genes(genes, 1_549, 3_098)
    assert [(g["start"], g["end"]) for g in resolved] == [(914, 1_390)]


def test_a_truncated_origin_spanning_gene_is_not_reported_intact():
    """A gene longer than the appended sequence runs off the right edge. It is kept as a
    partial, never reported as intact."""
    length = 20_000
    genes = [
        {"start": 3, "end": 5_600, "strand": 1, "partial": 1},
        {"start": 19_001, "end": 25_000, "strand": 1, "partial": 1},
    ]
    resolved = circular.resolve_origin_genes(genes, length, 25_000)
    assert len(resolved) == 1 and resolved[0]["partial"] == 1


def test_an_origin_spanning_gene_longer_than_half_the_molecule_is_called_whole():
    """A 1,206-bp gene on a 1,906-bp circle, the size of a small cryptic plasmid with one
    Rep gene: every rotation, including those that cut the gene, gives the whole
    401-residue protein."""
    genecall.configure(63)
    rng = random.Random(5)
    gene = "ATG" + "".join(rng.choice(SENSE) for _ in range(400)) + "TAA"
    molecule = "".join(rng.choice("ACGT") for _ in range(700)) + gene
    truth = proteins(molecule)
    assert 401 in {len(p) for p in truth}
    for offset in (800, 1_300, 1_800):
        assert proteins(rotate(molecule, offset)) == truth, offset


def test_inverted_terminal_repeat_is_linear():
    """An inverted terminal repeat marks a genuinely linear replicon with hairpin or
    protein-capped telomeres (Borrelia and Streptomyces linear plasmids); joining its ends
    would invent a gene across a junction that does not exist in the cell."""
    assert not circular.is_circular("inverted terminal repeat")
    assert circular.is_circular("direct terminal repeat")


def test_an_unknown_topology_is_treated_as_linear():
    """Treating a circular molecule as linear loses the genes that cross the origin;
    treating a linear molecule as circular invents a chimeric gene. Unknown is linear."""
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
    genecall.configure(90)
    # Seed 2 puts a gene across the junction, so the untrimmed record really does call a
    # different protein set: the test fails if the repeat is not removed.
    rng = random.Random(2)
    molecule = "".join(rng.choice("ACGT") for _ in range(6_000))
    repeat = 101                                   # not a multiple of 3
    record = molecule + molecule[:repeat]
    k = circular.terminal_repeat_length(record, 20)
    assert k == repeat
    truth = proteins(rotate(molecule, 2_345))
    assert proteins(record[:-k], "direct terminal repeat") == truth
    assert proteins(record, "direct terminal repeat") != truth


def test_a_gene_running_off_the_left_edge_is_partial_at_its_begin_only():
    """A linear record that starts inside a 300-codon open reading frame with no ATG: the
    call is partial at its begin (position 1), complete at its end, and says which."""
    genecall.configure(63)
    rng = random.Random(3)
    body = "".join(rng.choice([c for c in SENSE if c != "ATG"]) for _ in range(300))
    seq = body + "TAA" + "".join(rng.choice("ACGT") for _ in range(300))
    (gene,) = genecall.call_genes("p", seq, "linear")[1]
    assert (gene["start"], gene["partial"], gene["partial_begin"], gene["partial_end"]) == (
        1, 1, 1, 0)
