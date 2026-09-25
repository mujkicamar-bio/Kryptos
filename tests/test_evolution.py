"""S7: dN/dS - the strongest single piece of evidence that a dark ORF is a real protein.

It converts "nobody has named it" into "evolution is paying to keep it", which is a
completely different claim and the one a reviewer will ask for. FESNov required dN/dS < 0.5
and discarded 94.3% of its novel clusters on this and the other reality filters.

Nei-Gojobori counting is used rather than a codon model: it needs no tree, no optimiser and
no external process, so it runs over hundreds of thousands of small families and is fully
testable in-process.
"""
import pytest

from plasmidann.evolution import (
    back_translate,
    codon_differences,
    consensus,
    dnds,
    dnds_detail,
    synonymous_sites,
)

# --- site counting -------------------------------------------------------------------

def test_a_fourfold_degenerate_codon_has_one_synonymous_site():
    """GCx all encode alanine, and nothing else does, so exactly the third position is
    synonymous."""
    s, n = synonymous_sites("GCT")
    assert s == pytest.approx(1.0)
    assert s + n == pytest.approx(3.0)


def test_leucine_has_more_than_one_synonymous_site():
    """Leucine spans two codon families (CTx and TTA/TTG), so CTA is also synonymous at
    the FIRST position: CTA -> TTA is silent. Counting only the third position would
    understate synonymous sites and inflate dN/dS."""
    s, _ = synonymous_sites("CTA")
    assert s == pytest.approx(4 / 3)


def test_methionine_has_no_synonymous_sites():
    """ATG is the only codon for methionine: every change alters the amino acid."""
    s, n = synonymous_sites("ATG")
    assert s == pytest.approx(0.0)
    assert n == pytest.approx(3.0)


# --- pairwise differences ------------------------------------------------------------

def test_a_silent_third_position_change_is_synonymous():
    sd, nd = codon_differences("CTA", "CTG")   # Leu -> Leu
    assert (sd, nd) == (1, 0)


def test_an_amino_acid_change_is_nonsynonymous():
    sd, nd = codon_differences("ATG", "ATA")   # Met -> Ile
    assert (sd, nd) == (0, 1)


def test_identical_codons_differ_in_nothing():
    assert codon_differences("ATG", "ATG") == (0, 0)


# --- the ratio -----------------------------------------------------------------------

def test_identical_sequences_give_no_signal():
    """Zero divergence carries no information about selection, and must not be reported as
    dN/dS = 0, which would look like maximal purifying selection."""
    seq = "ATGCTAGCTAGCAAA"
    assert dnds(seq, seq) is None


def test_only_silent_changes_indicate_strong_purifying_selection():
    """Divergence is kept well below the Jukes-Cantor saturation point (pS = 0.75); at or
    beyond it the correction has no domain and the estimate is meaningless."""
    a = "CTA" * 10
    b = "CTA" * 7 + "CTG" + "CTT" + "CTC"    # 3 silent changes in 10 Leu codons
    ratio = dnds(a, b)
    assert ratio is not None and ratio == 0.0


def test_only_replacement_changes_give_a_ratio_above_one():
    a = "ATGATGATGATGATG"          # Met x5
    b = "ATAATAATAATAATA"          # Ile x5
    assert dnds(a, b) == float("inf")


def test_sequences_of_unequal_length_are_refused():
    with pytest.raises(ValueError):
        dnds("ATGATG", "ATG")


def test_a_sequence_that_is_not_a_whole_number_of_codons_is_refused():
    with pytest.raises(ValueError):
        dnds("ATGAT", "ATGAT")


def test_gapped_and_ambiguous_codons_are_skipped_not_counted():
    """Alignment gaps are missing data, not evidence of conservation. Counting a gapped
    column as identical would inflate apparent purifying selection precisely where the
    alignment is least trustworthy."""
    a = "CTA" * 4 + "---" + "CTA" * 5
    b = "CTA" * 2 + "CTG" + "CTT" + "---" + "CTA" * 5
    assert dnds(a, b) == 0.0


def test_sequences_too_diverged_to_correct_return_no_estimate():
    """Beyond pS = 0.75 the Jukes-Cantor correction has no domain. Returning a number
    there would be inventing one."""
    a = "CTA" * 5
    b = "CTG" * 3 + "CTT" + "CTC"     # every codon changed: saturated
    assert dnds(a, b) is None


# --- codon alignment from a protein alignment ----------------------------------------

def test_back_translation_places_gaps_on_codon_boundaries():
    """mafft aligns proteins; dN/dS needs codons. A gap in the protein alignment becomes
    exactly three nucleotide gaps, never one or two, or the reading frame is destroyed."""
    assert back_translate("M-K", "ATGAAA") == "ATG---AAA"


def test_back_translation_refuses_a_cds_that_is_too_short():
    with pytest.raises(ValueError):
        back_translate("MKV", "ATGAAA")


# --- distinguishing WHY there is no estimate -----------------------------------------

def test_no_divergence_is_reported_separately_from_a_failed_estimate():
    """A family of identical sequences has no divergence to measure. That is not the same
    as a family that was measured and found neutral, and it must not be scored as absence
    of evidence - clonal redundancy would otherwise silently penalise the most conserved
    families in the collection."""
    from plasmidann.evolution import dnds_detail

    value, status = dnds_detail("ATGCTAGCTAGCAAA", "ATGCTAGCTAGCAAA")
    assert value is None
    assert status == "NO_DIVERGENCE"


def test_saturation_is_reported_as_its_own_status():
    from plasmidann.evolution import dnds_detail

    value, status = dnds_detail("CTA" * 5, "CTG" * 3 + "CTT" + "CTC")
    assert value is None and status == "SATURATED"


def test_a_successful_estimate_is_reported_as_measured():
    from plasmidann.evolution import dnds_detail

    value, status = dnds_detail("CTA" * 10, "CTA" * 7 + "CTG" + "CTT" + "CTC")
    assert value == 0.0 and status == "MEASURED"


def test_too_little_usable_alignment_is_reported_as_its_own_status():
    from plasmidann.evolution import dnds_detail

    value, status = dnds_detail("CTA---", "CTG---")
    assert value is None and status == "TOO_SHORT"


# --- min_codons was declared in config and never applied ------------------------------

def test_a_pair_with_too_few_usable_codons_reports_no_estimate():
    """config/targets.yaml declares `min_codons: 20`, described as the point below which a
    dN/dS estimate is noise whatever it says. Nothing read it. The only floor actually
    applied was a hardcoded 3, so eight-codon fragments produced dN/dS values that were
    then used as `purifying_selection` - one of the four reality tests, and the strongest
    of them.

    Eight codons is above the absolute floor of 3 and below the declared 20."""
    a = "ATGAAAGTGCTGGCGACCACCCTG"          # 8 codons
    b = "ATGAAAGTACTGGCGACCACCCTA"          # same protein, two silent changes
    assert len(a) == len(b) == 24

    value, status = dnds_detail(a, b, min_codons=20)
    assert value is None
    assert status == "TOO_SHORT"


def test_the_same_pair_is_measurable_when_the_floor_allows_it():
    """The floor has to be the reason, not the alignment: at the default floor this pair
    yields an estimate, so the test above is measuring min_codons and nothing else."""
    a = "ATGAAAGTGCTGGCGACCACCCTG"
    b = "ATGAAAGTACTGGCGACCACCCTA"
    value, status = dnds_detail(a, b, min_codons=3)
    assert status == "MEASURED"
    assert value is not None


# --- the family-consensus re-check ----------------------------------------------------

def test_the_consensus_is_the_commonest_residue_per_column():
    """Pavlopoulos removed 6.5% of clusters this way: a family can be collectively
    recognisable while every member individually misses the per-sequence threshold. The
    consensus is what carries the family's shared signal, so it is what gets re-searched."""
    aln = {"a": "MKVL-AT", "b": "MKIL-AT", "c": "MKVLQAS"}
    assert consensus(aln) == "MKVLAT"


def test_a_majority_gap_column_is_dropped():
    """A column that is mostly gap is an insertion in one member, not part of the family's
    shared sequence. Keeping it would put a residue into the consensus that most members do
    not have, and the re-search would be of a sequence no member actually carries."""
    aln = {"a": "MK--VL", "b": "MK--VL", "c": "MKQQVL"}
    assert consensus(aln) == "MKVL"


def test_ties_are_broken_deterministically():
    """Two residues equally common must not resolve by dict ordering, or the consensus -
    and therefore the re-check verdict - changes between runs on identical input."""
    aln = {"a": "MA", "b": "MC"}
    assert consensus(aln) == consensus({"b": "MC", "a": "MA"})
    assert consensus(aln) == "MA"


def test_an_empty_alignment_has_no_consensus():
    assert consensus({}) == ""
