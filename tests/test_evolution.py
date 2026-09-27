"""Rule family_evolution: yn00 dN/dS, codon projection, consensus and RNAcode."""
import pytest
from conftest import requires

from plasmidann.evolution import back_translate, consensus, parse_yn00, yn00
from plasmidann.evolution_worker import rnacode

# --- yn00 -----------------------------------------------------------------------------

# The lines of a yn00 4.10.7 output file that parse_yn00 reads.
YN00_OUT = """YN00         aln.phy

ns =   3\tls =   9

(B) Yang & Nielsen (2000) method

seq. seq.     S       N        t   kappa   omega     dN +- SE    dS +- SE

   2    1     8.2    18.8   0.2686  4.6000  0.0000 -0.0000 +- 0.0000  0.2938 +- 0.2212
   3    1     5.1    21.9   0.1152  4.6000 99.0000 0.0473 +- 0.0480 -0.0000 +- 0.0000
   3    2    -nan    -nan     -nan  4.6000    -nan   -nan +-   -nan    -nan +-   -nan


(C) LWL85, LPB93 & LWLm methods

LWL85:  dS =  0.5748 dN =  0.0000 w = 0.0000 S =    7.0 N =   20.0
"""


def test_yn00_pairs_are_read_as_yn00_printed_them():
    """dS = 0 gives omega 99.0000 and an inestimable pair gives nan; both are kept."""
    codons, pairs = parse_yn00(YN00_OUT, ["a", "b", "c"])
    assert codons == 9
    assert pairs[0] == {"seq1": "b", "seq2": "a", "S": "8.2", "N": "18.8", "t": "0.2686",
                        "kappa": "4.6000", "omega": "0.0000", "dN": "-0.0000",
                        "dN_SE": "0.0000", "dS": "0.2938", "dS_SE": "0.2212"}
    assert (pairs[1]["omega"], pairs[1]["dN"], pairs[1]["dS"]) == ("99.0000", "0.0473",
                                                                   "-0.0000")
    assert pairs[2]["omega"] == "-nan"


def test_yn00_output_without_the_yang_nielsen_table_is_no_result():
    assert parse_yn00("ns =   2\tls =   9\n", ["a", "b"]) is None


ALN4 = {"m0": "ATGCTAGCTAGCAAACTAGCTAGCAAACTA",
        "m1": "ATGCTGGCTAGCAAACTAGCTAGCAAACTG",
        "m2": "ATGATAGCTAGCAAACTAGCTAGCAAACTA",
        "m3": "ATGCTAGCTAGC---CTAGCTAGCAAACTA"}


@requires("yn00")
def test_yn00_runs_every_pair_in_one_call(tmp_path):
    """Four members give six pairs; yn00 drops the gapped codon column from all of them."""
    codons, pairs = yn00(ALN4, tmp_path)
    assert codons == 9
    assert {(p["seq1"], p["seq2"]) for p in pairs} == {
        ("m1", "m0"), ("m2", "m0"), ("m2", "m1"), ("m3", "m0"), ("m3", "m1"), ("m3", "m2")}
    by = {(p["seq1"], p["seq2"]): p for p in pairs}
    assert by[("m1", "m0")]["dS"] == "0.2938"
    # m3 equals m0 once the gap column is dropped: yn00 reports omega 99 at dS = 0.
    assert by[("m3", "m0")]["omega"] == "99.0000"


@requires("yn00")
def test_a_table_4_gene_with_an_inner_tga_is_run_under_the_table_4_code(tmp_path):
    """Under the standard code yn00 exits with an error on the inner TGA."""
    aln = dict(ALN4, m1="ATGCTGGCTAGCAAATGAGCTAGCAAACTG")
    assert yn00(aln, tmp_path) is not None


def test_a_failing_yn00_is_no_result(tmp_path, monkeypatch):
    exe = tmp_path / "bin" / "yn00"
    exe.parent.mkdir()
    exe.write_text("#!/bin/sh\nexit 255\n")
    exe.chmod(0o755)
    monkeypatch.setenv("PATH", f"{exe.parent}:/usr/bin:/bin")
    assert yn00(ALN4, tmp_path / "w") is None


# --- codon alignment from a protein alignment ----------------------------------------

def test_back_translation_places_gaps_on_codon_boundaries():
    """mafft aligns proteins; dN/dS needs codons. A gap in the protein alignment becomes
    exactly three nucleotide gaps, never one or two, or the reading frame is destroyed."""
    assert back_translate("M-K", "ATGAAA") == "ATG---AAA"


def test_back_translation_refuses_a_cds_that_is_too_short():
    with pytest.raises(ValueError):
        back_translate("MKV", "ATGAAA")


# --- the family-consensus re-check ----------------------------------------------------

def test_the_consensus_is_the_commonest_residue_per_column():
    """The consensus carries the family's shared signal, so it is what gets re-searched."""
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


# --- RNAcode ----------------------------------------------------------------------------

def _stand_in_rnacode(tmp_path, monkeypatch, body):
    exe = tmp_path / "bin" / "RNAcode"
    exe.parent.mkdir()
    exe.write_text(f"#!/bin/sh\n{body}\n")
    exe.chmod(0o755)
    monkeypatch.setenv("PATH", f"{exe.parent}:/usr/bin:/bin")


ALN = {"m0": "ATGAAAGTG", "m1": "ATGAAGGTG"}


def test_rnacode_best_p_is_read_per_strand(tmp_path, monkeypatch):
    rows = ["1 + 1 30 1 30 m0 1 90 5.1 0.001", "2 - 2 30 1 30 m0 1 90 3.0 0.2",
            "3 + 3 30 1 30 m0 1 90 4.0 0.01"]
    _stand_in_rnacode(tmp_path, monkeypatch, "printf '" + "\\n".join(rows) + "\\n'")
    assert rnacode(ALN, tmp_path / "a.aln") == (0.001, 0.2, "MEASURED")


def test_rnacode_that_reports_nothing_is_no_signal(tmp_path, monkeypatch):
    _stand_in_rnacode(tmp_path, monkeypatch, "exit 0")
    assert rnacode(ALN, tmp_path / "a.aln") == (None, None, "NO_SIGNAL")


def test_rnacode_that_prints_an_error_and_exits_0_is_no_output(tmp_path, monkeypatch):
    _stand_in_rnacode(tmp_path, monkeypatch,
                      "echo 'ERROR: Unknown alignment file format'; exit 0")
    assert rnacode(ALN, tmp_path / "a.aln") == (None, None, "NO_OUTPUT")


def test_rnacode_that_exits_non_zero_is_no_output(tmp_path, monkeypatch):
    _stand_in_rnacode(tmp_path, monkeypatch, "exit 3")
    assert rnacode(ALN, tmp_path / "a.aln") == (None, None, "NO_OUTPUT")


def test_a_missing_rnacode_is_no_output(tmp_path, monkeypatch):
    monkeypatch.setenv("PATH", str(tmp_path))
    assert rnacode(ALN, tmp_path / "a.aln") == (None, None, "NO_OUTPUT")
