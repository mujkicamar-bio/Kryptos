"""Parsing HH-suite output for the PHROGs tier (spec sections 18 and 79).

WHY HH-SUITE AND NOT MMseqs2

PHROGs distributes HH-suite HHM profiles and documents HH-suite as the way to search them.
An earlier attempt converted them to an MMseqs2 profile database instead, to stay inside
section 79's exclusion of HMM-HMM search. That was abandoned: the conversion lost the PHROG
identifier - `mmseqs convertprofiledb` takes each profile's header from its NAME line, which
names the seed protein, and assigns database keys in an order that does not follow the
ffindex - and searching 38,880 converted profiles never completed on the test set.

Using the profiles as their authors distribute them is the decision recorded in the spec.
The PHROG identifier is then native: it is the ffindex key.

WHY THE .hhr FORMAT AND NOT -blasttab

Section 18 requires phrog_score, and HH-suite's PROBABILITY is the statistic PHROGs' own
work reports and filters on. -blasttab gives an E-value and a bit score but no probability,
so the native hit table is what carries the required field.
"""
import pytest

from plasmidann import phrogs


# A real hhsearch/hhblits hit table. The header block, the column header, the hits, then
# the per-hit alignments, which are not parsed.
HHR = """\
Query         2171fc6c84273f21a14035f8eb570b4c
Match_columns 480
No_of_seqs    1 out of 1
Neff          1
Searched_HMMs 38880
Date          Fri Sep 18 14:02:11 2026
Command       hhsearch -i query.a3m -d data/refs/phrogs/phrogs

 No Hit                             Prob E-value P-value  Score    SS Cols Query HMM  Template HMM
  1 phrog_1234 terminase large sub 100.0 1.2E-45 2.6E-50  350.2   0.0  250    5-260     1-255 (687)
  2 phrog_99 hypothetical protein   87.3   0.004 1.0E-07   45.1   0.0   60  100-160    12-72  (120)
  3 phrog_7 tail fiber              22.0     8.1   0.00021  20.3   0.0   30  300-330     1-31  (95)

No 1
>phrog_1234 terminase large subunit
Probab=100.00  E-value=1.2e-45  Score=350.24  Aligned_cols=250
"""


def test_the_hit_table_is_parsed_with_every_field_section_18_requires():
    hits = phrogs.parse_hhr(HHR)

    assert len(hits) == 3
    top = hits[0]
    assert top["query"] == "2171fc6c84273f21a14035f8eb570b4c"
    assert top["phrog_id"] == "phrog_1234"
    assert top["phrog_description"] == "terminase large sub"
    assert top["phrog_score"] == 350.2
    assert top["phrog_probability"] == 100.0
    assert top["phrog_evalue"] == "1.2E-45"
    assert top["query_start"] == 5
    assert top["query_end"] == 260


def test_query_coverage_is_measured_against_the_query_not_the_profile():
    """Section 18 asks for phrog_query_coverage. The Query HMM range is the span on OUR
    protein; the Template HMM range and its parenthesised length belong to the PHROG.
    Mixing them up produces a coverage above 1 whenever the profile is shorter than the
    protein, and nothing downstream would reject it."""
    hits = phrogs.parse_hhr(HHR)

    # 5-260 of a 480-column query.
    assert hits[0]["phrog_query_coverage"] == round(256 / 480, 4)


def test_a_search_with_no_hits_yields_no_rows():
    """hhsearch writes the header and the column line and nothing else when a query
    matches nothing. That is the commonest outcome for a dark protein and must read as
    zero rows rather than as a parse failure."""
    empty = "\n".join(HHR.splitlines()[:8]) + "\n"

    assert phrogs.parse_hhr(empty) == []


def test_the_probability_filter_is_applied_and_is_configuration():
    """PHROGs' own annotation work filters on HH-suite probability. The threshold is a
    parameter rather than a constant here, because it governs what enters the annotated
    set and therefore what the DARK set is the complement of."""
    hits = phrogs.parse_hhr(HHR, min_probability=80.0)

    assert [h["phrog_id"] for h in hits] == ["phrog_1234", "phrog_99"], (
        "the 22.0-probability hit is below the threshold and must not be reported")


def test_a_description_may_be_empty_and_the_id_still_parses():
    """Most PHROGs NAME lines are 'hypothetical protein' and some carry nothing at all.
    An absent description must not shift the columns, because every field after it -
    probability, E-value, score, the coordinates - would then be read from the wrong
    place."""
    text = ("Query         q1\nMatch_columns 100\n\n"
            " No Hit                             Prob E-value P-value  Score    SS Cols "
            "Query HMM  Template HMM\n"
            "  1 phrog_5                         99.9 1.0E-30 2.0E-35  200.0   0.0  90  "
            "  1-90      1-90  (90)\n")

    hits = phrogs.parse_hhr(text)

    assert len(hits) == 1
    assert hits[0]["phrog_id"] == "phrog_5"
    assert hits[0]["phrog_description"] == ""
    assert hits[0]["phrog_probability"] == 99.9
    assert hits[0]["query_start"] == 1


def test_the_ffindex_key_suffix_is_stripped_from_the_identifier():
    """ffindex_build keys entries by file name, so the database reports 'phrog_1234.hhm'.
    The identifier the spec asks for, and the key that joins to the PHROGs annotation
    table, is the stem."""
    text = ("Query         q1\nMatch_columns 100\n\n"
            " No Hit                             Prob E-value P-value  Score    SS Cols "
            "Query HMM  Template HMM\n"
            "  1 phrog_1234.hhm terminase        99.9 1.0E-30 2.0E-35  200.0   0.0  90  "
            "  1-90      1-90  (90)\n")

    assert phrogs.parse_hhr(text)[0]["phrog_id"] == "phrog_1234"


# ------------------------------------------------------------------------------------
# hhblits_omp writes its results as an ffindex rather than as files, so reading them back
# is part of the stage.
# ------------------------------------------------------------------------------------

def test_ffindex_entries_are_read_back_by_name(tmp_path):
    """ffindex is <name>\\t<offset>\\t<length>, and each entry is NUL-terminated with the
    recorded length INCLUDING that byte. Reading the recorded length verbatim therefore
    puts a stray NUL at the end of every entry, which for a .hhr means the final hit row
    never matches its pattern."""
    data = tmp_path / "res.ffdata"
    index = tmp_path / "res.ffindex"
    first, second = b"hello\0", b"world!\0"
    data.write_bytes(first + second)
    index.write_text(f"q1\t0\t{len(first)}\nq2\t{len(first)}\t{len(second)}\n")

    entries = dict(phrogs.read_ffindex(data, index))

    assert entries == {"q1": "hello", "q2": "world!"}


def test_an_empty_ffindex_reads_as_nothing(tmp_path):
    """A tier whose every query matched nothing still produces the pair of files."""
    data = tmp_path / "res.ffdata"
    index = tmp_path / "res.ffindex"
    data.write_bytes(b"")
    index.write_text("")

    assert list(phrogs.read_ffindex(data, index)) == []


# ------------------------------------------------------------------------------------
# Building the prefilter from what PHROGs ships.
#
# HH-suite's cstranslate builds the cs219 prefilter from an ALIGNMENT, and its 3.3 release
# accepts prf, seq, fas, a2m, a3m or ca3m - not HHM. PHROGs ships HHM only. But an HHM
# written by hhmake carries the alignment it was built from in its SEQ block, as a3m, so
# the documented input exists inside the file the tool refuses to read.
# ------------------------------------------------------------------------------------

HHM = """\
HHsearch 1.5
NAME  p205037 VI_04338
FAM   
FILE  phrog_10
LENG  120 match states, 300 columns in multiple alignment
SEQ
>Consensus
xxxMKLTxxxALL
>p205037 VI_04338
--MMKLTEKQALL
>other_seq
MKMMKLTeKQ-LL
#
NULL   3706	5728	4211	4064	4839	3729	4763	4308	4069	3323	5509	4640	4464	4937	4285	4423	3815	3783	6325	4665
HMM    A	C	D	E	F	G	H	I	K	L	M	N	P	Q	R	S	T	V	W	Y
//
"""


def test_the_alignment_is_extracted_from_the_hhm_seq_block():
    a3m = phrogs.a3m_from_hhm(HHM)

    assert a3m.startswith(">Consensus\n"), "the consensus is the master sequence"
    assert ">p205037 VI_04338\n--MMKLTEKQALL\n" in a3m
    assert ">other_seq\nMKMMKLTeKQ-LL\n" in a3m, "lowercase inserts must survive"
    assert "NULL" not in a3m and "HMM" not in a3m, "profile lines are not alignment"


def test_an_hhm_without_a_seq_block_is_an_error_not_an_empty_alignment():
    """An empty a3m entry gives cstranslate nothing to build a state sequence from, and
    hhblits then treats that profile as unsearchable while the hhm entry still exists.
    A profile silently missing from the prefilter is the failure the build must refuse."""
    with pytest.raises(ValueError):
        phrogs.a3m_from_hhm("HHsearch 1.5\nNAME  x\nLENG  5\n#\nHMM ...\n//\n")
