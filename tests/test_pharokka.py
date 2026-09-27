"""Parsing pharokka's protein-mode output for the phage tier.

Built against the real output of pharokka 1.10.1 on 200 proteins from the test set. A
protein with no hit carries 'No_MMseqs' in the `phrog` column, so a family hit is
recognised by an integer there; when both searches found the family, the stronger E-value
is reported.
"""
from plasmidann import pharokka

HEADER = ("ID\tlength\tphrog\tannot\tcategory\tmmseqs_phrog\tmmseqs_alnScore\t"
          "mmseqs_seqIdentity\tmmseqs_eVal\tpyhmmer_phrog\tpyhmmer_bitscore\t"
          "pyhmmer_evalue\tcolor\tvfdb_hit\tvfdb_alnScore\tvfdb_seqIdentity\tvfdb_eVal\t"
          "vfdb_species\tvfdb_short_name\tvfdb_description\tCARD_hit\tCARD_alnScore\t"
          "CARD_seqIdentity\tCARD_eVal\tCARD_species\tARO_Accession\tCARD_short_name\t"
          "Protein_Accession\tDNA_Accession\tAMR_Gene_Family\tDrug_Class\t"
          "Resistance_Mechanism\n")

BOTH = ("0cab68e6\t250\t164\tParA-like partition protein\tDNA, RNA and nucleotide "
        "metabolism\t164\t105\t0.326\t3.543e-27\t164\t146.517853\t2.5937970813855673e-42\t"
        "#a6cee3\tNone\tNone\tNone\tNone\tNone\tNone\tNone\tNone\tNone\tNone\tNone\tNone\t"
        "None\tNone\tNone\tNone\tNone\tNone\tNone\n")
HMM_ONLY = ("aaaa0001\t138\t14168\tnuclease\tDNA, RNA and nucleotide metabolism\t"
            "No_MMseqs\tNo_MMseqs\tNo_MMseqs\tNo_MMseqs\t14168\t66.215622\t6.66e-18\t"
            "#fb9a99\tNone\tNone\tNone\tNone\tNone\tNone\tNone\tNone\tNone\tNone\tNone\t"
            "None\tNone\tNone\tNone\tNone\tNone\tNone\tNone\n")
HYPOTHETICAL = ("bbbb0002\t90\t9928\thypothetical protein\tunknown function\t9928\t40\t"
                "0.3\t1e-8\t9928\t30.1\t1e-9\t#ffffff\tNone\tNone\tNone\tNone\tNone\tNone\t"
                "None\tNone\tNone\tNone\tNone\tNone\tNone\tNone\tNone\tNone\tNone\tNone\t"
                "None\n")
NO_HIT = ("cccc0003\t77\tNo_MMseqs\thypothetical protein\tunknown function\tNo_MMseqs\t"
          "No_MMseqs\tNo_MMseqs\tNo_MMseqs\tNo_PHROGs_HMM\tNo_PHROGs_HMM\tNo_PHROGs_HMM\t"
          "None\tNone\tNone\tNone\tNone\tNone\tNone\tNone\tNone\tNone\tNone\tNone\tNone\t"
          "None\tNone\tNone\tNone\tNone\tNone\tNone\n")
CARD = ("dddd0004\t286\tNo_MMseqs\thypothetical protein\tunknown function\tNo_MMseqs\t"
        "No_MMseqs\tNo_MMseqs\tNo_MMseqs\tNo_PHROGs_HMM\tNo_PHROGs_HMM\tNo_PHROGs_HMM\t"
        "None\tNone\tNone\tNone\tNone\tNone\tNone\tNone\t"
        "gb|AAP20891.1|ARO:3000873|TEM-1\t560\t0.99\t1e-200\tEscherichia coli\t"
        "ARO:3000873\tTEM-1\tAAP20891.1\tAY293345.1\tTEM beta-lactamase\t"
        "cephalosporin;penam\tantibiotic inactivation\n")


def test_a_family_hit_carries_family_label_category_and_length():
    hits = pharokka.parse_merged(HEADER + BOTH)

    assert len(hits) == 1
    h = hits[0]
    assert h["query"] == "0cab68e6"
    assert h["source"] == "pharokka"
    assert h["family_id"] == "phrog_164"
    assert h["label"] == "ParA-like partition protein"
    assert h["category"] == "DNA, RNA and nucleotide metabolism"
    assert h["query_length"] == 250


def test_the_stronger_of_the_two_evalues_is_reported():
    """The family was found by both searches. The assignment rests on the stronger
    evidence, and the cascade ranks hits by E-value, so reporting the weaker one would
    let a Swiss-Prot hit outrank a family assignment that pyhmmer found at 1e-42."""
    h = pharokka.parse_merged(HEADER + BOTH)[0]

    assert h["evalue"] == "2.5937970813855673e-42"


def test_a_pyhmmer_only_hit_is_a_hit():
    h = pharokka.parse_merged(HEADER + HMM_ONLY)[0]

    assert h["family_id"] == "phrog_14168"
    assert h["evalue"] == "6.66e-18"


def test_no_hit_yields_no_row_whatever_sentinel_the_column_holds():
    """'No_MMseqs' in the phrog column is pharokka's own sentinel for 'nothing at all', a
    consequence of its null-fill order. A parser keyed on 'No_PHROG' would report
    'No_MMseqs' as a family identifier."""
    assert pharokka.parse_merged(HEADER + NO_HIT) == []


def test_a_hypothetical_family_is_still_a_hit():
    """Most families are named 'hypothetical protein'. The hit is real - the protein
    belongs to a family - and it is reported. Whether the LABEL is informative is
    cascade.is_informative's decision, made identically at every tier, so that a family
    membership with no name does not make a dark protein annotated."""
    h = pharokka.parse_merged(HEADER + HYPOTHETICAL)[0]

    assert h["family_id"] == "phrog_9928"
    assert h["label"] == "hypothetical protein"


def test_card_and_vfdb_hits_are_separate_rows_with_their_own_source():
    """pharokka searches CARD and VFDB in the same run. They are different databases with
    different meanings - an AMR gene family is not a phage family - so each is its own
    row with its own source, never folded into the phage hit."""
    hits = pharokka.parse_merged(HEADER + CARD)

    assert [h["source"] for h in hits] == ["card"]
    assert hits[0]["label"] == "TEM beta-lactamase"
    assert hits[0]["family_id"] == "ARO:3000873"
    assert hits[0]["evalue"] == "1e-200"


def test_no_hit_carries_coordinates():
    """pharokka deletes its raw alignment tables on exit, so no coordinates survive. The
    parser says so rather than inventing a full-length span: start and end are absent,
    and the cascade treats the hit as a family-level assignment whose completeness is
    NOT_MEASURED."""
    h = pharokka.parse_merged(HEADER + BOTH)[0]

    assert "start" not in h and "end" not in h
