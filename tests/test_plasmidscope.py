"""PlasmidScope label transfer: the class, the join key, and the choice between occurrences."""
from plasmidann.dereplicate import _seq_id
from plasmidann.plasmidscope import annot_label, protein_seq_id, ps_class, reduce_rows


def _row(seq="MKVLA*", cog_cat="S", cog_id="-", ko="-", pfams="-", ec="-",
         source="Prodigal:2.6"):
    """A PlasmidScope ALL.protein_list row, with its column names."""
    return {"Sequence": seq, "COG_category": cog_cat, "COG_id": cog_id, "KEGG_ko": ko,
            "PFAMs": pfams, "EC_number": ec, "KEGG_Pathway": "-", "GOs": "-",
            "Orf Prediction Source": source}


def test_the_field_patterns_of_the_real_table():
    """The patterns counted in ALL.protein_list: nothing but category S; a COG in S only;
    a named Pfam, a KO or an EC with the product still 'hypothetical'."""
    assert ps_class(_row()) == "NONE"
    assert ps_class(_row(cog_id="COG4710")) == "UNKNOWN_ORTHOLOG"
    assert ps_class(_row(pfams="RHH_1")) == "ANNOTATED"
    assert ps_class(_row(ko="ko:K18918")) == "ANNOTATED"
    assert ps_class(_row(ec="2.7.7.7")) == "ANNOTATED"


def test_unknown_function_evidence_is_not_an_annotation():
    """The cascade's own bar: a label must name a function. A domain of unknown function,
    or a COG category letter with nothing else, names none - such a protein must be
    searched, not declared FUNCTIONAL and removed from the dark set."""
    assert ps_class(_row(pfams="DUF1845")) == "UNKNOWN_ORTHOLOG"
    assert ps_class(_row(pfams="DUF1845,UPF0102", cog_id="COG3befr")) == "UNKNOWN_ORTHOLOG"
    assert ps_class(_row(cog_cat="L", cog_id="COG2026")) == "UNKNOWN_ORTHOLOG"
    # One named domain beside a DUF is enough.
    assert ps_class(_row(pfams="DUF1845,Relaxase")) == "ANNOTATED"


def test_the_join_key_is_our_seq_id_with_or_without_the_stop():
    """Our ORFs are stored without '*'. A PlasmidScope sequence that keeps it must still
    hash to the same seq_id, or no protein would ever match."""
    assert protein_seq_id(_row(seq="MKVLA*")) == _seq_id("MKVLA")
    assert protein_seq_id(_row(seq="MKVLA")) == _seq_id("MKVLA")


def test_an_identical_sequence_takes_its_most_informative_occurrence():
    """The same protein on a Prokka-called and on a PGAP-deposited plasmid can receive
    different eggNOG results. The annotated occurrence wins, whichever comes first."""
    rows = [_row(), _row(pfams="RHH_1", source="Protein Homology"), _row(cog_id="COG1")]
    out = reduce_rows(rows, wanted={_seq_id("MKVLA")})
    rec = out[_seq_id("MKVLA")]
    assert rec["ps_class"] == "ANNOTATED"
    assert rec["pfams"] == "RHH_1"
    assert rec["orf_source"] == "Protein Homology"
    assert rec["n_orfs"] == 3


def test_proteins_we_do_not_have_are_dropped_and_placeholders_cleared():
    out = reduce_rows([_row(seq="MKVLA*"), _row(seq="MOTHER*")], wanted={_seq_id("MKVLA")})
    assert list(out) == [_seq_id("MKVLA")]
    assert out[_seq_id("MKVLA")]["kegg_ko"] == "", "the '-' placeholder became a term"


def test_the_label_is_the_first_named_pfam_then_ko_then_ec():
    assert annot_label({"pfams": "RHH_1,ParE_toxin", "kegg_ko": "ko:K1"}) == "RHH_1"
    assert annot_label({"pfams": "DUF1845,Relaxase", "kegg_ko": "ko:K1"}) == "Relaxase"
    assert annot_label({"pfams": "DUF1845", "kegg_ko": "ko:K1"}) == "ko:K1"
    assert annot_label({"pfams": "", "kegg_ko": "", "ec": "2.7.7.7"}) == "2.7.7.7"
