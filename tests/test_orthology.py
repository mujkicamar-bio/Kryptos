"""Parsing eggNOG-mapper's annotation table."""
from plasmidann.orthology import parse_annotations

HEADER = ("#query\tseed_ortholog\tevalue\tscore\teggNOG_OGs\tmax_annot_lvl\t"
          "COG_category\tDescription\tPreferred_name\tGOs\tEC\tKEGG_ko\tKEGG_Pathway\n")


def _rows(*lines):
    return HEADER + "".join(lines)


def test_a_protein_with_pathways_yields_them_as_a_list():
    text = _rows("P1\t83333.b0002\t1e-90\t500\tCOG0527@1\tBacteria\tE\t"
                 "Aspartokinase\tthrA\tGO:0004072\t2.7.2.4\tko:K00928\t"
                 "ko00260,ko01230,map00260\n")
    got = parse_annotations(text)
    assert got["P1"]["cog_category"] == "E"
    assert got["P1"]["kegg_pathways"] == ["ko00260", "ko01230", "map00260"]
    assert got["P1"]["preferred_name"] == "thrA"


def test_emappers_placeholder_dash_is_not_an_annotation():
    """eggNOG-mapper writes '-' for every field it has nothing for. Read literally that
    becomes a COG category named '-' and a pathway named '-', and both then aggregate into
    the neighbourhood composition as though they meant something."""
    text = _rows("P2\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\n")
    got = parse_annotations(text)
    assert got["P2"]["cog_category"] == ""
    assert got["P2"]["kegg_pathways"] == []
    assert got["P2"]["preferred_name"] == ""


def test_comment_lines_and_the_trailing_summary_are_skipped():
    """emapper writes '##' banner lines at the top and a '## Total time' line at the end.
    A parser that takes any non-header line produces a row keyed '## Total time'."""
    text = ("## emapper-2.1.12\n" + HEADER
            + "P3\t-\t-\t-\t-\t-\tS\tUncharacterised\t-\t-\t-\t-\t-\n"
            + "## Total time = 3 s\n")
    got = parse_annotations(text)
    assert list(got) == ["P3"]


def test_multiple_cog_categories_are_kept_whole():
    """Several single-letter COG categories are written adjacently, as 'EGP'. The table
    keeps them as one string; labels_from_orthology splits them."""
    text = _rows("P4\tx\t1e-5\t100\tCOG1\tBacteria\tEGP\tTransporter\t-\t-\t-\t-\t-\n")
    assert parse_annotations(text)["P4"]["cog_category"] == "EGP"


# --- every identifying column ----------------------------------------------------------

EMAPPER_FULL = "\n".join([
    "## emapper-2.1.12",
    "#query\tseed_ortholog\tevalue\tscore\teggNOG_OGs\tmax_annot_lvl\tCOG_category"
    "\tDescription\tPreferred_name\tGOs\tEC\tKEGG_ko\tKEGG_Pathway\tPFAMs",
    "p1\t83333.b0001\t1e-50\t200.0\tCOG5527@2,2QV1F@1224\t2|Bacteria\tL"
    "\tPlasmid replication initiator protein\trepA"
    "\tGO:0006270,GO:0003677\t2.7.7.7\tko:K02314\tko03030\tRepA_N,Bac_RepA_C",
    "p2\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-",
    "## Total time: 1 s",
])


def test_gene_symbol_pfams_go_ec_and_ko_are_all_parsed():
    """The gene symbol and the controlled-vocabulary columns, multi-valued ones as lists."""
    records = parse_annotations(EMAPPER_FULL)

    assert records["p1"]["preferred_name"] == "repA"
    assert records["p1"]["pfams"] == ["RepA_N", "Bac_RepA_C"]
    assert records["p1"]["gos"] == ["GO:0006270", "GO:0003677"]
    assert records["p1"]["ec"] == ["2.7.7.7"]
    assert records["p1"]["kegg_ko"] == ["ko:K02314"]


def test_the_placeholder_is_absence_in_every_list_column():
    """emapper writes a bare '-' for every field it has nothing for. Read literally that
    becomes a Pfam family named '-' and a GO term named '-', and both would aggregate as
    though they were real terms."""
    records = parse_annotations(EMAPPER_FULL)

    assert records["p2"]["pfams"] == []
    assert records["p2"]["gos"] == []
    assert records["p2"]["ec"] == []
    assert records["p2"]["kegg_ko"] == []
    assert records["p2"]["preferred_name"] == ""


def test_a_missing_column_yields_an_empty_list_not_a_failure():
    """A table without some columns parses, with those fields empty."""
    minimal = "\n".join([
        "#query\tCOG_category\tDescription",
        "p3\tL\tSome protein",
    ])

    records = parse_annotations(minimal)

    assert records["p3"]["pfams"] == []
    assert records["p3"]["gos"] == []
    assert records["p3"]["cog_category"] == "L"
