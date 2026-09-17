"""Collapsing free-text product names to a canonical form.

Pfam family names and COG identifiers are controlled and finite, so they can be enumerated.
Product names from nr and Swiss-Prot are neither: the same protein is written 'Plasmid
replication initiator protein RepA', 'replication initiation protein' and 'putative
replication initiator protein' by different submitters. Left as written, one concept
becomes dozens of categories, each too rare to test.

Normalisation is what makes the frequency ranking meaningful: once the surface forms
collapse, the head of the distribution is a bounded list a person can review.
"""
from plasmidann import normalise


def test_the_same_concept_written_three_ways_collapses_to_one_string():
    forms = ["Plasmid replication initiator protein",
             "plasmid replication initiator protein",
             "putative plasmid replication initiator protein"]

    assert len({normalise.product_name(f) for f in forms}) == 1


def test_hedging_words_are_stripped():
    """'putative', 'probable' and 'predicted' describe the submitter's confidence, not the
    protein. Keeping them would split every product into a confident and a hedged form."""
    assert normalise.product_name("putative relaxase") == "relaxase"
    assert normalise.product_name("probable relaxase") == "relaxase"
    assert normalise.product_name("predicted relaxase") == "relaxase"


def test_the_generic_protein_suffixes_are_stripped():
    """'MobA/MobL family protein' and 'MobA/MobL protein' are one product. The suffixes
    carry no functional information and are applied inconsistently by submitters."""
    assert normalise.product_name("MobA/MobL family protein") == "moba/mobl"
    assert normalise.product_name("MobA/MobL protein") == "moba/mobl"
    assert normalise.product_name("DUF1234 domain-containing protein") == "duf1234"


def test_whitespace_and_case_are_canonicalised():
    assert normalise.product_name("  Conjugal   Transfer  Protein  TraG ") == \
        "conjugal transfer protein trag"


def test_a_name_that_is_only_a_suffix_survives_whole():
    """Stripping must not empty a product. 'protein' alone is a real, if useless, product
    name, and returning '' would merge it with every parse failure in the collection."""
    assert normalise.product_name("protein") == "protein"
    assert normalise.product_name("hypothetical protein") == "hypothetical"


def test_an_empty_input_gives_an_empty_string():
    assert normalise.product_name("") == ""
    assert normalise.product_name(None) == ""
