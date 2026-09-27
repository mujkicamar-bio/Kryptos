"""The label vocabulary: what a tool said, with the kind of statement it made."""
import pytest

from plasmidann import labels


def test_a_pfam_hit_yields_the_family_the_accession_the_description_and_the_clan():
    """One hmmsearch hit gives the family, its description and its clan, all with the
    accession."""
    pfam = {"RepA_N": {"accession": "PF06970.19",
                       "description": "Replication initiator protein A (RepA) N-terminus",
                       "type": "Domain", "clan": "CL0123"}}
    row = {"tier": "T1", "source": "pfam", "label": "RepA_N",
           "target_accession": "PF06970.19",
           "informative": "True"}

    out = labels.labels_from_hit(row, pfam=pfam)
    by_kind = {d["kind"]: d["label"] for d in out}

    assert by_kind["pfam_family"] == "RepA_N"
    assert by_kind["pfam_description"] == "Replication initiator protein A (RepA) N-terminus"
    assert by_kind["pfam_clan"] == "CL0123"
    assert all(d["accession"] == "PF06970.19" for d in out)


def test_a_pfam_family_with_no_clan_emits_no_clan_label():
    """An empty clan is absence. A label reading '' would group every clanless family
    together - 16,201 of the 30,134 in Pfam 38.2 - which is the largest false group it is
    possible to create here."""
    pfam = {"CagE_TrbE_VirB": {
        "accession": "PF03135.20",
        "description": "CagE, TrbE, VirB family, component of type IV transporter system",
        "type": "Family", "clan": ""}}
    row = {"tier": "T1", "source": "pfam", "label": "CagE_TrbE_VirB",
           "target_accession": "PF03135.20",
           "informative": "True"}

    kinds = {d["kind"] for d in labels.labels_from_hit(row, pfam=pfam)}

    assert "pfam_clan" not in kinds
    assert "pfam_family" in kinds


def test_a_pfam_family_absent_from_the_release_still_yields_the_family_name():
    """A family name the installed release does not know is a real event - a database
    mismatch - and the name must survive so the mismatch is visible in the table rather
    than silently dropping the hit."""
    row = {"tier": "T2", "source": "pfam", "label": "Not_In_Pfam",
           "target_accession": "PF99999.1",
           "informative": "True"}

    out = labels.labels_from_hit(row, pfam={})

    assert [d["kind"] for d in out] == ["pfam_family"]


def test_a_diamond_hit_yields_the_product_name_and_the_accession():
    """The NCBI title is 'ACC RecName: Full=Toxin CcdB; AltName: ... [organism]' for
    swissprot and 'ACC product name [organism]' for nr - verified against the installed
    database. The product name is the label; the organism is not a functional statement
    and is not emitted as one."""
    row = {"tier": "T3", "source": "swissprot", "informative": "True",
           "target_accession": "P62554.1",
           "label": "P62554.1 RecName: Full=Toxin CcdB; AltName: Full=Protein LetD "
                    "[Escherichia coli K-12]"}

    out = labels.labels_from_hit(row)

    assert {d["kind"]: d["label"] for d in out} == {"swissprot_product": "Toxin CcdB"}
    assert not any("Escherichia coli K-12" in d["label"] for d in out)


def test_a_qualifier_before_recname_is_kept_and_recname_is_still_parsed():
    """'Q47718.1 PUTATIVE PSEUDOGENE: RecName: Full=...', seen in the installed swissprot:
    the qualifier stays in front of the recommended name."""
    title = ("Q47718.1 PUTATIVE PSEUDOGENE: RecName: Full=Putative transposase InsO for "
             "insertion sequence element IS911B [Escherichia coli K-12]")

    parsed = labels.parse_ncbi_title(title)

    assert parsed["product"] == ("PUTATIVE PSEUDOGENE: Putative transposase InsO for "
                                 "insertion sequence element IS911B")
    assert parsed["accession"] == "Q47718.1"


def test_an_nr_hit_yields_the_product_name_of_the_clusterednr_title():
    row = {"tier": "T4", "source": "nr", "informative": "True",
           "target_accession": "WP_000813620.1",
           "label": "WP_000813620.1 type II toxin-antitoxin system RelE/ParE family "
                    "toxin [Escherichia coli]"}

    by_kind = {d["kind"]: d["label"] for d in labels.labels_from_hit(row)}

    assert by_kind["nr_product"] == ("type II toxin-antitoxin system RelE/ParE family "
                                     "toxin")


def test_an_uninformative_hit_yields_no_label():
    """'hypothetical protein' names nothing. It is evidence that someone has seen the
    protein, which the cascade already records as dark evidence, and it must never enter
    the label table."""
    row = {"tier": "T4", "source": "nr", "informative": "False",
           "target_accession": "WP_2.1",
           "label": "WP_2.1 hypothetical protein [Escherichia coli]"}

    assert labels.labels_from_hit(row) == []


def test_an_unparsable_title_still_yields_the_whole_title_as_the_product():
    """nr titles are not uniform and this parser will meet shapes it was not shown. The
    honest fallback is the whole string: a dropped label is invisible, a strange label is
    not."""
    row = {"tier": "T4", "source": "nr", "informative": "True", "target_accession": "",
           "label": "something with no recognisable structure"}

    by_kind = {d["kind"]: d["label"] for d in labels.labels_from_hit(row)}

    assert by_kind["nr_product"] == "something with no recognisable structure"


def test_orthology_yields_the_gene_symbol_and_every_controlled_identifier():
    row = {"seq_id": "p1", "cog_category": "L", "preferred_name": "repA",
           "eggnog_description": "Plasmid replication initiator protein",
           "pfams": "RepA_N,Bac_RepA_C", "gos": "GO:0006270", "ec": "2.7.7.7",
           "kegg_ko": "ko:K02314", "kegg_pathways": "ko03030", "eggnog_ogs": "COG5527@2"}

    pairs = {(d["kind"], d["label"]) for d in labels.labels_from_orthology(row)}

    assert ("gene_symbol", "repA") in pairs
    assert ("cog_category", "L") in pairs
    assert ("cog_id", "COG5527") in pairs
    assert ("eggnog_pfam", "RepA_N") in pairs
    assert ("eggnog_pfam", "Bac_RepA_C") in pairs
    assert ("go", "GO:0006270") in pairs
    assert ("ec", "2.7.7.7") in pairs
    assert ("kegg_ko", "ko:K02314") in pairs


def test_a_multi_letter_cog_category_is_split_into_one_label_per_letter():
    """COG categories are written adjacently as 'LKV', and each letter is a separate
    category. Kept whole it is a category named 'LKV' that no COG release defines; split,
    it is three real ones. The orthology TABLE keeps the string whole for provenance; the
    VOCABULARY needs the individual categories."""
    row = {"seq_id": "p1", "cog_category": "LKV", "preferred_name": "",
           "eggnog_description": "", "pfams": "", "gos": "", "ec": "",
           "kegg_ko": "", "kegg_pathways": "", "eggnog_ogs": ""}

    cats = {d["label"] for d in labels.labels_from_orthology(row)
            if d["kind"] == "cog_category"}

    assert cats == {"L", "K", "V"}


def test_an_empty_orthology_row_yields_no_labels():
    row = {"seq_id": "p2", "cog_category": "", "preferred_name": "",
           "eggnog_description": "", "pfams": "", "gos": "", "ec": "",
           "kegg_ko": "", "kegg_pathways": "", "eggnog_ogs": ""}

    assert labels.labels_from_orthology(row) == []


def test_every_emitted_kind_is_declared():
    """protein_labels refuses a kind not in KINDS, so every kind emitted must be there."""
    pfam = {"RepA_N": {"accession": "PF1.1", "description": "d", "type": "Family",
                       "clan": "CL0123"}}
    emitted = set()
    emitted |= {d["kind"] for d in labels.labels_from_hit(
        {"tier": "T1", "source": "pfam", "label": "RepA_N", "target_accession": "PF1.1",
         "informative": "True"}, pfam=pfam)}
    emitted |= {d["kind"] for d in labels.labels_from_hit(
        {"tier": "T3", "source": "swissprot",
         "label": "P1.1 RecName: Full=Toxin CcdB [E. coli]",
         "target_accession": "P1.1", "informative": "True"})}
    emitted |= {d["kind"] for d in labels.labels_from_hit(
        {"tier": "T4", "source": "nr", "label": "W1.1 relaxase [E. coli]",
         "target_accession": "W1.1",
         "informative": "True"})}
    emitted |= {d["kind"] for d in labels.labels_from_orthology(
        {"seq_id": "p", "cog_category": "L", "preferred_name": "repA",
         "eggnog_description": "x", "pfams": "A", "gos": "GO:1", "ec": "1.1.1.1",
         "kegg_ko": "ko:K1", "kegg_pathways": "", "eggnog_ogs": "COG1@2"})}

    assert emitted <= labels.KINDS, f"undeclared kinds: {sorted(emitted - labels.KINDS)}"


def test_a_swissprot_hit_is_recognised_by_its_source_not_its_tier_id():
    row = {"tier": "T4", "source": "swissprot", "informative": "True",
           "target_accession": "P62554.1",
           "label": "P62554.1 RecName: Full=Toxin CcdB; AltName: Full=Protein LetD "
                    "[Escherichia coli K-12]"}

    by_kind = {d["kind"]: d["label"] for d in labels.labels_from_hit(row)}

    assert by_kind == {"swissprot_product": "Toxin CcdB"}, (
        "Swiss-Prot at T4 was read as nr: the label kind is being guessed from tier "
        "position instead of read from the row's source")


def test_a_pharokka_hit_yields_the_annotation_and_the_category():
    """Two labels from one hit, because they are two different statements: what the family
    is called, and which functional group the database puts it in. The category is the
    grouping PHROGs' curators made, carried as a label like any other."""
    row = {"tier": "T3", "source": "pharokka", "informative": "True",
           "target_accession": "phrog_164", "label": "ParA-like partition protein",
           "category": "DNA, RNA and nucleotide metabolism"}

    by_kind = {d["kind"]: d["label"] for d in labels.labels_from_hit(row)}

    assert by_kind == {"pharokka_annotation": "ParA-like partition protein",
                       "pharokka_category": "DNA, RNA and nucleotide metabolism"}


def test_card_and_vfdb_hits_from_the_pharokka_tier_have_their_own_kinds():
    card = {"tier": "T3", "source": "card", "informative": "True",
            "target_accession": "ARO:3000873", "label": "TEM beta-lactamase",
            "category": "antibiotic inactivation"}
    vfdb = {"tier": "T3", "source": "vfdb", "informative": "True",
            "target_accession": "VFG000001", "label": "type IV pilus", "category": ""}

    assert {d["kind"] for d in labels.labels_from_hit(card)} == {
        "card_gene_family", "card_mechanism"}
    assert {d["kind"] for d in labels.labels_from_hit(vfdb)} == {"vfdb_factor"}


def test_a_row_without_a_source_is_refused():
    """The source decides the label kind, so a row without one is refused rather than
    guessed from its tier id."""
    row = {"tier": "T1", "informative": "True", "target_accession": "PF00000",
           "label": "RepA_N"}
    with pytest.raises(ValueError):
        labels.labels_from_hit(row)


def test_ncbi_record_prefixes_stay_in_the_product_as_ncbi_writes_them():
    for title, product in (
            ("MBD3193859.1 MAG: ABC transporter permease [Clostridia bacterium]",
             "MAG: ABC transporter permease"),
            ("DAD1.1 TPA_asm: MobA/MobL family protein [Siphoviridae sp.]",
             "TPA_asm: MobA/MobL family protein"),
            ("WP_1.1 MULTISPECIES: relaxase [Bacillus]", "MULTISPECIES: relaxase")):
        assert labels.parse_ncbi_title(title)["product"] == product


def test_the_plasmid_label_database_kinds_are_declared():
    """protein_labels refuses an undeclared kind, so every kind the plasmid label
    databases emit must be in KINDS."""
    from plasmidann import labeldb

    assert labeldb.KINDS <= labels.KINDS
    # pharokka's CARD kinds and the direct CARD search stay distinct statements.
    assert {"card_gene_family", "card_amr_family"} <= labels.KINDS
