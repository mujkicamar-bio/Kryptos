"""Tests for the plasmid-specific label databases (plasmidann.labeldb) and their stage script.

Each rule is tested at its boundary, because every one of them is a threshold that a
methods section quotes: the identity and coverage tiers on BOTH sequences, CARD's Perfect
and Strict calls including a hit just below the curated cut-off, the AMRFinderPlus element
types.
"""
import json
import stat

import pytest
from conftest import FakeSnakemake, read_tsv, requires, run_script, write_fasta, write_tsv

from plasmidann import labeldb


# ------------------------------------------------------------------------------------
# Tiers (Islam et al. 2026): identity and coverage, coverage on query AND subject
# ------------------------------------------------------------------------------------
def test_tier_one_at_its_inclusive_boundaries():
    assert labeldb.tier(80.0, 90.0, 90.0) == 1


def test_tier_two_is_strictly_above_its_boundaries():
    assert labeldb.tier(79.9, 95.0, 95.0) == 2
    assert labeldb.tier(60.1, 70.1, 70.1) == 2
    assert labeldb.tier(60.0, 95.0, 95.0) is None
    assert labeldb.tier(95.0, 70.0, 95.0) is None


def test_tier_requires_coverage_on_the_subject_as_well_as_the_query():
    # A short query aligned in full to a long reference: query coverage alone would pass.
    assert labeldb.tier(99.0, 100.0, 40.0) is None
    # A long query carrying the whole reference as one domain: subject coverage alone.
    assert labeldb.tier(99.0, 40.0, 100.0) is None
    # Tier 1 on the query side but only tier 2 on the subject side is tier 2.
    assert labeldb.tier(95.0, 95.0, 85.0) == 2


def _hit(query, key, pident, qcov, scov, bitscore, **extra):
    return {"query": query, "key": key, "pident": pident, "qcov": qcov, "scov": scov,
            "bitscore": bitscore, **extra}


def test_best_tiered_hit_prefers_the_better_tier_over_the_higher_bitscore():
    hits = [_hit("q1", "1", 70.0, 80.0, 80.0, 900.0),   # tier 2, higher score
            _hit("q1", "2", 85.0, 95.0, 95.0, 300.0),   # tier 1
            _hit("q1", "3", 99.0, 95.0, 30.0, 950.0)]   # fails subject coverage
    best = labeldb.best_tiered_hits(hits)
    assert best["q1"][0] == 1 and best["q1"][1]["key"] == "2"


def test_best_tiered_hit_breaks_ties_by_bitscore_then_the_smaller_subject_number():
    hits = [_hit("q1", "10", 85.0, 95.0, 95.0, 300.0),
            _hit("q1", "2", 85.0, 95.0, 95.0, 300.0),
            _hit("q1", "3", 85.0, 95.0, 95.0, 310.0),
            _hit("q2", "26", 50.0, 95.0, 95.0, 999.0)]
    best = labeldb.best_tiered_hits(hits)
    assert best["q1"][1]["key"] == "3"
    assert "q2" not in best          # below both tiers: no label at all
    assert labeldb.best_tiered_hits(hits[:2])["q1"][1]["key"] == "2"
    assert labeldb.best_tiered_hits(hits[1::-1])["q1"][1]["key"] == "2"


# ------------------------------------------------------------------------------------
# CARD (Alcock et al. 2023): Perfect and Strict on protein homolog models only
# ------------------------------------------------------------------------------------
def test_card_perfect_is_full_identity_over_the_whole_reference():
    assert labeldb.card_call(100.0, 1, 286, 286, 580.0, 500.0) == "Perfect"
    # Perfect does not depend on the cut-off: the reference itself is resistance.
    assert labeldb.card_call(100.0, 1, 286, 286, 90.0, 500.0) == "Perfect"


def test_card_identity_over_part_of_the_reference_is_not_perfect():
    # 100% identical but starting at reference residue 5: Strict when above the cut-off.
    assert labeldb.card_call(100.0, 5, 286, 286, 560.0, 500.0) == "Strict"
    assert labeldb.card_call(100.0, 1, 280, 286, 560.0, 500.0) == "Strict"
    assert labeldb.card_call(99.6, 1, 286, 286, 559.0, 500.0) == "Strict"


def test_card_strict_is_bitscore_at_or_above_the_cut_off_and_loose_is_rejected():
    assert labeldb.card_call(90.0, 1, 286, 286, 500.0, 500.0) == "Strict"
    assert labeldb.card_call(90.0, 1, 286, 286, 499.9, 500.0) is None


def _card_json(path):
    def model(mid, name, mtype, seq, cutoff, families):
        cats = {str(i): {"category_aro_class_name": "AMR Gene Family",
                         "category_aro_name": f} for i, f in enumerate(families)}
        cats["99"] = {"category_aro_class_name": "Drug Class",
                      "category_aro_name": "cephalosporin"}
        entry = {"model_id": mid, "model_name": name, "model_type": mtype,
                 "ARO_accession": "300" + mid, "ARO_category": cats,
                 "model_sequences": {"sequence": {"1": {"protein_sequence": {
                     "accession": "X", "sequence": seq}}}}}
        if cutoff is not None:
            entry["model_param"] = {"blastp_bit_score": {"param_value": str(cutoff)}}
        return entry
    data = {
        "_version": "4.0.2", "_timestamp": "2026-08-11",
        "1": model("1", "CTX-M-15", "protein homolog model", "mvkkslrqft", 500,
                   ["CTX-M beta-lactamase"]),
        "2": model("2", "vanY gene in vanM cluster", "protein homolog model", "MKKLLW", 200,
                   ["glycopeptide resistance gene cluster", "vanY"]),
        "3": model("3", "Escherichia coli gyrA", "protein variant model", "MSDLAREIT", 800,
                   ["fluoroquinolone resistant gyrA"]),
        "4": model("4", "tetR", "protein overexpression model", "MARLNRE", 300,
                   ["major facilitator superfamily (MFS) antibiotic efflux pump"]),
    }
    path.write_text(json.dumps(data))
    return path


def test_card_models_keep_protein_homolog_models_only_with_their_cut_off(tmp_path):
    models = labeldb.card_models(_card_json(tmp_path / "card.json"))
    assert set(models) == {"1", "2"}
    assert models["1"]["cut_off"] == 500.0
    assert models["1"]["sequence"] == "MVKKSLRQFT"
    assert models["1"]["subject"] == "ARO:3001"
    assert models["1"]["labels"] == [("CTX-M beta-lactamase", "CTX-M-15")]
    # A model in several AMR gene families gives one label per family.
    assert models["2"]["labels"] == [("glycopeptide resistance gene cluster",
                                      "vanY gene in vanM cluster"),
                                     ("vanY", "vanY gene in vanM cluster")]


def test_card_best_hit_prefers_perfect_then_the_highest_strict_bitscore(tmp_path):
    models = labeldb.card_models(_card_json(tmp_path / "card.json"))
    hits = [
        # q1: a Strict hit with a higher score and a Perfect hit - Perfect wins.
        dict(query="q1", key="2", pident=90.0, sstart=1, send=6, slen=6, bitscore=900.0),
        dict(query="q1", key="1", pident=100.0, sstart=1, send=10, slen=10, bitscore=550.0),
        # q2: two Strict hits, the higher bitscore wins.
        dict(query="q2", key="1", pident=95.0, sstart=1, send=10, slen=10, bitscore=520.0),
        dict(query="q2", key="2", pident=95.0, sstart=1, send=6, slen=6, bitscore=600.0),
        # q3: Loose only (below the cut-off) - no call.
        dict(query="q3", key="1", pident=90.0, sstart=1, send=10, slen=10, bitscore=499.0),
        # q4: a hit to a model that is not a protein homolog model is never read.
        dict(query="q4", key="3", pident=100.0, sstart=1, send=9, slen=9, bitscore=999.0),
    ]
    best = labeldb.best_card_hits(hits, models)
    assert best["q1"][0] == "Perfect" and best["q1"][1]["key"] == "1"
    assert best["q2"][0] == "Strict" and best["q2"][1]["key"] == "2"
    assert "q3" not in best and "q4" not in best


# ------------------------------------------------------------------------------------
# Header and metadata parsers: (label, sub_label) per database entry
# ------------------------------------------------------------------------------------
def test_tadb_role_and_type_come_from_the_file_stem_kept_in_the_id():
    header = ("T28|type_II_T T28 WP_000916169.1 NC_000915:c946611-946345 [Helicobacter "
              "pylori 26695]")
    assert labeldb.parse_tadb(header) == ("T28", [("type II toxin", "")])
    assert labeldb.parse_tadb("AT1042|type_IV_AT AT1042 NP_416508.1 x [E. coli]") == \
        ("AT1042", [("type IV antitoxin", "")])
    assert labeldb.parse_tadb("RE6024|regulator RE6024 WP_1.1 x [S. pyogenes]") == \
        ("RE6024", [("regulator", "")])
    with pytest.raises(ValueError):
        labeldb.parse_tadb("T1|tadb T1 WP_1.1")


def test_bacmet_labels_are_the_compounds_verbatim_with_gene_and_class_as_sub_label():
    header = ("BAC0224|merA|sp|P00392|MERA_PSEAI Mercuric reductase OS=Pseudomonas "
              "aeruginosa GN=merA")
    compounds = {"BAC0224": "Mercury (Hg), Organo-mercury compounds [class: Organo-mercury]",
                 "BAC0001": "Triclosan [class: Phenolic compounds], Acriflavine [class: "
                            "Acridine], Proflavine [class: Acridine], Copper (Cu), "
                            "Wex-cide-128."}
    assert labeldb.parse_bacmet(header, compounds) == \
        ("BAC0224", [("Mercury (Hg)", "merA"),
                     ("Organo-mercury compounds", "merA; class=Organo-mercury")])
    # One label per compound, a class shared by two compounds included; the '.' that ends
    # the list is not part of the last compound.
    assert labeldb.parse_bacmet("BAC0001|abeM|tr|Q5FAM9|x", compounds) == \
        ("BAC0001", [("Acriflavine", "abeM; class=Acridine"), ("Copper (Cu)", "abeM"),
                     ("Proflavine", "abeM; class=Acridine"),
                     ("Triclosan", "abeM; class=Phenolic compounds"),
                     ("Wex-cide-128", "abeM")])


def test_bacmet_mapping_file_is_read_by_bacmet_id(tmp_path):
    path = tmp_path / "BacMet2_EXP.753.mapping.txt"
    write_tsv(path, ["BacMet_ID", "Gene_name", "Accession", "Organism", "Location",
                     "Compound"],
              [["BAC0224", "merA", "P00392", "P. aeruginosa", "Plasmid", "Mercury (Hg)"]])
    assert labeldb.read_bacmet_compounds(path) == {"BAC0224": "Mercury (Hg)"}


def test_oritdb_role_from_the_installed_id_and_family_from_the_original_header():
    assert labeldb.parse_oritdb(
        "relaxase_00001 TraI_RP4 CAA38336 MOBP id=1 [Escherichia coli HB101]") == \
        ("TraI_RP4", [("relaxase", "MOBP")])
    assert labeldb.parse_oritdb("t4cp_00001 TraG_RP4 CAA38334 VirD4/TraG id=1 [E. coli]") \
        == ("TraG_RP4", [("T4CP", "VirD4/TraG")])
    assert labeldb.parse_oritdb(
        "auxiliary_00001 TraJ_RP4 CAA38338 id=1 [Escherichia coli HB101]") == \
        ("TraJ_RP4", [("auxiliary protein", "")])
    # The token before 'id=' is the family, whatever precedes it; '_' is no family.
    assert labeldb.parse_oritdb("relaxase_00009 MobA_x WP_1 WP_000119405 MOBV id=9 [x]")[1] \
        == [("relaxase", "MOBV")]
    assert labeldb.parse_oritdb("relaxase_00010 MobA_y WP_2 _ id=9 [x]")[1] == \
        [("relaxase", "")]
    with pytest.raises(ValueError):
        labeldb.parse_oritdb("TraI_RP4 CAA38336 MOBP id=1")


def test_mobileog_major_category_with_minor_category_and_evidence():
    assert labeldb.parse_mobileog(
        "mobileOG_000001713|aq_aa05|O66401|integration/excision|NA|Multiple|Manual") == \
        ("mobileOG_000001713", [("integration/excision", "evidence=Manual")])
    assert labeldb.parse_mobileog(
        "mobileOG_000000002|cII_1|A0A0J8XRP9|phage|infection,regulation|Multiple|Homology") \
        == ("mobileOG_000000002", [("phage", "infection,regulation; evidence=Homology")])
    assert labeldb.parse_mobileog(
        "mobileOG_000278000|cI|WP_154180136.1|phage|replication,transfer|Plasmid RefSeq|"
        "Homology") == ("mobileOG_000278000",
                        [("phage", "replication,transfer; evidence=Homology")])


def test_mobileog_keyword_search_entries_are_kept_with_their_evidence():
    assert labeldb.parse_mobileog(
        "mobileOG_000500000|x|Q1|transfer|NA|Multiple|Keyword Search") == \
        ("mobileOG_000500000", [("transfer", "evidence=Keyword Search")])


def test_dbapis_family_with_gene_and_evidence():
    assert labeldb.parse_dbapis(
        "YP_009986712.1 gene=U56 family=APIS125 evidence=verified hypothetical protein "
        "JR326_gp056 [Escherichia phage ukendt]") == \
        ("YP_009986712.1", [("APIS125", "U56; evidence=verified")])
    assert labeldb.parse_dbapis("WP_1.1 gene=gp54 family=gp54 evidence=homolog")[1] == \
        [("gp54", "gp54; evidence=homolog")]
    # gene=NA is the installer's mark for a family without a verified seed: no gene.
    assert labeldb.parse_dbapis("WP_2.1 gene=NA family=APIS030 evidence=homolog")[1] == \
        [("APIS030", "evidence=homolog")]
    # An accession in two families gives one label per family.
    assert labeldb.parse_dbapis("MGV_7 gene=ArdA family=APIS003,APIS067 evidence=homolog")[1] \
        == [("APIS003", "ArdA; evidence=homolog"), ("APIS067", "ArdA; evidence=homolog")]
    with pytest.raises(ValueError):
        labeldb.parse_dbapis("WP_1.1 gene=gp54 family=gp54")


def test_acrdb_family_with_crispr_type_and_evidence():
    assert labeldb.parse_acrdb("anti_CRISPR0001 family=AcrIF1 type=I-F acc=YP_007392342.1 "
                               "evidence=Verified") == \
        ("anti_CRISPR0001", [("AcrIF1", "I-F; evidence=Verified")])
    assert labeldb.parse_acrdb("anti_CRISPR0900 family=AcrIIA11a.1 type=V-A,I-C acc=X "
                               "evidence=Putative")[1] == \
        [("AcrIIA11a.1", "V-A,I-C; evidence=Putative")]
    assert labeldb.parse_acrdb("anti_CRISPR0901 family=AcrIF2 type= acc=X "
                               "evidence=PLiterature")[1] == \
        [("AcrIF2", "evidence=PLiterature")]
    with pytest.raises(ValueError):
        labeldb.parse_acrdb("anti_CRISPR0002 family= type=I-F acc=X evidence=Putative")
    with pytest.raises(ValueError):
        labeldb.parse_acrdb("anti_CRISPR0003 family=AcrIF1 type=I-F acc=X")


def test_fasta_records_keep_the_whole_header_and_drop_whitespace_in_sequences(tmp_path):
    path = tmp_path / "x.faa"
    path.write_text(">MobA_x CAA1 MOBQ id=2 [B. sp.]\nMATK\nLLQ*\n>B y\nmk lw\n")
    assert list(labeldb.fasta_records(path)) == [("MobA_x CAA1 MOBQ id=2 [B. sp.]",
                                                  "MATKLLQ"), ("B y", "MKLW")]


def test_load_reference_numbers_entries_in_file_order(tmp_path):
    d = tmp_path / "tadb"
    write_fasta(d / "tadb.faa", [("AT1|type_II_AT AT1 WP_3.1 x [a]", "MSTQRE"),
                                 ("T1|type_II_T T1 WP_1.1 x [a]", "MKKLLVVLAAG")])
    entries = labeldb.load_reference("tadb", d)
    assert [(e["key"], e["subject"]) for e in entries] == [("0", "AT1"), ("1", "T1")]
    assert entries[1]["labels"] == [("type II toxin", "")]
    assert entries[1]["sequence"] == "MKKLLVVLAAG"


def test_load_reference_reads_the_bacmet_mapping_under_raw(tmp_path):
    d = tmp_path / "bacmet"
    write_fasta(d / "bacmet.faa", [("BAC0224|merA|sp|P00392|MERA x", "MTHLKITGMTCDS")])
    write_tsv(d / "raw" / "BacMet2_EXP.753.mapping.txt",
              ["BacMet_ID", "Gene_name", "Accession", "Organism", "Location", "Compound"],
              [["BAC0224", "merA", "P00392", "x", "Plasmid", "Mercury (Hg)"]])
    (entry,) = labeldb.load_reference("bacmet", d)
    assert entry["labels"] == [("Mercury (Hg)", "merA")]


def test_load_reference_keeps_every_mobileog_evidence_class(tmp_path):
    d = tmp_path / "mobileog"
    write_fasta(d / "mobileog.faa", [
        ("mobileOG_1|intI1|Q1|integration/excision|NA|Multiple|Manual", "MKTLLAAG"),
        ("mobileOG_2|x|Q2|transfer|NA|Multiple|Keyword Search", "MSTQREWW")])
    entries = labeldb.load_reference("mobileog", d)
    assert [(e["subject"], e["labels"][0][1]) for e in entries] == \
        [("mobileOG_1", "evidence=Manual"), ("mobileOG_2", "evidence=Keyword Search")]


# ------------------------------------------------------------------------------------
# AMRFinderPlus (Feldgarden et al. 2021): protein mode, --plus
# ------------------------------------------------------------------------------------
AMR_V4 = ["Protein id", "Contig id", "Start", "Stop", "Strand", "Element symbol",
          "Element name", "Scope", "Type", "Subtype", "Class", "Subclass", "Method",
          "Target length", "Reference sequence length", "% Coverage of reference",
          "% Identity to reference", "Alignment length", "Closest reference accession",
          "Closest reference name", "HMM accession", "HMM description"]


def _amr_row(pid, symbol, etype, subtype, method="EXACTP", cov="100.00", ident="100.00",
             acc="WP_000027057.1"):
    return [pid, "NA", "NA", "NA", "NA", symbol, "name", "core", etype, subtype, "CLASS",
            "SUBCLASS", method, "286", "286", cov, ident, "286", acc, "ref", "NA", "NA"]


def test_amrfinder_elements_become_labels_with_their_type_as_sub_label(tmp_path):
    path = tmp_path / "amrfinder.tsv"
    write_tsv(path, AMR_V4, [
        _amr_row("p1", "blaTEM-1", "AMR", "AMR"),
        _amr_row("p2", "merA", "STRESS", "METAL", method="BLASTP", ident="97.20"),
        _amr_row("p3", "qacEdelta1", "STRESS", "BIOCIDE"),
        _amr_row("p4", "asr", "STRESS", "ACID", method="HMM", cov="NA", ident="NA",
                 acc="NA"),
        _amr_row("p5", "iutA", "VIRULENCE", "VIRULENCE"),
    ])
    rows = {r["seq_id"]: r for r in labeldb.parse_amrfinder(path, "2026-08-01.1")}
    assert set(rows) == {"p1", "p2", "p3", "p4", "p5"}
    for r in rows.values():
        assert r["source"] == "amrfinder" and r["label_kind"] == "amrfinder_gene"
        assert r["database_version"] == "2026-08-01.1"
        assert list(r) == labeldb.COLUMNS
    assert rows["p1"]["label"] == "blaTEM-1" and rows["p1"]["sub_label"] == "AMR/AMR"
    assert rows["p2"]["sub_label"] == "STRESS/METAL" and rows["p2"]["tier"] == "BLASTP"
    assert rows["p2"]["pident"] == "97.20" and rows["p2"]["scov"] == "100.00"
    assert rows["p4"]["pident"] == "" and rows["p4"]["subject"] == ""


def test_amrfinder_unknown_element_type_is_refused(tmp_path):
    path = tmp_path / "amrfinder.tsv"
    write_tsv(path, AMR_V4, [_amr_row("p1", "x", "NEWTYPE", "NEWTYPE")])
    with pytest.raises(ValueError):
        labeldb.parse_amrfinder(path, "v")


def test_label_rows_carry_the_label_columns():
    entries = {"7": {"subject": "BAC0224", "labels": [("Mercury", "merA")], "cut_off": ""}}
    best = {"q1": (2, _hit("q1", "7", 75.5, 88.0, 91.0, 250.0))}
    (row,) = labeldb.label_rows("bacmet", best, entries, "2.0")
    assert list(row) == labeldb.COLUMNS
    assert row == {"seq_id": "q1", "source": "bacmet", "label_kind": "bacmet_compound",
                   "label": "Mercury", "sub_label": "merA", "tier": "2", "cut_off": "",
                   "pident": "75.5", "qcov": "88.0", "scov": "91.0", "bitscore": "250.0",
                   "subject": "BAC0224", "database_version": "2.0"}


# ------------------------------------------------------------------------------------
# The stage script, outside Snakemake
# ------------------------------------------------------------------------------------
# Two unrelated 120-residue proteins, so a DIAMOND search can tell them apart.
PROT_A = ("MSKIAVLGAGGWGTALAIVLAKRGHEVRLWGRDPEKIAALRAEGENRRYLPGVKLPEGLTLTADLAEALAGAD"
          "LILVAVPSQALRELLPQLAPHLKAGATVVSLAKGIEPGTLKLLSEIV")
PROT_B = ("MTTQRDYLVLGAGSAGLAAAWHLRQAGHRVTVLEAQDRVGGRCWTRDLDGAPVELGAQWIHGGDNALYRALCEE"
          "LGLPLRDRTRPWLRLHDGSGWHSVDEAALRAALEAALAAGRDTGLS")


def _refs(root, dbs=("tadb", "bacmet"), card=False):
    for db in dbs:
        (root / db).mkdir(parents=True)
        (root / db / "VERSION").write_text(f"{db}-v1\n")
        (root / db / "SOURCE").write_text("test fixture\n")
    if "tadb" in dbs:
        write_fasta(root / "tadb" / "tadb.faa", [("T28|type_II_T T28 WP_1.1 x [a]", PROT_A)])
    if "bacmet" in dbs:
        write_fasta(root / "bacmet" / "bacmet.faa",
                    [("BAC0224|merA|sp|P00392|MERA x", PROT_B)])
        write_tsv(root / "bacmet" / "raw" / "BacMet2_EXP.753.mapping.txt",
                  ["BacMet_ID", "Gene_name", "Accession", "Organism", "Location",
                   "Compound"],
                  [["BAC0224", "merA", "P00392", "x", "Plasmid", "Mercury (Hg)"]])
    if card:
        (root / "card").mkdir()
        (root / "card" / "VERSION").write_text("4.0.2\n")
        (root / "card" / "SOURCE").write_text("x\n")
        write_fasta(root / "card" / "protein_fasta_protein_homolog_model.fasta",
                    [("gb|X|ARO:3001|CTX-M-15 [x]", PROT_A)])
        data = json.loads(_card_json(root / "card" / "card.json").read_text())
        data["1"]["model_sequences"]["sequence"]["1"]["protein_sequence"]["sequence"] = PROT_A
        (root / "card" / "card.json").write_text(json.dumps(data))


def _fake_amrfinder(tmp_path, rows):
    """An executable standing in for amrfinder: checks its arguments, writes a table."""
    table = tmp_path / "canned.tsv"
    write_tsv(table, AMR_V4, rows)
    exe = tmp_path / "bin" / "amrfinder"
    exe.parent.mkdir()
    exe.write_text(
        "#!/bin/sh\n"
        f"echo \"$@\" > {tmp_path}/amrfinder_args\n"
        "while [ $# -gt 0 ]; do [ \"$1\" = -o ] && out=$2; shift; done\n"
        f"cp {table} \"$out\"\n")
    exe.chmod(exe.stat().st_mode | stat.S_IEXEC)
    db = tmp_path / "amrdb" / "2026-08-01.1"
    db.mkdir(parents=True)
    (db / "version.txt").write_text("2026-08-01.1\n")
    return str(exe), str(db)


def _snake(tmp_path, refs, amr, required=False, amr_required=False):
    faa = tmp_path / "unique_proteins.faa"
    write_fasta(faa, [("q1", PROT_A), ("q2", PROT_B), ("q3", PROT_A[:40])])
    out = tmp_path / "08_protein_labels"
    return FakeSnakemake(
        input={"faa": str(faa)},
        output={"tsv": str(out / "protein_labels_plasmid.tsv"),
                "status": str(out / "label_databases_status.tsv")},
        params={"labels": {"dir": str(refs), "required": required},
                "amrfinder": {"executable": amr[0], "database": amr[1],
                              "required": amr_required}},
        threads=2)


@requires("diamond")
def test_script_labels_proteins_and_records_absent_databases_as_not_run(tmp_path):
    refs = tmp_path / "refs"
    _refs(refs, card=True)
    amr = _fake_amrfinder(tmp_path, [_amr_row("q2", "merA", "STRESS", "METAL")])
    snake = _snake(tmp_path, refs, amr)
    run_script("label_databases.py", snake)
    rows = read_tsv(snake.output.tsv)
    got = {(r["seq_id"], r["source"], r["label"], r["tier"]) for r in rows}
    assert got == {("q1", "tadb", "type II toxin", "1"),
                   ("q2", "bacmet", "Mercury (Hg)", "1"),
                   ("q1", "card", "CTX-M beta-lactamase", "Perfect"),
                   ("q2", "amrfinder", "merA", "EXACTP")}
    # q3 is 40 residues of PROT_A: identical, but covering a third of the subject.
    assert not any(r["seq_id"] == "q3" for r in rows)
    assert all(list(r) == labeldb.COLUMNS for r in rows)
    card_row = next(r for r in rows if r["source"] == "card")
    assert card_row["cut_off"] == "500.0" and card_row["database_version"] == "4.0.2"
    assert next(r for r in rows if r["source"] == "tadb")["database_version"] == "tadb-v1"
    args = open(tmp_path / "amrfinder_args").read().split()
    assert args[:2] == ["-p", snake.input.faa] and "--plus" in args
    assert args[args.index("--database") + 1] == amr[1]
    status = {r["database"]: r for r in read_tsv(snake.output.status)}
    assert status["tadb"]["status"] == "SUCCESS" and status["tadb"]["n_proteins"] == "1"
    assert status["amrfinder"]["version"] == "2026-08-01.1"
    for db in ("oritdb", "mobileog", "dbapis", "acrdb"):
        assert status[db]["status"] == "NOT_RUN"


def test_script_halts_when_a_required_database_is_absent(tmp_path):
    refs = tmp_path / "refs"
    _refs(refs, dbs=("tadb",))
    amr = _fake_amrfinder(tmp_path, [])
    with pytest.raises(SystemExit) as err:
        run_script("label_databases.py", _snake(tmp_path, refs, amr, required=True))
    assert "bacmet" in str(err.value)


def test_script_halts_when_amrfinder_is_required_and_absent(tmp_path):
    refs = tmp_path / "refs"
    refs.mkdir()
    amr = (str(tmp_path / "no" / "amrfinder"), str(tmp_path / "nodb"))
    with pytest.raises(SystemExit) as err:
        run_script("label_databases.py", _snake(tmp_path, refs, amr, amr_required=True))
    assert "amrfinder" in str(err.value).lower()


def test_script_records_everything_not_run_when_nothing_is_installed(tmp_path):
    refs = tmp_path / "refs"
    refs.mkdir()
    amr = (str(tmp_path / "no" / "amrfinder"), str(tmp_path / "nodb"))
    snake = _snake(tmp_path, refs, amr)
    run_script("label_databases.py", snake)
    assert read_tsv(snake.output.tsv) == []
    assert open(snake.output.tsv).readline().rstrip("\n").split("\t") == labeldb.COLUMNS
    status = read_tsv(snake.output.status)
    assert {r["database"] for r in status} == set(labeldb.DATABASES) | {"amrfinder"}
    assert {r["status"] for r in status} == {"NOT_RUN"}
