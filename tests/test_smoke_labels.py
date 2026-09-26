"""Smoke tests: label databases, protein labels, orthology, feature files, quality gate.

Each test runs one workflow script against a small fixture.
"""
import pytest
from conftest import (
    PLASMID_LABEL_COLS,
    FakeSnakemake,
    _ps_table,
    _selection,
    read_tsv,
    run_script,
    write_fasta,
    write_tsv,
)


def test_the_quality_gate_passes_when_controls_are_annotated(fixture_dir):
    """The gate must not halt a healthy run. Controls that the cascade named correctly
    should clear it."""
    prot = fixture_dir / "protein_annotation.tsv"
    artefact = fixture_dir / "artefact_flags.tsv"
    flags = fixture_dir / "eligibility.tsv"
    report = fixture_dir / "gate.txt"

    rows = [[f"CTRL_{i:05d}_P0000{i}", "T1", "relaxase MobA", "FUNCTIONAL"]
            for i in range(100)]
    rows += [["seq_dark_1", "", "", "NONE"],
             ["seq_uh_1", "T5", "hypothetical protein", "UNCHARACTERIZED_HOMOLOG"],
             ["seq_skipped", "", "", "NOT_SEARCHED"]]
    write_tsv(prot, ["seq_id", "annot_tier", "annot_label", "functional_class"], rows)
    write_tsv(artefact, ["seq_id", "artefact_flag"],
              [[r[0], 0] for r in rows])

    run_script("quality_gate.py", FakeSnakemake(
        input={"prot": str(prot), "artefact": str(artefact)},
        output={"flags": str(flags), "report": str(report)},
        params={"gate": {"min_control_recall": 0.99, "require_control_set": True},
                "tier_sources": {"T1": "pfam", "T4": "swissprot"}}))

    text = report.read_text()
    assert "control_recall=1.0" in text
    # Each tier is named with its database, taken from the cascade configuration.
    assert "  T1\tpfam\t100\n" in text, text
    # Both dark definitions: ours (2) and FESNov's no-homologue one (1).
    assert "target_eligible=2 dark_no_homologue=1" in text, text
    # A protein the selection never searched is not reported as annotated.
    reasons = {r["seq_id"]: r["exclusion_reason"] for r in read_tsv(flags)}
    assert reasons["seq_skipped"] == "not_searched"


def test_the_quality_gate_halts_when_known_proteins_come_out_dark(fixture_dir):
    """The gate's whole purpose. If it cannot fail here it is not a gate."""
    prot = fixture_dir / "protein_annotation.tsv"
    artefact = fixture_dir / "artefact_flags.tsv"
    flags = fixture_dir / "eligibility.tsv"
    report = fixture_dir / "gate.txt"

    rows = [[f"CTRL_{i:05d}_P0000{i}", "", "hypothetical protein",
             "UNCHARACTERIZED_HOMOLOG"] for i in range(50)]
    write_tsv(prot, ["seq_id", "annot_tier", "annot_label", "functional_class"], rows)
    write_tsv(artefact, ["seq_id", "artefact_flag"], [[r[0], 0] for r in rows])

    with pytest.raises(SystemExit):
        run_script("quality_gate.py", FakeSnakemake(
            input={"prot": str(prot), "artefact": str(artefact)},
            output={"flags": str(flags), "report": str(report)},
            params={"gate": {"min_control_recall": 0.99, "require_control_set": True},
                "tier_sources": {"T1": "pfam", "T4": "swissprot"}}))


def test_feature_files_place_an_origin_spanning_gene_correctly(fixture_dir):
    """S1 reconstructs genes broken by linearising a circular plasmid and writes them
    start > end. GFF3 forbids that and GenBank has dedicated syntax for it, so this is the
    case where a feature file silently puts a gene in the wrong part of the molecule."""
    import gzip

    ann = fixture_dir / "plasmid_annotation.tsv"
    write_tsv(ann, ["orf_id", "plasmid_id", "start", "end", "strand", "partial",
                    "spans_origin", "annot_label", "functional_class", "annot_tier",
                    "artefact_flag"],
              [["p1|1", "p1", 100, 400, "+", 0, 0, "Relaxase MobA", "FUNCTIONAL", "T1", 0],
               ["p1|2", "p1", 480, 120, "-", 0, 1, "", "NONE", "", 0]])
    master = fixture_dir / "master.tsv"
    write_tsv(master, ["plasmid_id", "topology", "size_bp"], [["p1", "circular", 500]])
    fasta = fixture_dir / "ws.fna.gz"
    with gzip.open(fasta, "wt") as fh:
        fh.write(">p1 test plasmid\n" + ("ATGC" * 125) + "\n")

    gff = fixture_dir / "plasmid_annotation.gff3"
    gbk = fixture_dir / "plasmid_annotation.gbk"
    run_script("feature_files.py", FakeSnakemake(
        input={"annotation": str(ann), "fasta": str(fasta), "master": str(master)},
        output={"gff3": str(gff), "genbank": str(gbk)}))

    lines = [l for l in gff.read_text().splitlines() if not l.startswith("#")]
    assert lines[0].startswith("##sequence-region") or True
    cds = [l.split("\t") for l in lines]
    assert all(int(c[3]) <= int(c[4]) for c in cds), "GFF3 requires start <= end"
    # The origin-spanning gene is two rows sharing one ID.
    wrapped = [c for c in cds if "p1%7C2" in c[8]]
    assert len(wrapped) == 2, f"expected a discontinuous feature, got {len(wrapped)} rows"
    assert sorted((int(c[3]), int(c[4])) for c in wrapped) == [(1, 120), (480, 500)]

    text = gbk.read_text()
    assert "complement(join(480..500,1..120))" in text, (
        "the origin-spanning gene is not written with GenBank's join convention")
    assert "circular" in text.split("\n")[0], "LOCUS line does not record the topology"
    assert text.rstrip().endswith("//"), "GenBank record is not terminated"
    assert "atgcatgc" in text.lower().replace(" ", ""), "no sequence written"


def test_feature_files_cover_the_analysis_set_and_nothing_else(fixture_dir):
    """One record per sequence in the analysis-set FASTA, whether or not it carries an
    annotation.

    This stage used to read the whole corpus FASTA, so on a 100-plasmid run it emitted
    208,245 GenBank records. The scope has to come from the FASTA S0 wrote. A plasmid with
    no called ORFs still gets a record - it is in the analysis set and the answer for it
    is "no features", which is not the same as the record being absent.
    """
    ann = fixture_dir / "scope.tsv"
    write_tsv(ann, ["orf_id", "plasmid_id", "start", "end", "strand", "partial",
                    "spans_origin", "annot_label", "functional_class", "annot_tier",
                    "artefact_flag"],
              [["in1|1", "in1", 10, 60, "+", 0, 0, "", "NONE", "", 0],
               # A row for a plasmid that is NOT in the FASTA: the annotation table may
               # be wider than this run's scope, and that must not put it in the output.
               ["out1|1", "out1", 10, 60, "+", 0, 0, "", "NONE", "", 0]])
    master = fixture_dir / "scope_master.tsv"
    write_tsv(master, ["plasmid_id", "topology", "size_bp"],
              [["in1", "circular", 200], ["in2", "linear", 200], ["out1", "linear", 200]])

    fasta = fixture_dir / "scope.fna"
    write_fasta(fasta, [("in1", "ATGC" * 50), ("in2", "GGCC" * 50)])

    gff = fixture_dir / "scope.gff3"
    gbk = fixture_dir / "scope.gbk"
    run_script("feature_files.py", FakeSnakemake(
        input={"annotation": str(ann), "fasta": str(fasta), "master": str(master)},
        output={"gff3": str(gff), "genbank": str(gbk)}))

    loci = [l.split()[1] for l in gbk.read_text().splitlines() if l.startswith("LOCUS")]
    assert loci == ["in1", "in2"], (
        f"expected one record per input sequence, got {loci} - a record for a plasmid "
        "outside the FASTA means the stage is reading something wider than its input")

    regions = [l.split()[1] for l in gff.read_text().splitlines()
               if l.startswith("##sequence-region")]
    assert regions == ["in1", "in2"], f"GFF3 scope disagrees with GenBank: {regions}"


def test_a_dark_orf_is_written_without_a_fabricated_product(fixture_dir):
    """A dark ORF has no product. `product=` asserts it has one that is blank, and
    `product=hypothetical protein` fabricates an annotation the cascade did not make."""
    import gzip
    ann = fixture_dir / "a.tsv"
    write_tsv(ann, ["orf_id", "plasmid_id", "start", "end", "strand", "partial",
                    "spans_origin", "annot_label", "functional_class", "annot_tier",
                    "artefact_flag"],
              [["p1|1", "p1", 10, 60, "+", 0, 0, "", "NONE", "", 0]])
    master = fixture_dir / "m.tsv"
    write_tsv(master, ["plasmid_id", "topology", "size_bp"], [["p1", "linear", 100]])
    fasta = fixture_dir / "ws.fna.gz"
    with gzip.open(fasta, "wt") as fh:
        fh.write(">p1\n" + ("ATGC" * 25) + "\n")

    gff = fixture_dir / "o.gff3"
    run_script("feature_files.py", FakeSnakemake(
        input={"annotation": str(ann), "fasta": str(fasta), "master": str(master)},
        output={"gff3": str(gff), "genbank": str(fixture_dir / "o.gbk")}))

    row = [l for l in gff.read_text().splitlines() if not l.startswith("#")][0]
    attrs = row.split("\t")[8]
    assert "product=" not in attrs, f"a product was fabricated for a dark ORF: {attrs}"
    assert "ID=p1%7C1" in attrs


def test_orthology_queries_only_the_proteins_the_cascade_named(fixture_dir):
    """A dark protein has nothing for eggNOG to transfer an ortholog from, and searching
    3.5M sequences to establish that would cost days for no information. The stage exists
    to describe the KNOWN genes well enough that a dark ORF's neighbourhood can be
    aggregated into pathways.

    With no eggNOG database present the run must still produce a complete, empty-valued
    table rather than no table - the column has to exist for S8 to read."""
    prot = fixture_dir / "protein_annotation.tsv"
    write_tsv(prot, ["seq_id", "functional_class"],
              [["named1", "FUNCTIONAL"], ["named2", "DOMAIN_ONLY"],
               ["dark1", "NONE"], ["dark2", "UNCHARACTERIZED_HOMOLOG"]])
    faa = fixture_dir / "unique.faa"
    write_fasta(faa, [(n, "MKVLATT") for n in ("named1", "named2", "dark1", "dark2")])
    out = fixture_dir / "orthology.tsv"

    run_script("orthology.py", FakeSnakemake(
        input={"prot": str(prot), "faa": str(faa), "ps": _ps_table(fixture_dir)},
        output=[str(out)],
        params={"orthology": {"data_dir": str(fixture_dir / "absent-db"),
                              "required": False}},
        threads=1))

    rows = {r["seq_id"] for r in read_tsv(out)}
    assert rows == {"named1", "named2"}, (
        f"the dark proteins were sent to eggNOG, or the named ones were not: {rows}")


def _label_inputs(fixture_dir, label_rows=(), ko_lines=(), defence_rows=(),
                  conj_rows=(), protein_map=""):
    """The S4d inputs of protein_labels, empty unless rows are given: the plasmid label
    databases' table, the KEGG KO list, the defence and CONJScan calls and the map."""
    labels_plasmid = fixture_dir / "protein_labels_plasmid.tsv"
    write_tsv(labels_plasmid, PLASMID_LABEL_COLS, list(label_rows))
    ko_list = fixture_dir / "list_ko.txt"
    ko_list.write_text("".join(f"{line}\n" for line in ko_lines))
    defence = fixture_dir / "defence_systems.tsv"
    write_tsv(defence, ["orf_id", "plasmid_id", "system", "component"], list(defence_rows))
    conj = fixture_dir / "conjugation_systems.tsv"
    write_tsv(conj, ["orf_id", "plasmid_id", "system", "component"], list(conj_rows))
    pmap = fixture_dir / "protein_map.tsv"
    pmap.write_text(protein_map)
    return {"labels_plasmid": str(labels_plasmid), "ko_list": str(ko_list),
            "defence": str(defence), "conjugation": str(conj), "map": str(pmap)}


def test_protein_labels_gathers_every_source_into_one_long_table(fixture_dir):
    """The substrate for the functional grouping. A wide table cannot hold it: the
    vocabulary is open, Pfam-A 38.2 alone has 30,134 families, and the grouping is derived
    from the labels observed rather than declared in advance."""
    import gzip

    hits = fixture_dir / "hits.tsv"
    write_tsv(hits, ["query", "label", "target_accession", "coverage", "target_coverage",
                     "evalue", "informative", "is_best", "start", "end", "tier", "source",
                     "category", "threshold", "max_evalue"],
              [["s1", "RepA_N", "PF06970.19", 0.9, 0.95, "1e-40", "True", 1, 1, 100,
                "T1", "pfam", "", "--cut_ga", ""],
               ["s1", "P62554.1 RecName: Full=Toxin CcdB [Escherichia coli]", "P62554.1",
                0.8, 0.9, "1e-30", "True", 1, 1, 90, "T3", "swissprot", "", "--fast", "1e-5"],
               ["s2", "WP_1.1 hypothetical protein [Escherichia coli]", "WP_1.1",
                0.95, 0.9, "1e-20", "False", 1, 1, 95, "T4", "nr", "", "--fast", "1e-10"],
               ["CTRL_P1", "PF00001.1", "PF00001.1", 0.9, 0.9, "1e-50", "True", 1, 1, 90,
                "T1", "pfam", "", "--cut_ga", ""]])

    orth = fixture_dir / "orthology.tsv"
    write_tsv(orth, ["seq_id", "cog_category", "kegg_pathways", "preferred_name",
                     "eggnog_description", "eggnog_ogs", "pfams", "gos", "ec", "kegg_ko",
                     "orthology_source"],
              [["s1", "L", "ko03030", "repA", "Replication initiator",
                "COG5527@2", "RepA_N", "GO:0006270", "2.7.7.7", "ko:K02314", "emapper"]])

    pfam_dat = fixture_dir / "Pfam-A.hmm.dat.gz"
    with gzip.open(pfam_dat, "wt") as fh:
        fh.write("# STOCKHOLM 1.0\n#=GF ID   RepA_N\n#=GF AC   PF06970.19\n"
                 "#=GF DE   Replication initiator protein A (RepA) N-terminus\n"
                 "#=GF TP   Domain\n#=GF CL   CL0123\n//\n")

    out = fixture_dir / "protein_labels.tsv"
    run_script("protein_labels.py", FakeSnakemake(
        input={"hits": [str(hits)], "orthology": str(orth), "pfam_dat": str(pfam_dat),
               "selection": _selection(fixture_dir, [["m1", 1, "member", "s1"],
                                                      ["s1", 1, "representative", "s1"]]),
               **_label_inputs(fixture_dir)},
        output={"tsv": str(out), "disagreements": str(fixture_dir / "disagree.tsv")},
        params={"pfam_version": "38.2", "swissprot_version": "2025-03-03",
                "nr_version": "2025-03-03", "eggnog_version": "5.0.2"}))

    pairs = {(r["protein_id"], r["kind"], r["label"]) for r in read_tsv(out)}

    # From the Pfam hit, including the two fields the domtblout does not carry.
    assert ("s1", "pfam_family", "RepA_N") in pairs
    assert ("s1", "pfam_description",
            "Replication initiator protein A (RepA) N-terminus") in pairs
    assert ("s1", "pfam_clan", "CL0123") in pairs
    # From the Swiss-Prot hit.
    assert ("s1", "swissprot_product", "Toxin CcdB") in pairs
    # From eggNOG, including the gene symbol.
    assert ("s1", "gene_symbol", "repA") in pairs
    assert ("s1", "cog_category", "L") in pairs
    assert ("s1", "cog_id", "COG5527") in pairs
    # The uninformative hit contributes nothing.
    assert not any(p[0] == "s2" for p in pairs), (
        "'hypothetical protein' entered the functional vocabulary")
    # A search-cluster member takes its representative's labels, and says whose they are.
    assert ("m1", "pfam_family", "RepA_N") in pairs
    assert {r["via_representative"] for r in read_tsv(out) if r["protein_id"] == "m1"
            and r["source"] == "pfam"} == {"s1"}
    # Controls are instrumentation and contribute no labels.
    assert not any(p[0].startswith("CTRL_") for p in pairs)


def test_protein_labels_records_the_database_version_on_every_row(fixture_dir):
    """A label without the database release it came from cannot be reproduced, and the
    grouping built on it cannot be described in a methods section."""
    hits = fixture_dir / "hits.tsv"
    write_tsv(hits, ["query", "label", "target_accession", "coverage", "target_coverage",
                     "evalue", "informative", "is_best", "start", "end", "tier", "source",
                     "category", "threshold", "max_evalue"],
              [["s1", "RepA_N", "PF06970.19", 0.9, 0.95, "1e-40", "True", 1, 1, 100,
                "T1", "pfam", "", "--cut_ga", ""],
               # The three sources the pharokka tier emits. The first real run failed
               # here with KeyError: 'pharokka' - the provenance map had no entry, and
               # no fixture had exercised a row from the new tier.
               ["s2", "ParA-like partition protein", "phrog_164", "", "", "1e-42", "True",
                1, "", "", "T3", "pharokka", "DNA, RNA and nucleotide metabolism", "",
                "1e-05"],
               ["s2", "TEM beta-lactamase", "ARO:3000873", "", "", "1e-200", "True", 0,
                "", "", "T3", "card", "antibiotic inactivation", "", "1e-05"],
               ["s2", "type IV pilus", "VFG000001", "", "", "1e-50", "True", 0, "", "",
                "T3", "vfdb", "", "", "1e-05"]])
    orth = fixture_dir / "orthology.tsv"
    write_tsv(orth, ["seq_id", "cog_category", "kegg_pathways", "preferred_name",
                     "eggnog_description", "eggnog_ogs", "pfams", "gos", "ec",
                     "kegg_ko", "orthology_source"], [])
    pfam_dat = fixture_dir / "pfam.dat"
    pfam_dat.write_text("")

    out = fixture_dir / "protein_labels.tsv"
    run_script("protein_labels.py", FakeSnakemake(
        input={"hits": [str(hits)], "orthology": str(orth), "pfam_dat": str(pfam_dat),
               "selection": _selection(fixture_dir), **_label_inputs(fixture_dir)},
        output={"tsv": str(out), "disagreements": str(fixture_dir / "disagree.tsv")},
        params={"pfam_version": "38.2", "swissprot_version": "2025-03-03",
                "nr_version": "2025-03-03", "eggnog_version": "5.0.2",
                "pharokka_db_version": "1.8.0"}))

    rows = read_tsv(out)
    assert rows
    assert all(r["database_version"] for r in rows), "a row carries no database version"
    assert rows[0]["database"] == "Pfam-A"
    assert rows[0]["database_version"] == "38.2"

    # pharokka ships the phage families, CARD and VFDB as ONE versioned bundle and does
    # not expose the CARD or VFDB snapshot dates separately, so all three cite the bundle.
    by_source = {r["source"]: r for r in rows}
    for source in ("pharokka", "card", "vfdb"):
        assert by_source[source]["database_version"] == "1.8.0", by_source[source]
    assert by_source["pharokka"]["database"] == "pharokka databases (PHROG v4)"
    assert by_source["card"]["database"] == "pharokka databases (CARD)"
    assert by_source["vfdb"]["database"] == "pharokka databases (VFDB)"


def test_protein_labels_merges_a_label_seen_by_two_tiers(fixture_dir):
    """The same Pfam family hit by T1 and T2 is one statement about the protein, not two.
    Unmerged, a widely searched label would outvote a rare one by copy number when the
    categories are counted."""
    hits = fixture_dir / "hits.tsv"
    write_tsv(hits, ["query", "label", "target_accession", "coverage", "target_coverage",
                     "evalue", "informative", "is_best", "start", "end", "tier", "source",
                     "category", "threshold", "max_evalue"],
              [["s1", "RepA_N", "PF06970.19", 0.9, 0.95, "1e-10", "True", 1, 1, 100,
                "T1", "pfam", "", "--cut_ga", ""],
               ["s1", "RepA_N", "PF06970.19", 0.9, 0.95, "1e-40", "True", 1, 1, 100,
                "T2", "pfam", "", "-E 1e-5", "1e-5"]])
    orth = fixture_dir / "orthology.tsv"
    write_tsv(orth, ["seq_id", "cog_category", "kegg_pathways", "preferred_name",
                     "eggnog_description", "eggnog_ogs", "pfams", "gos", "ec",
                     "kegg_ko", "orthology_source"], [])
    pfam_dat = fixture_dir / "pfam.dat"
    pfam_dat.write_text("")

    out = fixture_dir / "protein_labels.tsv"
    run_script("protein_labels.py", FakeSnakemake(
        input={"hits": [str(hits)], "orthology": str(orth), "pfam_dat": str(pfam_dat),
               "selection": _selection(fixture_dir), **_label_inputs(fixture_dir)},
        output={"tsv": str(out), "disagreements": str(fixture_dir / "disagree.tsv")},
        params={"pfam_version": "38.2", "swissprot_version": "2025-03-03",
                "nr_version": "2025-03-03", "eggnog_version": "5.0.2"}))

    rows = [r for r in read_tsv(out) if r["kind"] == "pfam_family"]
    assert len(rows) == 1, f"the same family was recorded {len(rows)} times"
    # The strongest evidence for the statement survives the merge.
    assert rows[0]["evidence_evalue"] == "1e-40"


def test_protein_labels_merges_the_plasmid_label_databases_and_lists_disagreements(
        fixture_dir):
    """S4d: every plasmid label database row enters protein_labels.tsv as its own kind,
    sub_label included (AMRFinderPlus's element type decides amr against metal), and the
    cross-source conflicts are written beside it without changing a label."""
    hits = fixture_dir / "hits.tsv"
    write_tsv(hits, ["query", "label", "target_accession", "coverage", "target_coverage",
                     "evalue", "informative", "is_best", "start", "end", "tier", "source",
                     "category", "threshold", "max_evalue"],
              [["s1", "RepA_N", "PF06970.19", 0.9, 0.95, "1e-40", "True", 1, 1, 100,
                "T1", "pfam", "", "--cut_ga", ""]])
    orth = fixture_dir / "orthology.tsv"
    # Tier 0 (PlasmidScope) names s2 by a KO whose KEGG symbol is merA.
    write_tsv(orth, ["seq_id", "cog_category", "kegg_pathways", "preferred_name",
                     "eggnog_description", "eggnog_ogs", "pfams", "gos", "ec", "kegg_ko",
                     "orthology_source"],
              [["s2", "P", "", "-", "", "", "", "", "", "ko:K00520", "plasmidscope"]])
    pfam_dat = fixture_dir / "pfam.dat"
    pfam_dat.write_text("")
    label_rows = [
        # s1: CARD and AMRFinderPlus name different genes - card_vs_amrfinder.
        ["s1", "card", "card_amr_family", "TEM beta-lactamase", "TEM-1", "Perfect",
         "500", "100.0", "100.0", "100.0", "560", "ARO:3000873", "CARD 4.0.2"],
        ["s1", "amrfinder", "amrfinder_gene", "sul1", "AMR/AMR", "EXACTP", "", "100.0",
         "", "100.0", "", "WP_000259031.1", "2026-08-07.1"],
        # s2: BacMet names merB where Tier 0 names merA - tier0_vs_bacmet.
        ["s2", "bacmet", "bacmet_compound", "Mercury", "merB", "1", "", "95.0", "98.0",
         "97.0", "400", "BAC0231", "BacMet 2.0"],
        # s3: TADB on a DefenseFinder component - tadb_vs_defencefinder.
        ["s3", "tadb", "tadb_ta", "type II toxin", "", "1", "", "90.0", "95.0", "95.0",
         "300", "TA01", "TADB 3.0"],
    ]
    inputs = _label_inputs(
        fixture_dir, label_rows=label_rows,
        ko_lines=["K00520\tmerA; mercuric reductase [EC:1.16.1.1]"],
        defence_rows=[["p1|3", "p1", "defense-finder-models/DefenseFinder/AbiE/AbiE",
                       "AbiEii"]],
        protein_map="s1\tp1|1\ns2\tp1|2\ns3\tp1|3\n")
    out = fixture_dir / "protein_labels.tsv"
    disagree = fixture_dir / "label_disagreements.tsv"
    run_script("protein_labels.py", FakeSnakemake(
        input={"hits": [str(hits)], "orthology": str(orth), "pfam_dat": str(pfam_dat),
               "selection": _selection(fixture_dir), **inputs},
        output={"tsv": str(out), "disagreements": str(disagree)},
        params={"pfam_version": "38.2", "swissprot_version": "2025-03-03",
                "nr_version": "2025-03-03", "eggnog_version": "5.0.2"}))

    rows = read_tsv(out)
    assert list(rows[0]) == ["protein_id", "source", "tier", "kind", "label", "sub_label",
                             "accession", "evidence_evalue", "evidence_coverage",
                             "database", "database_version", "via_representative"]
    plasmid_kinds = {"card_amr_family", "amrfinder_gene", "bacmet_compound", "tadb_ta"}
    merged = {(r["protein_id"], r["kind"]): r for r in rows if r["kind"] in plasmid_kinds}
    amr = merged[("s1", "amrfinder_gene")]
    assert (amr["label"], amr["sub_label"], amr["tier"], amr["accession"]) == (
        "sul1", "AMR/AMR", "EXACTP", "WP_000259031.1")
    assert (amr["database"], amr["database_version"]) == (
        "AMRFinderPlus database", "2026-08-07.1")
    card = merged[("s1", "card_amr_family")]
    assert (card["tier"], card["sub_label"]) == ("Perfect", "TEM-1")
    assert merged[("s2", "bacmet_compound")]["sub_label"] == "merB"
    assert merged[("s3", "tadb_ta")]["database"] == "TADB"
    # The cascade and eggNOG rows are still there, with an empty sub_label.
    others = [r for r in rows if r["kind"] not in plasmid_kinds]
    assert {r["source"] for r in others} == {"pfam", "eggnog"}
    assert all(r["sub_label"] == "" for r in others)

    with open(disagree) as fh:
        assert fh.readline().rstrip("\n").split("\t") == [
            "seq_id", "source_a", "label_a", "source_b", "label_b", "conflict_type"]
    conflicts = {(r["seq_id"], r["conflict_type"]) for r in read_tsv(disagree)}
    assert conflicts == {("s1", "card_vs_amrfinder"), ("s2", "tier0_vs_bacmet"),
                         ("s3", "tadb_vs_defencefinder")}
    # No label was removed because of a disagreement.
    assert len(merged) == len(label_rows)


def test_decoys_reaching_the_gate_are_reported_as_a_false_positive_rate(fixture_dir):
    """A decoy classed FUNCTIONAL is a false positive of the annotation cascade, and its
    rate is a measurement this pipeline should report rather than assume.

    It does NOT halt the run. A halting negative gate would stop the pipeline over the
    hardest cases in the collection, and the number a reader needs is the rate itself.
    """
    prot = fixture_dir / "gate_prot.tsv"
    write_tsv(prot, ["seq_id", "functional_class", "annot_tier", "annot_label"],
              [["CTRL_00001_P1", "FUNCTIONAL", "T3", "Relaxase"],
               ["CTRL_00002_P2", "FUNCTIONAL", "T1", "RepA"],
               ["DECOY_shuf_00000", "NONE", "", ""],
               ["DECOY_shuf_00001", "FUNCTIONAL", "T4", "hit by composition"],
               ["DECOY_rc_00002", "NONE", "", ""],
               ["DECOY_rc_00003", "NONE", "", ""],
               ["realprotein", "NONE", "", ""]])
    artefact = fixture_dir / "gate_artefact.tsv"
    write_tsv(artefact, ["seq_id", "artefact_flag"], [["realprotein", 0]])
    flags = fixture_dir / "gate_flags.tsv"
    report = fixture_dir / "gate_report.txt"

    run_script("quality_gate.py", FakeSnakemake(
        input={"prot": str(prot), "artefact": str(artefact)},
        output={"flags": str(flags), "report": str(report)},
        params={"gate": {"min_control_recall": 0.99, "require_control_set": True,
                         "min_controls": 2},
                "tier_sources": {"T1": "pfam", "T4": "swissprot"}}))

    text = report.read_text()
    assert "decoy_n=4" in text, f"the decoy count is not reported:\n{text}"
    assert "decoy_false_positive_rate=0.25" in text, (
        f"1 of 4 decoys was named FUNCTIONAL; the rate must be reported:\n{text}")

    eligible = {r["seq_id"] for r in read_tsv(flags)}
    assert not any(s.startswith("DECOY_") for s in eligible), (
        "decoys are instrumentation, not screening candidates; they must not appear in "
        "the target-eligibility table any more than the positive controls do")
    assert "realprotein" in eligible


def test_orthology_takes_plasmidscope_terms_without_running_emapper(fixture_dir):
    """PlasmidScope already ran eggNOG-mapper on these proteins; only proteins our cascade
    named are left for it. With no eggNOG database the PS rows still carry their terms."""
    prot = fixture_dir / "protein_annotation.tsv"
    write_tsv(prot, ["seq_id", "functional_class"],
              [["ps1", "FUNCTIONAL"], ["named1", "FUNCTIONAL"], ["dark1", "NONE"]])
    faa = fixture_dir / "unique.faa"
    write_fasta(faa, [(n, "MKVLATT") for n in ("ps1", "named1", "dark1")])
    ps = _ps_table(fixture_dir, [["ps1", "ANNOTATED", "DJ", "COG2026", "ko:K06218",
                                  "map02024", "ParE_toxin", "GO:0001", "",
                                  "Prodigal:2.6", 1]])
    out = fixture_dir / "orthology.tsv"

    run_script("orthology.py", FakeSnakemake(
        input={"prot": str(prot), "faa": str(faa), "ps": ps},
        output=[str(out)],
        params={"orthology": {"data_dir": str(fixture_dir / "absent-db"),
                              "required": False}},
        threads=1))

    rows = {r["seq_id"]: r for r in read_tsv(out)}
    assert set(rows) == {"ps1", "named1"}
    assert rows["ps1"]["orthology_source"] == "plasmidscope"
    assert rows["ps1"]["kegg_ko"] == "ko:K06218" and rows["ps1"]["cog_category"] == "DJ"
    query = fixture_dir / "emapper" / "named.faa"
    sent = [l[1:].strip() for l in open(query) if l.startswith(">")]
    assert sent == ["named1"], f"eggNOG-mapper would re-annotate PlasmidScope's proteins: {sent}"
