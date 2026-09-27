"""Smoke tests: the final annotation report.

Each test runs one workflow script against a small fixture.
"""
from conftest import PLASMID_LABEL_COLS, FakeSnakemake, _is_table, run_script, write_tsv


def _report_fixture(fixture_dir):
    """Two dark families: F1 an ORPHAN with nothing measurable, F2 a family
    with evidence on every axis. Every table the report joins is built here so
    that both report tests see the same collection."""
    ann = fixture_dir / "plasmid_annotation.tsv"
    write_tsv(ann, ["orf_id", "plasmid_id", "start", "end", "strand", "annot_label",
                    "functional_class", "artefact_flag"],
              [["p1|1", "p1", 1, 90, "+", "Relaxase MobA", "FUNCTIONAL", 0],
               ["p1|2", "p1", 100, 200, "+", "", "NONE", 1],
               ["p1|3", "p1", 300, 400, "+", "", "NONE", 0]])
    pmap = fixture_dir / "map.tsv"
    pmap.write_text("S1\tp1|1\nS2\tp1|2\nS3\tp1|3\n")
    fams = fixture_dir / "families.tsv"
    write_tsv(fams, ["family_id", "representative", "n_members", "n_orfs", "n_plasmids",
                     "n_mob_clusters", "family_class", "members"],
              [["F1", "S2", 1, 1, 1, 1, "ORPHAN", "S2"],
               ["F2", "S3", 4, 9, 9, 3, "FAMILY", "S3"]])
    evo = fixture_dir / "evo.tsv"
    write_tsv(evo, ["family_id", "dnds_median", "dnds_status", "under_purifying_selection",
                    "rnacode_p", "rnacode_status", "coding_signal"],
              [["F1", "", "TOO_FEW_MEMBERS", "", "", "TOO_FEW_MEMBERS", ""],
               ["F2", 0.21, "MEASURED", 1, 0.002, "MEASURED", 1]])
    rec = fixture_dir / "recheck.tsv"
    write_tsv(rec, ["family_id", "consensus_hit", "consensus_label", "collectively_novel"],
              [["F1", 0, "", 1], ["F2", 0, "", 1]])
    ctx = fixture_dir / "ctx.tsv"
    write_tsv(ctx, ["family_id", "n_units", "cons_defence", "cons_integron",
                    "cons_is_element", "cons_annotated_neighbour",
                    "cons_operon_with_annotated", "cons_two_gene_operon", "cons_conj"],
              [["F2", 9, 0.8, 0.0, 0.0, 1.0, 0.6, 0.2, 0.3]])
    struct = fixture_dir / "struct.tsv"
    write_tsv(struct, ["seq_id", "target", "target_description", "evalue"], [])
    orth = fixture_dir / "orth.tsv"
    write_tsv(orth, ["seq_id", "cog_category", "kegg_pathways", "preferred_name",
                     "eggnog_description"],
              [["S1", "L", "ko03430", "mobA", "Relaxase"]])
    recur = fixture_dir / "recurrence.tsv"
    write_tsv(recur, ["family_id", "family_resolution", "representative",
                      "plasmid_occurrence_count", "unique_plasmid_count",
                      "independent_plasmid_cluster_count", "independent_cluster_status",
                      "host_count", "genus_count", "MOB_count",
                      "habitat_count", "database_source_count"],
              [["F1", "broad", "S2", 1, 1, 1, "SUCCESS", 1, 1, 1, 1, 1],
               # 40 gene copies on 9 records that are only 2 independent lineages.
               ["F2", "broad", "S3", 40, 9, 2, "SUCCESS", 3, 2, 3, 2, 1]])
    syn = fixture_dir / "synteny.tsv"
    measures = ["n_occurrences", "context_recurrence", "n_lineages",
                "n_lineages_discordant", "lineage_left_conservation",
                "lineage_right_conservation", "lineage_neighborhood_conservation",
                "lineage_operon_like_conservation", "lineage_synteny_conservation",
                "modal_left", "modal_right", "modal_synteny", "status"]
    f1 = [1, 1, 1, 0, "", "", "", "", "", "", "", "", "TOO_FEW_LINEAGES"]
    f2 = [9, 9, 3, 1, 0.9, 0.7, 0.8, 0.6, 0.7, "mobA", "repA", "mobA|repA", "SUCCESS"]
    write_tsv(syn, ["family_id", "level", "intermediate_family_ids",
                    "synteny_min_lineages", *measures, *(f"small_{m}" for m in measures)],
              # Stage 9 measures the close level too; the family table reads the
              # primary rows and the ORF table the close rows.
              [["close:S3", "close", "F2", 2, 4, 4, 2, 0, 1.0, 1.0, 1.0, 0.5, 1.0,
                "close:S1", "", "close:S1|", "SUCCESS", *f2],
               ["F1", "intermediate", "F1", 2, *f1, *f1],
               ["F2", "intermediate", "F2", 2, *f2, *f2]])
    rarity_tsv = fixture_dir / "family_rarity.tsv"
    write_tsv(rarity_tsv, ["family_id", "rarity_labels",
                           "independent_plasmid_cluster_count", "unique_plasmid_count",
                           "MOB_count", "host_count", "genus_count",
                           "rare_max_lineages", "widespread_min_lineages"],
              [["F1", "RARE,LINEAGE_SPECIFIC", 1, 1, 1, 1, 1, 3, 50],
               ["F2", "RARE,CROSS_MOB", 2, 9, 3, 3, 2, 3, 50]])

    return (ann, pmap, fams, evo, rec, ctx, struct, orth, recur, syn,
            rarity_tsv)


def _run_report(fixture_dir, *tables):
    """Drive annotation_report.py. With no tables passed, build the standard fixture."""
    if not tables:
        tables = _report_fixture(fixture_dir)
    ann, pmap, fams, evo, rec, ctx, struct, orth, recur, syn, rarity_tsv = tables
    # One IS element over p1|3 (300-400), none elsewhere.
    is_tsv = _is_table(fixture_dir,
                       [["p1", "p1|IS1", "IS3", "IS3_1", 280, 1500, "+", 1, "1e-50", ""]])
    registry = fixture_dir / "report_registry.tsv"
    write_tsv(registry, ["plasmid_id", "species", "genus", "predicted_host_range"],
              [["p1", "Escherichia coli", "Escherichia", "Enterobacterales"]])
    # Close clusters: S3 with S2 (measured by Stage 9), S1 alone (no close row).
    clusters_close = fixture_dir / "families_close_cluster.tsv"
    clusters_close.write_text("S1\tS1\nS3\tS3\nS3\tS2\n")
    # S1 is the relaxase: oriTDB and CARD label it, and an AMRFinderPlus VIRULENCE element
    # (no term type, so no column); CONJScan calls p1|1 the MOB of a MOB system.
    labels_plasmid = fixture_dir / "protein_labels_plasmid.tsv"
    write_tsv(labels_plasmid, PLASMID_LABEL_COLS,
              [["S1", "oritdb", "oritdb_role", "relaxase", "MOBP", "1", "", "95", "98",
                "97", "500", "TraI_RP4", "oriTDB 2.0"],
               ["S1", "card", "card_amr_family", "sulfonamide resistant sul, x", "sul1",
                "Strict", "300", "80", "100", "100", "400", "ARO:1", "CARD 4.0.2"],
               ["S1", "amrfinder", "amrfinder_gene", "sul1", "AMR/AMR", "EXACTP", "",
                "100", "", "100", "", "WP_1", "2026-08-07.1"],
               ["S1", "amrfinder", "amrfinder_gene", "iutA", "VIRULENCE/VIRULENCE",
                "BLASTP", "", "99", "", "100", "", "WP_2", "2026-08-07.1"]])
    conj = fixture_dir / "conjugation_systems.tsv"
    write_tsv(conj, ["orf_id", "plasmid_id", "system", "system_id", "component",
                     "hit_status", "sys_wholeness", "conjscan_version"],
              [["p1|1", "p1", "MOB", "p1_MOB_1", "T4SS_MOBP1", "mandatory", "1.000",
                "2.1.0"]])
    conj_class = fixture_dir / "conjugation_plasmid_class.tsv"
    write_tsv(conj_class, ["plasmid_id", "class"], [["p1", "pMOB"]])
    # S8g: F2's member S3 travels with S9 (q 0.01) and less clearly with S8; F1's member
    # S2 was in no reported pair.
    cooc = fixture_dir / "dark_cooccurrence.tsv"
    write_tsv(cooc, ["seq_a", "seq_b", "n_lineages_a", "n_lineages_b",
                     "n_lineages_together", "n_lineages_total", "fraction_of_a",
                     "fraction_of_b", "expected_together", "p_value", "q_value"],
              [["S3", "S9", 3, 2, 2, 100, 0.6667, 1.0, 0.06, 0.001, 0.01],
               ["S3", "S8", 3, 40, 2, 100, 0.6667, 0.05, 1.2, 0.3, 0.3]])
    out_ann = fixture_dir / "annotation_complete.csv"
    out_fam = fixture_dir / "dark_families_complete.csv"
    run_script("annotation_report.py", FakeSnakemake(
        input={"annotation": str(ann), "map": str(pmap), "families": str(fams),
               "evolution": str(evo), "recheck": str(rec),
               "context": str(ctx), "structure": str(struct), "orthology": str(orth),
               "recurrence": str(recur), "synteny": str(syn),
               "rarity": str(rarity_tsv), "is_elements": is_tsv,
               "registry": str(registry), "clusters_close": str(clusters_close),
               "labels_plasmid": str(labels_plasmid), "conjugation": str(conj),
               "conjugation_class": str(conj_class), "cooccurrence": str(cooc)},
        output={"annotation": str(out_ann), "families": str(out_fam)},
        params={"evolution": {"min_members_for_dnds": 3, "dnds_purifying_max": 0.5},
                "cooccurrence": {"min_lineages_together": 2, "fdr": 0.05}}))
    return out_ann, out_fam


def test_the_report_carries_every_orf_and_every_family(fixture_dir):
    """The pipeline's output is every annotation it could produce, in a form you can sort
    and filter yourself. Nothing is dropped for being artefactual, ORPHAN, or evidence-free,
    and nothing is ranked - selecting candidates is a decision made on this table, not one
    baked into a rule."""
    out_ann, out_fam = _run_report(fixture_dir)

    import csv as _csv
    orfs = list(_csv.DictReader(open(out_ann)))
    assert len(orfs) == 3, "an ORF was dropped from the complete annotation"
    by_orf = {r["orf_id"]: r for r in orfs}
    # The annotated ORF carries its orthology terms.
    assert by_orf["p1|1"]["kegg_pathways"] == "ko03430"
    # The ORF's own plasmid: observed host, and MOB-suite's range as a separate column.
    assert (by_orf["p1|1"]["host_species"], by_orf["p1|1"]["host_genus"],
            by_orf["p1|1"]["predicted_host_range"]) == (
        "Escherichia coli", "Escherichia", "Enterobacterales")
    # The artefact-flagged dark ORF is present, flagged, and joined to its family.
    assert by_orf["p1|2"]["artefact_flag"] == "1"
    assert by_orf["p1|2"]["family_id"] == "F1"
    # Family evidence travels down to the ORF row.
    assert by_orf["p1|3"]["dnds_median"] == "0.21"

    fam_rows = {r["family_id"]: r for r in _csv.DictReader(open(out_fam))}
    assert set(fam_rows) == {"F1", "F2"}, "an ORPHAN family was dropped"
    assert fam_rows["F1"]["dnds_status"] == "TOO_FEW_MEMBERS", (
        "absence of a measurement must be an explicit status, not a blank")
    assert fam_rows["F1"]["reality_n"] == "0"
    # purifying_selection fired, so is_family is entailed and does not count twice.
    assert fam_rows["F2"]["reality_n"] == "2"
    assert fam_rows["F2"]["reality_lines_implied"] == "is_family"
    # multi_lineage is read from Stage 6 lineages (2 for F2), not MOB-suite clusters.
    assert "multi_lineage" in fam_rows["F2"]["reality_lines"].split("+")
    # The column sets are the contract, asserted by equality, so an added, renamed or
    # removed column fails here. It also keeps any score or rank column out of the report.
    FAMILY_COLUMNS = [
        "family_id", "representative", "family_class",
        "n_members", "n_orfs", "n_plasmids", "n_mob_clusters",
        "dark_member_count", "annotated_member_count", "percentage_dark_in_family",
        "dark_only",
        # small plasmids and large ones
        "n_small_members", "n_large_members", "scope", "known_from",
        "reality_n", "reality_lines", "reality_lines_implied",
        "dnds_median", "dnds_min", "dnds_status", "under_purifying_selection", "n_pairs",
        "rnacode_p", "rnacode_p_antisense", "rnacode_status", "coding_signal",
        "small_n_aligned", "small_dnds_median", "small_dnds_min", "small_dnds_status",
        "small_under_purifying_selection", "small_n_pairs", "small_rnacode_p",
        "small_rnacode_p_antisense", "small_rnacode_status", "small_coding_signal",
        "consensus_status", "consensus_hit", "consensus_label", "collectively_novel",
        "darkness_state", "structural_match", "structural_description", "structure_evalue",
        "cons_defence", "cons_integron", "cons_is_element",
        "cons_annotated_neighbour", "cons_operon_with_annotated", "cons_two_gene_operon",
        "cons_conj",
        # Stage 7: seven counts, never collapsed into one.
        "plasmid_occurrence_count", "unique_plasmid_count",
        "independent_plasmid_cluster_count", "independent_cluster_status",
        "host_count", "genus_count", "n_plasmids_with_host",
        "n_plasmids_with_species", "host_count_status", "n_plasmids_with_predicted_range",
        "predicted_host_range_count",
        "predicted_host_ranges", "MOB_count", "habitat_count",
        "database_source_count",
        # Stage 9: six conservation measurements, kept apart because they fail apart,
        # counted over lineages; the primary-level (family) rows.
        "context_recurrence", "n_occurrences", "n_lineages", "n_lineages_discordant",
        "lineage_left_conservation", "lineage_right_conservation",
        "lineage_neighborhood_conservation", "lineage_operon_like_conservation",
        "lineage_synteny_conservation", "modal_left", "modal_right", "modal_synteny",
        "synteny_status", "synteny_min_lineages",
        "small_n_occurrences", "small_context_recurrence", "small_n_lineages",
        "small_n_lineages_discordant", "small_lineage_left_conservation",
        "small_lineage_right_conservation", "small_lineage_neighborhood_conservation",
        "small_lineage_operon_like_conservation", "small_lineage_synteny_conservation",
        "small_modal_left", "small_modal_right", "small_modal_synteny",
        "small_synteny_status",
        # S8g: dark sequences its members travel with.
        "cooccurrence_status", "n_cooccurring_partners", "top_cooccurring_partner",
        "top_cooccurring_partner_q", "top_cooccurring_partner_fraction",
        "cooccurrence_fdr", "cooccurrence_min_lineages",
        # Stage 14: descriptors, not a ranking.
        "rarity_labels",
        # Stage 15: dimensions counted, never scored.
        "evidence_dimensions_present", "evidence_dimension_count",
        "supporting_observations_count", "supporting_observations_are_not_independent",
    ]
    assert list(fam_rows["F2"]) == FAMILY_COLUMNS, (
        f"family table columns changed: {list(fam_rows['F2'])}")
    CARRIED_TO_ORFS = {"family_id", "scope", "reality_n", "reality_lines", "darkness_state",
                       "dnds_median", "dnds_status", "coding_signal", "collectively_novel",
                       "structural_match", "structural_description"}
    annotation_cols = {"orf_id", "plasmid_id", "start", "end", "strand", "annot_label",
                       "functional_class", "artefact_flag"}
    # seq_id is the dereplicated-protein key that joins an ORF to its sequence; it is
    # written beside the orthology terms and is neither annotation nor family evidence.
    join_and_orthology = {"seq_id", "cog_category", "kegg_pathways", "preferred_name",
                          "eggnog_description", "is_element", "host_species", "host_genus",
                          "predicted_host_range"}
    # The ORF's own evidence beside the family's: its protein's plasmid label database
    # labels by term type, its CONJScan call, and its close cluster's synteny.
    orf_evidence = ["amr_labels", "metal_labels", "ta_labels", "conj_role_labels",
                    "mge_labels", "antidefence_labels", "conj_system", "conj_component",
                    "plasmid_conjscan_class", "close_family_id", "close_n_lineages",
                    "close_lineage_synteny_conservation", "close_modal_synteny",
                    "close_synteny_status"]
    carried = set(by_orf["p1|3"]) - annotation_cols - join_and_orthology - set(orf_evidence)
    assert carried == CARRIED_TO_ORFS, (
        f"family evidence carried to the ORF table changed: {sorted(carried)}")
    columns = list(by_orf["p1|3"])
    assert columns[columns.index("predicted_host_range") + 1:][:len(orf_evidence)] == \
        orf_evidence, f"ORF evidence columns changed: {columns}"

    # --- the ORF's own evidence ------------------------------------------------------
    relaxase = by_orf["p1|1"]
    assert relaxase["conj_role_labels"] == "oritdb:relaxase"
    assert relaxase["amr_labels"] == "amrfinder:sul1; card:sulfonamide resistant sul, x"
    assert relaxase["metal_labels"] == "", "a VIRULENCE element has no term type"
    assert (relaxase["conj_system"], relaxase["conj_component"]) == ("MOB", "T4SS_MOBP1")
    assert {r["plasmid_conjscan_class"] for r in orfs} == {"pMOB"}
    assert by_orf["p1|3"]["conj_system"] == ""
    # S2 and S3 share the close cluster close:S3, which Stage 9 measured; S1's was not.
    for orf in ("p1|2", "p1|3"):
        assert (by_orf[orf]["close_family_id"], by_orf[orf]["close_n_lineages"],
                by_orf[orf]["close_lineage_synteny_conservation"],
                by_orf[orf]["close_modal_synteny"], by_orf[orf]["close_synteny_status"]) \
            == ("close:S3", "2", "1.0", "close:S1|", "SUCCESS")
    assert relaxase["close_family_id"] == relaxase["close_synteny_status"] == ""

    # --- Stage 7: the counts stay apart ---------------------------------------------
    # F2 is 40 gene copies on 9 records that are 2 lineages: three separate counts.
    assert fam_rows["F2"]["plasmid_occurrence_count"] == "40"
    assert fam_rows["F2"]["unique_plasmid_count"] == "9"
    assert fam_rows["F2"]["independent_plasmid_cluster_count"] == "2", (
        "the independent-lineage count is what a recurrence claim needs, and it is not "
        "the record count")

    # --- Stage 9: synteny, with its own status name ---------------------------------
    assert fam_rows["F2"]["lineage_synteny_conservation"] == "0.7"
    assert (fam_rows["F2"]["n_lineages"], fam_rows["F2"]["synteny_min_lineages"]) == (
        "3", "2")
    assert fam_rows["F2"]["cons_conj"] == "0.3"
    assert fam_rows["F1"]["synteny_status"] == "TOO_FEW_LINEAGES", (
        "one lineage is perfectly conserved with itself; that must read as a status, "
        "not as a conservation of 1.0")

    # --- S8g: partners, from the tested pairs --------------------------------------
    f2 = fam_rows["F2"]
    assert (f2["cooccurrence_status"], f2["n_cooccurring_partners"],
            f2["top_cooccurring_partner"], f2["top_cooccurring_partner_q"],
            f2["top_cooccurring_partner_fraction"], f2["cooccurrence_fdr"],
            f2["cooccurrence_min_lineages"]) == (
        "SUCCESS", "1", "S9", "0.01", "0.6667", "0.05", "2")
    f1 = fam_rows["F1"]
    assert (f1["cooccurrence_status"], f1["n_cooccurring_partners"],
            f1["top_cooccurring_partner"]) == ("TOO_FEW_LINEAGES", "0", "")

    # --- Stage 14: labels are descriptors -------------------------------------------
    assert fam_rows["F1"]["rarity_labels"] == "RARE,LINEAGE_SPECIFIC"

    # --- Stage 15: dimensions counted, never scored ----------------------------------
    dims = fam_rows["F2"]["evidence_dimensions_present"].split(",")
    assert "EVOLUTIONARY_CONSERVATION" in dims, "dnds_status is a measurement"
    assert "DISTRIBUTION" in dims, "independent_cluster_status is a measurement"
    assert "GENOMIC_CONTEXT" in dims, "the context rates are a measurement"
    assert fam_rows["F2"]["cons_defence"] == "0.8"
    assert fam_rows["F2"]["evidence_dimension_count"] == str(len(dims))
    assert fam_rows["F2"]["supporting_observations_are_not_independent"] == "1", (
        "the observation count must carry its own warning, because a column selected "
        "into a downstream ranking takes the warning with it")


def test_the_report_names_the_is_family_of_an_orf_inside_an_element(fixture_dir):
    import csv as _csv
    out_ann, _ = _run_report(fixture_dir)
    by_orf = {r["orf_id"]: r for r in _csv.DictReader(open(out_ann))}
    assert by_orf["p1|3"]["is_element"] == "IS3"
    assert by_orf["p1|1"]["is_element"] == ""
