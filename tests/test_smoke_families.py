"""Smoke tests: host registry, lineages, protein families, network, rarity, recurrence.

Each test runs one workflow script against a small fixture.
"""
import pytest
from conftest import FakeSnakemake, read_tsv, requires, run_script, write_fasta, write_tsv


def _recurrence_fixture(fixture_dir, lineage_rows):
    """One family on three plasmids; the caller decides how independent those are."""
    families = fixture_dir / "protein_families.tsv"
    write_tsv(families, ["family_id", "family_resolution", "representative", "members"],
              [["broad:s1", "broad", "s1", "s1,s2"]])
    mapping = fixture_dir / "protein_map.tsv"
    # s1 on two plasmids with two copies on one of them; s2 on a third.
    mapping.write_text("s1\tpl1|1,pl1|2,pl2|1\ns2\tpl3|1\n")
    registry = fixture_dir / "clonal_registry.tsv"
    write_tsv(registry, ["plasmid_id", "mob_cluster", "species", "genus",
                         "predicted_host_range", "topology", "size_bp", "hab_top"],
              [["pl1", "MOB_A", "Escherichia coli", "Escherichia", "Enterobacterales",
                "circular", 100, "Host-associated"],
               ["pl2", "MOB_A", "Escherichia coli", "Escherichia", "Enterobacterales",
                "circular", 100, "Host-associated"],
               ["pl3", "MOB_A", "", "", "Actinomycetota,Pseudomonadota", "circular", 100,
                "Host-associated"]])
    lineage = fixture_dir / "plasmid_lineage.tsv"
    write_tsv(lineage, ["plasmid_id", "plasmid_lineage_cluster"], lineage_rows)
    master = fixture_dir / "master.tsv"
    write_tsv(master, ["plasmid_id", "sources"],
              [["pl1", "PLSDB,IMG"], ["pl2", "PLSDB"], ["pl3", "PLSDB"]])

    out = fixture_dir / "recurrence.tsv"
    run_script("recurrence.py", FakeSnakemake(
        input={"families": str(families), "map": str(mapping),
               "registry": str(registry), "lineage": str(lineage),
               "master": str(master)},
        output={"tsv": str(out)}))
    return read_tsv(out)[0]


def test_recurrence_separates_occurrences_plasmids_and_lineages(fixture_dir):
    """Spec section 34.2: 'database record counts must never be treated as independent
    biological observations.' Four gene copies on three plasmid records that are all ONE
    lineage is one independent observation, and the three numbers must not agree."""
    row = _recurrence_fixture(fixture_dir,
                              [["pl1", "L1"], ["pl2", "L1"], ["pl3", "L1"]])

    assert row["plasmid_occurrence_count"] == "4", "gene copies miscounted"
    assert row["unique_plasmid_count"] == "3"
    assert row["independent_plasmid_cluster_count"] == "1", (
        "three redepositions of one lineage were counted as independent observations")
    assert row["independent_cluster_status"] == "SUCCESS"


def test_genuinely_independent_plasmids_are_counted_as_such(fixture_dir):
    """The other direction: the conservative count must not flatten real breadth."""
    row = _recurrence_fixture(fixture_dir,
                              [["pl1", "L1"], ["pl2", "L2"], ["pl3", "L3"]])

    assert row["independent_plasmid_cluster_count"] == "3"


def test_mob_breadth_is_not_evolutionary_independence(fixture_dir):
    """Section 33: MOB classification and sequence similarity are separate concepts. All
    three plasmids share one MOB cluster while being three lineages, so the two counts
    must be able to disagree in both directions."""
    row = _recurrence_fixture(fixture_dir,
                              [["pl1", "L1"], ["pl2", "L2"], ["pl3", "L3"]])

    assert row["MOB_count"] == "1"
    assert row["independent_plasmid_cluster_count"] == "3"


def test_host_counts_say_how_many_plasmids_had_a_host(fixture_dir):
    """Two of the three plasmids have a recorded host. host_count counts species over
    those; n_plasmids_with_host says how many that is, so 1 host from 2 of 3 plasmids is
    not read as 1 host from all three."""
    row = _recurrence_fixture(fixture_dir, [["pl1", "L1"]])
    assert row["host_count"] == "1"
    assert row["genus_count"] == "1"
    assert row["n_plasmids_with_host"] == "2"
    assert row["n_plasmids_with_species"] == "2"
    assert row["host_count_status"] == "SUCCESS"
    # The MOB-suite range is a separate measurement over EVERY plasmid, hosted or not, and
    # never enters the host counts.
    assert row["n_plasmids_with_predicted_range"] == "3"
    assert row["predicted_host_range_count"] == "2"
    assert row["predicted_host_ranges"] == (
        "Actinomycetota,Pseudomonadota;Enterobacterales")


def test_clonal_registry_takes_the_host_from_three_sources(fixture_dir):
    """PLSDB species first, then PlasmidScope's per-record host (IMG/PR's among them),
    then the GenBank/RefSeq organism. A name that is not an organism is no host."""
    master = fixture_dir / "master.tsv"
    write_tsv(master, ["plasmid_id", "mob_cluster", "plsdb_species", "mob_host_range",
                       "topology", "size_bp", "hab_top"],
              [["a", "M1", "Klebsiella_pneumoniae", "Enterobacterales", "circular", 5000,
                "Host-associated"],
               ["b", "", "", "", "circular", 5000, "Environmental"],
               ["c", "M2", "", "", "linear", 5000, "Environmental"],
               ["d", "M3", "", "Bacteroides", "circular", 5000, "Environmental"]])
    ids_file = fixture_dir / "ids.txt"
    ids_file.write_text("a\nb\nc\nd\n")
    prov = fixture_dir / "complete_provenance.tsv"
    write_tsv(prov, ["plasmid_id", "host"],
              [["a", "Escherichia coli"], ["b", "Acidipila rosea"],
               ["c", "human gut metagenome"], ["d", "-"]])
    ws = fixture_dir / "working_set.tsv"
    write_tsv(ws, ["plasmid_id", "lifestyle", "organism"],
              [["a", "isolate", ""], ["b", "metagenomic", ""],
               ["c", "isolate", "Acinetobacter sp. X1"], ["d", "metagenomic", ""]])
    out = fixture_dir / "registry.tsv"

    run_script("clonal_registry.py", FakeSnakemake(
        input={"master": str(master), "ids": str(ids_file), "ps_hosts": str(prov),
               "working_set": str(ws)},
        output=[str(out)]))

    rows = {r["plasmid_id"]: r for r in read_tsv(out)}
    assert (rows["a"]["species"], rows["a"]["host_source"]) == \
        ("Klebsiella pneumoniae", "plsdb")
    assert (rows["b"]["species"], rows["b"]["host_source"]) == \
        ("Acidipila rosea", "plasmidscope")
    assert (rows["c"]["species"], rows["c"]["genus"], rows["c"]["host_source"]) == \
        ("", "Acinetobacter", "organism")
    assert rows["d"]["genus"] == "" and rows["d"]["host_source"] == ""
    # MOB-suite's predicted range is its own column, never the host.
    assert rows["d"]["predicted_host_range"] == "Bacteroides"
    assert rows["b"]["lifestyle"] == "metagenomic"


def test_unmeasured_independence_is_not_reported_as_zero(fixture_dir):
    """A family whose plasmids are absent from the lineage table has not been measured.
    Reporting 0 would read as 'no independent lineages', a much stronger claim than 'not
    measured' (section 2.9)."""
    row = _recurrence_fixture(fixture_dir, [])

    assert row["independent_cluster_status"] == "NOT_RUN"


def test_database_sources_are_provenance_not_biology(fixture_dir):
    """Section 34.1 asks for these counts so a reader can see when a number is large for a
    database reason. They are reported and are never a denominator."""
    row = _recurrence_fixture(fixture_dir,
                              [["pl1", "L1"], ["pl2", "L1"], ["pl3", "L1"]])

    assert row["database_source_count"] == "2", "PLSDB and IMG were not both counted"


@requires("mmseqs")
def test_family_network_links_a_dark_cluster_to_an_annotated_relative(fixture_dir):
    """Two 50%-identity clusters of related sequences, one annotated and one dark, and an
    unrelated dark protein. The network must link the first two, mark only the unrelated
    one as unconnected, and say so in the summary."""
    import random
    rng = random.Random(3)
    aa = "ACDEFGHIKLMNPQRSTVWY"
    base = "".join(rng.choice(aa) for _ in range(180))
    # Resampling 55% of positions leaves ~47% identity: below the 50% clustering bar, so a
    # separate cluster, but still a significant full-length alignment.
    relative = "".join(c if rng.random() > 0.55 else rng.choice(aa) for c in base)
    loner = "".join(rng.choice(aa) for _ in range(180))
    reps = fixture_dir / "reps.fasta"
    write_fasta(reps, [("known", base), ("darkrel", relative), ("loner", loner)])
    clusters = fixture_dir / "clusters.tsv"
    clusters.write_text("known\tknown\ndarkrel\tdarkrel\nloner\tloner\n")
    prot = fixture_dir / "protein_annotation.tsv"
    write_tsv(prot, ["seq_id", "functional_class", "annot_label", "explained_fraction",
                     "annot_completeness"],
              [["known", "FUNCTIONAL", "Relaxase", 0.95, "FULL"],
               ["darkrel", "NONE", "", 0.0, "NONE"],
               ["loner", "NONE", "", 0.0, "NONE"]])
    dark_ids = fixture_dir / "dark_ids.txt"
    dark_ids.write_text("darkrel\nloner\n")
    # A cluster of reference proteins alone is in MMseqs2's files but not in the family
    # table, and must not become a node.
    write_fasta(reps, [("known", base), ("darkrel", relative), ("loner", loner),
                       ("refonly", relative[::-1])])
    clusters.write_text("known\tknown\ndarkrel\tdarkrel\nloner\tloner\n"
                        "refonly\trefonly\n")
    fams = fixture_dir / "protein_families.tsv"
    write_tsv(fams, ["family_id", "family_resolution", "representative", "family_size",
                     "family_class", "scope", "members"],
              [[f"{res}:{m}", res, m, 1, "ORPHAN", scope, m]
               for res in ("broad", "intermediate")
               for m, scope in (("known", "small_only_known"),
                                ("darkrel", "mixed_unknown"),
                                ("loner", "small_only_unknown"))])
    pmap = fixture_dir / "map.tsv"
    pmap.write_text("known\tp1|1\ndarkrel\tp2|1\nloner\tp3|1\n")
    out = {k: str(fixture_dir / f"network_{k}.tsv") for k in ("nodes", "edges", "summary")}

    run_script("family_network.py", FakeSnakemake(
        input={"families": str(fams), "reps": str(reps), "clusters": str(clusters),
               "prot": str(prot), "dark_ids": str(dark_ids), "map": str(pmap)},
        output=out,
        params={"network": {"min_cov": 0.5, "max_evalue": 1e-4, "max_out_edges": 4,
                            "dark_brightness": 0.05},
                "primary": "broad", "node_resolution": "intermediate", "seed": 1},
        threads=2))

    nodes = {r["node_id"]: r for r in read_tsv(out["nodes"])}
    assert set(nodes) == {"known", "darkrel", "loner"}
    assert nodes["darkrel"]["scope"] == "mixed_unknown"
    assert nodes["known"]["dark"] == "0" and nodes["known"]["label"] == "Relaxase"
    assert nodes["darkrel"]["dark"] == "1"
    assert nodes["darkrel"]["degree"] == "1" and nodes["loner"]["degree"] == "0"
    summary = {r["metric"]: r["value"] for r in read_tsv(out["summary"])}
    assert summary["dark_nodes_connected"] == "1"
    assert summary["dark_connected_to_bright"] == "1"
    assert summary["dark_family_singletons_with_edge"] == "1"


def test_dark_cooccurrence_writes_the_tested_pairs(fixture_dir):
    """S8g on the tables the rule reads. d1 and d2 share a plasmid in lineages A and B
    (A's two redeposited copies count once); d3 is on B's plasmid too but in one lineage
    only, so no pair with it is tested; d4 shares lineage C with d1 on another plasmid,
    which is not together."""
    fams = fixture_dir / "dark_families.tsv"
    write_tsv(fams, ["family_id", "members"],
              [["intermediate:d1", "s1,s1b"], ["intermediate:d2", "s2"],
               ["intermediate:d3", "s3"], ["intermediate:d4", "s4"]])
    mapping = fixture_dir / "protein_map.tsv"
    mapping.write_text("s1\tA1|1,A2|1\ns1b\tB1|3,C1|1\ns2\tA1|2,A2|2,B1|1\n"
                       "s3\tB1|2\ns4\tC2|1,D1|1\nknown\tA1|3\n")
    lineage = fixture_dir / "plasmid_lineage.tsv"
    write_tsv(lineage, ["plasmid_id", "plasmid_lineage_cluster"],
              [["A1", "A"], ["A2", "A"], ["B1", "B"], ["C1", "C"], ["C2", "C"],
               ["D1", "D"], ["E1", "E"], ["F1", "F"]])
    out = fixture_dir / "dark_cooccurrence.tsv"

    run_script("dark_cooccurrence.py", FakeSnakemake(
        input={"families": str(fams), "map": str(mapping), "lineage": str(lineage)},
        output=[str(out)],
        params={"cooccurrence": {"min_lineages_together": 2, "fdr": 0.05}}))

    rows = read_tsv(out)
    assert list(rows[0]) == ["family_a", "family_b", "n_lineages_a", "n_lineages_b",
                             "n_lineages_together", "n_lineages_total", "fraction_of_a",
                             "fraction_of_b", "expected_together", "p_value", "q_value",
                             "status"]
    (row,) = rows
    assert (row["family_a"], row["family_b"], row["n_lineages_a"], row["n_lineages_b"],
            row["n_lineages_together"], row["n_lineages_total"]) == (
        "intermediate:d1", "intermediate:d2", "3", "2", "2", "6")
    # P(X >= 2) with N=6, K=3, n=2: C(3,2)/C(6,2) = 3/15.
    assert float(row["p_value"]) == pytest.approx(0.2)
    assert (row["fraction_of_a"], row["fraction_of_b"]) == ("0.6667", "1.0")


def test_rarefaction_samples_every_small_plasmid(fixture_dir):
    """The x-axis is the small plasmids, with a dark family or without; a dark family's
    copies on large plasmids do not put those plasmids on the axis."""
    recurrence = fixture_dir / "recurrence.tsv"
    write_tsv(recurrence, ["family_id", "independent_plasmid_cluster_count"],
              [["broad:d", 2]])
    fams = fixture_dir / "dark_families.tsv"
    write_tsv(fams, ["family_id", "members", "small_members"], [["broad:d", "d", "d"]])
    mapping = fixture_dir / "protein_map.tsv"
    mapping.write_text("d\tS1|1,L1|4\n")
    small_ids = fixture_dir / "small_plasmids.txt"
    small_ids.write_text("S1\nS2\nS3\nS4\n")
    out = fixture_dir / "rarefaction.tsv"

    run_script("rarity.py", FakeSnakemake(
        input={"recurrence": str(recurrence), "dark_families": str(fams),
               "map": str(mapping), "small_ids": str(small_ids)},
        output={"rarity": str(fixture_dir / "rarity.tsv"), "rarefaction": str(out)},
        params={"rarity": {"rare_max_lineages": 3, "widely_conserved_min_lineages": 50,
                           "cross_min_hosts": 2, "cross_min_genera": 2,
                           "rarefaction_replicates": 5},
                "seed": 1}))

    final = read_tsv(out)[-1]
    assert final["n_plasmids"] == "4"
    assert final["mean_families"] == "1.0"
