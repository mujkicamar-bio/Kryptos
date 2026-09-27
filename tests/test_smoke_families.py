"""Smoke tests: host registry, lineages, protein families, network, rarity, recurrence.

Each test runs one workflow script against a small fixture.
"""
import pytest
from conftest import FakeSnakemake, read_tsv, requires, run_script, write_fasta, write_tsv


def _recurrence_fixture(fixture_dir, lineage_rows, mob="MOB_A"):
    """One family on three plasmids, all in MOB cluster `mob`; the caller decides how
    independent those are."""
    families = fixture_dir / "protein_families.tsv"
    write_tsv(families, ["family_id", "family_resolution", "representative", "members"],
              [["broad:s1", "broad", "s1", "s1,s2"]])
    mapping = fixture_dir / "protein_map.tsv"
    # s1 on two plasmids with two copies on one of them; s2 on a third.
    mapping.write_text("s1\tpl1|1,pl1|2,pl2|1\ns2\tpl3|1\n")
    registry = fixture_dir / "clonal_registry.tsv"
    write_tsv(registry, ["plasmid_id", "mob_cluster", "species", "genus",
                         "predicted_host_range", "topology", "size_bp", "hab_top"],
              [["pl1", mob, "Escherichia coli", "Escherichia", "Enterobacterales",
                "circular", 100, "Host-associated"],
               ["pl2", mob, "Escherichia coli", "Escherichia", "Enterobacterales",
                "circular", 100, "Host-associated"],
               ["pl3", mob, "", "", "Actinomycetota,Pseudomonadota", "circular", 100,
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
    """Four gene copies on three plasmid records of one lineage: three different counts."""
    row = _recurrence_fixture(fixture_dir,
                              [["pl1", "L1"], ["pl2", "L1"], ["pl3", "L1"]])

    assert row["plasmid_occurrence_count"] == "4", "gene copies miscounted"
    assert row["unique_plasmid_count"] == "3"
    assert row["independent_plasmid_cluster_count"] == "1", (
        "three redepositions of one lineage were counted as independent observations")
    assert row["independent_cluster_status"] == "SUCCESS"


def test_independent_lineages_are_counted_apart_from_mob_clusters(fixture_dir):
    """Three lineages sharing one MOB cluster: three independent observations, one MOB
    cluster."""
    row = _recurrence_fixture(fixture_dir,
                              [["pl1", "L1"], ["pl2", "L2"], ["pl3", "L3"]])

    assert row["independent_plasmid_cluster_count"] == "3"
    assert row["MOB_count"] == "1"


def test_plasmids_without_a_mob_cluster_give_no_mob_count(fixture_dir):
    """The registry leaves mob_cluster empty when MOB-suite assigned none, so a family on
    untyped plasmids has MOB_count 0 and carries neither MOB label."""
    row = _recurrence_fixture(fixture_dir, [["pl1", "L1"]], mob="")
    assert row["MOB_count"] == "0"


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
    # MOB-suite assigned b no cluster: the column stays empty, so no count treats it as one.
    assert rows["b"]["mob_cluster"] == "" and rows["a"]["mob_cluster"] == "M1"


@requires("mash")
def test_plasmid_lineage_links_near_identical_plasmids(fixture_dir):
    """Two plasmids differing at 1% of positions are one lineage, two identical small ones
    another; the thresholds that produced the lineages are on every row."""
    import random
    rng = random.Random(5)
    big = "".join(rng.choice("ACGT") for _ in range(5000))
    near = "".join(rng.choice("ACGT") if i % 100 == 0 else c for i, c in enumerate(big))
    small = "".join(rng.choice("ACGT") for _ in range(3000))
    fasta = fixture_dir / "analysis_set.fna"
    write_fasta(fasta, [("pA", big), ("pB", near), ("pS", small), ("pS2", small)])
    out = fixture_dir / "lineage" / "plasmid_lineage.tsv"

    run_script("plasmid_lineage.py", FakeSnakemake(
        input={"fasta": str(fasta)}, output={"tsv": str(out)},
        params={"lineage": {"max_distance": 0.05, "max_pvalue": 1e-10, "kmer": 21,
                            "sketch_size": 1000}},
        threads=1))

    rows = read_tsv(out)
    assert list(rows[0]) == ["plasmid_id", "plasmid_lineage_cluster", "lineage_max_distance",
                             "lineage_kmer", "lineage_sketch_size"]
    assert {r["plasmid_id"]: r["plasmid_lineage_cluster"] for r in rows} == {
        "pA": "pA", "pB": "pA", "pS": "pS", "pS2": "pS"}
    assert {(r["lineage_max_distance"], r["lineage_kmer"]) for r in rows} == {("0.05", "21")}


def test_unmeasured_independence_is_not_reported_as_zero(fixture_dir):
    """A family whose plasmids are absent from the lineage table has not been measured, so
    its status is NOT_RUN."""
    row = _recurrence_fixture(fixture_dir, [])

    assert row["independent_cluster_status"] == "NOT_RUN"


def test_database_sources_are_provenance_not_biology(fixture_dir):
    """Source databases are counted per family as provenance."""
    row = _recurrence_fixture(fixture_dir,
                              [["pl1", "L1"], ["pl2", "L1"], ["pl3", "L1"]])

    assert row["database_source_count"] == "2", "PLSDB and IMG were not both counted"


@requires("mmseqs")
def test_family_network_links_a_dark_cluster_to_an_annotated_relative(fixture_dir):
    """Two 50%-identity clusters of related sequences, one annotated and one dark, an
    unrelated dark protein and an unsearched one. The network must link the first two,
    mark the unrelated dark one as unconnected, leave the unsearched one neither dark nor
    bright, and say so in the summary."""
    import random
    rng = random.Random(3)
    aa = "ACDEFGHIKLMNPQRSTVWY"
    base = "".join(rng.choice(aa) for _ in range(180))
    # Resampling 55% of positions leaves ~47% identity: below the 50% clustering bar, so a
    # separate cluster, but still a significant full-length alignment.
    relative = "".join(c if rng.random() > 0.55 else rng.choice(aa) for c in base)
    loner = "".join(rng.choice(aa) for _ in range(180))
    # unsearched is an AntiFam-flagged protein the cascade never searched: its brightness
    # is unknown, so it is neither dark nor bright.
    unsearched = "".join(rng.choice(aa) for _ in range(180))
    # refonly is a cluster of large-plasmid proteins alone: in MMseqs2's files but not in
    # the family table, so it must not become a node.
    reps = fixture_dir / "reps.fasta"
    write_fasta(reps, [("known", base), ("darkrel", relative), ("loner", loner),
                       ("unsearched", unsearched), ("refonly", relative[::-1])])
    clusters = fixture_dir / "clusters.tsv"
    clusters.write_text("known\tknown\ndarkrel\tdarkrel\nloner\tloner\n"
                        "unsearched\tunsearched\nrefonly\trefonly\n")
    prot = fixture_dir / "protein_annotation.tsv"
    write_tsv(prot, ["seq_id", "functional_class", "annot_label", "explained_fraction",
                     "annot_completeness"],
              [["known", "FUNCTIONAL", "Relaxase", 0.95, "FULL"],
               ["darkrel", "NONE", "", 0.0, "NONE"],
               ["loner", "NONE", "", 0.0, "NONE"],
               ["unsearched", "NOT_SEARCHED", "", "", ""]])
    dark_ids = fixture_dir / "dark_ids.txt"
    dark_ids.write_text("darkrel\nloner\n")
    fams = fixture_dir / "protein_families.tsv"
    write_tsv(fams, ["family_id", "family_resolution", "representative", "family_size",
                     "family_class", "scope", "members"],
              [[f"{res}:{m}", res, m, 1, "ORPHAN", scope, m]
               for res in ("broad", "intermediate")
               for m, scope in (("known", "small_only_known"),
                                ("darkrel", "mixed_unknown"),
                                ("loner", "small_only_unknown"),
                                ("unsearched", "small_only_unknown"))])
    pmap = fixture_dir / "map.tsv"
    pmap.write_text("known\tp1|1\ndarkrel\tp2|1\nloner\tp3|1\nunsearched\tp4|1\n")
    out = {k: str(fixture_dir / f"network_{k}.tsv") for k in ("nodes", "edges", "summary")}

    run_script("family_network.py", FakeSnakemake(
        input={"families": str(fams), "reps": str(reps), "clusters": str(clusters),
               "prot": str(prot), "dark_ids": str(dark_ids), "map": str(pmap)},
        output=out,
        params={"network": {"min_cov": 0.5, "max_evalue": 1e-4, "max_out_edges": 4,
                            "dark_brightness": 0.05},
                "primary": "intermediate", "node_resolution": "intermediate", "seed": 1},
        threads=2))

    nodes = {r["node_id"]: r for r in read_tsv(out["nodes"])}
    assert set(nodes) == {"known", "darkrel", "loner", "unsearched"}
    assert (nodes["unsearched"]["brightness"], nodes["unsearched"]["dark"]) == ("", "")
    assert nodes["darkrel"]["scope"] == "mixed_unknown"
    assert nodes["known"]["dark"] == "0" and nodes["known"]["label"] == "Relaxase"
    assert nodes["known"]["family_id"] == "intermediate:known"
    # No member has an informative label, so the node has none.
    assert nodes["darkrel"]["dark"] == "1" and nodes["darkrel"]["label"] == ""
    assert nodes["darkrel"]["degree"] == "1" and nodes["loner"]["degree"] == "0"
    summary = {r["metric"]: r["value"] for r in read_tsv(out["summary"])}
    assert (summary["dark_nodes"], summary["nodes_not_measured"]) == ("2", "1")
    assert summary["dark_nodes_connected"] == "1"
    assert summary["dark_connected_to_bright"] == "1"
    assert summary["dark_family_singletons"] == "2"
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


def _run_rarity(fixture_dir, small_plasmids):
    """One dark family d, on small plasmid S1 and large plasmid L1."""
    recurrence = fixture_dir / "recurrence.tsv"
    write_tsv(recurrence, ["family_id", "independent_plasmid_cluster_count"],
              [["broad:d", 2]])
    fams = fixture_dir / "dark_families.tsv"
    write_tsv(fams, ["family_id", "members", "small_members"], [["broad:d", "d", "d"]])
    mapping = fixture_dir / "protein_map.tsv"
    mapping.write_text("d\tS1|1,L1|4\n")
    small_ids = fixture_dir / "small_plasmids.txt"
    small_ids.write_text("".join(f"{p}\n" for p in small_plasmids))
    out = fixture_dir / "rarefaction.tsv"

    run_script("rarity.py", FakeSnakemake(
        input={"recurrence": str(recurrence), "dark_families": str(fams),
               "map": str(mapping), "small_ids": str(small_ids)},
        output={"rarity": str(fixture_dir / "rarity.tsv"), "rarefaction": str(out)},
        params={"rarity": {"rare_max_lineages": 3, "widely_conserved_min_lineages": 50,
                           "cross_min_hosts": 2, "cross_min_genera": 2,
                           "rarefaction_replicates": 5},
                "seed": 1}))
    return read_tsv(out)


def test_rarefaction_samples_every_small_plasmid(fixture_dir):
    """The x-axis is the small plasmids, with a dark family or without; a dark family's
    copies on large plasmids do not put those plasmids on the axis."""
    final = _run_rarity(fixture_dir, ["S1", "S2", "S3", "S4"])[-1]

    assert final["n_plasmids"] == "4"
    assert final["mean_families"] == "1.0"


def test_family_rarity_writes_the_labels_and_the_thresholds_behind_them(fixture_dir):
    """broad:d is in two lineages: RARE (at most 3), not LINEAGE_SPECIFIC (exactly 1)."""
    _run_rarity(fixture_dir, ["S1"])
    (row,) = read_tsv(fixture_dir / "rarity.tsv")
    assert (row["family_id"], row["rarity_labels"],
            row["independent_plasmid_cluster_count"]) == ("broad:d", "RARE", "2")
    assert (row["rare_max_lineages"], row["widely_conserved_min_lineages"]) == ("3", "50")


def test_an_undefined_saturation_is_not_reported_as_flattened(fixture_dir, capsys):
    """Two small plasmids give two curve points, too few for a saturation value."""
    assert len(_run_rarity(fixture_dir, ["S1", "S2"])) == 2

    log = capsys.readouterr().out
    assert "saturation undefined" in log
    assert "flattened" not in log


def _run_families(fixture_dir, input, output, params, threads=2, classes=None,
                  small=None):
    """S2f then Stage 5, as the workflow runs them.

    functional_class defaults to NONE for the dark ids and FUNCTIONAL for everything else;
    every plasmid in the map is small unless `small` lists them.
    """
    input = dict(input)
    res = sorted(params["clustering"]["resolutions"])
    clusters = [str(fixture_dir / f"families_{r}_cluster.tsv") for r in res]
    run_script("protein_clustering.py", FakeSnakemake(
        input={"faa": input.pop("faa")},
        output={"reps": [str(fixture_dir / f"families_{r}_rep_seq.fasta") for r in res],
                "clusters": clusters},
        params=params, threads=threads))
    dark = set(open(input["dark_ids"]).read().split())
    seq_ids, plasmids = [], set()
    for line in open(input["map"]):
        sid, orfs = line.rstrip("\n").split("\t")
        seq_ids.append(sid)
        plasmids |= {o.rsplit("|", 1)[0] for o in orfs.split(",")}
    prot = fixture_dir / "protein_annotation.tsv"
    write_tsv(prot, ["seq_id", "functional_class"],
              [[sid, (classes or {}).get(sid, "NONE" if sid in dark else "FUNCTIONAL")]
               for sid in seq_ids])
    small_ids = fixture_dir / "small_plasmids.txt"
    small_ids.write_text("".join(f"{p}\n" for p in sorted(plasmids if small is None
                                                             else small)))
    run_script("protein_families.py", FakeSnakemake(
        input={**input, "clusters": clusters, "prot": str(prot),
               "small_ids": str(small_ids)},
        output=output, params=params, threads=threads))


@requires("mmseqs")
def test_a_conserved_protein_on_many_plasmids_is_not_reported_as_a_singleton(fixture_dir):
    """Clustering runs on unique sequences, so a protein identical on four plasmids is one
    member; n_orfs and the plasmid and lineage counts show how widely it is carried."""
    faa = fixture_dir / "unique_proteins.faa"
    write_fasta(faa, [("S1", "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEK"),
                      ("S2", "MQQTTLNRSDEIVWCAPGHKGGAFLNDVWRDNPHLAGCVLLTSDGKLLWQRRD")])
    pmap = fixture_dir / "protein_map.tsv"
    # S1 is one unique sequence carried by four plasmids; S2 by one.
    pmap.write_text("S1\tp1|1,p2|1,p3|1,p4|1\nS2\tp5|1\n")
    registry = fixture_dir / "clonal_registry.tsv"
    write_tsv(registry, ["plasmid_id", "mob_cluster", "species", "hab_top"],
              [["p1", "AA1", "E. coli", "H"], ["p2", "AA2", "E. coli", "H"],
               ["p3", "AA3", "E. coli", "H"], ["p4", "AA1", "E. coli", "H"],
               ["p5", "AA9", "E. coli", "H"]])
    lineage = fixture_dir / "plasmid_lineage.tsv"
    write_tsv(lineage, ["plasmid_id", "plasmid_lineage_cluster"],
              [["p1", "L1"], ["p2", "L2"], ["p3", "L3"], ["p4", "L1"], ["p5", "L9"]])
    dark_ids = fixture_dir / "dark_ids.txt"
    dark_ids.write_text("S1\nS2\n")
    out = fixture_dir / "protein_families.tsv"

    _run_families(fixture_dir,
        input={"faa": str(faa), "map": str(pmap),
               "registry": str(registry), "lineage": str(lineage),
               "dark_ids": str(dark_ids)},
        output={"families": str(out),
                "dark_families": str(fixture_dir / "dark_families.tsv")},
        params={"clustering": {
            "resolutions": {"broad": {"min_seq_id": 0.30, "coverage": 0.50}},
            "primary": "broad", "cov_mode": 0, "cluster_mode": 0}},
        threads=2)

    rows = {r["representative"]: r for r in read_tsv(out)}
    assert "S1" in rows, f"S1 did not survive clustering: {list(rows)}"
    assert rows["S1"]["family_size"] == "1", "S1 is one unique sequence"
    assert rows["S1"]["n_orfs"] == "4", (
        "the ORF count is missing, so a protein on four plasmids is indistinguishable "
        "from one seen once")
    assert rows["S1"]["family_plasmid_count"] == "4"
    assert rows["S1"]["family_MOB_count"] == "3"
    # p1 and p4 are one lineage: four plasmid records, three lineages.
    assert rows["S1"]["family_plasmid_lineage_count"] == "3"


@requires("mmseqs")
def test_protein_families_clusters_annotated_and_dark_together(fixture_dir):
    """A dark and an annotated protein cluster into one family, which is then 50% dark."""

    # Two near-identical proteins that must cluster together, one dark and one annotated,
    # plus an unrelated dark one that must not join them.
    shared = ("MKVLATTLLGAAFAASSALAQKKWLVRNGDTLSGIAQRYGVSVAQLQRWNHLSSDTIHPGQ"
              "KLRVGSDAPQAAPKAEPKVEAKPAAKPVAKPAAKPVAKPAAKPAAKPKAEEKPKAEEK")
    variant = shared.replace("SSDTIHPGQ", "SSDTIHPGK")
    other = ("MPQRSTVWYACDEFGHIKLMNPQRSTVWYACDEFGHIKLMNPQRSTVWYACDEFGHIKLMN"
             "PQRSTVWYACDEFGHIKLMNPQRSTVWYACDEFGHIKLMNPQRSTVWY")

    faa = fixture_dir / "unique_proteins.faa"
    write_fasta(faa, [("p_dark", shared), ("p_annot", variant), ("p_lone", other)])

    dark_ids = fixture_dir / "dark_ids.txt"
    dark_ids.write_text("p_dark\np_lone\n")

    mapping = fixture_dir / "protein_map.tsv"
    mapping.write_text("p_dark\tpl1|1\np_annot\tpl2|1\np_lone\tpl3|1\n")

    registry = fixture_dir / "clonal_registry.tsv"
    write_tsv(registry, ["plasmid_id", "mob_cluster", "species", "topology", "size_bp",
                         "hab_top"],
              [["pl1", "MOB_A", "Escherichia coli", "circular", 5000, "Host-associated"],
               ["pl2", "MOB_B", "Salmonella enterica", "circular", 6000, "Host-associated"],
               ["pl3", "MOB_A", "Escherichia coli", "linear", 7000, "Environmental"]])

    lineage_tsv = fixture_dir / "plasmid_lineage.tsv"
    write_tsv(lineage_tsv, ["plasmid_id", "plasmid_lineage_cluster"],
              [["pl1", "L1"], ["pl2", "L2"], ["pl3", "L1"]])

    families = fixture_dir / "protein_families.tsv"
    dark_families = fixture_dir / "dark_families.tsv"
    _run_families(fixture_dir,
        input={"faa": str(faa), "dark_ids": str(dark_ids),
               "map": str(mapping), "registry": str(registry),
               "lineage": str(lineage_tsv)},
        output={"families": str(families), "dark_families": str(dark_families)},
        params={"clustering": {
            "resolutions": {"broad": {"min_seq_id": 0.3, "coverage": 0.5}},
            "primary": "broad", "cov_mode": 0, "cluster_mode": 0}},
        threads=2)

    rows = read_tsv(families)
    assert rows, "no families written"

    mixed = [r for r in rows if int(r["family_size"]) > 1]
    assert mixed, "the two near-identical proteins did not cluster together"
    row = mixed[0]
    assert int(row["dark_member_count"]) == 1
    assert int(row["annotated_member_count"]) == 1
    assert float(row["percentage_dark_in_family"]) == 50.0
    assert row["dark_only"] == "0"
    assert int(row["family_plasmid_count"]) == 2
    assert int(row["family_host_count"]) == 2
    assert int(row["family_MOB_count"]) == 2


@requires("mmseqs")
def test_a_family_without_a_dark_member_is_not_dark_only(fixture_dir):
    """An AntiFam-flagged protein (NOT_SEARCHED) and one the dark set excludes (NONE, not
    dark) are unnamed but not dark, so neither family is dark-only. No member plasmid is in
    the lineage table, so the lineage count was not measured."""
    faa = fixture_dir / "unique_proteins.faa"
    write_fasta(faa, [("af", "MKVLATTLLGAAFAASSALAQKKWLVRNGDTLSGIAQRYGVSVAQLQRWNH"),
                      ("lc", "MPQRSTVWYACDEFGHIKLMNPQRSTVWYACDEFGHIKLMNPQRSTVWYAC")])
    dark_ids = fixture_dir / "dark_ids.txt"
    dark_ids.write_text("")
    mapping = fixture_dir / "protein_map.tsv"
    mapping.write_text("af\tS1|1\nlc\tS1|2\n")
    registry = fixture_dir / "clonal_registry.tsv"
    write_tsv(registry, ["plasmid_id", "mob_cluster", "species", "hab_top"],
              [["S1", "M1", "E. coli", "H"]])
    lineage_tsv = fixture_dir / "plasmid_lineage.tsv"
    write_tsv(lineage_tsv, ["plasmid_id", "plasmid_lineage_cluster"], [])
    families = fixture_dir / "protein_families.tsv"
    dark_families = fixture_dir / "dark_families.tsv"

    _run_families(fixture_dir,
        input={"faa": str(faa), "dark_ids": str(dark_ids), "map": str(mapping),
               "registry": str(registry), "lineage": str(lineage_tsv)},
        output={"families": str(families), "dark_families": str(dark_families)},
        params={"clustering": {
            "resolutions": {"broad": {"min_seq_id": 0.3, "coverage": 0.5}},
            "primary": "broad", "cov_mode": 0, "cluster_mode": 0}},
        classes={"af": "NOT_SEARCHED", "lc": "NONE"})

    rows = read_tsv(families)
    assert {m for r in rows for m in r["members"].split(",")} == {"af", "lc"}
    for row in rows:
        assert (row["dark_member_count"], row["dark_only"]) == ("0", "0"), row
        assert row["family_plasmid_lineage_status"] == "NOT_RUN"
    assert read_tsv(dark_families) == []


@requires("mmseqs")
def test_protein_families_calls_each_family_small_only_or_mixed_and_known_or_unknown(
        fixture_dir):
    """A family is written when it holds a small-plasmid protein. It is mixed when a member
    also occurs on a large plasmid, and known when any member is named - on either side.
    Counts cover ALL members; the small/large split is in its own columns."""
    import random
    rng = random.Random(7)
    aa = "ACDEFGHIKLMNPQRSTVWY"

    def seq():
        return "M" + "".join(rng.choice(aa) for _ in range(150))

    def variant(s):
        return "".join(c if i % 20 else rng.choice(aa) for i, c in enumerate(s))

    a, b, c, d, f = (seq() for _ in range(5))
    proteins = [("s_known", a), ("s_dark_rel", variant(a)),    # small_only_known
                ("s_alone", b),                               # small_only_unknown
                ("s_dark_c", c), ("l_known", variant(c)),     # mixed_known, known on large
                ("both", d),                                  # on a small AND a large one
                ("l_only", f), ("l_only2", variant(f))]       # no small member: not written
    plasmid = {"s_known": "S1", "s_dark_rel": "S2", "s_alone": "S3", "s_dark_c": "S4",
               "l_known": "L1", "l_only": "L2", "l_only2": "L3"}
    faa = fixture_dir / "unique_proteins.faa"
    write_fasta(faa, proteins)
    dark_ids = fixture_dir / "dark_ids.txt"
    dark_ids.write_text("s_dark_rel\ns_alone\ns_dark_c\nboth\n")
    mapping = fixture_dir / "protein_map.tsv"
    mapping.write_text("".join(f"{sid}\t{plasmid[sid]}|1\n" for sid in plasmid)
                       + "both\tS5|1,L4|1\n")
    plasmids = sorted(set(plasmid.values()) | {"S5", "L4"})
    registry = fixture_dir / "clonal_registry.tsv"
    write_tsv(registry, ["plasmid_id", "mob_cluster", "species", "hab_top"],
              [[p, "M", "E. coli", "H"] for p in plasmids])
    lineage_tsv = fixture_dir / "plasmid_lineage.tsv"
    write_tsv(lineage_tsv, ["plasmid_id", "plasmid_lineage_cluster"],
              [[p, "L"] for p in plasmids])
    families = fixture_dir / "protein_families.tsv"
    dark_families = fixture_dir / "dark_families.tsv"

    _run_families(fixture_dir,
        input={"faa": str(faa), "dark_ids": str(dark_ids), "map": str(mapping),
               "registry": str(registry), "lineage": str(lineage_tsv)},
        output={"families": str(families), "dark_families": str(dark_families)},
        params={"clustering": {
            "resolutions": {"broad": {"min_seq_id": 0.3, "coverage": 0.5}},
            "primary": "broad", "cov_mode": 0, "cluster_mode": 0}},
        classes={"l_only": "NOT_SEARCHED", "l_only2": "NOT_SEARCHED"},
        small=["S1", "S2", "S3", "S4", "S5"])

    by_member = {m: r for r in read_tsv(families) for m in r["members"].split(",")}
    assert "l_only" not in by_member, "a family without a small-plasmid member was written"
    assert by_member["s_known"]["scope"] == "small_only_known"
    assert by_member["s_known"]["known_from"] == "small"
    assert by_member["s_alone"]["scope"] == "small_only_unknown"
    mixed = by_member["s_dark_c"]
    assert (mixed["scope"], mixed["known_from"]) == ("mixed_known", "large")
    # Every member counts; the split is in its own columns.
    assert (mixed["family_size"], mixed["n_small_members"],
            mixed["n_large_members"]) == ("2", "1", "1")
    both = by_member["both"]
    assert both["scope"] == "mixed_unknown" and both["family_class"] == "ORPHAN"
    assert (both["n_small_members"], both["n_large_members"]) == ("1", "1")
    assert {r["scope"] for r in read_tsv(dark_families)} == {
        "small_only_known", "small_only_unknown", "mixed_known", "mixed_unknown"}


@requires("mmseqs")
def test_family_ids_are_content_derived_not_ordinal(fixture_dir):
    """family_id is <resolution>:<representative>, which does not depend on cluster order."""
    faa = fixture_dir / "unique_proteins.faa"
    write_fasta(faa, [("a", "MKVLATTLLGAAFAASSALAQKKWLVRNGDTLSGIAQRYGVSVAQLQRWNH"),
                      ("b", "MPQRSTVWYACDEFGHIKLMNPQRSTVWYACDEFGHIKLMNPQRSTVWYAC")])
    dark_ids = fixture_dir / "dark_ids.txt"
    dark_ids.write_text("a\n")
    mapping = fixture_dir / "protein_map.tsv"
    mapping.write_text("a\tpl1|1\nb\tpl2|1\n")
    registry = fixture_dir / "clonal_registry.tsv"
    write_tsv(registry, ["plasmid_id", "mob_cluster", "species", "topology", "size_bp",
                         "hab_top"], [["pl1", "M1", "E. coli", "circular", 100, "H"],
                                      ["pl2", "M2", "E. coli", "circular", 100, "H"]])

    lineage_tsv = fixture_dir / "plasmid_lineage.tsv"
    write_tsv(lineage_tsv, ["plasmid_id", "plasmid_lineage_cluster"],
              [["pl1", "L1"], ["pl2", "L2"]])

    families = fixture_dir / "protein_families.tsv"
    _run_families(fixture_dir,
        input={"faa": str(faa), "dark_ids": str(dark_ids),
               "map": str(mapping), "registry": str(registry),
               "lineage": str(lineage_tsv)},
        output={"families": str(families),
                "dark_families": str(fixture_dir / "dark_families.tsv")},
        params={"clustering": {
            "resolutions": {"broad": {"min_seq_id": 0.3, "coverage": 0.5}},
            "primary": "broad", "cov_mode": 0, "cluster_mode": 0}},
        threads=2)

    for row in read_tsv(families):
        assert row["family_id"] == f"broad:{row['representative']}", (
            f"family_id {row['family_id']!r} is not <resolution>:<representative>")


@requires("mmseqs")
def test_the_dark_family_representative_is_a_dark_protein(fixture_dir):
    """S8d searches the dark family's representative structurally, so it must be a dark
    member even when MMseqs2 picks an annotated one, and only dark members are listed."""
    shared = ("MKVLATTLLGAAFAASSALAQKKWLVRNGDTLSGIAQRYGVSVAQLQRWNHLSSDTIHPGQ"
              "KLRVGSDAPQAAPKAEPKVEAKPAAKPVAKPAAKPVAKPAAKPAAKPKAEEKPKAEEK")
    faa = fixture_dir / "unique_proteins.faa"
    # The longer sequence is the one MMseqs2 tends to pick as representative; make it the
    # ANNOTATED member so the test fails if the representative is taken unconditionally.
    write_fasta(faa, [("p_annot_long", shared + "AAAKPAAKPAAKPAAKPKAEEK"),
                      ("p_dark", shared)])
    dark_ids = fixture_dir / "dark_ids.txt"
    dark_ids.write_text("p_dark\n")
    mapping = fixture_dir / "protein_map.tsv"
    mapping.write_text("p_annot_long\tpl1|1\np_dark\tpl2|1\n")
    registry = fixture_dir / "clonal_registry.tsv"
    write_tsv(registry, ["plasmid_id", "mob_cluster", "species", "topology", "size_bp",
                         "hab_top"], [["pl1", "M1", "E. coli", "circular", 100, "H"],
                                      ["pl2", "M2", "E. coli", "circular", 100, "H"]])

    lineage_tsv = fixture_dir / "plasmid_lineage.tsv"
    write_tsv(lineage_tsv, ["plasmid_id", "plasmid_lineage_cluster"],
              [["pl1", "L1"], ["pl2", "L2"]])

    dark_families = fixture_dir / "dark_families.tsv"
    _run_families(fixture_dir,
        input={"faa": str(faa), "dark_ids": str(dark_ids),
               "map": str(mapping), "registry": str(registry),
               "lineage": str(lineage_tsv)},
        output={"families": str(fixture_dir / "protein_families.tsv"),
                "dark_families": str(dark_families)},
        params={"clustering": {
            "resolutions": {"broad": {"min_seq_id": 0.3, "coverage": 0.5}},
            "primary": "broad", "cov_mode": 0, "cluster_mode": 0}},
        threads=2)

    rows = read_tsv(dark_families)
    assert rows, "no dark families derived"
    for row in rows:
        assert row["representative"] == "p_dark", (
            f"the dark family's representative is {row['representative']!r}, which is not "
            "a dark protein - S8d would search the wrong sequence")
        assert "p_annot_long" not in row["members"], (
            "an annotated member leaked into the dark family's member list, which would "
            "widen every downstream evolution and context measurement")
