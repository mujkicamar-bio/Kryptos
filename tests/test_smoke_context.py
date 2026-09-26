"""Smoke tests: synteny, evolution, context features, structure, defence, mobile elements.

Each test runs one workflow script against a small fixture.
"""
import pathlib

import pytest
from conftest import (
    FakeSnakemake,
    _is_table,
    read_tsv,
    requires,
    run_script,
    write_fasta,
    write_tsv,
)


@requires("mafft")
def test_family_evolution_runs_and_reports_a_dnds_status(fixture_dir):
    """S7b crashed with NameError on the first family with >=3 members. It survived
    review because without mafft on PATH the script takes the ALIGNMENT_FAILED branch and
    exits 0 - so the bug only appears when the tool is actually present."""
    faa = fixture_dir / "dark.faa"
    cds = fixture_dir / "dark.fna"
    fams = fixture_dir / "families.tsv"
    out = fixture_dir / "family_evolution.tsv"

    # Three real, alignable members of one family, divergent only at silent sites.
    prot = "MKVLATTLLGAAFAASSALAQ"
    codons = {"M": "ATG", "K": "AAA", "V": "GTG", "L": "CTG", "A": "GCG", "T": "ACC",
              "G": "GGC", "F": "TTT", "S": "AGC", "Q": "CAG"}
    alt = {"L": "CTA", "A": "GCA", "T": "ACT", "S": "TCT", "G": "GGT"}
    base = "".join(codons[a] for a in prot)
    var1 = "".join(alt.get(a, codons[a]) if i % 4 == 0 else codons[a]
                   for i, a in enumerate(prot))
    var2 = "".join(alt.get(a, codons[a]) if i % 5 == 0 else codons[a]
                   for i, a in enumerate(prot))

    write_fasta(faa, [("m1", prot), ("m2", prot), ("m3", prot)])
    write_fasta(cds, [("m1", base), ("m2", var1), ("m3", var2)])
    write_tsv(fams, ["family_id", "representative", "n_members", "n_plasmids",
                     "n_mob_clusters", "family_class", "members", "small_members"],
              [["F0000001", "m1", 3, 3, 2, "FAMILY", "m1,m2,m3", "m1,m2"]])

    run_script("family_evolution.py", FakeSnakemake(
        input={"families": str(fams), "faa": str(faa), "cds": str(cds)},
        output={"tsv": str(out), "consensus": str(out.parent / "consensus.faa")},
        params={"evolution": {"min_codons": 20, "min_members_for_dnds": 3,
                              "dnds_purifying_max": 0.5, "rnacode_max_p": 0.05,
                              "max_members_aligned": 50}},
        threads=2))

    rows = read_tsv(out)
    assert len(rows) == 1
    # The small-plasmid members are measured on their own: two are too few for dN/dS.
    assert rows[0]["small_dnds_status"] == "TOO_FEW_MEMBERS"
    assert rows[0]["small_n_aligned"] == "0"
    # The family is alignable and divergent, so a status must have been reached. An empty
    # status means the script fell through an error branch and reported success.
    assert rows[0]["dnds_status"], "no dnds_status - the script took a silent error branch"
    assert rows[0]["dnds_status"] != "ALIGNMENT_FAILED", (
        "alignment failed with mafft available - the tool is not being invoked correctly")


def _run_context(fixture_dir, is_rows=(), genes=None, topology="linear"):
    """One plasmid: dark ORF pl1|1 in a two-gene directon with an annotated partner, and a
    dark ORF pl1|3 inside an integron cassette array. F1 is pl1|1's family."""
    ann = fixture_dir / "plasmid_annotation.tsv"
    write_tsv(ann, ["plasmid_id", "orf_id", "start", "end", "strand", "annot_label",
                    "functional_class"],
              genes or [["pl1", "pl1|1", 100, 400, "+", "", "NONE"],
                        ["pl1", "pl1|2", 430, 700, "+", "MobA_MobL", "FUNCTIONAL"],
                        ["pl1", "pl1|3", 3000, 3300, "+", "", "NONE"]])
    fam = fixture_dir / "dark_families.tsv"
    write_tsv(fam, ["family_id", "representative", "n_members", "n_plasmids",
                    "n_mob_clusters", "family_class", "members"],
              [["F1", "S1", 1, 1, 1, "ORPHAN", "S1"]])
    pmap = fixture_dir / "protein_map.tsv"
    pmap.write_text("S1\tpl1|1\nS2\tpl1|2\nS3\tpl1|3\n")
    defence = fixture_dir / "defence_systems.tsv"
    write_tsv(defence, ["orf_id", "system"], [])
    integrons = fixture_dir / "integrons.tsv"
    write_tsv(integrons, ["plasmid_id", "integron_id", "element", "start", "end",
                          "integron_type", "annotation", "type_elt"],
              [["pl1", "in1", "protein", 3000, 3300, "complete", "protein", "protein"]])
    master = fixture_dir / "context_master.tsv"
    write_tsv(master, ["plasmid_id", "size_bp", "topology"], [["pl1", 5000, topology]])
    lengths = fixture_dir / "plasmid_lengths.tsv"
    write_tsv(lengths, ["plasmid_id", "length_bp"], [["pl1", 5000]])
    conj = fixture_dir / "conjugation_systems.tsv"
    write_tsv(conj, ["orf_id", "plasmid_id", "system", "system_id", "component"], [])
    labels = fixture_dir / "protein_labels.tsv"
    write_tsv(labels, ["protein_id", "source", "tier", "kind", "label", "sub_label"], [])
    all_fams = fixture_dir / "protein_families.tsv"
    write_tsv(all_fams, ["family_id", "family_resolution", "representative", "members"],
              [["F1", "intermediate", "S1", "S1"]])
    lineage = fixture_dir / "plasmid_lineage.tsv"
    write_tsv(lineage, ["plasmid_id", "plasmid_lineage_cluster"], [["pl1", "pl1"]])
    fams_out = fixture_dir / "family_context.tsv"
    run_script("context_features.py", FakeSnakemake(
        input={"annotation": str(ann), "families": str(fam), "map": str(pmap),
               "defence": str(defence), "conjugation": str(conj),
               "integrons": str(integrons),
               "is_elements": _is_table(fixture_dir, is_rows), "master": str(master),
               "lengths": str(lengths), "labels": str(labels),
               "all_families": str(all_fams), "lineage": str(lineage)},
        output={"families": str(fams_out),
                "terms": str(fixture_dir / "family_context_terms.tsv")},
        params={"context": {"max_operon_gap": 100, "neighbourhood_window": 3},
                "primary": "intermediate"}))
    return read_tsv(fams_out)


def test_context_writes_one_row_of_rates_per_family(fixture_dir):
    """The unit is the plasmid, and an absent feature is a rate of 0, not a blank."""
    rows = _run_context(fixture_dir)

    assert len(rows) == 1
    row = rows[0]
    assert row["family_id"] == "F1" and row["n_units"] == "1"
    assert row["cons_annotated_neighbour"] == "1.0"
    assert row["cons_operon_with_annotated"] == "1.0"
    assert row["cons_two_gene_operon"] == "1.0"
    # pl1|3 is in the integron, not the family's ORF.
    assert row["cons_integron"] == "0.0"
    assert row["cons_defence"] == "0.0"
    assert row["cons_is_element"] == "0.0"
    assert row["cons_conj"] == "0.0"


FOLDSEEK_DB = "data/refs/foldseek/pdb"


PROSTT5 = "data/refs/foldseek/prostt5"


@pytest.mark.slow
@requires("foldseek")
@pytest.mark.skipif(not pathlib.Path(PROSTT5).exists(),
                    reason="ProstT5 model not downloaded")
def test_structure_search_reports_what_the_match_actually_is(fixture_dir):
    """Foldseek's `target` is a PDB accession - `12as-assembly1_A`. Every downstream use
    of structural evidence needs to know WHAT the fold is, not which entry it came from:
    the nucleic_acid_binding stratum tested for the word "nucle" in the accession and so
    could never be filled, and the hypothesis written into the synthesis order read
    `structural:12as-assembly1_A`, which tells a bench scientist nothing.

    Takes about 35 s because ProstT5 has to load; marked slow."""
    faa = fixture_dir / "dark_proteins.faa"
    write_fasta(faa, [("q1",
                       "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVK"
                       "ALPDAQFEVVHSLAKWKRQTLGQHDFSAGEGLYTHMKALRPDEDRLSPLHSVYVDQWDWE"
                       "RVMGDGERQFSTLKSTVEAIWAGIKATEAAVSEEFGLAPFLPDQIHFVHSQELLSRYPDL"
                       "DAKGRERAIAKDLGAVFLVGIGGKLSDGHRHDVRAPDYDDWSTPSELGHAGLNGDILVWN"
                       "PVLEDAFELSSMGIRVDADTLKHQLALTGDEDRLELEWHQALLRGEMPQTIGGGIGQSRL"
                       "TMLLLQLPHIGQVQAGVWPAAVRESVPSLL")])
    out = fixture_dir / "structure_hits.tsv"

    # scope: all, so the test searches the query it wrote rather than resolving family
    # representatives. The representative path is covered by its own test below.
    run_script("structure_search.py", FakeSnakemake(
        input={"faa": str(faa), "families": str(fixture_dir / "unused_families.tsv")},
        output=[str(out)],
        params={"structure": {"max_evalue": 1.0e-3, "min_plddt": 70, "required": True,
                              "scope": "all"},
                "target_db": FOLDSEEK_DB, "prostt5": PROSTT5},
        threads=4))

    rows = read_tsv(out)
    assert rows, "foldseek returned no structural match for a protein with an exact one"
    assert rows[0]["target_description"], (
        "no description recorded - only the accession, which names nothing")
    assert "SYNTHETASE" in rows[0]["target_description"].upper()


@requires("mafft", "RNAcode")
def test_family_evolution_reports_coding_potential(fixture_dir):
    """The design names RNAcode as one of two positive-evidence tests at S7: coding signal
    independent of the gene caller. `evolution.rnacode_max_p` was declared, the binary was
    installed, and no code invoked it - the columns were written empty for every family and
    the run reported success.

    RNAcode scores BOTH strands, which matters more here than anywhere: a shadow ORF is the
    reverse complement of a real gene, so its antisense signal should beat its sense signal.
    Recording only the sense P would throw away the one number that distinguishes the
    artefact class this stage exists to catch."""
    prot = "MKVLATTLLGAAFAASSALAQKKWLVRDGHIYQPLMNEATSGKLW"
    cod = {"M": "ATG", "K": "AAA", "V": "GTG", "L": "CTG", "A": "GCG", "T": "ACC",
           "G": "GGC", "F": "TTT", "S": "AGC", "Q": "CAG", "W": "TGG", "R": "CGT",
           "D": "GAT", "H": "CAT", "I": "ATT", "Y": "TAT", "P": "CCG", "N": "AAC",
           "E": "GAA"}
    alt = {"L": "CTA", "A": "GCA", "T": "ACT", "S": "TCT", "G": "GGT", "V": "GTA",
           "R": "CGC", "P": "CCA"}
    members = []
    for k in range(4):
        members.append((f"m{k}", "".join(
            alt.get(a, cod[a]) if (i + k) % 5 == 0 else cod[a]
            for i, a in enumerate(prot))))

    faa = fixture_dir / "dark.faa"
    cds = fixture_dir / "dark.fna"
    write_fasta(faa, [(n, prot) for n, _ in members])
    write_fasta(cds, members)
    fams = fixture_dir / "families.tsv"
    write_tsv(fams, ["family_id", "representative", "n_members", "n_orfs", "n_plasmids",
                     "n_mob_clusters", "family_class", "members"],
              [["F1", "m0", 4, 4, 4, 3, "FAMILY", ",".join(n for n, _ in members)]])
    out = fixture_dir / "family_evolution.tsv"

    run_script("family_evolution.py", FakeSnakemake(
        input={"families": str(fams), "faa": str(faa), "cds": str(cds)},
        output={"tsv": str(out), "consensus": str(out.parent / "consensus.faa")},
        params={"evolution": {"min_codons": 20, "min_members_for_dnds": 3,
                              "dnds_purifying_max": 0.5, "rnacode_max_p": 0.05,
                              "max_members_aligned": 50}},
        threads=1))

    row = read_tsv(out)[0]
    assert row["rnacode_p"], "rnacode_p is still empty - RNAcode was never invoked"
    assert float(row["rnacode_p"]) < 0.05
    assert row["coding_signal"] == "1"
    assert row["rnacode_p_antisense"], (
        "the antisense P is not recorded, so a shadow ORF cannot be distinguished by the "
        "one test designed to catch it")
    assert row["rnacode_status"] == "MEASURED"


@requires("mafft", "RNAcode")
def test_an_rnacode_failure_is_not_silently_a_pass(fixture_dir):
    """RNAcode prints `ERROR: Unknown alignment file format` to stdout and exits 0. Parsing
    its output without checking would read zero rows as "no coding signal" - which is
    evidence AGAINST a family - when in fact the tool never ran."""
    faa = fixture_dir / "dark.faa"
    cds = fixture_dir / "dark.fna"
    # Two members: below min_members_for_dnds, so no alignment is attempted at all.
    write_fasta(faa, [("m0", "MKVL"), ("m1", "MKVL")])
    write_fasta(cds, [("m0", "ATGAAAGTGCTG"), ("m1", "ATGAAAGTACTG")])
    fams = fixture_dir / "families.tsv"
    write_tsv(fams, ["family_id", "representative", "n_members", "n_orfs", "n_plasmids",
                     "n_mob_clusters", "family_class", "members"],
              [["F1", "m0", 2, 2, 2, 1, "FAMILY", "m0,m1"]])
    out = fixture_dir / "family_evolution.tsv"

    run_script("family_evolution.py", FakeSnakemake(
        input={"families": str(fams), "faa": str(faa), "cds": str(cds)},
        output={"tsv": str(out), "consensus": str(out.parent / "consensus.faa")},
        params={"evolution": {"min_codons": 20, "min_members_for_dnds": 3,
                              "dnds_purifying_max": 0.5, "rnacode_max_p": 0.05,
                              "max_members_aligned": 50}},
        threads=1))

    row = read_tsv(out)[0]
    assert row["coding_signal"] == "", "absence of a test must not read as a failed test"
    assert row["rnacode_status"] == "TOO_FEW_MEMBERS"


@requires("mafft")
def test_family_evolution_writes_a_consensus_per_family(fixture_dir):
    """The re-check needs one sequence per family carrying the family's shared signal.
    S7b already has the protein alignment in hand, so building it here costs nothing and
    avoids aligning every family a second time."""
    prot_a = "MKVLATTLLGAAFAASSALAQKKWLVRDGHIY"
    prot_b = "MKVLATTLLGAAFCASSALAQKKWLVRDGHIY"
    faa = fixture_dir / "dark.faa"
    cds = fixture_dir / "dark.fna"
    cod = {"M": "ATG", "K": "AAA", "V": "GTG", "L": "CTG", "A": "GCG", "T": "ACC",
           "G": "GGC", "F": "TTT", "S": "AGC", "Q": "CAG", "W": "TGG", "R": "CGT",
           "D": "GAT", "H": "CAT", "I": "ATT", "Y": "TAT", "C": "TGC"}
    write_fasta(faa, [("m1", prot_a), ("m2", prot_a), ("m3", prot_b)])
    write_fasta(cds, [(n, "".join(cod[a] for a in p))
                      for n, p in (("m1", prot_a), ("m2", prot_a), ("m3", prot_b))])
    fams = fixture_dir / "families.tsv"
    write_tsv(fams, ["family_id", "representative", "n_members", "n_orfs", "n_plasmids",
                     "n_mob_clusters", "family_class", "members"],
              [["F1", "m1", 3, 3, 3, 2, "FAMILY", "m1,m2,m3"]])
    out = fixture_dir / "family_evolution.tsv"
    cons = fixture_dir / "family_consensus.faa"

    run_script("family_evolution.py", FakeSnakemake(
        input={"families": str(fams), "faa": str(faa), "cds": str(cds)},
        output={"tsv": str(out), "consensus": str(cons)},
        params={"evolution": {"min_codons": 20, "min_members_for_dnds": 3,
                              "dnds_purifying_max": 0.5, "rnacode_max_p": 0.05,
                              "max_members_aligned": 50}},
        threads=1))

    text = cons.read_text()
    assert text.startswith(">F1"), f"no consensus written: {text[:80]!r}"
    seq = "".join(l for l in text.splitlines() if not l.startswith(">"))
    assert seq == prot_a, f"consensus is not the majority sequence: {seq}"


@requires("hmmsearch")
@pytest.mark.skipif(not pathlib.Path("data/refs/pfam/Pfam-A.hmm").exists(),
                    reason="Pfam-A not downloaded")
@pytest.mark.slow
def test_the_consensus_recheck_finds_a_family_that_is_collectively_recognisable(fixture_dir):
    """A family whose consensus hits Pfam is NOT collectively novel, however dark each
    member looked on its own. Pavlopoulos removed 6.5% of clusters this way, and those are
    exactly the clusters this pipeline would otherwise send to the bench as novel.

    This is a LABEL, not a filter - the family stays in the table with the verdict on it."""
    import subprocess
    # A real Pfam consensus sequence, so the re-check has something true to find.
    subprocess.run("hmmfetch --index data/refs/pfam/Pfam-A.hmm", shell=True,
                   capture_output=True)
    fetched = subprocess.run("hmmfetch data/refs/pfam/Pfam-A.hmm RepA_N", shell=True,
                             capture_output=True, text=True)
    hmm = fixture_dir / "one.hmm"
    hmm.write_text(fetched.stdout)
    emitted = subprocess.run(f"hmmemit -c {hmm}", shell=True, capture_output=True,
                             text=True).stdout
    repa = "".join(l.strip() for l in emitted.splitlines() if not l.startswith(">"))

    cons = fixture_dir / "family_consensus.faa"
    write_fasta(cons, [("F1", repa), ("F2", "MKEEDDKKEEDDKKEEDDKKEEQ" * 3)])
    out = fixture_dir / "consensus_recheck.tsv"
    fams = fixture_dir / "dark_families.tsv"
    write_tsv(fams, ["family_id"], [["F1"], ["F2"], ["F3"]])

    run_script("consensus_recheck.py", FakeSnakemake(
        input={"consensus": str(cons), "families": str(fams)},
        output=[str(out)],
        params={"db": "data/refs/pfam/Pfam-A.hmm", "args": "--cut_ga",
                "hmmer_z": 3497616},
        threads=4))

    rows = {r["family_id"]: r for r in read_tsv(out)}
    assert rows["F1"]["consensus_hit"] == "1", (
        "a family whose consensus is a Pfam family consensus was called novel")
    assert rows["F1"]["consensus_label"] == "RepA_N"
    assert rows["F2"]["consensus_hit"] == "0"
    assert rows["F2"]["collectively_novel"] == "1"
    # F3 had too few members for a consensus: in the table, and saying it was not tested.
    assert rows["F3"]["consensus_status"] == "NOT_RUN"
    assert rows["F3"]["collectively_novel"] == ""


def test_defence_systems_keeps_component_status_and_system_wholeness(fixture_dir):
    """A mandatory component of a complete system and a neutral component of a fragment
    are not the same evidence, and both arrived as the same row. MacSyFinder reports the
    distinction and it was discarded."""
    out = fixture_dir / "defence_systems.tsv"
    phase2 = out.parent / "phase2" / "run"
    phase2.mkdir(parents=True)
    write_tsv(phase2 / "best_solution.tsv",
              ["replicon", "hit_id", "gene_name", "hit_pos", "model_fqn", "sys_id",
               "sys_wholeness", "sys_score", "hit_gene_ref", "hit_status",
               "hit_i_eval", "hit_profile_cov"],
              [["p1", "GB1", "RM_Type_II_REase", 3,
                "defense-finder-models/Defense/RM_Type_II", "p1_RM_1",
                "1.000", "5.5", "RM_Type_II_REase", "mandatory", "1e-40", "0.95"],
               ["p1", "GB2", "RM_Type_II_MTase", 4,
                "defense-finder-models/Defense/RM_Type_II", "p1_RM_1",
                "1.000", "5.5", "RM_Type_II_MTase", "accessory", "1e-20", "0.60"]])

    mapping = fixture_dir / "map.tsv"
    write_tsv(mapping, ["gembase_id", "orf_id", "plasmid_id"],
              [["GB1", "p1|1", "p1"], ["GB2", "p1|2", "p1"]])
    faa = fixture_dir / "cand.faa"
    write_fasta(faa, [("GB1", "MKV"), ("GB2", "MKW")])

    # macsyfinder is not invoked: its output tree is pre-populated above, which is what the
    # parser under test reads. skip_run lets this test run without the model set installed.
    models = fixture_dir / "models"
    models.mkdir()
    (models / "placeholder").write_text("")

    run_script("defence_systems.py", FakeSnakemake(
        input={"faa": str(faa), "map": str(mapping)},
        output={"tsv": str(out)},
        params={"models_dir": str(models), "skip_run": True, "required": False},
        threads=1))

    rows = {r["orf_id"]: r for r in read_tsv(out)}
    assert rows["p1|1"]["hit_status"] == "mandatory"
    assert rows["p1|2"]["hit_status"] == "accessory"
    assert rows["p1|1"]["sys_wholeness"] == "1.000"
    assert rows["p1|1"]["hit_profile_cov"] == "0.95"
    assert rows["p1|1"]["system"] == "defense-finder-models/Defense/RM_Type_II"


def test_defence_records_not_run_when_the_models_are_absent(fixture_dir):
    """Spec section 7.2 separates NOT_RUN from NO_HIT, and this is where the distinction
    is earned. An empty defence table with no status reads as 'this collection carries no
    defence systems', which is a biological claim a run without the models has not made.

    The models are an OPTIONAL database: the primary deliverable is complete annotation of
    every ORF, which does not depend on them, so their absence must not halt the run."""
    faa = fixture_dir / "cand.faa"
    write_fasta(faa, [("GB1", "MKV")])
    mapping = fixture_dir / "map.tsv"
    write_tsv(mapping, ["gembase_id", "orf_id", "plasmid_id"], [["GB1", "p1|1", "p1"]])
    out = fixture_dir / "defence_systems.tsv"

    # The script ends with sys.exit(0) - a successful early return for a Snakemake script,
    # which the in-process harness sees as SystemExit. The code is asserted rather than
    # swallowed: exit 0 is the whole claim being made, that this is a clean skip and not a
    # failure.
    with pytest.raises(SystemExit) as exit_info:
        run_script("defence_systems.py", FakeSnakemake(
            input={"faa": str(faa), "map": str(mapping)},
            output={"tsv": str(out)},
            params={"models_dir": str(fixture_dir / "absent"), "required": False},
            threads=1))
    assert exit_info.value.code == 0, "a missing optional database exited non-zero"

    assert out.exists(), "no table written, so downstream stages cannot read the status"
    rows = read_tsv(out)
    assert rows == [], "rows were invented for a search that never ran"
    assert "status" in out.read_text().split("\n")[0], (
        "the table carries no status column, so absent-because-not-searched cannot be "
        "told from absent-because-searched")


def test_defence_halts_when_the_models_are_required_and_absent(fixture_dir):
    """The other direction: in a production run a silently missing defence axis is a
    defect, so required true must fail loudly rather than record NOT_RUN."""
    faa = fixture_dir / "cand.faa"
    write_fasta(faa, [("GB1", "MKV")])
    mapping = fixture_dir / "map.tsv"
    write_tsv(mapping, ["gembase_id", "orf_id", "plasmid_id"], [["GB1", "p1|1", "p1"]])

    with pytest.raises(SystemExit, match="not installed"):
        run_script("defence_systems.py", FakeSnakemake(
            input={"faa": str(faa), "map": str(mapping)},
            output={"tsv": str(fixture_dir / "out.tsv")},
            params={"models_dir": str(fixture_dir / "absent"), "required": True},
            threads=1))


def test_structure_search_restricts_to_family_representatives(fixture_dir):
    """Spec section 49 sets representative scale as the discovery-scale strategy and
    section 79 makes it a success criterion. ProstT5 is a transformer and the query count
    is the cost of the stage, so searching every dark protein rather than one per family is
    the difference between the largest job in the pipeline and a modest one.

    Foldseek is not invoked here: the assertion is on which sequences reach the query file,
    which is what the scope setting controls."""
    faa = fixture_dir / "dark_proteins.faa"
    write_fasta(faa, [("rep_a", "MKTAYIAKQRQISFVKSHFSRQ"),
                      ("member_a", "MKTAYIAKQRQISFVKSHFSRK"),
                      ("rep_b", "MQQTTLNRSDEIVWCAPGHKGG")])
    families = fixture_dir / "dark_families.tsv"
    write_tsv(families, ["family_id", "representative", "members"],
              [["broad:rep_a", "rep_a", "rep_a,member_a"],
               ["broad:rep_b", "rep_b", "rep_b"]])

    out = fixture_dir / "structure_hits.tsv"
    # foldseek will fail on the absent database; the query file is written before that, and
    # it is the only thing under test.
    try:
        run_script("structure_search.py", FakeSnakemake(
            input={"faa": str(faa), "families": str(families)},
            output=[str(out)],
            params={"structure": {"max_evalue": 1.0e-3, "scope": "representatives"},
                    "target_db": str(fixture_dir / "absent_db"),
                    "prostt5": str(fixture_dir / "absent_model")},
            threads=1))
    except Exception:
        pass

    query = fixture_dir / "structure_query.faa"
    assert query.exists(), "no representative query file was written"
    names = {l[1:].split()[0] for l in query.read_text().splitlines() if l.startswith(">")}
    assert names == {"rep_a", "rep_b"}, (
        f"the query set is {sorted(names)}; it must be one sequence per family, and "
        "member_a is a family member rather than a representative")


def test_an_operon_across_the_origin_of_a_circular_plasmid_counts(fixture_dir):
    """pl1|1 sits at the start of the record and its annotated partner pl1|2 at the end;
    on the 5,000 bp circle they are 49 nt apart and one operon."""
    genes = [["pl1", "pl1|1", 50, 400, "+", "", "NONE"],
             ["pl1", "pl1|3", 2000, 2300, "-", "", "NONE"],
             ["pl1", "pl1|2", 4000, 4950, "+", "MobA_MobL", "FUNCTIONAL"]]

    circular = _run_context(fixture_dir, genes=genes, topology="circular")[0]
    linear = _run_context(fixture_dir, genes=genes, topology="linear")[0]

    assert circular["cons_operon_with_annotated"] == "1.0"
    assert linear["cons_operon_with_annotated"] == "0.0"


def test_an_orf_inside_an_is_element_gets_is_element_context(fixture_dir):
    """An IS element is an island, like a defence system or an integron."""
    rows = _run_context(
        fixture_dir,
        is_rows=[["pl1", "pl1|IS1", "IS3", "IS3_1", 50, 450, "+", 1, "1e-50", ""]])
    assert rows[0]["cons_is_element"] == "1.0"
