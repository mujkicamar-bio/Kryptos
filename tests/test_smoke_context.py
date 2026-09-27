"""Smoke tests: synteny, evolution, context features, structure, defence, mobile elements.

Each test runs one workflow script against a small fixture.
"""
import csv
import gzip
import os
import pathlib
import sys

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


def _evolution_fixture(fixture_dir):
    """Three alignable members of one family of 21 codons, divergent only at silent sites;
    two of them are on small plasmids."""
    faa = fixture_dir / "dark.faa"
    cds = fixture_dir / "dark.fna"
    fams = fixture_dir / "families.tsv"
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
    return {"families": str(fams), "faa": str(faa), "cds": str(cds)}


def _run_evolution(fixture_dir, min_codons):
    out = fixture_dir / "family_evolution.tsv"
    pairs = fixture_dir / "family_yn00_pairs.tsv.gz"
    run_script("family_evolution.py", FakeSnakemake(
        input=_evolution_fixture(fixture_dir),
        output={"tsv": str(out), "consensus": str(fixture_dir / "consensus.faa"),
                "pairs": str(pairs)},
        params={"evolution": {"min_codons": min_codons, "min_members_for_dnds": 3,
                              "dnds_purifying_max": 0.5, "rnacode_max_p": 0.05,
                              "max_members_aligned": 50}},
        threads=2))
    with gzip.open(pairs, "rt") as fh:
        pair_rows = list(csv.DictReader(fh, delimiter="\t"))
    return read_tsv(out), pair_rows


@requires("mafft", "yn00")
def test_family_evolution_summarises_the_yn00_pairs(fixture_dir):
    """The family median is the median yn00 omega over its three pairs, every pair row is
    written, and the small-plasmid member set is measured on its own."""
    rows, pairs = _run_evolution(fixture_dir, min_codons=20)
    assert len(rows) == 1
    assert rows[0]["dnds_status"] == "MEASURED"
    assert rows[0]["n_pairs"] == "3"
    omegas = sorted(float(p["omega"]) for p in pairs)
    assert [p["member_set"] for p in pairs] == ["all"] * 3
    assert float(rows[0]["dnds_median"]) == pytest.approx(omegas[1], abs=1e-4)
    assert float(rows[0]["dnds_min"]) == pytest.approx(omegas[0], abs=1e-4)
    # Two small-plasmid members are too few for dN/dS.
    assert rows[0]["small_dnds_status"] == "TOO_FEW_MEMBERS"
    assert rows[0]["small_n_aligned"] == "0"


@requires("mafft", "yn00")
def test_an_alignment_shorter_than_min_codons_is_too_short(fixture_dir):
    """min_codons applies to the 21 codons yn00 used; the pairs are still written."""
    rows, pairs = _run_evolution(fixture_dir, min_codons=22)
    assert rows[0]["dnds_status"] == "TOO_SHORT"
    assert rows[0]["dnds_median"] == ""
    assert len(pairs) == 3


def _run_context(fixture_dir, is_rows=(), genes=None, topology="linear",
                 defence_rows=(), conj_rows=(), integron_rows=None):
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
    write_tsv(defence, ["orf_id", "system", "status"], list(defence_rows))
    integrons = fixture_dir / "integrons.tsv"
    write_tsv(integrons, ["plasmid_id", "integron_id", "element", "start", "end",
                          "integron_type", "annotation", "type_elt"],
              integron_rows or
              [["pl1", "in1", "protein", 3000, 3300, "complete", "protein", "protein"]])
    master = fixture_dir / "context_master.tsv"
    write_tsv(master, ["plasmid_id", "size_bp", "topology"], [["pl1", 5000, topology]])
    lengths = fixture_dir / "plasmid_lengths.tsv"
    write_tsv(lengths, ["plasmid_id", "length_bp"], [["pl1", 5000]])
    conj = fixture_dir / "conjugation_systems.tsv"
    write_tsv(conj, ["orf_id", "plasmid_id", "system", "system_id", "component", "status"],
              list(conj_rows))
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


def test_a_not_run_system_stage_gives_an_empty_rate_not_zero(fixture_dir):
    """NOT_RUN means nothing was searched, so the rate is missing, not a measured 0."""
    row = _run_context(fixture_dir, defence_rows=[["", "", "NOT_RUN"]],
                       conj_rows=[["", "", "", "", "", "NOT_RUN"]])[0]
    assert row["cons_defence"] == "" and row["cons_conj"] == ""
    assert row["cons_integron"] == "0.0"


def test_only_an_integron_with_a_cassette_array_is_an_island(fixture_dir):
    """pl1|1 (100..400) lies in an integron element: counted for complete and CALIN
    integrons, not for In0, an integrase without attC sites."""
    def rate(integron_type):
        row = ["pl1", "in1", "intI_1", 50, 450, integron_type, "intI", "protein"]
        return _run_context(fixture_dir, integron_rows=[row])[0]["cons_integron"]
    assert (rate("complete"), rate("CALIN"), rate("In0")) == ("1.0", "1.0", "0.0")


def test_an_orf_that_is_a_defence_component_gets_defence_context(fixture_dir):
    row = _run_context(fixture_dir, defence_rows=[["pl1|1", "Clover", "SUCCESS"]])[0]
    assert row["cons_defence"] == "1.0"


def test_an_orf_that_overlaps_a_component_without_being_one_has_no_system_context(
        fixture_dir):
    """pl1|1 overlaps the defence and conjugation component pl1|2 by four bases."""
    genes = [["pl1", "pl1|1", 100, 400, "+", "", "NONE"],
             ["pl1", "pl1|2", 397, 700, "+", "MobA_MobL", "FUNCTIONAL"]]
    row = _run_context(fixture_dir, genes=genes,
                       defence_rows=[["pl1|2", "Clover", "SUCCESS"]],
                       conj_rows=[["pl1|2", "pl1", "T4SS_typeF", "s1", "MOBF", "SUCCESS"]])[0]
    assert row["cons_defence"] == "0.0" and row["cons_conj"] == "0.0"


FOLDSEEK_DB = "data/refs/foldseek/pdb"


PROSTT5 = "data/refs/foldseek/prostt5"


@pytest.mark.slow
@requires("foldseek")
@pytest.mark.skipif(not pathlib.Path(PROSTT5).exists(),
                    reason="ProstT5 model not downloaded")
def test_structure_search_reports_what_the_match_actually_is(fixture_dir):
    """Foldseek's `target` is a PDB accession (12as-assembly1_A), which names nothing; the
    description (theader) says what the matched structure is.

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
    """RNAcode gives coding signal independent of the gene caller, on both strands: a
    shadow ORF, the reverse complement of a real gene, should score higher antisense."""
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
        output={"tsv": str(out), "consensus": str(out.parent / "consensus.faa"),
                "pairs": str(out.parent / "pairs.tsv.gz")},
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


def test_extract_cds_recovers_origin_spanning_and_minus_strand_genes(fixture_dir):
    fasta = fixture_dir / "analysis_set.fna"
    write_fasta(fasta, [("p1", "ATGAAACCCGGGTTTTAG"), ("p2", "AAAAAA")])
    ids = fixture_dir / "dark_ids.txt"
    ids.write_text("S1\nS2\nS3\n")
    pmap = fixture_dir / "protein_map.tsv"
    pmap.write_text("S1\tp1|1\nS2\tp1|2,p1|9\nS3\tp1|3\nS4\tp2|1\n")
    index = fixture_dir / "orf_index.tsv"
    write_tsv(index, ["orf_id", "plasmid_id", "start", "end", "strand", "spans_origin"],
              [["p1|1", "p1", 1, 9, "1", "0"],
               ["p1|2", "p1", 16, 3, "1", "1"],      # 16..18 then 1..3
               ["p1|3", "p1", 4, 9, "-1", "0"],
               ["p2|1", "p2", 1, 6, "1", "0"]])
    out = fixture_dir / "dark_cds.fna"
    run_script("extract_cds.py", FakeSnakemake(
        input={"ids": str(ids), "map": str(pmap), "index": str(index), "fasta": str(fasta)},
        output=[str(out)]))
    assert out.read_text() == ">S1\nATGAAACCC\n>S2\nTAGATG\n>S3\nGGGTTT\n"


@requires("mafft")
def test_family_evolution_writes_a_consensus_per_family(fixture_dir):
    """S7b writes one consensus per family from its protein alignment, for the S7c
    re-check."""
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
        output={"tsv": str(out), "consensus": str(cons),
                "pairs": str(out.parent / "pairs.tsv.gz")},
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
    """A family whose consensus hits Pfam is not collectively novel, however dark each
    member looked on its own. This is a label, not a filter: the family stays in the table."""
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


def test_a_consensus_that_hits_only_a_duf_or_upf_family_stays_collectively_novel(
        fixture_dir, monkeypatch):
    """The hit and its name are recorded, but a domain of unknown function names nothing.
    The stand-in hmmsearch writes a domain table with one hit per family."""
    rows = ["F1 - 60 RepA_N PF01051.1 90 1e-20 70 0 1 1 1e-21 1e-21",
            "F2 - 60 DUF1234 PF06776.1 90 1e-20 70 0 1 1 1e-21 1e-21",
            "F3 - 60 UPF0126 PF03458.1 90 1e-20 70 0 1 1 1e-21 1e-21"]
    domtbl = fixture_dir / "hits.domtbl"
    domtbl.write_text("".join(r + "\n" for r in rows))
    exe = fixture_dir / "bin" / "hmmsearch"
    exe.parent.mkdir()
    exe.write_text("#!/bin/sh\nwhile [ \"$1\" != --domtblout ]; do shift; done\n"
                   f"cp {domtbl} \"$2\"\n")
    exe.chmod(0o755)
    monkeypatch.setenv("PATH", f"{exe.parent}:{os.environ['PATH']}")
    cons = fixture_dir / "family_consensus.faa"
    write_fasta(cons, [(f, "MKV") for f in ("F1", "F2", "F3")])
    fams = fixture_dir / "dark_families.tsv"
    write_tsv(fams, ["family_id"], [["F1"], ["F2"], ["F3"]])
    out = fixture_dir / "consensus_recheck.tsv"
    run_script("consensus_recheck.py", FakeSnakemake(
        input={"consensus": str(cons), "families": str(fams)}, output=[str(out)],
        params={"db": "Pfam-A.hmm", "args": "--cut_ga", "hmmer_z": 1}))
    got = {r["family_id"]: (r["consensus_hit"], r["consensus_label"], r["collectively_novel"])
           for r in read_tsv(out)}
    assert got == {"F1": ("1", "RepA_N", "0"), "F2": ("1", "DUF1234", "1"),
                   "F3": ("1", "UPF0126", "1")}


def test_defence_systems_keeps_component_status_and_system_wholeness(fixture_dir):
    """A mandatory component of a complete system and a neutral component of a fragment
    are different evidence; MacSyFinder's hit_status and sys_wholeness say which."""
    out = fixture_dir / "defence_systems.tsv"
    phase2 = out.parent / "phase2" / "RM"
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
    """Both phases write one NOT_RUN row without an ORF, and exit 0: the models are an
    optional database, and an empty table would read as 'searched, no defence system'."""
    faa = fixture_dir / "cand.faa"
    write_fasta(faa, [("GB1", "MKV")])
    mapping = fixture_dir / "map.tsv"
    write_tsv(mapping, ["gembase_id", "orf_id", "plasmid_id"], [["GB1", "p1|1", "p1"]])
    absent = str(fixture_dir / "absent")
    phase1 = fixture_dir / "defence_components.tsv"
    phase2 = fixture_dir / "defence_systems.tsv"

    for script, inputs, out in (("defence_search.py", {"faa": str(faa)}, phase1),
                                ("defence_systems.py", {"faa": str(faa), "map": str(mapping)},
                                 phase2)):
        with pytest.raises(SystemExit) as exit_info:
            run_script(script, FakeSnakemake(
                input=inputs, output={"tsv": str(out)},
                params={"models_dir": absent, "required": False}, threads=1))
        assert exit_info.value.code == 0, f"{script}: a missing optional database failed"
        rows = read_tsv(out)
        assert [r["status"] for r in rows] == ["NOT_RUN"], script
        assert rows[0].get("seq_id", rows[0].get("orf_id")) == ""


def test_defence_gembase_writes_nothing_after_a_phase1_not_run(fixture_dir):
    components = fixture_dir / "defence_components.tsv"
    write_tsv(components, ["seq_id", "component", "model", "hit_evalue", "status"],
              [["", "", "", "", "NOT_RUN"]])
    pmap = fixture_dir / "protein_map.tsv"
    pmap.write_text("S1\tp1|1\n")
    index = fixture_dir / "orf_index.tsv"
    write_tsv(index, ["orf_id", "plasmid_id", "start", "end", "strand", "partial",
                      "spans_origin", "translation_table", "seq"],
              [["p1|1", "p1", 1, 300, 1, 0, 0, 11, "MKV"]])
    faa, gmap = fixture_dir / "cand.faa", fixture_dir / "gembase_map.tsv"
    run_script("defence_gembase.py", FakeSnakemake(
        input={"components": str(components), "map": str(pmap), "index": str(index)},
        output={"faa": str(faa), "map": str(gmap)}))
    assert faa.read_text() == ""
    assert read_tsv(gmap) == []


def test_defence_gembase_writes_every_orf_of_a_candidate_plasmid(fixture_dir):
    components = fixture_dir / "defence_components.tsv"
    write_tsv(components, ["seq_id", "component", "model", "hit_evalue", "status"],
              [["S2", "CloA", "Clover", "1e-9", "SUCCESS"]])
    pmap = fixture_dir / "protein_map.tsv"
    pmap.write_text("S1\tp_1|1\nS2\tp_1|2\nS3\tp2|1\n")
    index = fixture_dir / "orf_index.tsv"
    write_tsv(index, ["orf_id", "plasmid_id", "start", "end", "strand", "partial",
                      "spans_origin", "translation_table", "seq"],
              [["p_1|1", "p_1", 900, 1200, 1, 0, 0, 11, "MKA"],
               ["p_1|2", "p_1", 100, 400, 1, 0, 0, 11, "MKB"],
               ["p2|1", "p2", 1, 300, 1, 0, 0, 11, "MKC"]])
    faa, gmap = fixture_dir / "cand.faa", fixture_dir / "gembase_map.tsv"
    run_script("defence_gembase.py", FakeSnakemake(
        input={"components": str(components), "map": str(pmap), "index": str(index)},
        output={"faa": str(faa), "map": str(gmap)}))
    assert faa.read_text() == ">p-1_00001\nMKB\n>p-1_00002\nMKA\n"
    assert [(r["gembase_id"], r["orf_id"]) for r in read_tsv(gmap)] == [
        ("p-1_00001", "p_1|2"), ("p-1_00002", "p_1|1")]


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


# Stand-in for foldseek: records its arguments and query beside itself, since the query and
# output lie in the stage's scratch directory, and writes one hit per query.
FAKE_FOLDSEEK = """#!{python}
import pathlib, shutil, sys
args = sys.argv[1:]
query, out = args[1], args[3]
here = pathlib.Path(sys.argv[0]).parent
here.joinpath("foldseek_argv.txt").write_text("\\n".join(args))
shutil.copy(query, here / "query.faa")
names = [l[1:].split()[0] for l in open(query) if l.startswith(">")]
with open(out, "w") as fh:
    for n in names:
        fh.write(f"{n}\\t1abc_A\\tA SYNTHETASE\\t0.4\\t100\\t1e-10\\t80\\n")
"""


def _structure(fixture_dir, monkeypatch, required=True, with_refs=True):
    faa = fixture_dir / "dark_proteins.faa"
    write_fasta(faa, [("rep_a", "MKTAYIAKQRQISFVKSHFSRQ"),
                      ("member_a", "MKTAYIAKQRQISFVKSHFSRK"),
                      ("rep_b", "MQQTTLNRSDEIVWCAPGHKGG")])
    families = fixture_dir / "dark_families.tsv"
    write_tsv(families, ["family_id", "representative", "members"],
              [["broad:rep_a", "rep_a", "rep_a,member_a"],
               ["broad:rep_b", "rep_b", "rep_b"]])
    db, model = fixture_dir / "pdb", fixture_dir / "prostt5"
    if with_refs:
        db.write_text("")
        model.mkdir()
    exe = fixture_dir / "bin" / "foldseek"
    exe.parent.mkdir()
    exe.write_text(FAKE_FOLDSEEK.replace("{python}", sys.executable))
    exe.chmod(0o755)
    monkeypatch.setenv("PATH", f"{exe.parent}:{os.environ['PATH']}")
    out = fixture_dir / "structure_hits.tsv"
    run_script("structure_search.py", FakeSnakemake(
        input={"faa": str(faa), "families": str(families)}, output=[str(out)],
        params={"structure": {"max_evalue": 1.0e-3, "scope": "representatives",
                              "required": required},
                "target_db": str(db), "prostt5": str(model)},
        threads=1))
    return out


def test_structure_search_restricts_to_family_representatives(fixture_dir, monkeypatch):
    """The default scope searches one sequence per dark family, not every dark protein."""
    out = _structure(fixture_dir, monkeypatch)

    query = (fixture_dir / "bin" / "query.faa").read_text()
    names = {l[1:].split()[0] for l in query.splitlines() if l.startswith(">")}
    assert names == {"rep_a", "rep_b"}
    argv = (fixture_dir / "bin" / "foldseek_argv.txt").read_text().split("\n")
    assert argv[0] == "easy-search"
    assert argv[argv.index("--prostt5-model") + 1] == str(fixture_dir / "prostt5")
    rows = read_tsv(out)
    assert {r["seq_id"] for r in rows} == {"rep_a", "rep_b"}
    assert {(r["target_description"], r["status"]) for r in rows} == {
        ("A SYNTHETASE", "SUCCESS")}
    # The query and Foldseek's raw table are intermediates, removed with the scratch
    # directory on success.
    assert sorted(p.name for p in fixture_dir.iterdir() if p.is_file()) == [
        "dark_families.tsv", "dark_proteins.faa", "pdb", "structure_hits.tsv"]


def test_structure_search_records_not_run_when_optional_and_absent(fixture_dir, monkeypatch):
    """With structure.required false and the references absent, the stage exits 0 with one
    NOT_RUN row, so the run continues and the table does not read as 'no match'."""
    with pytest.raises(SystemExit) as done:
        _structure(fixture_dir, monkeypatch, required=False, with_refs=False)
    assert done.value.code == 0
    rows = read_tsv(fixture_dir / "structure_hits.tsv")
    assert [(r["seq_id"], r["status"]) for r in rows] == [("", "NOT_RUN")]
    assert not (fixture_dir / "bin" / "foldseek_argv.txt").exists()


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
    """An IS element is an island, like an integron."""
    rows = _run_context(
        fixture_dir,
        is_rows=[["pl1", "pl1|IS1", "IS3", "IS3_1", 50, 450, "+", 1, "1e-50", ""]])
    assert rows[0]["cons_is_element"] == "1.0"


# Stand-in for integron_finder: keeps its topology file beside itself and writes one In0
# integrase per replicon in the IntegronFinder 2.0.6 column order.
FAKE_INTEGRON_FINDER = """#!{python}
import pathlib, shutil, sys
args = sys.argv[1:]
here = pathlib.Path(sys.argv[0]).parent
shutil.copy(args[args.index("--topology-file") + 1], here / "topology.txt")
out = pathlib.Path(args[args.index("--outdir") + 1])
out.mkdir(parents=True)
with open(out / "chunk.integrons", "w") as fh:
    for line in open(args[-1]):
        if line.startswith(">"):
            rep = line[1:].split()[0]
            fh.write("\\t".join(["integron_01", rep, "intI_1", "10", "900", "1", "1e-50",
                                 "protein", "intI", "intersection_tyr_intI", "In0", "Yes",
                                 "NA", "circ"]) + "\\n")
"""


def test_integrons_passes_each_plasmid_its_registry_topology(fixture_dir, monkeypatch):
    exe = fixture_dir / "bin" / "integron_finder"
    exe.parent.mkdir()
    exe.write_text(FAKE_INTEGRON_FINDER.replace("{python}", sys.executable))
    exe.chmod(0o755)
    monkeypatch.setenv("PATH", f"{exe.parent}:{os.environ['PATH']}")
    fasta = fixture_dir / "analysis_set.fna"
    write_fasta(fasta, [("c1", "ACGT" * 10), ("l1", "ACGT" * 10), ("d1", "ACGT" * 10)])
    master = fixture_dir / "master.tsv"
    write_tsv(master, ["plasmid_id", "topology"],
              [["c1", "circular"], ["l1", "linear"], ["d1", "direct terminal repeat"]])
    out = fixture_dir / "out" / "integrons.tsv"
    out.parent.mkdir()
    run_script("integrons.py", FakeSnakemake(
        input={"fasta": str(fasta), "master": str(master)}, output=[str(out)], threads=1))
    # The format integron_finder.topology parses: "<replicon> <circ|lin>".
    assert (exe.parent / "topology.txt").read_text().splitlines() == [
        "c1 circ", "l1 lin", "d1 circ"]
    assert {r["plasmid_id"]: r["integron_type"] for r in read_tsv(out)} == {
        "c1": "In0", "l1": "In0", "d1": "In0"}
