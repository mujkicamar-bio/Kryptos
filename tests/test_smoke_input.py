"""Smoke tests: analysis set, gene calling, and dereplication.

Each test runs one workflow script against a small fixture.
"""
import csv
import pathlib

import pytest
from conftest import (
    FakeSnakemake,
    read_tsv,
    requires,
    run_script,
    write_fasta,
    write_tsv,
)

from plasmidann.fasta import iter_fasta


TAXDUMP = str(pathlib.Path(__file__).parent / "data" / "taxdump")


def _host_sources(fixture_dir, ps_hosts=(), organisms=()):
    """analysis_set's host inputs: PlasmidScope hosts [(plasmid_id, host)] and working-set
    rows [(plasmid_id, lifestyle, organism)]."""
    ps = fixture_dir / "ps_hosts.tsv"
    write_tsv(ps, ["plasmid_id", "host"], list(ps_hosts))
    ws = fixture_dir / "working_set.tsv"
    write_tsv(ws, ["plasmid_id", "lifestyle", "organism"], list(organisms))
    return {"ps_hosts": str(ps), "working_set": str(ws)}


def test_analysis_set_excludes_plasmids_with_a_eukaryotic_host(fixture_dir):
    """A yeast host is excluded and counted; a bacterial host, a genus name that is also an
    animal genus, and the sampled human of a metagenomic record are not."""
    ids = ["yeast", "coli", "bacillus", "human_gut"]
    master = fixture_dir / "master.tsv"
    write_tsv(master, ["plasmid_id", "hab_top", "size_bp", "plsdb_species"],
              [[p, "Unknown", 4, "Escherichia_coli" if p == "coli" else ""] for p in ids])
    fasta = fixture_dir / "in.fna"
    write_fasta(fasta, [(p, "ATGC") for p in ids])
    out = {k: str(fixture_dir / f"{k}.out")
           for k in ("ids", "fasta", "small_ids", "repeats", "lengths")}
    log = fixture_dir / "analysis_set.log"

    run_script("analysis_set.py", FakeSnakemake(
        input={"master": str(master), "fasta": str(fasta), **_host_sources(
            fixture_dir, ps_hosts=[("yeast", "Saccharomyces cerevisiae S288C"),
                                   ("bacillus", "Bacillus sp. X1")],
            organisms=[("human_gut", "metagenomic", "Homo sapiens")])},
        output=out, log=[str(log)],
        params={"exclude": [], "max_size_bp": 20000, "min_terminal_repeat_bp": 20,
                "taxdump": TAXDUMP}))

    assert open(out["ids"]).read().split() == ["coli", "bacillus", "human_gut"]
    assert "1 plasmids with a eukaryotic host excluded: Saccharomyces cerevisiae 1" \
        in log.read_text()


def test_a_scripts_output_reaches_its_declared_log(fixture_dir):
    """Snakemake does not redirect a script's stdout to the rule's `log:`; _ctx does, so
    each rule's diagnostics can be attributed to it."""
    master = fixture_dir / "master.tsv"
    write_tsv(master, ["plasmid_id", "hab_top", "size_bp"],
              [["p1", "Unknown", 4], ["p2", "Unknown", 4]])
    fasta = fixture_dir / "in.fna"
    write_fasta(fasta, [("p1", "ATGC"), ("p2", "GGCC")])
    log = fixture_dir / "logs" / "analysis_set.log"

    run_script("analysis_set.py", FakeSnakemake(
        input={"master": str(master), "fasta": str(fasta), **_host_sources(fixture_dir)},
        output={"ids": str(fixture_dir / "ids.txt"), "fasta": str(fixture_dir / "out.fna"),
                "small_ids": str(fixture_dir / "small.txt"),
                "repeats": str(fixture_dir / "repeats.tsv"),
                "lengths": str(fixture_dir / "lengths.tsv")},
        params={"exclude": [], "max_size_bp": 20000, "min_terminal_repeat_bp": 20,
                "taxdump": TAXDUMP},
        log=[str(log)]))

    assert log.exists(), "the declared log file was never created"
    assert "2 plasmids" in log.read_text(), (
        f"the script's own diagnostics did not reach its log: {log.read_text()!r}")


def test_analysis_set_keeps_every_plasmid_and_lists_the_small_ones(fixture_dir):
    """Every plasmid that passes the exclusion is in the analysis set; those below the
    size cut-off are also listed as small. 20,000 bp itself is not small."""
    master = fixture_dir / "master.tsv"
    write_tsv(master, ["plasmid_id", "hab_top", "size_bp"],
              [["small", "Unknown", 19999], ["edge", "Unknown", 20000],
               ["large", "Unknown", 90000], ["sim", "Simulated-artifact", 5000]])
    fasta = fixture_dir / "in.fna"
    write_fasta(fasta, [(p, "ATGC") for p in ("small", "edge", "large", "sim")])
    out = {"ids": str(fixture_dir / "ids.txt"), "fasta": str(fixture_dir / "out.fna"),
           "small_ids": str(fixture_dir / "small.txt"),
           "repeats": str(fixture_dir / "repeats.tsv"),
           "lengths": str(fixture_dir / "lengths.tsv")}

    run_script("analysis_set.py", FakeSnakemake(
        input={"master": str(master), "fasta": str(fasta), **_host_sources(fixture_dir)}, output=out,
        params={"exclude": ["Simulated-artifact"], "max_size_bp": 20000,
                "min_terminal_repeat_bp": 20,
                "taxdump": TAXDUMP}))

    assert open(out["ids"]).read().split() == ["small", "edge", "large"]
    assert [l[1:].strip() for l in open(out["fasta"]) if l.startswith(">")] == [
        "small", "edge", "large"]
    assert open(out["small_ids"]).read().split() == ["small"]


def test_analysis_set_lists_only_plasmids_with_a_sequence(fixture_dir):
    """An in-scope plasmid with no FASTA record would enter the small-plasmid denominators
    without any genes; a record given twice would contribute its genes twice."""
    master = fixture_dir / "master.tsv"
    write_tsv(master, ["plasmid_id", "hab_top", "size_bp"],
              [["p1", "Unknown", 4], ["p2", "Unknown", 4]])
    fasta = fixture_dir / "in.fna"
    out = {k: str(fixture_dir / f"{k}.out")
           for k in ("ids", "fasta", "small_ids", "repeats", "lengths")}
    snake = FakeSnakemake(
        input={"master": str(master), "fasta": str(fasta), **_host_sources(fixture_dir)}, output=out,
        params={"exclude": [], "max_size_bp": 20000, "min_terminal_repeat_bp": 20,
                "taxdump": TAXDUMP})

    write_fasta(fasta, [("p1", "ATGC")])
    run_script("analysis_set.py", snake)
    assert open(out["ids"]).read().split() == ["p1"]
    assert open(out["small_ids"]).read().split() == ["p1"]

    write_fasta(fasta, [("p1", "ATGC"), ("p1", "ATGC")])
    with pytest.raises(SystemExit, match="twice"):
        run_script("analysis_set.py", snake)


def test_a_circular_records_terminal_repeat_is_written_once(fixture_dir):
    """A 'direct terminal repeat' record starts with a copy of its own last bases. S0
    writes the molecule with the last copy removed; a linear record with the same ends is
    a genuinely linear molecule and is written as it is."""
    import random
    random.seed(1)
    core = "".join(random.choice("ACGT") for _ in range(600))
    repeat = "GATTACAGATTACAGATTACAGA"        # 23 bp: not a multiple of 3
    record = repeat + core + repeat
    master = fixture_dir / "master.tsv"
    write_tsv(master, ["plasmid_id", "hab_top", "size_bp", "topology"],
              [["dtr", "Unknown", len(record), "direct terminal repeat"],
               ["lin", "Unknown", len(record), "linear"]])
    fasta = fixture_dir / "in.fna"
    write_fasta(fasta, [("dtr", record), ("lin", record)])
    out = {"ids": str(fixture_dir / "ids.txt"), "fasta": str(fixture_dir / "out.fna"),
           "small_ids": str(fixture_dir / "small.txt"),
           "repeats": str(fixture_dir / "repeats.tsv"),
           "lengths": str(fixture_dir / "lengths.tsv")}

    run_script("analysis_set.py", FakeSnakemake(
        input={"master": str(master), "fasta": str(fasta), **_host_sources(fixture_dir)}, output=out,
        params={"exclude": [], "max_size_bp": 20000, "min_terminal_repeat_bp": 20,
                "taxdump": TAXDUMP}))

    written = dict(iter_fasta([out["fasta"]]))
    assert written["dtr"] == repeat + core
    assert written["lin"] == record
    rows = list(csv.DictReader(open(out["repeats"]), delimiter="\t"))
    assert [(r["plasmid_id"], r["repeat_bp"]) for r in rows] == [("dtr", "23")]
    # Gene coordinates refer to the written record, so its length is the one recorded.
    lengths = {r["plasmid_id"]: int(r["length_bp"])
               for r in csv.DictReader(open(out["lengths"]), delimiter="\t")}
    assert lengths == {"dtr": len(repeat + core), "lin": len(record)}


ANTIFAM = "data/refs/antifam/AntiFam.hmm"


# A sequence derived from the Spurious_ORF_67 consensus and degraded until it sits between
# the two thresholds: domain score 32.4 bits, which clears that family's curated GA of
# 19.2, and i-Evalue 1.5e-05 at Z=3,497,616, which does not clear a 1e-5 floor.
SHADOW_LIKE = "AIFSANLAARPTAAAAAMRRAEKAFSLAQGAFAAGAPAAFAAAAAYLLGK"


NOT_AN_ARTEFACT = ("MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVK"
                   "ALPDAQFEVVHSLAKWKR")


@pytest.mark.skipif(not pathlib.Path(ANTIFAM).exists(), reason="AntiFam not downloaded")
@requires("hmmsearch", "tantan")
def test_the_artefact_screen_uses_antifams_curated_thresholds(fixture_dir):
    """274 of AntiFam's 278 profiles carry a curated gathering threshold LOOSER than
    E=1e-5 at Z=3,497,616 - the median curated cut corresponds to E=7.3e-4, about 73x
    looser. A blanket E-value floor therefore overrides the curator on 98.6% of the
    artefact database, in the one screen whose job is to stop shadow ORFs entering the
    dark set with a perfect score.

    This is the same argument that already sets --cut_ga on T1, applied to the database
    where the consequence is worse: a missed artefact is not a missed annotation, it is a
    non-protein sent to the bench."""
    faa = fixture_dir / "unique_proteins.faa"
    write_fasta(faa, [("shadow_like", SHADOW_LIKE), ("real", NOT_AN_ARTEFACT)])
    out = fixture_dir / "artefact_flags.tsv"

    run_script("artefact_screen.py", FakeSnakemake(
        input={"faa": str(faa), "preflight": ""},
        output=[str(out)],
        params={"artefact": {"antifam_db": ANTIFAM, "antifam_args": "--cut_ga",
                             "antifam_max_evalue": None,
                             "max_low_complexity_fraction": 0.5},
                "hmmer_z": 3497616},
        threads=2))

    rows = {r["seq_id"]: r for r in read_tsv(out)}
    assert rows["shadow_like"]["artefact_flag"] == "1", (
        "a sequence above AntiFam's own gathering threshold was not flagged - the screen "
        "is still applying a blanket E-value floor over the curation")
    assert rows["shadow_like"]["antifam_family"] == "Spurious_ORF_67"
    assert rows["real"]["artefact_flag"] == "0", "a real protein must not be flagged"


def test_orf_call_reconstructs_a_gene_across_the_origin(fixture_dir):
    """orf_call drives the origin repair: a plasmid cut inside a gene yields the same
    proteins as the uncut plasmid, with the gene written start > end and spans_origin=1.
    The second pair carries a gene longer than half the molecule."""
    import random
    pytest.importorskip("pyrodigal")

    def planted(seed, n_codons, left, right, cut):
        rng = random.Random(seed)
        sense = [a + b + c for a in "ACGT" for b in "ACGT" for c in "ACGT"
                 if a + b + c not in ("TAA", "TAG", "TGA")]
        gene = "ATG" + "".join(rng.choice(sense) for _ in range(n_codons)) + "TAA"
        flank = "".join(rng.choice("ACGT") for _ in range(left + right))
        uncut = flank[:left] + gene + flank[left:]
        cut_at = left + cut                         # the new origin lies inside the gene
        return uncut, uncut[cut_at:] + uncut[:cut_at]

    short_uncut, short_cut = planted(11, 300, 1_000, 1_000, 453)
    long_uncut, long_cut = planted(5, 400, 700, 0, 100)
    records = [("uncut", short_uncut), ("cut", short_cut),
               ("long_uncut", long_uncut), ("long_cut", long_cut)]
    master = fixture_dir / "master.tsv"
    write_tsv(master, ["plasmid_id", "topology"], [[pid, "circular"] for pid, _ in records])
    fasta = fixture_dir / "analysis_set.fna"
    write_fasta(fasta, records)
    out = fixture_dir / "orfs.tsv"

    run_script("orf_call.py", FakeSnakemake(
        input={"fasta": str(fasta), "master": str(master)},
        output={"tsv": str(out)},
        params={"min_orf_aa": 20}))

    rows = read_tsv(out)
    proteins = {pid: sorted(r["seq"] for r in rows if r["plasmid_id"] == pid)
                for pid, _ in records}
    assert proteins["uncut"] == proteins["cut"], "protein set depends on the cut point"
    assert proteins["long_uncut"] == proteins["long_cut"], "long gene depends on the cut"

    for pid, n_aa in (("cut", 301), ("long_cut", 401)):
        gene = [r for r in rows if r["plasmid_id"] == pid and len(r["seq"]) == n_aa]
        assert len(gene) == 1, f"the planted gene was not called exactly once on {pid}"
        assert gene[0]["spans_origin"] == "1"
        assert int(gene[0]["start"]) > int(gene[0]["end"]), "wrapped end not applied"
        assert gene[0]["partial"] == "0"
    assert not any(r["partial"] == "1" for r in rows), "a left-edge stub survived"
    assert {r["translation_table"] for r in rows} <= {"11", "4"}, rows[0]


def test_plasmidscope_import_keeps_only_our_proteins(fixture_dir):
    """The whole ALL table is read; only proteins identical to one of ours are written,
    keyed by our seq_id, so every later stage can join on it."""
    from plasmidann.dereplicate import _seq_id
    ours, theirs = "MKVLATTLLG", "MQQQQQQQQQ"
    faa = fixture_dir / "unique.faa"
    write_fasta(faa, [(_seq_id(ours), ours)])
    header = ["Plasmid_ID", "Protein_ID", "Orf Prediction Source", "Product",
              "COG_category", "COG_id", "KEGG_ko", "KEGG_Pathway", "PFAMs", "GOs",
              "EC_number", "Sequence"]
    ps = fixture_dir / "ALL.protein_list.tsv"
    write_tsv(ps, header,
              [["p1", "X_1", "Prodigal:2.6", "hypothetical protein", "S", "-", "-", "-",
                "RHH_1", "-", "-", ours + "*"],
               ["p2", "X_2", "Prodigal:2.6", "hypothetical protein", "S", "-", "-", "-",
                "-", "-", "-", theirs + "*"]])
    out = fixture_dir / "plasmidscope_proteins.tsv"

    run_script("plasmidscope_import.py", FakeSnakemake(
        input={"faa": str(faa), "ps": str(ps)}, output=[str(out)]))

    rows = read_tsv(out)
    assert [r["seq_id"] for r in rows] == [_seq_id(ours)]
    assert rows[0]["ps_class"] == "ANNOTATED" and rows[0]["pfams"] == "RHH_1"


def test_make_test_set_refuses_a_master_without_hab_top(fixture_dir):
    """Without the column the locked exclusion cannot be applied, so the tool must fail
    rather than sample simulated records."""
    import subprocess
    import sys
    master = fixture_dir / "master.tsv"
    write_tsv(master, ["plasmid_id", "topology", "size_bp"], [["p1", "circular", 4000]])
    tool = pathlib.Path(__file__).resolve().parents[1] / "tools" / "make_test_set.py"
    r = subprocess.run([sys.executable, str(tool), "--master", str(master),
                        "--fasta", "unused.fna.gz", "--out", str(fixture_dir / "out.fna")],
                       capture_output=True, text=True)
    assert r.returncode != 0 and "hab_top" in r.stderr
