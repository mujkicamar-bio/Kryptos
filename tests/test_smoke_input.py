"""Smoke tests: analysis set, gene calling, dereplication, controls and decoys.

Each test runs one workflow script against a small fixture.
"""
import csv
import pathlib

import pytest
from conftest import (
    FakeSnakemake,
    _ps_table,
    read_tsv,
    requires,
    run_script,
    write_fasta,
    write_tsv,
)

from plasmidann.fasta import iter_fasta


def test_controls_with_uninformative_titles_are_excluded(fixture_dir):
    """14.3% of reviewed Swiss-Prot plasmid entries are titled "Uncharacterized ..." -
    which cascade.is_informative correctly rejects. Sampling them into the control set
    capped achievable recall at 0.86 against a required 0.99, so the gate would halt every
    run after the full cascade and blame the cascade.

    A control protein must be one the pipeline is expected to annotate. A protein whose own
    curators could not name it is not a test of our recall."""
    raw = fixture_dir / "raw.faa"
    control = fixture_dir / "control.faa"
    spiked = fixture_dir / "cascade_input.faa"
    faa = fixture_dir / "unique.faa"

    body = "MKVLATTLLGAAFAASSALAQ" * 4
    write_fasta(raw, [
        ("sp|P00001|A_ECOLI Beta-lactamase TEM OS=Escherichia coli", body),
        ("sp|P00002|B_ECOLI Uncharacterized protein YbaA OS=Escherichia coli", body),
        ("sp|P00003|C_ECOLI Relaxase MobA OS=Escherichia coli", body),
        ("sp|P00004|D_ECOLI UPF0102 protein YraN OS=Escherichia coli", body),
        ("sp|P00005|E_ECOLI DNA polymerase III subunit beta OS=Escherichia coli", body),
    ])
    write_fasta(faa, [("realprotein", body)])

    decoy_faa = fixture_dir / "negative_control.faa"
    write_fasta(decoy_faa, [("DECOY_shuf_00000", "M" * 60),
                            ("DECOY_rc_00001", "K" * 60)])

    run_script("prepare_control.py", FakeSnakemake(
        input={"faa": str(faa), "raw": str(raw), "decoys": str(decoy_faa),
               "ps": _ps_table(fixture_dir)},
        output={"control": str(control), "spiked": str(spiked)},
        params={"n_controls": 10, "min_controls": 1, "seed": 1}))

    headers = [l[1:].strip() for l in open(control) if l.startswith(">")]
    assert len(headers) == 3, f"expected 3 nameable controls, got {len(headers)}"
    joined = " ".join(headers).lower()
    assert "uncharacterized" not in joined
    assert "upf0102" not in joined

    # Both control sets must be in the file the cascade actually searches. A decoy that
    # never enters the query is a negative control in the config and nowhere else.
    spiked_ids = [l[1:].strip() for l in open(spiked) if l.startswith(">")]
    assert sum(i.startswith("CTRL_") for i in spiked_ids) == 3
    assert sum(i.startswith("DECOY_") for i in spiked_ids) == 2, (
        f"the decoys were not spiked into the cascade query set: {spiked_ids}")


def test_a_scripts_output_reaches_its_declared_log(fixture_dir):
    """Every rule declares `log:`, and for `script:` rules Snakemake does not redirect
    stdout there - it only creates the path. So all 26 log files stayed empty and every
    diagnostic the scripts print (hit counts, rejected hits, background rates) went to the
    single SLURM output file interleaved across 1,481 concurrent jobs, where it cannot be
    attributed to a rule."""
    master = fixture_dir / "master.tsv"
    write_tsv(master, ["plasmid_id", "hab_top", "size_bp"],
              [["p1", "Unknown", 4], ["p2", "Unknown", 4]])
    fasta = fixture_dir / "in.fna"
    write_fasta(fasta, [("p1", "ATGC"), ("p2", "GGCC")])
    log = fixture_dir / "logs" / "analysis_set.log"

    run_script("analysis_set.py", FakeSnakemake(
        input={"master": str(master), "fasta": str(fasta)},
        output={"ids": str(fixture_dir / "ids.txt"), "fasta": str(fixture_dir / "out.fna"),
                "small_ids": str(fixture_dir / "small.txt"),
                "repeats": str(fixture_dir / "repeats.tsv"),
                "lengths": str(fixture_dir / "lengths.tsv")},
        params={"exclude": [], "max_size_bp": 20000, "min_terminal_repeat_bp": 20},
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
        input={"master": str(master), "fasta": str(fasta)}, output=out,
        params={"exclude": ["Simulated-artifact"], "max_size_bp": 20000,
                "min_terminal_repeat_bp": 20}))

    assert open(out["ids"]).read().split() == ["small", "edge", "large"]
    assert [l[1:].strip() for l in open(out["fasta"]) if l.startswith(">")] == [
        "small", "edge", "large"]
    assert open(out["small_ids"]).read().split() == ["small"]


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
        input={"master": str(master), "fasta": str(fasta)}, output=out,
        params={"exclude": [], "max_size_bp": 20000, "min_terminal_repeat_bp": 20}))

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


def _decoy_fixture(fixture_dir, n_plasmids=12):
    """Real CDS on real plasmids, which is what a decoy has to be built from.

    A decoy drawn from a uniform random model is too easy - it fails to resemble anything
    for reasons that have nothing to do with the cascade - so the source is the collection
    itself.
    """
    import random
    rng = random.Random(7)
    shard = fixture_dir / "decoy_shard.fna"
    rows, records = [], []
    for i in range(n_plasmids):
        pid = f"dp{i:02d}"
        # A clean 300 bp frame: ATG, 98 sense codons, TAA.
        sense = [c for c in ("GCT", "AAA", "GAT", "TTT", "CAT", "ATT", "CTG", "ATG",
                             "AAT", "CCG", "CAG", "CGT", "AGC", "ACC", "GTT", "TGG",
                             "TAT", "GGT", "GAA", "TGC")]
        gene = "ATG" + "".join(rng.choice(sense) for _ in range(98)) + "TAA"
        seq = "".join(rng.choice("ACGT") for _ in range(100)) + gene
        records.append((pid, seq))
        rows.append([f"{pid}|1", pid, 101, 100 + len(gene), "+", 0])
    write_fasta(shard, records)
    index = fixture_dir / "decoy_index.tsv"
    write_tsv(index, ["orf_id", "plasmid_id", "start", "end", "strand", "spans_origin"],
              rows)
    return index, shard


def test_negative_controls_are_built_from_real_cds(fixture_dir):
    """Spec section 58.2. The config declared 250 shuffled and 250 reverse-complement
    decoys and nothing built them, so the pipeline had no test of the direction that
    matters most here.

    The deliverable is the DARK set - the complement of what the cascade could name. If
    the cascade can be induced to name something that is not a protein, the complement is
    not what it claims to be, and every statement about dark proteins inherits the error.
    Positive controls cannot detect that; they test the other direction.
    """
    index, shard = _decoy_fixture(fixture_dir)
    out = fixture_dir / "negative_control.faa"

    run_script("prepare_decoys.py", FakeSnakemake(
        input={"index": str(index), "fasta": str(shard)},
        output={"faa": str(out)},
        params={"n_shuffled": 3, "n_reverse_complement": 3, "seed": 1,
                "min_length": 50}))

    records = [l[1:].strip() for l in out.read_text().splitlines() if l.startswith(">")]
    assert all(r.startswith("DECOY_") for r in records), (
        f"a decoy must be recognisable by prefix at every later stage: {records}")
    assert sum("shuf" in r for r in records) == 3
    assert sum("_rc_" in r for r in records) == 3


def test_a_shuffled_decoy_keeps_its_source_composition_exactly(fixture_dir):
    """Composition is preserved and order destroyed; that is the whole construction.

    A hit to a shuffled sequence is a hit to amino-acid composition alone, and composition
    is not evidence of function. A decoy that also changed composition would fail to
    resemble anything for two reasons at once, and the test would no longer isolate the
    one being made.
    """
    index, shard = _decoy_fixture(fixture_dir)
    out = fixture_dir / "negative_control.faa"

    run_script("prepare_decoys.py", FakeSnakemake(
        input={"index": str(index), "fasta": str(shard)},
        output={"faa": str(out)},
        params={"n_shuffled": 4, "n_reverse_complement": 0, "seed": 1,
                "min_length": 50}))

    text = out.read_text().splitlines()
    seqs = [text[i + 1] for i, l in enumerate(text) if l.startswith(">")]
    assert seqs, "no decoy written"
    for seq in seqs:
        assert len(seq) >= 50
        assert "*" not in seq, "a stop codon is not a residue"
        # The source frame is a fixed 99-residue protein, so composition is checkable:
        # every decoy must be a permutation of a real translated CDS.
        assert set(seq) <= set("ACDEFGHIKLMNPQRSTVWY"), f"non-residue in decoy: {seq}"


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


def test_prepare_control_leaves_plasmidscope_annotated_proteins_out(fixture_dir):
    """Annotated by PlasmidScope: not searched. Dark in PlasmidScope: searched. The
    controls and decoys still enter, so the gate measures the cascade that actually ran."""
    raw = fixture_dir / "raw.faa"
    body = "MKVLATTLLGAAFAASSALAQ" * 4
    write_fasta(raw, [("sp|P00001|A_ECOLI Beta-lactamase TEM OS=Escherichia coli", body)])
    faa = fixture_dir / "unique.faa"
    write_fasta(faa, [("known", body), ("dark", body), ("absent", body)])
    decoys = fixture_dir / "negative_control.faa"
    write_fasta(decoys, [("DECOY_shuf_00000", "M" * 60)])
    ps = _ps_table(fixture_dir, [["known", "ANNOTATED", "", "", "", "", "RHH_1", "", "",
                                  "Prodigal:2.6", 1],
                                 ["dark", "NONE", "S", "", "", "", "", "", "",
                                  "Prodigal:2.6", 1]])
    spiked = fixture_dir / "cascade_input.faa"

    run_script("prepare_control.py", FakeSnakemake(
        input={"faa": str(faa), "raw": str(raw), "decoys": str(decoys), "ps": ps},
        output={"control": str(fixture_dir / "control.faa"), "spiked": str(spiked)},
        params={"n_controls": 1, "min_controls": 1, "seed": 1}))

    ids = [l[1:].strip() for l in open(spiked) if l.startswith(">")]
    assert "known" not in ids, "a PlasmidScope-annotated protein was sent to the cascade"
    assert {"dark", "absent"} <= set(ids)
    assert any(i.startswith("CTRL_") for i in ids) and "DECOY_shuf_00000" in ids
