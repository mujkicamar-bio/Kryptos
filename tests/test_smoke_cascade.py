"""Smoke tests: pre-flight, the annotation cascade search and its resolution.

Each test runs one workflow script against a small fixture.
"""
import pytest
from conftest import (
    FakeSnakemake,
    _ps_table,
    _selection,
    read_tsv,
    requires,
    run_script,
    write_fasta,
    write_tsv,
)

from plasmidann.tools import required_tools

# Pre-flight fails unless every registry tool of a run without structure and orthology is
# installed, so a test that expects it to pass needs them all.
every_tool = requires(*(t["name"] for t in required_tools(structure_required=False,
                                                          orthology_required=False)))

# The S2b artefact flags as artefact_screen writes them. Empty unless rows are given.
ARTEFACT_COLS = ["seq_id", "artefact_flag", "antifam_family", "antifam_ievalue",
                 "low_complexity_fraction", "artefact_reason"]


def _artefact_flags(fixture_dir, rows=()):
    path = fixture_dir / "artefact_flags.tsv"
    write_tsv(path, ARTEFACT_COLS, list(rows))
    return str(path)


def _s4d_params(fixture_dir, labels_required=False, amr_required=False,
                conj_required=False, conj_exe=None, conj_models=None, labels_dir=None):
    """Pre-flight params for the label databases, AMRFinderPlus and CONJScan. By default
    none is installed and none is required, so a test not about them sees them skipped."""
    return {"labels": {"dir": str(labels_dir or fixture_dir / "no_labels"),
                       "required": labels_required},
            "amrfinder": {"executable": str(fixture_dir / "no_amr" / "amrfinder"),
                          "database": str(fixture_dir / "no_amr_db"),
                          "required": amr_required},
            "conjugation": {"required": conj_required, "version": "2.1.0",
                            "exe": str(conj_exe or fixture_dir / "no_msf" / "macsyfinder")},
            "conjscan_models": str(conj_models or fixture_dir / "no_conjscan")}


@every_tool
def test_preflight_refuses_an_incomplete_diamond_index(fixture_dir):
    """A makedb killed part-way leaves a non-empty .dmnd. Existence is not completeness:
    the sequence count must match the one the database's metadata declares."""
    db = _tiny_diamond_db(fixture_dir)
    antifam = fixture_dir / "AntiFam.hmm"
    antifam.write_text("HMMER3/f\n")
    (fixture_dir / "AntiFam.hmm.h3i").write_text("")
    params = {"artefact": {"antifam_db": str(antifam)}, "structure": {"required": False},
              "orthology": {"required": False, "data_dir": ""},
              "foldseek_db": str(fixture_dir / "absent"),
              "prostt5": str(fixture_dir / "absent"), **_s4d_params(fixture_dir)}
    tier = {"id": "T5", "method": "diamond", "source": "nr", "db": str(db),
            "args": "", "max_evalue": 1e-5}

    with pytest.raises(SystemExit, match="holds 3 sequences, expected 4"):
        run_script("preflight.py", FakeSnakemake(
            output=[str(fixture_dir / "p1.tsv")],
            params={**params, "tiers": [{**tier, "expected_sequences": 4}]}))
    run_script("preflight.py", FakeSnakemake(
        output=[str(fixture_dir / "p2.tsv")],
        params={**params, "tiers": [{**tier, "expected_sequences": 3}]}))


def test_preflight_fails_when_a_downstream_tool_is_missing(fixture_dir, monkeypatch):
    """A missing tool of a later stage fails pre-flight, naming the tool and its stage.

    The harness prepends the project's bin to PATH, so the absence is injected at
    shutil.which."""
    import shutil
    db = fixture_dir / "fake.hmm"
    db.write_text("HMMER3/f\n")
    (fixture_dir / "fake.hmm.h3i").write_text("")
    antifam = fixture_dir / "AntiFam.hmm"
    antifam.write_text("HMMER3/f\n")
    (fixture_dir / "AntiFam.hmm.h3i").write_text("")

    real_which = shutil.which
    monkeypatch.setattr(shutil, "which",
                        lambda name, *a, **k: None if name == "mafft" else real_which(name))

    with pytest.raises(SystemExit) as exc:
        run_script("preflight.py", FakeSnakemake(
            output=[str(fixture_dir / "preflight.tsv")],
            params={"tiers": [{"id": "T1", "method": "hmmer", "source": "pfam", "db": str(db),
                               "args": "--cut_ga", "max_evalue": None}],
                    "artefact": {"antifam_db": str(antifam)},
                    "structure": {"required": False},
                    "foldseek_db": "", "prostt5": "", **_s4d_params(fixture_dir)}))
    assert "mafft" in str(exc.value), "a missing mafft must be named by pre-flight"
    assert "S7b" in str(exc.value), "pre-flight must say which stage the tool belongs to"


def _preflight_s4d(fixture_dir, **s4d):
    """Run pre-flight with every cascade check satisfied and the given S4d/S8f params."""
    db = fixture_dir / "fake.hmm"
    db.write_text("HMMER3/f\n")
    (fixture_dir / "fake.hmm.h3i").write_text("")
    out = fixture_dir / "preflight.tsv"
    run_script("preflight.py", FakeSnakemake(
        output=[str(out)],
        params={"tiers": [{"id": "T1", "method": "hmmer", "source": "pfam", "db": str(db),
                           "args": "--cut_ga", "max_evalue": None}],
                "artefact": {"antifam_db": str(db)},
                "structure": {"required": False},
                "orthology": {"required": False, "data_dir": ""},
                "foldseek_db": "", "prostt5": "", **_s4d_params(fixture_dir, **s4d)}))
    return out


def _label_dbs(root, omit=()):
    """A complete data/refs/labels tree, less the databases named in `omit`."""
    from plasmidann import labeldb
    for db in labeldb.DATABASES:
        if db in omit:
            continue
        d = root / db
        d.mkdir(parents=True)
        (d / "VERSION").write_text(f"{db} 1.0\n")
        (d / ("card.json" if db == "card" else f"{db}.faa")).write_text(">x\nMKV\n")
    (root / "bacmet" / "raw").mkdir(exist_ok=True)
    (root / "bacmet" / "raw" / "BacMet2_EXP.753.mapping.txt").write_text("x\n")
    return root


@every_tool
def test_preflight_fails_fast_when_a_required_label_database_is_missing(fixture_dir):
    labels = _label_dbs(fixture_dir / "labels", omit=("acrdb",))
    with pytest.raises(SystemExit, match="label database acrdb"):
        _preflight_s4d(fixture_dir, labels_dir=labels, labels_required=True)
    # Not required: the stage records it NOT_RUN, so pre-flight lets it through.
    _preflight_s4d(fixture_dir, labels_dir=labels, labels_required=False)


@every_tool
def test_preflight_refuses_an_installed_label_database_without_its_version(fixture_dir):
    """The stage halts on a missing VERSION whether the databases are required or not,
    so pre-flight must too - in seconds, not after the cascade."""
    labels = _label_dbs(fixture_dir / "labels")
    (labels / "tadb" / "VERSION").write_text("")
    with pytest.raises(SystemExit, match="tadb/VERSION is missing or empty"):
        _preflight_s4d(fixture_dir, labels_dir=labels, labels_required=False)


@every_tool
def test_preflight_fails_when_amrfinder_is_required_and_absent(fixture_dir):
    with pytest.raises(SystemExit, match="AMRFinderPlus executable"):
        _preflight_s4d(fixture_dir, amr_required=True)
    _preflight_s4d(fixture_dir, amr_required=False)


def _conjscan_install(fixture_dir, grammar, reported):
    """A CONJScan 2.1.0 tree of the given grammar, and a MacSyFinder that reports
    `reported` from --version, as 2.1.4 ('Macsyfinder 2.1.4') and 2.1.6 do."""
    models = fixture_dir / f"conjscan_{grammar}_{reported.split()[-1]}"
    definitions = models / "CONJScan" / "definitions" / "Plasmids"
    definitions.mkdir(parents=True)
    (models / "CONJScan" / "metadata.yml").write_text("name: CONJScan\nvers: 2.1.0\n")
    (definitions / "MOB.xml").write_text(
        f'<model inter_gene_max_space="500" min_mandatory_genes_required="1" '
        f'min_genes_required="1" vers="{grammar}">\n</model>\n')
    exe = models / "bin" / "macsyfinder"
    exe.parent.mkdir()
    exe.write_text(f"#!/bin/sh\necho '{reported} '\n")
    exe.chmod(0o755)
    return models, exe


@every_tool
def test_preflight_refuses_a_macsyfinder_too_old_for_the_conjscan_grammar(fixture_dir):
    """CONJScan 2.1.0 is written in grammar 2.1; MacSyFinder 2.1.4, which DefenseFinder
    pins, stops on it with a parse error. That must fail here, naming both versions, and
    also when conjugation is not required - the stage runs whenever both are installed."""
    models, exe = _conjscan_install(fixture_dir, "2.1", "Macsyfinder 2.1.4")
    with pytest.raises(SystemExit) as exc:
        _preflight_s4d(fixture_dir, conj_exe=exe, conj_models=models)
    assert "MacSyFinder >= 2.1.6" in str(exc.value)
    assert "MacSyFinder 2.1.4" in str(exc.value)

    models, exe = _conjscan_install(fixture_dir, "2.1", "MacSyFinder 2.1.6")
    out = _preflight_s4d(fixture_dir, conj_exe=exe, conj_models=models)
    assert "conjscan_macsyfinder\t2.1.6" in out.read_text()


@every_tool
def test_preflight_names_a_missing_conjscan_only_when_it_is_required(fixture_dir):
    with pytest.raises(SystemExit, match="CONJScan is not installed"):
        _preflight_s4d(fixture_dir, conj_required=True)
    out = _preflight_s4d(fixture_dir, conj_required=False)
    assert "conjscan\tnot installed" in out.read_text()


def _resolve_fixture(fixture_dir):
    """Two tiers. P1 is named by Pfam at T1; P2 is called hypothetical at T2. The spans
    come from the last tier alone, which writes forward everything it inherited."""
    hits = []
    for pid, label, tier in [("P1", "RepA_N", "T1"), ("P2", "hypothetical protein", "T2")]:
        h = fixture_dir / f"{tier}_hits.tsv"
        write_tsv(h, ["query", "label", "coverage", "target_coverage", "evalue",
                      "informative", "start", "end", "tier", "threshold", "max_evalue"],
                  [[pid, label, 0.9, 0.9, "1e-40", label == "RepA_N", 1, 90, tier,
                    "--cut_ga", ""]])
        hits.append(str(h))
    spans = fixture_dir / "T2_spans.tsv"
    write_tsv(spans, ["seq_id", "qlen", "intervals", "explained_fraction"],
              [["P1", 100, "1-90", 0.9], ["P2", 100, "", 0.0]])

    faa = fixture_dir / "query.faa"
    write_fasta(faa, [("P1", "M" * 100), ("P2", "K" * 100)])
    return hits, str(spans), faa


def test_cascade_resolve_takes_the_explained_fraction_from_the_last_tier(fixture_dir):
    """The last tier's spans hold every protein the cascade ever saw, because each tier
    writes forward what it inherited. A protein missing from them would read as explained
    fraction zero - 'nothing named it' - and the plasmid backbone would enter the
    screening pool."""
    hits, spans, faa = _resolve_fixture(fixture_dir)
    out = fixture_dir / "protein_annotation.tsv"

    run_script("cascade_resolve.py", FakeSnakemake(
        input={"hits": hits, "spans": spans, "faa": str(faa),
               "selection": _selection(fixture_dir), "ps": _ps_table(fixture_dir)},
        output=[str(out)],
        params={"thresholds": {"narrow_at": 0.9, "min_explained": 0.5, "min_coverage": 0.5,
                               "full_at": 0.8, "partial_at": 0.5},
                "tier_order": ["T1", "T2"]}))

    rows = {r["seq_id"]: r for r in read_tsv(out)}
    assert rows["P1"]["functional_class"] == "FUNCTIONAL"
    assert float(rows["P1"]["explained_fraction"]) == 0.9
    assert rows["P2"]["functional_class"] == "UNCHARACTERIZED_HOMOLOG"


def _tiny_diamond_db(fixture_dir):
    """A three-sequence DIAMOND database, built here so the test needs no shared state."""
    import subprocess
    ref = fixture_dir / "ref.faa"
    write_fasta(ref, [
        ("sp|P1|REPA Replication initiator protein RepA", "MKVLATTLLGAAFAASSALAQKKWLVRD"),
        ("sp|P2|HYP hypothetical protein", "MQQTTLNRSDEIVWCAPGHKGGAFLND"),
        ("sp|P3|TOX Toxin RelE", "MRTFEVRLSPQAEKDLDDIYDYIAQD"),
    ])
    db = fixture_dir / "ref.dmnd"
    subprocess.run(f"diamond makedb --in {ref} -d {db} --quiet",
                   shell=True, check=True)
    return db


@requires("diamond")
def test_a_diamond_tier_without_an_evalue_criterion_runs(fixture_dir):
    """A tier may declare max_evalue: null, meaning "the tool's own threshold decides".
    For a DIAMOND tier the flag must then be omitted: `--evalue None` is accepted by
    DIAMOND, returns zero hits and exits 0."""
    db = _tiny_diamond_db(fixture_dir)
    faa = fixture_dir / "q.faa"
    write_fasta(faa, [("Q1", "MKVLATTLLGAAFAASSALAQKKWLVRD")])
    out = fixture_dir / "t" / "hits.tsv"
    out.parent.mkdir()

    run_script("tier_search.py", FakeSnakemake(
        input={"faa": str(faa), "spans": [], "sweep": "",
               "preflight": ""},
        output={"hits": str(out), "unresolved": str(out.parent / "unresolved.faa"),
                "spans": str(out.parent / "spans.tsv")},
        params={"spec": {"id": "T3", "method": "diamond", "source": "swissprot", "db": str(db),
                         "args": "--very-sensitive", "max_evalue": None},
                "narrow_at": 0.9, "hmmer_z": 3497616, "max_target_seqs": 5},
        threads=1))

    rows = read_tsv(out)
    assert rows, "no hits written for a query identical to a database sequence"
    assert rows[0]["query"] == "Q1"
    # Identity, alignment length, bit score and subject length travel with a hit.
    assert float(rows[0]["identity"]) == 1.0
    assert rows[0]["align_length"] and rows[0]["bitscore"] and rows[0]["target_length"]


@requires("diamond")
def test_a_protein_pfam_or_swissprot_named_is_not_searched_in_nr(fixture_dir):
    """skip_if_named_by: a protein an earlier tier named from a listed source keeps that
    name and never reaches this tier; one it could not name is searched."""
    db = _tiny_diamond_db(fixture_dir)
    faa = fixture_dir / "q.faa"
    write_fasta(faa, [("NAMED", "MKVLATTLLGAAFAASSALAQKKWLVRD"),
                      ("OPEN", "MRTFEVRLSPQAEKDLDDIYDYIAQD")])
    earlier = fixture_dir / "T1_hits.tsv"
    write_tsv(earlier, ["query", "label", "informative", "source"],
              [["NAMED", "RepA_N", "True", "pfam"], ["OPEN", "DUF1", "False", "pfam"]])
    out = fixture_dir / "t" / "hits.tsv"
    out.parent.mkdir()

    run_script("tier_search.py", FakeSnakemake(
        input={"faa": str(faa), "spans": [], "sweep": "", "preflight": "",
               "named": [str(earlier)]},
        output={"hits": str(out), "unresolved": str(out.parent / "unresolved.faa"),
                "spans": str(out.parent / "spans.tsv")},
        params={"spec": {"id": "T5", "method": "diamond", "source": "nr", "db": str(db),
                         "args": "--very-sensitive", "max_evalue": 1e-5,
                         "skip_if_named_by": ["pfam", "swissprot"]},
                "narrow_at": 0.9, "hmmer_z": 3497616, "max_target_seqs": 5},
        threads=1))

    assert {r["query"] for r in read_tsv(out)} == {"OPEN"}


@requires("diamond")
def test_a_tier_carries_forward_what_it_could_not_explain(fixture_dir):
    """The cascade's central mechanic. A protein explained past narrow_at stops here; one
    that nothing named must reach the next tier, because a 'hypothetical protein' hit is a
    homolog, not an annotation."""
    db = _tiny_diamond_db(fixture_dir)
    faa = fixture_dir / "q.faa"
    write_fasta(faa, [("named", "MKVLATTLLGAAFAASSALAQKKWLVRD"),
                      ("unnamed", "MQQTTLNRSDEIVWCAPGHKGGAFLND")])
    out = fixture_dir / "t" / "hits.tsv"
    out.parent.mkdir()

    run_script("tier_search.py", FakeSnakemake(
        input={"faa": str(faa), "spans": [], "sweep": "", "preflight": ""},
        output={"hits": str(out), "unresolved": str(out.parent / "unresolved.faa"),
                "spans": str(out.parent / "spans.tsv")},
        params={"spec": {"id": "T3", "method": "diamond", "source": "swissprot", "db": str(db),
                         "args": "--very-sensitive", "max_evalue": 1.0e-5},
                "narrow_at": 0.9, "hmmer_z": 3497616, "max_target_seqs": 5},
        threads=1))

    carried = {l[1:].split()[0]
               for l in (out.parent / "unresolved.faa").read_text().splitlines()
               if l.startswith(">")}
    assert "unnamed" in carried, (
        "a protein whose only hit is 'hypothetical protein' was withheld from deeper "
        "tiers - an uninformative label must never satisfy the narrowing threshold")
    assert "named" not in carried, "a fully explained protein was searched again"


def test_every_informative_label_survives_into_the_resolved_row(fixture_dir):
    """Every informative label reaches the row, not only the winning one."""
    hits = fixture_dir / "hits.tsv"
    write_tsv(hits, ["query", "label", "coverage", "target_coverage", "evalue",
                     "informative", "is_best", "start", "end", "tier", "source",
                     "category", "threshold", "max_evalue"],
              [["P1", "RepA_N", 0.45, 0.9, "1e-20", True, 0, 1, 45, "T1", "pfam", "",
                "--cut_ga", ""],
               ["P1", "MULTISPECIES: replication protein [Enterobacteriaceae]", 0.95,
                0.9, "1e-90", True, 1, 1, 95, "T4", "nr", "", "--very-sensitive",
                "1e-05"]])
    spans = fixture_dir / "spans.tsv"
    write_tsv(spans, ["seq_id", "qlen", "intervals", "explained_fraction"],
              [["P1", 100, "1-95", 0.95]])
    faa = fixture_dir / "query.faa"
    write_fasta(faa, [("P1", "M" * 100)])
    prot = fixture_dir / "protein_annotation.tsv"

    run_script("cascade_resolve.py", FakeSnakemake(
        input={"hits": [str(hits)], "spans": str(spans), "faa": str(faa),
               "selection": _selection(fixture_dir),
               "ps": _ps_table(fixture_dir)},
        output=[str(prot)],
        params={"thresholds": {"narrow_at": 0.9, "min_explained": 0.5, "min_coverage": 0.5,
                               "full_at": 0.8, "partial_at": 0.5},
                "tier_order": ["T1", "T2", "T3", "T4"]}))

    row = read_tsv(prot)[0]
    assert "RepA_N" in row["informative_labels"], (
        "the curated Pfam family is absent from the resolved row")
    assert row["functional_class"] == "FUNCTIONAL", (
        "a fully explained replication initiator must be annotated, not dark")


@requires("diamond")
def test_every_domain_of_a_multi_domain_protein_is_recorded(fixture_dir):
    """Every informative hit of a tier is kept, not only the best one, so that
    n_informative_hits counts domains and explained_fraction can be recomputed from
    hits.tsv as a cross-check of spans.tsv."""
    import subprocess
    ref = fixture_dir / "ref.faa"
    # Two clearly distinct domains, joined into one query protein.
    dom_a = "MKVLATTLLGAAFAASSALAQKKWLVRDGHIYQPLMNE"
    dom_b = "MQQTTLNRSDEIVWCAPGHKGGAFLNDVWRDNPHLAGC"
    write_fasta(ref, [("sp|P1|A Replication initiator RepA_N domain", dom_a),
                      ("sp|P2|B Replication initiator RepA_C domain", dom_b)])
    db = fixture_dir / "ref.dmnd"
    subprocess.run(f"diamond makedb --in {ref} -d {db} --quiet", shell=True, check=True)

    faa = fixture_dir / "q.faa"
    write_fasta(faa, [("P1", dom_a + dom_b)])
    out = fixture_dir / "t" / "hits.tsv"
    out.parent.mkdir()

    run_script("tier_search.py", FakeSnakemake(
        input={"faa": str(faa), "spans": [], "sweep": "", "preflight": ""},
        output={"hits": str(out), "unresolved": str(out.parent / "unresolved.faa"),
                "spans": str(out.parent / "spans.tsv")},
        params={"spec": {"id": "T3", "method": "diamond", "source": "swissprot", "db": str(db),
                         "args": "--very-sensitive", "max_evalue": 1.0e-3},
                "narrow_at": 0.9, "hmmer_z": 3497616, "max_target_seqs": 5},
        threads=1))

    rows = [r for r in read_tsv(out) if r["informative"] == "True"]
    assert len(rows) >= 2, (
        f"only {len(rows)} informative hit(s) recorded for a two-domain protein")
    assert sum(r["is_best"] == "1" for r in rows) == 1, (
        "exactly one hit per protein per tier must be marked as the winning label")

    # The cross-check the docstring promises: spans recomputed from hits.tsv must agree.
    from plasmidann.cascade import explained_fraction
    spans = read_tsv(out.parent / "spans.tsv")[0]
    recomputed = explained_fraction(
        int(spans["qlen"]), [(int(r["start"]), int(r["end"])) for r in rows])
    assert recomputed == float(spans["explained_fraction"]), (
        f"hits.tsv recomputes to {recomputed} but spans.tsv says "
        f"{spans['explained_fraction']} - the stated cross-check cannot be performed")


@requires("diamond")
def test_tier_search_keeps_the_subject_accession(fixture_dir):
    """A DIAMOND title is free text; the accession is the only key a later join to UniProt
    or RefSeq can use. On the NCBI swissprot database used here the title is
    'P62554.1 RecName: Full=Toxin CcdB; ... [Escherichia coli]', with no gene symbol."""
    import subprocess

    subject = fixture_dir / "db.faa"
    prot = ("MQFKVYTYKRESRYRLFVDVQSDIIDTPGRRMVIPLASARLLSDKVSRELYPVVHIGDESW"
            "RMMTTDMASVPVSVIGEEVADLSHRENDIKNAINLMFWGI")
    write_fasta(subject, [("P62554.1", prot)])
    dmnd = fixture_dir / "db.dmnd"
    subprocess.run(f"diamond makedb --in {subject} -d {dmnd} --quiet",
                   shell=True, check=True)

    faa = fixture_dir / "q.faa"
    write_fasta(faa, [("q1", prot)])
    spans = fixture_dir / "spans_in.tsv"
    write_tsv(spans, ["seq_id", "qlen", "intervals", "explained_fraction"], [])
    sweep = fixture_dir / "sweep.txt"
    sweep.write_text("")
    preflight = fixture_dir / "preflight.tsv"
    write_tsv(preflight, ["check", "status"], [["stub", "SUCCESS"]])

    hits = fixture_dir / "hits.tsv"
    run_script("tier_search.py", FakeSnakemake(
        input={"faa": str(faa), "spans": str(spans), "sweep": str(sweep),
               "preflight": str(preflight)},
        output={"hits": str(hits), "unresolved": str(fixture_dir / "un.faa"),
                "spans": str(fixture_dir / "spans_out.tsv")},
        params={"spec": {"id": "T3", "method": "diamond", "source": "swissprot", "db": str(dmnd),
                         "args": "--fast", "max_evalue": 1e-5},
                "narrow_at": 0.9, "hmmer_z": 1000, "max_target_seqs": 5},
        threads=1))

    rows = read_tsv(hits)
    assert rows, "no hits written for a query identical to a database sequence"
    assert "target_accession" in rows[0], "hits.tsv has no target_accession column"
    assert rows[0]["target_accession"] == "P62554.1", (
        f"accession not captured, got {rows[0]['target_accession']!r}")
    assert "CcdB" in rows[0]["label"] or "P62554" in rows[0]["label"], (
        "the title is not kept as the label")


def test_a_family_level_hit_resolves_functional_with_completeness_not_measured(fixture_dir):
    """The pharokka tier reports a family, an annotation and an E-value and no coordinates.
    A protein named by it alone comes out FUNCTIONAL, and what depends on a measured span
    reads as not measured: completeness NOT_MEASURED, not NONE, and meets_min_explained
    empty, not 0."""
    hits = fixture_dir / "phage_hits.tsv"
    write_tsv(hits, ["query", "label", "target_accession", "coverage", "target_coverage",
                     "evalue", "informative", "is_best", "start", "end", "tier", "source",
                     "category", "threshold", "max_evalue"],
              [["P1", "ParA-like partition protein", "phrog_164", "", "", "1e-42", True,
                1, "", "", "T3", "pharokka", "DNA, RNA and nucleotide metabolism", "",
                "1e-05"]])
    spans = fixture_dir / "phage_spans.tsv"
    write_tsv(spans, ["seq_id", "qlen", "intervals", "explained_fraction"],
              [["P1", 250, "", 0.0]])
    faa = fixture_dir / "phage.faa"
    write_fasta(faa, [("P1", "M" * 250)])
    out = fixture_dir / "phage_annotation.tsv"

    run_script("cascade_resolve.py", FakeSnakemake(
        input={"hits": [str(hits)], "spans": str(spans), "faa": str(faa),
               "selection": _selection(fixture_dir),
               "ps": _ps_table(fixture_dir)},
        output=[str(out)],
        params={"thresholds": {"min_coverage": 0.5, "min_explained": 0.5,
                               "narrow_at": 0.9, "full_at": 0.8, "partial_at": 0.5},
                "tier_order": ["T1", "T2", "T3", "T4"]}))

    row = read_tsv(out)[0]
    assert row["functional_class"] == "FUNCTIONAL", row
    assert row["annot_tier"] == "T3"
    assert row["annot_completeness"] == "NOT_MEASURED", (
        f"completeness must say it was not measured, got {row['annot_completeness']!r}")
    assert row["explained_fraction"] == "0.0"
    assert row["meets_min_explained"] == "", (
        "an unmeasured explained fraction must not read as 'below min_explained'")


def test_cascade_resolve_adds_plasmidscope_rows_as_functional(fixture_dir):
    """Every protein must have a row, or annotate_plasmids would give its ORFs no class.
    The skipped ones are FUNCTIONAL, tier PS, with completeness NOT_MEASURED: eggNOG
    gives no span, so nothing measured how much of the protein is explained."""
    hits, spans, faa = _resolve_fixture(fixture_dir)
    ps = _ps_table(fixture_dir, [["P9", "ANNOTATED", "DJ", "COG2026", "ko:K06218", "",
                                  "ParE_toxin", "", "", "Prodigal:2.6", 1],
                                 ["P8", "NONE", "S", "", "", "", "", "", "",
                                  "Prodigal:2.6", 1]])
    out = fixture_dir / "protein_annotation.tsv"

    run_script("cascade_resolve.py", FakeSnakemake(
        input={"hits": hits, "spans": spans, "faa": str(faa),
               "selection": _selection(fixture_dir), "ps": ps},
        output=[str(out)],
        params={"thresholds": {"narrow_at": 0.9, "min_explained": 0.5, "min_coverage": 0.5,
                               "full_at": 0.8, "partial_at": 0.5},
                "tier_order": ["T1", "T2"]}))

    rows = {r["seq_id"]: r for r in read_tsv(out)}
    assert set(rows) == {"P1", "P2", "P9"}, (
        "a PlasmidScope-dark protein must come from the cascade, not from PlasmidScope")
    assert rows["P9"]["functional_class"] == "FUNCTIONAL"
    assert rows["P9"]["annot_tier"] == "PS" and rows["P9"]["annot_label"] == "ParE_toxin"
    assert rows["P9"]["annot_completeness"] == "NOT_MEASURED"
    assert rows["P9"]["explained_fraction"] == ""


def _check_hmmer_z(fixture_dir, hmmer_z):
    """Four unique proteins, of which PlasmidScope may annotate any number (they are
    counted all the same)."""
    unique = fixture_dir / "unique.faa"
    write_fasta(unique, [(f"u{i}", "MKV") for i in range(4)])
    run_script("check_hmmer_z.py", FakeSnakemake(
        input={"unique": str(unique)},
        output=[str(fixture_dir / "hmmer_z_checked.tsv")],
        params={"hmmer_z": hmmer_z}))


def test_hmmer_z_counts_every_unique_protein(fixture_dir):
    """-Z is 4, the unique proteins, whatever the cascade then searches. The searched-set
    count would halve E-values whenever PlasmidScope covers half the set."""
    _check_hmmer_z(fixture_dir, 4)
    with pytest.raises(ValueError, match="Set hmmer_z: 4"):
        _check_hmmer_z(fixture_dir, 6)


@requires("mmseqs")
def test_cascade_selection_searches_representatives_of_families_on_small_plasmids(
        fixture_dir):
    """The family of an unexplained small-plasmid protein is annotated on both
    sides; within it, proteins at >=90% identity over 80% of both lengths share one
    search. Tier 0 proteins and families of large-plasmid proteins alone are not
    searched."""
    import random
    rng = random.Random(11)
    aa = "ACDEFGHIKLMNPQRSTVWY"
    a = "M" + "".join(rng.choice(aa) for _ in range(200))
    near = "".join(c if i % 50 else rng.choice(aa.replace(c, "")) for i, c in enumerate(a))
    far = "".join(c if i % 3 else rng.choice(aa.replace(c, "")) for i, c in enumerate(a))
    other = "M" + "".join(rng.choice(aa) for _ in range(200))
    known = "M" + "".join(rng.choice(aa) for _ in range(200))
    faa = fixture_dir / "unique_proteins.faa"
    write_fasta(faa, [("s1", a), ("l_near", near), ("l_far", far), ("l_other", other),
                      ("s_known", known)])
    pmap = fixture_dir / "protein_map.tsv"
    pmap.write_text("s1\tS1|1\nl_near\tL1|1\nl_far\tL2|1\nl_other\tL3|1\n"
                    "s_known\tS2|1\n")
    small_ids = fixture_dir / "small_plasmids.txt"
    small_ids.write_text("S1\nS2\n")
    families = fixture_dir / "families_broad_cluster.tsv"
    families.write_text("s1\ts1\ns1\tl_near\ns1\tl_far\nl_other\tl_other\n"
                        "s_known\ts_known\n")
    ps = _ps_table(fixture_dir, [["s_known", "ANNOTATED", "", "", "", "", "RHH_1", "", "",
                                  "Prodigal:2.6", 1]])
    out = {"tsv": str(fixture_dir / "selection.tsv"),
           "faa": str(fixture_dir / "search_representatives.faa")}

    run_script("cascade_selection.py", FakeSnakemake(
        input={"faa": str(faa), "map": str(pmap), "small_ids": str(small_ids),
               "ps": ps, "families": str(families),
               "artefact": _artefact_flags(fixture_dir)},
        output=out,
        params={"search": {"min_seq_id": 0.9, "coverage": 0.8, "cov_mode": 0,
                           "cluster_mode": 2}},
        threads=2))

    rows = {r["seq_id"]: r for r in read_tsv(out["tsv"])}
    assert rows["s_known"]["role"] == "plasmidscope"
    assert rows["l_other"]["role"] == "not_selected"
    assert rows["l_far"]["role"] == "representative", "67% identity must not share a search"
    pair = {rows["s1"]["role"], rows["l_near"]["role"]}
    assert pair == {"representative", "member"}
    member = "s1" if rows["s1"]["role"] == "member" else "l_near"
    assert rows[member]["search_representative"] == ({"s1", "l_near"} - {member}).pop()
    assert (rows["s1"]["on_small"], rows["l_near"]["on_small"]) == ("1", "0")
    searched = {l[1:].strip() for l in open(out["faa"]) if l.startswith(">")}
    assert searched == {sid for sid, r in rows.items() if r["role"] == "representative"}


def test_cascade_resolve_gives_members_their_representatives_result(fixture_dir):
    """A search-cluster member takes its representative's row and says so; a protein no
    selected family holds is NOT_SEARCHED, which is neither dark nor annotated."""
    hits, spans, faa = _resolve_fixture(fixture_dir)
    out = fixture_dir / "protein_annotation.tsv"
    selection = _selection(fixture_dir, [["P1", 1, "representative", "P1"],
                                         ["M1", 0, "member", "P1"],
                                         ["X", 0, "not_selected", ""]])

    run_script("cascade_resolve.py", FakeSnakemake(
        input={"hits": hits, "spans": spans, "faa": str(faa), "selection": selection,
               "ps": _ps_table(fixture_dir)},
        output=[str(out)],
        params={"thresholds": {"narrow_at": 0.9, "min_explained": 0.5, "min_coverage": 0.5,
                               "full_at": 0.8, "partial_at": 0.5},
                "tier_order": ["T1", "T2"]}))

    rows = {r["seq_id"]: r for r in read_tsv(out)}
    assert rows["P1"]["annot_source"] == "self"
    member = rows["M1"]
    assert (member["annot_source"], member["annot_representative"]) == (
        "representative", "P1")
    assert (member["functional_class"], member["annot_label"]) == (
        rows["P1"]["functional_class"], rows["P1"]["annot_label"])
    assert (rows["X"]["functional_class"], rows["X"]["annot_source"]) == (
        "NOT_SEARCHED", "not_searched")


@requires("mmseqs")
def test_antifam_flagged_proteins_are_in_no_tier_query_set(fixture_dir):
    """A protein AntiFam flags skips every annotation tier. s_art would otherwise be
    searched (an unexplained small-plasmid protein, its own family's representative);
    here it is not, and it does not open its family for l_rel either. ps_art is annotated
    by Tier 0 and flagged: the flag wins. s_lc is flagged for low complexity only, which
    is no reason to skip, so it is searched."""
    import random
    rng = random.Random(5)
    aa = "ACDEFGHIKLMNPQRSTVWY"
    seqs = {n: "M" + "".join(rng.choice(aa) for _ in range(150))
            for n in ("s_art", "l_rel", "s_lc", "ps_art", "s_ok")}
    faa = fixture_dir / "unique_proteins.faa"
    write_fasta(faa, list(seqs.items()))
    pmap = fixture_dir / "protein_map.tsv"
    pmap.write_text("s_art\tS1|1\nl_rel\tL1|1\ns_lc\tS2|1\nps_art\tS3|1\ns_ok\tS3|2\n")
    small_ids = fixture_dir / "small_plasmids.txt"
    small_ids.write_text("S1\nS2\nS3\n")
    families = fixture_dir / "families_intermediate_cluster.tsv"
    families.write_text("s_art\ts_art\ns_art\tl_rel\ns_lc\ts_lc\nps_art\tps_art\n"
                        "s_ok\ts_ok\n")
    ps = _ps_table(fixture_dir, [["ps_art", "ANNOTATED", "", "", "", "", "RHH_1", "", "",
                                  "Prodigal:2.6", 1]])
    artefact = _artefact_flags(fixture_dir, [
        ["s_art", 1, "AntiFam_ANF00001", "1e-30", 0.0, "antifam"],
        ["ps_art", 1, "AntiFam_ANF00002", "1e-20", 0.6, "antifam,low_complexity"],
        ["s_lc", 1, "", "", 0.7, "low_complexity"],
        ["l_rel", 0, "", "", 0.0, ""], ["s_ok", 0, "", "", 0.0, ""]])
    out = {"tsv": str(fixture_dir / "selection.tsv"),
           "faa": str(fixture_dir / "search_representatives.faa")}

    run_script("cascade_selection.py", FakeSnakemake(
        input={"faa": str(faa), "map": str(pmap), "small_ids": str(small_ids),
               "ps": ps, "families": str(families), "artefact": artefact},
        output=out,
        params={"search": {"min_seq_id": 0.9, "coverage": 0.8, "cov_mode": 0,
                           "cluster_mode": 2}},
        threads=2))

    roles = {r["seq_id"]: r["role"] for r in read_tsv(out["tsv"])}
    assert roles == {"s_art": "artefact_antifam", "ps_art": "artefact_antifam",
                     "l_rel": "not_selected", "s_lc": "representative",
                     "s_ok": "representative"}
    searched = {l[1:].strip() for l in open(out["faa"]) if l.startswith(">")}
    assert searched == {"s_lc", "s_ok"}, "an AntiFam-flagged protein reached the cascade"


def test_an_antifam_skipped_protein_is_not_searched_and_never_dark(fixture_dir):
    """The skipped protein has a row, NOT_SEARCHED with annot_source
    artefact_antifam - not a PlasmidScope row, even when Tier 0 annotated it - and the
    target-eligibility step keeps it out of the dark set."""
    hits, spans, faa = _resolve_fixture(fixture_dir)
    out = fixture_dir / "protein_annotation.tsv"
    selection = _selection(fixture_dir, [["P1", 1, "representative", "P1"],
                                         ["A", 1, "artefact_antifam", ""],
                                         ["PA", 1, "artefact_antifam", ""]])
    ps = _ps_table(fixture_dir, [["PA", "ANNOTATED", "DJ", "COG2026", "ko:K06218", "",
                                  "ParE_toxin", "", "", "Prodigal:2.6", 1]])

    run_script("cascade_resolve.py", FakeSnakemake(
        input={"hits": hits, "spans": spans, "faa": str(faa), "selection": selection,
               "ps": ps},
        output=[str(out)],
        params={"thresholds": {"narrow_at": 0.7, "min_explained": 0.5, "min_coverage": 0.5,
                               "full_at": 0.8, "partial_at": 0.5},
                "tier_order": ["T1", "T2"]}))

    rows = read_tsv(out)
    assert [r["seq_id"] for r in rows].count("PA") == 1, "one row per protein"
    by_id = {r["seq_id"]: r for r in rows}
    for sid in ("A", "PA"):
        assert (by_id[sid]["functional_class"], by_id[sid]["annot_source"],
                by_id[sid]["annot_label"]) == ("NOT_SEARCHED", "artefact_antifam", "")

    artefact = _artefact_flags(fixture_dir, [
        ["A", 1, "AntiFam_ANF00001", "1e-30", 0.0, "antifam"],
        ["PA", 1, "AntiFam_ANF00002", "1e-20", 0.0, "antifam"]])
    flags = fixture_dir / "eligibility.tsv"
    run_script("target_eligibility.py", FakeSnakemake(
        input={"prot": str(out), "artefact": artefact},
        output={"flags": str(flags), "report": str(fixture_dir / "eligibility.txt")}))
    gate = {r["seq_id"]: r for r in read_tsv(flags)}
    for sid in ("A", "PA"):
        assert gate[sid]["target_eligible"] == "0"
        assert gate[sid]["exclusion_reason"] == "artefact,not_searched"


def test_a_dark_protein_named_only_by_a_span_less_family_has_dark_completeness_not_measured(
        fixture_dir):
    """A pharokka family named 'hypothetical protein' carries no span. Its dark coverage
    was not measured, so dark_completeness reads NOT_MEASURED, not NONE."""
    hits = fixture_dir / "phage_hits.tsv"
    write_tsv(hits, ["query", "label", "target_accession", "coverage", "target_coverage",
                     "evalue", "informative", "is_best", "start", "end", "tier", "source",
                     "category", "threshold", "max_evalue"],
              [["P1", "hypothetical protein", "phrog_9928", "", "", "1e-9", False, 0, "",
                "", "T3", "pharokka", "unknown function", "", "1e-05"]])
    spans = fixture_dir / "phage_spans.tsv"
    write_tsv(spans, ["seq_id", "qlen", "intervals", "explained_fraction"],
              [["P1", 90, "", 0.0]])
    faa = fixture_dir / "phage.faa"
    write_fasta(faa, [("P1", "M" * 90)])
    out = fixture_dir / "phage_annotation.tsv"

    run_script("cascade_resolve.py", FakeSnakemake(
        input={"hits": [str(hits)], "spans": str(spans), "faa": str(faa),
               "selection": _selection(fixture_dir), "ps": _ps_table(fixture_dir)},
        output=[str(out)],
        params={"thresholds": {"min_coverage": 0.5, "min_explained": 0.5,
                               "narrow_at": 0.9, "full_at": 0.8, "partial_at": 0.5},
                "tier_order": ["T1", "T2", "T3", "T4"]}))

    row = read_tsv(out)[0]
    assert row["functional_class"] == "UNCHARACTERIZED_HOMOLOG"
    assert row["dark_completeness"] == "NOT_MEASURED", row["dark_completeness"]


def _sweep_cohort(fixture_dir, ids, name):
    faa = fixture_dir / f"{name}.faa"
    write_fasta(faa, [(i, "MKV") for i in ids])
    out = fixture_dir / f"{name}_cohort.txt"
    run_script("sweep_cohort.py", FakeSnakemake(
        input={"faa": str(faa)}, output=[str(out)],
        params={"fraction": 0.2, "seed": 7}))
    return out.read_text().split()


def test_the_sweep_cohort_does_not_depend_on_the_fasta_order(fixture_dir):
    """The query FASTA is written in MMseqs2 representative order, which can change
    between versions and thread counts; the cohort must be the same set regardless."""
    ids = [f"p{i:03d}" for i in range(100)]
    forward = _sweep_cohort(fixture_dir, ids, "forward")
    assert len(forward) == 20
    assert forward == _sweep_cohort(fixture_dir, ids[::-1], "reverse")


def _run_tier(fixture_dir, faa, spec, sweep="", narrow_at=0.9, hmmer_z=1000, tag="t"):
    out = fixture_dir / tag / "hits.tsv"
    out.parent.mkdir()
    run_script("tier_search.py", FakeSnakemake(
        input={"faa": str(faa), "spans": [], "sweep": sweep, "preflight": ""},
        output={"hits": str(out), "unresolved": str(out.parent / "unresolved.faa"),
                "spans": str(out.parent / "spans.tsv")},
        params={"spec": spec, "narrow_at": narrow_at, "hmmer_z": hmmer_z,
                "max_target_seqs": 5},
        threads=1))
    carried = {l[1:].split()[0] for l in (out.parent / "unresolved.faa").read_text()
               .splitlines() if l.startswith(">")}
    return read_tsv(out), read_tsv(out.parent / "spans.tsv"), carried


@requires("diamond")
def test_a_sweep_cohort_protein_is_carried_forward_however_well_explained(fixture_dir):
    """Two proteins identical to a named database sequence are fully explained; only the
    one outside the cohort stops here."""
    db = _tiny_diamond_db(fixture_dir)
    faa = fixture_dir / "q.faa"
    write_fasta(faa, [("in_cohort", "MKVLATTLLGAAFAASSALAQKKWLVRD"),
                      ("outside", "MKVLATTLLGAAFAASSALAQKKWLVRD")])
    sweep = fixture_dir / "sweep.txt"
    sweep.write_text("in_cohort\n")

    _, spans, carried = _run_tier(
        fixture_dir, faa, sweep=str(sweep),
        spec={"id": "T4", "method": "diamond", "source": "swissprot", "db": str(db),
              "args": "--very-sensitive", "max_evalue": 1e-5})

    assert {r["seq_id"]: float(r["explained_fraction"]) for r in spans} == {
        "in_cohort": 1.0, "outside": 1.0}
    assert carried == {"in_cohort"}


@requires("hmmbuild", "hmmsearch")
def test_the_hmmer_tier_reads_domain_coordinates_accession_and_i_evalue(fixture_dir):
    """A 60-residue domain profile against a protein holding the domain at 31-90. The
    coordinates, family name and accession come from the right --domtblout columns, the
    i-Evalue scales with the pinned -Z, and max_evalue rejects the hit when it is weaker."""
    import random
    import subprocess
    rng = random.Random(3)
    aa = "ACDEFGHIKLMNPQRSTVWY"
    domain = "".join(rng.choice(aa) for _ in range(60))
    msa = fixture_dir / "fam.sto"
    rows = ["".join(c if rng.random() > 0.1 else rng.choice(aa) for c in domain)
            for _ in range(4)]
    msa.write_text("# STOCKHOLM 1.0\n#=GF ID TestFam\n#=GF AC PF99999.1\n"
                   + "".join(f"s{i} {r}\n" for i, r in enumerate(rows)) + "//\n")
    hmm = fixture_dir / "fam.hmm"
    subprocess.run(["hmmbuild", str(hmm), str(msa)], check=True,
                   stdout=subprocess.DEVNULL)
    faa = fixture_dir / "q.faa"
    left, right = ("".join(rng.choice(aa) for _ in range(30)) for _ in range(2))
    write_fasta(faa, [("P1", left + domain + right)])
    spec = {"id": "T1", "method": "hmmer", "source": "pfam", "db": str(hmm),
            "args": "--domE 10 -E 10", "max_evalue": None}

    hits, spans, _ = _run_tier(fixture_dir, faa, spec, hmmer_z=1, tag="z1")
    (hit,) = hits
    assert (hit["query"], hit["label"], hit["target_accession"]) == (
        "P1", "TestFam", "PF99999.1")
    assert abs(int(hit["start"]) - 31) <= 3 and abs(int(hit["end"]) - 90) <= 3, hit
    assert hit["target_length"] == "60" and spans[0]["qlen"] == "120"

    hits_z, _, _ = _run_tier(fixture_dir, faa, spec, hmmer_z=1000, tag="z1000")
    ratio = float(hits_z[0]["evalue"]) / float(hit["evalue"])
    assert 900 < ratio < 1100, f"the i-Evalue did not scale with -Z: ratio {ratio}"

    strict = {**spec, "max_evalue": float(hit["evalue"]) / 10}
    rejected, _, _ = _run_tier(fixture_dir, faa, strict, hmmer_z=1, tag="strict")
    assert rejected == []


def test_the_pharokka_tier_turns_the_merged_table_into_span_less_hits(fixture_dir):
    """A script standing in for pharokka writes its merged table: the family hit
    reaches hits.tsv with empty coordinates, the protein length comes from the table, and
    an unexplained protein is carried forward."""
    from test_pharokka import BOTH, HEADER
    table = fixture_dir / "merged.tsv"
    table.write_text(HEADER + BOTH)
    exe = fixture_dir / "bin" / "pharokka"
    exe.parent.mkdir()
    exe.write_text("#!/bin/sh\nwhile [ $# -gt 0 ]; do [ \"$1\" = -o ] && out=$2; shift; "
                   f"done\nmkdir -p \"$out\"\ncp {table} \"$out/tier_full_merged_output.tsv\"\n")
    exe.chmod(0o755)
    faa = fixture_dir / "q.faa"
    write_fasta(faa, [("0cab68e6", "M" * 250)])

    hits, spans, carried = _run_tier(
        fixture_dir, faa, {"id": "T3", "method": "pharokka", "source": "pharokka",
                           "db": str(fixture_dir), "exe": str(exe), "args": "",
                           "max_evalue": 1e-5})

    (hit,) = hits
    assert (hit["label"], hit["target_accession"], hit["source"]) == (
        "ParA-like partition protein", "phrog_164", "pharokka")
    assert (hit["start"], hit["end"], hit["informative"]) == ("", "", "True")
    assert spans == [{"seq_id": "0cab68e6", "qlen": "250", "intervals": "",
                      "explained_fraction": "0.0"}]
    assert carried == {"0cab68e6"}
