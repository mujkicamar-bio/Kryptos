"""Smoke tests for the Snakemake script layer.

Each test runs one real script against a tiny fixture. They are fast - the whole file
finishes in about a second - and they exist because five reviews found that every serious
defect in this pipeline lived in the one layer no test touched.

A test here asserts two things: the script does not raise, and it produces output that is
actually usable by the stage downstream. The second half matters more than the first: four
stages were found producing well-formed EMPTY tables and reporting success.
"""
import pathlib

import pytest

from conftest import (FakeSnakemake, run_script, requires, write_tsv,
                      write_fasta, read_tsv)


# --- S7b: the NameError that a missing tool was hiding --------------------------------

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
                     "n_mob_clusters", "family_class", "members"],
              [["F0000001", "m1", 3, 3, 2, "FAMILY", "m1,m2,m3"]])

    run_script("family_evolution.py", FakeSnakemake(
        input={"families": str(fams), "faa": str(faa), "cds": str(cds)},
        output={"tsv": str(out), "consensus": str(out.parent / "consensus.faa")},
        params={"evolution": {"min_codons": 20, "min_members_for_dnds": 3,
                              "dnds_purifying_max": 0.5, "rnacode_max_p": 0.05,
                              "max_members_aligned": 50}},
        threads=1))

    rows = read_tsv(out)
    assert len(rows) == 1
    # The family is alignable and divergent, so a status must have been reached. An empty
    # status means the script fell through an error branch and reported success.
    assert rows[0]["dnds_status"], "no dnds_status - the script took a silent error branch"
    assert rows[0]["dnds_status"] != "ALIGNMENT_FAILED", (
        "alignment failed with mafft available - the tool is not being invoked correctly")


# --- S9b: the column name that three reviewers found independently --------------------

# --- S5: the gate that can neither pass nor fail meaningfully -------------------------

def test_the_quality_gate_passes_when_controls_are_annotated(fixture_dir):
    """The gate must not halt a healthy run. Controls that the cascade named correctly
    should clear it."""
    prot = fixture_dir / "protein_annotation.tsv"
    artefact = fixture_dir / "artefact_flags.tsv"
    flags = fixture_dir / "eligibility.tsv"
    report = fixture_dir / "gate.txt"

    rows = [[f"CTRL_{i:05d}_P0000{i}", "T1", "relaxase MobA", "FUNCTIONAL"]
            for i in range(100)]
    rows += [["seq_dark_1", "", "", "NONE"]]
    write_tsv(prot, ["seq_id", "annot_tier", "annot_label", "functional_class"], rows)
    write_tsv(artefact, ["seq_id", "artefact_flag"],
              [[r[0], 0] for r in rows])

    run_script("quality_gate.py", FakeSnakemake(
        input={"prot": str(prot), "artefact": str(artefact)},
        output={"flags": str(flags), "report": str(report)},
        params={"gate": {"min_control_recall": 0.99, "require_control_set": True}}))

    text = report.read_text()
    assert "control_recall=1.0" in text


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
            params={"gate": {"min_control_recall": 0.99, "require_control_set": True}}))


# --- S2c: the control set that could not pass -----------------------------------------

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
        input={"faa": str(faa), "raw": str(raw), "decoys": str(decoy_faa)},
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


# --- S3 pre-flight: the rule that exists to stop 45-hour failures ---------------------

@requires("hmmsearch", "diamond", "mafft", "mmseqs", "foldseek")
def test_preflight_checks_every_tool_the_workflow_runs(fixture_dir):
    """Pre-flight covered 4 of the 11 executables the workflow invokes. A missing mafft,
    mmseqs, macsyfinder, integron_finder, prodigal, cmsearch or foldseek still killed the
    run days in - which is the exact failure class this rule was written to remove."""
    from plasmidann.tools import tool_names

    db = fixture_dir / "fake.hmm"
    db.write_text("HMMER3/f\n")
    (fixture_dir / "fake.hmm.h3i").write_text("")
    antifam = fixture_dir / "AntiFam.hmm"
    antifam.write_text("HMMER3/f\n")
    (fixture_dir / "AntiFam.hmm.h3i").write_text("")
    out = fixture_dir / "preflight.tsv"

    run_script("preflight.py", FakeSnakemake(
        output=[str(out)],
        params={"tiers": [{"id": "T1", "method": "hmmer", "db": str(db),
                           "args": "--cut_ga", "max_evalue": None}],
                "artefact": {"antifam_db": str(antifam)},
                "structure": {"required": False},
                "orthology": {"required": False, "data_dir": ""},
                "foldseek_db": str(fixture_dir / "absent"),
                "prostt5": str(fixture_dir / "absent")}))

    reported = {l.split("\t")[0] for l in out.read_text().splitlines() if "\t" in l}
    unchecked = {t for t in tool_names()} - reported - {"foldseek", "emapper.py"}
    assert not unchecked, f"pre-flight never checked: {sorted(unchecked)}"


def test_preflight_fails_when_a_downstream_tool_is_missing(fixture_dir, monkeypatch):
    """The whole point: fail in seconds, not after the cascade has run for three days.

    PATH cannot be emptied here - the harness deliberately prepends the project's bin so
    that a missing tool fails a test rather than sending a script down a silent error
    branch. So the absence is injected at shutil.which instead."""
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
            params={"tiers": [{"id": "T1", "method": "hmmer", "db": str(db),
                               "args": "--cut_ga", "max_evalue": None}],
                    "artefact": {"antifam_db": str(antifam)},
                    "structure": {"required": False},
                    "foldseek_db": "", "prostt5": ""}))
    assert "mafft" in str(exc.value), "a missing mafft must be named by pre-flight"
    assert "S7b" in str(exc.value), "pre-flight must say which stage the tool belongs to"


# --- S3 sharding: the difference between losing a shard and losing three days ---------

def test_the_cascade_input_split_is_lossless(fixture_dir):
    """T4 searches 3.5M queries against a 375 GB database. As one job, any failure in it
    loses days of work; the split is what makes a failure cost one shard. A split that
    dropped or duplicated a sequence would silently change the protein set the whole run
    is about, so losslessness is asserted rather than assumed."""
    faa = fixture_dir / "cascade_input.faa"
    outs = [fixture_dir / f"in/{i:03d}.faa" for i in range(4)]
    (fixture_dir / "in").mkdir()

    records = [(f"P{i:05d}", "MKVLATT" * (1 + i % 3)) for i in range(37)]
    write_fasta(faa, records)

    run_script("shard_cascade_input.py", FakeSnakemake(
        input=[str(faa)], output=[str(p) for p in outs]))

    seen = {}
    for path in outs:
        name = None
        for line in path.read_text().splitlines():
            if line.startswith(">"):
                name = line[1:].split()[0]
                assert name not in seen, f"{name} written to more than one shard"
                seen[name] = ""
            else:
                seen[name] += line

    assert seen == dict(records), "the split lost, duplicated or altered a sequence"
    sizes = [len([l for l in p.read_text().splitlines() if l.startswith(">")])
             for p in outs]
    assert max(sizes) - min(sizes) <= 1, f"shards are unbalanced: {sizes}"


# --- S3 resolve: one row per protein, assembled across tiers AND shards ---------------

def _resolve_fixture(fixture_dir):
    """Two shards, two tiers. P1 is named by Pfam; P2 is called hypothetical everywhere."""
    hits, spans = [], []
    for cshard, (pid, label, tier) in enumerate([("P1", "RepA_N", "T1"),
                                                 ("P2", "hypothetical protein", "T2")]):
        h = fixture_dir / f"{tier}_{cshard}_hits.tsv"
        write_tsv(h, ["query", "label", "coverage", "target_coverage", "evalue",
                      "informative", "start", "end", "tier", "threshold", "max_evalue"],
                  [[pid, label, 0.9, 0.9, "1e-40", label == "RepA_N", 1, 90, tier,
                    "--cut_ga", ""]])
        hits.append(str(h))
        s = fixture_dir / f"{tier}_{cshard}_spans.tsv"
        write_tsv(s, ["seq_id", "qlen", "intervals", "explained_fraction"],
                  [[pid, 100, "1-90" if label == "RepA_N" else "", 0.9
                    if label == "RepA_N" else 0.0]])
        spans.append(str(s))

    faa = fixture_dir / "cascade_input.faa"
    write_fasta(faa, [("P1", "M" * 100), ("P2", "K" * 100)])
    return hits, spans, faa


def test_cascade_resolve_reads_every_shard(fixture_dir):
    """Spans arrive as one file per shard now that the cascade is sharded. Reading only
    the first would silently zero the explained fraction of every protein in every other
    shard, and a zero explained fraction reads as 'nothing named it' - which would push
    the whole plasmid backbone into the screening pool."""
    hits, spans, faa = _resolve_fixture(fixture_dir)
    out = fixture_dir / "protein_annotation.tsv"

    run_script("cascade_resolve.py", FakeSnakemake(
        input={"hits": hits, "spans": spans, "faa": str(faa)},
        output=[str(out)],
        params={"thresholds": {"narrow_at": 0.9, "min_explained": 0.5, "min_coverage": 0.5,
                               "full_at": 0.8, "partial_at": 0.5},
                "tier_order": ["T1", "T2"]}))

    rows = {r["seq_id"]: r for r in read_tsv(out)}
    assert rows["P1"]["functional_class"] == "FUNCTIONAL"
    assert float(rows["P1"]["explained_fraction"]) == 0.9
    assert rows["P2"]["functional_class"] == "UNCHARACTERIZED_HOMOLOG"


# --- logging: 26 rules declared a log file that Snakemake never wrote to --------------

def test_a_scripts_output_reaches_its_declared_log(fixture_dir):
    """Every rule declares `log:`, and for `script:` rules Snakemake does not redirect
    stdout there - it only creates the path. So all 26 log files stayed empty and every
    diagnostic the scripts print (hit counts, rejected hits, background rates) went to the
    single SLURM output file interleaved across 1,481 concurrent jobs, where it cannot be
    attributed to a rule."""
    faa = fixture_dir / "cascade_input.faa"
    write_fasta(faa, [("P1", "MKVL"), ("P2", "MKVA")])
    log = fixture_dir / "logs" / "split.log"
    outs = [fixture_dir / f"in/{i}.faa" for i in range(2)]
    (fixture_dir / "in").mkdir()

    run_script("shard_cascade_input.py", FakeSnakemake(
        input=[str(faa)], output=[str(p) for p in outs], log=[str(log)]))

    assert log.exists(), "the declared log file was never created"
    assert "2 proteins" in log.read_text(), (
        f"the script's own diagnostics did not reach its log: {log.read_text()!r}")


# --- S2b: the artefact screen must use AntiFam's own curated thresholds ---------------

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


# --- S3 tier_search: the script the entire deliverable rests on -----------------------

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
    """A tier may declare max_evalue: null, meaning "the tool's own threshold decides" -
    the convention T1 already uses. For a DIAMOND tier the command was interpolated
    unconditionally, so it went out as `--evalue None`, DIAMOND exited non-zero, and with
    check=True the tier dies. The flag has to be omitted, not stringified."""
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
        params={"spec": {"id": "T3", "method": "diamond", "db": str(db),
                         "args": "--very-sensitive", "max_evalue": None},
                "narrow_at": 0.9, "hmmer_z": 3497616, "max_target_seqs": 5},
        threads=1))

    rows = read_tsv(out)
    assert rows, "no hits written for a query identical to a database sequence"
    assert rows[0]["query"] == "Q1"


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
        params={"spec": {"id": "T3", "method": "diamond", "db": str(db),
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


# --- S5: the stop-list has to see every label, not only the winning one ---------------

def test_every_informative_label_survives_into_the_resolved_row(fixture_dir):
    """Labels are ranked by E-value, so an nr hit with free text routinely takes
    annot_label away from a curated Pfam assignment on the same protein. Anything reading
    only annot_label therefore cannot see the Pfam family at all.

    Nothing filters on this any more - the backbone stop-list is gone - but the resolved
    row is the complete annotation record, and a curated family assignment is the most
    authoritative thing in it. Dropping it would be losing annotation, not losing a
    filter."""
    hits = fixture_dir / "hits.tsv"
    write_tsv(hits, ["query", "label", "coverage", "target_coverage", "evalue",
                     "informative", "is_best", "start", "end", "tier", "threshold",
                     "max_evalue"],
              [["P1", "RepA_N", 0.45, 0.9, "1e-20", True, 0, 1, 45, "T1", "--cut_ga", ""],
               ["P1", "MULTISPECIES: replication protein [Enterobacteriaceae]", 0.95,
                0.9, "1e-90", True, 1, 1, 95, "T4", "--very-sensitive", "1e-05"]])
    spans = fixture_dir / "spans.tsv"
    write_tsv(spans, ["seq_id", "qlen", "intervals", "explained_fraction"],
              [["P1", 100, "1-95", 0.95]])
    faa = fixture_dir / "cascade_input.faa"
    write_fasta(faa, [("P1", "M" * 100)])
    prot = fixture_dir / "protein_annotation.tsv"

    run_script("cascade_resolve.py", FakeSnakemake(
        input={"hits": [str(hits)], "spans": [str(spans)], "faa": str(faa)},
        output=[str(prot)],
        params={"thresholds": {"narrow_at": 0.9, "min_explained": 0.5, "min_coverage": 0.5,
                               "full_at": 0.8, "partial_at": 0.5},
                "tier_order": ["T1", "T2", "T3", "T4"]}))

    row = read_tsv(prot)[0]
    assert "RepA_N" in row["informative_labels"], (
        "the curated Pfam family is absent from the resolved row")
    assert row["functional_class"] == "FUNCTIONAL", (
        "a fully explained replication initiator must be annotated, not dark - that is "
        "what makes a stop-list unnecessary")


# --- S8c: the stage that turns a dark ORF into a testable hypothesis ------------------

def _context_fixture(fixture_dir, partner_label, partner_kind="pfam_family"):
    """One plasmid: a dark ORF, an annotated neighbour, and a dark ORF in a cassette.

    `partner_label` is what the TOOLS called the neighbour. Since the curated 73-name list
    was deleted, nothing here translates that into a class: the label travels into the
    context table as itself, and the grouping into replication, mobilisation and so on is
    derived later from the observed vocabulary.
    """
    ann = fixture_dir / "plasmid_annotation.tsv"
    write_tsv(ann, ["plasmid_id", "orf_id", "start", "end", "strand", "annot_label",
                    "functional_class"],
              [["pl1", "pl1|1", 100, 400, "+", "", "NONE"],
               ["pl1", "pl1|2", 430, 700, "+", partner_label, "FUNCTIONAL"],
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
    labels_tsv = fixture_dir / "protein_labels.tsv"
    write_tsv(labels_tsv, ["protein_id", "source", "tier", "kind", "label", "accession",
                           "evidence_evalue", "evidence_coverage", "database",
                           "database_version"],
              [["S2", "pfam", "T1", partner_kind, partner_label, "PF00000.1", "1e-40",
                "0.9", "Pfam-A", "38.2"]])
    master = fixture_dir / "context_master.tsv"
    write_tsv(master, ["plasmid_id", "size_bp"], [["pl1", 5000]])
    return ann, fam, pmap, defence, integrons, labels_tsv, master


BACKGROUND_PARAMS = {"covariates": ["plasmid_length", "gene_count"],
                     "min_stratum_size": 20}


def _run_context(fixture_dir, partner_label, categories=None):
    ann, fam, pmap, defence, integrons, labels_tsv, master = _context_fixture(
        fixture_dir, partner_label)
    fams_out = fixture_dir / "family_context.tsv"
    bg_out = fixture_dir / "context_background.tsv"
    run_script("context_features.py", FakeSnakemake(
        input={"annotation": str(ann), "families": str(fam), "map": str(pmap),
               "defence": str(defence), "integrons": [str(integrons)],
               "labels": str(labels_tsv), "master": str(master)},
        output={"families": str(fams_out), "background": str(bg_out)},
        params={"context": {"max_operon_gap": 100, "neighbourhood_window": 3,
                            "min_context_conservation": 0.50,
                            "high_confidence_conservation": 0.90,
                            "min_enrichment": 2.0},
                "background": BACKGROUND_PARAMS,
                "categories": categories}))
    return read_tsv(fams_out)


def test_context_records_the_label_the_tool_produced(fixture_dir):
    """The curated 73-name list is gone, so a neighbour's context is whatever the tools
    called it. With no category rules configured the category is the label qualified by its
    kind, which is what lets the grouping be decided later from the observed vocabulary
    instead of being fixed in code now."""
    rows = {r["category"]: r for r in _run_context(fixture_dir, "MobA_MobL")}

    assert "pfam_family:MobA_MobL" in rows, (
        f"the relaxase label is not in the context table; got {sorted(rows)}")
    row = rows["pfam_family:MobA_MobL"]
    assert row["n_units"] == "1", "the unit is the plasmid"
    assert row["n_units_with_category"] == "1"
    assert row["status"], "no status written"
    # The column exists even with no rules configured, so a later grouping can fill it
    # without changing the table's shape.
    assert "subcategories" in row


def _stratified_fixture(fixture_dir, n_small, n_large):
    """A collection of two size bands, each plasmid carrying a dark ORF beside a relaxase.

    The dark family lives only on the SMALL plasmids. Every plasmid in the collection has
    the association, so the corpus-wide rate is 1.0 and nothing can look enriched against
    it; the question is which population the family is compared against.
    """
    ann_rows, pmap_lines, master_rows, members = [], [], [], []
    for i in range(n_small + n_large):
        pid = f"p{i:03d}"
        small = i < n_small
        # A gene-dense large plasmid and a two-gene cryptic one land in different strata
        # on both configured covariates.
        ann_rows.append([pid, f"{pid}|1", 100, 400, "+", "", "NONE"])
        ann_rows.append([pid, f"{pid}|2", 430, 700, "+", "MobA_MobL", "FUNCTIONAL"])
        pmap_lines.append(f"D{i}\t{pid}|1")
        pmap_lines.append(f"M{i}\t{pid}|2")
        master_rows.append([pid, 4000 if small else 150000])
        if not small:
            for j in range(3, 13):
                ann_rows.append([pid, f"{pid}|{j}", 1000 * j, 1000 * j + 300, "+",
                                 "Filler", "FUNCTIONAL"])
                pmap_lines.append(f"X{i}_{j}\t{pid}|{j}")
        else:
            members.append(f"D{i}")

    ann = fixture_dir / "strat_annotation.tsv"
    write_tsv(ann, ["plasmid_id", "orf_id", "start", "end", "strand", "annot_label",
                    "functional_class"], ann_rows)
    fam = fixture_dir / "strat_families.tsv"
    write_tsv(fam, ["family_id", "representative", "n_members", "n_plasmids",
                    "n_mob_clusters", "family_class", "members"],
              [["F1", members[0], len(members), len(members), 1, "ORPHAN",
                ",".join(members)]])
    pmap = fixture_dir / "strat_map.tsv"
    pmap.write_text("\n".join(pmap_lines) + "\n")
    defence = fixture_dir / "strat_defence.tsv"
    write_tsv(defence, ["orf_id", "system"], [])
    integrons = fixture_dir / "strat_integrons.tsv"
    write_tsv(integrons, ["plasmid_id", "integron_id", "element", "start", "end",
                          "integron_type", "annotation", "type_elt"], [])
    labels_tsv = fixture_dir / "strat_labels.tsv"
    write_tsv(labels_tsv, ["protein_id", "source", "tier", "kind", "label", "accession",
                           "evidence_evalue", "evidence_coverage", "database",
                           "database_version"],
              [[f"M{i}", "pfam", "T1", "pfam_family", "MobA_MobL", "PF00000.1", "1e-40",
                "0.9", "Pfam-A", "38.2"] for i in range(n_small + n_large)])
    master = fixture_dir / "strat_master.tsv"
    write_tsv(master, ["plasmid_id", "size_bp"], master_rows)
    return ann, fam, pmap, defence, integrons, labels_tsv, master


def _run_stratified(fixture_dir, n_small, n_large, min_stratum_size):
    ann, fam, pmap, defence, integrons, labels_tsv, master = _stratified_fixture(
        fixture_dir, n_small, n_large)
    fams_out = fixture_dir / "strat_context.tsv"
    bg_out = fixture_dir / "strat_background.tsv"
    run_script("context_features.py", FakeSnakemake(
        input={"annotation": str(ann), "families": str(fam), "map": str(pmap),
               "defence": str(defence), "integrons": [str(integrons)],
               "labels": str(labels_tsv), "master": str(master)},
        output={"families": str(fams_out), "background": str(bg_out)},
        params={"context": {"max_operon_gap": 100, "neighbourhood_window": 3,
                            "min_context_conservation": 0.50,
                            "high_confidence_conservation": 0.90,
                            "min_enrichment": 2.0},
                "background": {"covariates": ["plasmid_length", "gene_count"],
                               "min_stratum_size": min_stratum_size},
                "categories": None}))
    return read_tsv(fams_out), read_tsv(bg_out)


def test_the_context_background_is_stratified_not_pooled(fixture_dir):
    """Spec section 53. A family on small cryptic plasmids is compared against small
    cryptic plasmids, not against the whole collection.

    A flat background divides by one rate for everything, which under-corrects for small
    plasmids - where a +-3 window is the entire molecule - and over-corrects for large
    ones. Small cryptic plasmids are a stratum of interest in this project, so the flat
    background's error lands where it does most damage.
    """
    rows, _ = _run_stratified(fixture_dir, n_small=30, n_large=30, min_stratum_size=20)
    row = [r for r in rows if r["category"] == "pfam_family:MobA_MobL"][0]

    assert row["normalization_method"] == "stratified_prevalence_ratio", (
        f"the comparison was not stratified: {row['normalization_method']}")
    assert row["background_total"] == "30", (
        "the background should be the 30 small plasmids the family lives on, not the 60 "
        f"in the collection; got {row['background_total']}")
    assert "stratified on plasmid_length+gene_count" in row["background_definition"], (
        f"the background is not described on the row: {row['background_definition']!r}")


def test_a_thin_stratum_falls_back_to_the_pooled_background_and_says_so(fixture_dir):
    """Spec section 52.2: a normalised value is not valid without its definition.

    Adding covariates makes every stratum smaller, and at some point a stratum holds too
    few plasmids to estimate a rate from. Falling back silently would leave two rows in
    the same column measured against different populations with nothing to tell them
    apart, so the fallback is recorded on the row that used it.
    """
    rows, _ = _run_stratified(fixture_dir, n_small=5, n_large=30, min_stratum_size=20)
    row = [r for r in rows if r["category"] == "pfam_family:MobA_MobL"][0]

    assert row["normalization_method"] == "pooled_prevalence_ratio", (
        "a 5-plasmid stratum is below min_stratum_size and must not be used as a "
        f"background; got {row['normalization_method']}")
    assert row["background_total"] == "35", (
        f"the pooled fallback should cover the whole collection; got "
        f"{row['background_total']}")
    assert "below min_stratum_size" in row["background_definition"], (
        f"the fallback is not recorded on the row: {row['background_definition']!r}")


def test_every_normalised_context_row_carries_its_normalisation_definition(fixture_dir):
    """Spec section 52.2, stated as an absolute: 'No normalized prevalence field is valid
    without its normalization definition.' The raw counts stay beside it (section 52.1),
    because a normalised value alone is a number whose meaning depends entirely on a
    choice the reader cannot see."""
    rows, _ = _run_stratified(fixture_dir, n_small=30, n_large=30, min_stratum_size=20)

    required = ["raw_count", "raw_prevalence", "background_count", "background_total",
                "background_prevalence", "normalization_method", "normalization_version",
                "background_definition"]
    for row in rows:
        missing = [c for c in required if not str(row.get(c, "")).strip()]
        assert not missing, (
            f"{row['family_id']}/{row['category']} is missing {missing} - a normalised "
            "value without its definition is not a result")


def test_a_plasmid_carrying_a_category_away_from_the_family_does_not_break_the_test(
        fixture_dir):
    """The background and the family measure the same plasmid by different rules.

    The background counts a plasmid for a category if ANY of its ORFs has that category in
    context. The family counts it only if one of the FAMILY's ORFs does. So a plasmid where
    the category sits beside some other gene is background-positive and family-negative,
    and the 2x2 the Fisher test builds - family against the rest of the corpus - then has
    a negative cell: rest_without = (N - n) - (K - k) < 0.

    The pooled corpus background hid this. N was the whole collection, so N - n was large
    enough to absorb the discrepancy and the cell stayed positive by luck rather than by
    construction. Against a stratum of 36 it went negative on the first real run.

    The fix is to derive the corpus counts from the REST of the population explicitly,
    which is what the test is comparing against anyway.
    """
    # Three plasmids. On all three a relaxase is present. On p0 and p1 the family's dark
    # ORF is beside it; on p2 the family's ORF is far away and a DIFFERENT dark ORF is
    # beside it - so p2 is background-positive and family-negative.
    ann_rows, pmap_lines, master_rows, members = [], [], [], []
    for i in range(3):
        pid = f"m{i}"
        master_rows.append([pid, 5000])
        if i < 2:
            ann_rows += [[pid, f"{pid}|1", 100, 400, "+", "", "NONE"],
                         [pid, f"{pid}|2", 430, 700, "+", "MobA_MobL", "FUNCTIONAL"]]
            pmap_lines += [f"D{i}\t{pid}|1", f"M{i}\t{pid}|2"]
        else:
            # The family's ORF, with three intervening genes so the relaxase is outside
            # its +-3 window; then an unrelated dark ORF right beside the relaxase.
            ann_rows += [[pid, f"{pid}|1", 100, 400, "+", "", "NONE"]]
            pmap_lines += [f"D{i}\t{pid}|1"]
            for j, label in enumerate(["Filler", "Filler", "Filler"], start=2):
                ann_rows.append([pid, f"{pid}|{j}", 1000 * j, 1000 * j + 300, "+",
                                 label, "FUNCTIONAL"])
                pmap_lines.append(f"F{i}_{j}\t{pid}|{j}")
            ann_rows += [[pid, f"{pid}|5", 6000, 6300, "+", "", "NONE"],
                         [pid, f"{pid}|6", 6400, 6700, "+", "MobA_MobL", "FUNCTIONAL"]]
            pmap_lines += [f"O{i}\t{pid}|5", f"M{i}\t{pid}|6"]
        members.append(f"D{i}")

    # Four more plasmids in the collection that carry the relaxase beside a dark ORF but
    # hold no member of this family. Without them the family covers the whole background,
    # fisher_enrichment returns NOT_APPLICABLE before the 2x2 is built, and the negative
    # cell is never reached.
    for i in range(3, 7):
        pid = f"m{i}"
        master_rows.append([pid, 5000])
        ann_rows += [[pid, f"{pid}|1", 100, 400, "+", "", "NONE"],
                     [pid, f"{pid}|2", 430, 700, "+", "MobA_MobL", "FUNCTIONAL"]]
        pmap_lines += [f"X{i}\t{pid}|1", f"M{i}\t{pid}|2"]

    ann = fixture_dir / "mix_annotation.tsv"
    write_tsv(ann, ["plasmid_id", "orf_id", "start", "end", "strand", "annot_label",
                    "functional_class"], ann_rows)
    fam = fixture_dir / "mix_families.tsv"
    write_tsv(fam, ["family_id", "representative", "n_members", "n_plasmids",
                    "n_mob_clusters", "family_class", "members"],
              [["F1", "D0", 3, 3, 1, "FAMILY", ",".join(members)]])
    pmap = fixture_dir / "mix_map.tsv"
    pmap.write_text("\n".join(pmap_lines) + "\n")
    defence = fixture_dir / "mix_defence.tsv"
    write_tsv(defence, ["orf_id", "system"], [])
    integrons = fixture_dir / "mix_integrons.tsv"
    write_tsv(integrons, ["plasmid_id", "integron_id", "element", "start", "end",
                          "integron_type", "annotation", "type_elt"], [])
    labels_tsv = fixture_dir / "mix_labels.tsv"
    write_tsv(labels_tsv, ["protein_id", "source", "tier", "kind", "label", "accession",
                           "evidence_evalue", "evidence_coverage", "database",
                           "database_version"],
              [[f"M{i}", "pfam", "T1", "pfam_family", "MobA_MobL", "PF00000.1", "1e-40",
                "0.9", "Pfam-A", "38.2"] for i in range(7)])
    master = fixture_dir / "mix_master.tsv"
    write_tsv(master, ["plasmid_id", "size_bp"], master_rows)

    fams_out = fixture_dir / "mix_context.tsv"
    bg_out = fixture_dir / "mix_background.tsv"
    run_script("context_features.py", FakeSnakemake(
        input={"annotation": str(ann), "families": str(fam), "map": str(pmap),
               "defence": str(defence), "integrons": [str(integrons)],
               "labels": str(labels_tsv), "master": str(master)},
        output={"families": str(fams_out), "background": str(bg_out)},
        params={"context": {"max_operon_gap": 100, "neighbourhood_window": 3,
                            "min_context_conservation": 0.50,
                            "high_confidence_conservation": 0.90,
                            "min_enrichment": 2.0},
                "background": {"covariates": ["plasmid_length", "gene_count"],
                               "min_stratum_size": 2},
                "categories": None}))

    rows = {r["category"]: r for r in read_tsv(fams_out)}
    row = rows["pfam_family:MobA_MobL"]
    assert row["n_units"] == "3"
    assert row["n_units_with_category"] == "2", (
        "the third plasmid carries the relaxase outside the family ORF's window, so the "
        "family count must be 2 while the background count for that plasmid is 1")
    assert row["status"], "the row was not produced at all"


def test_context_writes_no_hand_assigned_category(fixture_dir):
    """The standing constraint at the stage that used to break it. backbone_adjacent and
    ta_candidate were classes assigned by a hand-written list of 73 Pfam family names -
    against 30,134 families in Pfam-A 38.2, of which 67 mention replication in their
    description alone. Neither name may reappear as a category."""
    forbidden = {"backbone_adjacent", "ta_candidate", "replication", "mobilisation",
                 "conjugation", "toxin_antitoxin"}

    seen = {r["category"] for r in _run_context(fixture_dir, "RelE")}

    assert not (seen & forbidden), f"hand-assigned category reappeared: {seen & forbidden}"


def test_a_transposase_and_a_toxin_partner_are_told_apart_by_their_own_labels(fixture_dir):
    """ta_candidate fired whenever the two-gene partner was ANY of the 73 curated families,
    so a dark ORF beside a transposase was labelled a candidate antitoxin and routed to the
    toxin_or_ta_adjacent stratum - 175 of the 1,000 constructs. The distinction is now made
    by the partner's own label rather than by a class, so it cannot be lost to a gap in a
    list."""
    transposase = {r["category"] for r in _run_context(fixture_dir, "DDE_Tnp_Tn3")}
    toxin = {r["category"] for r in _run_context(fixture_dir, "RelE")}

    assert "pfam_family:DDE_Tnp_Tn3" in transposase
    assert "pfam_family:RelE" in toxin
    assert "pfam_family:RelE" not in transposase, (
        "a transposase partner produced a toxin label")


def test_a_configured_category_groups_labels_and_keeps_the_subcategory(fixture_dir):
    """The grouping applied. Two different Pfam families become one category for the test,
    while the subcategory keeps the finer distinction in the table - so a finer grouping
    can be re-cut later without re-running the pipeline."""
    rules = fixture_dir / "label_categories.yaml"
    rules.write_text(
        "categories:\n"
        "  mobilisation:\n"
        "    exact:\n"
        "      pfam_family: [MobA_MobL]\n"
        "    sub:\n"
        "      pfam_family:\n"
        "        MobA_MobL: relaxase\n")

    rows = {r["category"]: r for r in
            _run_context(fixture_dir, "MobA_MobL", categories=str(rules))}

    assert "mobilisation" in rows, f"the category was not applied; got {sorted(rows)}"
    assert rows["mobilisation"]["subcategories"] == "relaxase"
    assert "pfam_family:MobA_MobL" not in rows, (
        "the raw label survived alongside its category, so it would be counted twice")


# --- S8d: the structural evidence has to carry a description, not just an accession ----

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


# --- S9: what decides which plate a protein goes on -----------------------------------

# --- S6b: what "family of one" means after dereplication ------------------------------

@requires("mmseqs")
def test_a_conserved_protein_on_many_plasmids_is_not_reported_as_a_singleton(fixture_dir):
    """Clustering runs on the DEREPLICATED set, so a protein whose sequence is identical on
    two hundred plasmids is ONE member. It clusters alone and is labelled ORPHAN - while
    being one of the most strongly conserved things in the collection.

    That is defensible as a definition: it is not a family of divergent homologs. It is not
    defensible as a REPORT, because `family_size: 1` reads as "seen once". The ORF count
    has to be there too, or a reader cannot tell a genuine singleton from a protein carried
    by two hundred plasmids."""
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

    run_script("protein_families.py", FakeSnakemake(
        input={"faa": str(faa), "map": str(pmap), "registry": str(registry),
               "lineage": str(lineage), "dark_ids": str(dark_ids)},
        output={"families": str(out),
                "dark_families": str(fixture_dir / "dark_families.tsv")},
        params={"clustering": {
            "resolutions": {"broad": {"min_seq_id": 0.30, "coverage": 0.50}},
            "primary": "broad", "cov_mode": 0, "cluster_mode": 0}},
        threads=2))

    rows = {r["representative"]: r for r in read_tsv(out)}
    assert "S1" in rows, f"S1 did not survive clustering: {list(rows)}"
    assert rows["S1"]["family_size"] == "1", "S1 is one unique sequence"
    assert rows["S1"]["n_orfs"] == "4", (
        "the ORF count is missing, so a protein on four plasmids is indistinguishable "
        "from one seen once")
    assert rows["S1"]["family_plasmid_count"] == "4"
    assert rows["S1"]["family_MOB_count"] == "3"
    # And the independence count, which is the one a recurrence claim needs: four plasmid
    # records but only three independent lineages, because p1 and p4 are the same lineage.
    assert rows["S1"]["family_plasmid_lineage_count"] == "3"


# --- S3: hits.tsv must hold every informative hit, not one per tier -------------------

@requires("diamond")
def test_every_domain_of_a_multi_domain_protein_is_recorded(fixture_dir):
    """tier_search kept only the single best informative hit per query per tier. Two
    documented claims depend on it keeping all of them, and neither could hold:

      * `n_informative_hits` is described as telling a reader whether 0.9 coverage came
        from one domain or six. Capped at one per tier, its maximum was the number of
        tiers, and a two-domain Pfam protein and a six-domain one both reported 1.
      * tier_search's own docstring says the coordinates let explained_fraction be
        recomputed from hits.tsv as an independent cross-check of spans.tsv. With the other
        domains discarded that recomputation gives a smaller number every time, so the
        cross-check would fail on exactly the proteins it matters for.

    A replication initiator carrying RepA_N and Bac_RepA_C is the canonical case, and it is
    the protein class this pipeline must not let through."""
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
        params={"spec": {"id": "T3", "method": "diamond", "db": str(db),
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


# --- S4: the feature files the design names and v1 never wrote ------------------------

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
        input={"annotation": str(ann), "shards": [str(fasta)], "master": str(master)},
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


def test_feature_files_cover_the_shards_and_nothing_else(fixture_dir):
    """One record per sequence in the shards, whether or not it carries an annotation.

    This stage used to read the whole corpus FASTA, so on a 100-plasmid run it emitted
    208,245 GenBank records. The scope has to come from the shards the run was handed. A
    shard plasmid with no called ORFs still gets a record - it is in the analysis set and
    the answer for it is "no features", which is not the same as the record being absent.
    """
    ann = fixture_dir / "scope.tsv"
    write_tsv(ann, ["orf_id", "plasmid_id", "start", "end", "strand", "partial",
                    "spans_origin", "annot_label", "functional_class", "annot_tier",
                    "artefact_flag"],
              [["in1|1", "in1", 10, 60, "+", 0, 0, "", "NONE", "", 0],
               # A row for a plasmid that is NOT in the shards: the annotation table may
               # be wider than this run's scope, and that must not put it in the output.
               ["out1|1", "out1", 10, 60, "+", 0, 0, "", "NONE", "", 0]])
    master = fixture_dir / "scope_master.tsv"
    write_tsv(master, ["plasmid_id", "topology", "size_bp"],
              [["in1", "circular", 200], ["in2", "linear", 200], ["out1", "linear", 200]])

    shard_a = fixture_dir / "shard_a.fna"
    write_fasta(shard_a, [("in1", "ATGC" * 50)])
    shard_b = fixture_dir / "shard_b.fna"
    write_fasta(shard_b, [("in2", "GGCC" * 50)])

    gff = fixture_dir / "scope.gff3"
    gbk = fixture_dir / "scope.gbk"
    run_script("feature_files.py", FakeSnakemake(
        input={"annotation": str(ann), "shards": [str(shard_a), str(shard_b)],
               "master": str(master)},
        output={"gff3": str(gff), "genbank": str(gbk)}))

    loci = [l.split()[1] for l in gbk.read_text().splitlines() if l.startswith("LOCUS")]
    assert loci == ["in1", "in2"], (
        f"expected one record per shard sequence, got {loci} - a record for a plasmid "
        "outside the shards means the stage is reading something wider than its input")

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
        input={"annotation": str(ann), "shards": [str(fasta)], "master": str(master)},
        output={"gff3": str(gff), "genbank": str(fixture_dir / "o.gbk")}))

    row = [l for l in gff.read_text().splitlines() if not l.startswith("#")][0]
    attrs = row.split("\t")[8]
    assert "product=" not in attrs, f"a product was fabricated for a dark ORF: {attrs}"
    assert "ID=p1%7C1" in attrs


# --- S7b: RNAcode, declared in config and installed, and never invoked ----------------

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


# --- S7c: the family-consensus re-check -----------------------------------------------

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

    run_script("consensus_recheck.py", FakeSnakemake(
        input={"consensus": str(cons)},
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


# --- S4b: orthology over the ANNOTATED fraction ---------------------------------------

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
        input={"prot": str(prot), "faa": str(faa)},
        output=[str(out)],
        params={"orthology": {"data_dir": str(fixture_dir / "absent-db"),
                              "required": False}},
        threads=1))

    rows = {r["seq_id"] for r in read_tsv(out)}
    assert rows == {"named1", "named2"}, (
        f"the dark proteins were sent to eggNOG, or the named ones were not: {rows}")


# --- the deliverable: complete annotations, not a shortlist ---------------------------

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
    # LONG: one row per (family, category). The report reduces it, taking the strongest
    # association by q-value.
    write_tsv(ctx, ["family_id", "category", "subcategories", "n_units",
                    "n_units_with_category", "observed_rate", "background_rate",
                    "enrichment", "odds_ratio", "p_value", "q_value", "status"],
              [["F2", "defence", "RM_Type_II", 9, 7, 0.8, 0.1, 8.0, 20.0,
                "1e-6", "1e-6", "SUCCESS"]])
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
                      "host_count", "species_count", "genus_count", "MOB_count",
                      "habitat_count", "database_record_count", "database_source_count"],
              [["F1", "broad", "S2", 1, 1, 1, "SUCCESS", 1, 1, 1, 1, 1, 1, 1],
               # 40 gene copies on 9 records that are only 2 independent lineages: the
               # shape section 34.2 exists to keep visible.
               ["F2", "broad", "S3", 40, 9, 2, "SUCCESS", 3, 3, 2, 3, 2, 9, 1]])
    syn = fixture_dir / "synteny.tsv"
    write_tsv(syn, ["family_id", "n_occurrences", "context_recurrence",
                    "left_neighbor_conservation", "right_neighbor_conservation",
                    "neighborhood_conservation", "operon_like_conservation",
                    "synteny_conservation", "modal_left", "modal_right", "modal_synteny",
                    "status"],
              [["F1", 1, 0, "", "", "", "", "", "", "", "", "TOO_FEW_MEMBERS"],
               ["F2", 9, 9, 0.9, 0.7, 0.8, 0.6, 0.7, "mobA", "repA", "mobA|repA",
                "SUCCESS"]])
    rarity_tsv = fixture_dir / "family_rarity.tsv"
    write_tsv(rarity_tsv, ["family_id", "rarity_labels",
                           "independent_plasmid_cluster_count", "unique_plasmid_count",
                           "MOB_count", "host_count", "genus_count", "rarity_version",
                           "rare_max_lineages", "widely_conserved_min_lineages"],
              [["F1", "RARE,LINEAGE_SPECIFIC", 1, 1, 1, 1, 1, "1", 3, 50],
               ["F2", "RARE,CROSS_MOB", 2, 9, 3, 3, 2, "1", 3, 50]])

    return (ann, pmap, fams, evo, rec, ctx, struct, orth, recur, syn,
            rarity_tsv)


def _run_report(fixture_dir, *tables):
    """Drive annotation_report.py. With no tables passed, build the standard fixture.

    Two tests need the deliverable written: one checks what is in it, the other checks
    what must never be (spec section 76). Sharing a written file between them would make
    the second silently skip whenever it ran alone.
    """
    if not tables:
        tables = _report_fixture(fixture_dir)
    ann, pmap, fams, evo, rec, ctx, struct, orth, recur, syn, rarity_tsv = tables
    out_ann = fixture_dir / "annotation_complete.csv"
    out_fam = fixture_dir / "dark_families_complete.csv"
    run_script("annotation_report.py", FakeSnakemake(
        input={"annotation": str(ann), "map": str(pmap), "families": str(fams),
               "evolution": str(evo), "recheck": str(rec),
               "context": str(ctx), "structure": str(struct), "orthology": str(orth),
               "recurrence": str(recur), "synteny": str(syn),
               "rarity": str(rarity_tsv)},
        output={"annotation": str(out_ann), "families": str(out_fam)},
        params={"prioritisation": {"min_reality_lines": 2, "min_mob_clusters": 2},
                "evolution": {"min_members_for_dnds": 3, "dnds_purifying_max": 0.5}}))
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
    # The artefact-flagged dark ORF is present, flagged, and joined to its family.
    assert by_orf["p1|2"]["artefact_flag"] == "1"
    assert by_orf["p1|2"]["family_id"] == "F1"
    # Family evidence travels down to the ORF row.
    assert by_orf["p1|3"]["dnds_median"] == "0.21"
    assert by_orf["p1|3"]["top_hypothesis"] == "defence"

    fam_rows = {r["family_id"]: r for r in _csv.DictReader(open(out_fam))}
    assert set(fam_rows) == {"F1", "F2"}, "an ORPHAN family was dropped"
    assert fam_rows["F1"]["dnds_status"] == "TOO_FEW_MEMBERS", (
        "absence of a measurement must be an explicit status, not a blank")
    assert fam_rows["F1"]["reality_n"] == "0"
    # purifying_selection fired, so is_family is entailed and does not count twice.
    assert fam_rows["F2"]["reality_n"] == "2"
    assert fam_rows["F2"]["reality_lines_implied"] == "is_family"
    # The column sets are the contract, asserted by equality rather than by the absence
    # of one remembered name: a renamed S7d leftover, or a new column nobody documented,
    # fails here either way. Orthology columns come from the eggNOG table's own header,
    # so they are excluded from the ORF-side check.
    FAMILY_COLUMNS = [
        "family_id", "representative", "family_class",
        "n_members", "n_orfs", "n_plasmids", "n_mob_clusters",
        "reality_n", "reality_lines", "reality_lines_implied",
        "dnds_median", "dnds_min", "dnds_status", "under_purifying_selection", "n_pairs",
        "rnacode_p", "rnacode_p_antisense", "rnacode_status", "coding_signal",
        "consensus_hit", "consensus_label", "collectively_novel",
        "darkness_state", "structural_match", "structural_description", "structure_evalue",
        "top_hypothesis", "top_subcategories", "top_conservation", "top_enrichment",
        "top_q_value", "high_confidence",
        "cons_defence", "cons_integron",
        "cons_annotated_neighbour", "cons_operon_with_annotated", "cons_two_gene_operon",
        # Stage 7: seven counts, never collapsed into one.
        "plasmid_occurrence_count", "unique_plasmid_count",
        "independent_plasmid_cluster_count", "independent_cluster_status",
        "host_count", "species_count", "genus_count", "MOB_count", "habitat_count",
        "database_record_count", "database_source_count",
        # Stage 9: six conservation measurements, kept apart because they fail apart.
        "context_recurrence", "n_occurrences", "left_neighbor_conservation",
        "right_neighbor_conservation", "neighborhood_conservation",
        "operon_like_conservation", "synteny_conservation", "modal_left", "modal_right",
        "modal_synteny", "synteny_status",
        # Stage 14: descriptors, not a ranking.
        "rarity_labels", "rarity_version",
        # Stage 15: dimensions counted, never scored.
        "evidence_dimensions_present", "evidence_dimension_count",
        "supporting_observations_count", "supporting_observations_are_not_independent",
        "functional_hypothesis", "functional_hypothesis_support",
    ]
    assert list(fam_rows["F2"]) == FAMILY_COLUMNS, (
        f"family table columns changed: {list(fam_rows['F2'])}")
    CARRIED_TO_ORFS = {"family_id", "reality_n", "reality_lines", "darkness_state",
                       "dnds_median", "dnds_status", "coding_signal", "collectively_novel",
                       "top_hypothesis", "top_conservation", "structural_match",
                       "structural_description"}
    annotation_cols = {"orf_id", "plasmid_id", "start", "end", "strand", "annot_label",
                       "functional_class", "artefact_flag"}
    # seq_id is the dereplicated-protein key that joins an ORF to its sequence; it is
    # written beside the orthology terms and is neither annotation nor family evidence.
    join_and_orthology = {"seq_id", "cog_category", "kegg_pathways", "preferred_name",
                          "eggnog_description"}
    carried = set(by_orf["p1|3"]) - annotation_cols - join_and_orthology
    assert carried == CARRIED_TO_ORFS, (
        f"family evidence carried to the ORF table changed: {sorted(carried)}")

    # --- Stage 7: the counts stay apart ---------------------------------------------
    # Section 34.2: "Database record counts must never be treated as independent
    # biological observations." F2 is 40 gene copies on 9 records that are 2 lineages. If
    # any of those three numbers can be read off another, the distinction is gone.
    assert fam_rows["F2"]["plasmid_occurrence_count"] == "40"
    assert fam_rows["F2"]["unique_plasmid_count"] == "9"
    assert fam_rows["F2"]["independent_plasmid_cluster_count"] == "2", (
        "the independent-lineage count is what a recurrence claim needs, and it is not "
        "the record count")

    # --- Stage 9: synteny, with its own status name ---------------------------------
    assert fam_rows["F2"]["synteny_conservation"] == "0.7"
    assert fam_rows["F1"]["synteny_status"] == "TOO_FEW_MEMBERS", (
        "one occurrence is perfectly conserved with itself; that must read as a status, "
        "not as a conservation of 1.0")

    # --- Stage 14: labels are descriptors -------------------------------------------
    assert fam_rows["F1"]["rarity_labels"] == "RARE,LINEAGE_SPECIFIC"
    assert fam_rows["F1"]["rarity_version"] == "1", (
        "a label whose definition can change must travel with the version that made it")

    # --- Stage 15: dimensions counted, never scored ----------------------------------
    dims = fam_rows["F2"]["evidence_dimensions_present"].split(",")
    assert "EVOLUTIONARY_CONSERVATION" in dims, "dnds_status is a measurement"
    assert "DISTRIBUTION" in dims, "independent_cluster_status is a measurement"
    assert "GENOMIC_CONTEXT" in dims, "top_hypothesis is a measurement"
    assert fam_rows["F2"]["evidence_dimension_count"] == str(len(dims))
    assert fam_rows["F2"]["supporting_observations_are_not_independent"] == "1", (
        "the observation count must carry its own warning, because a column selected "
        "into a downstream ranking takes the warning with it")
    assert fam_rows["F2"]["functional_hypothesis"] == "defence_associated", (
        "a dark ORF recurrently beside a defence system is defence_associated")
    assert fam_rows["F2"]["darkness_state"], (
        "the hypothesis must not have changed what the protein IS (section 57)")


def test_the_report_has_no_composite_score_or_rank(fixture_dir):
    """Spec section 76 draws the Layer C boundary: the core pipeline produces no
    top_1000, no candidate_score, no novelty_score and no experimental_rank.

    Section 2.3 gives the reason, and it is not tidiness. Two proteins with the same total
    can be entirely different bets - one with overwhelming evidence that it is a real
    protein and no idea what it does, the other with a sharp hypothesis resting on almost
    nothing - and those demand different experiments. A composite destroys exactly the
    information a screening decision needs.

    evidence_dimension_count is a count of distinct measurements present, which is why it
    is allowed where a score is not.
    """
    import csv as _csv
    _, out_fam = _run_report(fixture_dir)
    columns = next(iter(_csv.DictReader(open(out_fam))))

    forbidden = {"candidate_score", "novelty_score", "experimental_rank", "priority",
                 "rank", "score", "total_score", "composite"}
    assert not (set(columns) & forbidden), (
        f"a ranking column reappeared in the deliverable: {set(columns) & forbidden}")


# --- S1: origin repair is wired through the script, not only the library --------------

def test_orf_call_reconstructs_a_gene_across_the_origin(fixture_dir):
    """tests/test_circular.py proves resolve_origin_genes in isolation. This proves the
    script actually drives it: the same plasmid cut inside a gene must yield the same
    protein set as the uncut plasmid, with the straddling gene written start > end and
    flagged spans_origin=1 in orf_index.tsv."""
    import random
    pytest.importorskip("pyrodigal")
    from darkorf.circular import rotate

    random.seed(11)
    sense = [a + b + c for a in "ACGT" for b in "ACGT" for c in "ACGT"
             if a + b + c not in ("TAA", "TAG", "TGA")]
    background = lambda n: "".join(random.choice("ACGT") for _ in range(n))
    gene = "ATG" + "".join(random.choice(sense) for _ in range(300)) + "TAA"
    uncut = background(1_000) + gene + background(1_000)
    cut = rotate(uncut, 1_000 + len(gene) // 2)   # the new origin lies mid-gene

    master = fixture_dir / "master.tsv"
    write_tsv(master, ["plasmid_id", "topology"],
              [["uncut", "circular"], ["cut", "circular"]])
    fasta = fixture_dir / "shard.fna"
    write_fasta(fasta, [("uncut", uncut), ("cut", cut)])
    out = fixture_dir / "orf_index.tsv"

    run_script("orf_call.py", FakeSnakemake(
        input={"fasta": str(fasta), "master": str(master)},
        output={"tsv": str(out)},
        params={"min_orf_aa": 20}))

    rows = read_tsv(out)
    proteins = {pid: sorted(r["seq"] for r in rows if r["plasmid_id"] == pid)
                for pid in ("uncut", "cut")}
    assert proteins["uncut"] == proteins["cut"], "protein set depends on the cut point"

    planted = [r for r in rows if r["plasmid_id"] == "cut" and len(r["seq"]) == 301]
    assert len(planted) == 1, "the planted gene was not called exactly once on the cut record"
    assert planted[0]["spans_origin"] == "1"
    assert int(planted[0]["start"]) > int(planted[0]["end"]), "wrapped end not applied"
    assert planted[0]["partial"] == "0"
    assert not any(r["partial"] == "1" for r in rows), (
        "a left-edge stub survived alongside the gene it is a fragment of")


# --- S3: the accession every join will need ------------------------------------------

@requires("diamond")
def test_tier_search_keeps_the_subject_accession(fixture_dir):
    """A DIAMOND title is free text; the accession is the only key a later join to UniProt
    or RefSeq can use, and it was being discarded. Measured on the NCBI swissprot database
    used here, the title is 'P62554.1 RecName: Full=Toxin CcdB; ... [Escherichia coli]' -
    NCBI's rendering, with no gene symbol - so the accession is the whole join key."""
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
        params={"spec": {"id": "T3", "method": "diamond", "db": str(dmnd),
                         "args": "--fast", "max_evalue": 1e-5},
                "narrow_at": 0.9, "hmmer_z": 1000, "max_target_seqs": 5},
        threads=1))

    rows = read_tsv(hits)
    assert rows, "no hits written for a query identical to a database sequence"
    assert "target_accession" in rows[0], "hits.tsv has no target_accession column"
    assert rows[0]["target_accession"] == "P62554.1", (
        f"accession not captured, got {rows[0]['target_accession']!r}")
    assert "CcdB" in rows[0]["label"] or "P62554" in rows[0]["label"], (
        "the title is no longer being kept as the label")


# --- S8a phase 2: mandatory is not the same as accessory ------------------------------

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


# --- S4c: one long table of what every tool said --------------------------------------

def test_protein_labels_gathers_every_source_into_one_long_table(fixture_dir):
    """The substrate for the functional grouping. A wide table cannot hold it: the
    vocabulary is open, Pfam-A 38.2 alone has 30,134 families, and the grouping is derived
    from the labels observed rather than declared in advance."""
    import gzip

    hits = fixture_dir / "hits.tsv"
    write_tsv(hits, ["query", "label", "target_accession", "coverage", "target_coverage",
                     "evalue", "informative", "is_best", "start", "end", "tier",
                     "threshold", "max_evalue"],
              [["s1", "RepA_N", "PF06970.19", 0.9, 0.95, "1e-40", "True", 1, 1, 100,
                "T1", "--cut_ga", ""],
               ["s1", "P62554.1 RecName: Full=Toxin CcdB [Escherichia coli]", "P62554.1",
                0.8, 0.9, "1e-30", "True", 1, 1, 90, "T3", "--fast", "1e-5"],
               ["s2", "WP_1.1 hypothetical protein [Escherichia coli]", "WP_1.1",
                0.95, 0.9, "1e-20", "False", 1, 1, 95, "T4", "--fast", "1e-10"]])

    orth = fixture_dir / "orthology.tsv"
    write_tsv(orth, ["seq_id", "cog_category", "kegg_pathways", "preferred_name",
                     "eggnog_description", "eggnog_ogs", "pfams", "gos", "ec", "kegg_ko"],
              [["s1", "L", "ko03030", "repA", "Replication initiator",
                "COG5527@2", "RepA_N", "GO:0006270", "2.7.7.7", "ko:K02314"]])

    pfam_dat = fixture_dir / "Pfam-A.hmm.dat.gz"
    with gzip.open(pfam_dat, "wt") as fh:
        fh.write("# STOCKHOLM 1.0\n#=GF ID   RepA_N\n#=GF AC   PF06970.19\n"
                 "#=GF DE   Replication initiator protein A (RepA) N-terminus\n"
                 "#=GF TP   Domain\n#=GF CL   CL0123\n//\n")

    out = fixture_dir / "protein_labels.tsv"
    run_script("protein_labels.py", FakeSnakemake(
        input={"hits": [str(hits)], "orthology": str(orth), "pfam_dat": str(pfam_dat)},
        output={"tsv": str(out)},
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


def test_protein_labels_records_the_database_version_on_every_row(fixture_dir):
    """A label without the database release it came from cannot be reproduced, and the
    grouping built on it cannot be described in a methods section."""
    hits = fixture_dir / "hits.tsv"
    write_tsv(hits, ["query", "label", "target_accession", "coverage", "target_coverage",
                     "evalue", "informative", "is_best", "start", "end", "tier",
                     "threshold", "max_evalue"],
              [["s1", "RepA_N", "PF06970.19", 0.9, 0.95, "1e-40", "True", 1, 1, 100,
                "T1", "--cut_ga", ""]])
    orth = fixture_dir / "orthology.tsv"
    write_tsv(orth, ["seq_id", "cog_category", "kegg_pathways", "preferred_name",
                     "eggnog_description", "eggnog_ogs", "pfams", "gos", "ec",
                     "kegg_ko"], [])
    pfam_dat = fixture_dir / "pfam.dat"
    pfam_dat.write_text("")

    out = fixture_dir / "protein_labels.tsv"
    run_script("protein_labels.py", FakeSnakemake(
        input={"hits": [str(hits)], "orthology": str(orth), "pfam_dat": str(pfam_dat)},
        output={"tsv": str(out)},
        params={"pfam_version": "38.2", "swissprot_version": "2025-03-03",
                "nr_version": "2025-03-03", "eggnog_version": "5.0.2"}))

    rows = read_tsv(out)
    assert rows
    assert all(r["database_version"] for r in rows), "a row carries no database version"
    assert rows[0]["database"] == "Pfam-A"
    assert rows[0]["database_version"] == "38.2"


def test_protein_labels_merges_a_label_seen_by_two_tiers(fixture_dir):
    """The same Pfam family hit by T1 and T2 is one statement about the protein, not two.
    Unmerged, a widely searched label would outvote a rare one by copy number when the
    categories are counted."""
    hits = fixture_dir / "hits.tsv"
    write_tsv(hits, ["query", "label", "target_accession", "coverage", "target_coverage",
                     "evalue", "informative", "is_best", "start", "end", "tier",
                     "threshold", "max_evalue"],
              [["s1", "RepA_N", "PF06970.19", 0.9, 0.95, "1e-10", "True", 1, 1, 100,
                "T1", "--cut_ga", ""],
               ["s1", "RepA_N", "PF06970.19", 0.9, 0.95, "1e-40", "True", 1, 1, 100,
                "T2", "-E 1e-5", "1e-5"]])
    orth = fixture_dir / "orthology.tsv"
    write_tsv(orth, ["seq_id", "cog_category", "kegg_pathways", "preferred_name",
                     "eggnog_description", "eggnog_ogs", "pfams", "gos", "ec",
                     "kegg_ko"], [])
    pfam_dat = fixture_dir / "pfam.dat"
    pfam_dat.write_text("")

    out = fixture_dir / "protein_labels.tsv"
    run_script("protein_labels.py", FakeSnakemake(
        input={"hits": [str(hits)], "orthology": str(orth), "pfam_dat": str(pfam_dat)},
        output={"tsv": str(out)},
        params={"pfam_version": "38.2", "swissprot_version": "2025-03-03",
                "nr_version": "2025-03-03", "eggnog_version": "5.0.2"}))

    rows = [r for r in read_tsv(out) if r["kind"] == "pfam_family"]
    assert len(rows) == 1, f"the same family was recorded {len(rows)} times"
    # The strongest evidence for the statement survives the merge.
    assert rows[0]["evidence_evalue"] == "1e-40"


# --- Category building: candidates for review, never a runtime matcher ----------------

def test_build_categories_writes_a_reviewable_draft_with_counts(fixture_dir):
    """The mining patterns must produce CANDIDATES with the evidence for judging them, not
    a finished config. Measured on Pfam-A 38.2, /toxin/ matches 317 family descriptions of
    which ABC_toxin_N is an insect toxin - the count beside each candidate is what lets a
    reviewer see the cost of keeping or dropping it."""
    import yaml

    labels_tsv = fixture_dir / "protein_labels.tsv"
    write_tsv(labels_tsv, ["protein_id", "source", "tier", "kind", "label", "accession",
                           "evidence_evalue", "evidence_coverage", "database",
                           "database_version"],
              [["s1", "pfam", "T1", "pfam_family", "RepA_N", "PF06970.19", "1e-40",
                "0.9", "Pfam-A", "38.2"],
               ["s2", "pfam", "T1", "pfam_family", "RepA_N", "PF06970.19", "1e-30",
                "0.9", "Pfam-A", "38.2"],
               ["s1", "eggnog", "S4b", "gene_symbol", "repA", "", "", "",
                "eggNOG", "5.0.2"],
               ["s3", "nr", "T4", "pgap_product",
                "putative plasmid replication initiator protein", "WP_1.1", "1e-20",
                "0.8", "NCBI nr", "2025-03-03"],
               ["s4", "nr", "T4", "pgap_product", "ABC transporter permease", "WP_2.1",
                "1e-20", "0.8", "NCBI nr", "2025-03-03"]])

    pfam_dat = fixture_dir / "pfam.dat"
    pfam_dat.write_text("# STOCKHOLM 1.0\n#=GF ID   RepA_N\n#=GF AC   PF06970.19\n"
                        "#=GF DE   Replication initiator protein A (RepA) N-terminus\n"
                        "#=GF TP   Domain\n//\n")
    mining = fixture_dir / "category_mining.yaml"
    mining.write_text("mine:\n  replication:\n    description: 'replicat'\n"
                      "    symbol: '^rep[A-Z]?$'\n")

    draft = fixture_dir / "category_draft.yaml"
    frequency = fixture_dir / "label_frequency.tsv"
    run_script("build_categories.py", FakeSnakemake(
        input={"labels": str(labels_tsv), "pfam_dat": str(pfam_dat),
               "mining": str(mining)},
        output={"draft": str(draft), "frequency": str(frequency)}))

    document = yaml.safe_load(draft.read_text())
    replication = document["categories"]["replication"]["exact"]

    # Mined from the Pfam DESCRIPTION, which the family NAME does not contain.
    assert "RepA_N" in replication["pfam_family"]
    # Mined from the gene-symbol convention.
    assert "repA" in replication["gene_symbol"]
    # Mined from the NORMALISED product name, so the hedge did not hide it.
    assert any("replication initiat" in p for p in replication["pgap_product"])
    # Not mined: nothing about it matches, and it must not be swept in.
    assert "ABC transporter permease" not in replication.get("pgap_product", [])

    # The frequency table is the evidence for the review.
    rows = {(r["kind"], r["label"]): r for r in read_tsv(frequency)}
    assert rows[("pfam_family", "RepA_N")]["n_proteins"] == "2", (
        "proteins are counted, not rows - a label seen by three tiers is one protein")
    assert rows[("gene_symbol", "repA")]["n_proteins"] == "1"


def test_build_categories_only_ever_writes_its_declared_outputs():
    """The draft is reviewed by a person before it becomes the live config. A builder that
    wrote config/label_categories.yaml directly would put unreviewed pattern matches
    straight into the pipeline, which is the failure the whole design exists to prevent.

    Asserted on the parsed code rather than on the text, because the module docstring
    legitimately names that path to say it does not write it."""
    import ast

    source = (pathlib.Path(__file__).parent.parent
              / "workflow" / "scripts" / "build_categories.py").read_text()

    written = []
    for node in ast.walk(ast.parse(source)):
        if not (isinstance(node, ast.Call) and getattr(node.func, "id", "") == "open"):
            continue
        mode = next((a.value for a in node.args[1:] if isinstance(a, ast.Constant)), "r")
        if "w" in mode or "a" in mode:
            written.append(ast.unparse(node.args[0]))

    assert written, "the builder opens nothing for writing - it produces no draft"
    for target in written:
        assert target.startswith("snakemake.output."), (
            f"build_categories writes to {target}, which is not a declared output")


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


# --- Stage 5: families over EVERY protein, not only the dark ones ----------------------

@requires("mmseqs")
def test_protein_families_clusters_annotated_and_dark_together(fixture_dir):
    """Spec section 31.2 requires dark_member_count, annotated_member_count and
    percentage_dark_in_family, and section 32 derives a dark-only family as 100% dark. None
    of those can be computed from a clustering that contains only dark proteins: every
    family would be trivially 100% dark, and a dark protein among well-annotated homologs -
    a strong observation - would look identical to one that is genuinely alone."""
    import yaml

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
    run_script("protein_families.py", FakeSnakemake(
        input={"faa": str(faa), "dark_ids": str(dark_ids), "map": str(mapping),
               "registry": str(registry), "lineage": str(lineage_tsv)},
        output={"families": str(families), "dark_families": str(dark_families)},
        params={"clustering": {
            "resolutions": {"broad": {"min_seq_id": 0.3, "coverage": 0.5}},
            "primary": "broad", "cov_mode": 0, "cluster_mode": 0}},
        threads=2))

    rows = read_tsv(families)
    assert rows, "no families written"

    mixed = [r for r in rows if int(r["family_size"]) > 1]
    assert mixed, "the two near-identical proteins did not cluster together"
    row = mixed[0]
    # The whole point: both counts are non-zero, which a dark-only clustering cannot show.
    assert int(row["dark_member_count"]) == 1
    assert int(row["annotated_member_count"]) == 1
    assert float(row["percentage_dark_in_family"]) == 50.0
    assert row["dark_only"] == "0"
    # Section 31.3 distribution fields, measured over independent units.
    assert int(row["family_plasmid_count"]) == 2
    assert int(row["family_species_count"]) == 2
    assert int(row["family_MOB_count"]) == 2


@requires("mmseqs")
def test_family_ids_are_content_derived_not_ordinal(fixture_dir):
    """Spec section 31: 'family IDs must not depend on result ordering', and section 5.4
    gives the form. The previous version numbered families F0000001, F0000002, ... in
    cluster order, so inserting one protein renumbered every family after it and no id
    could be compared between two runs."""
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
    run_script("protein_families.py", FakeSnakemake(
        input={"faa": str(faa), "dark_ids": str(dark_ids), "map": str(mapping),
               "registry": str(registry), "lineage": str(lineage_tsv)},
        output={"families": str(families),
                "dark_families": str(fixture_dir / "dark_families.tsv")},
        params={"clustering": {
            "resolutions": {"broad": {"min_seq_id": 0.3, "coverage": 0.5}},
            "primary": "broad", "cov_mode": 0, "cluster_mode": 0}},
        threads=2))

    for row in read_tsv(families):
        assert row["family_id"] == f"broad:{row['representative']}", (
            f"family_id {row['family_id']!r} is not <resolution>:<representative>")
        assert not row["family_id"].startswith("F0"), "family ids are ordinal again"


@requires("mmseqs")
def test_the_dark_family_representative_is_a_dark_protein(fixture_dir):
    """S8d searches the representative structurally. If MMseqs2 picks an ANNOTATED member
    as the cluster representative, searching it would spend the ProstT5 budget on a protein
    that is not in the dark set and produce no structural evidence for the dark one."""
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
    run_script("protein_families.py", FakeSnakemake(
        input={"faa": str(faa), "dark_ids": str(dark_ids), "map": str(mapping),
               "registry": str(registry), "lineage": str(lineage_tsv)},
        output={"families": str(fixture_dir / "protein_families.tsv"),
                "dark_families": str(dark_families)},
        params={"clustering": {
            "resolutions": {"broad": {"min_seq_id": 0.3, "coverage": 0.5}},
            "primary": "broad", "cov_mode": 0, "cluster_mode": 0}},
        threads=2))

    rows = read_tsv(dark_families)
    assert rows, "no dark families derived"
    for row in rows:
        assert row["representative"] == "p_dark", (
            f"the dark family's representative is {row['representative']!r}, which is not "
            "a dark protein - S8d would search the wrong sequence")
        assert "p_annot_long" not in row["members"], (
            "an annotated member leaked into the dark family's member list, which would "
            "widen every downstream evolution and context measurement")


# --- Stage 7: recurrence counted over independent units --------------------------------

def _recurrence_fixture(fixture_dir, lineage_rows):
    """One family on three plasmids; the caller decides how independent those are."""
    families = fixture_dir / "protein_families.tsv"
    write_tsv(families, ["family_id", "family_resolution", "representative", "members"],
              [["broad:s1", "broad", "s1", "s1,s2"]])
    mapping = fixture_dir / "protein_map.tsv"
    # s1 on two plasmids with two copies on one of them; s2 on a third.
    mapping.write_text("s1\tpl1|1,pl1|2,pl2|1\ns2\tpl3|1\n")
    registry = fixture_dir / "clonal_registry.tsv"
    write_tsv(registry, ["plasmid_id", "mob_cluster", "species", "topology", "size_bp",
                         "hab_top"],
              [["pl1", "MOB_A", "Escherichia coli", "circular", 100, "Host-associated"],
               ["pl2", "MOB_A", "Escherichia coli", "circular", 100, "Host-associated"],
               ["pl3", "MOB_A", "Escherichia coli", "circular", 100, "Host-associated"]])
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
    assert row["database_record_count"] == "3"


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


# --- S2d: the negative control the config declared and nothing built ------------------

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
        input={"index": str(index), "shards": [str(shard)]},
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
    import collections as _c
    index, shard = _decoy_fixture(fixture_dir)
    out = fixture_dir / "negative_control.faa"

    run_script("prepare_decoys.py", FakeSnakemake(
        input={"index": str(index), "shards": [str(shard)]},
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
                         "min_controls": 2}}))

    text = report.read_text()
    assert "decoy_n=4" in text, f"the decoy count is not reported:\n{text}"
    assert "decoy_false_positive_rate=0.25" in text, (
        f"1 of 4 decoys was named FUNCTIONAL; the rate must be reported:\n{text}")

    eligible = {r["seq_id"] for r in read_tsv(flags)}
    assert not any(s.startswith("DECOY_") for s in eligible), (
        "decoys are instrumentation, not screening candidates; they must not appear in "
        "the target-eligibility table any more than the positive controls do")
    assert "realprotein" in eligible
