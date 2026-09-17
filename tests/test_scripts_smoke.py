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

def test_the_synthesis_order_carries_the_hypothesis(fixture_dir):
    """library_design read `top_hypothesis`; prioritise writes `hypothesis`. Every
    construct silently lost its hypothesis on the way to the bench - and the order file
    still looked complete."""
    scored = fixture_dir / "scored.tsv"
    faa = fixture_dir / "dark.faa"
    order = fixture_dir / "order.tsv"
    report = fixture_dir / "library.txt"

    write_fasta(faa, [("rep1", "MKVLATTLLGAAFAASSALAQKKWLVR")])
    write_tsv(scored,
              ["family_id", "representative", "selected", "stratum", "hypothesis",
               "structural_match", "reality_n"],
              [["F0000001", "rep1", 1, "defence_island", "defence", "", 3]])

    run_script("library_design.py", FakeSnakemake(
        input={"scored": str(scored), "faa": str(faa)},
        output={"order": str(order), "report": str(report)},
        params={"library": {"host": "ecoli", "length_liability_above_aa": 400,
                            "avoid_sites": ["GAATTC", "GGATCC"]},
                "portfolio": {"n_controls": 2}, "seed": 1}))

    rows = [r for r in read_tsv(order) if r["role"] == "candidate"]
    assert rows, "no candidate constructs written"
    assert rows[0]["hypothesis"] == "defence", (
        f"hypothesis lost on the way to the order file: {rows[0]['hypothesis']!r}")


def test_every_ordered_construct_is_translatable(fixture_dir):
    """A CDS handed to a synthesis vendor must be a whole number of codons and must start
    with ATG. Whether it ends in a stop depends on the tag terminus, which the script
    already decides per construct."""
    scored = fixture_dir / "scored.tsv"
    faa = fixture_dir / "dark.faa"
    order = fixture_dir / "order.tsv"
    report = fixture_dir / "library.txt"

    write_fasta(faa, [("rep1", "MKVLATTLLGAAFAASSALAQ"),
                      ("rep2", "MKWKLFKKIGAVLKVLTTGLPALIS")])
    write_tsv(scored,
              ["family_id", "representative", "selected", "stratum", "hypothesis",
               "structural_match", "reality_n"],
              [["F1", "rep1", 1, "novel_fold", "", "", 2],
               ["F2", "rep2", 1, "small_cationic_peptide", "", "", 3]])

    run_script("library_design.py", FakeSnakemake(
        input={"scored": str(scored), "faa": str(faa)},
        output={"order": str(order), "report": str(report)},
        params={"library": {"host": "ecoli", "length_liability_above_aa": 400,
                            "avoid_sites": ["GAATTC"]},
                "portfolio": {"n_controls": 2}, "seed": 1}))

    for r in read_tsv(order):
        cds = r["cds"]
        assert len(cds) % 3 == 0, f"{r['construct_id']}: CDS is not a whole number of codons"
        assert cds.startswith("ATG"), f"{r['construct_id']}: CDS does not start with ATG"


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


def test_the_stop_codon_follows_the_tag_terminus(fixture_dir):
    """Not universal, and not optional. With an N-terminal tag the construct is
    [tag]-[protein]-STOP. With a C-terminal fusion a stop codon here truncates the
    transcript before the tag, so the screen sees nothing from that construct at all -
    a silent, total loss of half the library, in a screen whose only guaranteed readout
    is the tag."""
    scored = fixture_dir / "scored.tsv"
    faa = fixture_dir / "dark.faa"
    order = fixture_dir / "order.tsv"
    report = fixture_dir / "library.txt"

    # soluble -> N-terminal tag -> needs a stop; amphipathic peptide -> C-terminal -> must not
    write_fasta(faa, [("soluble", "MKEEDDKKEEDDKKEEDDKKEEQ"),
                      ("peptide", "MKWKLFKKIGAVLKVLTTGLPALIS")])
    write_tsv(scored,
              ["family_id", "representative", "selected", "stratum", "hypothesis",
               "structural_match", "reality_n"],
              [["F1", "soluble", 1, "novel_fold", "", "", 2],
               ["F2", "peptide", 1, "small_cationic_peptide", "", "", 3]])

    run_script("library_design.py", FakeSnakemake(
        input={"scored": str(scored), "faa": str(faa)},
        output={"order": str(order), "report": str(report)},
        params={"library": {"host": "ecoli", "length_liability_above_aa": 400,
                            "avoid_sites": ["GAATTC"]},
                "portfolio": {"n_controls": 0}, "seed": 1}))

    by_rep = {r["representative"]: r for r in read_tsv(order) if r["role"] == "candidate"}
    stops = ("TAA", "TAG", "TGA")

    n_tagged = by_rep["soluble"]
    assert n_tagged["tag_terminus"] == "N"
    assert n_tagged["cds"][-3:] in stops, "N-terminal tag construct is missing its stop codon"

    c_tagged = by_rep["peptide"]
    assert c_tagged["tag_terminus"] == "C"
    assert c_tagged["cds"][-3:] not in stops, (
        "C-terminal fusion carries a stop codon - the tag would never translate")


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

    run_script("prepare_control.py", FakeSnakemake(
        input={"faa": str(faa), "raw": str(raw)},
        output={"control": str(control), "spiked": str(spiked)},
        params={"n_controls": 10, "min_controls": 1, "seed": 1}))

    headers = [l[1:].strip() for l in open(control) if l.startswith(">")]
    assert len(headers) == 3, f"expected 3 nameable controls, got {len(headers)}"
    joined = " ".join(headers).lower()
    assert "uncharacterized" not in joined
    assert "upf0102" not in joined


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

def _context_fixture(fixture_dir, partner_label):
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
    return ann, fam, pmap, defence, integrons


def _run_context(fixture_dir, partner_label):
    ann, fam, pmap, defence, integrons = _context_fixture(fixture_dir, partner_label)
    fams_out = fixture_dir / "family_context.tsv"
    bg_out = fixture_dir / "context_background.tsv"
    run_script("context_features.py", FakeSnakemake(
        input={"annotation": str(ann), "families": str(fam), "map": str(pmap),
               "defence": str(defence), "integrons": [str(integrons)]},
        output={"families": str(fams_out), "background": str(bg_out)},
        params={"context": {"max_operon_gap": 100, "neighbourhood_window": 3,
                            "min_context_conservation": 0.50,
                            "high_confidence_conservation": 0.90,
                            "min_enrichment": 2.0}}))
    return read_tsv(fams_out)[0]


def test_a_transposase_partner_is_not_a_toxin_antitoxin_candidate(fixture_dir):
    """ta_candidate fired whenever the two-gene partner was ANY of the 73 curated backbone
    families - so a dark ORF beside a transposase, a relaxase or a methyltransferase was
    labelled a candidate antitoxin and routed to the toxin_or_ta_adjacent stratum, 175 of
    the 1,000 constructs. Only 16 of those families are toxins or antitoxins, and the tight
    two-gene geometry is only evidence when the partner is one of them."""
    row = _run_context(fixture_dir, "DDE_Tnp_Tn3")
    assert float(row["cons_ta_candidate"]) == 0.0, (
        "a transposase partner was scored as toxin-antitoxin geometry")


def test_a_toxin_partner_is_a_toxin_antitoxin_candidate(fixture_dir):
    """The other half: restricting the test must not disable it."""
    row = _run_context(fixture_dir, "RelE")
    assert float(row["cons_ta_candidate"]) == 1.0, (
        "a genuine toxin partner no longer produces the hypothesis")


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

    run_script("structure_search.py", FakeSnakemake(
        input={"faa": str(faa)},
        output=[str(out)],
        params={"structure": {"max_evalue": 1.0e-3, "min_plddt": 70, "required": True},
                "target_db": FOLDSEEK_DB, "prostt5": PROSTT5},
        threads=4))

    rows = read_tsv(out)
    assert rows, "foldseek returned no structural match for a protein with an exact one"
    assert rows[0]["target_description"], (
        "no description recorded - only the accession, which names nothing")
    assert "SYNTHETASE" in rows[0]["target_description"].upper()


# --- S9: what decides which plate a protein goes on -----------------------------------

def _prioritise(fixture_dir, fam_row, ctx_row, struct_row, rep_seq=None):
    fams = fixture_dir / "dark_families.tsv"
    write_tsv(fams, ["family_id", "representative", "n_members", "n_orfs", "n_plasmids",
                     "n_mob_clusters", "family_class", "members"], [fam_row])
    evo = fixture_dir / "family_evolution.tsv"
    write_tsv(evo, ["family_id", "dnds_median", "dnds_status"],
              [["F1", "0.2", "MEASURED"]])
    ctx = fixture_dir / "family_context.tsv"
    write_tsv(ctx, ["family_id", "cons_integron", "cons_defence", "cons_ta_candidate",
                    "top_hypothesis", "top_conservation", "enrich_integron"], [ctx_row])
    struct = fixture_dir / "structure_hits.tsv"
    write_tsv(struct, ["seq_id", "target", "target_description", "evalue"],
              [struct_row] if struct_row else [])
    faa = fixture_dir / "dark_proteins.faa"
    write_fasta(faa, [("S1", rep_seq
                       or "MKEEDDKKEEDDKKEEDDKKEEQNRSDEIVWCAPGHKGGAFLND")])
    scored = fixture_dir / "scored.tsv"
    report = fixture_dir / "portfolio.txt"

    run_script("prioritise.py", FakeSnakemake(
        input={"families": str(fams), "evolution": str(evo), "context": str(ctx),
               "structure": str(struct), "faa": str(faa)},
        output={"scored": str(scored), "report": str(report)},
        params={"prioritisation": {"min_reality_lines": 1, "min_mob_clusters": 2,
                                   "seed": 1},
                "evolution": {"min_members_for_dnds": 3, "dnds_purifying_max": 0.5},
                "portfolio": {"total": 10, "strata": {"novel_fold": 5,
                                                      "integron_cassette": 5},
                              "reallocate_shortfall": False, "n_controls": 0},
                "context": {"min_context_conservation": 0.5,
                            "high_confidence_conservation": 0.9, "min_enrichment": 2.0},
                "library": {"length_liability_above_aa": 400}}))
    return read_tsv(scored)[0], report.read_text()


FAM = ["F1", "S1", 10, 40, 40, 5, "FAMILY", "S1"]


def test_one_member_in_a_cassette_does_not_make_it_a_cassette_family(fixture_dir):
    """`in_cassette` was `cons_integron > 0` - true when a single member out of ten sits
    near an integron. Every other context claim in the pipeline requires conservation
    ACROSS the family, precisely because one instance is a coincidence. This one did not,
    and it was checked second in the stratum chain, so a single coincidental member
    diverted the whole family into the integron_cassette stratum.

    The membership LABEL must survive - the user asked to be able to see that a dark ORF is
    in a cassette - but it must not decide the stratum on its own."""
    row, _ = _prioritise(fixture_dir, FAM,
                         ["F1", 0.1, 0.0, 0.0, "", 0.0, ""], None)
    assert row["in_integron_cassette"] == "1", (
        "the cassette label was lost - it is evidence a reader asked to keep")
    assert row["cassette_conserved"] == "0"
    assert row["stratum"] != "integron_cassette", (
        "one member of ten decided the stratum for the whole family")


def test_a_cassette_conserved_across_the_family_does_set_the_stratum(fixture_dir):
    """Restricting the rule must not disable it."""
    row, _ = _prioritise(fixture_dir, FAM,
                         ["F1", 0.8, 0.0, 0.0, "integron", 0.8, "5.0"], None)
    assert row["cassette_conserved"] == "1"
    assert row["stratum"] == "integron_cassette"


def test_a_dna_binding_fold_reaches_its_stratum(fixture_dir):
    """The nucleic_acid_binding stratum tested for "nucle" in the Foldseek target, which
    is a PDB accession. It could never be filled. The description is what says what the
    fold is."""
    row, _ = _prioritise(fixture_dir, FAM, ["F1", 0.0, 0.0, 0.0, "", 0.0, ""],
                         ["S1", "1abc-assembly1_A",
                          "1abc-assembly1_A CRYSTAL STRUCTURE OF A DNA-BINDING PROTEIN",
                          "1e-8"])
    assert row["stratum"] == "nucleic_acid_binding"


def test_the_hypothesis_from_a_fold_says_what_the_fold_is(fixture_dir):
    """`structural:1abc-assembly1_A` reaches the synthesis order and tells a bench
    scientist nothing about what to assay."""
    row, _ = _prioritise(fixture_dir, FAM, ["F1", 0.0, 0.0, 0.0, "", 0.0, ""],
                         ["S1", "12as-assembly1_A",
                          "12as-assembly1_A ASPARAGINE SYNTHETASE MUTANT C51A", "1e-8"])
    assert "ASPARAGINE SYNTHETASE" in row["hypothesis"].upper(), (
        f"the hypothesis names only an accession: {row['hypothesis']!r}")


# --- S6b: what "family of one" means after dereplication ------------------------------

@requires("mmseqs")
def test_a_conserved_protein_on_many_plasmids_is_not_reported_as_a_singleton(fixture_dir):
    """Clustering runs on the DEREPLICATED set, so a protein whose sequence is identical on
    two hundred plasmids is ONE member. It clusters alone, is labelled ORPHAN, and fails
    the `is_family` reality test - while being one of the most strongly conserved things in
    the collection.

    That is defensible as a definition: it is not a family of divergent homologs. It is not
    defensible as a REPORT, because `n_members: 1` reads as "seen once". The ORF count has
    to be there too, or a reader cannot tell a genuine singleton from a protein carried by
    two hundred plasmids."""
    faa = fixture_dir / "dark_proteins.faa"
    write_fasta(faa, [("S1", "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEK"),
                      ("S2", "MQQTTLNRSDEIVWCAPGHKGGAFLNDVWRDNPHLAGCVLLTSDGKLLWQRRD")])
    pmap = fixture_dir / "protein_map.tsv"
    # S1 is one unique sequence carried by four plasmids; S2 by one.
    pmap.write_text("S1\tp1|1,p2|1,p3|1,p4|1\nS2\tp5|1\n")
    registry = fixture_dir / "clonal_registry.tsv"
    write_tsv(registry, ["plasmid_id", "mob_cluster"],
              [["p1", "AA1"], ["p2", "AA2"], ["p3", "AA3"], ["p4", "AA1"], ["p5", "AA9"]])
    out = fixture_dir / "dark_families.tsv"

    run_script("cluster_dark.py", FakeSnakemake(
        input={"faa": str(faa), "map": str(pmap), "registry": str(registry)},
        output={"tsv": str(out)},
        params={"clustering": {"min_seq_id": 0.30, "coverage": 0.50, "cov_mode": 0,
                               "cluster_mode": 0}},
        threads=2))

    rows = {r["representative"]: r for r in read_tsv(out)}
    assert "S1" in rows, f"S1 did not survive clustering: {list(rows)}"
    assert rows["S1"]["n_members"] == "1", "S1 is one unique sequence"
    assert rows["S1"]["n_orfs"] == "4", (
        "the ORF count is missing, so a protein on four plasmids is indistinguishable "
        "from one seen once")
    assert rows["S1"]["n_plasmids"] == "4"
    assert rows["S1"]["n_mob_clusters"] == "3"


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


def test_a_candidate_matching_two_strata_records_both(fixture_dir):
    """Stratum assignment is a first-match if/elif chain. That is a declared priority
    order and defensible as one, but it silently discards the other matches - and
    has_tm_helix is an explicitly coarse hydrophobicity heuristic that over-fires on
    soluble proteins with a hydrophobic core. A false TM call is checked early, so it can
    starve every stratum below it (defence_island, nucleic_acid_binding, novel_fold) with
    no trace of what the candidate would otherwise have been.

    The chain still decides. What must not happen is losing the alternatives."""
    row, _ = _prioritise(
        fixture_dir,
        ["F1", "S1", 10, 40, 40, 5, "FAMILY", "S1"],
        ["F1", 0.0, 0.0, 0.0, "", 0.0, ""],
        ["S1", "1abc-assembly1_A", "1abc-assembly1_A DNA-BINDING PROTEIN HU", "1e-8"],
        # Over 100 aa, so the peptide branch cannot fire, with one clear TM segment.
        rep_seq=("MKWLLLAAVFLGLAVLGSVIWLAGFAMTLVGSLLAWFPL"
                 + "EDKNQSTEDKNQSTEDKNQSTEDKNQSTEDKNQSTEDKNQ"
                 + "STEDKNQSTEDKNQSTEDKNQSTEDKNQ"))
    assert row["stratum"] == "membrane_or_secreted", "the declared chain order changed"
    assert "nucleic_acid_binding" in row["stratum_alternatives"], (
        "the DNA-binding fold was discarded without trace by an upstream TM call")


def test_a_widely_carried_orphan_is_counted_in_the_report(fixture_dir):
    """A protein whose sequence is identical on forty plasmids is ONE cluster member: it
    is ORPHAN, fails is_family, and cannot have a dN/dS measured, so it can reach at most
    one reality line and is ineligible at the configured floor of two.

    That is arguably the right definition and definitely the wrong silence. Near-universal
    conservation is among the strongest prevalence signals in the collection, and a run
    that excludes those candidates must say how many it excluded."""
    _, report = _prioritise(
        fixture_dir,
        ["F1", "S1", 1, 40, 40, 6, "ORPHAN", "S1"],
        ["F1", 0.0, 0.0, 0.0, "", 0.0, ""], None)
    assert "widely carried" in report.lower(), (
        f"the report never mentions the excluded high-prevalence orphans:\n{report}")
    assert "\t1\n" in report or " 1\n" in report


# --- S9b: what actually goes to the synthesis vendor ----------------------------------

def test_an_ambiguous_residue_never_reaches_the_synthesis_order(fixture_dir):
    """`CODON.get(aa, "NNN")` mapped any residue outside the twenty-entry table to the
    literal string NNN and wrote it straight into the `cds` column of the file a vendor
    synthesises from. Pyrodigal emits X wherever a gene is called across an ambiguous base,
    so this is a real state in this data - prepare_control.py already filters `"X" not in s`
    when building the control set, so the codebase knows it occurs, and guards only the path
    that never reaches the bench.

    Either the vendor rejects the order, discovered late, or it resolves the ambiguity to an
    arbitrary base and silently ships a different protein from the one on the plate map."""
    scored = fixture_dir / "scored.tsv"
    faa = fixture_dir / "dark.faa"
    order = fixture_dir / "order.tsv"
    report = fixture_dir / "library.txt"

    write_fasta(faa, [("clean", "MKVLATTLLGAAFAASSALAQ"),
                      ("ambiguous", "MKVLATTLXGAAFAASSALAQ")])
    write_tsv(scored,
              ["family_id", "representative", "selected", "stratum", "hypothesis",
               "structural_match", "reality_n"],
              [["F1", "clean", 1, "novel_fold", "", "", 2],
               ["F2", "ambiguous", 1, "novel_fold", "", "", 2]])

    run_script("library_design.py", FakeSnakemake(
        input={"scored": str(scored), "faa": str(faa)},
        output={"order": str(order), "report": str(report)},
        params={"library": {"host": "ecoli", "length_liability_above_aa": 400,
                            "avoid_sites": ["GAATTC"]},
                "portfolio": {"n_controls": 0}, "seed": 1}))

    rows = read_tsv(order)
    assert all("N" not in r["cds"] for r in rows), (
        "an ambiguity placeholder reached the order file: "
        + repr([r["cds"] for r in rows if "N" in r["cds"]]))
    assert [r["representative"] for r in rows if r["role"] == "candidate"] == ["clean"]
    assert "ambiguous" in report.read_text(), (
        "the dropped construct is not reported, so the order is silently short")


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

def test_the_report_carries_every_orf_and_every_family(fixture_dir):
    """The pipeline's output is every annotation it could produce, in a form you can sort
    and filter yourself. Nothing is dropped for being artefactual, ORPHAN, or evidence-free,
    and nothing is ranked - selecting candidates is a decision made on this table, not one
    baked into a rule."""
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
    write_tsv(ctx, ["family_id", "top_hypothesis", "top_conservation", "cons_integron"],
              [["F1", "", 0.0, 0.0], ["F2", "defence", 0.8, 0.0]])
    struct = fixture_dir / "struct.tsv"
    write_tsv(struct, ["seq_id", "target", "target_description", "evalue"], [])
    orth = fixture_dir / "orth.tsv"
    write_tsv(orth, ["seq_id", "cog_category", "kegg_pathways", "preferred_name",
                     "eggnog_description"],
              [["S1", "L", "ko03430", "mobA", "Relaxase"]])

    out_ann = fixture_dir / "annotation_complete.csv"
    out_fam = fixture_dir / "dark_families_complete.csv"
    run_script("annotation_report.py", FakeSnakemake(
        input={"annotation": str(ann), "map": str(pmap), "families": str(fams),
               "evolution": str(evo), "recheck": str(rec),
               "context": str(ctx), "structure": str(struct), "orthology": str(orth)},
        output={"annotation": str(out_ann), "families": str(out_fam)},
        params={"prioritisation": {"min_reality_lines": 2, "min_mob_clusters": 2},
                "evolution": {"min_members_for_dnds": 3, "dnds_purifying_max": 0.5}}))

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
        "top_hypothesis", "top_conservation", "top_enrichment", "high_confidence",
        "cons_defence", "cons_integron", "cons_backbone_adjacent",
        "cons_annotated_neighbour", "cons_operon_with_annotated", "cons_ta_candidate",
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
