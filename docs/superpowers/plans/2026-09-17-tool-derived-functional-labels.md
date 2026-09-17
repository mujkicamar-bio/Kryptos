# Tool-Derived Functional Labels Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the hand-curated 73-name plasmid-backbone list with functional labels taken verbatim from the tools the pipeline already runs, captured in one long table, and make the genomic-context significance test run over groupings of those labels rather than over hand-picked features.

**Architecture:** Three layers, strictly separated. (1) **Capture** - each search stage keeps every identifying field its tool emits, including the ones currently discarded. (2) **Vocabulary** - one long table, `protein_labels.tsv`, with one row per (protein, source, label), where a label always carries the kind of thing it is (`pfam_family`, `cog_category`, `gene_symbol`, `pgap_product`, `macsy_component`, ...). (3) **Grouping** - a pluggable label to category map whose default is the identity, so the enrichment machinery is built and tested now while the biological grouping into replication / mobilisation / conjugation is derived later from the observed vocabulary. Nothing in this plan hand-assigns a biological role.

**Tech Stack:** Python 3.11, Snakemake, pytest, hmmsearch (Pfam-A 38.2), DIAMOND (NCBI swissprot + nr), eggNOG-mapper v2, MacSyFinder 2.1.4, scipy (new dependency, for `fisher_exact`).

## Global Constraints

- Every threshold and parameter needs a published source, or a measured justification for deviating. Record it in `docs/PARAMETER_PROVENANCE.md`.
- Standard international English in all code, comments, docstrings, commit messages and reports. No slang, no abbreviations.
- Surgical changes only. Every changed line must trace to this plan. Do not reformat or improve adjacent code.
- No label may be assigned a biological role by hand anywhere in this plan. Roles are derived from tool vocabularies in a later, separate piece of work.
- Every label row carries its source: tool, database, database version, and tier where applicable. A label without provenance is not admissible.
- Nullable numeric fields carry a companion `<field>_status` column drawn from `darkorf.status`. Never a bare empty cell.
- Pfam release is 38.2 (`data/refs/pfam/Pfam.version.gz`); 30,134 families. The DIAMOND databases are the 2025-03-03 NCBI snapshot at `/sw/data/diamond_databases/Blast/latest/`.
- Tests run with `envs/plasmidann/bin/python -m pytest` from the repository root. `pytest.ini` already sets `pythonpath = src` and `testpaths = tests`.
- Commit after every task. Commit messages end with:
  `Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>`

---

## Background: what each tool actually emits

Verified on this machine, not recalled. Anything marked *discarded* is emitted by the tool and thrown away by the current code.

**hmmsearch `--domtblout`, Pfam-A 38.2** (tiers T1, T2; `workflow/scripts/tier_search.py`)

| field | index | content | status |
|---|---|---|---|
| query name | 3 | Pfam family name, e.g. `RepA_N` | kept as `label` |
| query accession | 4 | `PF06970.16` | **discarded** |
| i-Evalue, ali from/to, tlen, qlen | 12, 17, 18, 2, 5 | statistics and coordinates | kept |
| description of target | 22+ | the *protein's* description, empty for our FASTA | not useful |

The Pfam **description**, **type** (Family/Domain/Repeat/Motif) and **clan** are not in the domtblout at all. They are in `data/refs/pfam/Pfam-A.hmm.dat.gz`, already on disk, as `#=GF DE`, `#=GF TP`, `#=GF CL`. This is why Task 1 exists.

**DIAMOND `--outfmt 6`** (tiers T3 NCBI swissprot, T4 nr)

Current format string: `qseqid stitle qcovhsp scovhsp evalue qstart qend qlen`.

Measured title format for T3, from a real search of CcdB against `swissprot.dmnd`:

```
P62554.1 RecName: Full=Toxin CcdB; AltName: Full=Cytotoxic protein CcdB; ... [Escherichia coli K-12]
```

This is the **NCBI** rendering of Swiss-Prot, not UniProt's. There is no `OS=`/`GN=`/`PE=` structure and **no gene symbol**. What it does give is `RecName: Full=<protein name>`, the `AltName` synonyms, the organism in brackets, and `sseqid` = `P62554.1`, a UniProt accession that a later UniProt join can use for keywords and gene names. `sseqid` is currently **discarded**.

T4 nr titles are NCBI protein titles (`<accession> <product name> [<organism>]`, sometimes prefixed `MULTISPECIES:`). On `WP_` accessions the product name comes from PGAP's controlled product vocabulary, which is the largest single source of role-bearing names in this project. The parser written in Task 5 must tolerate a title that does not match the expected shape rather than assume it; the exact nr title distribution is verified on first real output, not assumed here.

**eggNOG-mapper v2 `.emapper.annotations`** (S4b; `src/plasmidann/orthology.py`)

| column | content | status |
|---|---|---|
| `COG_category` | one or more letters: `L` replication/recombination/repair, `D` partitioning, `V` defence, `U` secretion | kept |
| `Preferred_name` | **gene symbol**: `repA`, `traG`, `mobA`, `ccdB`, `parB`, `tnpA` | kept |
| `Description` | orthologous-group description | kept |
| `eggNOG_OGs` | OG identifiers per taxonomic level, including COG ids such as `COG5527@2` | kept, unparsed |
| `PFAMs` | Pfam names eggNOG assigns to the group | **discarded** |
| `GOs` | GO identifiers | **discarded** |
| `EC` | EC numbers | **discarded** |
| `KEGG_ko` | KO identifiers | **discarded** |
| `KEGG_Pathway` | pathway identifiers | kept |

Gene symbols are the axis the later grouping work will lean on hardest, because symbols are systematic: `rep*`, `tra*`, `trb*`, `mob*`, `par*`, `tnp*`, `ccd*`, `rel*`, `vap*`, `hig*`.

**MacSyFinder 2.1.4 `best_solution.tsv`** (S8a phase 2; `workflow/scripts/defence_systems.py`)

| column | content | status |
|---|---|---|
| `gene_name` | component name from the model, e.g. `RM_Type_II_REase`, `Cas1` | kept as `component` |
| `model_fqn` | `defense-finder-models/.../RM_Type_I` - the system's own name | kept as `system` |
| `sys_id` | system instance | kept |
| `hit_status` | `mandatory` / `accessory` / `neutral` in that model | **discarded** |
| `sys_wholeness` | how complete the called system is, 0-1 | **discarded** |
| `hit_gene_ref`, `hit_profile_cov` | which profile matched, and its coverage | **discarded** |

`macsydata available` reports **CONJScan 2.1.0 installed but never invoked**. It is the native source for conjugation and mobilisation calls (MOB relaxase classes, T4CP, MPF types). Running it is deliberately **out of scope here** and recorded as a decision at the end of this plan, because this plan commits only to capturing what already runs.

**IntegronFinder** (S8b): `annotation` (`intI`, `attC`, `attI`, `protein`), `type_elt`, `model`, integron `type`. Mostly kept.

**Measured evidence that a GO-first approach fails here.** `pfam2go` (2026-07-06) maps only 4 of 16 curated replication families, 1 of 16 conjugation families, 0 of 11 mobilisation families, and exactly 1 family to a toxin-antitoxin-specific term. In a sample of 8 reviewed Swiss-Prot toxin-antitoxin plasmid entries, 8 of 8 carry the UniProt keyword and only 2 of 8 carry `GO:0110001`. GO is therefore not the grouping vocabulary; it is captured as one label kind among several and nothing depends on it.

---

## File Structure

**Create**

| path | responsibility |
|---|---|
| `src/plasmidann/pfam_meta.py` | Parse `Pfam-A.hmm.dat.gz` into family name to {accession, description, type, clan}. |
| `src/plasmidann/labels.py` | The label vocabulary: kinds, extraction of labels from each tool's raw fields, normalisation. |
| `src/plasmidann/categories.py` | Pluggable label to category mapping, loaded from config. Identity by default. |
| `src/plasmidann/enrich.py` | Fisher exact test and Benjamini-Hochberg correction over independent units. |
| `src/plasmidann/controls.py` | `control_recall`, moved out of `backbone.py` unchanged. |
| `config/label_categories.yaml` | The category rules. Ships empty, with the schema documented in place. |
| `workflow/scripts/protein_labels.py` | S4c: build `protein_labels.tsv` from every stage's output. |
| `tests/test_pfam_meta.py`, `tests/test_labels.py`, `tests/test_categories.py`, `tests/test_enrich.py`, `tests/test_controls.py` | One test file per new module. |

**Modify**

| path | change |
|---|---|
| `workflow/scripts/tier_search.py:60-175` | Capture the Pfam accession and the DIAMOND `sseqid`; one new `hits.tsv` column. |
| `src/plasmidann/orthology.py:60-82` | Parse `PFAMs`, `GOs`, `EC`, `KEGG_ko`. |
| `workflow/scripts/orthology.py:77-90` | Write the four new columns. |
| `workflow/scripts/defence_systems.py:50-62` | Capture `hit_status`, `sys_wholeness`, `hit_gene_ref`, `hit_profile_cov`. |
| `workflow/scripts/context_features.py` | Open-vocabulary context in long format, with significance. |
| `workflow/scripts/quality_gate.py:35` | Import `control_recall` from `plasmidann.controls`. |
| `workflow/rules/s3_cascade.smk`, `workflow/rules/s5_s9_targets.smk` | New S4c rule; new context inputs and params. |
| `src/darkorf/schemas.py:11-52` | Declare the `protein_labels` table. |
| `workflow/envs/plasmidann.yaml` | Add `scipy`. |
| `config/config.yaml`, `config/targets.yaml` | Reference paths and versions; suspend the TA stratum. |
| `docs/PARAMETER_PROVENANCE.md`, `docs/PIPELINE_CODE.md` | Document the new parameters and modules. |

**Delete**

| path | reason |
|---|---|
| `src/plasmidann/backbone.py` | The 73-name list is not citable, has measured gaps, and its role assignments have no source. `control_recall` moves to `controls.py` first. |
| `tests/test_backbone.py` | Replaced by `tests/test_controls.py` plus the new label tests. |

---

### Task 1: Pfam family metadata

The Pfam description, type and clan are the first citable label source and are not in the search output. This task makes them available.

**Files:**
- Create: `src/plasmidann/pfam_meta.py`
- Test: `tests/test_pfam_meta.py`

**Interfaces:**
- Consumes: nothing.
- Produces: `parse_pfam_dat(text) -> dict[str, dict]` keyed by family name, each value `{"accession": str, "description": str, "type": str, "clan": str}`; `load(path) -> dict[str, dict]` which accepts a plain or gzipped path.

- [ ] **Step 1: Write the failing test**

Create `tests/test_pfam_meta.py`:

```python
"""Pfam family metadata, which the search output does not carry.

hmmsearch --domtblout gives the family NAME and accession. The description, the type
(Family, Domain, Repeat, Motif) and the clan are in Pfam-A.hmm.dat, and they are what
makes a Pfam hit interpretable as anything other than a string. Pfam 38.2 ships 30,134
families; the parser must handle every one without a special case.
"""
import gzip

from plasmidann import pfam_meta

SAMPLE = """# STOCKHOLM 1.0
#=GF ID   RepA_N
#=GF AC   PF06970.16
#=GF DE   Replication initiator protein A (RepA) N-terminus
#=GF GA   25.5; 25.5;
#=GF TP   Family
#=GF CL   CL0123
//
# STOCKHOLM 1.0
#=GF ID   MobA_MobL
#=GF AC   PF03389.20
#=GF DE   MobA/MobL family protein
#=GF GA   21.7; 21.7;
#=GF TP   Family
//
"""


def test_a_family_with_a_clan_is_parsed_completely():
    meta = pfam_meta.parse_pfam_dat(SAMPLE)

    assert meta["RepA_N"] == {
        "accession": "PF06970.16",
        "description": "Replication initiator protein A (RepA) N-terminus",
        "type": "Family",
        "clan": "CL0123",
    }


def test_a_family_without_a_clan_gets_an_empty_clan_not_a_missing_key():
    """Most families belong to no clan. A missing key would make every consumer guard for
    it, and one that forgot would raise on the common case rather than the rare one."""
    meta = pfam_meta.parse_pfam_dat(SAMPLE)

    assert meta["MobA_MobL"]["clan"] == ""
    assert meta["MobA_MobL"]["description"] == "MobA/MobL family protein"


def test_every_family_in_the_sample_is_present():
    assert set(pfam_meta.parse_pfam_dat(SAMPLE)) == {"RepA_N", "MobA_MobL"}


def test_a_gzipped_file_is_read_transparently(tmp_path):
    """Pfam ships the file gzipped and it is stored gzipped, so requiring the caller to
    decompress it first would mean every caller decompressing a 30,134-entry file."""
    path = tmp_path / "Pfam-A.hmm.dat.gz"
    with gzip.open(path, "wt") as fh:
        fh.write(SAMPLE)

    assert pfam_meta.load(path)["RepA_N"]["accession"] == "PF06970.16"


def test_an_uncompressed_file_is_also_read(tmp_path):
    path = tmp_path / "Pfam-A.hmm.dat"
    path.write_text(SAMPLE)

    assert pfam_meta.load(path)["MobA_MobL"]["type"] == "Family"
```

- [ ] **Step 2: Run test to verify it fails**

Run: `envs/plasmidann/bin/python -m pytest tests/test_pfam_meta.py -v`
Expected: FAIL, `ModuleNotFoundError: No module named 'plasmidann.pfam_meta'`

- [ ] **Step 3: Write minimal implementation**

Create `src/plasmidann/pfam_meta.py`:

```python
"""Pfam family metadata, from the release file rather than from the search output.

hmmsearch --domtblout carries the family name (field 3) and the family accession
(field 4), and nothing else about the family. Its "description of target" column describes
the PROTEIN, not the profile, and is empty for our FASTA.

The description, the type and the clan come from Pfam-A.hmm.dat, which ships with the
release and is already on disk at data/refs/pfam/Pfam-A.hmm.dat.gz for Pfam 38.2. They are
the difference between a hit reading 'RepA_N' and a hit reading 'Replication initiator
protein A (RepA) N-terminus, Family, clan CL0123'.

The clan matters beyond readability: clans group families that are homologous but too
divergent to align as one, so two proteins hitting different families of one clan share an
origin. A grouping built on family names alone would treat them as unrelated.
"""
import gzip
import pathlib

# Stockholm '#=GF <tag> <value>' tags this module keeps, mapped to output key.
_TAGS = {
    "ID": "name",
    "AC": "accession",
    "DE": "description",
    "TP": "type",
    "CL": "clan",
}


def parse_pfam_dat(text):
    """Parse Pfam-A.hmm.dat contents into {family name: metadata}.

    Records are separated by '//'. A family with no '#=GF CL' line belongs to no clan and
    gets an empty string, not a missing key: most families have no clan, so a missing key
    would push a guard into every caller and fail on the common case.
    """
    families = {}
    record = {}
    for line in text.splitlines():
        if line.startswith("//"):
            name = record.pop("name", "")
            if name:
                families[name] = {
                    "accession": record.get("accession", ""),
                    "description": record.get("description", ""),
                    "type": record.get("type", ""),
                    "clan": record.get("clan", ""),
                }
            record = {}
            continue
        if not line.startswith("#=GF "):
            continue
        parts = line[len("#=GF "):].split(None, 1)
        if len(parts) != 2:
            continue
        tag, value = parts
        key = _TAGS.get(tag)
        if key:
            record[key] = value.strip()
    return families


def load(path):
    """Parse Pfam-A.hmm.dat from a path, gzipped or not.

    Decided by the suffix rather than by sniffing the magic bytes: the release file name is
    fixed, and a wrong guess here would be reported as a parse failure on a 30,134-entry
    file rather than as the wrong file being passed.
    """
    path = pathlib.Path(path)
    if path.suffix == ".gz":
        with gzip.open(path, "rt") as fh:
            return parse_pfam_dat(fh.read())
    return parse_pfam_dat(path.read_text())
```

- [ ] **Step 4: Run test to verify it passes**

Run: `envs/plasmidann/bin/python -m pytest tests/test_pfam_meta.py -v`
Expected: PASS, 5 tests

- [ ] **Step 5: Verify against the real Pfam release on disk**

Run:
```bash
envs/plasmidann/bin/python -c "
import sys; sys.path.insert(0, 'src')
from plasmidann import pfam_meta
m = pfam_meta.load('data/refs/pfam/Pfam-A.hmm.dat.gz')
print('families:', len(m))
print('with clan:', sum(1 for v in m.values() if v['clan']))
print('with description:', sum(1 for v in m.values() if v['description']))
for n in ('RepA_N', 'MobA_MobL', 'CagE_TrbE_VirB', 'PIN'):
    print(n, m.get(n))
"
```
Expected: `families: 30134`, every family carrying a description, roughly half carrying a clan, and each of the four named families printing a populated record. If `families` is not 30134, the parser is dropping records and the task is not done.

- [ ] **Step 6: Commit**

```bash
git add src/plasmidann/pfam_meta.py tests/test_pfam_meta.py
git commit -m "$(cat <<'EOF'
feat(labels): read Pfam family description, type and clan from the release file

hmmsearch --domtblout carries only the family name and accession; its
description column describes the protein, not the profile. The description,
type and clan come from Pfam-A.hmm.dat, which is already on disk for Pfam
38.2. The clan is load-bearing rather than cosmetic: clans group families
that are homologous but too divergent to align together, so a grouping built
on family names alone treats them as unrelated.

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>
EOF
)"
```

---

### Task 2: Capture the target accession in the cascade hits

Every DIAMOND hit currently loses the subject accession, which is the only key a later UniProt or RefSeq join can use. Every hmmsearch hit loses the Pfam accession, which is the stable identifier across Pfam releases.

**Files:**
- Modify: `workflow/scripts/tier_search.py:60-175`
- Test: `tests/test_scripts_smoke.py` (add one test)

**Interfaces:**
- Consumes: `plasmidann.cascade` unchanged.
- Produces: `hits.tsv` gains one column, `target_accession`, after `label`. Full column order becomes: `query, label, target_accession, coverage, target_coverage, evalue, informative, is_best, start, end, tier, threshold, max_evalue`.

- [ ] **Step 1: Write the failing test**

Append to `tests/test_scripts_smoke.py`:

```python
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
```

- [ ] **Step 2: Run test to verify it fails**

Run: `PATH=envs/plasmidann/bin:$PATH envs/plasmidann/bin/python -m pytest tests/test_scripts_smoke.py -k subject_accession -v`
Expected: FAIL with `AssertionError: hits.tsv has no target_accession column`

- [ ] **Step 3: Write minimal implementation**

In `workflow/scripts/tier_search.py`, change the `record` signature:

```python
def record(q, label, qcov, tcov, ev, start, end, tlen, accession=""):
```

Inside `record`, in the informative branch, replace the `hit = {...}` assignment with:

```python
        hit = {"query": q, "label": label, "target_accession": accession,
               "coverage": round(qcov, 4),
               "target_coverage": round(tcov, 4), "evalue": ev,
               "informative": True, "start": start, "end": end}
```

In the uninformative branch replace the appended dict with:

```python
        unnamed.setdefault(q, []).append(
            {"query": q, "label": label, "target_accession": accession,
             "coverage": round(qcov, 4),
             "target_coverage": round(tcov, 4), "evalue": ev, "informative": False,
             "start": start, "end": end})
```

In the hmmer branch, replace the `record(...)` call, adding the comment because the index is not self-evident:

```python
            # f[4] is the QUERY accession, which for hmmsearch is the Pfam accession
            # (PF06970.16) - hmmsearch searches profiles against sequences, so the "query"
            # is the profile. f[1] is the target accession, which is the protein's and is
            # '-' for our FASTA. The accession is kept because it is stable across Pfam
            # releases while a family NAME can be changed by a curator.
            record(f[0], f[3],
                   (b - a + 1) / tlen if tlen else 0.0,
                   (int(f[16]) - int(f[15]) + 1) / hlen if hlen else 0.0,
                   f[12], a, b, tlen, accession=f[4])
```

In the DIAMOND branch, add `sseqid` to the format string and unpack it:

```python
        # sseqid is the subject ACCESSION and stitle is free text. Both are kept: the title
        # is the only human-readable identification at T4, and the accession is the only
        # key that joins to UniProt keywords or to RefSeq. Measured on the NCBI swissprot
        # database in use here, the title has no gene symbol, so the accession carries the
        # whole join.
        evalue_flag = "" if max_evalue is None else f"--evalue {max_evalue} "
        cmd = (f"diamond blastp -q {faa} -d {spec['db']} -o {raw} {spec['args']} "
               f"{evalue_flag}--threads {snakemake.threads} "
               f"--max-target-seqs {max_target_seqs} "
               f"--outfmt 6 qseqid sseqid stitle qcovhsp scovhsp evalue qstart qend qlen "
               f"--quiet")
        subprocess.run(cmd, shell=True, check=True)
        for line in open(raw):
            q, sid, title, qc, tc, ev, qs, qe, ql = line.rstrip("\n").split("\t")
            record(q, title, float(qc) / 100, float(tc) / 100, ev,
                   int(qs), int(qe), int(ql), accession=sid)
```

Add the column to `cols`:

```python
cols = ["query", "label", "target_accession", "coverage", "target_coverage", "evalue",
        "informative", "is_best", "start", "end", "tier", "threshold", "max_evalue"]
```

- [ ] **Step 4: Run the test and the whole suite**

Run:
```bash
PATH=envs/plasmidann/bin:$PATH envs/plasmidann/bin/python -m pytest \
  tests/test_scripts_smoke.py -k "subject_accession or cascade or tier" -v
```
Expected: the new test PASSES; no previously passing test fails.

Then the whole suite, because `cascade_resolve.py` reads `hits.tsv`:
```bash
PATH=envs/plasmidann/bin:$PATH envs/plasmidann/bin/python -m pytest -q
```
Expected: no new failures. `csv.DictReader` is keyed by name, so an added column is inert; if anything fails it is reading by position and that is a real defect to fix here.

- [ ] **Step 5: Commit**

```bash
git add workflow/scripts/tier_search.py tests/test_scripts_smoke.py
git commit -m "$(cat <<'EOF'
feat(labels): keep the subject accession on every cascade hit

DIAMOND hits kept only the free-text title and hmmsearch hits kept only the
family name, so the one stable join key each tier produces was discarded.
The Swiss-Prot database in use is NCBI's rendering, whose title has no gene
symbol, so the accession carries the whole join to UniProt keywords. The
Pfam accession is stable across releases while a family name can be changed
by a curator.

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>
EOF
)"
```

---

### Task 3: Capture the full eggNOG column set

eggNOG-mapper already writes gene symbols, Pfam assignments, GO identifiers, EC numbers and KO identifiers on every named protein. Four of those columns are parsed and discarded.

**Files:**
- Modify: `src/plasmidann/orthology.py:60-82`
- Modify: `workflow/scripts/orthology.py:77-90`
- Test: `tests/test_orthology.py` (add tests)

**Interfaces:**
- Consumes: `parse_annotations(text)` as it exists.
- Produces: each record gains `pfams` (list), `gos` (list), `ec` (list), `kegg_ko` (list). `orthology.tsv` gains four comma-joined columns: `pfams`, `gos`, `ec`, `kegg_ko`.

- [ ] **Step 1: Write the failing test**

Append to `tests/test_orthology.py`:

```python
# --- every identifying column, not only the two the first consumer needed -------------

EMAPPER_FULL = "\n".join([
    "## emapper-2.1.12",
    "#query\tseed_ortholog\tevalue\tscore\teggNOG_OGs\tmax_annot_lvl\tCOG_category"
    "\tDescription\tPreferred_name\tGOs\tEC\tKEGG_ko\tKEGG_Pathway\tPFAMs",
    "p1\t83333.b0001\t1e-50\t200.0\tCOG5527@2,2QV1F@1224\t2|Bacteria\tL"
    "\tPlasmid replication initiator protein\trepA"
    "\tGO:0006270,GO:0003677\t2.7.7.7\tko:K02314\tko03030\tRepA_N,Bac_RepA_C",
    "p2\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-",
    "## Total time: 1 s",
])


def test_gene_symbol_pfams_go_ec_and_ko_are_all_parsed():
    """These four columns were parsed and thrown away. The gene symbol is the axis a
    functional grouping leans on hardest, because symbols are systematic: rep*, tra*,
    trb*, mob*, par*, tnp*, ccd*."""
    records = parse_annotations(EMAPPER_FULL)

    assert records["p1"]["preferred_name"] == "repA"
    assert records["p1"]["pfams"] == ["RepA_N", "Bac_RepA_C"]
    assert records["p1"]["gos"] == ["GO:0006270", "GO:0003677"]
    assert records["p1"]["ec"] == ["2.7.7.7"]
    assert records["p1"]["kegg_ko"] == ["ko:K02314"]


def test_the_placeholder_is_absence_in_every_new_column_too():
    """emapper writes a bare '-' for every field it has nothing for. Read literally that
    becomes a Pfam family named '-' and a GO term named '-', and both would aggregate as
    though they were real terms."""
    records = parse_annotations(EMAPPER_FULL)

    assert records["p2"]["pfams"] == []
    assert records["p2"]["gos"] == []
    assert records["p2"]["ec"] == []
    assert records["p2"]["kegg_ko"] == []
    assert records["p2"]["preferred_name"] == ""


def test_a_missing_column_yields_an_empty_list_not_a_failure():
    """An older emapper writes fewer columns. A KeyError here would fail the whole stage
    for a column nothing load-bearing depends on."""
    minimal = "\n".join([
        "#query\tCOG_category\tDescription",
        "p3\tL\tSome protein",
    ])

    records = parse_annotations(minimal)

    assert records["p3"]["pfams"] == []
    assert records["p3"]["gos"] == []
    assert records["p3"]["cog_category"] == "L"
```

- [ ] **Step 2: Run test to verify it fails**

Run: `envs/plasmidann/bin/python -m pytest tests/test_orthology.py -v`
Expected: FAIL with `KeyError: 'pfams'`

- [ ] **Step 3: Write minimal implementation**

In `src/plasmidann/orthology.py`, add a list-splitting helper directly below `_clean`:

```python
def _clean_list(value):
    """Parse a comma-separated emapper field into a list, absence as an empty list.

    emapper separates multi-valued fields with commas and writes a bare '-' when it has
    nothing. An empty list is the honest representation of absence; a list containing '-'
    would aggregate as a real term, which is the placeholder trap in the module docstring
    applied to four more columns.
    """
    cleaned = _clean(value)
    return [v for v in (p.strip() for p in cleaned.split(",")) if v] if cleaned else []
```

In `parse_annotations`, replace the `records[query] = {...}` block with:

```python
        records[query] = {
            # Several single-letter categories are written adjacently, as "EGP". They are
            # kept whole: splitting them into characters invents three annotations, and
            # keeping only the first discards two.
            "cog_category": _clean(row.get("COG_category")),
            "kegg_pathways": [p for p in pathways.split(",") if p] if pathways else [],
            "preferred_name": _clean(row.get("Preferred_name")),
            "description": _clean(row.get("Description")),
            "eggnog_ogs": _clean(row.get("eggNOG_OGs")),
            # Four columns emapper always writes and this parser used to discard. They are
            # the controlled-vocabulary identifiers on the annotated fraction: the Pfam
            # families eggNOG assigns to the orthologous group, GO terms, EC numbers and
            # KEGG orthologs. A functional grouping has to be derived from vocabularies
            # like these rather than written by hand, so discarding them made the grouping
            # impossible to build.
            "pfams": _clean_list(row.get("PFAMs")),
            "gos": _clean_list(row.get("GOs")),
            "ec": _clean_list(row.get("EC")),
            "kegg_ko": _clean_list(row.get("KEGG_ko")),
        }
```

In `workflow/scripts/orthology.py`, replace the output block:

```python
cols = ["seq_id", "cog_category", "kegg_pathways", "preferred_name", "eggnog_description",
        "eggnog_ogs", "pfams", "gos", "ec", "kegg_ko"]
with open(snakemake.output[0], "w", newline="") as out:
    w = csv.DictWriter(out, fieldnames=cols, delimiter="\t")
    w.writeheader()
    for sid in sorted(named):
        r = records.get(sid, {})
        w.writerow({"seq_id": sid,
                    "cog_category": r.get("cog_category", ""),
                    "kegg_pathways": ",".join(r.get("kegg_pathways", [])),
                    "preferred_name": r.get("preferred_name", ""),
                    "eggnog_description": r.get("description", ""),
                    "eggnog_ogs": r.get("eggnog_ogs", ""),
                    "pfams": ",".join(r.get("pfams", [])),
                    "gos": ",".join(r.get("gos", [])),
                    "ec": ",".join(r.get("ec", [])),
                    "kegg_ko": ",".join(r.get("kegg_ko", []))})

n_kegg = sum(1 for r in records.values() if r["kegg_pathways"])
n_symbol = sum(1 for r in records.values() if r["preferred_name"])
print(f"orthology: queried={n_query} annotated={len(records)} with_kegg={n_kegg} "
      f"with_gene_symbol={n_symbol}")
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `envs/plasmidann/bin/python -m pytest tests/test_orthology.py tests/test_scripts_smoke.py -q`
Expected: PASS, no new failures

- [ ] **Step 5: Commit**

```bash
git add src/plasmidann/orthology.py workflow/scripts/orthology.py tests/test_orthology.py
git commit -m "$(cat <<'EOF'
feat(labels): keep eggNOG gene symbols, Pfam, GO, EC and KO columns

emapper writes these on every named protein and the parser discarded four of
them. They are the controlled-vocabulary identifiers on the annotated
fraction, and a functional grouping has to be derived from vocabularies like
these rather than written by hand. The gene symbol matters most, because
symbols are systematic where free-text descriptions are not.

The '-' placeholder is parsed as an empty list in each new column: a list
containing '-' would aggregate as a real term.

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>
EOF
)"
```

---

### Task 4: Capture the MacSyFinder hit detail

MacSyFinder says whether a component is mandatory or accessory to its model, and how complete the called system is. Both are discarded, and both change how much a system call is worth.

**Files:**
- Modify: `workflow/scripts/defence_systems.py:22-62`
- Test: `tests/test_scripts_smoke.py` (add one test)

**Interfaces:**
- Consumes: MacSyFinder `best_solution.tsv`.
- Produces: `defence_systems.tsv` gains `hit_status`, `sys_wholeness`, `hit_gene_ref`, `hit_profile_cov`. Existing columns keep their names and meaning.

- [ ] **Step 1: Write the failing test**

Append to `tests/test_scripts_smoke.py`:

```python
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
    run_script("defence_systems.py", FakeSnakemake(
        input={"faa": str(faa), "map": str(mapping)},
        output={"tsv": str(out)},
        params={"models_dir": str(fixture_dir / "models"), "skip_run": True},
        threads=1))

    rows = {r["orf_id"]: r for r in read_tsv(out)}
    assert rows["p1|1"]["hit_status"] == "mandatory"
    assert rows["p1|2"]["hit_status"] == "accessory"
    assert rows["p1|1"]["sys_wholeness"] == "1.000"
    assert rows["p1|1"]["hit_profile_cov"] == "0.95"
    assert rows["p1|1"]["system"] == "defense-finder-models/Defense/RM_Type_II"
```

- [ ] **Step 2: Run test to verify it fails**

Run: `envs/plasmidann/bin/python -m pytest tests/test_scripts_smoke.py -k component_status -v`
Expected: FAIL, because `macsyfinder` is invoked unconditionally

- [ ] **Step 3: Write minimal implementation**

In `workflow/scripts/defence_systems.py`, replace the `subprocess.run(...)` block with:

```python
# skip_run exists for the parser test, which pre-populates the output tree. It is never
# set by the workflow: a missing MacSyFinder run in production must fail, not be skipped.
if not getattr(snakemake.params, "skip_run", False):
    subprocess.run(
        f"macsyfinder --models-dir {snakemake.params.models_dir} "
        f"--models defense-finder-models all "
        f"--sequence-db {snakemake.input.faa} "
        f"--db-type gembase --replicon-topology circular "
        f"--worker {snakemake.threads} --out-dir {outdir} --mute",
        shell=True, check=True)
```

Replace the `rows.append({...})` block with:

```python
            rows.append({"orf_id": src["orf_id"], "plasmid_id": src["plasmid_id"],
                         "gembase_id": gid, "system": rec.get("model_fqn", ""),
                         "system_id": rec.get("sys_id", ""),
                         "component": rec.get("gene_name", ""),
                         "hit_evalue": rec.get("hit_i_eval", ""),
                         # A mandatory component of a complete system and a neutral
                         # component of a fragment are different evidence and arrived as
                         # identical rows. hit_status is the model's own declaration of
                         # which it is, and sys_wholeness is how much of the model was
                         # found - both are MacSyFinder's numbers, not our
                         # reinterpretation of them.
                         "hit_status": rec.get("hit_status", ""),
                         "sys_wholeness": rec.get("sys_wholeness", ""),
                         "hit_gene_ref": rec.get("hit_gene_ref", ""),
                         "hit_profile_cov": rec.get("hit_profile_cov", "")})
```

Replace the writer's `fieldnames` with:

```python
    w = csv.DictWriter(out, fieldnames=["orf_id", "plasmid_id", "gembase_id", "system",
                                        "system_id", "component", "hit_evalue",
                                        "hit_status", "sys_wholeness", "hit_gene_ref",
                                        "hit_profile_cov"],
                       delimiter="\t")
```

Do not change the `outdir` line: the script already globs `outdir.rglob("best_solution.tsv")` under `<output parent>/phase2`, which is where the test writes. Read it to confirm rather than editing it.

- [ ] **Step 4: Run tests to verify they pass**

Run: `envs/plasmidann/bin/python -m pytest tests/test_scripts_smoke.py -k "component_status or defence" -q`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add workflow/scripts/defence_systems.py tests/test_scripts_smoke.py
git commit -m "$(cat <<'EOF'
feat(labels): keep MacSyFinder component status and system wholeness

A mandatory component of a complete system and a neutral component of a
fragment are different evidence and arrived as identical rows. hit_status,
sys_wholeness, hit_gene_ref and hit_profile_cov are MacSyFinder's own
declarations and were discarded.

skip_run is added for the parser test, which pre-populates the output tree.
The workflow never sets it: a missing run in production must fail.

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>
EOF
)"
```

---

### Task 5: The label vocabulary

One module that knows what kinds of label exist and how to extract each kind from a stage's output row. This is the layer that makes the later grouping possible, and the layer that must never assign a biological role.

**Files:**
- Create: `src/plasmidann/labels.py`
- Test: `tests/test_labels.py`

**Interfaces:**
- Consumes: `plasmidann.pfam_meta` from Task 1.
- Produces:
  - `KINDS` - a frozenset of the declared label kinds.
  - `labels_from_hit(row, pfam=None) -> list[dict]` where `row` is one `hits.tsv` row; each returned dict is `{"kind": str, "label": str, "accession": str}`.
  - `labels_from_orthology(row) -> list[dict]` for one `orthology.tsv` row.
  - `labels_from_defence(row) -> list[dict]` for one `defence_systems.tsv` row.
  - `labels_from_integron(row) -> list[dict]` for one integron row.
  - `parse_ncbi_title(title) -> dict` with keys `accession`, `product`, `organism`.

- [ ] **Step 1: Write the failing test**

Create `tests/test_labels.py`:

```python
"""The label vocabulary: what a tool said, with the kind of statement it made.

A functional grouping of plasmid proteins into replication, mobilisation and conjugation
cannot be written by hand. Measured on the curated list this replaces: 73 names, of which
pfam2go maps 4 of 16 replication families, 1 of 16 conjugation families and 0 of 11
mobilisation families. Pfam alone holds roughly 30 T4SS families against the 15 listed. A
hand list cannot reach the scope.

So every label is captured verbatim from the tool, tagged with its KIND, and grouped later
from the observed vocabulary. The kind is not decoration: 'repA' as a gene symbol and
'RepA_N' as a Pfam family are different statements with different reliability, and a
grouping that cannot tell them apart cannot be audited.
"""
from plasmidann import labels


def test_a_pfam_hit_yields_the_family_the_accession_the_description_and_the_clan():
    """One hmmsearch hit is four statements, and three of them were being discarded."""
    pfam = {"RepA_N": {"accession": "PF06970.16",
                       "description": "Replication initiator protein A (RepA) N-terminus",
                       "type": "Family", "clan": "CL0123"}}
    row = {"tier": "T1", "label": "RepA_N", "target_accession": "PF06970.16",
           "informative": "True"}

    out = labels.labels_from_hit(row, pfam=pfam)
    by_kind = {d["kind"]: d["label"] for d in out}

    assert by_kind["pfam_family"] == "RepA_N"
    assert by_kind["pfam_description"] == "Replication initiator protein A (RepA) N-terminus"
    assert by_kind["pfam_clan"] == "CL0123"
    assert all(d["accession"] == "PF06970.16" for d in out)


def test_a_pfam_family_with_no_clan_emits_no_clan_label():
    """An empty clan is absence. A label reading '' would group every clanless family
    together, which is the largest false group it is possible to create here."""
    pfam = {"MobA_MobL": {"accession": "PF03389.20",
                          "description": "MobA/MobL family protein",
                          "type": "Family", "clan": ""}}
    row = {"tier": "T1", "label": "MobA_MobL", "target_accession": "PF03389.20",
           "informative": "True"}

    kinds = {d["kind"] for d in labels.labels_from_hit(row, pfam=pfam)}

    assert "pfam_clan" not in kinds
    assert "pfam_family" in kinds


def test_a_pfam_family_absent_from_the_release_still_yields_the_family_name():
    """A family name the installed release does not know is a real event - a database
    mismatch - and the name must survive so the mismatch is visible in the table rather
    than silently dropping the hit."""
    row = {"tier": "T2", "label": "Not_In_Pfam", "target_accession": "PF99999.1",
           "informative": "True"}

    out = labels.labels_from_hit(row, pfam={})

    assert [d["kind"] for d in out] == ["pfam_family"]


def test_a_diamond_hit_yields_the_product_name_and_the_accession():
    """The NCBI title is 'ACC RecName: Full=Toxin CcdB; AltName: ... [organism]' for
    swissprot and 'ACC product name [organism]' for nr. The product name is the label; the
    organism is not a functional statement and is not emitted as one."""
    row = {"tier": "T3", "informative": "True",
           "target_accession": "P62554.1",
           "label": "P62554.1 RecName: Full=Toxin CcdB; AltName: Full=Protein LetD "
                    "[Escherichia coli K-12]"}

    by_kind = {d["kind"]: d["label"] for d in labels.labels_from_hit(row)}

    assert by_kind["swissprot_product"] == "Toxin CcdB"
    assert "organism" not in by_kind


def test_an_nr_hit_yields_the_pgap_product_name():
    """On WP_ accessions the product name comes from PGAP's controlled vocabulary, which
    will be the largest single source of role-bearing names in the collection."""
    row = {"tier": "T4", "informative": "True",
           "target_accession": "WP_000813620.1",
           "label": "WP_000813620.1 type II toxin-antitoxin system RelE/ParE family "
                    "toxin [Escherichia coli]"}

    by_kind = {d["kind"]: d["label"] for d in labels.labels_from_hit(row)}

    assert by_kind["pgap_product"] == ("type II toxin-antitoxin system RelE/ParE family "
                                       "toxin")


def test_the_multispecies_prefix_is_removed_from_an_nr_product():
    """'MULTISPECIES: relaxase' and 'relaxase' are the same product, and keeping the
    prefix would split every widespread protein into two labels - exactly the proteins a
    grouping most needs to see as one."""
    row = {"tier": "T4", "informative": "True", "target_accession": "WP_1.1",
           "label": "WP_1.1 MULTISPECIES: conjugal transfer protein TraG "
                    "[Enterobacteriaceae]"}

    by_kind = {d["kind"]: d["label"] for d in labels.labels_from_hit(row)}

    assert by_kind["pgap_product"] == "conjugal transfer protein TraG"


def test_an_uninformative_hit_yields_no_label():
    """'hypothetical protein' names nothing. It is evidence that someone has seen the
    protein, which the cascade already records, and it must never enter the functional
    vocabulary - a 'hypothetical protein' category would be the most enriched context
    feature in the run."""
    row = {"tier": "T4", "informative": "False", "target_accession": "WP_2.1",
           "label": "WP_2.1 hypothetical protein [Escherichia coli]"}

    assert labels.labels_from_hit(row) == []


def test_an_unparsable_title_still_yields_the_whole_title_as_the_product():
    """nr titles are not uniform and this parser will meet shapes it was not shown. The
    honest fallback is the whole string: a dropped label is invisible, a strange label is
    not."""
    row = {"tier": "T4", "informative": "True", "target_accession": "",
           "label": "something with no recognisable structure"}

    by_kind = {d["kind"]: d["label"] for d in labels.labels_from_hit(row)}

    assert by_kind["pgap_product"] == "something with no recognisable structure"


def test_orthology_yields_the_gene_symbol_and_every_controlled_identifier():
    row = {"seq_id": "p1", "cog_category": "L", "preferred_name": "repA",
           "eggnog_description": "Plasmid replication initiator protein",
           "pfams": "RepA_N,Bac_RepA_C", "gos": "GO:0006270", "ec": "2.7.7.7",
           "kegg_ko": "ko:K02314", "kegg_pathways": "ko03030", "eggnog_ogs": "COG5527@2"}

    pairs = {(d["kind"], d["label"]) for d in labels.labels_from_orthology(row)}

    assert ("gene_symbol", "repA") in pairs
    assert ("cog_category", "L") in pairs
    assert ("cog_id", "COG5527") in pairs
    assert ("eggnog_pfam", "RepA_N") in pairs
    assert ("eggnog_pfam", "Bac_RepA_C") in pairs
    assert ("go", "GO:0006270") in pairs
    assert ("ec", "2.7.7.7") in pairs
    assert ("kegg_ko", "ko:K02314") in pairs


def test_a_multi_letter_cog_category_is_split_into_one_label_per_letter():
    """COG categories are written adjacently as 'LKV', and each letter is a separate
    category. Kept whole it is a category named 'LKV' that no COG release defines; split,
    it is three real ones. The orthology TABLE keeps the string whole for provenance; the
    VOCABULARY needs the individual categories."""
    row = {"seq_id": "p1", "cog_category": "LKV", "preferred_name": "",
           "eggnog_description": "", "pfams": "", "gos": "", "ec": "",
           "kegg_ko": "", "kegg_pathways": "", "eggnog_ogs": ""}

    cats = {d["label"] for d in labels.labels_from_orthology(row)
            if d["kind"] == "cog_category"}

    assert cats == {"L", "K", "V"}


def test_an_empty_orthology_row_yields_no_labels():
    row = {"seq_id": "p2", "cog_category": "", "preferred_name": "",
           "eggnog_description": "", "pfams": "", "gos": "", "ec": "",
           "kegg_ko": "", "kegg_pathways": "", "eggnog_ogs": ""}

    assert labels.labels_from_orthology(row) == []


def test_a_defence_row_yields_the_system_name_and_the_component():
    """model_fqn is a path; the system name is its last element. The full path is kept as
    the accession so the model set that made the call stays visible."""
    row = {"orf_id": "p1|1", "system": "defense-finder-models/Defense/RM_Type_II",
           "component": "RM_Type_II_REase", "hit_status": "mandatory",
           "sys_wholeness": "1.000"}

    by_kind = {d["kind"]: d["label"] for d in labels.labels_from_defence(row)}

    assert by_kind["macsy_system"] == "RM_Type_II"
    assert by_kind["macsy_component"] == "RM_Type_II_REase"


def test_an_integron_row_yields_the_element_annotation():
    row = {"orf_id": "p1|4", "annotation": "intI", "integron_type": "complete"}

    by_kind = {d["kind"]: d["label"] for d in labels.labels_from_integron(row)}

    assert by_kind["integron_element"] == "intI"
    assert by_kind["integron_type"] == "complete"


def test_every_emitted_kind_is_declared():
    """An undeclared kind is a label nothing downstream knows how to group, and it would
    be discovered as a missing category rather than as a bug here."""
    pfam = {"RepA_N": {"accession": "PF1.1", "description": "d", "type": "Family",
                       "clan": "CL0123"}}
    emitted = set()
    emitted |= {d["kind"] for d in labels.labels_from_hit(
        {"tier": "T1", "label": "RepA_N", "target_accession": "PF1.1",
         "informative": "True"}, pfam=pfam)}
    emitted |= {d["kind"] for d in labels.labels_from_hit(
        {"tier": "T3", "label": "P1.1 RecName: Full=Toxin CcdB [E. coli]",
         "target_accession": "P1.1", "informative": "True"})}
    emitted |= {d["kind"] for d in labels.labels_from_hit(
        {"tier": "T4", "label": "W1.1 relaxase [E. coli]", "target_accession": "W1.1",
         "informative": "True"})}
    emitted |= {d["kind"] for d in labels.labels_from_orthology(
        {"seq_id": "p", "cog_category": "L", "preferred_name": "repA",
         "eggnog_description": "x", "pfams": "A", "gos": "GO:1", "ec": "1.1.1.1",
         "kegg_ko": "ko:K1", "kegg_pathways": "", "eggnog_ogs": "COG1@2"})}
    emitted |= {d["kind"] for d in labels.labels_from_defence(
        {"orf_id": "o", "system": "m/Defense/X", "component": "X_c",
         "hit_status": "mandatory", "sys_wholeness": "1"})}
    emitted |= {d["kind"] for d in labels.labels_from_integron(
        {"orf_id": "o", "annotation": "intI", "integron_type": "complete"})}

    assert emitted <= labels.KINDS, f"undeclared kinds: {sorted(emitted - labels.KINDS)}"


def test_no_extraction_function_assigns_a_biological_role():
    """The standing constraint, asserted rather than trusted. Nothing in this module may
    decide that RepA_N is 'replication': that grouping is derived later from the observed
    vocabulary, and a role appearing here would be a hand-curated list growing back."""
    forbidden = {"replication", "mobilisation", "conjugation", "partition",
                 "transposition", "toxin_antitoxin", "restriction_modification",
                 "backbone"}

    assert not (labels.KINDS & forbidden)
```

- [ ] **Step 2: Run test to verify it fails**

Run: `envs/plasmidann/bin/python -m pytest tests/test_labels.py -v`
Expected: FAIL, `ModuleNotFoundError: No module named 'plasmidann.labels'`

- [ ] **Step 3: Write minimal implementation**

Create `src/plasmidann/labels.py`:

```python
"""The functional-label vocabulary: what each tool said, and what kind of statement it is.

WHY THIS REPLACES A CURATED LIST

An earlier version of this pipeline classified plasmid proteins into replication,
mobilisation, conjugation, partition, transposition, toxin-antitoxin and
restriction-modification using a hand-written list of 73 Pfam family names. It could not
work, and the measurements say why:

  * Pfam holds roughly 30 T4SS families; the list named 15. It named no MobB and no MobD.
  * pfam2go, the published Pfam-to-GO mapping, covers 4 of the 16 listed replication
    families, 1 of 16 conjugation families and 0 of 11 mobilisation families - so the gaps
    are in the curated mapping too, not only in the list.
  * Nine of the 73 names did not exist in Pfam-A at all, so those entries had never once
    matched anything, and no output could have revealed it.
  * The role assignments had no source. They were one person's reading.

A hand list cannot reach the scope of a 30,134-family database, and its gaps are silent.
So no role is assigned here. Every label is taken verbatim from the tool that produced it,
tagged with the KIND of statement it is, and the grouping into biological categories is
derived later from the vocabulary actually observed in the data.

WHY THE KIND IS LOAD-BEARING

'repA' as a gene symbol from eggNOG and 'RepA_N' as a Pfam family from hmmsearch are
different statements about a protein, from different evidence, with different reliability.
A vocabulary that could not tell them apart could not be audited, and a category built from
both would be impossible to describe in a methods section. The kind travels with every
label for that reason, and downstream code groups on (kind, label), never on the label
alone.

WHAT IS DELIBERATELY NOT A LABEL

  * Organism names. A functional grouping is not a taxonomy, and an organism in a title is
    the source of the reference sequence, not a statement about the query.
  * Uninformative descriptions. 'hypothetical protein' names nothing; the cascade already
    records it as dark evidence. Admitted here it would become the most frequent, and
    therefore most apparently enriched, category in the collection.
  * The empty string in any field. Absence is absence, and a label of '' would group every
    protein missing that field into one enormous false category.
"""

# Every kind of label this module can emit. A kind not listed here is a statement nothing
# downstream knows how to group, so emitting one is a bug rather than a new feature.
#
# Deliberately absent: any biological role name. Roles are derived from these kinds later;
# a role appearing in this set would be the curated list growing back.
KINDS = frozenset({
    # hmmsearch against Pfam-A, tiers T1 and T2
    "pfam_family",        # RepA_N            the family name, as tier_search records it
    "pfam_description",   # 'Replication initiator protein A (RepA) N-terminus'
    "pfam_clan",          # CL0123            families too divergent to align as one
    # DIAMOND, tier T3 (NCBI rendering of Swiss-Prot) and tier T4 (nr)
    "swissprot_product",  # 'Toxin CcdB'      from 'RecName: Full=...'
    "pgap_product",       # 'conjugal transfer protein TraG'
    # eggNOG-mapper, S4b
    "gene_symbol",        # repA, traG, mobA  the most systematic axis available
    "cog_category",       # L, D, V           one label per letter
    "cog_id",             # COG5527
    "eggnog_pfam",        # the Pfam families eggNOG assigns to the orthologous group
    "eggnog_description",
    "go",                 # GO:0006270
    "ec",                 # 2.7.7.7
    "kegg_ko",            # ko:K02314
    # MacSyFinder, S8a phase 2
    "macsy_system",       # RM_Type_II        the model's own system name
    "macsy_component",    # RM_Type_II_REase
    # IntegronFinder, S8b
    "integron_element",   # intI, attC, attI
    "integron_type",      # complete, In0, CALIN
})

# Tiers whose label is a Pfam family name rather than a sequence title. Taken from the tier
# method in config/cascade.yaml: the hmmer tiers search Pfam-A and report family names.
_PFAM_TIERS = frozenset({"T1", "T2"})

# The tier whose database is NCBI's rendering of Swiss-Prot. Its titles carry
# 'RecName: Full=<name>;' rather than the UniProt 'OS=/GN=/PE=' structure, and notably no
# gene symbol - verified against the installed database, not assumed.
_SWISSPROT_TIERS = frozenset({"T3"})

# NCBI marks a title shared by several organisms with this prefix. 'MULTISPECIES: relaxase'
# and 'relaxase' are the same product, so keeping the prefix would split every widespread
# protein into two labels - exactly the proteins a grouping most needs to see as one.
_MULTISPECIES = "MULTISPECIES:"


def _informative(row):
    """Whether a hits.tsv row was judged informative by the cascade.

    The column is written by csv as the string 'True' or 'False'. Comparing the raw value
    to a boolean would make every row falsy and silently empty the vocabulary, so the
    accepted spellings are tested explicitly.
    """
    return row.get("informative") in (True, "True", "true", "1", 1)


def parse_ncbi_title(title):
    """Split an NCBI protein title into accession, product name and organism.

    Two shapes are handled, both verified against the databases in use:

        swissprot   P62554.1 RecName: Full=Toxin CcdB; AltName: Full=Protein LetD [E. coli]
        nr          WP_000813620.1 type II toxin-antitoxin system RelE/ParE toxin [E. coli]

    A title that matches neither returns the whole string as the product. That is the
    honest fallback: nr titles are not uniform, this parser will meet shapes it was not
    shown, and a dropped label is invisible while a strange one is not.
    """
    text = (title or "").strip()
    if not text:
        return {"accession": "", "product": "", "organism": ""}

    organism = ""
    if text.endswith("]") and "[" in text:
        head, _, tail = text.rpartition("[")
        organism = tail[:-1].strip()
        text = head.strip()

    accession, product = "", text
    first, _, rest = text.partition(" ")
    # An accession is the leading token when it looks like one: NCBI and UniProt
    # accessions carry a version suffix and no spaces. Testing for the '.' rather than
    # matching a pattern per database keeps this from needing a rule per accession style.
    if rest and "." in first:
        accession, product = first, rest.strip()

    if product.startswith(_MULTISPECIES):
        product = product[len(_MULTISPECIES):].strip()

    if product.startswith("RecName:"):
        # 'RecName: Full=Toxin CcdB; AltName: Full=Protein LetD' - the recommended name is
        # the first Full=, and the AltName synonyms are dropped: they are the same protein
        # under other names, so admitting them would multiply one statement into several.
        body = product[len("RecName:"):].strip()
        head = body.split(";", 1)[0].strip()
        product = head[len("Full="):].strip() if head.startswith("Full=") else head

    return {"accession": accession, "product": product, "organism": organism}


def labels_from_hit(row, pfam=None):
    """Labels from one hits.tsv row.

    `pfam` is the mapping from plasmidann.pfam_meta.load; when it is None or lacks the
    family, only the family name is emitted. A family the installed release does not know
    is a database mismatch, and the name has to survive so the mismatch appears in the
    table rather than silently removing the hit.
    """
    if not _informative(row):
        return []

    tier = row.get("tier", "")
    label = (row.get("label") or "").strip()
    accession = (row.get("target_accession") or "").strip()
    if not label:
        return []

    out = []
    if tier in _PFAM_TIERS:
        out.append({"kind": "pfam_family", "label": label, "accession": accession})
        meta = (pfam or {}).get(label)
        if meta:
            if meta.get("description"):
                out.append({"kind": "pfam_description", "label": meta["description"],
                            "accession": accession})
            if meta.get("clan"):
                out.append({"kind": "pfam_clan", "label": meta["clan"],
                            "accession": accession})
        return out

    parsed = parse_ncbi_title(label)
    if not parsed["product"]:
        return []
    kind = "swissprot_product" if tier in _SWISSPROT_TIERS else "pgap_product"
    return [{"kind": kind, "label": parsed["product"],
             "accession": accession or parsed["accession"]}]


def _split(value):
    """Comma-separated table field to a list, with blanks removed."""
    return [v for v in (p.strip() for p in (value or "").split(",")) if v]


def labels_from_orthology(row):
    """Labels from one orthology.tsv row.

    Every controlled-vocabulary identifier eggNOG assigns, plus the gene symbol. The gene
    symbol is singled out because it is the most systematic axis available anywhere in this
    pipeline: symbols are assigned in families (rep*, tra*, trb*, mob*, par*, tnp*, ccd*,
    rel*, vap*, hig*), which free-text descriptions are not.
    """
    out = []
    symbol = (row.get("preferred_name") or "").strip()
    if symbol:
        out.append({"kind": "gene_symbol", "label": symbol, "accession": ""})

    # Adjacent letters, as 'LKV', are three separate COG categories. The orthology table
    # keeps the string whole for provenance; the vocabulary needs each category, because
    # 'LKV' is not a category any COG release defines.
    for letter in (row.get("cog_category") or "").strip():
        if letter.isalpha():
            out.append({"kind": "cog_category", "label": letter, "accession": ""})

    # eggNOG_OGs is 'COG5527@2,2QV1F@1224' - the group identifier at each taxonomic level.
    # The COG identifiers are the citable ones; the numeric eggNOG groups have no published
    # functional description, so admitting them would add labels nothing can interpret.
    for group in _split(row.get("eggnog_ogs")):
        name = group.split("@", 1)[0]
        if name.startswith("COG"):
            out.append({"kind": "cog_id", "label": name, "accession": ""})

    for family in _split(row.get("pfams")):
        out.append({"kind": "eggnog_pfam", "label": family, "accession": ""})
    for term in _split(row.get("gos")):
        out.append({"kind": "go", "label": term, "accession": ""})
    for number in _split(row.get("ec")):
        out.append({"kind": "ec", "label": number, "accession": ""})
    for ko in _split(row.get("kegg_ko")):
        out.append({"kind": "kegg_ko", "label": ko, "accession": ""})

    description = (row.get("eggnog_description") or "").strip()
    if description:
        out.append({"kind": "eggnog_description", "label": description, "accession": ""})
    return out


def labels_from_defence(row):
    """Labels from one defence_systems.tsv row.

    `system` is MacSyFinder's model_fqn, a path such as
    'defense-finder-models/Defense/RM_Type_II'. The system name is its last element; the
    full path is kept as the accession so the model set that made the call stays visible in
    the table, which is what makes the call traceable to its publication.
    """
    out = []
    fqn = (row.get("system") or "").strip()
    if fqn:
        out.append({"kind": "macsy_system", "label": fqn.rsplit("/", 1)[-1],
                    "accession": fqn})
    component = (row.get("component") or "").strip()
    if component:
        out.append({"kind": "macsy_component", "label": component, "accession": fqn})
    return out


def labels_from_integron(row):
    """Labels from one integron row: the element type and the integron class."""
    out = []
    element = (row.get("annotation") or "").strip()
    if element:
        out.append({"kind": "integron_element", "label": element, "accession": ""})
    integron_type = (row.get("integron_type") or "").strip()
    if integron_type:
        out.append({"kind": "integron_type", "label": integron_type, "accession": ""})
    return out
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `envs/plasmidann/bin/python -m pytest tests/test_labels.py -v`
Expected: PASS, 15 tests

- [ ] **Step 5: Commit**

```bash
git add src/plasmidann/labels.py tests/test_labels.py
git commit -m "$(cat <<'EOF'
feat(labels): add the tool-derived label vocabulary

Every functional label is taken verbatim from the tool that produced it and
tagged with the kind of statement it is, so a biological grouping can be
derived later from the vocabulary actually observed rather than written by
hand. No role is assigned here, and a test asserts that no kind is a role
name, so the curated list cannot grow back.

The kind is load-bearing: 'repA' as a gene symbol and 'RepA_N' as a Pfam
family are different statements with different reliability, and downstream
code groups on (kind, label) rather than on the label alone.

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>
EOF
)"
```

---

### Task 6: The protein_labels table

One long table, one row per (protein, source, label). This is the substrate the later grouping reads, and the only place a reader has to look to see what the tools said about a protein.

**Files:**
- Create: `workflow/scripts/protein_labels.py`
- Modify: `workflow/rules/s3_cascade.smk` (add rule `protein_labels` after rule `orthology`)
- Modify: `config/config.yaml` (reference paths and database versions)
- Modify: `src/darkorf/schemas.py:11-52` (declare the table)
- Test: `tests/test_scripts_smoke.py` (add two tests), `tests/test_schemas.py` (add one test)

**Interfaces:**
- Consumes: `plasmidann.labels` and `plasmidann.pfam_meta` from Tasks 1 and 5.
- Produces: `results/s4c/protein_labels.tsv` with columns
  `protein_id, source, tier, kind, label, accession, evidence_evalue, evidence_coverage, database, database_version`.
  One row per distinct (protein_id, source, kind, label); the best supporting statistics are kept when a label is seen more than once.

- [ ] **Step 1: Write the failing tests**

Append to `tests/test_scripts_smoke.py`:

```python
# --- S4c: one long table of what every tool said --------------------------------------

def test_protein_labels_gathers_every_source_into_one_long_table(fixture_dir):
    """The substrate for the functional grouping. A wide table cannot hold it: the
    vocabulary is open, Pfam alone has 30,134 families, and the grouping is derived from
    the labels observed rather than declared in advance."""
    import gzip

    hits = fixture_dir / "hits.tsv"
    write_tsv(hits, ["query", "label", "target_accession", "coverage", "target_coverage",
                     "evalue", "informative", "is_best", "start", "end", "tier",
                     "threshold", "max_evalue"],
              [["s1", "RepA_N", "PF06970.16", 0.9, 0.95, "1e-40", "True", 1, 1, 100,
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
        fh.write("# STOCKHOLM 1.0\n#=GF ID   RepA_N\n#=GF AC   PF06970.16\n"
                 "#=GF DE   Replication initiator protein A (RepA) N-terminus\n"
                 "#=GF TP   Family\n#=GF CL   CL0123\n//\n")

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
              [["s1", "RepA_N", "PF06970.16", 0.9, 0.95, "1e-40", "True", 1, 1, 100,
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
```

Append to `tests/test_schemas.py`:

```python
def test_the_protein_labels_table_is_declared():
    """The substrate for the functional grouping is a declared table, not a loose file.
    An undeclared column here is a label kind nothing downstream can group."""
    cols = schemas.columns("protein_labels")

    for required in ("protein_id", "source", "kind", "label", "database",
                     "database_version"):
        assert required in cols, f"protein_labels has no {required} column"
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `envs/plasmidann/bin/python -m pytest tests/test_scripts_smoke.py -k protein_labels tests/test_schemas.py -v`
Expected: FAIL - script not found, and `KeyError: "no schema declared for table 'protein_labels'"`

- [ ] **Step 3: Declare the table**

In `src/darkorf/schemas.py`, add inside `TABLES`, after the `"proteins"` entry:

```python
    # Spec section 62.3 addendum - one row per (protein, source, label). LONG, not wide:
    # the label vocabulary is open (Pfam-A 38.2 alone has 30,134 families, and nr product
    # names are unbounded), so it cannot be columns. This table is the substrate from which
    # functional categories - replication, mobilisation, conjugation - are derived; the
    # hand-curated family list it replaces covered 15 of roughly 30 T4SS families and
    # assigned roles with no source.
    "protein_labels": [
        "protein_id", "source", "tier", "kind", "label", "accession",
        "evidence_evalue", "evidence_coverage", "database", "database_version",
    ],
```

- [ ] **Step 4: Write the script**

Create `workflow/scripts/protein_labels.py`:

```python
"""S4c: one long table of every functional label every tool produced.

WHY THIS TABLE EXISTS

The question it serves is "is this protein replication, mobilisation, or conjugation", and
the answer cannot come from a hand-written list. The list this replaces held 73 Pfam family
names against roughly 30 T4SS families in Pfam alone, named no MobB and no MobD, included
nine names that do not exist in Pfam-A at all, and assigned every role with no source.

So the labels are collected verbatim, with their kind and their provenance, and the
grouping into biological categories is derived afterwards from the vocabulary observed
here. That order matters: a category built from the labels the data actually contains can be
described in a methods section, and a category built from recollection cannot.

WHY IT IS LONG AND NOT WIDE

The vocabulary is open. Pfam-A 38.2 has 30,134 families, nr product names are unbounded,
and a protein carries a different number of labels from every source. One row per
(protein, source, kind, label) is the only shape that holds that without a column per
family, and it is the shape a grouping step reads naturally: select the distinct labels of
one kind, decide their categories, join back.

DEDUPLICATION

A label seen several times for one protein - the same Pfam family hit by two tiers, the
same product name from ten nr subjects - is one statement, not ten. Rows are keyed on
(protein_id, source, kind, label) and the best supporting statistics are kept, because an
unmerged table would let a widespread label outvote a rare one purely by copy number when
the categories are counted.
"""
import _ctx  # noqa: F401
import csv
import sys

from plasmidann import labels, pfam_meta

# Which database and version produced each source, for the provenance columns. A label
# without its database release cannot be reproduced, and a category built on it cannot be
# described in a paper.
_DATABASE = {
    "pfam": ("Pfam-A", "pfam_version"),
    "swissprot": ("NCBI swissprot", "swissprot_version"),
    "nr": ("NCBI nr", "nr_version"),
    "eggnog": ("eggNOG", "eggnog_version"),
}

# Tier to source. The hmmer tiers search Pfam-A; T3 is NCBI's rendering of Swiss-Prot and
# T4 is nr. Taken from config/cascade.yaml rather than inferred from the label text.
_TIER_SOURCE = {"T1": "pfam", "T2": "pfam", "T3": "swissprot", "T4": "nr"}

params = snakemake.params
pfam = pfam_meta.load(snakemake.input.pfam_dat)

# (protein_id, source, kind, label) -> row. Identity is fixed on first sight and later
# sightings only improve the statistics, so a label's evidence is the strongest seen rather
# than the last row read.
rows = {}


def _as_float(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return float("inf")


def add(protein_id, source, tier, entry, evalue="", coverage=""):
    """Admit one label, merging it with any previous sighting of the same statement."""
    if entry["kind"] not in labels.KINDS:
        # An undeclared kind is a statement nothing downstream can group. Failing here
        # names the kind; discovering it later looks like a missing category.
        raise SystemExit(
            f"protein_labels: undeclared label kind {entry['kind']!r} from source "
            f"{source!r}. Add it to plasmidann.labels.KINDS or stop emitting it.")
    database, version_key = _DATABASE[source]
    key = (protein_id, source, entry["kind"], entry["label"])
    existing = rows.get(key)
    if existing is None:
        rows[key] = {
            "protein_id": protein_id, "source": source, "tier": tier,
            "kind": entry["kind"], "label": entry["label"],
            "accession": entry.get("accession", ""),
            "evidence_evalue": evalue, "evidence_coverage": coverage,
            "database": database, "database_version": str(params.get(version_key, "")),
        }
        return
    if _as_float(evalue) < _as_float(existing["evidence_evalue"]):
        existing["evidence_evalue"] = evalue
        existing["evidence_coverage"] = coverage
        existing["accession"] = entry.get("accession", "") or existing["accession"]


# ------------------------------------------------------------------------------------
# Cascade hits: Pfam families from the hmmer tiers, product names from the DIAMOND tiers.
# ------------------------------------------------------------------------------------
n_hits = 0
for path in snakemake.input.hits:
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            n_hits += 1
            source = _TIER_SOURCE.get(row.get("tier", ""))
            if source is None:
                # A tier this script does not know about would otherwise be silently
                # dropped, which is how a whole database's labels go missing.
                raise SystemExit(
                    f"protein_labels: tier {row.get('tier')!r} has no source mapping. "
                    "Add it to _TIER_SOURCE when a tier is added to config/cascade.yaml.")
            for entry in labels.labels_from_hit(row, pfam=pfam):
                add(row["query"], source, row.get("tier", ""), entry,
                    evalue=row.get("evalue", ""), coverage=row.get("coverage", ""))

# ------------------------------------------------------------------------------------
# Orthology: gene symbols and every controlled identifier eggNOG assigns.
# ------------------------------------------------------------------------------------
n_orth = 0
with open(snakemake.input.orthology, newline="") as fh:
    for row in csv.DictReader(fh, delimiter="\t"):
        n_orth += 1
        for entry in labels.labels_from_orthology(row):
            add(row["seq_id"], "eggnog", "S4b", entry)

cols = ["protein_id", "source", "tier", "kind", "label", "accession",
        "evidence_evalue", "evidence_coverage", "database", "database_version"]
with open(snakemake.output.tsv, "w", newline="") as out:
    writer = csv.DictWriter(out, fieldnames=cols, delimiter="\t")
    writer.writeheader()
    for key in sorted(rows):
        writer.writerow(rows[key])

by_kind = {}
for (_, _, kind, _) in rows:
    by_kind[kind] = by_kind.get(kind, 0) + 1
distinct_labels = len({(k, l) for (_, _, k, l) in rows})

print(f"protein_labels: read {n_hits} hits and {n_orth} orthology rows -> "
      f"{len(rows)} label rows on {len({k[0] for k in rows})} proteins, "
      f"{distinct_labels} distinct (kind, label) pairs")
for kind in sorted(by_kind):
    print(f"  {kind:<22} {by_kind[kind]}")

# An empty table here is always a bug: the cascade named a large fraction of the
# collection, and every named protein carries at least one label by construction.
if not rows:
    sys.exit("protein_labels: no labels extracted from any source - the functional "
             "grouping has no substrate. Check that hits.tsv carries a tier column and "
             "that orthology.tsv is not empty.")
```

- [ ] **Step 5: Wire the rule and the config**

In `workflow/rules/s3_cascade.smk`, append after rule `orthology`:

```python
rule protein_labels:
    """S4c: every functional label every tool produced, in one long table.

    The substrate for grouping proteins into replication, mobilisation and conjugation.
    Long rather than wide because the vocabulary is open: Pfam-A 38.2 alone has 30,134
    families. Nothing here assigns a biological role; the grouping is derived from the
    labels observed, which is what makes it describable in a methods section.
    """
    input:
        hits=expand(f"{OUT}/s3/{{tier}}/{{cshard}}/hits.tsv",
                    tier=TIER_IDS, cshard=CASCADE_SHARDS),
        orthology=f"{OUT}/s4b/orthology.tsv",
        pfam_dat=config["pfam_dat"],
    output:
        tsv=f"{OUT}/s4c/protein_labels.tsv",
    params:
        pfam_version=config["pfam_version"],
        swissprot_version=config["swissprot_version"],
        nr_version=config["nr_version"],
        eggnog_version=config["eggnog_version"],
    resources:
        mem_mb=32000,
        runtime=240,
    log:
        f"{OUT}/logs/s4c/protein_labels.log",
    conda:
        "../envs/plasmidann.yaml"
    script:
        "../scripts/protein_labels.py"
```

In `config/config.yaml`, add alongside the other reference paths:

```yaml
  # Pfam family metadata: the description, type and clan that --domtblout does not carry.
  # Release read from data/refs/pfam/Pfam.version.gz, which reports 38.2, 30,134 families.
  pfam_dat: data/refs/pfam/Pfam-A.hmm.dat.gz
  pfam_version: "38.2"
  # The NCBI BLAST database snapshot behind the DIAMOND indexes. The path
  # /sw/data/diamond_databases/Blast/latest resolves to .../20250303-003013.
  swissprot_version: "2025-03-03"
  nr_version: "2025-03-03"
  eggnog_version: "5.0.2"
```

- [ ] **Step 6: Run tests to verify they pass**

Run: `envs/plasmidann/bin/python -m pytest tests/test_scripts_smoke.py -k protein_labels tests/test_schemas.py tests/test_workflow_wiring.py -v`
Expected: PASS. If `tests/test_workflow_wiring.py` checks that every config key a rule reads is declared, it will catch a missing key here; add the key rather than relaxing the test.

- [ ] **Step 7: Commit**

```bash
git add workflow/scripts/protein_labels.py workflow/rules/s3_cascade.smk \
        config/config.yaml src/darkorf/schemas.py \
        tests/test_scripts_smoke.py tests/test_schemas.py
git commit -m "$(cat <<'EOF'
feat(labels): add S4c protein_labels, one long table of every tool label

The substrate for grouping proteins into replication, mobilisation and
conjugation. Long rather than wide because the vocabulary is open: Pfam-A
38.2 has 30,134 families and nr product names are unbounded, so a column per
family is not possible and a hand-picked subset is what failed before.

Every row carries its database and release, because a label without its
source cannot be reproduced and a category built on it cannot be described
in a methods section. Labels repeated across tiers are merged, so a
widespread label cannot outvote a rare one by copy number.

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>
EOF
)"
```

---

### Task 7: The pluggable label-to-category map

The biological grouping is deferred, but the seam it plugs into is built and tested now. The default is the identity: each (kind, label) is its own category. That makes every consumer work today and makes the later grouping a configuration change rather than a rewrite.

**Files:**
- Create: `src/plasmidann/categories.py`
- Create: `config/label_categories.yaml`
- Test: `tests/test_categories.py`

**Interfaces:**
- Consumes: `plasmidann.labels.KINDS` from Task 5.
- Produces: `load_rules(path) -> dict`; `load_rules_from_text(text) -> dict`; `assign(kind, label, rules) -> list[tuple[str, str]]`; `CategoryMap(rules)` with `.categories_for(kind, label) -> list[tuple[str, str]]` and `.is_identity -> bool`.
- Every returned pair is `(category, subcategory)`. The category is the FLAT axis the significance test runs on; the subcategory is carried alongside and never tested, so a finer grouping can be re-cut from the table without re-running the pipeline.

- [ ] **Step 1: Write the failing test**

Create `tests/test_categories.py`:

```python
"""Grouping labels into functional categories - the seam, not the grouping.

The grouping itself (replication, mobilisation, conjugation, ...) is derived from the label
vocabulary after a full annotation run, because that is the only way it can cover the
scope: Pfam alone holds roughly 30 T4SS families, and the hand list this replaces named 15.

What is built here is the mechanism and its default. With no rules configured, each
(kind, label) is its own category, so the enrichment machinery downstream is exercised on
real data from the first run and the grouping becomes a configuration change rather than a
rewrite.
"""
import pytest

from plasmidann import categories


def test_with_no_rules_a_label_is_its_own_category_qualified_by_kind():
    """The identity default. The kind qualifies the category name because 'repA' the gene
    symbol and 'RepA_N' the Pfam family are different statements, and a category that
    merged them could not be audited."""
    rules = categories.load_rules(None)

    assert categories.assign("pfam_family", "RepA_N", rules) == \
        [("pfam_family:RepA_N", "")]
    assert categories.assign("gene_symbol", "repA", rules) == [("gene_symbol:repA", "")]


def test_the_identity_map_reports_itself_as_the_identity():
    """A run whose categories are still the identity has not had the grouping applied, and
    a report must be able to say so rather than presenting per-label enrichment as if it
    were per-category."""
    assert categories.CategoryMap(categories.load_rules(None)).is_identity is True


def test_an_exact_rule_assigns_a_category():
    rules = categories.load_rules_from_text("""
categories:
  replication:
    exact:
      pfam_family: [RepA_N, Rep_3]
      gene_symbol: [repA]
""")

    assert categories.assign("pfam_family", "RepA_N", rules) == [("replication", "")]
    assert categories.assign("gene_symbol", "repA", rules) == [("replication", "")]


def test_a_label_matching_no_rule_falls_back_to_itself():
    """A label with no category must stay visible. Dropping it would make the vocabulary
    shrink silently as rules are added, which is how the previous list's nine dead entries
    went unnoticed for a whole version."""
    rules = categories.load_rules_from_text("""
categories:
  replication:
    exact:
      pfam_family: [RepA_N]
""")

    assert categories.assign("pfam_family", "MobA_MobL", rules) == \
        [("pfam_family:MobA_MobL", "")]


def test_a_label_may_belong_to_several_categories():
    """A relaxase is both mobilisation and, on a conjugative plasmid, part of transfer. A
    map that forced one category per label would decide that by rule order rather than by
    biology."""
    rules = categories.load_rules_from_text("""
categories:
  mobilisation:
    exact:
      pfam_family: [MobA_MobL]
  transfer:
    exact:
      pfam_family: [MobA_MobL]
""")

    assert sorted(categories.assign("pfam_family", "MobA_MobL", rules)) == \
        [("mobilisation", ""), ("transfer", "")]


def test_a_rule_naming_an_undeclared_kind_is_an_error():
    """A rule on a kind nothing emits can never fire, and it would look like a category
    with no members rather than a typo. Nine such entries in the list this replaces had
    never once matched anything."""
    with pytest.raises(ValueError, match="unknown label kind"):
        categories.load_rules_from_text("""
categories:
  replication:
    exact:
      pfam_familly: [RepA_N]
""")


def test_matching_is_exact_and_case_sensitive():
    """Substring matching would be a trap: 'PIN' and 'rve' are Pfam families whose names
    appear inside many unrelated descriptions, and a mislabelled protein sends someone to
    the bench to test the wrong thing."""
    rules = categories.load_rules_from_text("""
categories:
  toxin_antitoxin:
    exact:
      pfam_family: [PIN]
""")

    assert categories.assign("pfam_family", "PIN", rules) == [("toxin_antitoxin", "")]
    assert categories.assign("pfam_family", "PIN_2", rules) == \
        [("pfam_family:PIN_2", "")]
    assert categories.assign("pfam_family", "pin", rules) == [("pfam_family:pin", "")]


def test_the_map_memoises_repeated_lookups():
    """The table has one row per protein per label; the same (kind, label) is looked up
    millions of times. A per-call scan of every rule would be the same quadratic class
    already fixed twice in this pipeline, in tier_search and in context_features."""
    cmap = categories.CategoryMap(categories.load_rules_from_text("""
categories:
  replication:
    exact:
      pfam_family: [RepA_N]
"""))

    first = cmap.categories_for("pfam_family", "RepA_N")
    second = cmap.categories_for("pfam_family", "RepA_N")

    assert first == second == [("replication", "")]
    assert cmap.categories_for("pfam_family", "RepA_N") is second


def test_a_subcategory_is_carried_but_does_not_change_the_category():
    """Flat for the test, finer for the record. A subcategory lets 'replication/initiator'
    and 'replication/control' be told apart in the table without splitting the category the
    significance test runs on - more categories means more tests and weaker FDR power, and
    re-cutting a finer grouping must not require re-running the pipeline."""
    rules = categories.load_rules_from_text("""
categories:
  replication:
    exact:
      pfam_family: [Rop]
    sub:
      pfam_family:
        RepA_N: initiator
        Rop: control
""")

    assert categories.assign("pfam_family", "Rop", rules) == [("replication", "control")]


def test_a_label_with_a_subcategory_but_no_exact_rule_is_not_categorised():
    """A sub entry is a refinement of a category membership, not a way to declare one. If
    it created membership by itself, the sub block would become a second, undocumented
    place where categories are assigned."""
    rules = categories.load_rules_from_text("""
categories:
  replication:
    sub:
      pfam_family:
        RepA_N: initiator
""")

    assert categories.assign("pfam_family", "RepA_N", rules) == \
        [("pfam_family:RepA_N", "")]


def test_the_shipped_config_file_loads_and_is_still_the_identity():
    """config/label_categories.yaml ships with no rules, because the grouping is derived
    from a full annotation run. This test fails the day rules are added without updating
    it, which is the point at which a reviewer should look at them."""
    rules = categories.load_rules("config/label_categories.yaml")

    assert categories.CategoryMap(rules).is_identity is True, (
        "rules have been added to config/label_categories.yaml - update this test and "
        "record the provenance of each category in docs/PARAMETER_PROVENANCE.md")
```

- [ ] **Step 2: Run test to verify it fails**

Run: `envs/plasmidann/bin/python -m pytest tests/test_categories.py -v`
Expected: FAIL, `ModuleNotFoundError: No module named 'plasmidann.categories'`

- [ ] **Step 3: Write minimal implementation**

Create `src/plasmidann/categories.py`:

```python
"""Grouping tool labels into functional categories.

WHAT THIS IS AND IS NOT

This module is the mechanism for grouping labels, and its default is to do no grouping at
all. The biological categories - replication, mobilisation, conjugation, partition,
transposition, toxin-antitoxin - are NOT defined here. They are derived from the label
vocabulary after a full annotation run, from config/label_categories.yaml, and each one
needs a recorded justification like every other parameter in this pipeline.

The reason for that order is measured. The hand-written list this replaces covered 15 of
roughly 30 T4SS families in Pfam, named no MobB and no MobD, and contained nine family
names that do not exist in Pfam-A at all - entries that had never matched anything in a
whole version of the pipeline, invisibly. Rules written against the labels the data actually
contains can be checked against the table; rules written from recollection cannot.

THE IDENTITY DEFAULT

With no rules configured, every (kind, label) is its own category, named 'kind:label'. Two
consequences, both wanted:

  * The enrichment machinery downstream runs on real data from the first annotation run,
    rather than waiting for a grouping that does not exist yet.
  * Applying the grouping later is a configuration change, not a rewrite.

The kind qualifies the category name because 'repA' as a gene symbol and 'RepA_N' as a Pfam
family are different statements about a protein. A category that merged them could not be
audited, and a reader could not tell which evidence produced it.

FALLBACK, NOT EXCLUSION

A label matching no rule falls back to being its own category. It is never dropped. A
vocabulary that shrank silently as rules were added is how the previous list's nine dead
entries survived unnoticed, and it is the one failure mode this design exists to prevent.
"""
import pathlib

import yaml

from plasmidann.labels import KINDS


def load_rules_from_text(text):
    """Parse category rules from YAML text.

    Shape:

        categories:
          replication:
            exact:
              pfam_family: [RepA_N, Rep_3]
              gene_symbol: [repA]

    A rule naming a kind that plasmidann.labels does not emit raises: such a rule can never
    fire, and it would present as a category with no members rather than as the typo it is.
    """
    document = yaml.safe_load(text) or {}
    declared = document.get("categories") or {}
    index = {}
    for category, spec in declared.items():
        spec = spec or {}
        # 'sub' refines a membership that 'exact' declares; it never creates one, so that
        # there is exactly one place in the file where a category gains a member.
        subs = {}
        for kind, mapping in (spec.get("sub") or {}).items():
            _check_kind(kind, category)
            for value, name in (mapping or {}).items():
                subs[(kind, value)] = str(name)
        for kind, values in (spec.get("exact") or {}).items():
            _check_kind(kind, category)
            for value in values or []:
                index.setdefault((kind, value), []).append(
                    (category, subs.get((kind, value), "")))
    return {"exact": index}


def _check_kind(kind, category):
    """A rule on a kind nothing emits can never fire.

    It would present as a category with no members rather than as the typo it is, which is
    exactly how nine entries of the list this replaces survived a whole version.
    """
    if kind not in KINDS:
        raise ValueError(
            f"unknown label kind {kind!r} in category {category!r}; "
            f"known kinds: {', '.join(sorted(KINDS))}")


def load_rules(path):
    """Load category rules from a path, or the empty (identity) rules when path is None."""
    if path is None:
        return {"exact": {}}
    return load_rules_from_text(pathlib.Path(path).read_text())


def assign(kind, label, rules):
    """The (category, subcategory) pairs of one label.

    Exact, case-sensitive matching. Substring matching would be a trap: 'PIN' and 'rve' are
    Pfam family names that appear inside many unrelated descriptions, and a mislabelled
    protein sends someone to the bench to test the wrong thing.

    A label may belong to several categories, and they are all returned. A relaxase is
    mobilisation and, on a conjugative plasmid, part of transfer; forcing one category per
    label would decide that by rule order rather than by biology.

    The subcategory is a finer name carried alongside the category and never tested. The
    significance test runs on the FLAT category, because more categories means more tests
    and weaker FDR power; the subcategory is in the table so a finer grouping can be re-cut
    from the output without re-running the pipeline.
    """
    matched = rules["exact"].get((kind, label))
    return list(matched) if matched else [(f"{kind}:{label}", "")]


class CategoryMap:
    """A memoised view of the rules, for the per-row lookups.

    protein_labels has one row per protein per label, so the same (kind, label) is resolved
    millions of times. Scanning the rules per call would be the same quadratic class already
    fixed twice in this pipeline, in tier_search and in context_features.
    """

    def __init__(self, rules):
        self._rules = rules
        self._cache = {}

    @property
    def is_identity(self):
        """Whether no grouping is configured.

        A report has to be able to say that a run's categories are still per-label, rather
        than presenting per-label enrichment as though it were per-category.
        """
        return not self._rules["exact"]

    def categories_for(self, kind, label):
        key = (kind, label)
        if key not in self._cache:
            self._cache[key] = assign(kind, label, self._rules)
        return self._cache[key]
```

Create `config/label_categories.yaml`:

```yaml
# Grouping tool labels into functional categories.
#
# THIS FILE SHIPS EMPTY, DELIBERATELY.
#
# The categories this project needs - replication, mobilisation, conjugation, partition,
# transposition, toxin-antitoxin, restriction-modification - are derived from the label
# vocabulary in results/s4c/protein_labels.tsv after a full annotation run. They are not
# written from recollection, because that was tried and measured:
#
#   * the hand-written list this replaces held 73 Pfam family names, against roughly 30
#     T4SS families in Pfam alone;
#   * it named no MobB and no MobD;
#   * nine of its 73 names do not exist in Pfam-A at all, so those entries had never once
#     matched anything, through a whole version of the pipeline, invisibly;
#   * pfam2go, the published Pfam-to-GO mapping, covers 4 of its 16 replication families,
#     1 of 16 conjugation families and 0 of 11 mobilisation families - so the gap is in the
#     curated mappings too, not only in one person's reading.
#
# Until rules are added, every (kind, label) is its own category and the enrichment
# machinery runs per label. That is a valid run: it reports which individual labels a dark
# protein's neighbourhood is enriched for, which is a weaker but honest statement.
#
# WHEN RULES ARE ADDED
#
# Each category needs, in docs/PARAMETER_PROVENANCE.md:
#   1. the query over protein_labels.tsv that produced its member list;
#   2. the published source for the grouping (a database's own category vocabulary, a
#      review, or a nomenclature paper) - not a judgement;
#   3. the count of labels and proteins it covers, and what it leaves unassigned.
#
# SHAPE
#
# categories:
#   <category name>:
#     exact:
#       <label kind>: [<label>, <label>, ...]
#     sub:                      # optional, and never creates membership
#       <label kind>:
#         <label>: <subcategory name>
#
# The CATEGORY is the flat axis the significance test runs on. The SUBCATEGORY is carried
# beside it in the output and never tested, so 'replication/initiator' and
# 'replication/control' can be told apart in the table while the test still sees one
# 'replication' category. More categories means more tests and weaker FDR power, and a
# finer grouping must be re-cuttable from the output without re-running the pipeline.
#
# A label listed only under 'sub' is NOT categorised: 'sub' refines a membership that
# 'exact' declares, so that there is exactly one place in this file where a category gains
# a member.
#
# Matching is exact and case-sensitive. Label kinds are those declared in
# plasmidann.labels.KINDS; a rule on any other kind is an error, because it could never
# fire and would present as an empty category rather than as a typo.
#
# A label may appear in several categories. A relaxase is mobilisation and, on a
# conjugative plasmid, part of transfer, and forcing one category per label would decide
# that by rule order.

categories: {}
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `envs/plasmidann/bin/python -m pytest tests/test_categories.py -v`
Expected: PASS, 9 tests

- [ ] **Step 5: Commit**

```bash
git add src/plasmidann/categories.py config/label_categories.yaml tests/test_categories.py
git commit -m "$(cat <<'EOF'
feat(labels): add the pluggable label-to-category map, identity by default

The mechanism for grouping labels into functional categories, with no
grouping configured. The categories themselves are derived from the observed
vocabulary after a full annotation run, because the alternative was measured:
the list this replaces covered 15 of roughly 30 T4SS families and contained
nine names absent from Pfam-A that had never matched anything.

The identity default means the enrichment machinery runs on real data from
the first run and the grouping becomes a configuration change. A label
matching no rule falls back to being its own category and is never dropped:
a vocabulary that shrinks silently as rules are added is how the previous
list's dead entries survived a whole version.

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>
EOF
)"
```

---

### Task 8: Significance over independent units

The current context statistic is a ratio of two rates with no test attached, computed over family members. Members of one family are homologs on plasmids that may be near-identical, so the effective sample size is not the member count. This task adds a test, fixes the unit, and corrects for the number of categories tested.

**Files:**
- Create: `src/plasmidann/enrich.py`
- Modify: `workflow/envs/plasmidann.yaml` (add `scipy`)
- Modify: `docs/PARAMETER_PROVENANCE.md`
- Test: `tests/test_enrich.py`

**Interfaces:**
- Consumes: `darkorf.status` for the status vocabulary.
- Produces:
  - `fisher_enrichment(k, n, K, N, min_units=2) -> dict` with keys `observed_rate`, `background_rate`, `enrichment`, `odds_ratio`, `p_value`, `status`, where `k` of `n` family units carry the category and `K` of `N` corpus units do.
  - `benjamini_hochberg(p_values) -> list` returning q-values in input order, non-numeric entries passed through unchanged.

- [ ] **Step 1: Add scipy to the environments**

In `workflow/envs/plasmidann.yaml`, add to the dependency list, matching the file's existing ordering and comment style:

```yaml
  # scipy for fisher_exact. The context enrichment was a bare ratio of two rates with no
  # test attached, so a family with three members at 100% and one with three hundred at
  # 90% were indistinguishable in the output.
  - scipy
```

Install it into the locally built environment the tests use:

```bash
envs/plasmidann/bin/python -m pip install --no-input scipy
envs/plasmidann/bin/python -c "import scipy; print('scipy', scipy.__version__)"
```
Expected: a version prints. If pip is unavailable in that environment, use
`conda install -p envs/plasmidann -y scipy` instead.

- [ ] **Step 2: Write the failing test**

Create `tests/test_enrich.py`:

```python
"""Significance for the genomic-context association, over independent units.

Three defects in the statistic this replaces:

  1. No test. Enrichment was observed rate over background rate, so a family with three
     members all beside a relaxase scored identically to one with three hundred.
  2. The wrong unit. The rate was over family MEMBERS, which are homologs, often on
     near-identical plasmids. A protein family present on forty copies of one sequenced
     plasmid is one observation, not forty.
  3. No correction. Six hand-picked features needed none; the open label vocabulary means
     thousands of tests per family, where an uncorrected p-value is meaningless.
"""
import pytest

from plasmidann import enrich


def test_a_strong_association_is_significant():
    """Eighteen of twenty independent plasmids carrying this family have the category in
    the neighbourhood, against 5% of the corpus."""
    result = enrich.fisher_enrichment(k=18, n=20, K=500, N=10_000)

    assert float(result["p_value"]) < 1e-10
    assert result["enrichment"] > 3
    assert result["status"] == "SUCCESS"


def test_an_association_at_the_background_rate_is_not_significant():
    """The small-plasmid trap: on a six-gene plasmid a plus or minus three neighbourhood is
    the whole molecule, so everything co-occurs with everything. Observed equal to
    background must come out unremarkable, whatever the absolute rate."""
    result = enrich.fisher_enrichment(k=5, n=10, K=5_000, N=10_000)

    assert float(result["p_value"]) > 0.5
    assert result["enrichment"] == pytest.approx(1.0, abs=0.01)


def test_a_single_unit_cannot_be_tested():
    """One plasmid is an anecdote. Returning a p-value for n=1 would let a family found
    once rank alongside one found on two hundred independent plasmids."""
    result = enrich.fisher_enrichment(k=1, n=1, K=100, N=10_000)

    assert result["status"] == "TOO_FEW_MEMBERS"
    assert result["p_value"] == ""


def test_a_category_absent_from_the_corpus_is_reported_not_divided_by_zero():
    """A zero background with a non-zero observation is the honest infinite enrichment,
    and the caller decides how to rank it."""
    result = enrich.fisher_enrichment(k=3, n=10, K=0, N=10_000)

    assert result["enrichment"] == float("inf")


def test_no_observation_and_no_background_is_no_signal():
    result = enrich.fisher_enrichment(k=0, n=10, K=0, N=10_000)

    assert result["enrichment"] == 1.0


def test_an_impossible_count_is_an_error_not_a_silent_result():
    """k > n means the caller's units disagree with its counts, which would produce a
    plausible-looking p-value from nonsense."""
    with pytest.raises(ValueError):
        enrich.fisher_enrichment(k=11, n=10, K=100, N=10_000)
    with pytest.raises(ValueError):
        enrich.fisher_enrichment(k=1, n=10, K=100, N=10)


def test_the_family_is_tested_against_the_rest_of_the_corpus(): 
    """The two-by-two table compares the family's plasmids against the plasmids that are
    NOT the family's. Including them in both margins would test the family against a
    background it contributes to, which shrinks any real effect - and for a family that
    covers most of the corpus it would shrink it to nothing."""
    result = enrich.fisher_enrichment(k=10, n=10, K=10, N=1_000)

    assert float(result["p_value"]) < 1e-10, (
        "a category unique to this family came out unremarkable, so the family is being "
        "tested against a background that includes it")


def test_benjamini_hochberg_preserves_input_order():
    """The caller joins q-values back to rows by position, so a sorted return would
    silently attach every q-value to the wrong category."""
    q = enrich.benjamini_hochberg([0.5, 0.001, 0.04, 0.2])

    assert len(q) == 4
    assert q[1] == min(q), "the smallest p-value did not get the smallest q-value"


def test_benjamini_hochberg_matches_the_published_procedure():
    """Benjamini and Hochberg 1995, J R Stat Soc B 57:289. With m=4, the q-value of the
    i-th smallest p is min over j >= i of (m/j) * p_j, which is monotone by construction."""
    q = enrich.benjamini_hochberg([0.01, 0.02, 0.03, 0.04])

    assert q == pytest.approx([0.04, 0.04, 0.04, 0.04])


def test_benjamini_hochberg_ignores_untested_entries():
    """A category that could not be tested has no p-value, and treating an empty string as
    zero would give it the strongest q-value in the table."""
    q = enrich.benjamini_hochberg([0.01, "", 0.02])

    assert q[1] == ""
    assert q[0] != ""


def test_no_q_value_exceeds_one():
    q = enrich.benjamini_hochberg([0.9, 0.95, 0.99])

    assert all(v <= 1.0 for v in q)
```

- [ ] **Step 3: Run test to verify it fails**

Run: `envs/plasmidann/bin/python -m pytest tests/test_enrich.py -v`
Expected: FAIL, `ModuleNotFoundError: No module named 'plasmidann.enrich'`

- [ ] **Step 4: Write minimal implementation**

Create `src/plasmidann/enrich.py`:

```python
"""Significance for a genomic-context association.

WHAT WAS WRONG WITH THE PREVIOUS STATISTIC

Context enrichment was observed rate divided by background rate, and nothing else. Three
consequences:

  * No test. A family with three members all sitting beside a relaxase produced the same
    number as a family with three hundred. The output could not distinguish a coincidence
    from a result, and the ranking that selects 1,000 constructs for synthesis read that
    number.

  * The wrong unit. The rate was over family MEMBERS. Members of a protein family are
    homologs, and they sit on plasmids that are frequently near-identical - the same
    clinical plasmid sequenced forty times is forty members and one observation. Counting
    members treats sequencing effort as evidence.

  * No correction. Six hand-picked features needed none. The label vocabulary is open, so
    a family is now tested against every category present in its neighbourhoods, which is
    thousands of tests, where an uncorrected p-value means nothing.

THE UNIT

The unit of observation is the PLASMID, not the family member. A family is counted as
carrying a category if at least one of its members on that plasmid has the category in its
plus or minus three neighbourhood, and the background is the fraction of all plasmids where
any ORF has it. This removes copy-number inflation within a plasmid.

It does NOT remove clonal redundancy between plasmids: forty independent depositions of the
same clinical plasmid remain forty units. That is a separate correction, and the registry
for it already exists in workflow/scripts/clonal_registry.py; applying it here is recorded
as follow-up rather than silently half-done.

THE TEST

Fisher's exact test on the two-by-two table of (family plasmids, rest of the corpus) against
(category present, absent). Exact rather than chi-squared because many categories are rare,
and the expected counts in a rare cell are far below the five that the chi-squared
approximation needs.

Multiple testing is controlled by the Benjamini-Hochberg procedure (Benjamini and Hochberg
1995, J R Stat Soc B 57:289), which controls the false discovery rate. FDR rather than
family-wise error because the purpose is to rank candidates for a 1,000-construct screen,
where a false positive costs a well and a false negative costs a discovery.
"""
from scipy.stats import fisher_exact

from darkorf import status


def fisher_enrichment(k, n, K, N, min_units=2):
    """Test whether a category is over-represented in one family's neighbourhoods.

    `k` of `n` plasmids carrying the family have the category nearby; `K` of `N` plasmids in
    the corpus have it. Returns the rates, the enrichment ratio, the odds ratio, a p-value
    and a status from the pipeline's single status vocabulary.

    A family on fewer than `min_units` plasmids is not tested. One plasmid is an anecdote,
    and returning a p-value for it would let a family seen once rank beside one seen on two
    hundred independent plasmids - which is the error the family-level statistic exists to
    prevent.
    """
    if not 0 <= k <= n:
        raise ValueError(f"family counts are impossible: k={k}, n={n}")
    if not 0 <= K <= N:
        raise ValueError(f"corpus counts are impossible: K={K}, N={N}")

    observed = k / n if n else 0.0
    background = K / N if N else 0.0
    if background == 0:
        # An observation against a zero background is infinitely enriched, which is the
        # honest answer; zero over zero is no signal and no surprise.
        ratio = float("inf") if observed else 1.0
    else:
        ratio = observed / background

    result = {
        "observed_rate": round(observed, 4),
        "background_rate": round(background, 6),
        "enrichment": ratio,
        "odds_ratio": "",
        "p_value": "",
        "status": status.SUCCESS,
    }

    if n < min_units:
        result["status"] = status.TOO_FEW_MEMBERS
        return result
    if N - n <= 0:
        # The family covers the whole corpus, so there is no background to test against.
        result["status"] = status.NOT_APPLICABLE
        return result

    # The family's plasmids against the REST of the corpus. Testing against a background
    # the family contributes to shrinks any real effect, and for a family covering most of
    # the corpus it shrinks it to nothing.
    rest_with = K - k
    rest_without = (N - n) - rest_with
    if rest_with < 0 or rest_without < 0:
        raise ValueError(
            f"corpus counts are inconsistent with the family counts: k={k}, n={n}, "
            f"K={K}, N={N} - the family's plasmids are not a subset of the corpus")

    odds, p_value = fisher_exact([[k, n - k], [rest_with, rest_without]],
                                 alternative="greater")
    result["odds_ratio"] = round(float(odds), 4) if odds == odds else ""
    result["p_value"] = f"{p_value:.6g}"
    return result


def benjamini_hochberg(p_values):
    """False-discovery-rate q-values, in the input order.

    Benjamini and Hochberg 1995, J R Stat Soc B 57:289. FDR rather than a family-wise
    correction because these values rank candidates for a 1,000-construct screen: a false
    positive costs one well, a false negative costs a discovery.

    Entries that are not numbers - a category that could not be tested - pass through
    unchanged. Treating an empty value as zero would give an untested category the
    strongest q-value in the table.

    The order is preserved because the caller joins the result back to its rows by
    position; returning a sorted list would attach every q-value to the wrong category.
    """
    indexed = []
    for i, value in enumerate(p_values):
        try:
            indexed.append((float(value), i))
        except (TypeError, ValueError):
            continue

    out = list(p_values)
    m = len(indexed)
    if not m:
        return out

    indexed.sort()
    # Walk from the largest p-value down, keeping the running minimum. This enforces the
    # monotonicity the procedure requires: the q-value of the i-th smallest p is the
    # minimum over j >= i of (m / j) * p_j.
    running = 1.0
    for rank in range(m, 0, -1):
        p, i = indexed[rank - 1]
        running = min(running, p * m / rank)
        out[i] = round(min(running, 1.0), 6)
    return out
```

- [ ] **Step 5: Run tests to verify they pass**

Run: `envs/plasmidann/bin/python -m pytest tests/test_enrich.py -v`
Expected: PASS, 11 tests

- [ ] **Step 6: Record the parameters**

Append to `docs/PARAMETER_PROVENANCE.md`:

```markdown
## Context enrichment significance

| parameter | value | source |
|---|---|---|
| test | Fisher's exact, one-sided (`alternative="greater"`) | Exact rather than chi-squared because many categories are rare and expected cell counts fall well below the 5 the chi-squared approximation requires. One-sided because the hypothesis is over-representation; a category a dark family avoids is not a screening hypothesis. |
| multiple-testing correction | Benjamini-Hochberg FDR | Benjamini and Hochberg 1995, J R Stat Soc B 57:289. FDR rather than family-wise error because these values rank candidates for a 1,000-construct screen: a false positive costs a well, a false negative costs a discovery. |
| unit of observation | plasmid | A family's members are homologs and its plasmids are frequently near-identical. Counting members makes sequencing effort look like evidence. Clonal redundancy BETWEEN plasmids is not corrected here; `workflow/scripts/clonal_registry.py` holds the registry for it and the correction is outstanding. |
| `min_units` | 2 | One plasmid is an anecdote. Below two units no test is performed and the status is `TOO_FEW_MEMBERS`, so an untested family cannot rank beside a tested one. |
```

- [ ] **Step 7: Commit**

```bash
git add src/plasmidann/enrich.py tests/test_enrich.py \
        workflow/envs/plasmidann.yaml docs/PARAMETER_PROVENANCE.md
git commit -m "$(cat <<'EOF'
feat(context): test context enrichment, per plasmid, with FDR control

The previous statistic was observed rate over background rate with no test,
so a family with three members all beside a relaxase scored the same as one
with three hundred, and the construct ranking read that number.

Three fixes: Fisher's exact test, because many categories are rare and the
chi-squared approximation does not hold there; the plasmid as the unit
rather than the family member, because members are homologs and counting
them makes sequencing effort look like evidence; and Benjamini-Hochberg
correction, because the open label vocabulary means thousands of tests per
family rather than six.

Clonal redundancy between plasmids is still uncorrected and is recorded as
outstanding in docs/PARAMETER_PROVENANCE.md rather than half-applied.

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>
EOF
)"
```

---

### Task 9: Retire the curated list

The list can only be deleted once nothing needs it. `control_recall` moves out first, then the context stage is rewritten onto the label vocabulary, then the file goes.

**Files:**
- Create: `src/plasmidann/controls.py`, `tests/test_controls.py`
- Delete: `src/plasmidann/backbone.py`, `tests/test_backbone.py`
- Modify: `workflow/scripts/quality_gate.py:25-35`
- Modify: `workflow/scripts/context_features.py` (whole file)
- Modify: `workflow/rules/s5_s9_targets.smk` (rule `context_features` inputs and params)
- Modify: `workflow/scripts/prioritise.py:150-170`, `config/targets.yaml:196-203`
- Modify: `workflow/scripts/annotation_report.py:84-85`, `tests/test_context.py`
- Test: `tests/test_scripts_smoke.py` (replace the two context tests)

**Interfaces:**
- Consumes: `plasmidann.labels`, `plasmidann.categories`, `plasmidann.enrich` from Tasks 5, 7, 8; `results/s4c/protein_labels.tsv` from Task 6.
- Produces:
  - `plasmidann.controls.control_recall(rows) -> float`, identical behaviour to the current function.
  - `results/s8/family_context.tsv` in LONG format, columns:
    `family_id, category, subcategories, n_units, n_units_with_category, observed_rate, background_rate, enrichment, odds_ratio, p_value, q_value, status`.
    The test runs on `category` alone; `subcategories` is the sorted distinct set observed
    for that (family, category), carried so a finer grouping can be re-cut without a rerun.
  - `results/s8/context_background.tsv`, columns: `category, n_units_with_category, n_units, background_rate`.

- [ ] **Step 1: Move `control_recall`, with its tests**

Create `src/plasmidann/controls.py`:

```python
"""The positive control for the annotation cascade (success criterion SC2).

Following ECLIPSE (Bioinformatics 42:8), which recovered 99.2-100% of 246 virulence, 42
AMR and 75 essential genes as annotated. Anything KNOWN that emerges dark is a recall
failure in the cascade, and the run should stop rather than produce a target list built on
a broken annotation step. Its absence from v1 was the single point of unanimous reviewer
criticism.

The control set is independent by construction: reviewed Swiss-Prot proteins of known
function, spiked into the query set at S2c and carried through every tier exactly as a real
protein is. An earlier version drew the control from proteins the cascade itself had
labelled, which cannot detect the failure that matters, because a protein the cascade
MISSED never enters a self-drawn control set.
"""


def control_recall(rows):
    """Fraction of a known-function control set that the cascade called FUNCTIONAL.

    The gate for success criterion SC2. Raises on an empty control set rather than
    returning 1.0 or 0.0: an empty control silently passing is exactly how a gate stops
    being a gate, and it is the failure mode that would reproduce v1's missing control.
    """
    if not rows:
        raise ValueError(
            "the positive control set is empty - a gate that cannot fail is not a gate. "
            "Check that the control ids resolved against the annotation table.")
    n = sum(1 for r in rows if r.get("functional_class") == "FUNCTIONAL")
    return round(n / len(rows), 4)
```

Create `tests/test_controls.py`, moving the three control tests from `tests/test_backbone.py` and changing only the import:

```python
"""The positive control for the annotation cascade (SC2).

Moved out of tests/test_backbone.py when the curated backbone family list was deleted. The
control is unrelated to that list and always was: it draws on reviewed Swiss-Prot proteins
spiked into the query set, which is what makes it able to detect a protein the cascade
missed.
"""
import pytest

from plasmidann.controls import control_recall


def test_control_recall_is_the_fraction_correctly_called_functional():
    """ECLIPSE recovered 99.2-100% of 246 virulence, 42 AMR and 75 essential genes as
    annotated. Anything known that comes out dark is a recall failure in the cascade."""
    rows = [{"functional_class": "FUNCTIONAL"}] * 99 + [{"functional_class": "NONE"}]

    assert control_recall(rows) == 0.99


def test_a_control_set_where_everything_is_dark_scores_zero():
    rows = [{"functional_class": "UNCHARACTERIZED_HOMOLOG"}] * 10

    assert control_recall(rows) == 0.0


def test_an_empty_control_set_is_an_error_not_a_pass():
    """Silently passing on an empty control is how a gate stops being a gate."""
    with pytest.raises(ValueError):
        control_recall([])
```

In `workflow/scripts/quality_gate.py`, change the import on line 35 to:

```python
from plasmidann.controls import control_recall
```

and replace the module-docstring sentence beginning "The curated family table survives in plasmidann.backbone" with:

```
The curated family table is gone. It covered 15 of roughly 30 T4SS families in Pfam, named
no MobB and no MobD, and nine of its 73 names did not exist in Pfam-A at all - entries that
had never once matched anything. Functional labels now come from the tools themselves, in
results/s4c/protein_labels.tsv, and are grouped into categories from the observed
vocabulary rather than by hand.
```

Run: `envs/plasmidann/bin/python -m pytest tests/test_controls.py tests/test_scripts_smoke.py -k "control or quality" -q`
Expected: PASS

- [ ] **Step 2: Commit the move**

```bash
git add src/plasmidann/controls.py tests/test_controls.py workflow/scripts/quality_gate.py
git commit -m "$(cat <<'EOF'
refactor(controls): move control_recall out of the backbone module

The positive control is unrelated to the curated family list and always was:
it draws on reviewed Swiss-Prot proteins spiked into the query set, which is
what lets it detect a protein the cascade missed. Moving it first makes the
list deletable.

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>
EOF
)"
```

- [ ] **Step 3: Write the failing tests for the new context stage**

Find the two existing context tests with `grep -n "backbone_adjacent\|ta_candidate" tests/test_scripts_smoke.py` and replace them with:

```python
# --- S8c: open-vocabulary context, in long format -------------------------------------

def test_context_features_records_the_labels_the_tools_produced(fixture_dir):
    """The curated 73-name list is gone, so a neighbour's context is whatever the tools
    called it. A dark ORF beside a relaxase must come out associated with the relaxase's
    own label, not with a hand-assigned 'backbone' class - the class is derived later from
    the observed vocabulary, and the raw label has to be in the table for that to be
    possible."""
    ann = fixture_dir / "plasmid_annotation.tsv"
    write_tsv(ann, ["orf_id", "plasmid_id", "start", "end", "strand", "partial",
                    "spans_origin", "annot_label", "functional_class", "annot_tier",
                    "artefact_flag"],
              # A dark ORF beside a relaxase on each of two plasmids, and a third plasmid
              # with no relaxase so the background is not 1.0.
              [["p1|1", "p1", 100, 400, "+", 0, 0, "MobA_MobL", "FUNCTIONAL", "T1", 0],
               ["p1|2", "p1", 500, 700, "+", 0, 0, "", "NONE", "", 0],
               ["p2|1", "p2", 100, 400, "+", 0, 0, "MobA_MobL", "FUNCTIONAL", "T1", 0],
               ["p2|2", "p2", 500, 700, "+", 0, 0, "", "NONE", "", 0],
               ["p3|1", "p3", 100, 400, "+", 0, 0, "Aspartokinase", "FUNCTIONAL", "T1", 0],
               ["p3|2", "p3", 500, 700, "+", 0, 0, "", "NONE", "", 0]])

    labels_tsv = fixture_dir / "protein_labels.tsv"
    write_tsv(labels_tsv, ["protein_id", "source", "tier", "kind", "label", "accession",
                           "evidence_evalue", "evidence_coverage", "database",
                           "database_version"],
              [["s_mob", "pfam", "T1", "pfam_family", "MobA_MobL", "PF03389.20",
                "1e-40", "0.9", "Pfam-A", "38.2"],
               ["s_ask", "pfam", "T1", "pfam_family", "Aspartokinase", "PF00696.1",
                "1e-40", "0.9", "Pfam-A", "38.2"]])

    mapping = fixture_dir / "protein_map.tsv"
    mapping.write_text("s_mob\tp1|1,p2|1\ns_ask\tp3|1\ns_dark\tp1|2,p2|2,p3|2\n")

    families = fixture_dir / "dark_families.tsv"
    write_tsv(families, ["family_id", "representative", "n_members", "n_orfs",
                         "n_plasmids", "n_mob_clusters", "family_class", "members"],
              [["F0000001", "s_dark", 1, 3, 3, 2, "FAMILY", "s_dark"]])

    defence = fixture_dir / "defence_systems.tsv"
    write_tsv(defence, ["orf_id", "plasmid_id", "gembase_id", "system", "system_id",
                        "component", "hit_evalue", "hit_status", "sys_wholeness",
                        "hit_gene_ref", "hit_profile_cov"], [])
    integrons = fixture_dir / "integrons.tsv"
    write_tsv(integrons, ["plasmid_id", "start", "end", "integron_type", "annotation",
                          "type_elt"], [])

    out = fixture_dir / "family_context.tsv"
    run_script("context_features.py", FakeSnakemake(
        input={"annotation": str(ann), "families": str(families), "map": str(mapping),
               "defence": str(defence), "integrons": [str(integrons)],
               "labels": str(labels_tsv)},
        output={"families": str(out),
                "background": str(fixture_dir / "context_background.tsv")},
        params={"context": {"max_operon_gap": 100, "neighbourhood_window": 3,
                            "min_context_conservation": 0.50, "max_q_value": 0.05},
                "categories": None}))

    rows = {r["category"]: r for r in read_tsv(out) if r["family_id"] == "F0000001"}

    # With no category rules configured, the category is the label qualified by its kind.
    assert "pfam_family:MobA_MobL" in rows, (
        f"the relaxase label is not in the context table; got {sorted(rows)}")
    row = rows["pfam_family:MobA_MobL"]
    assert row["n_units"] == "3", "the unit is the plasmid, so n_units should be 3"
    assert row["n_units_with_category"] == "2"
    assert float(row["enrichment"]) > 1.0
    assert row["q_value"], "no q-value written, so nothing corrected for multiple testing"
    # The column exists even with no rules configured, so a later grouping can fill it
    # without changing the table's shape.
    assert "subcategories" in row


def test_context_features_writes_no_hand_assigned_category(fixture_dir):
    """The standing constraint at the stage that used to break it. 'backbone_adjacent' and
    'ta_candidate' were hand-assigned classes built on a 73-name list, and neither may
    reappear as a category name."""
    ann = fixture_dir / "plasmid_annotation.tsv"
    write_tsv(ann, ["orf_id", "plasmid_id", "start", "end", "strand", "partial",
                    "spans_origin", "annot_label", "functional_class", "annot_tier",
                    "artefact_flag"],
              [["p1|1", "p1", 100, 400, "+", 0, 0, "RelB", "FUNCTIONAL", "T1", 0],
               ["p1|2", "p1", 420, 600, "+", 0, 0, "", "NONE", "", 0]])
    labels_tsv = fixture_dir / "protein_labels.tsv"
    write_tsv(labels_tsv, ["protein_id", "source", "tier", "kind", "label", "accession",
                           "evidence_evalue", "evidence_coverage", "database",
                           "database_version"],
              [["s_rel", "pfam", "T1", "pfam_family", "RelB", "PF04221.1", "1e-30",
                "0.9", "Pfam-A", "38.2"]])
    mapping = fixture_dir / "protein_map.tsv"
    mapping.write_text("s_rel\tp1|1\ns_dark\tp1|2\n")
    families = fixture_dir / "dark_families.tsv"
    write_tsv(families, ["family_id", "representative", "n_members", "n_orfs",
                         "n_plasmids", "n_mob_clusters", "family_class", "members"],
              [["F0000001", "s_dark", 1, 1, 1, 1, "FAMILY", "s_dark"]])
    defence = fixture_dir / "defence_systems.tsv"
    write_tsv(defence, ["orf_id", "plasmid_id", "gembase_id", "system", "system_id",
                        "component", "hit_evalue", "hit_status", "sys_wholeness",
                        "hit_gene_ref", "hit_profile_cov"], [])
    integrons = fixture_dir / "integrons.tsv"
    write_tsv(integrons, ["plasmid_id", "start", "end", "integron_type", "annotation",
                          "type_elt"], [])

    out = fixture_dir / "family_context.tsv"
    run_script("context_features.py", FakeSnakemake(
        input={"annotation": str(ann), "families": str(families), "map": str(mapping),
               "defence": str(defence), "integrons": [str(integrons)],
               "labels": str(labels_tsv)},
        output={"families": str(out),
                "background": str(fixture_dir / "context_background.tsv")},
        params={"context": {"max_operon_gap": 100, "neighbourhood_window": 3,
                            "min_context_conservation": 0.50, "max_q_value": 0.05},
                "categories": None}))

    forbidden = {"backbone_adjacent", "ta_candidate", "replication", "mobilisation",
                 "conjugation", "toxin_antitoxin"}
    seen = {r["category"] for r in read_tsv(out)}

    assert not (seen & forbidden), f"hand-assigned category reappeared: {seen & forbidden}"
```

- [ ] **Step 4: Run tests to verify they fail**

Run: `envs/plasmidann/bin/python -m pytest tests/test_scripts_smoke.py -k "context" -v`
Expected: FAIL, `KeyError: 'labels'` and `KeyError: 'category'`

- [ ] **Step 5: Rewrite the context stage**

Replace `workflow/scripts/context_features.py` in full:

```python
"""S8c: genomic context per ORF, aggregated to families and tested against a background.

WHAT CHANGED AND WHY

This stage used to describe a dark ORF's neighbourhood with six hand-picked features, two
of which - backbone_adjacent and ta_candidate - were decided by a hand-written list of 73
Pfam family names. That list covered 15 of roughly 30 T4SS families in Pfam, named no MobB
and no MobD, and nine of its names did not exist in Pfam-A at all, so those entries had
never once matched anything and no output could have revealed it.

A neighbourhood is now described by the labels the tools actually produced for the
neighbours, read from results/s4c/protein_labels.tsv. Those labels are grouped into
categories through plasmidann.categories, which with no rules configured makes each
(kind, label) its own category. The biological grouping - replication, mobilisation,
conjugation - is derived later from the observed vocabulary and applied here as a
configuration change.

THE TABLE IS LONG, NOT WIDE

Six features fitted in six column pairs. An open vocabulary does not: there are 30,134 Pfam
families alone. One row per (family, category) is the only shape that holds it, and it is
the shape the later grouping reads naturally.

THE STATISTICAL TRAP THIS STAGE EXISTS TO AVOID

On a 5 kb cryptic plasmid carrying six genes, a plus or minus three neighbourhood IS the
entire plasmid. Everything co-occurs with everything, and a raw co-occurrence frequency
would rank the smallest plasmids as the most informative when they are the least - which
would be catastrophic here, because small cryptic plasmids are a stratum of interest. So
every association is reported as enrichment over a corpus-wide background, and now also
with a Fisher exact p-value and a Benjamini-Hochberg q-value.

THE UNIT IS THE PLASMID

Not the family member. Members of a family are homologs on plasmids that are frequently
near-identical, so counting members makes sequencing effort look like evidence. Clonal
redundancy between distinct plasmids remains uncorrected; see
docs/PARAMETER_PROVENANCE.md.
"""
import _ctx  # noqa: F401
import collections
import csv

from plasmidann import categories, enrich
from plasmidann.context import directons, neighbourhood, overlapping_islands

cfg = snakemake.params.context
cmap = categories.CategoryMap(categories.load_rules(snakemake.params.categories))

# ------------------------------------------------------------------------------------
# Every ORF on every plasmid, with whatever the cascade named it.
# ------------------------------------------------------------------------------------
by_plasmid = collections.defaultdict(list)
class_of = {}
with open(snakemake.input.annotation, newline="") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        by_plasmid[r["plasmid_id"]].append({
            "orf_id": r["orf_id"], "start": int(r["start"]), "end": int(r["end"]),
            "strand": 1 if r["strand"] in ("1", "+") else -1})
        class_of[r["orf_id"]] = r.get("functional_class") or "NONE"

# ------------------------------------------------------------------------------------
# Which ORFs correspond to which unique protein. BOTH directions, once.
#
# The reverse direction matters: the per-family loop below previously rescanned the whole
# 9.3M-entry forward map per family, which is O(families x ORFs) - measured at 0.37-0.39 s
# per family, projecting to 10 hours to 3.6 days single-core for ~50k families.
# ------------------------------------------------------------------------------------
seq_of_orf = {}
orfs_of_seq = collections.defaultdict(list)
with open(snakemake.input.map) as fh:
    for line in fh:
        sid, ids = line.rstrip("\n").split("\t")
        for oid in ids.split(","):
            seq_of_orf[oid] = sid
            orfs_of_seq[sid].append(oid)

# ------------------------------------------------------------------------------------
# The categories of every unique protein, from the labels the tools produced.
#
# Read per protein, not per ORF: a protein identical on forty plasmids was labelled once.
# ------------------------------------------------------------------------------------
# Each entry is a (category, subcategory) pair. The category is the flat axis the
# significance test runs on; the subcategory travels with it and is never tested.
categories_of_seq = collections.defaultdict(set)
n_label_rows = 0
with open(snakemake.input.labels, newline="") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        n_label_rows += 1
        categories_of_seq[r["protein_id"]] |= set(
            cmap.categories_for(r["kind"], r["label"]))

# ------------------------------------------------------------------------------------
# Islands: defence systems and integron cassette arrays, as intervals per plasmid.
#
# These are system-level calls rather than per-protein labels, so they stay intervals: a
# dark ORF INSIDE a defence system is a different statement from one merely beside a
# defence component.
# ------------------------------------------------------------------------------------
islands = collections.defaultdict(list)
defence_orfs = {}
with open(snakemake.input.defence, newline="") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        if r.get("orf_id"):
            # The system's own name from MacSyFinder's model, not a class of our invention.
            defence_orfs[r["orf_id"]] = (r.get("system") or "").rsplit("/", 1)[-1]

for f in snakemake.input.integrons:
    with open(f, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            try:
                islands[r["plasmid_id"]].append(
                    {"name": "integron", "start": int(r["start"]), "end": int(r["end"]),
                     "detail": r.get("integron_type", "")})
            except (ValueError, KeyError):
                continue

for pid, genes in by_plasmid.items():
    for g in genes:
        system = defence_orfs.get(g["orf_id"])
        if system:
            islands[pid].append({"name": f"defence_system:{system}",
                                 "start": g["start"], "end": g["end"],
                                 "detail": system})

# ------------------------------------------------------------------------------------
# Per-ORF context: the categories of its neighbours, its directon partners, and the
# islands it sits inside.
# ------------------------------------------------------------------------------------
context_of = {}
for pid, genes in by_plasmid.items():
    units = directons(genes, max_gap=cfg["max_operon_gap"])
    unit_of = {oid: i for i, unit in enumerate(units) for oid in unit}
    plasmid_islands = islands.get(pid, [])

    for g in genes:
        oid = g["orf_id"]
        ctx = set()

        # EVERY island, not the first one found. Defence intervals are appended after
        # integron intervals, so returning the first made a dark ORF inside a defence
        # system that also sat in a cassette array label as 'integron' only.
        for island in overlapping_islands(g, plasmid_islands):
            ctx.add((island["name"], ""))

        # What the neighbours are, in the tools' own words. No hand-assigned class: a
        # neighbour labelled MobA_MobL contributes the category of MobA_MobL, which with no
        # rules configured is 'pfam_family:MobA_MobL'.
        for n in neighbourhood(genes, oid, window=cfg["neighbourhood_window"]):
            ctx |= categories_of_seq.get(seq_of_orf.get(n, ""), set())
            if class_of.get(n) == "FUNCTIONAL":
                ctx.add(("annotated_neighbour", ""))

        # Directon membership with at least one annotated partner: the dark ORF is
        # predicted to be co-transcribed with something we understand, which is a stronger
        # claim than adjacency.
        unit = units[unit_of[oid]] if oid in unit_of else [oid]
        partners = [p for p in unit if p != oid]
        if partners and any(class_of.get(p) == "FUNCTIONAL" for p in partners):
            ctx.add(("operon_with_annotated", ""))
        # The geometry of a tight two-gene unit, recorded as geometry. It used to be read
        # as a toxin-antitoxin call whenever the partner was on the 73-name list; whether
        # the partner's category makes it one is a question for the grouping, not for this
        # stage, and the partner's categories are already in ctx for that to be answered.
        if len(unit) == 2:
            ctx.add(("two_gene_operon", ""))

        context_of[oid] = ctx

# ------------------------------------------------------------------------------------
# Background, per PLASMID. A plasmid counts for a category if any of its ORFs has that
# category in context. Counting ORFs instead would let one gene-dense plasmid dominate.
# ------------------------------------------------------------------------------------
plasmids_with = collections.Counter()
for pid, genes in by_plasmid.items():
    present = set()
    for g in genes:
        present |= context_of.get(g["orf_id"], set())
    # Counted on the CATEGORY alone. Counting (category, subcategory) pairs would make the
    # background depend on how finely the categories happen to be subdivided, so adding a
    # subcategory would silently change every enrichment value in the table.
    for category in {c for c, _ in present}:
        plasmids_with[category] += 1

n_plasmids = len(by_plasmid)

with open(snakemake.output.background, "w", newline="") as out:
    w = csv.writer(out, delimiter="\t")
    w.writerow(["category", "n_units_with_category", "n_units", "background_rate"])
    for category, count in sorted(plasmids_with.items()):
        w.writerow([category, count, n_plasmids,
                    round(count / n_plasmids, 6) if n_plasmids else 0.0])

# ------------------------------------------------------------------------------------
# Per family: one row per category present in any member's context, tested against the
# background, with q-values corrected across the categories tested for that family.
# ------------------------------------------------------------------------------------
cols = ["family_id", "category", "subcategories", "n_units", "n_units_with_category",
        "observed_rate", "background_rate", "enrichment", "odds_ratio", "p_value",
        "q_value", "status"]

n_rows = n_families = 0
with open(snakemake.output.families, "w", newline="") as out:
    w = csv.DictWriter(out, fieldnames=cols, delimiter="\t")
    w.writeheader()
    with open(snakemake.input.families, newline="") as fh:
        for fam in csv.DictReader(fh, delimiter="\t"):
            n_families += 1
            seq_members = fam["members"].split(",")
            orfs = [o for m in seq_members for o in orfs_of_seq.get(m, ())]

            # The unit is the plasmid. orf_id is '<plasmid_id>|<ordinal>'.
            per_plasmid = collections.defaultdict(set)
            for oid in orfs:
                per_plasmid[oid.rsplit("|", 1)[0]] |= context_of.get(oid, set())
            n_units = len(per_plasmid)
            if not n_units:
                continue

            # Counts on the flat category; the subcategories seen for it are collected
            # beside the count and reported, never tested.
            counts = collections.Counter()
            subcategories = collections.defaultdict(set)
            for present in per_plasmid.values():
                for category, subcategory in present:
                    if subcategory:
                        subcategories[category].add(subcategory)
                for category in {c for c, _ in present}:
                    counts[category] += 1

            tested = []
            for category, k in sorted(counts.items()):
                tested.append((category, k, enrich.fisher_enrichment(
                    k=k, n=n_units,
                    K=plasmids_with.get(category, 0), N=n_plasmids)))

            # Corrected across the categories tested for THIS family. A family tested
            # against four hundred categories and one tested against three are not
            # comparable without it, and the open vocabulary makes that the normal case.
            q_values = enrich.benjamini_hochberg([r["p_value"] for _, _, r in tested])

            for (category, k, result), q in zip(tested, q_values):
                n_rows += 1
                w.writerow({
                    "family_id": fam["family_id"], "category": category,
                    "subcategories": ",".join(sorted(subcategories.get(category, ()))),
                    "n_units": n_units, "n_units_with_category": k,
                    "observed_rate": result["observed_rate"],
                    "background_rate": result["background_rate"],
                    "enrichment": result["enrichment"],
                    "odds_ratio": result["odds_ratio"],
                    "p_value": result["p_value"], "q_value": q,
                    "status": result["status"]})

grouping = "per label (no category rules configured)" if cmap.is_identity else "per category"
print(f"context: {n_label_rows} label rows, {len(plasmids_with)} categories over "
      f"{n_plasmids} plasmids, {n_rows} association rows for {n_families} families, "
      f"grouping is {grouping}")
```

- [ ] **Step 6: Rewire the rule**

In `workflow/rules/s5_s9_targets.smk`, in rule `context_features`, replace the `input:` and `params:` blocks with:

```python
    input:
        annotation=f"{OUT}/s4/plasmid_annotation.tsv",
        families=f"{OUT}/s6/dark_families.tsv",
        map=f"{OUT}/s2/protein_map.tsv",
        defence=f"{OUT}/s8/defence_systems.tsv",
        integrons=expand(f"{OUT}/s8/integrons/{{shard}}.tsv", shard=SHARDS),
        labels=f"{OUT}/s4c/protein_labels.tsv",
    output:
        families=f"{OUT}/s8/family_context.tsv",
        background=f"{OUT}/s8/context_background.tsv",
    params:
        context=targets["context"],
        # None means no grouping: every (kind, label) is its own category. The biological
        # grouping is derived from results/s4c/protein_labels.tsv after a full annotation
        # run and enabled by pointing this at config/label_categories.yaml.
        categories=config.get("label_categories"),
```

- [ ] **Step 7: Rewire the downstream consumers**

`workflow/scripts/prioritise.py` reads a wide context table that no longer exists. Replace the block that loads `snakemake.input.context` with a reduction over the long table:

```python
# The context table is LONG: one row per (family, category). The strongest association per
# family is selected here rather than written by S8c, because "strongest" depends on the
# q-value threshold, which is a prioritisation parameter and is swept.
context = {}
with open(snakemake.input.context, newline="") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        if r["status"] != "SUCCESS":
            continue
        try:
            q = float(r["q_value"])
        except (TypeError, ValueError):
            continue
        if q > cfg["max_q_value"]:
            continue
        if float(r["observed_rate"]) < cfg["min_context_conservation"]:
            continue
        best = context.get(r["family_id"])
        if best is None or q < best["q_value"]:
            context[r["family_id"]] = {
                "top_hypothesis": r["category"], "q_value": q,
                "top_conservation": float(r["observed_rate"]),
                "top_enrichment": r["enrichment"]}
```

Replace the `ta_candidate` stratum rule (currently `prioritise.py:164-165`) with:

```python
    # The toxin_or_ta_adjacent stratum was filled by a 'ta_candidate' hypothesis that fired
    # when a two-gene partner appeared on a 73-name curated list. The list is gone and the
    # replacement - a toxin_antitoxin CATEGORY - does not exist until the grouping is
    # derived from results/s4c/protein_labels.tsv. Until then the stratum is unfillable and
    # is declared suspended in config/targets.yaml rather than being filled by a rule that
    # no longer has a definition behind it.
```

In `config/targets.yaml`, replace the `portfolio:` block with:

```yaml
portfolio:
  total: 1000
  strata:
    membrane_or_secreted: 200      # localisation - informative for every member
    defence_island: 175            # localisation; full readout needs a phage arm
    integron_cassette: 175         # see below - real genes by construction
    small_cationic_peptide: 125    # lysis, membrane disruption
    nucleic_acid_binding: 75       # nucleoid association, filamentation
    novel_fold: 250                # open discovery arm, plus the 175 released below
  #
  # toxin_or_ta_adjacent (175) is SUSPENDED, not removed. It was filled by a hypothesis
  # that fired when a dark ORF's two-gene partner appeared on a hand-written list of 73
  # Pfam family names. That list is deleted: it covered 15 of roughly 30 T4SS families in
  # Pfam, named no MobB and no MobD, and nine of its names did not exist in Pfam-A at all.
  # Its replacement is a toxin_antitoxin CATEGORY derived from the observed label
  # vocabulary in results/s4c/protein_labels.tsv, which does not exist until a full
  # annotation run has been made. The quota is parked in novel_fold rather than left
  # unfilled, so the portfolio still totals 1,000 and the suspension is visible here
  # rather than appearing as a short order at the bench.
```

In `workflow/scripts/annotation_report.py:84-85`, the `cons_defence`, `cons_integron`, `cons_backbone_adjacent`, `cons_annotated_neighbour`, `cons_operon_with_annotated`, `cons_ta_candidate` column names no longer exist. Replace that column list with the long-table equivalent the report needs: the family's top category, its q-value and its enrichment, read from `family_context.tsv` with the same reduction as `prioritise.py`. Keep the report's other columns unchanged.

- [ ] **Step 8: Delete the list**

```bash
git rm src/plasmidann/backbone.py tests/test_backbone.py
```

Confirm nothing still imports it:
```bash
grep -rn "plasmidann.backbone\|from plasmidann import backbone\|is_backbone\|backbone_names_for_role\|BACKBONE_FAMILIES\|backbone_adjacent\|ta_candidate" \
  --include=*.py --include=*.smk --include=*.yaml src workflow tests config
```
Expected: no output. Any hit is a consumer that must be rewired before the deletion is complete.

- [ ] **Step 9: Run the whole suite**

Run: `PATH=envs/plasmidann/bin:$PATH envs/plasmidann/bin/python -m pytest -q`
Expected: PASS. Failures in `tests/test_context.py` are expected where it asserts the wide-table shape; update those assertions to the long shape rather than deleting the tests.

- [ ] **Step 10: Update the documentation**

In `docs/PIPELINE_CODE.md`, replace the section describing the backbone list with three new sections:

```markdown
## 6. `src/plasmidann/labels.py` - the tool-derived label vocabulary

Every functional label, verbatim from the tool that produced it, tagged with the KIND of
statement it is: `pfam_family`, `pfam_description`, `pfam_clan`, `swissprot_product`,
`pgap_product`, `gene_symbol`, `cog_category`, `cog_id`, `eggnog_pfam`, `go`, `ec`,
`kegg_ko`, `macsy_system`, `macsy_component`, `integron_element`, `integron_type`.

This replaces a hand-written list of 73 Pfam family names that assigned each one a
biological role. The list could not work, and the measurements say why: Pfam holds roughly
30 T4SS families and the list named 15; it named no MobB and no MobD; **nine of its 73
names do not exist in Pfam-A at all**, so those entries had never once matched anything,
through a whole version, invisibly; and pfam2go covers 4 of its 16 replication families, 1
of 16 conjugation families and 0 of 11 mobilisation families, so the gap is in the curated
mappings too.

No role is assigned in this module, and `tests/test_labels.py` asserts that no label kind
is a role name, so the list cannot grow back.

## 7. `src/plasmidann/categories.py` - the grouping seam

Groups `(kind, label)` pairs into functional categories from `config/label_categories.yaml`.
That file **ships empty**: the categories are derived from the vocabulary observed in
`results/s4c/protein_labels.tsv` after a full annotation run, which is the only way they can
cover the scope. With no rules, every `(kind, label)` is its own category, so the enrichment
machinery runs from the first run and applying the grouping later is a configuration change.

A label matching no rule falls back to being its own category and is never dropped. A
vocabulary that shrank silently as rules were added is how the previous list's nine dead
entries survived.

## 8. `src/plasmidann/enrich.py` - significance

Fisher's exact test, one-sided, on the two-by-two table of (family plasmids, rest of corpus)
against (category present, absent), with Benjamini-Hochberg FDR control across the
categories tested per family.

The unit is the **plasmid**, not the family member. Members of a family are homologs on
plasmids that are frequently near-identical, so counting members makes sequencing effort
look like evidence. Clonal redundancy between distinct plasmids is **not** corrected;
`workflow/scripts/clonal_registry.py` holds the registry and the correction is outstanding.
```

- [ ] **Step 11: Commit**

```bash
git add -A src/plasmidann workflow/scripts/context_features.py \
        workflow/scripts/prioritise.py workflow/scripts/annotation_report.py \
        workflow/rules/s5_s9_targets.smk config/targets.yaml tests docs/PIPELINE_CODE.md
git commit -m "$(cat <<'EOF'
feat(context)!: describe neighbourhoods with tool labels, delete the curated list

src/plasmidann/backbone.py is deleted. Its 73 hand-written Pfam family names
covered 15 of roughly 30 T4SS families in Pfam, named no MobB and no MobD,
and nine of them do not exist in Pfam-A at all - entries that had never once
matched anything, through a whole version, invisibly. The role assignments
had no source.

A neighbourhood is now described by the labels the tools produced for the
neighbours, read from results/s4c/protein_labels.tsv and grouped through
plasmidann.categories, which with no rules configured makes each (kind,
label) its own category. family_context.tsv becomes long: one row per
(family, category), because the vocabulary is open.

BREAKING: family_context.tsv changes from wide cons_*/enrich_* columns to
long format, and the toxin_or_ta_adjacent stratum is suspended with its
quota parked in novel_fold. The stratum was filled by a hypothesis that
fired on the deleted list; its replacement is a toxin_antitoxin category
derived from the observed vocabulary, which needs a full annotation run.

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>
EOF
)"
```

---

### Task 10: The category builder

Run after the first full annotation run. Turns the observed vocabulary into a **draft** `label_categories.yaml` with counts beside every candidate, so the review is bounded work on real data rather than recollection.

**Files:**
- Create: `src/plasmidann/normalise.py`, `tests/test_normalise.py`
- Create: `workflow/scripts/build_categories.py`
- Create: `config/category_mining.yaml`
- Test: `tests/test_scripts_smoke.py` (add one test)

**Interfaces:**
- Consumes: `results/s4c/protein_labels.tsv` from Task 6; `plasmidann.pfam_meta` from Task 1.
- Produces:
  - `normalise.product_name(text) -> str`
  - `build_categories.py` writes `results/s4c/category_draft.yaml` and `results/s4c/label_frequency.tsv`.

**Why a draft and not the final file.** Measured on Pfam-A 38.2: the description field matches `/replicat/` in **67** families where the deleted list named 16, `/conjug/` in **42** where it named 15, and `/transposas/` in **76** where it named 15. So mining gets the scope a hand list never had. In the same measurement `/toxin/` matches **317** families, of which `ABC_toxin_N` is an insect toxin and `ADPRTs_Tse2` is a T6SS effector, neither being a toxin-antitoxin system. So mining also produces false positives that need a person. The pattern therefore **generates candidates and never matches at runtime**: what ships in `config/label_categories.yaml` is the reviewed, enumerated list, and `categories.py` keeps its exact matching.

- [ ] **Step 1: Write the failing test for normalisation**

Create `tests/test_normalise.py`:

```python
"""Collapsing free-text product names to a canonical form.

Pfam family names and COG identifiers are controlled and finite, so they can be enumerated.
Product names from nr and Swiss-Prot are neither: the same protein is written 'Plasmid
replication initiator protein RepA', 'replication initiation protein' and 'putative
replication initiator protein' by different submitters. Left as written, one concept
becomes dozens of categories, each too rare to test.

Normalisation is what makes the frequency ranking meaningful: once the surface forms
collapse, the head of the distribution is a bounded list a person can review.
"""
from plasmidann import normalise


def test_the_same_concept_written_three_ways_collapses_to_one_string():
    forms = ["Plasmid replication initiator protein",
             "plasmid replication initiator protein",
             "putative plasmid replication initiator protein"]

    assert len({normalise.product_name(f) for f in forms}) == 1


def test_hedging_words_are_stripped():
    """'putative', 'probable' and 'predicted' describe the submitter's confidence, not the
    protein. Keeping them would split every product into a confident and a hedged form."""
    assert normalise.product_name("putative relaxase") == "relaxase"
    assert normalise.product_name("probable relaxase") == "relaxase"
    assert normalise.product_name("predicted relaxase") == "relaxase"


def test_the_generic_protein_suffixes_are_stripped():
    """'MobA/MobL family protein' and 'MobA/MobL protein' are one product. The suffixes
    carry no functional information and are applied inconsistently by submitters."""
    assert normalise.product_name("MobA/MobL family protein") == "moba/mobl"
    assert normalise.product_name("MobA/MobL protein") == "moba/mobl"
    assert normalise.product_name("DUF1234 domain-containing protein") == "duf1234"


def test_whitespace_and_case_are_canonicalised():
    assert normalise.product_name("  Conjugal   Transfer  Protein  TraG ") == \
        "conjugal transfer trag"


def test_a_name_that_is_only_a_suffix_survives_whole():
    """Stripping must not empty a product. 'protein' alone is a real, if useless, product
    name, and returning '' would merge it with every parse failure in the collection."""
    assert normalise.product_name("protein") == "protein"
    assert normalise.product_name("hypothetical protein") == "hypothetical protein"


def test_an_empty_input_gives_an_empty_string():
    assert normalise.product_name("") == ""
    assert normalise.product_name(None) == ""
```

- [ ] **Step 2: Run test to verify it fails**

Run: `envs/plasmidann/bin/python -m pytest tests/test_normalise.py -v`
Expected: FAIL, `ModuleNotFoundError: No module named 'plasmidann.normalise'`

- [ ] **Step 3: Write the normaliser**

Create `src/plasmidann/normalise.py`:

```python
"""Collapsing free-text product names to a canonical form.

WHY ONLY THE FREE-TEXT KINDS NEED THIS

Pfam family names, COG identifiers, GO terms, EC numbers and KO identifiers are controlled
vocabularies: a finite set of exact strings that can be enumerated and reviewed. Product
names from nr and Swiss-Prot are not. The same protein is deposited as 'Plasmid replication
initiator protein RepA', 'replication initiation protein' and 'putative replication
initiator protein' by three different submitters.

Left as written, one concept becomes dozens of categories, each with too few members to
test, and the frequency ranking that makes the review tractable does not work. Normalised,
the surface forms collapse and the head of the distribution is a bounded list.

WHAT IS DELIBERATELY NOT DONE HERE

No stemming, no synonym dictionary, no fuzzy matching. Those would merge products that a
reviewer should see separately, and the merge would be invisible in the output. This
function only removes tokens that carry no functional information: the submitter's hedging,
the generic suffixes, and inconsistent whitespace and case.
"""
import re

# Confidence hedges. They describe the submitter's certainty, not the protein, and keeping
# them splits every product into a confident and a hedged form.
_PREFIXES = ("putative ", "probable ", "predicted ", "possible ", "conserved ")

# Generic suffixes applied inconsistently across submissions. 'MobA/MobL family protein'
# and 'MobA/MobL protein' are one product.
_SUFFIXES = (" domain-containing protein", " domain containing protein",
             " family protein", " superfamily protein", " containing protein",
             " like protein", " protein")


def product_name(text):
    """Canonical form of a free-text product name.

    Lower-cased, hedges and generic suffixes removed, whitespace collapsed. Stripping never
    empties the string: a product that is nothing but a suffix is returned whole, because
    an empty result would merge it with every parse failure in the collection.
    """
    value = re.sub(r"\s+", " ", (text or "").strip().lower())
    if not value:
        return ""

    changed = True
    while changed:
        changed = False
        for prefix in _PREFIXES:
            if value.startswith(prefix) and len(value) > len(prefix):
                value, changed = value[len(prefix):].strip(), True
        for suffix in _SUFFIXES:
            if value.endswith(suffix) and len(value) > len(suffix):
                value, changed = value[:-len(suffix)].strip(), True
    return value or re.sub(r"\s+", " ", text.strip().lower())
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `envs/plasmidann/bin/python -m pytest tests/test_normalise.py -v`
Expected: PASS, 6 tests

- [ ] **Step 5: Write the mining configuration**

Create `config/category_mining.yaml`:

```yaml
# Patterns that GENERATE candidate category members for review.
#
# These patterns never run in the pipeline. build_categories.py applies them once, to the
# observed vocabulary in results/s4c/protein_labels.tsv, and writes a DRAFT config with the
# protein count beside every candidate. A person reviews that draft, removes what does not
# belong, and the reviewed result is what ships in config/label_categories.yaml, where
# matching is exact.
#
# Why generate rather than match at runtime, measured on Pfam-A 38.2:
#
#   /replicat/    matches  67 family descriptions   (the deleted hand list named 16)
#   /conjug/      matches  42                        (it named 15)
#   /transposas/  matches  76                        (it named 15)
#   /toxin/       matches 317  - and ABC_toxin_N is an insect toxin, ADPRTs_Tse2 is a T6SS
#                                effector; neither is a toxin-antitoxin system
#
# So the pattern supplies the scope a hand list cannot reach, and the review removes what
# the pattern cannot judge. Running the pattern at query time would keep the scope and lose
# the review, which is how 'PIN' matching inside an unrelated description sends someone to
# the bench to test the wrong thing.
#
# `search` fields, per category:
#   description   regular expression over the label's DESCRIPTION (Pfam DE, eggNOG
#                 description, normalised product name), case-insensitive
#   symbol        regular expression over a gene symbol; bacterial symbols are assigned in
#                 families, which makes them the most systematic axis available
#   cog_category  exact COG single-letter categories
#
# Every candidate is reported with the number of proteins carrying it, so the reviewer can
# see what including or excluding it costs before deciding.

mine:
  replication:
    description: 'replicat|repA|initiator protein'
    symbol: '^rep[A-Z]?$'
    cog_category: [L]
  mobilisation:
    description: 'relaxase|mobilis|mobiliz|MOB[A-Z]|nickase|oriT'
    symbol: '^(mob|nik|tra[IJK])[A-Z]?$'
  conjugation:
    description: 'conjug|type IV secretion|T4SS|pilus assembly|mating pair'
    symbol: '^(tra|trb|trw|vir)[A-Z]$'
  partition:
    description: 'partition|plasmid stability|segregation|centromere'
    symbol: '^(par|sop|stb)[A-C]$'
    cog_category: [D]
  transposition:
    description: 'transposas|integrase|recombinase|resolvase|insertion sequence'
    symbol: '^(tnp|ins|int|res)[A-Z]?$'
  toxin_antitoxin:
    description: 'toxin-antitoxin|antitoxin|addiction module|post-segregational'
    symbol: '^(ccd|rel|vap|hig|maz|par[DE]|hic|phd|doc|yef)[A-Z]?$'
  restriction_modification:
    description: 'restriction|modification methyl|methyltransferase|DNA methylase'
  defence:
    description: 'anti-phage|antiviral|abortive infection|CRISPR|Cas[0-9]'
```

- [ ] **Step 6: Write the failing test for the builder**

Append to `tests/test_scripts_smoke.py`:

```python
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
              [["s1", "pfam", "T1", "pfam_family", "RepA_N", "PF06970.16", "1e-40",
                "0.9", "Pfam-A", "38.2"],
               ["s2", "pfam", "T1", "pfam_family", "RepA_N", "PF06970.16", "1e-30",
                "0.9", "Pfam-A", "38.2"],
               ["s1", "eggnog", "S4b", "gene_symbol", "repA", "", "", "",
                "eggNOG", "5.0.2"],
               ["s3", "nr", "T4", "pgap_product",
                "putative plasmid replication initiator protein", "WP_1.1", "1e-20",
                "0.8", "NCBI nr", "2025-03-03"],
               ["s4", "nr", "T4", "pgap_product", "ABC transporter permease", "WP_2.1",
                "1e-20", "0.8", "NCBI nr", "2025-03-03"]])

    pfam_dat = fixture_dir / "pfam.dat"
    pfam_dat.write_text("# STOCKHOLM 1.0\n#=GF ID   RepA_N\n#=GF AC   PF06970.16\n"
                        "#=GF DE   Replication initiator protein A (RepA) N-terminus\n"
                        "#=GF TP   Family\n//\n")
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

    # Mined from the Pfam DESCRIPTION, which the family name alone does not contain.
    assert "RepA_N" in replication["pfam_family"]
    # Mined from the gene-symbol convention.
    assert "repA" in replication["gene_symbol"]
    # Mined from the NORMALISED product name, so the hedge did not hide it.
    assert any("replication initiat" in p for p in replication["pgap_product"])
    # Not mined: nothing about it matches, and it must not be swept in.
    assert "ABC transporter permease" not in replication.get("pgap_product", [])

    # The frequency table is the evidence for the review.
    rows = {(r["kind"], r["label"]): r for r in read_tsv(frequency)}
    assert rows[("pfam_family", "RepA_N")]["n_proteins"] == "2"
    assert rows[("gene_symbol", "repA")]["n_proteins"] == "1"


def test_build_categories_never_writes_the_live_config(fixture_dir):
    """The draft is reviewed by a person before it becomes the live config. A builder that
    wrote config/label_categories.yaml directly would put unreviewed pattern matches into
    the pipeline, which is the failure the whole design exists to prevent."""
    source = (pathlib.Path(__file__).parent.parent
              / "workflow" / "scripts" / "build_categories.py").read_text()

    assert "config/label_categories.yaml" not in source.replace(
        "# config/label_categories.yaml", ""), (
        "the builder references the live config path outside a comment")
```

- [ ] **Step 7: Write the builder**

Create `workflow/scripts/build_categories.py`:

```python
"""Build a DRAFT category configuration from the observed label vocabulary.

WHEN THIS RUNS

After a full annotation run, once results/s4c/protein_labels.tsv exists. Never as part of
the pipeline: it is a one-off that produces something for a person to read.

WHAT IT DOES AND WHY IT STOPS WHERE IT DOES

It applies the patterns in config/category_mining.yaml to the DESCRIPTIONS of the observed
labels and writes two files:

  category_draft.yaml    every candidate member of every category, with the number of
                         proteins carrying it, in the shape config/label_categories.yaml
                         expects
  label_frequency.tsv    the full observed vocabulary, per kind, ranked by protein count

It does NOT write config/label_categories.yaml. The draft is reviewed first, and the
reviewed result is what ships. Measured on Pfam-A 38.2, the reason for that separation:

  /replicat/    matches  67 family descriptions   (the deleted hand list named 16)
  /conjug/      matches  42                        (it named 15)
  /transposas/  matches  76                        (it named 15)
  /toxin/       matches 317  - including ABC_toxin_N, an insect toxin, and ADPRTs_Tse2, a
                               T6SS effector; neither is a toxin-antitoxin system

So the pattern reaches a scope no hand list can, and the review removes what the pattern
cannot judge. Applying the pattern at query time instead would keep the scope and lose the
review - and a substring match against a live vocabulary is how a mislabelled neighbourhood
sends someone to the bench to test the wrong thing.

The protein count beside each candidate is the evidence for the review: it says what
keeping or dropping that candidate costs, before the decision is made.
"""
import _ctx  # noqa: F401
import collections
import csv
import re

import yaml

from plasmidann import normalise, pfam_meta

# Label kinds whose text is free-form and must be normalised before it is counted or
# matched. The controlled kinds are exact strings and are used as written.
_FREE_TEXT = {"pgap_product", "swissprot_product", "eggnog_description",
              "pfam_description"}

mining = (yaml.safe_load(open(snakemake.input.mining).read()) or {}).get("mine") or {}
pfam = pfam_meta.load(snakemake.input.pfam_dat)

# ------------------------------------------------------------------------------------
# The observed vocabulary: how many distinct proteins carry each (kind, label).
#
# Proteins, not rows: a label seen on one protein through three tiers is one protein's
# worth of evidence, and counting rows would make widely searched labels look commoner.
# ------------------------------------------------------------------------------------
proteins = collections.defaultdict(set)
with open(snakemake.input.labels, newline="") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        label = r["label"]
        if r["kind"] in _FREE_TEXT:
            label = normalise.product_name(label)
        if label:
            proteins[(r["kind"], label)].add(r["protein_id"])

with open(snakemake.output.frequency, "w", newline="") as out:
    w = csv.writer(out, delimiter="\t")
    w.writerow(["kind", "label", "n_proteins", "description"])
    for (kind, label), ids in sorted(proteins.items(),
                                     key=lambda kv: (-len(kv[1]), kv[0])):
        description = pfam.get(label, {}).get("description", "") \
            if kind == "pfam_family" else ""
        w.writerow([kind, label, len(ids), description])


def describe(kind, label):
    """The text a pattern is matched against for one label.

    For a Pfam family this is the release description, which is where the functional words
    live - 'RepA_N' contains no matchable word, 'Replication initiator protein A (RepA)
    N-terminus' does. For everything else the label is its own description.
    """
    if kind == "pfam_family":
        return pfam.get(label, {}).get("description", "") or label
    return label


# ------------------------------------------------------------------------------------
# Candidates per category.
# ------------------------------------------------------------------------------------
draft = {}
for category, spec in mining.items():
    spec = spec or {}
    description_pattern = re.compile(spec["description"], re.I) \
        if spec.get("description") else None
    symbol_pattern = re.compile(spec["symbol"]) if spec.get("symbol") else None
    cog = set(spec.get("cog_category") or [])

    by_kind = collections.defaultdict(list)
    for (kind, label), ids in proteins.items():
        matched = False
        if kind == "cog_category":
            matched = label in cog
        elif kind == "gene_symbol":
            matched = bool(symbol_pattern and symbol_pattern.search(label))
        elif description_pattern:
            matched = bool(description_pattern.search(describe(kind, label)))
        if matched:
            by_kind[kind].append((len(ids), label))

    if by_kind:
        draft[category] = {
            "exact": {kind: [label for _, label in sorted(values, reverse=True)]
                      for kind, values in sorted(by_kind.items())},
            "_candidate_counts": {
                kind: {label: n for n, label in sorted(values, reverse=True)}
                for kind, values in sorted(by_kind.items())},
        }

header = (
    "# DRAFT - generated by workflow/scripts/build_categories.py. NOT the live config.\n"
    "#\n"
    "# Every entry below is a CANDIDATE produced by a pattern in\n"
    "# config/category_mining.yaml. Patterns reach a scope no hand-written list can, and\n"
    "# they also match things that do not belong: measured on Pfam-A 38.2, /toxin/ hits\n"
    "# 317 family descriptions of which ABC_toxin_N is an insect toxin and ADPRTs_Tse2 is\n"
    "# a T6SS effector.\n"
    "#\n"
    "# Review this file, delete what does not belong, delete the _candidate_counts blocks,\n"
    "# and copy the result into config/label_categories.yaml. Record for each category, in\n"
    "# docs/PARAMETER_PROVENANCE.md: the pattern that produced it, what was removed and\n"
    "# why, and the number of proteins it ends up covering.\n"
    "#\n"
    "# _candidate_counts gives the number of distinct proteins carrying each candidate, so\n"
    "# the cost of keeping or dropping one is visible before the decision is made.\n")

with open(snakemake.output.draft, "w") as out:
    out.write(header)
    yaml.safe_dump({"categories": draft}, out, default_flow_style=False, sort_keys=True)

print(f"build_categories: {len(proteins)} distinct (kind, label) pairs over "
      f"{len({p for ids in proteins.values() for p in ids})} proteins")
for category in sorted(draft):
    counts = draft[category]["_candidate_counts"]
    total = sum(sum(v.values()) for v in counts.values())
    per_kind = ", ".join(f"{k}={len(v)}" for k, v in sorted(counts.items()))
    print(f"  {category:<26} {sum(len(v) for v in counts.values()):>5} candidates "
          f"({per_kind}), {total} protein assignments")
```

- [ ] **Step 8: Run tests to verify they pass**

Run: `envs/plasmidann/bin/python -m pytest tests/test_normalise.py tests/test_scripts_smoke.py -k "normalise or build_categories" -v`
Expected: PASS

- [ ] **Step 9: Commit**

```bash
git add src/plasmidann/normalise.py workflow/scripts/build_categories.py \
        config/category_mining.yaml tests/test_normalise.py tests/test_scripts_smoke.py
git commit -m "$(cat <<'EOF'
feat(labels): add the category builder, which drafts rules for review

Turns the observed vocabulary into a draft config with the protein count
beside every candidate, so the grouping is decided on real data instead of
recollection. Mining reaches the scope a hand list cannot - measured on
Pfam-A 38.2, /replicat/ matches 67 family descriptions where the deleted
list named 16, /conjug/ 42 where it named 15 - and it also matches things
that do not belong, since /toxin/ hits 317 including an insect toxin and a
T6SS effector.

So the pattern generates candidates and never matches at runtime. What ships
in config/label_categories.yaml is the reviewed, enumerated list, and
categories.py keeps exact matching. A test asserts the builder never writes
the live config.

Free-text product names are normalised first, because the same protein is
deposited under several surface forms and unnormalised each becomes its own
category, too rare to test.

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>
EOF
)"
```

---

## What comes after this plan

Task 10 produces the draft; a person turns it into the live config. The sequence, once a
full annotation run exists:

1. Run `build_categories.py`. It writes `label_frequency.tsv` (the whole observed
   vocabulary, ranked by protein count) and `category_draft.yaml` (candidates per category,
   with counts).
2. **Review the draft.** Bounded work: roughly 67 lines for replication, 42 for
   conjugation, 76 for transposition, on the Pfam axis. Delete what does not belong -
   `CDT1` is a eukaryotic licensing factor, `ABC_toxin_N` is an insect toxin - and add a
   `sub:` entry where a finer name is worth keeping.
3. Copy the reviewed result into `config/label_categories.yaml`, point
   `config.label_categories` at it, and update
   `tests/test_categories.py::test_the_shipped_config_file_loads_and_is_still_the_identity`.
4. Record each category in `docs/PARAMETER_PROVENANCE.md`: the pattern that produced it,
   what was removed and why, and the proteins it covers.
5. Re-run S8c only. Nothing upstream changes, because the labels are already captured.

## Decisions recorded, not taken

These were considered and deliberately left out. Each is separate work with its own justification.

**CONJScan 2.1.0.** Installed in the `plasmidann` environment (`macsydata available` lists it) and never invoked. It is the native source for conjugation and mobilisation calls: MOB relaxase classes, T4CP and MPF types, resolved as *systems* on genomic order rather than as isolated domain hits, from published models (Cury et al. 2017; Abby et al. 2016). Adding it is one extra `--models CONJScan` run over the same input as the existing defence stage. Out of scope here because this plan commits only to capturing what already runs, and a plan that both adds a tool and rewires the label layer could not be reviewed cleanly. **Recommended as the next piece of work**, because conjugation is the category the label vocabulary covers least well: Pfam family names for Tra and Trb proteins are numerous and inconsistent, and only 218 reviewed Swiss-Prot entries carry the `Conjugation` keyword.

**UniProt keyword join.** `KW-0235 DNA replication` (7,114 reviewed entries), `KW-1277 Toxin-antitoxin system` (674), `KW-0614 Plasmid` (4,095) are curated per protein and applied systematically: measured on 8 reviewed toxin-antitoxin plasmid entries, 8 of 8 carry the keyword and only 2 of 8 carry `GO:0110001`. Task 2 captures the accession that makes this join possible. The join itself needs a pinned UniProt release and a reference file, and belongs with the grouping work that would consume it.

**pfam2go and a GO slim.** Measured and rejected as a primary grouping source: 4 of 16 replication families, 1 of 16 conjugation families, 0 of 11 mobilisation families. GO is captured as one label kind among many and nothing depends on it.

**PHROGs tier.** Specified as a fifth cascade tier and disabled at `config/cascade.yaml:137`. Each PHROG carries a native functional category among nine, which would make it the only tier with built-in categories. Enabling it is a cascade change, not a label change.

**Bakta.** Would supply PGAP product names, COG categories and GO on every protein from one `bakta_proteins` run. Explicitly deferred by the user. The nr `WP_` product names captured in Task 5 are the same PGAP vocabulary reached through a tool that already runs.

**Clonal redundancy between plasmids.** The Fisher test's unit is the plasmid, which removes copy-number inflation within a plasmid but not forty independent depositions of one clinical plasmid. `workflow/scripts/clonal_registry.py` holds the registry. Recorded as outstanding in `docs/PARAMETER_PROVENANCE.md`.

## Self-Review

**Spec coverage.** The request was: delete the curated list (Task 9); use only labels the tools provide (Tasks 2, 3, 4, 5); capture labels for the plus or minus three neighbourhood whether or not they group under one class (Task 9, Step 5 - the neighbourhood takes each neighbour's categories, which default to the raw label); run significance on categories (Task 8, applied in Task 9); defer the grouping into replication and mobilisation (Task 7 - the seam ships with no rules, and `config/label_categories.yaml` documents what adding them requires).

**Type consistency.** `labels_from_*` returns `list[dict]` with keys `kind`, `label`, `accession` throughout Tasks 5, 6 and 9. `CategoryMap.categories_for(kind, label) -> list[str]` is used identically in Tasks 7 and 9. `fisher_enrichment` returns the same key set everywhere it is consumed. `protein_labels.tsv` columns are declared once in `schemas.py` (Task 6, Step 3), written once (Step 4) and read once (Task 9, Step 5), with the same names.

**Known risk, stated rather than hidden.** Task 9 changes `family_context.tsv` from wide to long, and three consumers read it: `prioritise.py`, `annotation_report.py:84-85` and `tests/test_context.py`. Steps 7 and 9 instruct rewiring all three, and Step 8's grep is the check. A worker who treats a failure in those files as unrelated will leave the pipeline broken; the grep must come back empty before Task 9 is complete.

**Deliberate omission.** `labels_from_defence` and `labels_from_integron` are written and tested in Task 5 but not called by `protein_labels.py` in Task 6, because defence and integron calls are per ORF rather than per unique protein and enter the context stage as *intervals* (Task 9, Step 5), which is the stronger statement. They exist for the grouping work, which will need the MacSyFinder system vocabulary. This is the one piece of the plan that builds something before its consumer exists, and it is here because the alternative is a second pass over `labels.py` for two six-line functions.
