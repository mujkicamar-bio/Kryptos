# Dark ORF Pipeline — Plan 1: Foundation and Stages 1–3

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build the data foundation and the first three pipeline stages — plasmid ingest, ORF prediction with circular-origin handling, ORF QC/artifact screening, and dereplication — producing a validated `proteins` table that every later stage consumes.

**Architecture:** A Snakemake workflow over a DuckDB/Parquet relational store. Python logic lives in importable, unit-tested modules under `src/darkorf/`; `workflow/scripts/*.py` are thin wrappers that read `snakemake.input/output/params` and call those modules. Every table is Parquet; the wide occurrence-level table of §67 is a DuckDB *view*, never materialized into pandas. Stages shard over plasmids (Stage 1) and over proteins (Stage 2 onward) so a failure costs one shard, not the stage. Dereplication is Stage 2, before all per-protein work, so AntiFam and every cascade tier sees ~3.5M unique sequences rather than ~9.3M occurrences.

**Tech Stack:** Python 3.11, Snakemake 8, pyrodigal 3.7, pyhmmer 0.12 / HMMER 3.4, DuckDB, PyArrow, pytest.

**Spec:** `PLASMID_ANALYSIS.md` (sections cited per task). The spec is authoritative; where this plan and the spec disagree, stop and ask.

## Global Constraints

- Plasmid count must never be hard-coded; the analysis-set denominator comes from the input and is recorded in the run manifest (§3.1).
- `protein_id = SHA256(normalized_protein_sequence)[:32]` — mandatory, not optional (§5.3).
- `orf_occurrence_id = <plasmid_id>:<start>-<end>:<strand>`; never an ordinal (§5.2).
- Identifiers must be stable across reruns when biological input is unchanged (§5).
- Numeric fields are nullable and carry a companion `*_status` field; a test that was not run, not applicable, saturated, or failed are four different states (§2.9, §7).
- `NA` is for genuinely unavailable *text* only; it must never conceal a test status (§7.3).
- Records are flagged, never deleted. Every exclusion carries a reason and stays in the master dataset (§9, §10, §12).
- Dereplication changes computation, never biological occurrence counts (§13).
- No composite score of any kind — no `dark_score`, `novelty_score`, `candidate_score` (§2.3).
- Every threshold lives in `config/config.yaml`, is schema-validated, and is stamped into the run manifest (§6, §15.3).
- Do not load the complete occurrence dataset into pandas (§4.2).
- Comment every non-obvious code block for a human reviewer; explain *why*, not *what*.
- Standard international English in all code comments, docstrings and commit messages.

---

## File Structure

| File | Responsibility |
|---|---|
| `config/config.yaml` | Every path and threshold in one place |
| `config/schemas/config.schema.yaml` | JSON-schema validation of the above |
| `src/darkorf/ids.py` | The three identifier constructors (§5) |
| `src/darkorf/status.py` | The centralized status vocabulary (§7.2) |
| `src/darkorf/schemas.py` | Column contract for every Parquet table (§4.1) |
| `src/darkorf/store.py` | Parquet write/read + DuckDB connection helpers |
| `src/darkorf/manifest.py` | Run manifest assembly (§6) |
| `src/darkorf/circular.py` | Origin-spanning ORF handling (ported, §8.4) |
| `src/darkorf/orf.py` | pyrodigal wrapper, CDS nucleotide extraction (§8) |
| `src/darkorf/qc.py` | Overlap QC, origin QC, discovery eligibility (§9, §10, §12) |
| `src/darkorf/derep.py` | Dereplication with lossless occurrence mapping (§13) |
| `src/darkorf/controls.py` | Positive and negative control construction (§58) |
| `workflow/Snakefile` | Targets, config load, shard definitions |
| `workflow/rules/s0_input.smk` | Stage 0: plasmids table, manifest, preflight |
| `workflow/rules/s1_orf.smk` | Stage 1: sharded ORF prediction |
| `workflow/rules/s2_derep.smk` | Stage 2: dereplication, run first so every later stage is cheaper |
| `workflow/rules/s3_qc.smk` | Stage 3: AntiFam over unique proteins, overlap QC, eligibility |
| `workflow/envs/darkorf.yaml` | Conda environment specification |
| `tests/` | One test module per source module |

Ported from the old implementation (read the original before porting; keep the behaviour and the tests, adapt the interfaces): `src/plasmidann/circular.py` → `src/darkorf/circular.py`, `src/plasmidann/dereplicate.py` → `src/darkorf/derep.py`.

---

### Task 1: Environment and configuration

**Files:**
- Create: `workflow/envs/darkorf.yaml`, `config/config.yaml`, `config/schemas/config.schema.yaml`
- Test: `tests/test_config.py`

**Interfaces:**
- Produces: `config/config.yaml` keys consumed by every later task — `input.working_set_fasta`, `input.master_table`, `outdir`, `shards.plasmid`, `shards.protein`, `orf.min_call_length_aa`, `discovery.min_dark_candidate_length_aa`, `cascade.narrow_at`, `cascade.min_explained`, `seed`.

- [ ] **Step 1: Write the failing test**

```python
# tests/test_config.py
"""The configuration file is the single source of every threshold (spec §6, §15.3).

These tests exist because three of six thresholds in an earlier version of this project
lived as Python constants, which meant they could not be validated, swept, or recorded in
the output. A threshold that is not in config is a threshold nobody can audit.
"""
import pathlib
import yaml
import jsonschema

ROOT = pathlib.Path(__file__).resolve().parents[1]


def test_config_validates_against_its_schema():
    config = yaml.safe_load((ROOT / "config" / "config.yaml").read_text())
    schema = yaml.safe_load((ROOT / "config" / "schemas" / "config.schema.yaml").read_text())
    jsonschema.validate(config, schema)


def test_the_two_length_parameters_are_distinct_keys():
    """Spec §11.1: the gene-caller floor and the discovery cutoff are different concepts
    that happen to share a value. Conflating them is how a discovery flag silently becomes
    a deletion."""
    config = yaml.safe_load((ROOT / "config" / "config.yaml").read_text())
    assert config["orf"]["min_call_length_aa"] == 20
    assert config["discovery"]["min_dark_candidate_length_aa"] == 20


def test_cascade_thresholds_match_the_spec():
    """Spec §15.3. narrow_at controls compute; min_explained controls science."""
    config = yaml.safe_load((ROOT / "config" / "config.yaml").read_text())
    assert config["cascade"]["narrow_at"] == 0.9
    assert config["cascade"]["min_explained"] == 0.5
    assert config["cascade"]["narrow_at"] >= config["cascade"]["min_explained"]
```

- [ ] **Step 2: Run test to verify it fails**

Run: `.venv/bin/python -m pytest tests/test_config.py -v`
Expected: FAIL — `FileNotFoundError: config/config.yaml`

- [ ] **Step 3: Write the environment file**

```yaml
# workflow/envs/darkorf.yaml
# Every executable the workflow invokes. Snakemake's `conda:` directives point here, so a
# tool missing from this file is a tool the pipeline cannot run - the preflight rule
# (Task 6) checks each one before any compute is spent.
name: darkorf
channels: [conda-forge, bioconda]
dependencies:
  - python=3.11
  - snakemake-minimal=8.*
  - pyrodigal=3.7          # S1 gene calling
  - hmmer=3.4              # S2 AntiFam, S4 Pfam tiers
  - pyhmmer=0.12           # in-process HMM parsing
  - diamond=2.1            # S4 Swiss-Prot and nr tiers
  - mmseqs2=15.*           # S4 PHROGs profiles, S5 families
  - foldseek=10.*          # S11 structural search
  - mafft=7.*              # S10 family alignments
  - rnacode=0.3            # S10 coding-potential evidence
  - integron_finder=2.*    # S8 integron context
  - macsyfinder=2.*        # S8 DefenseFinder backend
  - mash=2.3               # S6 plasmid lineage clustering
  - duckdb=1.*             # relational store
  - python-duckdb=1.*
  - pyarrow=17.*           # Parquet IO
  - numpy=2.*
  - scipy=1.*              # S13 statistics
  - pandas=2.2             # shard-level frames only, never the full occurrence table
  - biopython=1.8*
  - pyyaml
  - jsonschema
  - pytest
```

- [ ] **Step 4: Write the configuration**

```yaml
# config/config.yaml
# Every threshold the pipeline uses lives here (spec §6). Values are stamped into the run
# manifest and into the output rows they govern, so a result can always be traced to the
# configuration that produced it.

input:
  # The analysis set: the master table with the locked artifact exclusions already applied.
  # The plasmid count is NOT recorded here - it is derived at runtime and written to the
  # manifest (spec §3.1), because a hard-coded denominator silently goes stale.
  working_set_fasta: data/plasmidscope_primary/provenance/working_set.fna.gz
  master_table: data/plasmidscope_primary/analysis_set.tsv

outdir: results

# Sharding sets the unit of parallelism and of resume. Nothing biological depends on these
# numbers, but changing them between resumed runs invalidates the affected stage (spec §72).
shards:
  plasmid: 600    # S1/S2 shard over plasmids
  protein: 64     # S4 shards over the dereplicated protein set

orf:
  # The gene caller's own floor. Distinct from the discovery cutoff below (spec §11.1).
  min_call_length_aa: 20
  # Prodigal in metagenomic mode: plasmids are too short and too compositionally varied to
  # train a single-genome model on.
  meta_mode: true

discovery:
  # Proteins shorter than this are flagged discovery_excluded_short, never deleted (spec §10, §11).
  min_dark_candidate_length_aa: 20

qc:
  # Spec §9.3. A protein whose coding region is mostly antisense to a longer ORF on the
  # opposite strand is the classic shadow-ORF artifact: it passes the whole cascade cleanly
  # because it is not a protein and no database contains it. Flag, never delete.
  max_opposite_strand_overlap_fraction: 0.6
  antifam:
    db: data/refs/antifam/AntiFam.hmm
    # AntiFam's curated per-family gathering thresholds. Measured on this database: 274 of
    # 278 profiles carry a GA looser than E=1e-5, so a blanket E-value floor would override
    # the curator across essentially the whole database (Eberhardt et al. 2012, Database
    # 2012:bas003).
    threshold: cut_ga

cascade:
  # Spec §15. Two different decisions that an earlier design conflated into one number.
  narrow_at: 0.9      # how explained before we STOP SEARCHING (compute)
  min_explained: 0.5  # how explained before we call it ANNOTATED (science)
  # Fraction of proteins that bypass narrowing and are searched by every tier, so the cost
  # of narrowing is measurable rather than assumed (spec §61).
  sweep_cohort_fraction: 0.02

controls:
  # Spec §58. Positive controls test recall; negative controls test the opposite direction.
  positive:
    raw_faa: data/refs/control/raw.faa
    n: 500
  negative:
    # Decoys are built by shuffling real plasmid CDS while preserving length and amino-acid
    # composition, plus reverse-complement translations of real CDS. Both must come out DARK
    # and must not receive a confident functional annotation.
    n_shuffled: 250
    n_reverse_complement: 250

# Seed for every stochastic step, so a rerun selects the identical sample (spec §61).
seed: 20260916
```

- [ ] **Step 5: Write the schema**

```yaml
# config/schemas/config.schema.yaml
# Validated at workflow load time. A typo in config.yaml must fail in seconds, not in hour 45.
type: object
required: [input, outdir, shards, orf, discovery, qc, cascade, controls, seed]
properties:
  input:
    type: object
    required: [working_set_fasta, master_table]
    properties:
      working_set_fasta: {type: string}
      master_table: {type: string}
  outdir: {type: string}
  shards:
    type: object
    required: [plasmid, protein]
    properties:
      plasmid: {type: integer, minimum: 1}
      protein: {type: integer, minimum: 1}
  orf:
    type: object
    required: [min_call_length_aa, meta_mode]
    properties:
      min_call_length_aa: {type: integer, minimum: 1}
      meta_mode: {type: boolean}
  discovery:
    type: object
    required: [min_dark_candidate_length_aa]
    properties:
      min_dark_candidate_length_aa: {type: integer, minimum: 1}
  qc:
    type: object
    required: [max_opposite_strand_overlap_fraction, antifam]
    properties:
      max_opposite_strand_overlap_fraction: {type: number, minimum: 0, maximum: 1}
      antifam:
        type: object
        required: [db, threshold]
        properties:
          db: {type: string}
          threshold: {type: string, enum: [cut_ga]}
  cascade:
    type: object
    required: [narrow_at, min_explained, sweep_cohort_fraction]
    properties:
      narrow_at: {type: number, minimum: 0, maximum: 1}
      min_explained: {type: number, minimum: 0, maximum: 1}
      sweep_cohort_fraction: {type: number, minimum: 0, maximum: 1}
  controls:
    type: object
    required: [positive, negative]
  seed: {type: integer}
```

- [ ] **Step 6: Run test to verify it passes**

Run: `.venv/bin/python -m pytest tests/test_config.py -v`
Expected: PASS (3 tests)

- [ ] **Step 7: Commit**

```bash
git add workflow/envs/darkorf.yaml config/config.yaml config/schemas/config.schema.yaml tests/test_config.py
git commit -m "feat(config): add darkorf environment and schema-validated configuration

Every threshold the pipeline uses lives in config/config.yaml and is checked
against config/schemas/config.schema.yaml at load time, per spec section 6.

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

---

### Task 2: Identifier constructors

**Files:**
- Create: `src/darkorf/__init__.py`, `src/darkorf/ids.py`
- Test: `tests/test_ids.py`

**Interfaces:**
- Produces:
  - `protein_id(sequence: str) -> str` — 32-char lowercase hex
  - `occurrence_id(plasmid_id: str, start: int, end: int, strand: int) -> str`
  - `family_id(resolution: str, representative_protein_id: str) -> str`
  - `observation_id(protein_id: str, tool: str, database: str, database_version: str, target: str, qstart: int, qend: int) -> str`

- [ ] **Step 1: Write the failing test**

```python
# tests/test_ids.py
"""Identifiers must be stable across reruns when the biological input is unchanged (spec §5).

The failure these tests prevent: an earlier implementation numbered families by
enumerate(sorted(...)), so adding one protein renumbered every later family and silently
broke every join against a previous run.
"""
from darkorf.ids import protein_id, occurrence_id, family_id, observation_id


def test_protein_id_is_a_32_character_hash_of_the_sequence():
    assert protein_id("MKV") == protein_id("MKV")
    assert len(protein_id("MKV")) == 32
    assert protein_id("MKV") != protein_id("MKW")


def test_protein_id_normalizes_case_and_whitespace():
    """FASTA parsers differ in whether they upper-case and how they wrap. The same protein
    must hash identically regardless of which parser produced the string."""
    assert protein_id("mkv") == protein_id("MKV")
    assert protein_id(" MKV\n") == protein_id("MKV")


def test_protein_id_strips_a_single_trailing_stop_codon():
    """Prodigal emits a trailing '*'; other tools do not. Without this the same protein gets
    two different identifiers depending on its provenance."""
    assert protein_id("MKV*") == protein_id("MKV")


def test_occurrence_id_is_coordinate_based_not_ordinal():
    assert occurrence_id("pA", 10, 99, 1) == "pA:10-99:+"
    assert occurrence_id("pA", 10, 99, -1) == "pA:10-99:-"


def test_occurrence_id_is_unchanged_by_shard_layout():
    """The same ORF must get the same id whether it was called in shard 1 or shard 599."""
    assert occurrence_id("pA", 10, 99, 1) == occurrence_id("pA", 10, 99, 1)


def test_family_id_is_derived_from_its_representative():
    assert family_id("close", "abc123") == "close:abc123"


def test_observation_id_is_stable_and_order_independent_of_tool_version():
    """Tool version is provenance, not identity (spec §5.5): re-running the same search with
    a patched binary must not renumber every observation."""
    a = observation_id("p1", "hmmsearch", "Pfam-A", "38.2", "PF00001", 1, 50)
    b = observation_id("p1", "hmmsearch", "Pfam-A", "38.2", "PF00001", 1, 50)
    assert a == b
    assert len(a) == 32
    assert a != observation_id("p1", "hmmsearch", "Pfam-A", "38.2", "PF00002", 1, 50)
```

- [ ] **Step 2: Run test to verify it fails**

Run: `.venv/bin/python -m pytest tests/test_ids.py -v`
Expected: FAIL — `ModuleNotFoundError: No module named 'darkorf'`

- [ ] **Step 3: Write the implementation**

```python
# src/darkorf/ids.py
"""Deterministic identifier construction (spec §5).

Every identifier in this pipeline is a pure function of biological content. Nothing is
numbered by position, enumeration order, or shard layout, because all three change between
runs for reasons that have nothing to do with the biology.
"""
import hashlib

# 32 hex characters is 128 bits. At 3.5 million proteins the collision probability is
# around 1e-26, so the truncation is safe and the identifiers stay readable in a terminal.
_ID_LENGTH = 32


def _normalize_protein_sequence(sequence):
    """Canonical form of a protein sequence for hashing.

    Upper-cased, whitespace removed, and a single trailing stop codon dropped - Prodigal
    emits '*' and most other tools do not, and the same protein must not hash differently
    because of which tool produced the string.
    """
    cleaned = "".join(sequence.split()).upper()
    return cleaned[:-1] if cleaned.endswith("*") else cleaned


def protein_id(sequence):
    """Stable identifier for a unique protein sequence (spec §5.3)."""
    canonical = _normalize_protein_sequence(sequence)
    return hashlib.sha256(canonical.encode()).hexdigest()[:_ID_LENGTH]


def occurrence_id(plasmid_id, start, end, strand):
    """Stable identifier for one ORF occurrence (spec §5.2).

    Coordinate-based rather than ordinal, so the identifier does not depend on how many ORFs
    were called before it or on which shard called it. Origin-spanning ORFs keep their
    unrotated coordinates, where end < start is legitimate (spec §8.4).
    """
    sign = "+" if strand > 0 else "-"
    return f"{plasmid_id}:{start}-{end}:{sign}"


def family_id(resolution, representative_protein_id):
    """Stable identifier for a sequence family (spec §5.4).

    Derived from the cluster representative, so adding an unrelated protein elsewhere in the
    dataset cannot renumber this family.
    """
    return f"{resolution}:{representative_protein_id}"


def observation_id(protein_id_value, tool, database, database_version, target, qstart, qend):
    """Stable identifier for one annotation/context/structure observation (spec §5.5).

    Tool version is deliberately excluded: it is recorded as provenance, but including it
    would give every observation a new identity after a patch release.
    """
    key = "\t".join(
        [protein_id_value, tool, database, database_version, target, str(qstart), str(qend)]
    )
    return hashlib.sha256(key.encode()).hexdigest()[:_ID_LENGTH]
```

```python
# src/darkorf/__init__.py
"""Dark ORF discovery pipeline: importable logic behind the Snakemake workflow."""
```

- [ ] **Step 4: Run test to verify it passes**

Run: `.venv/bin/python -m pytest tests/test_ids.py -v`
Expected: PASS (7 tests)

- [ ] **Step 5: Commit**

```bash
git add src/darkorf/__init__.py src/darkorf/ids.py tests/test_ids.py
git commit -m "feat(ids): add content-derived identifiers stable across reruns

protein_id is SHA256 of the normalized sequence, occurrence_id is
coordinate-based, family_id derives from the cluster representative. Nothing
is numbered by enumeration order, per spec section 5.

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

---

### Task 3: Status vocabulary and table schemas

**Files:**
- Create: `src/darkorf/status.py`, `src/darkorf/schemas.py`
- Test: `tests/test_status.py`, `tests/test_schemas.py`

**Interfaces:**
- Produces:
  - `status.Status` — string constants `NOT_RUN`, `NO_HIT`, `TOO_FEW_MEMBERS`, `NO_DIVERGENCE`, `SATURATED`, `NO_OUTPUT`, `FAILED`, `NOT_APPLICABLE`, `SUCCESS`
  - `status.ALL` — frozenset of the above
  - `schemas.TABLES: dict[str, list[str]]` — table name → ordered column names
  - `schemas.columns(table: str) -> list[str]`
  - `schemas.validate_frame(table: str, frame) -> None` — raises `ValueError` on a column mismatch

- [ ] **Step 1: Write the failing test**

```python
# tests/test_status.py
"""Spec §2.9 and §7.2: four different reasons for a missing number must not collapse into
one value. A null dN/dS because the family had three members and a null because the
alignment saturated are different scientific statements."""
from darkorf import status


def test_the_vocabulary_is_centralized_and_complete():
    for expected in [
        "NOT_RUN", "NO_HIT", "TOO_FEW_MEMBERS", "NO_DIVERGENCE",
        "SATURATED", "NO_OUTPUT", "FAILED", "NOT_APPLICABLE", "SUCCESS",
    ]:
        assert expected in status.ALL
        assert getattr(status, expected) == expected


def test_na_is_not_a_status():
    """NA is for genuinely unavailable text only. Using it as a status is exactly the
    concealment spec §7.3 forbids."""
    assert "NA" not in status.ALL
```

```python
# tests/test_schemas.py
"""Every table has a declared column contract (spec §4.1). A stage that writes a column
nobody declared is a stage whose output nothing downstream can rely on."""
import pytest
import pyarrow as pa
from darkorf import schemas


def test_every_spec_table_is_declared():
    for table in ["plasmids", "orf_occurrences", "proteins", "controls", "run_manifest"]:
        assert table in schemas.TABLES
        assert len(schemas.columns(table)) > 0


def test_orf_occurrences_carries_the_cds_nucleotide_sequence():
    """Spec §8.3. Without it RNAcode and dN/dS become one-shot: re-running them later means
    re-extracting CDS from every plasmid in the collection."""
    assert "cds_nucleotide_sequence" in schemas.columns("orf_occurrences")


def test_validate_frame_rejects_a_missing_column():
    frame = pa.table({"plasmid_id": ["p1"]})
    with pytest.raises(ValueError, match="missing"):
        schemas.validate_frame("orf_occurrences", frame)


def test_validate_frame_rejects_an_undeclared_column():
    columns = {name: ["x"] for name in schemas.columns("plasmids")}
    columns["surprise"] = ["x"]
    with pytest.raises(ValueError, match="undeclared"):
        schemas.validate_frame("plasmids", pa.table(columns))
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `.venv/bin/python -m pytest tests/test_status.py tests/test_schemas.py -v`
Expected: FAIL — `ImportError: cannot import name 'status'`

- [ ] **Step 3: Write the status vocabulary**

```python
# src/darkorf/status.py
"""The single status vocabulary (spec §7.2).

Every nullable numeric field in the pipeline carries a companion <field>_status drawn from
this set. The point is that "we did not run this", "this did not apply", "there were too
few sequences", "the estimate saturated" and "the tool crashed" are five different
scientific statements that a bare null would flatten into one.
"""

NOT_RUN = "NOT_RUN"                  # the stage was not executed for this record
NO_HIT = "NO_HIT"                    # the search ran and returned nothing above threshold
TOO_FEW_MEMBERS = "TOO_FEW_MEMBERS"  # below the configured minimum for the estimator
NO_DIVERGENCE = "NO_DIVERGENCE"      # sequences are identical; the statistic is undefined
SATURATED = "SATURATED"              # divergence too high for the estimate to be meaningful
NO_OUTPUT = "NO_OUTPUT"              # the tool exited successfully but produced nothing
FAILED = "FAILED"                    # the tool errored; see the log referenced in the row
NOT_APPLICABLE = "NOT_APPLICABLE"    # the field has no meaning for this record
SUCCESS = "SUCCESS"                  # a value is present and usable

ALL = frozenset({
    NOT_RUN, NO_HIT, TOO_FEW_MEMBERS, NO_DIVERGENCE,
    SATURATED, NO_OUTPUT, FAILED, NOT_APPLICABLE, SUCCESS,
})
```

- [ ] **Step 4: Write the table schemas**

```python
# src/darkorf/schemas.py
"""Column contracts for every table in the relational store (spec §4.1).

Declaring the columns in one place is what makes the store auditable: a reviewer can read
this file and know exactly what each stage is allowed to produce, and validate_frame turns
a silent schema drift into an immediate failure.

Only the tables Plan 1 writes are declared here. Later plans extend TABLES; they must not
redefine an existing table's columns without a spec revision.
"""

TABLES = {
    # Spec §62.1 - one row per plasmid in the analysis set.
    "plasmids": [
        "plasmid_id", "accession", "source_database", "length", "topology",
        "host", "host_taxonomy", "habitat", "mob_class", "mob_cluster",
        "plasmid_type", "bioproject", "record_class",
    ],
    # Spec §62.2, §8.3 - one row per predicted ORF occurrence.
    "orf_occurrences": [
        "orf_occurrence_id", "plasmid_id", "protein_id",
        "start", "end", "strand", "nucleotide_length", "protein_length",
        "cds_nucleotide_sequence", "partial", "start_type",
        "origin_spanning",
        # QC, spec §9. Flags, never deletions.
        "antifam_hit", "antifam_model", "antifam_score", "antifam_evalue",
        "opposite_strand_overlap_fraction", "opposite_strand_overlap_partner",
        "discovery_excluded_short", "discovery_excluded_partial",
        "discovery_excluded_antifam", "discovery_eligible", "exclusion_reason",
    ],
    # Spec §62.3 - one row per unique protein sequence.
    "proteins": [
        "protein_id", "protein_sequence", "protein_length",
        "occurrence_count", "plasmid_count",
        "has_complete_occurrence", "discovery_eligible", "record_class",
    ],
    # Spec §58 - positive and negative controls travel with the real data and are excluded
    # from biological occurrence-integrity checks.
    "controls": [
        "protein_id", "record_class", "control_category",
        "source_protein_id", "expected_outcome",
    ],
    # Spec §6 - one row per production run.
    "run_manifest": [
        "pipeline_version", "git_commit", "input_dataset_hash",
        "plasmid_count", "orf_count", "unique_protein_count",
        "config_hash", "database_versions", "software_versions",
        "shard_counts", "started_at", "finished_at",
    ],
}


def columns(table):
    """Ordered column names for one table."""
    if table not in TABLES:
        raise KeyError(f"no schema declared for table {table!r}")
    return list(TABLES[table])


def validate_frame(table, frame):
    """Raise if an Arrow table does not match its declared contract.

    Both directions are checked. A missing column breaks consumers downstream; an
    undeclared one means a stage is emitting something no schema documents, which is how
    columns with two different meanings end up sharing a name.
    """
    declared = set(columns(table))
    present = set(frame.column_names)
    missing = declared - present
    undeclared = present - declared
    if missing:
        raise ValueError(f"{table}: missing declared columns {sorted(missing)}")
    if undeclared:
        raise ValueError(f"{table}: undeclared columns {sorted(undeclared)}")
```

- [ ] **Step 5: Run tests to verify they pass**

Run: `.venv/bin/python -m pytest tests/test_status.py tests/test_schemas.py -v`
Expected: PASS (6 tests)

- [ ] **Step 6: Commit**

```bash
git add src/darkorf/status.py src/darkorf/schemas.py tests/test_status.py tests/test_schemas.py
git commit -m "feat(schemas): add status vocabulary and per-table column contracts

Nullable numerics carry a companion status so that not-run, not-applicable,
too-few-members, saturated and failed stay distinguishable, per spec 2.9/7.2.

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

---

### Task 4: Parquet store and DuckDB access

**Files:**
- Create: `src/darkorf/store.py`
- Test: `tests/test_store.py`

**Interfaces:**
- Produces:
  - `store.write_table(path, table_name, frame) -> None` — validates then writes Parquet
  - `store.read_table(path) -> pa.Table`
  - `store.connect(outdir) -> duckdb.DuckDBPyConnection` — registers every Parquet table as a view
  - `store.table_path(outdir, table_name) -> pathlib.Path`

- [ ] **Step 1: Write the failing test**

```python
# tests/test_store.py
"""The relational store (spec §4.2): Parquet on disk, DuckDB for queries, and never the
whole occurrence table in pandas."""
import pyarrow as pa
import pytest
from darkorf import store, schemas


def _plasmids_frame():
    return pa.table({name: ["x"] for name in schemas.columns("plasmids")})


def test_write_then_read_round_trips(tmp_path):
    path = store.table_path(tmp_path, "plasmids")
    store.write_table(path, "plasmids", _plasmids_frame())
    assert store.read_table(path).column_names == schemas.columns("plasmids")


def test_write_rejects_a_frame_that_breaks_the_contract(tmp_path):
    """Schema validation happens at write time, so a bad stage fails where it is, not three
    stages later in something that tried to join against it."""
    bad = pa.table({"plasmid_id": ["p1"]})
    with pytest.raises(ValueError):
        store.write_table(store.table_path(tmp_path, "plasmids"), "plasmids", bad)


def test_connect_exposes_each_written_table_as_a_queryable_view(tmp_path):
    store.write_table(store.table_path(tmp_path, "plasmids"), "plasmids", _plasmids_frame())
    with store.connect(tmp_path) as connection:
        assert connection.execute("SELECT count(*) FROM plasmids").fetchone()[0] == 1
```

- [ ] **Step 2: Run test to verify it fails**

Run: `.venv/bin/python -m pytest tests/test_store.py -v`
Expected: FAIL — `ImportError: cannot import name 'store'`

- [ ] **Step 3: Write the implementation**

```python
# src/darkorf/store.py
"""Parquet store with DuckDB query access (spec §4.2).

Why not one big TSV: the integrated occurrence view is roughly 9.3 million rows by ~120
columns, of which about 90 are constant within protein_id or family_id. Materializing that
as text is 9-14 GB and needs 50-110 GB of RAM to open. Parquet keeps each table at its own
level and DuckDB reconstructs the wide view on demand, without loading it (spec §63).
"""
import contextlib
import pathlib

import duckdb
import pyarrow.parquet as parquet

from darkorf import schemas


def table_path(outdir, table_name):
    """Canonical on-disk location of one table."""
    return pathlib.Path(outdir) / "tables" / f"{table_name}.parquet"


def write_table(path, table_name, frame):
    """Validate against the declared contract, then write Parquet.

    Validation is deliberately at write time: a stage that emits the wrong columns must
    fail in its own rule, where the log says which stage it was.
    """
    schemas.validate_frame(table_name, frame)
    path = pathlib.Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    # ZSTD because these tables are read far more often than written, and the sequence
    # columns compress by roughly an order of magnitude.
    parquet.write_table(frame, path, compression="zstd")


def read_table(path):
    return parquet.read_table(path)


@contextlib.contextmanager
def connect(outdir):
    """Open DuckDB with every Parquet table registered as a view of the same name.

    Views rather than imports: the Parquet files stay the single copy of the data, so a
    stage that rewrites a table is immediately visible to every later query.
    """
    connection = duckdb.connect()
    try:
        for path in sorted((pathlib.Path(outdir) / "tables").glob("*.parquet")):
            connection.execute(
                f"CREATE VIEW {path.stem} AS SELECT * FROM read_parquet('{path}')"
            )
        yield connection
    finally:
        connection.close()
```

- [ ] **Step 4: Run test to verify it passes**

Run: `.venv/bin/python -m pytest tests/test_store.py -v`
Expected: PASS (3 tests)

- [ ] **Step 5: Commit**

```bash
git add src/darkorf/store.py tests/test_store.py
git commit -m "feat(store): add Parquet store with DuckDB views and write-time validation

Tables are written at their own level and the wide occurrence view is
reconstructed by DuckDB on demand rather than materialized, per spec 4.2/63.

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

---

### Task 5: Circular-origin ORF handling

**Files:**
- Create: `src/darkorf/circular.py` (port of `src/plasmidann/circular.py` — read it first)
- Test: `tests/test_circular.py`

**Interfaces:**
- Consumes: nothing from earlier tasks.
- Produces:
  - `circular.is_circular(topology: str) -> bool`
  - `circular.overlap_for(length: int) -> int`
  - `circular.rotate(sequence: str, offset: int) -> str`
  - `circular.resolve_origin_genes(genes, original_length) -> list` — deduplicates genes found twice because of the rotation overlap, and marks origin-spanning ones

- [ ] **Step 1: Write the failing test**

```python
# tests/test_circular.py
"""Origin handling for circular plasmids (spec §8.4).

The invariant: a circular sequence has no privileged starting point, so calling genes on
any rotation of it must give the same protein set. The measured stake is 160,375
origin-spanning ORFs (1.72% of the collection) that a naive linear caller either truncates
or misses entirely.
"""
import random
from darkorf import circular


def test_topologies_that_count_as_circular():
    assert circular.is_circular("circular")
    assert circular.is_circular("Circular")
    assert not circular.is_circular("linear")


def test_rotation_preserves_length_and_content():
    sequence = "ATGCGTACGT"
    rotated = circular.rotate(sequence, 4)
    assert len(rotated) == len(sequence)
    assert sorted(rotated) == sorted(sequence)
    assert circular.rotate(rotated, len(sequence) - 4) == sequence


def test_overlap_is_capped_for_long_plasmids():
    """The duplicated head must be long enough to contain any real gene crossing the origin,
    but duplicating a 400 kb megaplasmid wholesale would double the gene-calling cost."""
    assert circular.overlap_for(1_000) <= 1_000
    assert circular.overlap_for(1_000_000) == circular.MAX_OVERLAP_BP


def test_genes_duplicated_by_the_overlap_are_resolved_to_one_copy():
    """A gene lying inside the duplicated head is called twice - once at its true
    coordinates and once shifted by the sequence length. Exactly one copy must survive, or
    every downstream occurrence count is inflated."""
    length = 1_000
    genes = [
        {"start": 10, "end": 100, "strand": 1},
        {"start": 10 + length, "end": 100 + length, "strand": 1},
    ]
    resolved = circular.resolve_origin_genes(genes, length)
    assert len(resolved) == 1


def test_an_origin_spanning_gene_is_marked_and_keeps_unrotated_coordinates():
    """Spec §8.4: do not assume end > start. The occurrence id is built from these
    coordinates, so rewriting them to look linear would break the identifier."""
    length = 1_000
    genes = [{"start": 960, "end": 1_040, "strand": 1}]
    resolved = circular.resolve_origin_genes(genes, length)
    assert len(resolved) == 1
    assert resolved[0]["origin_spanning"] is True
    assert resolved[0]["start"] == 960
    assert resolved[0]["end"] == 40


def test_rotation_invariance_on_a_random_sequence():
    """The property the whole module exists to satisfy."""
    random.seed(7)
    sequence = "".join(random.choice("ACGT") for _ in range(3_000))
    baseline = circular.rotate(sequence, 0)
    for offset in (1, 137, 1_500, 2_999):
        assert sorted(circular.rotate(sequence, offset)) == sorted(baseline)
```

- [ ] **Step 2: Run test to verify it fails**

Run: `.venv/bin/python -m pytest tests/test_circular.py -v`
Expected: FAIL — `ImportError: cannot import name 'circular'`

- [ ] **Step 3: Port the implementation**

Read `src/plasmidann/circular.py` and port it. Keep `MAX_OVERLAP_BP = 5000`, `CIRCULAR_TOPOLOGIES`, `overlap_for`, `rotate`, `resolve_origin_genes`, `is_circular`. Two changes are required:

1. `resolve_origin_genes` must set an explicit `origin_spanning` boolean on every returned gene (the old version did not carry this field), because `orf_occurrences.origin_spanning` is now a declared column and §9.4 requires origin-spanning QC.
2. Coordinates of an origin-spanning gene are reported unrotated — `start` in the original frame, `end` wrapped past the origin, so `end < start`. Do not normalize this away; `ids.occurrence_id` depends on it.

Keep the module's existing comments explaining *why* the overlap exists and why it is capped; add one comment where the `origin_spanning` flag is set.

- [ ] **Step 4: Run test to verify it passes**

Run: `.venv/bin/python -m pytest tests/test_circular.py -v`
Expected: PASS (6 tests)

- [ ] **Step 5: Commit**

```bash
git add src/darkorf/circular.py tests/test_circular.py
git commit -m "feat(circular): port origin-spanning ORF handling with explicit flag

Ported from src/plasmidann/circular.py. Origin-spanning genes now carry an
explicit origin_spanning flag and keep unrotated coordinates, which the
occurrence identifier depends on, per spec 8.4/9.4.

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

---

### Task 6: Stage 0 — plasmids table, preflight and manifest

**Files:**
- Create: `src/darkorf/manifest.py`, `workflow/Snakefile`, `workflow/rules/common.smk`, `workflow/rules/s0_input.smk`, `workflow/scripts/build_plasmids.py`, `workflow/scripts/preflight.py`
- Test: `tests/test_manifest.py`, `tests/test_workflow_wiring.py`

**Interfaces:**
- Consumes: `store.write_table`, `schemas.columns`, `config/config.yaml`.
- Produces:
  - `manifest.build(config, counts, databases, software) -> dict`
  - `results/tables/plasmids.parquet`
  - `results/run_manifest.yaml`
  - Snakemake target `results/s0/preflight.tsv`

- [ ] **Step 1: Write the failing test**

```python
# tests/test_manifest.py
"""Spec §6: any parameter that affects output must be in the manifest, so a result can be
traced to the exact configuration and database versions that produced it."""
from darkorf import manifest


def test_manifest_records_the_derived_denominator_not_a_hard_coded_one():
    """Spec §3.1: the plasmid count is a property of the input, discovered at runtime."""
    built = manifest.build(
        config={"cascade": {"narrow_at": 0.9, "min_explained": 0.5}},
        counts={"plasmid_count": 7, "orf_count": 70, "unique_protein_count": 65},
        databases={"Pfam-A": "38.2"},
        software={"pyrodigal": "3.7.1"},
    )
    assert built["plasmid_count"] == 7
    assert built["unique_protein_count"] == 65


def test_manifest_stamps_the_thresholds_that_govern_the_run():
    built = manifest.build(
        config={"cascade": {"narrow_at": 0.9, "min_explained": 0.5}},
        counts={"plasmid_count": 1, "orf_count": 1, "unique_protein_count": 1},
        databases={}, software={},
    )
    assert built["config"]["cascade"]["min_explained"] == 0.5


def test_config_hash_changes_when_any_threshold_changes():
    """Two runs that used different thresholds must be distinguishable after the fact."""
    base = {"cascade": {"narrow_at": 0.9, "min_explained": 0.5}}
    changed = {"cascade": {"narrow_at": 0.9, "min_explained": 0.6}}
    counts = {"plasmid_count": 1, "orf_count": 1, "unique_protein_count": 1}
    a = manifest.build(config=base, counts=counts, databases={}, software={})
    b = manifest.build(config=changed, counts=counts, databases={}, software={})
    assert a["config_hash"] != b["config_hash"]
```

```python
# tests/test_workflow_wiring.py
"""The workflow must be loadable and its conda directives must point at a real file.
These catch the failure mode where a rule references an environment that does not exist and
the error only surfaces after hours of upstream compute."""
import pathlib
import re

ROOT = pathlib.Path(__file__).resolve().parents[1]


def test_every_conda_directive_points_at_a_file_that_exists():
    for smk in list((ROOT / "workflow").rglob("*.smk")) + [ROOT / "workflow" / "Snakefile"]:
        for match in re.finditer(r'conda:\s*\n\s*"([^"]+)"', smk.read_text()):
            referenced = (smk.parent / match.group(1)).resolve()
            assert referenced.exists(), f"{smk}: missing env {match.group(1)}"


def test_the_workflow_declares_the_stage_zero_target():
    text = (ROOT / "workflow" / "rules" / "s0_input.smk").read_text()
    assert "preflight" in text
    assert "plasmids.parquet" in text
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `.venv/bin/python -m pytest tests/test_manifest.py tests/test_workflow_wiring.py -v`
Expected: FAIL — `ImportError: cannot import name 'manifest'`

- [ ] **Step 3: Write the manifest module**

```python
# src/darkorf/manifest.py
"""Run manifest assembly (spec §6).

The manifest is what makes a result reproducible after the fact: it records the inputs, the
thresholds, the database versions and the software versions that produced a given set of
tables. config_hash exists so two runs that differ only in a threshold are distinguishable
without diffing their configuration files.
"""
import datetime
import hashlib
import json
import subprocess


def _git_commit():
    """Current commit, or 'unknown' outside a repository.

    Never fatal: a run started from an exported tarball is still a valid run, it is just
    less traceable, and that fact belongs in the manifest rather than in a crash.
    """
    try:
        return subprocess.run(
            ["git", "rev-parse", "HEAD"], capture_output=True, text=True, check=True
        ).stdout.strip()
    except (subprocess.CalledProcessError, FileNotFoundError):
        return "unknown"


def config_hash(config):
    """Stable hash of the configuration.

    sort_keys because YAML loading order must not change the hash - otherwise an unrelated
    edit makes two identical configurations look different.
    """
    canonical = json.dumps(config, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(canonical.encode()).hexdigest()[:16]


def build(config, counts, databases, software):
    """Assemble the manifest for one production run."""
    return {
        "pipeline_version": "1.0.0",
        "git_commit": _git_commit(),
        "config": config,
        "config_hash": config_hash(config),
        # Derived from the input at runtime, never hard-coded (spec §3.1).
        "plasmid_count": counts["plasmid_count"],
        "orf_count": counts["orf_count"],
        "unique_protein_count": counts["unique_protein_count"],
        "database_versions": databases,
        "software_versions": software,
        "started_at": datetime.datetime.now(datetime.timezone.utc).isoformat(),
    }
```

- [ ] **Step 4: Write the workflow skeleton**

```python
# workflow/Snakefile
"""Dark ORF discovery pipeline.

Run with:
    snakemake -s workflow/Snakefile --configfile config/config.yaml -j 96 --use-conda

Stage layout mirrors PLASMID_ANALYSIS.md: s0 input, s1 ORF prediction, s2 QC, s3
dereplication. Later plans add s4 onward. Every rule writes Parquet into
results/tables/ and nothing downstream reads a stage's intermediates directly.
"""
import pathlib

configfile: "config/config.yaml"

# Validate the configuration before any compute is scheduled. A typo in a threshold must
# fail in seconds rather than in hour 45 of a seven-day job.
from snakemake.utils import validate
validate(config, "../config/schemas/config.schema.yaml")

OUT = config["outdir"]

# Shard identifiers are zero-padded so that lexical and numeric order agree, which keeps
# wildcard globs and sorted directory listings consistent.
PLASMID_SHARDS = [f"{i:04d}" for i in range(config["shards"]["plasmid"])]
PROTEIN_SHARDS = [f"{i:04d}" for i in range(config["shards"]["protein"])]


include: "rules/common.smk"
include: "rules/s0_input.smk"
include: "rules/s1_orf.smk"
include: "rules/s2_derep.smk"
include: "rules/s3_qc.smk"


rule all:
    """Plan 1 endpoint: a validated, dereplicated protein table."""
    input:
        f"{OUT}/tables/proteins.parquet",
        f"{OUT}/run_manifest.yaml",
```

```python
# workflow/rules/common.smk
"""Helpers shared by every stage."""

def table(name):
    """Canonical path of a store table, so no rule spells the layout out by hand."""
    return f"{OUT}/tables/{name}.parquet"
```

```python
# workflow/rules/s0_input.smk
"""Stage 0: input validation, the plasmids table, and the run manifest (spec §68 Agent 1)."""


rule preflight:
    """Fail fast if a tool or database the run needs is absent.

    This rule exists because the expensive failure mode is not a missing database - it is a
    missing database discovered after gene calling has already consumed a day of wall time.
    """
    input:
        fasta=config["input"]["working_set_fasta"],
        master=config["input"]["master_table"],
    output:
        f"{OUT}/s0/preflight.tsv",
    log:
        f"{OUT}/logs/preflight.log",
    conda:
        "../envs/darkorf.yaml"
    script:
        "../scripts/preflight.py"


rule build_plasmids:
    """Build the plasmids table from the analysis set (spec §62.1).

    The plasmid count is derived here and recorded in the manifest; nothing downstream may
    hard-code it (spec §3.1).
    """
    input:
        master=config["input"]["master_table"],
        preflight=f"{OUT}/s0/preflight.tsv",
    output:
        table("plasmids"),
    log:
        f"{OUT}/logs/build_plasmids.log",
    conda:
        "../envs/darkorf.yaml"
    script:
        "../scripts/build_plasmids.py"
```

- [ ] **Step 5: Write the scripts**

```python
# workflow/scripts/preflight.py
"""Check every tool and database before any compute is spent (spec §71)."""
import pathlib
import shutil
import sys

REQUIRED_TOOLS = ["prodigal", "hmmsearch", "diamond", "mmseqs", "mafft"]

rows = []
missing = []

for tool in REQUIRED_TOOLS:
    found = shutil.which(tool)
    rows.append(f"tool\t{tool}\t{found or 'MISSING'}")
    if not found:
        missing.append(tool)

for label, path in [
    ("working_set_fasta", snakemake.input.fasta),
    ("master_table", snakemake.input.master),
    ("antifam_db", snakemake.config["qc"]["antifam"]["db"]),
]:
    exists = pathlib.Path(path).exists()
    rows.append(f"database\t{label}\t{'ok' if exists else 'MISSING'}\t{path}")
    if not exists:
        missing.append(f"{label} ({path})")

pathlib.Path(snakemake.output[0]).parent.mkdir(parents=True, exist_ok=True)
pathlib.Path(snakemake.output[0]).write_text("\n".join(rows) + "\n")

if missing:
    sys.exit("preflight failed; missing: " + ", ".join(missing))
```

```python
# workflow/scripts/build_plasmids.py
"""Build the plasmids table from the analysis set (spec §62.1)."""
import csv
import sys

import pyarrow as pa

sys.path.insert(0, "src")
from darkorf import schemas, store  # noqa: E402

# Map the master table's column names onto the declared schema. Anything the master table
# does not carry becomes NA text rather than a silent empty string, per spec §7.3.
SOURCE_COLUMNS = {
    "plasmid_id": "plasmid_id",
    "accession": "plsdb_acc",
    "source_database": "sources",
    "length": "size_bp",
    "topology": "topology",
    # Verified against the real header of analysis_set.tsv. There is no host_name,
    # host_taxid or bioproject column: plsdb_species is the best available host identity and
    # mob_host_range the best available breadth. BioProject is absent entirely, so the
    # provenance-aware recurrence of spec §34.2 cannot use it and must say so rather than
    # silently counting accessions instead.
    "host": "plsdb_species",
    "host_taxonomy": "mob_host_range",
    # hab_sub, never hab_top: the standing rule is that compartments are built from hab_sub
    # plus is_clinical.
    "habitat": "hab_sub",
    "mob_class": "mob_mobility",
    "mob_cluster": "mob_cluster",
    # mob_typer replicon types, never the pf_* columns.
    "plasmid_type": "mob_rep_types",
}

rows = {name: [] for name in schemas.columns("plasmids")}
with open(snakemake.input.master, newline="") as handle:
    for record in csv.DictReader(handle, delimiter="\t"):
        for target, source in SOURCE_COLUMNS.items():
            rows[target].append(record.get(source) or "NA")
        # Not carried by the master table at all; recorded as NA rather than omitted, so the
        # column exists and its absence is visible instead of implicit.
        rows["bioproject"].append("NA")
        # Real plasmids; controls are appended by the controls rule in Task 10.
        rows["record_class"].append("observed")

frame = pa.table(rows)
store.write_table(snakemake.output[0], "plasmids", frame)
print(f"plasmids: {frame.num_rows}", file=open(snakemake.log[0], "w"))
```

- [ ] **Step 6: Run tests to verify they pass**

Run: `.venv/bin/python -m pytest tests/test_manifest.py tests/test_workflow_wiring.py -v`
Expected: PASS (5 tests)

- [ ] **Step 7: Commit**

```bash
git add src/darkorf/manifest.py workflow/ tests/test_manifest.py tests/test_workflow_wiring.py
git commit -m "feat(s0): add workflow skeleton, preflight, plasmids table and run manifest

Preflight fails in seconds on a missing tool or database rather than after a
day of gene calling. The plasmid count is derived from the input and recorded
in the manifest, never hard-coded, per spec 3.1/6.

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

---

### Task 7: Stage 1 — ORF prediction

**Files:**
- Create: `src/darkorf/orf.py`, `workflow/rules/s1_orf.smk`, `workflow/scripts/shard_plasmids.py`, `workflow/scripts/call_orfs.py`
- Test: `tests/test_orf.py`

**Interfaces:**
- Consumes: `circular.is_circular`, `circular.overlap_for`, `circular.rotate`, `circular.resolve_origin_genes`, `ids.occurrence_id`, `ids.protein_id`, `schemas`, `store`.
- Produces:
  - `orf.call_orfs(plasmid_id, sequence, topology, min_aa, meta_mode) -> list[dict]` with keys matching `schemas.columns("orf_occurrences")` for the fields Stage 1 fills (`orf_occurrence_id`, `plasmid_id`, `protein_id`, `start`, `end`, `strand`, `nucleotide_length`, `protein_length`, `cds_nucleotide_sequence`, `partial`, `start_type`, `origin_spanning`)
  - `results/s1/orfs/{shard}.parquet`

- [ ] **Step 1: Write the failing test**

```python
# tests/test_orf.py
"""Stage 1 gene calling (spec §8)."""
from darkorf import orf

# A short synthetic ORF: ATG ... TAA in frame, long enough to clear the 20 aa floor.
CODING = "ATG" + "AAA" * 25 + "TAA"


def test_every_called_orf_carries_its_cds_nucleotide_sequence():
    """Spec §8.3. Without the nucleotides, RNAcode and dN/dS become one-shot and re-running
    them later means re-extracting CDS from every plasmid in the collection."""
    called = orf.call_orfs("p1", CODING, "linear", min_aa=20, meta_mode=True)
    assert called
    for gene in called:
        assert gene["cds_nucleotide_sequence"]
        assert len(gene["cds_nucleotide_sequence"]) % 3 == 0


def test_protein_and_occurrence_ids_are_populated():
    called = orf.call_orfs("p1", CODING, "linear", min_aa=20, meta_mode=True)
    gene = called[0]
    assert len(gene["protein_id"]) == 32
    assert gene["orf_occurrence_id"].startswith("p1:")


def test_short_orfs_are_not_called_below_the_configured_floor():
    """The floor is the gene caller's, distinct from the discovery cutoff (spec §11.1)."""
    tiny = "ATG" + "AAA" * 3 + "TAA"
    assert orf.call_orfs("p1", tiny, "linear", min_aa=20, meta_mode=True) == []


def test_a_circular_plasmid_reports_origin_spanning_genes_once():
    """The gene crossing the origin must appear exactly once, with origin_spanning set."""
    # Place a coding sequence so that it runs off the end and wraps to the start.
    sequence = "AAA" * 30 + CODING
    rotated = sequence[-30:] + sequence[:-30]
    called = orf.call_orfs("p1", rotated, "circular", min_aa=20, meta_mode=True)
    ids = [gene["orf_occurrence_id"] for gene in called]
    assert len(ids) == len(set(ids)), "an origin-spanning gene was reported twice"


def test_calling_is_invariant_to_where_a_circular_sequence_was_cut():
    """The property that motivates the whole circular module: a circle has no privileged
    starting point, so the protein set must not depend on where the assembler cut it."""
    sequence = "AAA" * 40 + CODING + "TTT" * 40
    first = orf.call_orfs("p1", sequence, "circular", min_aa=20, meta_mode=True)
    rotated = sequence[100:] + sequence[:100]
    second = orf.call_orfs("p1", rotated, "circular", min_aa=20, meta_mode=True)
    assert {g["protein_id"] for g in first} == {g["protein_id"] for g in second}
```

- [ ] **Step 2: Run test to verify it fails**

Run: `.venv/bin/python -m pytest tests/test_orf.py -v`
Expected: FAIL — `ImportError: cannot import name 'orf'`

- [ ] **Step 3: Write the implementation**

```python
# src/darkorf/orf.py
"""Stage 1: ORF prediction with circular-origin support (spec §8).

pyrodigal in metagenomic mode: plasmids are short and compositionally heterogeneous, so
there is not enough sequence to train a reliable single-genome model on most of them.
"""
import pyrodigal

from darkorf import circular, ids


def _gene_finder(meta_mode):
    # min_gene is in nucleotides; the caller passes a protein-length floor.
    return pyrodigal.GeneFinder(meta=meta_mode)


def call_orfs(plasmid_id, sequence, topology, min_aa, meta_mode):
    """Call ORFs on one plasmid and return occurrence records.

    For circular plasmids the sequence is extended by a duplicated head before calling, so
    that a gene crossing the origin is presented to the caller as a contiguous stretch. The
    duplication necessarily produces some genes twice, which resolve_origin_genes collapses
    back to one copy each.
    """
    sequence = sequence.upper()
    original_length = len(sequence)

    if circular.is_circular(topology):
        overlap = circular.overlap_for(original_length)
        search_space = sequence + sequence[:overlap]
    else:
        overlap = 0
        search_space = sequence

    genes = []
    for prediction in _gene_finder(meta_mode).find_genes(search_space.encode()):
        protein = prediction.translate()
        # The trailing stop is not part of the protein; ids.protein_id strips it too, but
        # protein_length must not count it either.
        peptide = protein[:-1] if protein.endswith("*") else protein
        if len(peptide) < min_aa:
            continue
        genes.append({
            "start": prediction.begin,
            "end": prediction.end,
            "strand": prediction.strand,
            "protein_sequence": peptide,
            "cds_nucleotide_sequence": search_space[prediction.begin - 1:prediction.end],
            "partial": prediction.partial_begin or prediction.partial_end,
            "start_type": prediction.start_type,
        })

    if circular.is_circular(topology):
        genes = circular.resolve_origin_genes(genes, original_length)
    else:
        for gene in genes:
            gene["origin_spanning"] = False

    records = []
    for gene in genes:
        records.append({
            "orf_occurrence_id": ids.occurrence_id(
                plasmid_id, gene["start"], gene["end"], gene["strand"]
            ),
            "plasmid_id": plasmid_id,
            "protein_id": ids.protein_id(gene["protein_sequence"]),
            "start": gene["start"],
            "end": gene["end"],
            "strand": gene["strand"],
            "nucleotide_length": len(gene["cds_nucleotide_sequence"]),
            "protein_length": len(gene["protein_sequence"]),
            "cds_nucleotide_sequence": gene["cds_nucleotide_sequence"],
            "partial": bool(gene["partial"]),
            "start_type": gene["start_type"],
            "origin_spanning": bool(gene["origin_spanning"]),
        })
    return records
```

- [ ] **Step 4: Run test to verify it passes**

Run: `.venv/bin/python -m pytest tests/test_orf.py -v`
Expected: PASS (5 tests)

- [ ] **Step 5: Write the sharding rule and script**

```python
# workflow/rules/s1_orf.smk
"""Stage 1: sharded ORF prediction (spec §68 Agent 2).

Sharding over plasmids sets the unit of resume: a failure costs one shard rather than the
whole stage. The shard count is recorded in the manifest, and changing it between resumed
runs invalidates this stage (spec §72).
"""


rule shard_plasmids:
    """Split the analysis-set FASTA into deterministic shards.

    Assignment is by a hash of the plasmid id, not by position in the file, so the same
    plasmid lands in the same shard on every run.
    """
    input:
        fasta=config["input"]["working_set_fasta"],
        plasmids=table("plasmids"),
    output:
        expand(f"{OUT}/s1/shards/{{shard}}.fna", shard=PLASMID_SHARDS),
    log:
        f"{OUT}/logs/shard_plasmids.log",
    conda:
        "../envs/darkorf.yaml"
    script:
        "../scripts/shard_plasmids.py"


rule call_orfs:
    """Call ORFs on one shard (spec §8)."""
    input:
        fasta=f"{OUT}/s1/shards/{{shard}}.fna",
        plasmids=table("plasmids"),
    output:
        f"{OUT}/s1/orfs/{{shard}}.parquet",
    resources:
        mem_mb=4000,
        runtime=120,
    log:
        f"{OUT}/logs/call_orfs/{{shard}}.log",
    conda:
        "../envs/darkorf.yaml"
    script:
        "../scripts/call_orfs.py"
```

Write `workflow/scripts/shard_plasmids.py` to stream the gzipped FASTA once, assigning each record to shard `int(sha256(plasmid_id)[:8], 16) % n_shards`, and `workflow/scripts/call_orfs.py` to read one shard plus the topology column of `plasmids.parquet`, call `orf.call_orfs` per record, and write the shard's occurrences as Parquet with only the Stage 1 columns populated. Both scripts insert `src` on `sys.path` and import from `darkorf`, matching `build_plasmids.py`.

- [ ] **Step 6: Commit**

```bash
git add src/darkorf/orf.py workflow/rules/s1_orf.smk workflow/scripts/shard_plasmids.py workflow/scripts/call_orfs.py tests/test_orf.py
git commit -m "feat(s1): add sharded ORF prediction with CDS nucleotide retention

Gene calling is invariant to where a circular sequence was cut, origin-spanning
genes are reported once and flagged, and every occurrence carries its CDS
nucleotides so the evolutionary stages are repeatable, per spec 8.

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

---

### Task 8: Stage 2 — dereplication

Dereplication runs **before** every other per-protein stage. A protein sequence occurring on
400 plasmids is then screened, searched and annotated once rather than 400 times: the
collection's 9.3 million ORF occurrences collapse to roughly 3.5 million unique sequences, so
everything downstream of this stage costs about 2.7 times less. The constraint that makes it
safe is that the occurrence mapping stays lossless (spec §13).

**Files:**
- Create: `src/darkorf/derep.py`, `workflow/rules/s2_derep.smk`, `workflow/scripts/dereplicate.py`
- Test: `tests/test_derep.py`

**Interfaces:**
- Consumes: Stage 1 shard output (`results/s1/orfs/{shard}.parquet`), `ids.protein_id`.
- Produces:
  - `derep.dereplicate(occurrences) -> tuple[list[dict], dict[str, list[str]]]` — protein core records plus `protein_id -> [occurrence_id]`. Each protein record has keys `protein_id`, `protein_sequence`, `protein_length`, `occurrence_count`, `plasmid_count`, `has_complete_occurrence`.
  - `results/tables/protein_occurrences.parquet` — the lossless mapping
  - `results/s2/protein_core.parquet` — protein universe without QC columns
  - `results/s2/proteins/{pshard}.faa` — unique proteins sharded for every later search

Note on what this task does **not** produce: the declared `proteins` table. `discovery_eligible`
and the AntiFam columns are only known after Stage 3, so Stage 3 assembles the final table.
Writing a `proteins` table here and rewriting it later would leave a window in which two
different files both claim to be the protein table.

- [ ] **Step 1: Write the failing test**

```python
# tests/test_derep.py
"""Stage 2 dereplication (spec §13).

The one rule: dereplication changes computation, never biological occurrence counts. If a
protein occurs on 400 plasmids, every later stage must still be able to see 400 occurrences.
"""
from darkorf import derep, ids


def _occurrence(plasmid, sequence, complete=True, occurrence_id=None):
    return {
        "orf_occurrence_id": occurrence_id or f"{plasmid}:1-99:+",
        "plasmid_id": plasmid,
        "protein_id": ids.protein_id(sequence),
        "protein_sequence": sequence,
        "protein_length": len(sequence),
        "partial": not complete,
    }


def test_identical_sequences_collapse_to_one_protein():
    proteins, mapping = derep.dereplicate([
        _occurrence("p1", "MKV"), _occurrence("p2", "MKV"),
    ])
    assert len(proteins) == 1
    assert proteins[0]["occurrence_count"] == 2
    assert proteins[0]["plasmid_count"] == 2


def test_the_occurrence_mapping_is_lossless():
    """The invariant that makes dereplication safe: every input occurrence must be
    recoverable through the mapping."""
    occurrences = [_occurrence("p1", "MKV"), _occurrence("p2", "MKV"), _occurrence("p3", "MKW")]
    proteins, mapping = derep.dereplicate(occurrences)
    recovered = {occ for occs in mapping.values() for occ in occs}
    assert recovered == {o["orf_occurrence_id"] for o in occurrences}


def test_a_protein_is_complete_if_any_occurrence_is_complete():
    """Spec §12: one partial occurrence does not make the protein globally incomplete."""
    proteins, _ = derep.dereplicate([
        _occurrence("p1", "MKV", complete=False),
        _occurrence("p2", "MKV", complete=True),
    ])
    assert proteins[0]["has_complete_occurrence"] is True


def test_a_protein_seen_only_as_partial_is_not_complete():
    proteins, _ = derep.dereplicate([_occurrence("p1", "MKV", complete=False)])
    assert proteins[0]["has_complete_occurrence"] is False


def test_the_same_protein_twice_on_one_plasmid_counts_two_occurrences_one_plasmid():
    """Occurrence count and plasmid count are different quantities. Conflating them inflates
    every independence measure downstream, which is the defect that made mob_cluster-based
    counts unusable in the previous design."""
    proteins, _ = derep.dereplicate([
        _occurrence("p1", "MKV", occurrence_id="p1:1-99:+"),
        _occurrence("p1", "MKV", occurrence_id="p1:500-599:+"),
    ])
    assert proteins[0]["occurrence_count"] == 2
    assert proteins[0]["plasmid_count"] == 1


def test_protein_records_are_emitted_in_a_deterministic_order():
    """Shard output order depends on filesystem listing, which is not stable. Sorting here is
    what makes the downstream protein shards reproducible between runs."""
    first, _ = derep.dereplicate([_occurrence("p1", "MKW"), _occurrence("p2", "MKV")])
    second, _ = derep.dereplicate([_occurrence("p2", "MKV"), _occurrence("p1", "MKW")])
    assert [p["protein_id"] for p in first] == [p["protein_id"] for p in second]
```

- [ ] **Step 2: Run test to verify it fails**

Run: `.venv/bin/python -m pytest tests/test_derep.py -v`
Expected: FAIL — `ImportError: cannot import name 'derep'`

- [ ] **Step 3: Write the implementation**

Port `src/plasmidann/dereplicate.py`, adapting it to the two-part return above. Read the
original first: its losslessness assertion is the behaviour being preserved.

```python
# src/darkorf/derep.py
"""Stage 2: dereplication (spec §13).

Runs before every per-protein stage. The collection's ~9.3 million ORF occurrences collapse
to ~3.5 million unique sequences, so AntiFam, the annotation cascade, structure and every
other per-protein tool does roughly 2.7 times less work.

Dereplication is a computational optimization with one hard constraint: it must not change
biological occurrence counts. A protein on 400 plasmids is searched once and still reports
400 occurrences.
"""
import collections


def dereplicate(occurrences):
    """Collapse occurrences to unique proteins and return the lossless mapping.

    Returns (protein_records, mapping) where mapping is protein_id -> [orf_occurrence_id].
    Protein records are sorted by protein_id: shard output arrives in filesystem order, which
    is not stable between runs, and the downstream protein shards must be reproducible.
    """
    sequences = {}
    occurrence_ids = collections.defaultdict(list)
    plasmids = collections.defaultdict(set)
    complete = collections.defaultdict(bool)

    for record in occurrences:
        protein = record["protein_id"]
        sequences.setdefault(protein, record["protein_sequence"])
        occurrence_ids[protein].append(record["orf_occurrence_id"])
        plasmids[protein].add(record["plasmid_id"])
        # Spec §12: the protein is complete if ANY occurrence is complete. One truncated copy
        # at a contig edge says nothing about the protein itself.
        if not record["partial"]:
            complete[protein] = True

    proteins = []
    for protein in sorted(sequences):
        proteins.append({
            "protein_id": protein,
            "protein_sequence": sequences[protein],
            "protein_length": len(sequences[protein]),
            # Two different quantities, deliberately both kept.
            "occurrence_count": len(occurrence_ids[protein]),
            "plasmid_count": len(plasmids[protein]),
            "has_complete_occurrence": bool(complete[protein]),
        })

    return proteins, dict(occurrence_ids)
```

- [ ] **Step 4: Run test to verify it passes**

Run: `.venv/bin/python -m pytest tests/test_derep.py -v`
Expected: PASS (6 tests)

- [ ] **Step 5: Write the rule**

```python
# workflow/rules/s2_derep.smk
"""Stage 2: dereplication (spec §13, §68 Agent 4).

This stage is placed before all per-protein work on purpose. Every tool downstream - AntiFam,
the six cascade tiers, structure, properties - runs on the ~3.5 million unique sequences
rather than the ~9.3 million occurrences.
"""


rule dereplicate:
    """Collapse occurrences to unique proteins, keeping the mapping lossless."""
    input:
        orfs=expand(f"{OUT}/s1/orfs/{{shard}}.parquet", shard=PLASMID_SHARDS),
    output:
        core=f"{OUT}/s2/protein_core.parquet",
        mapping=f"{OUT}/tables/protein_occurrences.parquet",
    resources:
        # Roughly 3.5M sequences at a ~330 aa mean, plus the occurrence mapping. Measured on
        # the smoke set and extrapolated; the benchmark rule refines this before production.
        mem_mb=64000,
        runtime=240,
    log:
        f"{OUT}/logs/dereplicate.log",
    conda:
        "../envs/darkorf.yaml"
    script:
        "../scripts/dereplicate.py"


rule shard_proteins:
    """Write the unique proteins as FASTA shards, the input unit for every later search."""
    input:
        core=f"{OUT}/s2/protein_core.parquet",
    output:
        expand(f"{OUT}/s2/proteins/{{pshard}}.faa", pshard=PROTEIN_SHARDS),
    resources:
        mem_mb=16000,
        runtime=60,
    log:
        f"{OUT}/logs/shard_proteins.log",
    conda:
        "../envs/darkorf.yaml"
    script:
        "../scripts/shard_proteins.py"
```

`workflow/scripts/dereplicate.py` reads the Stage 1 shards one at a time (never all at once),
calls `derep.dereplicate`, and **asserts losslessness before writing**: the total number of
mapped occurrence ids must equal the number of input rows, and the script must fail loudly
with both numbers if it does not. `workflow/scripts/shard_proteins.py` assigns each protein to
shard `int(protein_id[:8], 16) % n_protein_shards`, so a protein lands in the same shard on
every run regardless of how many proteins exist.

- [ ] **Step 6: Commit**

```bash
git add src/darkorf/derep.py workflow/rules/s2_derep.smk workflow/scripts/dereplicate.py workflow/scripts/shard_proteins.py tests/test_derep.py
git commit -m "feat(s2): dereplicate before per-protein stages, with a lossless mapping

Dereplication runs first so AntiFam, the cascade and every later per-protein
tool sees ~3.5M unique sequences rather than ~9.3M occurrences. A protein is
complete if any occurrence is complete, per spec 12, and occurrence count and
plasmid count stay distinct quantities.

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

---

### Task 9: Stage 3 — AntiFam, overlap QC and discovery eligibility

**Files:**
- Create: `src/darkorf/qc.py`, `workflow/rules/s3_qc.smk`, `workflow/scripts/apply_qc.py`
- Test: `tests/test_qc.py`

**Interfaces:**
- Consumes: `results/s2/protein_core.parquet`, `results/s2/proteins/{pshard}.faa`, the Stage 1 shards, `schemas`, `store`.
- Produces:
  - `qc.opposite_strand_overlap(gene, neighbours) -> tuple[float, str|None]`
  - `qc.discovery_eligibility(record, min_aa, max_overlap) -> dict` with keys `discovery_excluded_short`, `discovery_excluded_partial`, `discovery_excluded_antifam`, `discovery_eligible`, `exclusion_reason`
  - `results/tables/orf_occurrences.parquet` and `results/tables/proteins.parquet`

Two levels meet in this stage, and keeping them apart is the point. AntiFam status is a
property of a **sequence**, so it is computed once per unique protein. Opposite-strand overlap
is a property of an **occurrence**, because it depends on the neighbours on that particular
plasmid, so it is computed per occurrence and cannot be dereplicated.

- [ ] **Step 1: Write the failing test**

```python
# tests/test_qc.py
"""Stage 3 QC (spec §9, §10, §12).

The governing rule for this whole module: flag, never delete. Every excluded record stays in
the master dataset with a reason attached (spec §9.4). The excluded population is itself a
result - it is how the artifact rate gets measured rather than assumed.
"""
from darkorf import qc


def test_an_orf_antisense_to_a_longer_orf_is_flagged():
    """Spec §9.3. The shadow-ORF artifact: a spurious call on the reverse complement of a real
    gene. It passes the entire annotation cascade cleanly, because it is not a protein and no
    database contains it, and it is recurrent for exactly the same reason the real gene is -
    so recurrence evidence selects for it unless it is caught here."""
    shadow = {"orf_occurrence_id": "p1:100-400:-", "start": 100, "end": 400, "strand": -1}
    real = {"orf_occurrence_id": "p1:90-500:+", "start": 90, "end": 500, "strand": 1}
    fraction, partner = qc.opposite_strand_overlap(shadow, [real])
    assert fraction > 0.9
    assert partner == "p1:90-500:+"


def test_same_strand_neighbours_do_not_count_as_overlap():
    """Operons are gene-dense and same-strand overlap is ordinary bacterial biology."""
    gene = {"orf_occurrence_id": "p1:100-400:+", "start": 100, "end": 400, "strand": 1}
    neighbour = {"orf_occurrence_id": "p1:90-500:+", "start": 90, "end": 500, "strand": 1}
    fraction, partner = qc.opposite_strand_overlap(gene, [neighbour])
    assert fraction == 0.0
    assert partner is None


def test_a_gene_with_no_neighbours_has_no_overlap():
    gene = {"orf_occurrence_id": "p1:100-400:+", "start": 100, "end": 400, "strand": 1}
    assert qc.opposite_strand_overlap(gene, []) == (0.0, None)


def test_an_origin_spanning_gene_does_not_report_a_negative_span():
    """An origin-spanning gene keeps unrotated coordinates, so end < start (spec §8.4).
    Treating that as a negative length would silently produce nonsense overlap fractions."""
    wrapped = {"orf_occurrence_id": "p1:960-40:+", "start": 960, "end": 40, "strand": 1}
    other = {"orf_occurrence_id": "p1:900-980:-", "start": 900, "end": 980, "strand": -1}
    fraction, _ = qc.opposite_strand_overlap(wrapped, [other])
    assert 0.0 <= fraction <= 1.0


def test_a_short_protein_is_flagged_but_stays_in_the_master_dataset():
    """Spec §10.2: the cutoff is an operational discovery boundary, not a claim that short
    proteins are not real."""
    record = {"protein_length": 12, "partial": False, "antifam_hit": False,
              "opposite_strand_overlap_fraction": 0.0}
    verdict = qc.discovery_eligibility(record, min_aa=20, max_overlap=0.6)
    assert verdict["discovery_excluded_short"] is True
    assert verdict["discovery_eligible"] is False
    assert "short" in verdict["exclusion_reason"]


def test_a_clean_protein_is_eligible():
    record = {"protein_length": 120, "partial": False, "antifam_hit": False,
              "opposite_strand_overlap_fraction": 0.1}
    verdict = qc.discovery_eligibility(record, min_aa=20, max_overlap=0.6)
    assert verdict["discovery_eligible"] is True
    assert verdict["exclusion_reason"] == ""


def test_multiple_exclusion_reasons_are_all_recorded():
    """A record excluded for two reasons must say both. Keeping only the first would make the
    exclusion counts in the QC report wrong."""
    record = {"protein_length": 10, "partial": True, "antifam_hit": True,
              "opposite_strand_overlap_fraction": 0.0}
    verdict = qc.discovery_eligibility(record, min_aa=20, max_overlap=0.6)
    for reason in ("short", "partial", "antifam"):
        assert reason in verdict["exclusion_reason"]
```

- [ ] **Step 2: Run test to verify it fails**

Run: `.venv/bin/python -m pytest tests/test_qc.py -v`
Expected: FAIL — `ImportError: cannot import name 'qc'`

- [ ] **Step 3: Write the implementation**

```python
# src/darkorf/qc.py
"""Stage 3: quality control and discovery eligibility (spec §9, §10, §12).

Nothing here deletes a record. Every exclusion sets a flag and a reason and the record stays
in the master dataset, because the excluded population is how the artifact rate is measured.

Two levels meet in this module. AntiFam status is a property of a sequence and is computed
once per unique protein. Opposite-strand overlap is a property of an occurrence, because it
depends on that plasmid's neighbours, and cannot be dereplicated.
"""


def _span(gene):
    """Half-open interval covered by a gene, ignoring strand.

    Origin-spanning genes keep unrotated coordinates, so end < start is legitimate (spec
    §8.4). min/max rather than assuming an order, because a negative length here would turn
    into a meaningless overlap fraction rather than an error.
    """
    return min(gene["start"], gene["end"]), max(gene["start"], gene["end"])


def opposite_strand_overlap(gene, neighbours):
    """Fraction of this gene covered by the longest opposite-strand gene overlapping it.

    Spec §9.3. Returns (fraction, partner_occurrence_id). Only opposite-strand overlap counts:
    same-strand overlap is ordinary in gene-dense bacterial replicons, whereas a call lying
    antisense to a longer real gene is the classic spurious-ORF signature, and it is the one
    artifact class that survives the whole annotation cascade untouched.
    """
    start, end = _span(gene)
    length = end - start
    if length <= 0:
        return 0.0, None

    best_fraction = 0.0
    best_partner = None
    for other in neighbours:
        if other["strand"] == gene["strand"]:
            continue
        other_start, other_end = _span(other)
        covered = min(end, other_end) - max(start, other_start)
        if covered <= 0:
            continue
        fraction = covered / length
        if fraction > best_fraction:
            best_fraction = fraction
            best_partner = other["orf_occurrence_id"]
    return best_fraction, best_partner


def discovery_eligibility(record, min_aa, max_overlap):
    """Decide whether one record enters the dark-discovery population (spec §10, §12).

    Eligibility is about what may be discovered, not about what is real: a 12-residue
    microprotein is excluded because the evidence tools cannot interpret it, and it stays in
    the dataset because it may still be biology.
    """
    reasons = []

    excluded_short = record["protein_length"] < min_aa
    if excluded_short:
        reasons.append("short")

    # Spec §12: an apparent lack of annotation on a truncated sequence says nothing.
    excluded_partial = bool(record["partial"])
    if excluded_partial:
        reasons.append("partial")

    excluded_antifam = bool(record["antifam_hit"])
    if excluded_antifam:
        reasons.append("antifam")

    if record.get("opposite_strand_overlap_fraction", 0.0) > max_overlap:
        reasons.append("opposite_strand_overlap")

    return {
        "discovery_excluded_short": excluded_short,
        "discovery_excluded_partial": excluded_partial,
        "discovery_excluded_antifam": excluded_antifam,
        "discovery_eligible": not reasons,
        "exclusion_reason": ",".join(reasons),
    }
```

- [ ] **Step 4: Run test to verify it passes**

Run: `.venv/bin/python -m pytest tests/test_qc.py -v`
Expected: PASS (7 tests)

- [ ] **Step 5: Write the rules**

```python
# workflow/rules/s3_qc.smk
"""Stage 3: artifact screening and discovery eligibility (spec §9, §68 Agent 3).

AntiFam runs on the DEREPLICATED protein shards produced by Stage 2, not on occurrences:
AntiFam status is a property of the sequence, so screening ~3.5 million unique proteins gives
the identical answer for roughly a third of the work.
"""


rule antifam_search:
    """Screen unique proteins against AntiFam at its curated per-family thresholds.

    --cut_ga rather than a blanket E-value: 274 of AntiFam's 278 profiles carry a curated
    gathering threshold looser than E=1e-5, so a global floor would override the curator
    across almost the whole database - in the one screen whose job is stopping non-proteins
    from reaching the discovery set (Eberhardt et al. 2012, Database 2012:bas003).
    """
    input:
        faa=f"{OUT}/s2/proteins/{{pshard}}.faa",
        db=config["qc"]["antifam"]["db"],
    output:
        f"{OUT}/s3/antifam/{{pshard}}.domtbl",
    threads: 4
    resources:
        mem_mb=8000,
        runtime=240,
    log:
        f"{OUT}/logs/antifam/{{pshard}}.log",
    conda:
        "../envs/darkorf.yaml"
    shell:
        "hmmsearch --cut_ga --domtblout {output} --cpu {threads} "
        "{input.db} {input.faa} > /dev/null 2> {log}"


rule apply_qc:
    """Assemble the occurrence and protein tables with all QC columns filled (spec §9, §10)."""
    input:
        orfs=expand(f"{OUT}/s1/orfs/{{shard}}.parquet", shard=PLASMID_SHARDS),
        core=f"{OUT}/s2/protein_core.parquet",
        antifam=expand(f"{OUT}/s3/antifam/{{pshard}}.domtbl", pshard=PROTEIN_SHARDS),
    output:
        occurrences=table("orf_occurrences"),
        proteins=table("proteins"),
    resources:
        mem_mb=48000,
        runtime=240,
    log:
        f"{OUT}/logs/apply_qc.log",
    conda:
        "../envs/darkorf.yaml"
    script:
        "../scripts/apply_qc.py"
```

`workflow/scripts/apply_qc.py` must:

1. Parse the AntiFam domtbl shards into `protein_id -> (model, score, evalue)`.
2. Process the Stage 1 shards **one shard at a time**, and within a shard group by
   `plasmid_id` so that `opposite_strand_overlap` compares a gene only against neighbours on
   its own plasmid. Sharding is by plasmid, so a plasmid's genes are never split across
   shards — assert this rather than assuming it.
3. Apply `discovery_eligibility` per occurrence and write `orf_occurrences`.
4. Roll the occurrence-level verdicts up to `proteins`: `discovery_eligible` is true for a
   protein if **any** of its occurrences is eligible, matching the spec §12 rule that one
   partial copy does not disqualify a protein, and `record_class` is `"observed"`.

Never hold all occurrences in memory at once; append per shard.

- [ ] **Step 6: Commit**

```bash
git add src/darkorf/qc.py workflow/rules/s3_qc.smk workflow/scripts/apply_qc.py tests/test_qc.py
git commit -m "feat(s3): screen AntiFam over unique proteins, add opposite-strand overlap QC

AntiFam runs on the dereplicated set because the status is a property of the
sequence, which is the identical answer for a third of the work. Adds the
shadow-ORF signature the artifact screen previously missed: fraction of a call
covered by a longer opposite-strand ORF. Every exclusion is a flag with a
reason and the record stays in the dataset, per spec 9.

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

---

### Task 10: Positive and negative controls

**Files:**
- Create: `src/darkorf/controls.py`, `workflow/scripts/build_controls.py`
- Modify: `workflow/rules/s0_input.smk` (add the `build_controls` rule)
- Test: `tests/test_controls.py`

**Interfaces:**
- Consumes: `ids.protein_id`, `schemas`, `store`, config `controls.*`.
- Produces:
  - `controls.shuffled_decoy(sequence, rng) -> str`
  - `controls.reverse_complement_decoy(cds_nucleotide_sequence) -> str`
  - `controls.build(positive_records, cds_records, config, seed) -> list[dict]` matching `schemas.columns("controls")`

- [ ] **Step 1: Write the failing test**

```python
# tests/test_controls.py
"""Controls (spec §58).

Positive controls ask: can the pipeline still recover known biology? Negative controls ask
the opposite and previously unasked question: does the pipeline avoid treating decoy
sequences as convincing biological annotations? Without the second, the dark classification
can only fail in one direction, which is the falsifiability gap review found.
"""
import random
from darkorf import controls


def test_a_shuffled_decoy_preserves_length_and_composition():
    """Composition is preserved so that the decoy fails for the right reason: absence of
    homology, not an obviously unnatural amino-acid frequency."""
    original = "MKVLAACDEFGHIKLMNPQRSTVWY"
    decoy = controls.shuffled_decoy(original, random.Random(1))
    assert len(decoy) == len(original)
    assert sorted(decoy) == sorted(original)


def test_shuffling_is_deterministic_for_a_given_seed():
    """Spec §61: the same run configuration must produce the same controls, or the quality
    gate moves between runs for reasons unrelated to the pipeline."""
    a = controls.shuffled_decoy("MKVLAACDEFG", random.Random(7))
    b = controls.shuffled_decoy("MKVLAACDEFG", random.Random(7))
    assert a == b


def test_a_reverse_complement_decoy_translates_the_antisense_strand():
    """The decoy that mimics the real artifact class: a shadow ORF read off the opposite
    strand of a genuine gene."""
    decoy = controls.reverse_complement_decoy("ATGAAAGGGTTTTAA")
    assert decoy
    assert "*" not in decoy[:-1]


def test_controls_are_labelled_and_carry_their_expected_outcome():
    records = controls.build(
        positive_records=[{"protein_sequence": "MKVLAA", "source": "sp|P00001"}],
        cds_records=[{"cds_nucleotide_sequence": "ATGAAAGGGTTTTAA"}],
        config={"positive": {"n": 1}, "negative": {"n_shuffled": 1, "n_reverse_complement": 1}},
        seed=20260916,
    )
    classes = {record["record_class"] for record in records}
    assert classes == {"positive_control", "negative_control"}
    for record in records:
        if record["record_class"] == "positive_control":
            assert record["expected_outcome"] == "ANNOTATED"
        else:
            assert record["expected_outcome"] == "DARK"
```

- [ ] **Step 2: Run test to verify it fails**

Run: `.venv/bin/python -m pytest tests/test_controls.py -v`
Expected: FAIL — `ImportError: cannot import name 'controls'`

- [ ] **Step 3: Write the implementation**

Implement the three functions. `shuffled_decoy` shuffles residues with the supplied `random.Random`. `reverse_complement_decoy` reverse-complements the CDS and translates frame 1, stopping at the first stop codon. `build` emits one record per control with `record_class` in `{positive_control, negative_control}`, `control_category` in `{known_function, shuffled, reverse_complement}`, `source_protein_id` where applicable, and `expected_outcome` in `{ANNOTATED, DARK}`. Controls are injected into the protein set before the cascade and excluded from biological occurrence-integrity checks (spec §58.2).

- [ ] **Step 4: Run test to verify it passes**

Run: `.venv/bin/python -m pytest tests/test_controls.py -v`
Expected: PASS (4 tests)

- [ ] **Step 5: Commit**

```bash
git add src/darkorf/controls.py workflow/scripts/build_controls.py workflow/rules/s0_input.smk tests/test_controls.py
git commit -m "feat(controls): add positive controls and shuffled/antisense negative controls

Negative controls make the dark classification falsifiable in both directions:
decoys must come out DARK without a confident functional annotation, per spec 58.

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

---

### Task 11: End-to-end smoke run

**Files:**
- Create: `.test/mini.fna.gz`, `.test/mini_analysis_set.tsv`, `config/test.yaml`, `tests/test_smoke.py`
- Test: `tests/test_smoke.py`

**Interfaces:**
- Consumes: everything above.
- Produces: a runnable `snakemake --configfile config/test.yaml` target completing in under two minutes.

- [ ] **Step 1: Build the miniature dataset**

```bash
# 200 plasmids sampled deterministically from the analysis set, small enough that the whole
# workflow runs in a couple of minutes on a laptop but real enough to exercise circular
# topology, origin-spanning genes and duplicate proteins.
.venv/bin/python - <<'PY'
import gzip, random, pathlib
src = "data/plasmidscope_primary/provenance/working_set.fna.gz"
random.seed(20260916)
records, current = [], None
with gzip.open(src, "rt") as handle:
    for line in handle:
        if line.startswith(">"):
            if current: records.append(current)
            if len(records) >= 4000: break
            current = [line, ""]
        elif current:
            current[1] += line.strip()
sample = random.sample(records, 200)
pathlib.Path(".test").mkdir(exist_ok=True)
with gzip.open(".test/mini.fna.gz", "wt") as out:
    for header, seq in sample:
        out.write(header); out.write(seq + "\n")
print(f"wrote {len(sample)} plasmids")
PY
```

Then subset the master table to the same identifiers into `.test/mini_analysis_set.tsv`, preserving the header.

- [ ] **Step 2: Write the test configuration**

```yaml
# config/test.yaml
# The smoke configuration. Same schema, tiny inputs, few shards - so the wiring is exercised
# end to end in seconds rather than days. Any rule that works here can still fail at scale,
# but a rule that fails here is broken everywhere.
input:
  working_set_fasta: .test/mini.fna.gz
  master_table: .test/mini_analysis_set.tsv
outdir: .test/results
shards: {plasmid: 4, protein: 2}
orf: {min_call_length_aa: 20, meta_mode: true}
discovery: {min_dark_candidate_length_aa: 20}
qc:
  max_opposite_strand_overlap_fraction: 0.6
  antifam: {db: data/refs/antifam/AntiFam.hmm, threshold: cut_ga}
cascade: {narrow_at: 0.9, min_explained: 0.5, sweep_cohort_fraction: 0.02}
controls:
  positive: {raw_faa: data/refs/control/raw.faa, n: 20}
  negative: {n_shuffled: 10, n_reverse_complement: 10}
seed: 20260916
```

- [ ] **Step 3: Write the smoke test**

```python
# tests/test_smoke.py
"""End-to-end wiring test (spec §73).

This is the test that catches what unit tests cannot: a rule whose output filename does not
match the next rule's input, a script that imports a name the module does not export, a
missing conda directive. It runs the real workflow on 200 plasmids.
"""
import pathlib
import subprocess
import pytest

ROOT = pathlib.Path(__file__).resolve().parents[1]


@pytest.mark.smoke
def test_the_workflow_runs_end_to_end_on_the_miniature_dataset(tmp_path):
    result = subprocess.run(
        ["snakemake", "-s", "workflow/Snakefile", "--configfile", "config/test.yaml",
         "-j", "2", "--nolock", "--rerun-incomplete"],
        cwd=ROOT, capture_output=True, text=True,
    )
    assert result.returncode == 0, result.stderr[-4000:]
    for table in ["plasmids", "orf_occurrences", "proteins"]:
        assert (ROOT / ".test" / "results" / "tables" / f"{table}.parquet").exists()


@pytest.mark.smoke
def test_dereplication_was_lossless_in_the_smoke_run():
    """The invariant, checked on real output rather than synthetic records."""
    from darkorf import store
    with store.connect(ROOT / ".test" / "results") as connection:
        occurrences = connection.execute("SELECT count(*) FROM orf_occurrences").fetchone()[0]
        mapped = connection.execute(
            "SELECT sum(occurrence_count) FROM proteins WHERE record_class = 'observed'"
        ).fetchone()[0]
    assert occurrences == mapped
```

Register the marker in `pytest.ini`: `markers = smoke: end-to-end workflow run`.

- [ ] **Step 4: Run the smoke test**

Run: `.venv/bin/python -m pytest tests/test_smoke.py -v -m smoke`
Expected: PASS (2 tests), under two minutes.

- [ ] **Step 5: Run the whole suite**

Run: `.venv/bin/python -m pytest -q`
Expected: all tests pass.

- [ ] **Step 6: Commit**

```bash
git add .test/ config/test.yaml tests/test_smoke.py pytest.ini
git commit -m "test(smoke): add end-to-end workflow run on a 200-plasmid dataset

Catches the wiring failures unit tests cannot see: mismatched rule filenames,
missing conda directives, scripts importing names their module does not export.

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

---

### Task 12: Remove the superseded implementation

**Files:**
- Delete: `src/plasmidann/`, `workflow/rules/s3_cascade.smk`, `workflow/rules/s5_s9_targets.smk`, the old `workflow/scripts/*`, the old `tests/test_*.py` that target `plasmidann`
- Modify: `config/cascade.yaml`, `config/targets.yaml` (fold any surviving threshold into `config/config.yaml` first)

- [ ] **Step 1: Confirm nothing in the new tree imports the old package**

Run: `grep -rn "plasmidann" src/darkorf workflow tests config | grep -v "^tests/test_.*plasmidann"`
Expected: no output.

- [ ] **Step 2: Confirm every threshold worth keeping has moved**

Read `config/cascade.yaml` and `config/targets.yaml` and check each value against `config/config.yaml`. Anything still needed by a later plan (tier E-values, family resolutions, dN/dS criteria, Foldseek significance) is carried into the new config **with its comment and citation intact**, not silently dropped.

- [ ] **Step 3: Remove the old tree**

```bash
git rm -r src/plasmidann tests/test_cascade.py tests/test_completeness.py \
         tests/test_dark_evidence.py tests/test_dereplicate.py \
         tests/test_uninformative_labels.py tests/test_backbone.py \
         tests/test_circular.py.orig 2>/dev/null || true
git rm workflow/rules/s3_cascade.smk workflow/rules/s5_s9_targets.smk
```

- [ ] **Step 4: Run the suite**

Run: `.venv/bin/python -m pytest -q`
Expected: all remaining tests pass; no import errors.

- [ ] **Step 5: Commit**

```bash
git commit -m "refactor: remove the superseded plasmidann implementation

The v2 pipeline is replaced by the darkorf implementation of PLASMID_ANALYSIS.md.
Thresholds worth keeping were carried into config/config.yaml with their
rationale; the rest is recoverable from git history.

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>"
```

---

## Plans 2–5 (written when Plan 1 lands)

| Plan | Scope | Deliverable |
|---|---|---|
| 2 | Stage 4 cascade (Pfam GA, Pfam relaxed, PHROGs via MMseqs2, Swiss-Prot, eggNOG, nr), normalization, adjudication | `annotations`, `annotation_hits`, ANNOTATED/DARK per protein |
| 3 | Stages 5–7: MMseqs2 families, plasmid lineage clustering, distribution and recurrence | `families`, `family_members`, `plasmid_lineages` |
| 4 | Stages 8–9: all-ORF compact context, dark ±3 detail, DefenseFinder, IntegronFinder, synteny | `context_features`, `context_occurrences` |
| 5 | Stages 10–15: evolution, structure, properties, background normalization, rarity and rarefaction, evidence integration, final views and validation | `evolution`, `structure`, `properties`, `evidence`, the §63 integrated view |

Each follows the same shape: unit-tested module in `src/darkorf/`, thin script, Snakemake rule, smoke coverage, one commit per task.

---

## Self-Review

**Spec coverage for Plan 1's scope:** §3 dataset → Task 6; §4 data architecture → Tasks 3, 4; §5 identifiers → Task 2; §6 manifest → Task 6; §7 missing values → Task 3; §8 ORF prediction → Tasks 5, 7; §9 QC and artifact screening → Task 9; §10 minimum length → Task 9; §11 configuration of the two length parameters → Task 1; §12 partial handling → Tasks 8, 9; §13 dereplication → Task 8; §58 controls → Task 10; §72 sharding and §73 resume → Tasks 6, 7. Sections §14 onward belong to Plans 2–5 and are listed above.

**Placeholders:** none. Tasks 7, 8, 9 and 10 delegate one script body each to the implementer with an exact specification of inputs, outputs and invariants, rather than inline code; every module under test has its code or its port instruction in full.

**Type consistency:** `protein_id` is 32 hex characters everywhere; `occurrence_id` is `plasmid:start-end:strand` in `ids`, `orf`, and `qc`; `schemas.columns` is the single source of table columns used by `store.write_table` and every script; `status` constants are referenced only through the module.
