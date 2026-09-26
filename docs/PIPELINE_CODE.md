# Pipeline code: a walkthrough

What each file does, why it does it that way, and where the design decisions live.

Companion documents:

| document | covers |
|---|---|
| `plans/2026-09-10-pipeline-v2-design.md` | the architecture and the scientific rationale |
| `docs/annotation_statistics.md` | `-Z`, `--cut_ga`, the two Pfam tiers, E-value placement |
| this file | the code itself: layout, contracts, invariants, how to run it |

---

## 1. How to run it

```bash
# Everything the pipeline invokes lives in this environment. All eleven executables.
export PATH=envs/plasmidann/bin:$PATH

snakemake -s workflow/Snakefile -n            # plan only: 1,481 jobs
snakemake -s workflow/Snakefile --lint        # static checks
snakemake -s workflow/Snakefile -j 96         # run

.venv/bin/python -m pytest -q -m "not slow"   # 232 tests, ~9 s
.venv/bin/python -m pytest -q                 # adds the Foldseek integration test (~35 s)
```

On the cluster, submit `workflow/run_pipeline.sbatch`. It puts `envs/plasmidann/bin` on
PATH and runs `preflight` on its own first, so a missing tool or database fails in seconds
rather than after 96 cores have been spent on gene calling. Re-submitting the same script
resumes: every stage is one job, `--rerun-incomplete` discards the partial output of the
job that was killed, and everything finished before it is kept.

**One environment, not one per rule.** `workflow/envs/plasmidann.yaml` declares everything.
Per-rule environments were tried and removed for a specific reason: pre-flight cannot check
an environment it is not running in, so under per-rule environments it would report on its
own environment and say nothing about the one S7b will get. The split also caused the
failure it was meant to prevent — the downstream environment declared `integron_finder`
without `prodigal` or `infernal`, which IntegronFinder shells out to, so the stage failed at
run time with both tools visibly installed on the machine.

The first rule to execute is always `preflight`. It walks
`plasmidann.tools.REQUIRED_TOOLS` — every executable the workflow invokes, including the two
IntegronFinder invokes indirectly — plus every configured database.
`tests/test_tools_registry.py` scans `workflow/scripts/` for subprocess calls and fails if a
script runs anything the registry does not name, so the list cannot fall behind the code.

---

## 2. Layout

```
config/
  config.yaml              paths, ORF length floor, random seed
  cascade.yaml             EVERY search threshold, and the tier list
  schemas/                 JSON schemas both configs are validated against

src/plasmidann/            pure logic - no I/O, no Snakemake, fully unit-tested
  cascade.py               labels, significance, classification, coverage, thresholds
  circular.py              origin repair for circular plasmids
  dereplicate.py           exact-identity collapse
  orfindex.py              stable orf_id assignment
  tools.py                 every external executable and the stage that needs it
  features.py              GFF3 and GenBank conventions, origin-spanning genes included
  orthology.py             eggNOG-mapper output, and its "-" placeholder trap
  labels.py                the tool-derived label vocabulary, and its kinds
  labeldb.py               S4d plasmid label databases: tiers, CARD calls, disagreements
  conjscan.py              S8f CONJScan output and the plasmid mobility class
  context_terms.py         S8c context terms per family, counted over lineages
  synteny.py               Stage 9 conservation, a fractional vote over lineages
  selection.py             S2s which proteins the cascade annotates
  cooccurrence.py          S8g dark families on the same plasmid, per lineage
  normalise.py             collapsing free-text product names
  pfam_meta.py             Pfam description, type and clan from the release
  controls.py              the positive control (SC2)
  context.py               directons, neighbourhoods, island overlap
  evolution.py             Nei-Gojobori dN/dS with status codes
  evidence.py              reality tests and darkness state (description, no ranking)

workflow/
  Snakefile                S0-S2b rules
  rules/
    common.smk             tier chaining helpers
    annotation_cascade.smk         S3-S4 rules
    evidence.smk      S5-S9b rules
  scripts/                 thin I/O wrappers around src/plasmidann
    _ctx.py                shared preamble: import path, and log capture
  envs/plasmidann.yaml     the one environment
  run_pipeline.sbatch      cluster submission: one job runs every stage
  structure_gpu.sbatch     optional: S8d on a GPU (STRUCTURE_ON_GPU=1)
  envs/conjscan.yaml       MacSyFinder 2.1.6 for CONJScan, built as envs/conjscan
  envs/amrfinder.yaml      AMRFinderPlus 4.2.7, built as envs/amrfinder
  bench_nr.sbatch          measures the one cost that is still unmeasured

tools/                     outside the DAG
  download_label_dbs.py    installs the label databases into data/refs/labels/<db>/
  install_tool_envs.py     builds envs/conjscan and envs/amrfinder, installs their data
  calibrate_context.py     post-run: calibrates the context terms on known families

tests/                     one file per concern, 559 tests
  conftest.py              the harness that runs a workflow script against a fixture
docs/, plans/              rationale
```

The split between `src/plasmidann/` and `workflow/scripts/` is deliberate and load-bearing:
**every decision lives in `src/`, where it can be unit-tested without Snakemake, a
cluster, or a database.** The scripts parse files, call tools and write files. When a
script contains a judgement, that judgement has escaped its test coverage.

---

## 3. Data flow and file contracts

```
config: master_table, fasta (one file, may be the whole working set)
   |
   |  S0  analysis_set.py      exclude simulated and lab artifacts, ids AND sequence
   v
results/01_analysis_set/analysis_set.txt                        143,503 plasmid ids
results/01_analysis_set/analysis_set.fna                        the in-scope sequences
   |
   |  S1  orf_call.py          pyrodigal + circular-origin repair, over a process pool
   v
results/02_orf_calling/orfs.tsv            plasmid_id start end strand partial spans_origin seq
   |
   |  S1  orf_index.py         assign orf_id once, over the complete set
   v
results/02_orf_calling/orf_index.tsv       + orf_id                    9,317,050 rows
   |
   |  S2  dereplicate.py       SHA-256 on the exact sequence
   v
results/03_dereplication/unique_proteins.faa                     3,497,616 sequences
results/03_dereplication/protein_map.tsv     seq_id -> orf_id,orf_id,...
   |
   +---> S2z check_hmmer_z.py    -Z confirmed before anything searches with it
   +---> S2b artefact_screen.py  AntiFam + tantan
   |     results/04_orf_qc/artefact_flags.tsv
   |
   |  S3  sweep_cohort.py      2% that bypass narrowing
   |      tier_search.py       once per tier, chained
   v
results/05_annotation_cascade/T{n}/hits.tsv       one row per hit, carrying its own thresholds
results/05_annotation_cascade/T{n}/spans.tsv      cumulative merged informative spans
results/05_annotation_cascade/T{n}/unresolved.faa the next tier's input
   |
   |  S3  cascade_resolve.py   one row per unique protein
   v
results/05_annotation_cascade/protein_annotation.tsv
   |
   |  S4  annotate_plasmids.py join back out to every ORF
   v
results/06_annotation_tables/plasmid_annotation.tsv                  the primary deliverable
```

### Invariants held across the flow

| invariant | where enforced |
|---|---|
| No ORF is lost by dereplication | assertion in `dereplicate.py`; test over 9.29M ids |
| `orf_id` is assigned once and never renumbered | `assign_orf_ids` raises if already indexed |
| Gene calling is independent of where a circle was cut | `tests/test_circular.py`, SC6 |
| A tier never receives a hit for an id it did not query | `narrow` and `narrow_by_explained` raise |
| Every output row carries the thresholds that produced it | `thr_*` columns |
| An empty output is a failure, not a result | assertions in `orf_call.py`, `annotate_plasmids.py` |

---

## 4. `src/plasmidann/cascade.py` - the decisions

Three questions, answered per protein.

### Did anything *name* it?

`UNINFORMATIVE` is a verbose regex matching labels that record an observation rather than
a function: `hypothetical protein`, `uncharacterized`, `DUF1234`, `UPF0102`, `predicted
protein`, `unnamed protein product`, and others.

Recall matters far more than precision here, and the asymmetry is the reason the pattern
is aggressive: **a missed pattern silently promotes an unknown protein to `FUNCTIONAL` and
removes it from the screening set forever.** A false positive merely keeps a named protein
in the pool a little longer, where later evidence demotes it. Recall is measured against a
25-case labelled set and CI fails if it drops.

`is_informative(None)` returns `False`. v1 returned `True`, because
`not UNINFORMATIVE.search(None or "")` evaluates to `True` - so any protein whose label
failed to parse was treated as functionally annotated. Failing toward "we do not know" is
the safe direction for a discovery pipeline.

### Is the hit real enough to count?

`passes_significance(evalue, max_evalue)` runs at parse time, before a hit can contribute a
span, a label or a row.

The quadrant that matters:

| | strong E-value | weak E-value |
|---|---|---|
| **high coverage** | `FUNCTIONAL`, correct | **a spurious long alignment silently removes a genuine dark protein** |
| **low coverage** | `DOMAIN_ONLY`, correct | should never have been recorded |

The top-right cell is the worst error this project can make, because no downstream stage
can recover from it. Hence: **the E-value gates whether a hit exists; coverage decides the
class among the survivors.** It is not a second classification axis - the class drives one
binary decision, and a `FUNCTIONAL_WEAK` class would either behave like `FUNCTIONAL`
(decoration) or like `DOMAIN_ONLY` (where it belongs).

`max_evalue: null` on T1 is deliberate, not an oversight. Pfam gathering thresholds are
per-family bit-score cutoffs set by each family's curator, and for some short families GA
is looser than any global cut. A blanket floor would override curation and degrade the
highest-quality signal in the cascade.

### What class, and how deep?

```python
classify(hits, explained, min_coverage, tier_order)
```

`explained` is the **merged** informative coverage of the whole protein. This is the fix
for the defect that would have put the plasmid backbone on your bench: a replication
initiator carrying `RepA_N` (0.40) and `Bac_RepA_C` (0.52) is 92% explained between them,
but classified per hit both domains fall below `min_coverage` and the protein came out
`DOMAIN_ONLY` - which is target-eligible. Measured on 17% of all Pfam-hit proteins, with
T4SS ATPases and Tn3 transposases among the named examples.

The **label** comes from the strongest E-value, not the widest alignment. A longer
alignment is not a better identification; 15.4% of labels change, and the case that settled
it was `ABC_membrane` at E=1e-23 being chosen over `Peptidase_C39` at E=6.5e-40.

`tier_order` is passed in rather than read from a constant, so a four-tier cascade reports
depths on a four-tier scale. v1 indexed into a hardcoded eight-tier list.

### Coverage, over two disjoint populations

`explained_fraction(length, intervals)` merges overlapping intervals rather than summing
them - a protein hit twice by the same family at 8-337 and 263-387 is 97% explained, not
168%. It is used twice, with different inputs:

```python
explained_fraction(length, informative_spans(hits))    # explained_fraction
explained_fraction(length, uninformative_spans(hits))  # dark_covered_fraction
```

`dark_covered_fraction` is the column that makes the dark set legible. `annot_completeness`
is constant `NONE` there by construction (measured 71/71 and 703/703), because it is built
from informative spans and a dark protein has none. The dark analogue distinguishes:

* 95% covered by "hypothetical protein" across three databases - a real, conserved,
  full-length protein nobody has named. **Strong target.**
* one 20-aa "hypothetical" fragment hit - weak, possibly spurious. **Weak target.**

v1 kept the labels of uninformative hits and discarded their coordinates, so these two were
indistinguishable in the output.

### Thresholds

No function reads a module-level threshold. All arrive as arguments from
`config/cascade.yaml`. `check_thresholds(cfg)` validates the relationships a JSON schema
cannot express:

* `narrow_at >= min_explained` - otherwise the cascade stops searching proteins it then
  reports as unexplained, and because the deeper tiers never ran, the contradiction cannot
  be investigated.
* `full_at >= partial_at` - otherwise the completeness bands overlap and `FULL` is
  unreachable.

`check_hmmer_z(declared, actual)` guards a subtler drift. `hmmer_z` is pinned in config so
that E-values are comparable across tiers and runs, but it is *derived* from the analysis
set - it is the number of unique protein sequences. Anything that changes the ORF set
changes it, and S1 origin repair changes the ORF set substantially by reconstructing
~160,000 genes the linearisation had split. A stale value would silently rescale every
E-value in the run, which is exactly the failure `-Z` exists to prevent. The check runs in
`check_hmmer_z.py`, directly after dereplication; the artefact screen and every tier depend
on it, so nothing searches with an unconfirmed -Z, and it raises with the correct number in
the message. (It used to run in `sweep_cohort.py`, hours later and after the artefact
screen had already searched.)

---

## 5. `src/darkorf/circular.py` - origin repair

A plasmid is a circle; a FASTA record is a line. Cutting the circle breaks any gene
spanning the cut into two fragments, one at each end.

```
    original    [1 .................................... L]
    extended    [1 .................................... L][1 ... overlap]
                                              ^gene now contiguous^
```

Measured on the current analysis set: **160,375 partial ORFs (1.72% of 9,317,050), 1.12
per plasmid, 88,600 of them starting at coordinate <= 3.** 134,748 of 143,503 plasmids
(93.9%) are effectively circular.

These are the *opposite* of what the artefact screen catches. AntiFam finds things that are
not genes; these are real genes that are not whole - and they need the opposite treatment,
repair rather than exclusion. A truncated protein aligns to only part of a domain, so it
drifts toward `DOMAIN_ONLY` or out of annotation entirely and can enter the screening pool
as a novel dark protein. **Synthesising half a protein guarantees a dead well.**

`resolve_origin_genes` maps extended coordinates back, in four cases:

| case | action |
|---|---|
| `end <= L` | ordinary gene, unchanged |
| `start <= 3` and partial | truncated at the record start: the head fragment of a gene called intact across the cut; drop |
| `start <= L < end` | crossed the cut: wrap `end`, set `origin_spanning=True`, clear `partial` |
| `start > L` | wholly in the appended tail - a duplicate of one already kept; drop |

A gene with `origin_spanning=True` runs `start..L` then `1..end`, the GenBank `join()`
convention, so `start > end` and a naive `end - start` is negative. `orf_call.py` writes
the flag to `orf_index.tsv` as the integer column `spans_origin`. Any writer that ignores
this flag has a bug.

**The invariant (SC6):** calling genes on any rotation of a circular sequence must give the
same protein set. If it does not, the coordinate system is contributing biology, which it
must never do. `tests/test_circular.py` asserts this over four offsets.

One measured caution, recorded so nobody "fixes" it: **32.2% of circular plasmids carry a
duplicated protein sequence.** That is genuine multi-copy IS and transposase biology, not
an artefact. DTR records show only 0.9%, confirming their terminal repeats are already
trimmed.

---

## 6. The workflow scripts

### `orf_call.py` (S1)

Calls genes with pyrodigal in meta mode, applying origin repair where topology says the
molecule is closed. Unknown topology is treated as linear - extending a genuinely linear
molecule would fabricate a junction and could invent a chimeric gene across the two ends.

Two v1 defects addressed:

* **The 0-byte FASTA.** v1 declared an output `faa`, opened it, and never wrote to it. All
  600 files were empty, and the DAG reported success because `dereplicate` declared them as
  inputs but only ever read `orf_index.tsv`. The fix **removes the output** rather than
  populating it - protein sequence already travels in the `seq` column, and the FASTA that
  downstream stages consume is `results/03_dereplication/unique_proteins.faa`. Deleting an unused output
  removes the class of failure; populating it would have removed only the symptom.
* **The length floor off by one.** pyrodigal's `min_gene` counts the stop codon, so
  `min_aa * 3` set the real floor one residue below the declared value. Now `(min_aa + 1) * 3`.
  This changes the ORF count by ~0.14%; it is corrected because a declared threshold that
  does not mean what it says cannot be reasoned about.

An assertion fails the rule if no genes were called - an empty output is always a bug.

### `artefact_screen.py` (S2b)

`hmmsearch` against AntiFam plus `tantan` for low complexity, producing `artefact_flag`
with the evidence for it.

AntiFam is Pfam's companion database, curated as a blocklist of the artefact families that
researchers kept independently rediscovering and reporting as novel conserved hypothetical
proteins: shadow ORFs on the reverse-complement strand of real genes, translated rRNA and
tRNA, repeat-derived ORFs.

**These matter disproportionately here.** Our selection criterion is *nothing named it*,
and a shadow ORF is by construction something nothing named - because it is not a protein
and no database contains it. It passes every tier cleanly. Worse, these artefacts are
*conserved*, because the real feature underneath them is conserved, so they would also
survive the multi-lineage and purifying-selection tests planned for S7. They look like
ideal candidates all the way to the plate. Plasmids are high-yield for them: gene-dense,
GC-skewed, saturated with mobile elements.

Since 2026-09-25 an AntiFam flag also skips the cascade: `cascade_selection` reads
`artefact_flags.tsv`, and a protein with an `antifam_family` gets the role
`artefact_antifam`, so it is in no tier's query set, Tier 0 included. A low-complexity flag
skips nothing. The label databases (S4d), DefenseFinder and CONJScan still read every
unique protein. The negative-control decoys are not in `unique_proteins.faa`, so AntiFam
never screens them and they are searched by every tier: the decoy false-positive rate
measures the cascade thresholds alone, not AntiFam and the cascade together.

Nothing is deleted (P5). A flagged protein stays in every table and count; the exclusion
happens at target selection (`exclusion_reason` `artefact`, or `artefact,not_searched`
for an AntiFam-flagged protein) and stays reversible.

### `protein_clustering.py` (S2f)

Every unique protein into families at the three resolutions, **before** the cascade. A
family is a sequence cluster and needs nothing else, so it can come first, and the
selection below is made on it. Every protein is clustered - small and large plasmids,
annotated or not - so no member is lost to an earlier filter. Stage 5 later reads the same
`families_<resolution>_cluster.tsv` files and adds the annotation.

### `cascade_selection.py` (S2s)

What the cascade annotates, and what it actually searches (spec section 13.3). The rule is
`plasmidann.selection.select(clusters, on_small, skipped)`: proteins that are not skipped,
in families (the primary, intermediate clustering) holding a small-plasmid protein that is
not skipped either. Skipped means annotated by Tier 0 or, since 2026-09-25, flagged by
AntiFam (it reads `04_orf_qc/artefact_flags.tsv`, so it depends on `artefact_screen`); an
AntiFam-flagged protein therefore does not open its family for searching. The selected
proteins are clustered again at 90% identity over 80% of BOTH lengths, and only the
representatives go into `cascade_input.faa`. `selection.tsv` gives every unique protein
its role, in this precedence: `artefact_antifam`, `plasmidscope`, `representative`,
`member` or `not_selected`. An AntiFam flag overrides a Tier 0 annotation of the same
protein.

Coverage of both, not of the member only as in the families: the representative's result is
copied to the member, and a member that is a fragment of a longer representative would
receive a domain it does not have.

### `preflight.py` (S3)

Confirms every configured tool and database exists, and that HMM libraries are pressed.
A dependency of every tier.

v1's most expensive failure was a 45-hour job that died because DIAMOND was not on PATH -
the module load line had been written for HMMER and never updated when the tier list
changed. **The same failure occurred twice.** This rule costs about a second.

### `sweep_cohort.py` (S3)

Selects a deterministic random 2% of proteins that bypass narrowing and are searched by
every tier.

The cascade is self-narrowing, which is a compute optimisation that costs the
counterfactual: for a narrowed protein there is no T4 result, so "what would a different
threshold have given?" needs a complete re-run, and each run is weeks. Measured on v1,
`min_explained` 0.3 -> 0.8 moved the deep-tier set by 76%, and the chosen value sat exactly
at the 25th percentile of the observed distribution - the densest possible place to put a
hard cut, and the one place it could not be checked.

For this cohort the counterfactual exists, so the cost of narrowing at `narrow_at` (0.7
since 2026-09-25) can be reported with a confidence interval. ~70,000 sequences, ~2% of
the compute. The cohort is drawn from `cascade_input.faa`, so AntiFam-flagged proteins,
which are not in it, are excluded automatically.

### `tier_search.py` (S3)

Runs one tier and hands on the residue. Four fixes, marked `FIX` in place:

1. **The O(n²).** v1 built `set(ids)` *inside* a dict comprehension, reconstructing it once
   per id: 13.0 s at n=20,000, projecting to **1.6-4.6 days at n=3.49M** - after the search
   had finished, with nothing to show for it. The set is now hoisted.
2. **The significance gate.** Nothing enters until it clears the tier's `max_evalue`.
3. **`-Z` and `--domZ`.** Pinned to 3,497,616. Without this, an E-value means something
   different on every tier, because each tier's input size depends on what the previous
   tier left. See `docs/annotation_statistics.md` §3.
4. **Best hit by significance, not width**, and the recorded statistic is the **i-Evalue**
   (field 13), not the full-sequence E-value (field 7) - the i-Evalue is what governs an
   individual domain.

Uninformative hits now keep their coordinates, which is what makes `dark_covered_fraction`
possible downstream.

The `--domtblout` column indices are documented in place. All were independently verified
correct during review; they are the one surface that was already sound.

### `cascade_resolve.py` (S3)

Joins every tier's hits into one row per unique protein. Calls `check_thresholds` first, so
an incoherent config fails before a table nobody can interpret is produced.

`min_explained` is applied **here**, post hoc, as a reported flag rather than a filter -
which is what makes it sweepable: the search narrowed on `narrow_at` (0.7), so every
protein explained below 0.7 was seen by every tier, and `min_explained` can be swept up to
0.7 without a re-run.

Every unique protein gets a row, and `annot_source` says where it came from: `self`
(searched), `representative` (a 90% member; the row is its representative's, named in
`annot_representative`), `plasmidscope` (Tier 0), `not_searched` (functional class
`NOT_SEARCHED`, neither dark nor annotated) or, since 2026-09-25, `artefact_antifam`
(AntiFam-flagged, functional class `NOT_SEARCHED`; a PlasmidScope annotation of the same
protein is not used).

### `annotate_plasmids.py` (S4)

A join, not a decision. Expands per-protein annotation back over every ORF sharing that
sequence, and attaches artefact flags. Asserts that no ORF lost its mapping.

---

## 7. Tests

559 tests (measured 2026-09-25; 556 without the 3 slow tool integration tests), one file
per concern; the table lists the cascade core. They are written to read as specifications: each name is a
sentence about behaviour, and each docstring says why the behaviour matters, usually with
the measurement that motivated it.

| file | concern |
|---|---|
| `test_cascade.py` | classification, including the RepA multi-domain case |
| `test_significance.py` | the E-value gate, including `max_evalue: null` for T1 |
| `test_dark_coverage.py` | `uninformative_spans`, `n_dark_databases` |
| `test_thresholds.py` | coherence invariants, configurable completeness bands |
| `test_circular.py` | origin repair and the SC6 rotation invariant |
| `test_uninformative_labels.py` | 25 labelled cases; recall regression guard |
| `test_dark_evidence.py` | the evidence ladder; the MULTISPECIES rung is removed (ClusteredNR titles rarely carry it) |
| `test_explained.py`, `test_narrow.py` | coverage merging, tier narrowing |
| `test_dereplicate.py`, `test_orf_index.py` | losslessness, stable ids |

Every change here was written test-first. Where a fix went in ahead of its test - the
`MULTISPECIES` anchor and `is_informative(None)`, both part of a module rewrite - the test
was afterwards verified to discriminate by running it against the v1 pattern, which fails
it.

---

## 8. `src/plasmidann/` - the S5-S9 decisions

Same rule as the cascade: every judgement lives in `src/`, unit-tested without Snakemake, a
cluster or a database. The scripts parse, call tools and write.

### `labels.py` (S4c) - the tool-derived label vocabulary

Every functional label, verbatim from the tool that produced it, tagged with the KIND of
statement it is: `pfam_family`, `pfam_description`, `pfam_clan`, `swissprot_product`,
`pgap_product`, `gene_symbol`, `cog_category`, `cog_id`, `eggnog_pfam`, `go`, `ec`,
`kegg_ko`, `macsy_system`, `macsy_component`, `integron_element`, `integron_type`, and
since 2026-09-25 the kinds of `labeldb.KINDS` (S4d): `card_amr_family`, `amrfinder_gene`,
`tadb_ta`, `bacmet_compound`, `oritdb_role`, `mobileog_category`, `dbapis_family`,
`acrdb_family`. `protein_labels.tsv` carries a `sub_label` column for them: the detail
each source gives (BacMet gene, oriTDB family, CARD model name, AMRFinderPlus element
type/subtype, ...) and, where a database has several evidence classes, the class of the
reference entry that gave the label.

This replaces a hand-written list of 73 Pfam family names that assigned each one a
biological role. The list could not work, and the measurements say why. Pfam-A 38.2 holds
**30,134 families**, of which **67 mention replication** in their description and **42
mention conjugation**; the list named 16 and 15. It named no MobB and no MobD. **Nine of
its 73 names do not exist in Pfam-A at all**, so those entries had never once matched
anything, through a whole version, invisibly. And `pfam2go` covers 4 of its 16 replication
families, 1 of 16 conjugation families and 0 of 11 mobilisation families, so the gap is in
the published mappings too, not only in one person's reading.

No role is assigned in this module, and `tests/test_labels.py` asserts that no label kind
is a role name, so the list cannot grow back.

Organism names and uninformative titles are not admitted as functional labels. Both stay in
`hits.tsv`, so nothing is lost from the record; a `hypothetical protein` category would
simply be the most frequent, and therefore most apparently enriched, feature in the
collection.

### `labeldb.py` (S4d) - plasmid label databases, and where two disagree

When a DIAMOND hit against TADB, BacMet, oriTDB, mobileOG-db, dbAPIS or Anti-CRISPRdb
becomes a label, when a CARD protein homolog hit is Perfect or Strict, how an AMRFinderPlus
table becomes label rows, and which label pairs conflict. Each source is its own label kind
(`KIND`), so a statement from one is never read as a statement from another.

* **Tiers** (`tier`): >= 80% identity and >= 90% coverage, or > 60% and > 70%, the tiers of
  Islam et al. 2026 (PlasAnn). Coverage is the smaller of query and subject coverage: query
  coverage alone lets a fragment carry a full-length label, subject coverage alone lets a
  multidomain protein carry one domain's label. PlasAnn's database and labels are not used.
* **Every target** (`DIAMOND_TIERED_ARGS`): `--max-target-seqs 0`, pre-filtered at the tier-2
  minimum. A target limit ranks by bit score, not tier, and on the test set dropped tier-1
  hits even at 1,000 targets.
* **CARD** (`card_call`, `card_models`): Perfect is 100% identity over the full reference,
  Strict a bit score at or above the model's curated cut-off from `card.json`, as RGI
  tests it. Protein homolog models only.
* **AMRFinderPlus** (`parse_amrfinder`): the method is the tier, element type/subtype the
  sub_label; an element type the module does not know raises rather than being guessed.
* **Header parsers** (`parse_tadb`, `parse_bacmet`, `parse_oritdb`, `parse_mobileog`,
  `parse_dbapis`, `parse_acrdb`) read the FASTA headers `tools/download_label_dbs.py` writes;
  the evidence class goes into sub_label (`evidence=...`). References are keyed by their
  position in the file, because oriTDB entry names repeat across roles and a DIAMOND
  subject id stops at the first space.
* **Disagreements** (`disagreements`): conflict types `tier0_vs_<source>`,
  `card_vs_amrfinder`, `bacmet_vs_amrfinder`, `card_vs_bacmet`, `tadb_vs_defencefinder`,
  `oritdb_vs_conjscan`. The function reads labels and returns new rows; **it changes no
  label**. Gene names match up to case, punctuation and a prefix of at least three
  characters (`same_gene`); the Tier 0 symbols are the KEGG gene symbols of PlasmidScope's
  KOs (`read_ko_symbols`, `tier0_symbols`), the only use of KEGG.
* **Term prefixes** (`term_prefix`): amr, metal, ta, conj_role, mge, antidefence.

### `conjscan.py` (S8f) - CONJScan output and the mobility class

Reads MacSyFinder's `best_solution.tsv` (`read_best_solution`), the installed CONJScan
version from `metadata.yml` (`installed_version`), and assigns each plasmid the class of
Coluzzi et al. 2022 (`plasmid_class`): pCONJ for a T4SS_type model, else pdCONJ for
dCONJ_type, else pMOB for MOB, else pMOBless. An unknown model type raises, so a model
renamed in a later release stops the stage rather than falling through to pMOBless.

### `context_terms.py` (S8c) - what a family's context holds

A term is `<type>:<name>`: amr, metal, ta, conj_role, mge, antidefence from the label
kinds (`label_term`; AMRFinderPlus by element type), defence and conj from the systems
(`system_term`). KEGG is not a source. `orf_term_sources` applies the two neighbour rules:
a gene label counts from a neighbour within +-3 genes in the ORF's directon (same strand,
gaps <= 100 nt, the FESNov rule); a system counts when the ORF is a component or a
component lies within +-3 genes on either strand. `family_term_rows` drops tandem
paralogues (a neighbour in the focal family), counts each term over Stage 6 lineages and
reports TOO_FEW_LINEAGES below two. `window_covers_plasmid` marks an ORF whose window
already holds every other gene of its plasmid.

### `selection.py` (S2s) - which proteins the cascade annotates

`select(clusters, on_small, skipped)` returns the proteins to annotate: every member of a
primary family that is not skipped, when at least one of those unskipped members is on a
small plasmid. `skipped` is the union of the Tier 0 annotated proteins and, since
2026-09-25, the AntiFam-flagged proteins, so an artefact is neither searched nor opens its
family for searching.

### `cooccurrence.py` (S8g) - dark families that travel together

`cooccurrence(plasmid_families, lineage_of, n_lineages, min_lineages_together)` tests every
pair of dark families that share a plasmid in at least `min_lineages_together` (2) Stage 6
lineages. Together means a member ORF of each on the same plasmid; a lineage counts once
however many of its plasmids carry both. The test is the hypergeometric upper tail
(`hypergeom_sf`) over N lineages, in log space with `math.lgamma` because scipy is not a
declared dependency and the tails fall below the smallest double; it sums the tail away from
the mode and stops once a term no longer changes the sum. Because k counts shared plasmids,
the test is conservative relative to lineage-level overlap. `benjamini_hochberg` gives the
step-up q-values across the tested pairs (Benjamini & Hochberg 1995), which is the tested
set only, since the minimum depends on the outcome. A family below the minimum is left out
of the enumeration, which is exact and keeps the pair count bounded. `family_partners`
gives, per family, the number of partners at q <= `fdr` and the best partner (lowest q,
named whether or not significant) with the fraction of this family's lineages in which
the two share a plasmid. Significance means linkage or co-transfer, not function; with
N of the order of 1e5 lineages, K = n = k = 2 gives p = 1/C(N, 2).

### `synteny.py` (Stage 9) - conservation counted over lineages

The six section-42 measurements (`conservation`), each a fractional vote over Stage 6
lineages (`_lineage_modal`): a lineage has weight 1, split equally over its non-empty
values; the conservation is the top score over the voting lineages. Scores are exact
fractions, so ten copies worth 1/10 are one vote. `n_lineages_discordant` counts lineages
whose copies disagree on the left-right arrangement. Fewer lineages than `min_lineages` (2)
give TOO_FEW_LINEAGES and no value. The column names carry the `lineage_` prefix so that
they cannot be mistaken for the occurrence statistic they replace.

### `controls.py` (S5) - the positive control

`control_recall()` raises on an empty control set rather than returning 1.0. An empty
control silently passing is exactly how a gate stops being a gate.

It lived in `backbone.py` and was never part of that list: it draws on reviewed Swiss-Prot
proteins spiked into the query set, which is what lets it detect a protein the cascade
MISSED - something a self-drawn control set can never do.

### `context.py` (S8)

`directons()` groups consecutive same-strand genes separated by at most 100 nt into
putative transcriptional units. On a circular plasmid the record's last and first units
merge when the strand matches and the gap across the origin is within the same limit; the
length it needs comes from S0's `01_analysis_set/plasmid_lengths.tsv` (after terminal-repeat
removal, which `size_bp` does not reflect). Membership is a far stronger claim than adjacency: the
genes are predicted to be co-transcribed, so a dark ORF inside an otherwise annotated
operon inherits that operon's hypothesis.

`neighbourhood()` truncates at the ends of a record rather than wrapping. Wrapping is
correct for a circular molecule but must be the caller's explicit decision, not an accident
of negative list indexing - which is how it would happen silently in Python.

**On a 5 kb cryptic plasmid carrying six genes, a ±3 neighbourhood IS the whole plasmid.**
S8c reports context as descriptive per-family rates (`cons_*`) with no background
correction - the enrichment test, its stratified background and the label-category layer
were removed on 2026-09-25 - so a high neighbour rate on small plasmids is expected by
construction. The context terms (`context_terms.py`) report the same situation per
family as `window_covers_plasmid_fraction`.

### `evolution.py` (S7)

Nei-Gojobori dN/dS by counting, and since 2026-09-14 the pipeline's only dN/dS estimate.
It needs no tree, no optimiser and no external process, so it runs over hundreds of
thousands of small families in-process and is fully testable. The codon model that used to
sit beside it is gone; see S7d below.

This is the only filter that catches the artefact class nothing else does: a shadow ORF on
the reverse-complement strand of a real gene is conserved, multi-species, and passes every
absence-based test - but the selection acting on that locus is acting on the gene on the
*other* strand.

Three deliberate return-value decisions:

* **Identical sequences return `None`, not `0.0`.** Reporting them as dN/dS = 0 would look
  like maximal purifying selection and would promote every clonal duplicate to the top of
  the list.
* **Gapped and ambiguous codons are skipped.** Counting a gapped column as identical would
  inflate apparent conservation precisely where the alignment is least trustworthy.
* **Saturation returns `None`.** Beyond pS = 0.75 the Jukes-Cantor correction has no
  domain; returning a number there would be inventing one.

`back_translate()` projects a protein alignment onto codons - every protein gap becomes
exactly three nucleotide gaps, never one or two, or the reading frame is destroyed.

### `evidence.py` (S9) - reality lines and darkness state

What used to be `peptide.py` and `targets.py` - the peptide heuristics, eligibility,
strata, the lexicographic ranking and the portfolio - was removed with Layer C on
2026-09-17 (spec section 76). What remains describes a family and selects nothing.

`reality_lines()` counts and names the four independent lines of evidence that a family is
a real protein (`purifying_selection`, `multi_lineage` over Stage 6 lineages, `is_family`,
`folds`); `IMPLIED_BY` records that `purifying_selection` entails `is_family`, so the two
count once. Every threshold arrives from config through `reality_thresholds()`.
`darkness_state()` returns one of two named states, DARK_FOLD_KNOWN or DARK_NO_STRUCTURE.

---

## 9. The S5-S9 scripts

| script | stage | what it does |
|---|---|---|
| `clonal_registry.py` | S0b | MOB cluster per plasmid, so "independent occurrence" means something |
| `quality_gate.py` | S5 | positive control (**halts the run**) and artefact flags |
| `dark_set.py` | S6a | eligible ∩ not-partial-only |
| `cluster_dark.py` | S6b | MMseqs2 at 30%/50%; singletons kept as `ORPHAN` |
| `extract_cds.py` | S7a | nucleotide recovery, honouring `spans_origin` and strand |
| `family_evolution.py` | S7b | mafft → codon projection → dN/dS per family |
| `label_databases.py` | S4d | every unique protein against TADB, BacMet, oriTDB, CARD, mobileOG-db, dbAPIS and Anti-CRISPRdb (DIAMOND) and AMRFinderPlus (`--plus`, from `envs/amrfinder`); writes `08_protein_labels/protein_labels_plasmid.tsv` and `label_databases_status.tsv` (SUCCESS, NO_HIT or NOT_RUN per database, so an absent database never reads as "nothing found") |
| `protein_labels.py` | S4c | merges every source into `protein_labels.tsv`, the S4d rows as their own kinds with `sub_label`; writes `label_disagreements.tsv` (`labeldb.disagreements`), which changes no label and is kept for later analysis |
| `defence_systems.py` | S8a | DefenseFinder, as ONE MacSyFinder database with `--worker` set to the threads: per-core chunks made the calls depend on `-c`, because HMMER's i-evalue scales with the database size |
| `conjugation_systems.py` | S8f | CONJScan 2.1.0 (Plasmids models) over every ORF of every plasmid as one MacSyFinder database, run by the MacSyFinder 2.1.6 of `envs/conjscan`; writes `conjugation_systems.tsv` and `conjugation_plasmid_class.tsv` (pCONJ, pdCONJ, pMOB, pMOBless); NOT_RUN writes header-only tables |
| `integrons.py` | S8b | IntegronFinder, `--local-max` for CALIN elements |
| `context_features.py` | S8c | directons, islands, neighbours; one row of `cons_*` rates per family (`cons_conj` included), and `family_context_terms.tsv`: the context terms per family counted over lineages, for the dark families and, as the calibration benchmark, every known family at the primary resolution |
| `dark_cooccurrence.py` | S8g | reads `dark_families.tsv` (primary-resolution dark families, through their dark members), `protein_map.tsv` and `plasmid_lineage.tsv`; writes `12_context_and_structure/dark_cooccurrence.tsv`, one row per TESTED pair sorted by p (family_a, family_b, n_lineages_a, n_lineages_b, n_lineages_together, n_lineages_total, fraction_of_a, fraction_of_b, expected_together, p_value, q_value, status); settings `targets.yaml` `cooccurrence` (`min_lineages_together` 2, `fdr` 0.05). On the 100-plasmid test set: 198 dark families on 55 plasmids, 93 lineages, 6 pairs tested, all at q <= 0.05 |
| `synteny.py` | S9a | Stage 9 at every level of `synteny.levels` (close, intermediate): one row per cluster holding a dark small-plasmid member, measured over lineages; the primary-level set must equal `dark_families.tsv` |
| `structure_search.py` | S8d | Foldseek with ProstT5 |
| `annotation_report.py` | final | the two CSV deliverables (section 14). Since 2026-09-25 the ORF table adds the S4d labels by term type (`amr_labels`, `metal_labels`, `ta_labels`, `conj_role_labels`, `mge_labels`, `antidefence_labels`), `conj_system`, `conj_component`, `plasmid_conjscan_class` and the close-level synteny of the ORF's close cluster (`close_*`); the family table carries the intermediate-level `lineage_*` synteny columns and `cons_conj`, and, after the synteny block, the S8g co-occurrence columns (`cooccurrence_status` SUCCESS or TOO_FEW_LINEAGES, `n_cooccurring_partners`, `top_cooccurring_partner`, `top_cooccurring_partner_q`, `top_cooccurring_partner_fraction`, `cooccurrence_fdr`, `cooccurrence_min_lineages`); a co-occurrence status counts towards the GENOMIC_CONTEXT dimension |

Two design notes worth carrying:

**`integrons.py` uses `--local-max`.** That is the sensitive mode, which finds attC sites
beyond those adjacent to a detected integrase - CALIN elements, cassette arrays whose
integrase has been lost. Those are common on plasmids and their cassettes are exactly the
uncharacterised genes we are hunting. A dark ORF in a cassette array is a real gene *by
construction*: it carries an attC site and has been physically excised, mobilised,
re-integrated and retained.

**`structure_search.py` degrades rather than fails.** If Foldseek or the ProstT5 model is
unavailable it writes an empty, well-formed table and S9 scores structure as absent. One
evidence group out of four should not take down a run that has already spent days on the
cascade.

### Tools outside the DAG

| tool | when | what it does |
|---|---|---|
| `tools/download_label_dbs.py` | before the first run | installs TADB, BacMet, oriTDB, CARD, mobileOG-db, dbAPIS, Anti-CRISPRdb and the KEGG KO list into `data/refs/labels/<db>/`, each with `raw/`, the built FASTA, `VERSION`, `SOURCE` (URL, citation, licence, evidence rule, download date) and `MANIFEST.sha256`. Pinned downloads are checked by SHA-256 (CARD and KEGG are taken as the current release); `--local DIR` uses a local copy when its checksum matches. Idempotent |
| `tools/install_tool_envs.py` | before the first run | builds `envs/conjscan` and `envs/amrfinder` from `workflow/envs/conjscan.yaml` and `amrfinder.yaml`, installs the CONJScan 2.1.0 models into `data/refs/conjscan` and the AMRFinderPlus database into `data/refs/amrfinder` (`latest` link). Rebuilds only what is missing or at another version |
| `tools/calibrate_context.py` | after the run | reads `family_context_terms.tsv`, `protein_labels.tsv`, `protein_families.tsv`, the protein map and the two system tables; for every term measures precision on the known families (benchmark: more than half of the members carry the term; negatives: labelled families none of whose members carries it) and writes the conservation at which precision reaches 50% and 90% (FESNov's two confidence levels), a precision curve and a summary. `MIN_FAMILIES = 10`: 90% precision cannot be observed on fewer than 10 families (one false positive in 9 is 89%), so a term with fewer benchmark or negative families is UNCALIBRATED and gets no threshold |

---

## 10. Running it

```bash
export PATH=envs/plasmidann/bin:$PATH

snakemake -s workflow/Snakefile -n -c 2             # plan: 47 jobs (measured 2026-09-25)
snakemake -s workflow/Snakefile annotate_only -j 96 # stop after the annotated plasmidome
snakemake -s workflow/Snakefile -j 96               # the whole project

sbatch workflow/run_pipeline.sbatch                 # on the cluster
```

`annotate_only` exists because S4 is a deliverable in its own right, and because the
quality gate at S5 may legitimately halt a run that produced a perfectly good annotation.

**One submission runs everything** (since 2026-09-25): the cascade, S4d, DefenseFinder,
CONJScan, IntegronFinder, ISEScan, synteny, the context terms, the structure search - on
CPU, inside the same job - and the report. `STRUCTURE_ON_GPU=1 sbatch
workflow/run_pipeline.sbatch` restores the GPU split: the job stops before
`structure_search`, `workflow/structure_gpu.sbatch` runs it, and a second submission
finishes. AMRFinderPlus and the MacSyFinder that runs CONJScan are named by path in config
(`amrfinder.executable`, `conjugation.exe`) and checked there by pre-flight, including
that the MacSyFinder can read CONJScan's model grammar.

---

## 11. One job per stage, and why

There are no shards. Every stage is a single job over one file, and a stage that can use
more than one core - the cascade tiers, gene calling, IntegronFinder - is given every core
the run has (`threads: workflow.cores`).

The pipeline was sharded twice over until 2026-09-22: 600 plasmid shards for gene calling
and IntegronFinder, and 64 protein shards for the cascade. The protein shards were removed
because of one measurement. A search against a streamed database has a fixed cost per
invocation - DIAMOND reads and indexes the whole of nr each time it runs - and
`workflow/bench_nr.sbatch` found that 2,000 queries at 16 threads did not finish in 12
hours. Sixty-four shards were sixty-four such passes; one job amortises the pass over every
query at once. The measured hmmer side never needed sharding: T1 is 65.8 s fixed + 0.0365
s/protein at 4 threads, and hmmsearch takes `--cpu`.

The plasmid shards went with them so that the pipeline takes ONE input file. The analysis
scope is then no longer "whatever files were handed in": S0 filters the configured FASTA
by the master table's locked exclusion and writes `analysis_set.fna`, and every stage that
needs sequence reads that. The configured FASTA may therefore be the whole working set,
simulated plasmids included.

What was given up is resumability granularity. A failure loses the stage, not one shard of
it; `--rerun-incomplete` still keeps every stage that finished. The stage that pays most is
IntegronFinder, which walks replicons one at a time and threads only its HMM searches.

`tier_query` and `tier_spans` in `common.smk` resolve to the previous tier's output, so
explained spans accumulate tier by tier and `cascade_resolve` reads every tier's hit file
and the last tier's span file.

---

## 12. Logging

Every rule declares `log:`. For `script:` rules Snakemake creates that path and does **not**
redirect the script's stdout into it — that only happens for `shell:` and `run:`. So all the
log files stayed empty while the diagnostics the scripts print went to one SLURM output
file, interleaved across up to 1,481 concurrent jobs and impossible to attribute to a rule.

`workflow/scripts/_ctx.py` fixes this once for every script. It finds the injected
`snakemake` object by walking back up the frame stack from its own import, then tees stdout
and stderr into the declared log. Teeing rather than moving: the log gets the rule's own
record and the SLURM file keeps a live progress trace. A script that declares no log is
unaffected.

---

## 13. The stages added after the first review round

Each was named in the design, declared in config, and implemented by nothing.

| stage | rule | what it produces |
|---|---|---|
| **S4** feature files | `feature_files` | GFF3 and GenBank beside the TSV. The work is the 160,375 origin-spanning ORFs: GFF3 forbids `start > end` and needs a discontinuous feature sharing one ID, GenBank wants `complement(join(a,b))` with the complement outside the join. `plasmidann.features` owns both. |
| **S4b** orthology | `orthology` | COG and KEGG terms for the proteins the cascade NAMED. Not a dark-hunting tier: S8 asks what a dark ORF's neighbours do, free text cannot be aggregated into pathways, and FESNov's neighbourhood metric is defined over KEGG membership. |
| **S7b** coding potential | inside `family_evolution` (parallel since 2026-09-24; the per-family work is `plasmidann.evolution_worker`, and every measurement is also reported over the small-plasmid members, prefix `small_`) | RNAcode over the codon alignment the rule already built. Both strands: for a shadow ORF the ANTISENSE signal should be the stronger one, and that is the artefact class nothing else here catches. |
| **S7c** consensus re-check | `consensus_recheck` | The family consensus searched back against Pfam at curated GA. A family can be collectively recognisable while every member individually misses the cut; Pavlopoulos removed 6.5% of clusters this way. |
| ~~**S7d** codon model~~ | ~~`busted_confirm`~~ | **Removed 2026-09-14.** HyPhy BUSTED fitted to a per-family FastTree tree. See below. |

**Why S7d was removed.** BUSTED needs a tree, and FastTree was the only thing in the pipeline that
built one. Measured per family before removal: 7.2 s at 4 sequences, 66.1 s at 10, 141.7 s
at 20, 853.2 s at 50 - against 0.40 s to 12.94 s for everything S7b does to the same family.
The cost is set by the number of sequences and almost not at all by alignment length
(10 sequences of 300 codons: 69.3 s, against 66.1 s for 120). On the measured family-size
distribution that is 204-681 wall-hours against a declared `runtime` of 48 hours.

Two levers were weighed. Lowering `max_members_aligned` from 50 to 10 saves about half,
because the cost sits in the many families of 4-19 members rather than the few large ones,
and still leaves roughly eight days. Parallelising - the rule already declared 16 threads
and used one - would have brought a corrected stage to roughly 13-43 wall-hours. Removal was chosen
with both on the table. Two defects would have had to be fixed under either lever: the rule
fed FastTree *unaligned* CDS, which exits 1 with `Wrong number of characters`, and the
extracted CDS carry their terminal stop codon, which HyPhy refuses outright.

**Absence is a status, never a blank** still holds across the remaining stages:
`TOO_FEW_MEMBERS`, `NO_DIVERGENCE`, `SATURATED`, `NO_OUTPUT` are distinct from a value,
because a blank reads as a failed test rather than a test that never ran.

RNAcode is worth one specific warning. Handed FASTA it prints `ERROR: Unknown alignment file
format` on stdout and **exits 0**, so its return code cannot be trusted; the wrapper writes
Clustal W and checks what came back rather than what the process said.

---

## 14. The deliverable

`rule all` produces **complete annotation, not a shortlist**:

```
results/15_report/annotation_complete.csv       one row per ORF, 9.3M rows
results/15_report/dark_families_complete.csv    one row per dark family
results/08_protein_labels/protein_labels.tsv    every label, per protein
results/08_protein_labels/label_disagreements.tsv   cross-source conflicts, for review
results/12_context_and_structure/family_context_terms.tsv   context terms per family
results/12_context_and_structure/conjugation_plasmid_class.tsv   mobility class
results/12_context_and_structure/dark_cooccurrence.tsv   dark families on one plasmid
results/06_annotation_tables/plasmid_annotation.gff3 / .gbk   the feature files
results/09_quality_gate/quality_gate.txt                 the run-halting control
```

Nothing in those tables is filtered or ranked. Every ORF appears, artefact-flagged ones
included; every dark family appears, ORPHANs and evidence-free families included.

One caveat belongs here rather than only in the source. `under_purifying_selection` and
`dnds_median` come from a single estimator, Nei-Gojobori counting, and since 2026-09-14
nothing in the pipeline can contradict it; the codon model that could was removed (section
13). Read those two columns together with `dnds_status` and `n_pairs`, which carry the whole
of the qualification a second method would otherwise have supplied.

The family
table is the selection surface - `reality_n`, the named evidence lines, dN/dS and its status,
RNAcode on both strands, the consensus re-check, the context hypothesis
and its conservation, the structural match and its description - all side by side, so that
choosing candidates is a decision made ON the table rather than one baked into a rule.

CSV rather than TSV because these are the files that get opened in a spreadsheet, and
`csv.writer` quotes the free-text fields that a DIAMOND `stitle` routinely fills with commas.

The stratified 1,000-construct portfolio and its synthesis order were removed with Layer C
on 2026-09-17 (spec section 76); the 1,000 proteins are chosen by hand from these tables.

---

## 15. What is still missing

| gap | consequence | fix |
|---|---|---|
| **ColabFold/ESMFold** not installed | no pLDDT; `structure.min_plddt` was removed from config because nothing could honour it (`docs/PARAMETER_PROVENANCE.md`) | only needed to confirm shortlisted novel folds; Foldseek + ProstT5 does the screening pass without it |
| **nr tier not benchmarked** | the walltime of the nr tier is unknown and it dominates the run. The one measurement so far: 2,000 queries at 16 threads with default `-b 2 -c 4` did not finish in 12 hours, which is the fixed cost of one pass over the database | `workflow/bench_nr.sbatch`, with a 48 h limit, `-c 1` and a large `-b`, on the 96-core node |
| **no independent check on dN/dS** | Nei-Gojobori counting is the only selection estimate; the codon model that could contradict it was removed on 2026-09-14 (section 13). `purifying_selection` is one of four reality lines and the cheapest to fire | a confirmatory codon model on a chosen shortlist, outside the DAG, if one is ever wanted; the two input defects recorded in section 13 must be fixed first |
| **8 parameters with no source** (count of 2026-09-25) | researcher degrees of freedom, listed in `docs/PARAMETER_PROVENANCE.md` | each needs a citation or an explicit methods defence with measured sensitivity |

None of these block a run. Each degrades one evidence group or leaves one number undefended,
and each is recorded as absent rather than silently assumed. The Foldseek target database and
the ProstT5 model, previously on this list, are downloaded and exercised by
`tests/test_scripts_smoke.py::test_structure_search_reports_what_the_match_actually_is`.
