# Review summary — `plans/2026-09-10-pipeline-v2-design.md`

Consolidates three independent reviews of the Dark ORF Discovery Pipeline specification
(2,304 lines), commissioned 2026-09-16.

| Report | Lens | Findings |
|---|---|---|
| [`scientific_review.md`](scientific_review.md) | Scientific validity, thresholds, feasibility at scale | 57 across 5 categories |
| [`engineering_review.md`](engineering_review.md) | Data model, DAG, IDs, resumability, spec-to-code gap | ~50 across 7 categories |
| [`devils_advocate_review.md`](devils_advocate_review.md) | What fails after the pipeline runs correctly | 13 findings + 7 answered questions |

Each reviewer read the spec, `docs/annotation_statistics.md`, `docs/PARAMETER_PROVENANCE.md`,
`docs/PIPELINE_CODE.md`, `config/*.yaml` and the existing implementation in `src/plasmidann/`,
`workflow/` and `tests/`. Scale figures throughout: 143,504 plasmids, ~9,317,050 ORF occurrences,
~3,497,616 unique proteins, one 96-core node, 600 GB RAM, 7-day wall limit = **16,128 core-hours
per job**.

---

## Verdict

The specification is sound in its *principles* — evidence preservation, no premature collapse into
a score, prioritisation deferred downstream — and unimplementable in its *current form*. Three
root causes account for nearly every finding:

1. **It costs 2–11× the available compute.** §10 runs seven evidence sources over all 3.5M
   proteins where the existing design searches deeply only what remains unexplained.
2. **Its central output has no definition and no test that could falsify it.** "Dark" is defined
   as "according to the configured adjudication rules", and those rules are never given; §32/§39
   validate column presence and tautologies.
3. **It silently discards measured decisions this repository already made.** Every parameter in
   `docs/annotation_statistics.md` — the document the spec names as its companion — is absent.

None of the three is a reason to abandon the spec. All three are reasons not to implement it
verbatim.

---

## 1. Compute: the constraint that reshapes the design

Budget is 16,128 core-hours. §10 as written costs **36,000–180,000 core-hours (15–78 days of the
full node)**. Adding §19's all-dark structural stage puts the floor near 90,000.

| Stage | Core-hours | Wall (96 cores) | Verdict |
|---|---:|---:|---|
| T1 Pfam GA, `hmmsearch` (measured) | ~142 | 1.5 h | fine |
| T1 as §10.1 specifies, `hmmscan` | 1,400–2,800 | 15–30 h | 10–20× waste, and breaks `-Z` |
| T2 Pfam relaxed, full set | 140–2,800 | 1.5–30 h | fine |
| T3 InterProScan (not installed) | 3,000–20,000 | 1.3–8.7 d | must be subset |
| T4 eggNOG-mapper (DBs not installed) | 1,500–4,000 | 16–42 h | feasible if installed |
| T5 DIAMOND Swiss-Prot | 20–60 | <1 h | fine |
| **T6 HH-suite (not installed)** | **20,000–100,000** | **9–43 d** | **impossible** |
| **T7 DIAMOND nr, all 3.5M** | **10,000–50,000** | **4–22 d** | **at/over the limit** |
| §19 ProstT5, dark representatives (~600k) | ~12,000 | ~5 d | at the limit, CPU-only |
| **§19 ProstT5, all dark (~2.6M) as specified** | **~52,000** | **~22 d** | **impossible** |
| Everything downstream (MMseqs2, DefenseFinder, IntegronFinder, RNAcode, dN/dS, Foldseek, properties) | 2,500–9,000 | — | fits comfortably |

Three specific conclusions:

- **HH-suite (§10.6) must be struck or made shortlist-only.** It needs a 100–270 GB database that
  isn't on disk, and §10.6 and §23.2 *both* declare its output non-independent of Pfam and DIAMOND.
  The most expensive tool in the hub buys evidence the spec forbids counting. At 10k shortlisted
  proteins it costs ~60–170 core-hours.
- **InterProScan (§10.3) has the same problem.** Novel plasmid ORFs miss the precalculated UniParc
  lookup by construction, so all 3.5M take the full local calculation — across Gene3D, SUPERFAMILY,
  PANTHER and CDD, which §10.3 itself says are non-independent of Pfam (already T1/T2).
- **§19's "all dark proteins, not family representatives" costs 4.3× for nothing.** Members within
  a 30%-identity cluster are exactly the set ProstT5 maps to near-identical 3Di strings. The
  cheaper fix than restricting the stage is a GPU: the whole dark set becomes ~26 hours.

Unmeasured and blocking: DIAMOND streams the entire 350 GB nr per invocation, so 64 shards pay
that cost 64 times — 6.2 h of pure I/O at 1 GB/s, ~31 h at a realistic 200 MB/s shared filesystem,
before a single alignment. `workflow/bench_nr.sbatch` exists to measure this and has never been
run. No schedule for §10 is credible until it has.

---

## 2. The outcome variable is undefined, and nothing can falsify it

- **§13.2/§13.3 are circular.** ANNOTATED and DARK are both "according to the configured
  adjudication rules"; the spec never states them. The previous design had a concrete, measured
  one (`min_explained: 0.5`, defended in `docs/annotation_statistics.md`). As written, no "X% of
  the plasmidome is dark" claim is supportable.
- **§32 and §39 cannot fail.** They check column presence and tautological inequalities
  (`0 <= percentage_dark_in_family <= 100`). No external referent, no recall control, no negative
  set. The headline count is uncontradictable by construction.
- **The run-halting positive-control gate is gone.** The repo has one: 500 spiked reviewed
  Swiss-Prot proteins that must come back FUNCTIONAL, halting below 0.99 recall
  (`prepare_control.py`, `quality_gate.py`). It is the only instrument that can detect an
  annotation-recall failure, and building §29's agents from the spec alone deletes it.
- **§13.5 `dark_evidence_level` is circular in the DAG as well as the science.** Levels 3–5 need
  family (§15), context (§16) and recurrence (§14), all of which §37 places downstream of the
  classification that assigns the level. It also imposes an order the data do not have, and is
  used for nothing.
- **Context labels have no null.** `defence_associated = TRUE` cannot be false as specified; no
  permutation background is computed.

---

## 3. Evidence independence is overstated

- **§23.3 counts bins, not observations.** `distribution`, `evolutionary_conservation`,
  `genomic_context` and `orf_qc` are all monotone in family size — one big family scores high on
  four "independent" categories for one reason. §23.2 blocks only within-homology double-counting.
  Structural evidence is not independent either: ProstT5 is a sequence model. Reviewers put the
  true count at ~3 separable sources, not 9. The repo's `targets.IMPLIED_BY` entailment map, which
  handled exactly this, is absent from the spec.
- **The independence unit is wrong.** §14, §17, §18.1 and §22 all rest on `mob_cluster`, which is
  not a lineage: cluster AA379 alone is 21,817 plasmids — 15.2% of the dataset and 24.4% of
  everything under 10 kb. A family spread over 400 unrelated cryptic plasmids scores as one
  lineage. §32.6's `MOB_count <= plasmid_count` check cannot detect this. A plasmid ANI/mash
  clustering stage is needed as the independence unit.
- **§14's recurrence measures submission sociology, not biology.** Deposit counts are confounded
  by PGAP/MULTISPECIES propagation, shared upstream records across PlasmidScope/PLSDB/IMG-PR, and
  clonal outbreak over-sampling. No BioProject or ANI dereplication is specified anywhere.
- **§14 is also uncomputable as specified**: `max_target_seqs: 5` caps
  `hypothetical_deposit_count` at 5, and "database record" / "independent source" are never defined.

---

## 4. Artefacts: an estimated 30–50% of the "dark" set

§7.1 freezes the artefact screen at AntiFam's 278 profiles and explicitly forbids additions.
Missed classes, all of them *recurrent by construction* and therefore actively selected for by
§14's recurrence evidence:

- shadow/antisense ORFs on the reverse complement of conserved genes
- **160,375 measured origin-spanning ORFs (1.72%)** that §6 mishandles in both directions
- host chromosomal contamination
- IS-derived ORFs, frameshift-split genes, phage genes

The only artefact evidence the spec could obtain but doesn't retain: nucleotide sequence,
opposite-strand overlap fraction, and RNAcode run on both strands. The existing implementation
already does the last of these.

---

## 5. Engineering defects

- **§25's deliverable is 9–14 GB** (9.3M rows × ~120 columns at ~1.0–1.5 kB/row), needs 50–110 GB
  RAM to load in pandas, and exceeds every spreadsheet row limit. ~90 of those columns are constant
  within `protein_id` or `family_id` — the duplication *is* the file size. Parquet or DuckDB, with
  the flat table reproduced as a view, is the answer.
- **`family_id` is not stable across reruns**, violating §4.2. `cluster_dark.py:86` numbers
  families `F{i:07d}` from `enumerate(sorted(...))`, so adding one dark protein renumbers every
  later family. IDs must be content-derived (hash of representative or member set).
- **nr does not fit a 7-day wall**: 64 shards × 48 h declared at 24 concurrent jobs = 3 waves ≈ 6
  days for that tier alone.
- **Four expensive stages have no shard wildcard** — eggNOG-mapper, Foldseek/ProstT5, family
  evolution, MMseqs2. `family_evolution` also declares `threads: 16` while running `mafft
  --thread 1` serially: the defect class that already got a stage removed once.
- **§5.1/Rule 4 are violated today and unbudgeted when fixed.** `tier_search.py` discards raw
  hmmsearch/DIAMOND output to a temp dir; persisting it is 15–30 GB before Rule 3's rejected-hit
  rows.
- **§9's stage number contradicts its input**: the normalization layer is Stage 4, before the
  Stage 5 hub whose output strings it normalizes.
- **25+ config keys are undefined or literally written `configurable`**, including all shard and
  thread counts; §31's `min_dark_length_aa: 30` contradicts §7.3/§7.4's 20 aa and `config.yaml`'s
  `min_orf_aa: 20`; §2.6 and §13.4 are missing entirely, orphaning the `DARK_HYPOTHETICAL` /
  `DARK_DUF` sublabels. §30 Rule 1 ("do not invent thresholds") is therefore unsatisfiable.
- **§10.1 specifies `hmmscan`, which inverts `-Z` semantics.** `hmmer_z: 3497616` is a sequence
  count; under `hmmscan` it denotes profiles, inflating every Pfam E-value by ~2 orders of
  magnitude.
- **§7.2/§31 specify AntiFam with no threshold**, reopening a defect a measurement and a regression
  test already closed (274 of 278 curated GA cuts are looser than 1e-5).

---

## 6. What already exists

Of roughly 30 stages in the spec: **7 exist and match**, **~16 exist but disagree** in scope or
columns, **7 are genuinely new**.

| | Stages |
|---|---|
| **Have** | dereplication (losslessness asserted + tested), Pfam GA + relaxed, Swiss-Prot, nr, DefenseFinder, IntegronFinder, RNAcode + dN/dS |
| **Differ** | plasmid table, ORF prediction (7 cols vs 11), AntiFam/QC (<20 aa is a caller floor, not a flag), eggNOG (named proteins only, unsharded), thresholds (two-way not three-way), adjudication (no `primary_annotation*`, no hierarchy, no transitive guard), dark classification (no 0–5 level), families (dark set only, one resolution, 8 of 19 columns), context (±3 only, no per-ORF neighbour table), synteny (none of the six §17.2 measures), distribution (2 of 7 counts), structure (no model id/confidence), properties (charge/GRAVY/TM only), normalization (feature background only), evidence integration (4 reality lines vs 9 categories), final tables (CSV, ~45 cols) |
| **New** | §9 normalization layer, §10.3 InterPro, §10.6 HMM-HMM, §14 recurrence, §18.2 conservation, §22 rarity labels, §27 plasmid summary |
| **In repo, absent from spec — keep** | preflight, circular-origin repair, spiked control gate, feature files, consensus re-check, clonal registry, 2% sweep cohort, opt-in S9 portfolio |

On §28's proposed `dark_orf_pipeline/` tree, the engineering reviewer's verdict: *"rename is pure
churn; only `raw/intermediate/normalized/final` is worth adopting."*

---

## 7. Two decisions that are cheap now and expensive later

Both change what the *first* run must record, so they cannot be deferred:

1. **§6.3 discards nucleotide sequence.** RNAcode (§18.3) and dN/dS (§18.4) both need codon-level
   data. Without it, those stages are one-shot: re-running them means re-extracting CDS from all
   143,504 plasmids. Keeping nucleotides costs disk and nothing else.
2. **§16.1 extracts context for dark ORFs only.** Every context enrichment claim then lacks a
   background — "dark ORFs are defence-associated" is meaningless without the rate for non-dark
   ORFs. Fixing it later is a re-run of all of Stage 8 over ~9.3M ORFs.

---

## 8. The ~1,000 figure is a budget number

Nothing in the spec establishes that the dataset supports ~1,000 distinct, experimentally tractable
dark families. This project's own earlier measurement gives ~5,028 dark-only families with ≥10
members and 15,262 cross-lineage. The rarefaction curve of dark families against plasmids sampled
would settle the question, is nearly free to compute, and is required by neither §26, §38 nor §39.

---

## 9. Recommended actions before implementation

In the order the reviewers put them:

1. **Restore the run-halting control gate** and add a held-out arm. Without it §39 is not a
   definition of success.
2. **Move §10.6 (HH-suite) and §10.3 (InterPro) behind adjudication, and restore the
   narrowing/reporting split.** This is where the budget for everything else comes from.
3. **Make three columns mandatory in §6.3/§25**: nucleotide sequence, opposite-strand overlap
   fraction, RNAcode on both strands. These are the only artefact evidence the pipeline can obtain.
4. **Extract neighbourhoods for all ~9.3M ORFs**, not only dark ones, and compute a permutation
   null during the run.
5. **Add plasmid ANI/mash clustering** and use it — not `mob_cluster` — as the independence unit in
   §14, §17 and §18.1.
6. **Require the rarefaction curve** as a §38 deliverable.
7. **Delete `dark_evidence_level` (§13.5)** or replace it with the counts it is a function of; drop
   `independent_evidence_categories` (§23.3) or ship the joint-indicator table beside it.
8. **Define "dark"** with a concrete, measured threshold, and state it in config where it can be
   validated, swept and stamped into output rows.
9. **Run `workflow/bench_nr.sbatch`** before committing to any §10 schedule.

Items 1, 2 and 3 already exist, working and tested, in `src/plasmidann/` and `workflow/scripts/`.

---

## 10. Open decisions requiring the project lead

| # | Decision | Options |
|---|---|---|
| D1 | Build strategy | Restructure existing tested code in place · greenfield `dark_orf_pipeline/` · greenfield porting tested modules |
| D2 | Search strategy | Cascade with permissive narrowing (measured, fits budget) · hub on a reduced source list · hybrid (hub for cheap sources, cascade for expensive) |
| D3 | Compute cuts | Strike HH-suite · structure on representatives only · drop InterPro · eggNOG on the annotated fraction only · provision a GPU |
| D4 | Output format | Parquet + DuckDB views · flat TSV as §25 specifies · SQLite plus TSV exports |
| D5 | Definition of dark | Restore `min_explained` · define a new adjudication rule · other |

---

*Produced by three independent review agents (Claude Opus 5) run 2026-09-16 against
`plans/2026-09-10-pipeline-v2-design.md` at commit-time working tree. Full reports:
`scientific_review.md`, `engineering_review.md`, `devils_advocate_review.md` in this directory.*
