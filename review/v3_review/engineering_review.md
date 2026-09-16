# Engineering review: `plans/2026-09-10-pipeline-v2-design.md`

Scope: implementability as software. Scientific validity and biological threshold choice are
another reviewer's. All repo references are `path:line` at the reviewed commit; spec
references are section numbers from the design document.

Baseline: the repository already contains a working, tested Snakemake implementation
(`workflow/Snakefile`, `workflow/rules/*.smk`, 34 scripts, `src/plasmidann/*`, ~248 test
functions). Section 6 is the stage-by-stage inventory of what that implementation already
covers, what it covers differently, and what the spec adds.

---

## 1. Data model and output format

**1.1 — `final/annotation_complete.tsv` at one row per ORF occurrence is 9-14 GB and is not
a working format.** §25 declares ~120 columns over 9,317,050 rows. Column mix is ~40 numeric
flags (~6 B), ~20 identifiers (`orf_occurrence_id`, 32-hex `protein_id`, targets; ~15 B) and
~15 free-text fields (`diamond_nr_target`, `primary_annotation`, `structural_description`,
`hypothesis_support_summary`, `liabilities`; 40-120 B), giving ~1.0-1.5 kB/row, i.e. **9-14 GB
uncompressed**. Loading it with pandas at default dtypes costs 50-110 GB RAM — feasible once
on a 600 GB node, never twice concurrently. It also exceeds every spreadsheet row limit
(1,048,576), which invalidates the stated rationale for the current CSV form
(`docs/PIPELINE_CODE.md` §14).

**1.2 — The format should be Parquet (or DuckDB over Parquet), not TSV.** The specific reason
is not compression: it is that ~90 of the 120 columns are constant within a `protein_id` or a
`family_id`, so dictionary/RLE encoding collapses exactly the columns that dominate the size,
and predicate pushdown plus column projection lets the common query ("give me the dark families
with a defence hypothesis") read <1% of the file. Expect 0.5-1.5 GB and second-scale queries.
Keep one TSV/CSV export of the *family* table (§26, ~10^5 rows) for spreadsheet use; the
occurrence table cannot serve that purpose at any format.

**1.3 — One row per occurrence duplicates three levels of data.** Protein-level columns
(§25.5-25.9 Pfam/HMM-HMM/Swiss-Prot/nr/eggNOG/InterPro, §25.13 structure, §25.14 properties,
§25.4 dark state — roughly 60 of the 120) are duplicated 9,317,050/3,497,616 = **2.66x**;
family-level columns (§25.10-25.12, ~30 columns, the widest text among them) are duplicated
once per occurrence, i.e. ~19x per family at 5x10^5 families. Storing each fact at its own
level (§4.1 already names the tables) is a 3-4x size reduction before any encoding.

**1.4 — Duplicated aggregates make wrong answers the default.** `family_size`,
`plasmid_count`, `host_count`, `percentage_dark_in_family` (§25.10, §25.12) on an occurrence
row mean any `groupby(family_id)` count, sum or mean is a count of ORFs, not of members —
silently, with no error. The spec's own validation §32.5 (`dark_member_count <= family_size`)
cannot be evaluated on this table without first deduplicating to the protein level. §26 already
specifies the normalized `dark_family_members.tsv`; the aggregates belong there and in
`dark_families_complete.tsv` only, with the occurrence table carrying `family_id` and nothing
else from the family level.

**1.5 — §4 demands a normalized relational model and §25/§38 deliver ten flat files.** Nothing
in the spec bridges them: no stage produces the eleven logical tables of §4.1, no format is
named for them, and §38's "optional" raw/intermediate deliverables are the only place they could
live. Decide once: the internal model is SQLite/DuckDB with the §4.1 tables as real tables, and
§38's TSVs are generated views. Otherwise every agent in §29 invents its own join key.

**1.6 — §5.3's single `NA` token is incompatible with numeric columns in TSV and with §30
Rule 3.** Writing `NA` into `dnds_value`, `pfam_evalue` and `foldseek_score` makes every numeric
column a string on read-back; and `NA` cannot distinguish "test not run" from "test ran, no
result", which §30 Rule 3 and the existing implementation both require (`TOO_FEW_MEMBERS`,
`NO_DIVERGENCE`, `SATURATED`, `NO_OUTPUT` in `workflow/scripts/family_evolution.py:126-140`).
Rule: value columns are typed and nullable, and every value column has a sibling `*_status`
column; `NA` appears only in text columns.

---

## 2. Stage dependency and ordering defects

**2.1 — `dark_evidence_level` is circular.** §13.5 assigns it at Stage 6, but levels 3-5 are
defined as "independently recurrent", "conserved dark family" and "strongly conserved/recurrent
dark protein with independent biological context" — family membership comes from §15 (Stage 7),
context from §16 (Stage 8), recurrence from §14. §37's diagram makes this explicit and worse:
MMseqs2, context and evolution all hang *below* the ANNOTATED/DARK split. The level cannot be
computed where the spec computes it. Fix: split into `dark_status` (Stage 6, sequence evidence
only) and `dark_evidence_level` (recomputed in §23 evidence integration, after families,
context and evolution).

**2.2 — §14's recurrence evidence cannot be computed from what §10 retains.**
`hypothetical_deposit_count`, `independent_database_count`, `independent_source_count` need the
full hit list per protein, but §10.7 specifies only best-hit columns and the implementation caps
reporting at `max_target_seqs: 5` (`config/cascade.yaml:65`). With five nr targets per query,
`hypothetical_deposit_count` saturates at 5 and is meaningless. The spec also never defines what
a "database record" or an "independent source" is. Either raise the target cap for the nr tier
and persist the full m8 (see 4.2 for the disk cost), or delete the fields.

**2.3 — §9's stage number contradicts its input.** §9 is "Stage 4 — Annotation Normalization
Layer", placed before §10's Stage 5 annotation hub, and it normalizes *database annotation
strings* — which do not exist until the hub has run. §37's diagram places normalization after
the hub, correctly. The stage must be renumbered and re-placed, or it will be built by Agent 5
(§29) against an input Agent 4 has not yet produced.

**2.4 — §16's neighbour fields require an ORF-level annotation table that no stage before §25
produces.** `neighbor_annotation` and `neighbor_function_class` (§16.3) are properties of
*annotated* neighbours, while §37 routes context down the DARK branch only. The existing
implementation solves this with an intermediate join at S4 (`workflow/rules/s3_cascade.smk:159`,
`workflow/scripts/context_features.py:41-49` reads `s4/plasmid_annotation.tsv`); the spec has no
equivalent stage between adjudication and the final table. Add one: an ORF-level annotation
table produced immediately after §12, consumed by §16 and §21.

**2.5 — §16.4's context categories have no producer.** The fourteen categories (relaxase,
mobilization, replication, partition, toxin_antitoxin, defense, …) are a *functional role*
vocabulary; §9's vocabulary is a *characterization status* vocabulary
(`FUNCTIONALLY_CHARACTERIZED`, `HYPOTHETICAL_PROTEIN`, …). Nothing in the spec maps a free-text
Pfam or nr label to a role. The repo has a curated map (`src/plasmidann/backbone.py`, ~75 Pfam
families across seven roles, matched by exact family name); the spec needs to either adopt it or
name its own source, or §16.7's seven `*_associated` labels are unimplementable.

**2.6 — §17 forces context into two stages, and the spec describes one.** Per-ORF context (§16)
must precede family construction's consumers, while synteny conservation (§17) and
`dominant_context` (§26) are aggregations *over family members* and must follow it. §15's "run
MMseqs2 after the full annotation phase or at the configured point where all protein sequences
are available" leaves this undetermined. Declare the order: per-ORF context → families →
family-level context aggregation → evidence integration.

**2.7 — §18 needs nucleotides and §6 stores only protein.** §6.3's ORF record has
`protein_sequence` and coordinates; §18.3's RNAcode and §18.4's dN/dS need codon alignments. No
stage in the spec recovers CDS. The repo has one (`workflow/scripts/extract_cds.py`, re-reading
the genome FASTA and honouring `spans_origin`); the spec must add it, including the wrap-around
case, or S7 has no input.

**2.8 — Four stages cannot be sharded and each is a single point of failure late in the run.**
(a) §15 MMseqs2 over 3.5M sequences — one process, declared 64 GB / 24 h
(`workflow/rules/s5_s9_targets.smk:78-97`); (b) §10.4 eggNOG-mapper, unsharded, declared 48 h
over the protein set (`workflow/rules/s3_cascade.smk:179-202`); (c) §12 adjudication and §25
final assembly, global joins over every tool output (`workflow/scripts/annotation_report.py:50-68`
holds six tables in RAM); (d) §21 background normalization, which is global by definition.
Everything else in §37 is shardable and must be sharded (see 4.3).

**2.9 — §37's diagram contradicts §7.3 and §30 Rule 3.** The diagram routes AntiFam-positive
records off the main line ("──> flagged / discovery excluded") and puts "partial / <20 aa
filtering" as an inline stage, so the flow below it carries survivors only. §7.3 and Rule 3 say
these records stay in the master dataset with an exclusion reason. Implemented literally, the
diagram deletes records the text requires. Redraw it as a column-setting stage with no branch.

**2.10 — §19 places structural analysis on "all dark proteins", §25.13 reports it per
occurrence, and no stage reconciles the two populations.** Dark status is a protein property;
the exclusions that define the discovery universe (§7.3: partial, <20 aa) are *occurrence*
properties. A protein can be partial in one occurrence and complete in another. The repo
resolves this explicitly (`workflow/scripts/dark_set.py:21-41`: partial-only proteins excluded,
proteins with at least one complete occurrence kept); the spec does not define it at all.

---

## 3. Identifier and join integrity

**3.1 — `family_id` as specified is not stable and the current implementation's is not either.**
`workflow/scripts/cluster_dark.py:86` assigns `F{i:07d}` from `enumerate(sorted(members))`, so
adding or removing one dark protein renumbers every subsequent family, and §15.5's three
resolutions would produce three mutually unstable namespaces. Derive it from content:
`family_id = <resolution>:<representative protein_id>`, or a hash of the sorted member set. This
is the single most damaging instability in the spec, because §26's family table, §25.10's
carried columns and any external analysis all key on it.

**3.2 — `annotation_id`, `context_id`, `structure_id` (§4.2) are named and never defined.**
Nothing in the spec or the repo emits them. The natural implementation — a row number in a tool
output — is unstable under any re-shard, because a hit's row index depends on the shard's
composition. Define each as a content hash over the fields that identify the observation:
`sha256(protein_id, tool, db_name, db_version, target_accession, query_start, query_end)`,
truncated to 16-32 hex. Tool version must NOT be in the hash, or every tool upgrade renames
every annotation.

**3.3 — `orf_occurrence_id` survives a re-shard but not a gene-caller change.**
`src/plasmidann/orfindex.py:10-15` numbers ORFs `{plasmid_id}|{n}` in coordinate order over the
complete set, so shard count does not affect it (the sort key is `plasmid_id, start, end`) — but
inserting or losing one ORF on a plasmid renumbers every later ORF on it, and `min_orf_aa` or a
pyrodigal upgrade does exactly that. Use coordinates: `{plasmid_id}:{start}-{end}:{strand}`.
Ordinal ids are what §4.2's "do not use array indices as biological identifiers" is warning
against, one level up.

**3.4 — `protein_id` is the one id that is already right; the spec should say so.**
`src/plasmidann/dereplicate.py:13-17` uses `sha256(seq)[:32]`, which is invariant to shard count,
input order, tool version and partial rerun, and the 32-hex width is justified in place against
the birthday bound. §4.2 should mandate content-hash derivation explicitly rather than
"recommended", because every join in §25 goes through it.

**3.5 — Cascade shard membership is positional and a change to `n_cascade_shards` silently
invalidates resumed state.** `workflow/scripts/shard_cascade_input.py:33-37` assigns round-robin
over file order, and explained spans accumulate *within* a shard
(`workflow/rules/common.smk:6-30`). Change the shard count between runs and every protein moves,
while the old shard outputs remain on disk under paths nothing now requests. Record the shard
count in a run-manifest file that is an input to every tier, so a mismatch fails at scheduling
time rather than producing a table assembled from two partitions.

**3.6 — The sweep cohort is stable only if the protein set is.** `workflow/scripts/sweep_cohort.py:47`
draws `rng.sample(ids, k)` from the id list; any change to the ORF set re-draws the whole cohort
and destroys cross-run comparability. Make membership a property of the id, not of the sample:
`int(protein_id[:8], 16) / 16**8 < fraction`. Same fix applies to any future subsampling.

**3.7 — Spiked control records break §32.1 as written.** `workflow/scripts/prepare_control.py`
injects `CTRL_`-prefixed sequences into the cascade query set; they are proteins with no ORF
occurrence, so "every unique protein maps to >=1 ORF occurrence" (§32.1) fails by construction.
The spec has no concept of instrumentation records. Add a `record_class` column
(`observed` | `control`) and scope every §32 check to `observed`.

---

## 4. Resumability and failure handling

**4.1 — The nr tier does not fit a 7-day wall on one node, and the spec never budgets it.**
64 cascade shards at 4 threads each is 24 concurrent jobs on 96 cores, i.e. three waves; the
declared per-shard runtime for an nr tier is 48 h (`workflow/rules/s3_cascade.smk:121`), giving
up to 6 days for T4 alone before T1-T3, S6-S9 and the final joins. The nr fixed cost per
invocation is explicitly unmeasured (`config/config.yaml:24-26`, `workflow/bench_nr.sbatch`), and
sharding multiplies it 64x. Run the bench before committing to the shard count, and submit the
cascade as a SLURM job array with per-shard resume rather than one 7-day allocation.

**4.2 — §5.1 and §30 Rule 4 are violated today and have an unbudgeted disk cost when fixed.**
`workflow/scripts/tier_search.py:113,121,160` writes hmmsearch `--domtblout` and DIAMOND m8 into
`tempfile.mkdtemp()` and discards them. Preserving raw output as §5.1 requires means ~17.5M m8
lines per DIAMOND tier at `max_target_seqs: 5` (~2-4 GB/tier) plus Pfam domtbl at both tiers plus
AntiFam — on the order of 15-30 GB of raw, before §2.2's larger target cap. Declare a
`raw/` retention policy with an explicit budget and a gzip step, or state which raw outputs are
deliberately not kept.

**4.3 — Four expensive stages have no shard wildcard and lose everything on one failure.**
Required sharding granularity, in order of exposure:
 * eggNOG-mapper (§10.4): shard by protein, 64 shards, same partition as the cascade.
 * Foldseek/ProstT5 (§19): shard by dark protein, 32-64 shards; `workflow/scripts/structure_search.py`
   is one job over the whole dark set today (`workflow/rules/s5_s9_targets.smk:264-285`).
 * Family evolution (§18): shard by family block, 64-256 shards. The rule declares `threads: 16`
   (`workflow/rules/s5_s9_targets.smk:134`) and the script runs one family at a time with
   `mafft --thread 1` (`workflow/scripts/family_evolution.py:157`) — 16 cores idle for the whole
   stage. This is the same defect class that got S7d removed (`docs/PIPELINE_CODE.md` §13).
 * MMseqs2 clustering (§15): genuinely unshardable. Give it its own job with
   `--split-memory-limit`, a persistent tmp directory on scratch, and an explicit restart note.

**4.4 — §30 Rule 3 requires a rejected-hit record that no stage writes.** Every filter must emit
`filter_reason` and `filter_stage`. `workflow/scripts/tier_search.py:56,70` counts significance
rejections into `n_rejected` and only prints it; the rejected hits vanish. Honouring Rule 3 means
a `rejected_hits` table per tier-shard, which is where most of the extra volume in 4.2 comes
from. Decide explicitly: keep rejections as rows, or amend Rule 3 to "counted and logged".

**4.5 — The last job in the run is the largest single point of failure.** §25's assembly reads
every table into memory (`workflow/scripts/annotation_report.py:50-68`) and streams 9.3M rows
out. A failure there, on day six, re-does nothing but costs a resubmission and a full re-read.
Emit it per plasmid-shard and concatenate, or write Parquet row groups per shard — either makes
it resumable and removes the 50-110 GB read-back problem from 1.1.

**4.6 — Resume correctness depends on `--rerun-incomplete` being used every time.**
`docs/PIPELINE_CODE.md` §1 says so; `workflow/run_pipeline.sbatch` is the only enforcement. A
killed job's partial TSV is otherwise indistinguishable from a complete one, since no stage
writes a completion sentinel or a row count. Add a `.done` sidecar carrying the expected row
count, and check it on read — this is cheap and it converts a silent truncation into an error.

---

## 5. Configuration and validation

**5.1 — Config keys the spec requires but leaves undefined.** Literally written as
"configurable"/"conditional" in §31: `pfam.relaxed_evalue`, `mmseqs.family_thresholds.close`,
`.intermediate`, `.broad`, `evolution.dnds: conditional` (a predicate, not a value; the predicate
is never given — §18.4's "enough independent homologous sequences" and "adequate alignment
quality" are undefined quantities). Required by the text and absent from §31 entirely: InterPro
(§10.3), eggNOG (§10.4), DIAMOND Swiss-Prot and nr thresholds and target caps (§10.5, §10.7),
HMM-HMM probability/coverage cuts (§10.6), annotation-transfer criteria (§12.4),
`dark_evidence_level` boundaries (§13.5), "convincing structural match" (§19.4), RNAcode p
(§18.3), §17's conservation thresholds, §21's normalization method, §22's rarity band
boundaries, §23's independence rules, database paths and versions (§28 names `databases.yaml`
and §31 shows none), and **shard counts and thread counts — a 9.3M-ORF spec with no parallelism
configuration at all**.

**5.2 — §31 contradicts §7.3 on the length floor, and the repo has a third value.** §7.3 and
§7.4 say 20 aa; §31 says `min_dark_length_aa: 30`; `config/config.yaml:37` sets `min_orf_aa: 20`
as a *gene-caller* floor, which means ORFs below it are never called and therefore cannot be
"labelled rather than deleted" as §7.4 demands. Three numbers and two different semantics. Pick
one number, and separate `orf.min_length_aa` (calling) from `discovery.min_length_aa` (labelling)
so §7.4 is satisfiable.

**5.3 — §32 checks that cannot be written as tests, because the quantity is never defined.**
 * §32.3 "primary_annotation exists only when evidence criteria are met" — criteria undefined (§12.2).
 * §32.3 "hypothetical labels are not interpreted as functional assignments" — no predicate given.
 * §32.2 "translated ORFs have expected lengths" — "expected" undefined, and it is false as stated
   for origin-spanning genes where `end < start` (`src/plasmidann/circular.py`, 160,375 genes).
 * §32.6 "genus_count <= host_count where definitions require it" — the spec itself defers the
   definition; the check is untestable until the units are fixed.
 * §32.5 `percentage_dark_in_family` — requires a family-level dark definition when dark is a
   protein property and the discovery exclusions are occurrence properties (see 2.10).
 * §32.4 "excluded records remain recoverable" — no definition of recoverable, no artefact named.

**5.4 — §30 Rule 10 (record software, version, database, database version/date, parameters) has
no implementation surface.** `workflow/scripts/preflight.py` records resolved paths and database
*sizes* only; no versions, no dates. Every output carries `thr_*` threshold columns
(`workflow/scripts/tier_search.py:188-196`) but nothing else from Rule 10. Add one
`run_manifest.tsv` written by preflight — tool name, `--version` output, database path, mtime,
size, config hash, git commit — and make it an input to every rule so it cannot go stale.

**5.5 — §31's example config would not validate against the repo's own schemas.**
`config/schemas/cascade.schema.yaml` sets `additionalProperties: false` and requires
`narrow_at`, `min_explained`, `min_coverage`, `full_at`, `partial_at`, `hmmer_z`,
`max_target_seqs`, `sweep_cohort_fraction`, `tiers`, `artefact_screen` — none of which appear in
§31, and §31's keys (`orf.exclude_partial`, `annotation.retain_all_hits`, `output.missing_value`)
appear in no schema. Adopting §31 verbatim means rewriting both schemas and the coherence checks
in `src/plasmidann/cascade.py::check_thresholds` and `plasmidann.targets::check_reality_config`,
which are wired into the Snakefile at load time (`workflow/Snakefile:29-56`).

---

## 6. Spec-to-code inventory

Legend: **HAVE** = implemented and consistent with the spec; **DIFFERS** = implemented, but
disagrees with the spec as written; **NEW** = nothing in the repo implements it.

| Spec | Stage | Status | Where / what disagrees |
|---|---|---|---|
| §3, §4.1 `plasmids` | metadata ingest, `plasmids.tsv` | DIFFERS | `workflow/scripts/analysis_set.py` emits an id list only (8 lines); `workflow/scripts/clonal_registry.py:30-32` emits 6 columns. No `raw_*`/`normalized_*` pairs (§3.3), no metadata validation report, no `plasmids_complete.tsv`. |
| §6 | ORF prediction | DIFFERS | `workflow/scripts/orf_call.py` + `src/plasmidann/circular.py` (origin repair, not in the spec at all). Columns written are 7 (`orf_call.py:52`); spec §6.3 asks for 11 incl. `protein_id_pre_derep`, `start_type`, `nucleotide_length`, `protein_length`. |
| §7.2 AntiFam | artefact screen | DIFFERS | `workflow/scripts/artefact_screen.py:123` writes `seq_id, artefact_flag, antifam_family, antifam_ievalue, …` plus tantan low-complexity. Spec asks for `antifam_score`, `antifam_description`; spec has no low-complexity screen. |
| §7.3-7.4 | discovery exclusions | DIFFERS | `<20 aa` is a gene-caller floor, not a flag (`config/config.yaml:37`; see 5.2). `discovery_excluded_short` / `discovery_excluded_partial` columns do not exist; partial-only exclusion is implicit in `workflow/scripts/dark_set.py:21-41`. |
| §8 | dereplication | HAVE | `workflow/scripts/dereplicate.py`, `src/plasmidann/dereplicate.py`; losslessness asserted (`dereplicate.py:9`), tested (`tests/test_dereplicate.py`). |
| §9 | annotation normalization layer | NEW | Nothing implements the 10-class vocabulary, `normalization_rule`, `normalization_version`, or `normalization_dictionary.tsv`. Nearest is the binary `UNINFORMATIVE` regex in `src/plasmidann/cascade.py` (+ `tests/test_uninformative_labels.py`, 25 labelled cases). |
| §10.1 T1 Pfam GA | annotation | HAVE | `config/cascade.yaml:86-93` (T1, `--cut_ga`), `workflow/scripts/tier_search.py:120-152`. |
| §10.2 T2 Pfam relaxed | annotation | HAVE | `config/cascade.yaml:95-116` (`-E 1e-5 --domE 1e-5 --incdomE 1e-5`). |
| §10.3 T3 InterPro | annotation | NEW | No InterProScan anywhere; not in `src/plasmidann/tools.py::REQUIRED_TOOLS`. |
| §10.4 T4 eggNOG | annotation | DIFFERS | `workflow/scripts/orthology.py` runs eggNOG-mapper on **only the proteins the cascade named** (`workflow/rules/s3_cascade.smk:179-202`, `tests/test_scripts_smoke.py:1167`). Spec says the whole ORF protein set. Also unsharded. |
| §10.5 T5 Swiss-Prot | annotation | HAVE | `config/cascade.yaml:118-125` (tier T3). Column names differ (`annot_label`, `annot_evalue` vs `diamond_swissprot_*`). |
| §10.6 T6 HMM-HMM | annotation | NEW | No HH-suite; `docs/PIPELINE_CODE.md` §15 lists a missing phage tier, not this. |
| §10.7 T7 broad DB | annotation | HAVE | `config/cascade.yaml:127-135` (tier T4, nr). Capped at `max_target_seqs: 5` — see 2.2. |
| §11 | search/evidence/reporting thresholds | DIFFERS | Two-way split implemented (`narrow_at` vs `min_explained`, `config/cascade.yaml:32-33`, applied at `workflow/scripts/cascade_resolve.py:118`), not the spec's three-way split. |
| §12 | adjudication | DIFFERS | `src/plasmidann/cascade.py::classify/completeness` produce `functional_class`, `homology_depth`, `annot_completeness`. No `primary_annotation` / `_source` / `_confidence`, no §12.1 8-level hierarchy, no §12.4 transitive safeguard. |
| §13 | dark classification | DIFFERS | Dark = `functional_class in ("UNCHARACTERIZED_HOMOLOG","NONE")` (`workflow/scripts/quality_gate.py:103`). `dark_evidence` ladder exists (`src/plasmidann/cascade.py`, `tests/test_dark_evidence.py`); the numeric 0-5 `dark_evidence_level` and its sublabels do not. |
| §14 | recurrence evidence | NEW | Nearest is `n_dark_databases` (`cascade.py`, `tests/test_dark_coverage.py`) and `dark_covered_fraction`. None of the five §14 fields exist, and see 2.2 for why they are not computable as specified. |
| §15 | MMseqs2 families | DIFFERS | `workflow/scripts/cluster_dark.py` clusters **only the dark set**, spec §15.1 says all ORFs. Single resolution (`config/targets.yaml:50-54`) vs §15.5's three. Family columns are 8 (`cluster_dark.py:69`); §15.3 asks for 11 more (`host_count`, `genus_count`, `habitat_count`, `percentage_dark_in_family`, `max_annotation_strength`, …). `family_id` unstable (3.1). |
| §16.1-16.3 | neighbourhoods | DIFFERS | `src/plasmidann/context.py` + `workflow/scripts/context_features.py`; window is ±3 (`config/targets.yaml:109`), spec wants ±10 retained and ±3 summarized. No per-ORF neighbour table with the seven §16.3 fields — only family-level conservation/enrichment. |
| §16.5 DefenseFinder | context | HAVE | `workflow/scripts/defence_search.py` → `defence_gembase.py` → `defence_systems.py` (MacSyFinder driven directly with replicon topology); `src/plasmidann/defence.py`, `tests/test_defence.py`. |
| §16.6 IntegronFinder | context | HAVE | `workflow/scripts/integrons.py` (`--local-max`), sharded over 600 plasmid shards. |
| §16.7 | context labels | DIFFERS | Six features produced (`context_features.py:147`: defence, integron, backbone_adjacent, annotated_neighbour, operon_with_annotated, ta_candidate). Spec names seven different ones (mobilization, replication, partition …) and provides no role map — see 2.5. |
| §17 | synteny conservation | DIFFERS | `context_conservation` + `enrichment` over a corpus background (`src/plasmidann/context.py`). None of `left_neighbor_conservation`, `right_neighbor_conservation`, `operon_like_conservation`, `context_recurrence`, `synteny_conservation` exist. |
| §18.1 | distribution counts | DIFFERS | Only `n_plasmids`, `n_mob_clusters`, `n_orfs` (`cluster_dark.py:69`). `host_count`, `species_count`, `genus_count`, `lineage_count`, `habitat_count` are absent, though the master table carries the source columns. |
| §18.2 | sequence conservation | NEW | `mean_pairwise_identity`, `identity_distribution`, `alignment_coverage`, `conserved_positions`, `conserved_motifs` — none exist; the alignment they would come from is already built (`family_evolution.py:157`). |
| §18.3 RNAcode | evolution | HAVE | `workflow/scripts/family_evolution.py:95-124`, both strands, Clustal W workaround for RNAcode's exit-0-on-error. |
| §18.4 dN/dS | evolution | HAVE | `src/plasmidann/evolution.py` (Nei-Gojobori), `tests/test_evolution.py`; status codes instead of blanks. Spec's `selection_model`, `selection_significance` absent (HyPhy removed, `docs/PIPELINE_CODE.md` §13). |
| §19 | structure | DIFFERS | `workflow/scripts/structure_search.py` (Foldseek + ProstT5). No `structure_model_id`, `structure_confidence`, `model_length` — ProstT5 emits no pLDDT, documented at `config/targets.yaml:124-139`. State names differ from §19.3. Unsharded. |
| §20 | properties / liabilities | DIFFERS | `src/plasmidann/peptide.py`: net charge, GRAVY, TM-window heuristic only. Missing `molecular_weight`, `isoelectric_point`, `disorder_fraction`, `signal_peptide`, `coiled_coil`, `secondary_structure_summary`. `low_complexity_fraction` exists only as a boolean flag in `artefact_screen.py`. |
| §21 | distribution / background normalization | DIFFERS | Only per-feature corpus background rates (`context_features.py:147-193`). No per-gene or per-plasmid normalization, no size-matched backgrounds, no `raw_prevalence`/`normalized_prevalence` pair. |
| §22 | rarity vs conservation labels | NEW | No `rare` / `lineage_specific` / `cross_MOB` / `widely_conserved` labels anywhere. |
| §23 | evidence integration | DIFFERS | `src/plasmidann/targets.py::reality_lines` counts 4 named boolean lines and reports which fired. Spec wants 9 named categories, `evidence_count_total`, `independent_evidence_categories`, and explicit non-independence handling (§23.2). |
| §24 | functional hypothesis layer | DIFFERS | `top_hypothesis`, `top_conservation`, `top_enrichment`, `high_confidence` (`context_features.py:152-186`). No `hypothesis_support_summary` free-text trace. |
| §25 | `annotation_complete` | DIFFERS | `workflow/scripts/annotation_report.py:131-152` — CSV not TSV, ~45 columns not ~120, no per-source column blocks (§25.5-25.9). |
| §26 | `dark_families_complete` | DIFFERS | `annotation_report.py:73-111`, 34 columns. `dark_family_members.tsv` does not exist (the member list is a comma-joined cell, `cluster_dark.py:90`). |
| §27 | `plasmid_annotation_summary` | NEW | Nothing aggregates to the plasmid. |
| §28 | repo layout `dark_orf_pipeline/` | DIFFERS | Current layout is `src/plasmidann/` + `workflow/` + `config/` + `tests/`, with a deliberate and documented split (`docs/PIPELINE_CODE.md` §2). §28 adds `agents/`, `scripts/`, `data/{raw,intermediate,normalized}`, `results/{...}`. A rename touches every `conda:`/`script:` path and every test import — see the note below the table. |
| §29 | 12-agent decomposition | NEW | No equivalent; contracts between stages are file paths in `workflow/rules/*.smk`. |
| §30 Rule 9 | schema + validation test per table | DIFFERS | Tests are per-concern, not per-table (`tests/`, ~248 functions). Only the two config files have schemas (`config/schemas/`); no output table has one. |
| §30 Rule 10 | provenance stamping | DIFFERS | `thr_*` columns only (`tier_search.py:188`); `preflight.tsv` records paths and db sizes, not versions/dates. See 5.4. |
| §32 | validation stage | DIFFERS | Inline assertions (`dereplicate.py:9`, `annotate_plasmids.py:77-78`, `orf_call.py:104`) plus the run-halting S5 control gate. No standalone validation stage, no `pipeline_qc_report.tsv`. |
| — | **positive control / quality gate** | **HAVE, and the spec deletes it** | `workflow/scripts/prepare_control.py` + `quality_gate.py` + `config/targets.yaml:10-47`: 500 reviewed Swiss-Prot proteins spiked into the query set, run through the identical code path, run halts below 0.99 recall. **The spec has no positive control anywhere.** Implementing §29 from the spec alone removes the only test that can detect an annotation-recall failure. |
| — | preflight | HAVE, unspecified | `workflow/scripts/preflight.py` + `src/plasmidann/tools.py` + `tests/test_tools_registry.py`. Not in the spec; keep it. |
| — | circular-origin repair | HAVE, unspecified | `src/plasmidann/circular.py`, 160,375 genes affected. The spec never mentions topology; without this stage those genes enter the dark set as fragments. |
| — | feature files (GFF3/GenBank) | HAVE, unspecified | `workflow/scripts/feature_files.py`, `src/plasmidann/features.py`. |
| — | consensus re-check | HAVE, unspecified | `workflow/scripts/consensus_recheck.py` (family consensus re-searched at Pfam GA). |
| — | S9 prioritisation / library design | HAVE, and §1.2/§36 forbid it in v1 | `workflow/scripts/prioritise.py`, `library_design.py`, `src/plasmidann/targets.py`, `config/targets.yaml:146-238`. Opt-in via `snakemake portfolio` (`workflow/Snakefile:84-89`), so it already sits outside `rule all` — keep it out, do not delete it. |

**On §28.** The rename to `dark_orf_pipeline/src/dark_orf/{io,orf,antifam,…}` is pure churn
against CLAUDE.md's surgical-changes rule: it invalidates every `script:` path in three `.smk`
files, the `_ctx.py` import mechanism (`workflow/scripts/_ctx.py`), `pytest.ini`'s
`pythonpath = src`, and the script-harness in `tests/conftest.py:22-24`, while changing no
behaviour. The one part worth adopting is the `data/{raw,intermediate,normalized}` and
`results/final/` separation that §5.1 actually requires; that can be done by adding directories
under the existing `results/` tree.

---

## 7. Testability of the §38 deliverables

**7.1 — `final/pipeline_qc_report.tsv` has no schema and cannot be tested at all.** §38 names it
and §32 lists the checks, but no section defines the report's columns, the representation of a
passed vs failed check, or whether a failure halts. Without that, Rule 9 cannot be satisfied for
the very artefact that implements Rule 9. Define it as one row per check:
`check_id, stage, scope, expected, observed, status, severity`, and make a `FAIL` at
`severity=halt` exit non-zero.

**7.2 — `plasmid_annotation_summary.tsv` cannot have a sum test, because §27's counts do not
partition.** §27 lists `annotated_gene_count`, `uncharacterized_gene_count` and `dark_gene_count`
side by side, while §2.5 and §13.3 define uncharacterized as a *subset* of dark. So
`annotated + uncharacterized + dark != gene_count` and the natural invariant is unwritable.
Either rename to `dark_hypothetical_gene_count` as a declared subset with its own
`<=` assertion, or drop one column.

**7.3 — `annotation_complete.tsv` is testable for structure but not for its headline column.**
Writable: row count equals occurrence count; `orf_occurrence_id` unique; every `protein_id`
present in the protein table; every `family_id` present in the family table; every family-level
column identical across all rows sharing a `family_id` (this is the test that catches the
denormalization bugs in 1.4); missing-value token is exactly `NA` in text columns.
Not writable: `primary_annotation_confidence` has no defined scale or domain (§12.2), so no
assertion about it exists beyond "non-empty".

**7.4 — `dark_families_complete.tsv`: three of the §26 columns are untestable as specified.**
Writable: §32.5's three inequalities; `family_id` unique; member counts equal the row count in
`dark_family_members.tsv`; counts monotone across §15.5's nested resolutions
(`close <= intermediate <= broad`). Not writable: `max_annotation_strength` (no scale defined),
`dark_evidence_level_summary` (no aggregation function defined), `context_conservation` (§17
names six measurements and defines none of them numerically).

**7.5 — `protein_annotations_complete.tsv` and `orf_occurrences_complete.tsv` are fully
testable, with one exception.** Assert: every protein maps to >=1 occurrence and every occurrence
to exactly one plasmid (§32.1); sequence alphabet valid (§32.2); no duplicate id with a different
sequence (§32.2, which for a hash-derived `protein_id` is a collision check and should be
asserted, not assumed — `src/plasmidann/dereplicate.py:5-12` explains why). The exception is
§32.2's "translated ORFs have expected lengths": false for the 160,375 origin-spanning genes
where `end < start` unless a `spans_origin` flag is specified, which §6.3 does not do.

**7.6 — `normalization_dictionary.tsv` is the best-specified deliverable and should be built
first.** Assert: every distinct `raw_annotation` observed in the corpus maps to exactly one
class; every class is in §9.2's closed vocabulary; every row carries a resolvable
`normalization_rule` and `normalization_version`; and — the test that matters — recall against a
labelled set, as `tests/test_uninformative_labels.py` already does for the binary case. Extend
that 25-case set to the 10-class vocabulary and make CI fail on a recall drop.

**7.7 — `plasmids_complete.tsv` is testable; add the pairing assertion.** Row count equals input
FASTA record count; `plasmid_id` unique; every `normalized_<field>` has its `raw_<field>` sibling
present (§3.3) — this one assertion prevents the normalization layer from overwriting source
metadata, which is §30 Rule 5.

**7.8 — No deliverable in §38 can test annotation correctness, and the spec provides no gold
set.** §39's success criteria are all structural ("every ORF has a stable identity", "AntiFam
status is available"). The only correctness instrument in the project is the spiked Swiss-Prot
control gate that the spec omits (see the inventory's last-but-five row). Restore it as a
§32 requirement: N reviewed proteins of known function, spiked pre-cascade, traversing the
identical code path, with a declared minimum recall that halts the run.
