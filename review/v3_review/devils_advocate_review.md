# Devil's advocate review — v3 spec, "what happens after it works"

Target: `plans/2026-09-10-pipeline-v2-design.md` as it now stands (2,304 lines, the "Dark ORF
Discovery Pipeline — Agent Implementation Specification").

Premise of this review, as briefed: **assume the pipeline builds, runs inside the budget, and emits
every table in §38 without a single failed assertion.** Nothing below is a bug report. Two other
reviewers hold scientific validity/thresholds and software engineering/DAG; I do not touch those.
My question is only: given a clean run, why is the result not worth what it cost?

Numbers are taken from `docs/annotation_statistics.md`, `docs/PARAMETER_PROVENANCE.md`,
`docs/PIPELINE_CODE.md`, `config/*.yaml`, `src/plasmidann/*`, and from
`reports/dark_orf_clustering.md` as it existed at `b095493^` (deleted in `b095493`). Nothing was
executed.

---

## 0. The finding that reframes all the others

**The spec is a regression against code that already exists in this repository and is better than
it.** `docs/annotation_statistics.md` names this file as its companion, and every measured decision
in that document — `-Z 3497616`, `--domE`/i-Evalue, AntiFam `--cut_ga`, the T1/T2 authority split —
is absent from the spec it is a companion to. `config/cascade.yaml`, `config/targets.yaml` and
`docs/PARAMETER_PROVENANCE.md` describe a pipeline with a narrowing threshold separated from a
reporting threshold, a 2% sweep cohort that preserves the counterfactual, a run-halting spiked
positive control at `min_control_recall: 0.99`, circular-origin repair for 160,375 measured ORFs,
RNAcode on both strands, a family-consensus re-check, and a clonal registry. **None of these appear
anywhere in the 2,304 lines under review.** The previous devil's advocate review listed ten of these
in its §8 as "what I tried to break and could not — do not touch these in the course of fixing the
rest." The spec has removed every one of them. A referee comparing this document to the repository's
own `docs/` will conclude the project got less rigorous over time, and they will be right.

---

## Findings, most damaging first

### 1. Nothing in the spec can fail. §39 is a software checklist wearing a science label.

All 22 bullets of §39 ("Definition of Success") are presence-of-column or internal-consistency
statements: "AntiFam status is available", "structural status is available where attempted", "every
final row is auditable". §32's validation is the same: §32.5 checks `0 <= percentage_dark_in_family
<= 100`; §32.6 checks `MOB_count <= plasmid_count`. These are tautologies over the code's own
arithmetic. There is **no external referent anywhere in the spec** — no held-out set of known plasmid
proteins that must come back ANNOTATED, no negative set that must come back excluded, no recall
measurement of any kind. The repo already has this: `config/targets.yaml` declares a spiked
Swiss-Prot control at `min_control_recall: 0.99` (cited to ECLIPSE, 99.2–100%) that **halts the
run**. The spec deletes it. Consequence after a clean run: the headline number "N dark ORFs" is
uncontradictable, because DARK is defined in §2.5/§13.3 purely as "not resolved by the workflow we
chose to run".

### 2. `independent_evidence_categories` (§23.3) counts bins, not independent observations.

§23.2 correctly forbids counting Pfam + InterPro + HHsearch + DIAMOND + MMseqs as five observations
of one homology. It stops there. Of §23.1's nine categories, `distribution`,
`evolutionary_conservation`, `genomic_context` and `orf_qc` are all monotone in the same latent
variable — family size / occurrence count. The previous review measured the consequence in the
implemented version: `purifying_selection ⊆ is_family` deterministically (dN/dS needs
`min_members_for_dnds: 3`, which *is* the family condition), and 79.2% of families with ≥10 members
span ≥2 MOB lineages against 43.8% of families with ≥2. `docs/PARAMETER_PROVENANCE.md` records the
repair — `targets.IMPLIED_BY`, an explicit entailment map, plus `check_reality_config` refusing a
config where the two thresholds diverge. **§23.3 reintroduces the defect with no entailment map at
all.** The fix costs nothing: emit the 2^9 joint indicator table and the pairwise correlation matrix
across categories. §32 requires neither, so the number will be reported as "independent" and no
observation can show it is not.

### 3. `dark_evidence_level` (§13.5) is a scale calibrated after seeing the data, used for nothing.

Six levels, 0–5, of which §13.5 says "the exact numerical boundaries should be treated as
configuration and revised after inspecting the dataset" and "do not use this level as the final
experimental score." A quantity whose boundaries are set post hoc, whose use is forbidden, and which
has no external referent, cannot be wrong — and cannot be useful. Read the semantics: level 2
("repeated sequence/homology"), 3 ("independently recurrent"), 4 ("conserved dark family"), 5
("strongly conserved/recurrent with independent biological context") are four re-labellings of
family size and occurrence count. That is one measurement presented as a ladder. §26 then propagates
it to the family table as `dark_evidence_level_summary`, where the same confound is aggregated.
Either delete it or replace it with the raw counts it is a function of; keeping both is how a
reviewer finds a circular column in your methods table.

### 4. The deferral is not neutral: §16.1 and §19.1 make later selection impossible without a re-run.

§1.2/§30 Rule 8/§35 forbid any ranking, which is defensible. But §16.1 scopes context extraction to
"**for every dark ORF**", and §21.5 rules that "exact statistical tests belong to a later analysis
module". Those two together are fatal. An enrichment statistic needs a background — the neighbourhood
composition of *non-dark* ORFs on size-matched plasmids — and the spec never computes it. §21.2's
denominators are plasmid-level (`plasmid_gene_count`, `host_taxonomy`, `MOB_cluster`, `habitat`,
`source_database`); none of them is a per-ORF neighbourhood. So the downstream prioritisation in §36,
which explicitly wants "dark proteins with conserved genomic context", must re-run Stage 8 over the
remaining ~7–8M annotated ORFs to obtain the denominator. That is the single most expensive stage in
the spec, re-run in full, because a scoping decision in §16.1 saved a fraction of it.

### 5. §6.3 discards nucleotide sequence, so every evolutionary measurement is one-shot and unrepeatable.

§6.3's required ORF fields end at `protein_sequence`. No nucleotide sequence, no CDS. §18.3 (RNAcode)
and §18.4 (dN/dS) both need codon alignments, and both are hedged — "where sufficient homologous
nucleotide sequence is available", "only when... the evolutionary comparison is appropriate",
`dnds: conditional` in §31. So whether a family has an evolutionary measurement is decided at run
time by an unstated rule, and a family judged "not appropriate" has a blank. To obtain it later you
must re-parse 143,504 plasmids and re-extract CDS — and you must do it *correctly across the origin*,
which §6 has no mechanism for (see §7 below). This directly contradicts §1.3's promise that
"reporting thresholds and experimental-selection rules can be changed later without rerunning the
entire analysis". The repo's `extract_cds.py` and `plasmidann.circular` already solve this; the spec
does not know they exist.

### 6. The compute bottleneck is the one stage whose output §23.2 forbids counting.

§37's diagram fans all seven annotation tools out of one hub **in parallel, with no cascade and no
narrowing**; §11.2 says "run broad/permissive searches where computationally practical". So every
tool sees all 3,497,616 unique proteins. Budget: 96 cores × 168 h = **16,128 core-hours**.
- §10.6 HH-suite needs an HHblits MSA per query. At an optimistic 10 core-s/query that is **9,715
  core-hours (60% of the entire budget)**; at a realistic 40–160 core-s/query, **39,000–155,000
  core-hours, i.e. 2.4×–9.6× the whole allocation for one tool**.
- §10.3 InterProScan at 200–750 seq/core-hour is **4,660–17,500 core-hours (29%–108%)**.
- §10.7 DIAMOND vs nr un-narrowed: `docs/PIPELINE_CODE.md` §15 records this tier as "dominates the
  run" and **still unbenchmarked** (`workflow/bench_nr.sbatch` needs an allocation).
And §10.6 and §23.2 both state that HMM-HMM, Pfam, InterPro, DIAMOND and MMseqs "may all derive from
one underlying homology relationship" and must not count separately. **The largest consumer of the
budget buys evidence the spec itself declares non-independent.**

**Cut:** run §10.6 and §10.3 only on the residue after adjudication, not before it — §37 has the
order backwards. For a protein already FUNCTIONAL at Pfam GA, the HH-HMM result changes nothing by
the spec's own §10.6. Restore the repo's `narrow_at: 0.9` / `min_explained: 0.5` split plus the 2%
`sweep_cohort_fraction` so the narrowing's cost is still measurable. **Spend it on:** (i) RNAcode on
both strands for the whole dark set (the only instrument that separates a shadow ORF from a protein,
and the antisense signal is the diagnostic one — §18.3 never mentions strand); (ii) neighbourhood
extraction for all 9,317,050 ORFs so §21's enrichments have a background; (iii) plasmid-level
mash/ANI clustering over 143,504 plasmids (hours, not days) to give every count in §14 and §18.1 a
real independence denominator.

### 7. §7.1 freezes the artefact screen at 278 profiles and forbids improving it.

"Version 1 uses only Prodigal/Pyrodigal output properties and AntiFam. **No additional artifact
detector should be introduced without a deliberate pipeline revision.**" AntiFam is 278 profiles and
is by construction a blocklist of artefact families somebody already rediscovered and reported — the
spec's own companion doc says so. Artefact classes that screen misses, with magnitude for this
corpus:

| class | why AntiFam+partial+<20aa misses it | magnitude |
|---|---|---|
| shadow/antisense ORFs over real genes | not in the 278; conserved, so every recurrence test fires | plausibly **5–15% of the dark set**; detectable **for free** from §6.3 coordinates (opposite-strand overlap fraction), never computed |
| origin-spanning genes on circular plasmids | §6 calls Prodigal on linear FASTA; no `topology` handling | **160,375 ORFs measured on this exact set (1.72%, 1.12/plasmid)**. §7.3 excludes partials, so these real genes are deleted *and* their chance-complete fragments enter the dark set as short novel ORFs — both error directions at once |
| host chromosomal / chromid contamination | no plasmid verification anywhere; §3.3 keeps whatever the source claimed | 5–15% typical for aggregated public plasmid collections; these ORFs are the **most** conserved, so they dominate any conservation-based selection |
| database redundancy / re-deposition | §14 converts re-deposition **into evidence of reality** (`independent_database_count`), and §3.2's four sources overlap by design | one `mob_cluster` (AA379) holds **21,817 plasmids, 15.2% of the corpus** |
| IS/transposon-derived ORFs | multi-copy and multi-backbone by construction; §16.4 has no IS/transposon class | unquantified; every distribution measure in §18.1 fires on them |
| frameshift/assembly-split genes | both halves get fresh start/stop and are "complete" and unnamed | detectable from data the spec already has (adjacent same-strand ORFs whose concatenation hits one target); never done |
| phage / phage-plasmid genes | no phage database; `PIPELINE_CODE` §15 records the PHROGs download failing on every mirror | dark for want of a database, not for want of characterisation |

Composite estimate: **30–50% of the "dark ORFs" this pipeline reports are not novel plasmid protein
biology**, and the recurrence-weighted subset that selection will favour is *worse*, not better,
because classes 1, 3, 4 and 5 are all recurrent by construction. §34.3 states the opposite in prose
— repeated hypothetical deposition is "high-confidence evidence of protein existence" — with no
mechanism distinguishing a recurrent protein from a recurrently copied mis-annotation.

### 8. There is no evidence for ~1,000, and the one measurement that would settle it is never required.

The spec never states a target family count, never estimates how many distinct dark families this
dataset yields, and §36 defers composition to "after inspecting the complete dataset". The project's
own deleted measurement (`reports/dark_orf_clustering.md` at `b095493^`, one subset only) gives:
378,552 dark ORFs → **92,752 families**; 34,817 with ≥2 members; **5,028 dark-only families with ≥10
members**; 15,262 spanning ≥2 MOB lineages. So the candidate pool is 10⁴–10⁵ and **1,000 is a 1–10%
sampling ratio set by the synthesis budget, not by the data.** The opposite risk is never tested:
after the artefact classes in finding 7 and a real independence unit (finding 10), the count of
families that are certainly proteins, genuinely uncharacterised, and tractable in *E. coli* could be
closer to 10². **Where it would become visible:** the family-size distribution of §15.3 and the table
of §26. **Does the spec check?** No — §26 has no summary, §38 lists no distribution figure, §32.5's
family checks are inequalities that cannot fail, and §39 does not ask for it. Add now, at zero cost:
a required rarefaction curve (dark families discovered vs plasmids sampled) at each of §15.5's three
clustering resolutions. It answers both "is 1,000 the right order" and "is this space saturating or
is it database redundancy" — and it is the first figure a referee will ask for.

### 9. §15.4's "especially interesting" set is dominated by singletons, by construction.

"Families with `percentage_dark_in_family = 100` are especially interesting." Every dark ORPHAN — a
family of one — is at 100% by arithmetic. §32.5 validates `0 <= percentage <= 100` and
`dark_member_count <= family_size`, neither of which can fail. On the measured subset, 15.3% of dark
ORFs (57,935) sit in families of one; 91,159 of 92,752 families are dark-only. So the flag selects
~14,000+ singletons and calls them the interesting stratum. §22.1 argues the right principle ("a
protein in many plasmids is not automatically more interesting than one in seven") and §15.4 then
violates it in the opposite direction. The field is not wrong; it is uninformative, and a referee
will read "100% dark family" in a figure legend as a claim rather than a tautology.

### 10. Every "independence" count in §14, §18.1 and §22 rests on a lineage variable that is not a lineage.

`MOB_count` (§18.1, §25.12), `cross_MOB` (§22), `independent_plasmid_count` /
`independent_source_count` (§14) are the spec's entire defence against clonal pseudo-replication.
Measured on this corpus: 7,040 distinct `mob_cluster` values over 143,503 plasmids, **median cluster
size 6**, and **AA379 alone = 21,817 plasmids (15.2%), containing 24.4% of every plasmid under 10 kb**
— median 4,665 bp, 6 CDS, 99.1% non-mobilizable, 97.0% with no typed replicon. AA379 is MOB-suite's
residual bin for small untypable plasmids: a category, not a clade, spanning environmental and
clinical ecologies. A dark family on 400 unrelated small cryptic plasmids scores **one** lineage. That
population is this project's own declared subject (`project-small-cryptic-plasmidome`). §32.6 checks
`MOB_count <= plasmid_count`, which cannot detect this. Nothing in §21's normalisation layer mentions
it. Selection later cannot repair it from the table — it would need per-plasmid ANI, a stage that
does not exist.

### 11. Context and recurrence labels have no null, so `defence_associated = TRUE` cannot be false.

§16.7 creates seven boolean association labels and §24 turns them into `functional_hypothesis`
values. §5.4/§16.7 correctly call them hypotheses. But §21.5 rules that "exact statistical tests
belong to a later analysis module", so in the delivered tables these fire on **any** adjacency, with
no background, no threshold, no confidence interval, and no multiple-testing control across what
will be 10⁵–10⁶ families × 7 features = **10⁶–10⁷ comparisons**. On a defence-carrying plasmid most
ORFs are within ±3 of a defence gene, so the label is near-universal there and absent elsewhere —
it measures plasmid content, not ORF biology. §24's worked example ("observed in 17 independent
plasmids") states a count with no denominator and no expectation. Add now: a permutation null
(shuffle gene labels within a plasmid, or draw size-matched ORFs) computed during the run, when the
neighbourhoods are in memory; it is unrecoverable afterwards without the re-run in finding 4.

### 12. §19's structural states carry a confidence field that this method cannot fill.

§19.1 requires `structure_confidence`; §19.3/§19.4 define `DARK_NO_STRUCTURE` vs `DARK_FOLD_KNOWN`.
`docs/PARAMETER_PROVENANCE.md` records `structure.min_plddt` as **removed** with the reason: ProstT5
does not predict a structure, it translates sequence directly into 3Di, so there is no pLDDT to
threshold — and while the parameter was declared-but-unread, "a family whose structure was too poor
to trust looked identical to a confident novel fold". **§19.1 reintroduces exactly that field.** It
will be `NA` for the whole set, and `DARK_NO_STRUCTURE` will mean "the search returned nothing",
which is also what a non-protein returns. §34.5 forbids calling it a novel fold, correctly — but the
open-discovery arm of the eventual library is then selected on a negative result with no confidence
axis, which is the class most enriched for the artefacts of finding 7. §19.1's "all dark proteins,
not only family representatives" also multiplies this stage ~4–5× for near-zero gain: two members of
a 30%-identity family return near-identical 3Di strings.

### 13. §31's configuration block contains the string `configurable` where the decisions belong.

`relaxed_evalue: configurable`; `family_thresholds: close/intermediate/broad: configurable`;
`dnds: conditional`. `docs/PARAMETER_PROVENANCE.md` exists precisely because this project has "already
been damaged twice this way — the unswept cascade thresholds in v1, and the invented scoring
coefficients in v2", and it counts 21 unreferenced parameters as a known weakness. The spec adds
more, un-numbered, and §30 Rule 1 ("do not invent biological thresholds not defined in this document
or configuration") then hands the choice to whichever agent implements the stage. §31 also
contradicts §7.3/§7.4 outright: the text excludes `protein_length < 20`, the config says
`min_dark_length_aa: 30`, and `config/config.yaml` says `min_orf_aa: 20`. FESNov's validated
antimicrobial peptide was **36 residues**; a 30-aa floor is not obviously safe and is nowhere
defended.

---

## The seven questions, answered directly

### Q1 — The deferral problem

**The evidence table is not sufficient.** Deferring selection is correct in principle (§1.2, §35);
the damage is that it was used to justify scoping measurements down. Three things selection will
need are not in the table and cannot be added without re-running a full stage:

1. **A background for every context statistic.** §16.1 scopes neighbourhood extraction to dark ORFs
   only; §21.5 defers all tests. Enrichment needs non-dark neighbourhoods on size-matched plasmids.
   **Re-run: the whole of Stage 8 over ~7–8M annotated ORFs.**
2. **Nucleotide sequence per ORF occurrence.** §6.3 stops at `protein_sequence`. RNAcode (§18.3) and
   dN/dS (§18.4) are both hedged to "where appropriate" / `conditional`, so coverage is partial and
   the gaps are not recoverable. **Re-run: CDS extraction over 143,504 plasmids, which additionally
   requires origin-aware coordinates the spec does not produce.**
3. **An independence unit.** §14 and §18.1 count database records and MOB clusters; finding 10 shows
   `mob_cluster` lumps 15.2% of the corpus. **Re-run: a plasmid-level ANI/mash clustering stage that
   does not exist in the spec at all.**

Two further gaps are cheap but must be decided now, not later: an opposite-strand-overlap column
(free from §6.3 coordinates) and RNAcode's antisense-strand result (§18.3 never mentions strand).
Both are single columns; both are the only evidence that separates a shadow ORF from a protein.

### Q2 — What the reviewer attacks

**(a) "How do you know these are proteins?"** — The central attack, and §7.1 makes it
**unanswerable**: it caps the artefact screen at AntiFam's 278 curated profiles and forbids adding
detectors. Nothing in the spec distinguishes a recurrent protein from a recurrently copied
mis-annotation, and §34.3 asserts the opposite without a mechanism. **Add now:** RNAcode on both
strands for every family with ≥3 CDS (`rnacode_max_p: 0.05`, FESNov, already in
`config/targets.yaml`), an opposite-strand-overlap fraction per ORF, and the family-consensus
re-check that Pavlopoulos used to remove 6.5% of clusters (implemented in this repo as
`consensus_recheck.py`, absent from the spec).

**(b) "What is your false-discovery rate?"** — **Unanswerable.** §32 and §39 contain no external
validation of any kind. **Add now:** the spiked reviewed-Swiss-Prot positive control at
`min_control_recall: 0.99`, which this repo already implements and halts on, plus the second
held-out arm the previous review demanded (proteins masked from the database being searched), because
a Swiss-Prot protein spiked into a Swiss-Prot search hits itself at 100% identity and the gate then
tests nothing.

**(c) "Your dark set is database sampling, not biology."** — **Answerable, but only if measured
now.** §21 lists the confounders and §21.5 defers every test that would address them; §14 actively
converts redundancy into evidence. **Add now:** plasmid-level ANI/mash de-duplication, cross-database
accession reconciliation across §3.2's four overlapping sources, and the rarefaction curve of
finding 8. All three are cheap during the run and unrecoverable after it.

### Q3 — Unfalsifiable output

Eight places, each already argued above: §13.3/§2.5 `DARK` (finding 1); §23.3
`independent_evidence_categories` (2); §13.5 `dark_evidence_level` (3); §14's
`hypothetical_deposit_count` / `independent_database_count` / `hypothetical_recurrence_level` (7 —
no null, and a shadow ORF's count rises with the gene beneath it); §16.7/§24's association labels
(11 — no background, no correction); §15.4's 100%-dark families (9 — arithmetic); §19.3/§19.4
`DARK_NO_STRUCTURE` (12 — no confidence axis); and §39 in its entirety (1 — a correctly implemented
but scientifically worthless run passes every bullet). §32.5 and §32.6, the only numeric checks in
the document, are inequalities over the code's own arithmetic and cannot fail.

### Q4 — The real bottleneck

Finding 6. §10.6 (HH-suite) is the single largest consumer — **9,715 core-hours at the most
optimistic per-query cost, 39,000–155,000 realistically, against a 16,128 core-hour budget** — with
§10.3 (InterProScan, 4,660–17,500) and an un-narrowed §10.7 (nr, unbenchmarked and stated to
dominate) close behind. Information yield is *inversely* proportional: §10.6 and §23.2 both declare
these outputs non-independent of Pfam and DIAMOND. **Cut:** move §10.6 and §10.3 after adjudication
so they run on the dark residue rather than all 3.5M, and restore the narrowing/reporting threshold
split with the 2% sweep cohort. **Spend on:** both-strand RNAcode over the dark set, neighbourhood
extraction over all 9.3M ORFs, and plasmid-level ANI. Those three turn three unanswerable referee
questions into answerable ones.

### Q5 — Garbage in

Finding 7, with the table. Headline: **30–50% of the reported dark ORFs are artefact or non-plasmid
in origin**, and the recurrent subset — the one any selection will prefer — is enriched for them,
not depleted, because shadow ORFs, chromosomal contamination, re-deposition and IS-derived ORFs are
all recurrent by construction. The two largest single contributions are quantified on this exact
dataset: **160,375 origin-spanning ORFs (1.72%)** that §6 cannot handle in either direction, and
**21,817 plasmids (15.2%) in one `mob_cluster`** that §14 and §18.1 will count as independent
occurrences.

### Q6 — The 1,000 number

Finding 8. No evidence in the spec; the number comes from the synthesis budget. The dataset's own
prior measurement puts the substantial dark family space at **~5,000 dark-only families with ≥10
members** and ~15,000 spanning ≥2 lineages, so 10³ is the right *order* but for reasons the spec
never states and never checks. It becomes visible in §15.3's family-size distribution and §26's
table; §26, §32, §38 and §39 all decline to require it. The rarefaction curve is the fix and costs
nothing.

### Q7 — Reintroduced defects

| previous review § | defect | where it returns |
|---|---|---|
| §1.1 | the run-halting positive-control recall gate | **deleted entirely** — §32, §39 have no external validation |
| §4.1 | evidence lines treated as independent when they are nested | §23.3, with no `IMPLIED_BY` entailment map |
| §4.2, §4.4 | shadow ORFs pass every absence-based test | §7.1 forbids new detectors; §18.3 RNAcode never specifies strand |
| §4.4 | AntiFam's 278 profiles as the sole artefact screen | §7.1, made explicit policy |
| §4.5 | dereplication before family definition inverts family measures | §8 + §15.3; `family_size` never defined over occurrences vs unique sequences |
| §5.4 | pseudo-replicated conservation denominator | §17.2; no clonal/independence unit anywhere |
| §5.6 | no multiple-testing control, no null | §21.5 defers all tests out of the pipeline |
| §5.7 | declared context classes never implemented | §16.4 still has no AMR-island or IS/transposon class |
| §6 | `mob_cluster` is not a lineage (AA379 = 15.2%) | §18.1, §22, §32.6 |
| §2.5 | "novel fold" means "the search failed" | §19.3/§19.4, plus §19.1 reintroducing `structure_confidence`, the field `PARAMETER_PROVENANCE` removed as unfillable |
| §8.1 (praised) | `narrow_at` vs `min_explained` + 2% sweep cohort | absent; §11's three threshold classes omit the narrowing threshold |
| §8.2 (praised) | `-Z`/`--domZ` pinned to 3,497,616 | absent; §11.2 promises re-thresholding on E-values that are not comparable without it |
| §8.2 (praised) | `--domE`/i-Evalue (7.2% of T2 resolutions depend on insignificant domains) | §10.1 records `pfam_domain_evalue`, no i-Evalue, no `--incdomE` |
| §8.6 (praised) | circular-origin repair, 160,375 ORFs | absent; §6 has no topology handling |
| §8.10 (praised) | pre-flight tool/database check | absent; §29's Validation Agent runs at the end, §32 "before producing final outputs" — the 45-hour DIAMOND-not-on-PATH failure, which occurred twice, is re-enabled |
| — | AntiFam `--cut_ga` (274/278 curated thresholds looser than 1e-5) | §7.2 records `antifam_evalue` only |

---

## What I would do before writing a line of code

1. Put back the run-halting control gate and add the held-out arm. Without it, §39 is not a
   definition of success.
2. Move §10.6 and §10.3 behind adjudication and restore the narrowing/reporting split. That is
   where the budget for everything below comes from.
3. Make three columns mandatory in §6.3/§25: nucleotide sequence, opposite-strand overlap fraction,
   and RNAcode on both strands. They are the only artefact evidence the spec can obtain.
4. Extract neighbourhoods for all 9.3M ORFs, not only dark ones, and compute a permutation null
   during the run.
5. Add a plasmid ANI/mash clustering stage and use it, not `mob_cluster`, as the independence unit
   in §14, §17 and §18.1.
6. Require one figure in §38: the rarefaction curve of dark families against plasmids sampled.
7. Delete `dark_evidence_level` (§13.5) or replace it with the counts it is a function of, and
   either drop `independent_evidence_categories` (§23.3) or ship the joint-indicator table beside it.

Every one of items 1, 2 and 3 already exists, working and tested, in `src/plasmidann/` and
`workflow/scripts/`. The spec's most serious problem is not what it gets wrong; it is what it
forgot this repository already knows.
