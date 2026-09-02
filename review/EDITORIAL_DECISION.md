# Editorial decision — internal rigour audit of the plasmid analysis programme

**Date:** 2026-09-02
**Scope:** whole project, both arcs (resource + AMR/One Health; small cryptic + dark plasmidome)
**Standard applied:** internal rigour ("is this true and reproducible"), not any venue's taste
**Panel:** 5 independent reviewers + an editor's workflow audit + an independent numerical recompute

| input | file |
|---|---|
| Editor-in-Chief | `review/reports/eic_review.md` |
| R1 — methodology & statistics | `review/reports/methodology_review.md` |
| R2 — plasmid biology | `review/reports/domain_review.md` |
| R3 — protein dark matter | `review/reports/perspective_review.md` |
| Devil's advocate | `review/reports/devils_advocate_review.md` |
| Editor's workflow audit | `review/reports/workflow_audit.md` |
| Independent recompute | `review/VERIFICATION_PACK.md` (+ `review/verify/*.py`) |
| Corroboration analysis | `review/CORROBORATION_MATRIX.md` |

Reviewers worked blind to one another. Nothing below is asserted that does not trace to a named
report; where seats disagree, the disagreement is shown and arbitrated rather than averaged.

---

## DECISION

**Major revision — of the claims, not of the work.**

The pipeline is sound. An independent recompute reproduced **47 of 52** headline quantities exactly,
and the two initially-flagged discrepancies that were mine have been retracted. Three of the five
seats tried hard to break the data layer and could not: all six master-table joins are strictly 1:1,
the dedup is byte-identical across all 147 duplicated plasmids, the CARD zero-fill is benign, and the
CDS-indexing fix holds with zero join loss.

**The two objections most likely to be fatal were both tested independently by two seats and both
failed.** This is the most important result of the review and it should be stated first:

- **"It's metagenomic assembly artefact."** DA F12 — *fails decisively*. Stratified by provenance the
  Pfam cross-tab is invariant (isolate-only 59.5% vs IMG-PR 56.8% in the "neither can name" cell);
  isolate-only families give **74.9% recurrence** and 663 Pfam-negative substantial families with
  *higher* lineage spread (80.5%). A phage-hallmark screen the project never ran finds **<1%** of the
  71,414 carry a structural hallmark.
- **"The dark ORFs are spurious gene calls."** PER Q1–Q7 and DA F11 — *fails*, from two directions.
  Against 66,013 uncalled ORFs from the same replicons, 83.0% of dark Pfam-negative ORFs score >0 on
  codon LLR vs 16.5% of uncalled; composition is near-identical to PlasAnn-named proteins (L1 0.099)
  and far from random-codon expectation (0.303); dark ORFs overlap neighbours *less* than named CDS.

So the dark proteome is real protein, on real plasmids. What does not survive is **how much of it is
novel, and what its recurrence means.**

---

## Adjudication of Devil's Advocate CRITICAL findings

*(Required: every DA CRITICAL is adjudicated visibly. Both are validated; both block the current
framing.)*

### DA F1 — the project's own clonal-redundancy control was never applied to the recurrence claim
**VALIDATED. This is the review's most consequential finding.**

The project uses one-plasmid-per-MOB-cluster dereplication in `cryptic_plasmids.ipynb` and
`defense_systems.ipynb` — the deck even calls it "the same MOB-cluster control used throughout Act
II". `grep` confirms it appears in **no dark-family notebook**. Applied (20 draws):

| scenario | ORFs | families | recurrence | ≥10 members | Pfam-negative |
|---|---:|---:|---:|---:|---:|
| pooled (published) | 378,552 | 92,752 | **84.7%** | 5,623 | **3,343** |
| cap 50 / cluster | 206,608 | 56,043 | 83.2% | 3,537 | 2,046 |
| cap 10 / cluster | 125,203 | 46,938 | **75.7%** | 2,014 | **1,040** |
| 1 per cluster | 32,763 | 23,514 | **39.2%** | 159 | **40** |
| random, ORF-matched | 32,553 | — | 61.3% | — | — |

The ORF-matched control is load-bearing: at equal ORF count, removing clonal structure costs **22
points**, so this is not sample size. The 93-family core drops to ~30 families still on ≥10 carriers,
with a **median of 7 dereplicated carriers** — families selected for being "on ≥100 distinct
plasmids" were measuring deposits.

**The DA's own fair counter is accepted:** one-per-cluster is too harsh (AA379 alone absorbs 26% of
the subset, and DA F5 shows it is a junk drawer spanning 236–19,992 bp, not a lineage). **The
defensible statement is the cap-K middle: ~75.7% recurrence and ~1,040 Pfam-negative substantial
families at a cap of 10 per cluster.** Not 84.7% / 3,343, and not 39.2% / 40.

### DA F2 — recurrence is a function of sampling depth, and the null was never computed
**VALIDATED, and independently corroborated by a second seat.**

Rarefying carriers: 39.2% at 1,000 plasmids → 68.6% at 10,000 → 84.7% at 63,993. "84.7% recurrent,
therefore recurrent-families-world rather than idiosyncratic-junk-world" is therefore not a
discriminating test — it reads where sampling stopped. **PER F3, working blind, built the null the DA
says is missing** and found core families sit at median obs/exp **0.18** against random expectation;
calibrated against Pfam-named families (0.312 vs 0.197, p=2×10⁻⁷²), the candidate-novel families
travel **less** than the known backbone.

Two seats, two methods, one conclusion: **the recurrence argument as written does not support the
functional inference drawn from it.** The sentence "conservation across 224 lineages is not the
signature of neutral junk" cannot stand.

---

## Consensus findings (3+ seats, reached independently)

### CF1 — "dark = novel" is inflated at every link (4 seats, 4 routes)
DOM F1 (`phmmer` vs the project's *own* named proteins: 22.8% of dark-only ORFs have a named homolog
at E<1e-5, best 2.6e-204 → dark-only **85.9% → ~66%**) · PER F1 (`--cut_ga` → `-E 1e-3` on the same
Pfam-A: **948/3,343, 28.4%** of the deliverable and **54.8%** of the core recover a name) · DA F4
(2,000 residue families vs a composition-shuffled null: 4% at Bonferroni to **27%** at E≤1e-3) ·
MET C4 (scope regex recall failure: "47.3% inside PlasAnn's categories" → **≥54.8%**).

Every route recovers the same domains: **Rep_1, RepC, Relaxase, Mob_Pre, NikA-like, HTH, TrfA, Rop.**
The plasmid backbone, unlabelled. This is the project's own §4.3 insight recurring one level deeper
than the project took it.

**And `grep` confirms Pfam-A is the only external database ever searched** (DA F4) — no DIAMOND or
BLAST against UniRef/NR anywhere in the project. The novelty claim is bounded by one database at a
family-membership threshold that was never meant to be a homology detector.

### CF2 — a systematic "count labels, not entities" defect across four resource columns (3 seats)
`pf_n_inc` counts alleles while `pf_inc_families` is computed and discarded (pack §7.2, EIC F2,
MET C1) — and is documented backwards in `master_table_data_dictionary.md:129` as "count of distinct
Inc **families**… the paper's definition". `card_multidrug` unions drug-class label strings, so
**4,402 (21.2%) of the 20,741 "multidrug" plasmids carry exactly one ARO** (MET C2). `plsdb_n_inc`
inherits the same (MET C3). The locked mob_typer fallback does not escape it: 23.0% of its 18,685
multireplicon calls are IncF-only and **60.4% contain an uncalibrated `rep_cluster_NNNN` token**
(DOM F4).

Compounding it, MET D1: the PlasmidFinder rule is documented as "identical criteria to CGE" but CGE's
default is 95% identity, not 80%. **Multireplicon spans 8,169–24,583 across two undocumented
choices**, with no sensitivity published — under the question the project was built to answer.

### CF3 — corrections are made, then not propagated (3 seats)
EIC F1/F3/F8, PER F2, WF W1/W2/W4. Sharpest instance is **PER F2**: the flagship §R2 table in
`dark_orf_clustering.md` cites five family representatives that **do not exist in the shipped
clustering** — including `GenBank_CP048556.1|10`, the source of the "224 lineages" headline — in a
report whose text states "every number above is post-fix". Meanwhile
`amr_onehealth_methodology.md:248` still asserts "Mobility rivals provenance", and three documents
still request a threshold sweep that was run on 2026-09-01.

### CF4 — the dispersal paradox is circular and its headline runs the wrong way (2 seats, 4 findings)
DOM F3: `mob_orit` is populated for **zero** of 79,388 non-mobilizable plasmids *by construction* —
mob_typer calls them non-mobilizable *because* it found no oriT — while PlasAnn independently flags an
oriT on **3,070 (3.9%)**, a floor given PlasAnn's 38% oriT sensitivity. `orit_db_folder/` has been
indexed since 10 July and never used. MET D2: `Unlabelled` and `Geography only` are counted as
habitats, and the conjugative arm carries 1.8× more of them. MET D3 + DA F7: the K-ladder the report
itself publishes declines monotonically (0.929→0.800), so "near-parity" is contradicted by the
project's own table; the isolate-only stratum is the *confounder-controlled* analysis and reporting it
as a caveat under the pooled result inverts the evidential hierarchy. MET D4: the rank-biserial effect
size is published with **inverted sign**.

Note the direction: **MET D2 and D3 are corrections against the project's own interest that the
project did not make.**

---

## Disagreements between seats, arbitrated

**D-1 · What should lead the paper?**
EIC §3 says lead with the catalogue and explicitly *not* the shortlist. DA F14 says the opposite — the
five UniRef90 survivors, which "withstood everything I threw at them" (14–42 carriers *after* full
clonal dereplication, dominant-lineage share 10–30%), should occupy the abstract.

**Arbitration: both are arguing about which count leads, and the count is the wrong headline.** The
EIC's recommendation rests on its F4 result that the catalogue is robust to *re-clustering*
(Jaccard 0.928) — a real and valuable finding, but on an axis nobody attacked. DA F1 and PER F3/F4
attack a different axis — clonal and sampling independence — that the EIC could not see, and on that
axis 3,343 does not hold. Conversely a five-family paper is too thin to carry the work.

**What survives all five seats is the decomposition with its positive control**, and DOM F7 names the
right frame: this is not a novelty claim, it is **the first calibrated measurement of the plasmid
annotation gap**. To a plasmid biologist "most genes on small cryptic plasmids are unannotated" is the
field's standing complaint; what is new here is the *magnitude at scale*, the *positive control* that
makes the number interpretable (NAMED 76.7% vs DARKREP 26.3%), and the *partition* into ~29%
annotation failure vs a residue. Framed as novelty it draws "we knew that"; framed as calibration it
is publishable. The catalogue then appears at its defensible size, and the five survivors appear as
the worked shortlist.

**D-2 · Is the catalogue robust?**
EIC F4 (robust, Jaccard 0.928 under `--cluster-reassign`) vs PER F4 (−22% of families under FESNov's
coverage rule; clusterings not nested) and DA F1 (3,343 → 1,040 at cap-10 dereplication).
**Arbitration: not in conflict — different axes, and both must be reported together.** The catalogue
is stable to *how you cluster* and unstable to *how you define novelty and how you handle clonality*.
Saying only the first, as the project currently does implicitly, is the more misleading half.

**D-3 · How well calibrated are the claims?**
DOM scores interpretive calibration **73** (the project's best axis) while EIC scores claim–evidence
calibration **55** and MET scores methods-prose fidelity **58**. **Arbitration: not in conflict —
the corpus is bimodal.** The newest documents reason about evidence exceptionally well; the older
documents and the code comments describe it inaccurately. EIC's formulation: *"the programme doubts
the wrong things."*

---

## Scope challenge that must be answered

DA's "So what?" verdict: **large non-conjugative plasmids are 73.3% dark and hold 3.67M dark ORFs.
The studied set is 6.8% of the dataset's own dark proteome.** Until a matched comparison shows that
small-plasmid dark families are *distinctive*, the compartment choice is convenience, not a finding.
This is answerable with data already on disk (`data/fam_3300022589/big_dark_orfs.faa`, 5,135,708
proteins, already extracted) and it is T7 of the project's own reclustering plan.

---

## What is genuinely strong (unanimous across seats)

1. **The Pfam positive control.** All four substantive seats independently said the same thing: most
   published dark-matter surveys omit it, and it is the only reason "26.3%" means anything.
2. **Publishing the CDS-indexing bug** — a defect that *inflated the project's own novelty claim* —
   with the corrupted cross-tab printed beside the corrected one. EIC: "the strongest single
   credibility signal in this repository."
3. **`reassessment_round1.md`** — killing your own title claim with matched-size predicted
   probabilities, and refuting one reviewer objection *in the paper's favour*.
4. **`clustering_threshold_sensitivity.md`** — DOM: "the best document in the repository". PER: the
   self-initiated 89.12% compliance disclosure is "rare and admirable".
5. **The dbAPIS adjudication** catching `AcrIIA21` firing on RepA_N replication initiators, and
   recognising that concordance of *artefacts* between two tools is itself evidence.
6. **Definitional discipline** where it was hardest: the habitat taxonomy locked before analysis, the
   One Health compartments built from `hab_sub`+`is_clinical` by importing the map so notebooks
   cannot drift, the simulated-community exclusion filtered on `hab_top` for a stated reason.

PER's summary judgement is the fairest one-line verdict on the programme: *"the instincts of a careful
lab using out-of-date tools."*

---

## Revision roadmap, ranked

Items 1–3 are **editing, not science** (~2–3 days) and nothing should be written up before they are
done. Items 4–6 are the analysis that decides how large the paper is.

**1 · Reconcile the claim record, and put the repo under version control.** (EIC F1/F3/F6/F8, PER F2,
WF W1–W4) Supersession headers on `amr_onehealth_methodology.md` and `amr_onehealth_NAR.md`; forward
pointers from `dark_orf_clustering.md`, `pfam_dark_validation.md` and `dark_plasmidome.ipynb` to
`clustering_threshold_sensitivity.md`; **regenerate the §R2 table, whose five representatives do not
exist in the current output**; rebuild `PROJECT_OVERVIEW.md` to cover both arcs with an explicit
claim-status table (live / superseded / withdrawn, and by what). `git init` in the same pass, plus
checksums on the frozen `dark30_*`/`mix30_*` files so the "non-negotiable constraint" becomes
enforceable. Execute and save `notebooks_my/*.ipynb` **with outputs** — every Arc-2 number currently
exists only as hand-typed prose. Export the `genesis` and `panaroo` env specs.

**2 · Fix the resource layer.** (CF2) Correct `master_table_data_dictionary.md:129`; add `pf_n_fam` to
`harvest_typing.py` (already computed, thrown away) and rebuild; publish the multireplicon count as a
**range across the identity threshold (8,169–24,583) with the 95%-CGE value as primary**, or drop
PlasmidFinder multireplicon in favour of a cleaned mob-based figure that excludes `rep_cluster_NNNN`
tokens; fix `card_multidrug` to count distinct AROs; widen the `pfam_verdict.py` scope regex.

**3 · Re-derive every dark-plasmidome headline under clonal control.** (DA F1/F2, PER F3) Adopt
**cap-10-per-MOB-cluster** as the standard, report the pooled and one-per-cluster values beside it,
and restate recurrence as **~76% (capped)** rather than 84.7%. Recompute the catalogue at that
control (~1,040 substantial Pfam-negative families). Add PER's sampling null and report obs/exp for
every family in the deliverable. Retire "spans ≥2 lineages" for DA F3's median-standard version
(no lineage >50%: **24.6%**, not 72.9%).

**4 · Finish the homology cascade before writing.** (CF1; PER §X1, Tier 0–1) UniRef90/NR by DIAMOND —
never run, and no database is on disk; a sub-GA Pfam pass with a composition-shuffled null; eggNOG.
PER sizes Tiers 0–1 at ~1 day of CPU. This decides the size of the residue, and every seat expects it
to shrink further.

**5 · Break the dispersal circularity.** (DOM F3) Run the relaxed oriT search against
`orit_db_folder/`, indexed since July and named by the project itself as the first follow-up. Then
restate §3 per MET D3: not "near-parity" but a modest, monotone 8–20% breadth deficit that widens with
sampling depth, led by the isolate-only stratum.

**6 · Name the ColE1 cassette, and answer the scope challenge.** (DOM's dedicated section; DA's
"So what?") ORF-B's best in-corpus homolog is `traD` at E=1.7×10⁻³⁶ and ORF-A's is a mislabelled
`sugE` at E=9.5×10⁻⁷⁶ (141 identical residues); with gene order and the invariant 177 aa these map
onto ColE1 **MbeB/MbeC**, characterised since Boyd, Archer & Sherratt 1989. **The named result is
stronger than the anonymous one**: two of the four genes of the canonical ColE1 mob operon are
unnamed by both PlasAnn and Pfam across 336 plasmids. Separately, run the matched comparison against
the large-plasmid dark proteome already extracted.

**Deferred, correctly:** the structure test. PER F7 establishes that FlashFold's absence is not a
scheduling choice — `data/flashfold_run/install.log` is a **disk-quota failure**, and foldseek,
hhblits and colabfold are installed nowhere on the system. Fix the infrastructure or drop the claim
that structure is the next step; do not keep deferring to a test that cannot currently be run.

---

## Consolidated scores

| dimension | EIC | MET | DOM | PER | panel range |
|---|---:|---:|---:|---:|---|
| internal / statistical / biological validity | 72 | 66 | 56 | — | 56–72 |
| external validity | 58 | — | — | — | 58 |
| reproducibility | 68 | 65 | — | — | 65–68 |
| claim–evidence calibration | 55 | 58 | 73 | — | 55–73 |
| contribution | 74 | — | 61 | — | 61–74 |
| documentation / literature / usability | 63 | — | 34 | 25 | 25–63 |
| novelty-claim validity | — | — | — | 35 | 35 |
| code correctness | — | 72 | — | — | 72 |
| methods currency | — | — | — | 40 | 40 |

**Lowest scores cluster on one axis:** novelty-claim validity 35, catalogue usability 25, literature
engagement 34, methods currency 40. Those four are the same problem seen from four seats — the
science is measured well and then described in the vocabulary of a field whose current standards it
does not meet.

**Disposition.** Arc 1 is closed; declare it so. Arc 2 is a real paper that is currently claiming the
wrong thing at the wrong size. The measurement is sound, the artefact objections failed, and the
deliverable exists — but the headline must move from *novelty* to *calibration*, the counts must be
re-derived under clonal control, and the homology cascade must be finished before anything is written.
