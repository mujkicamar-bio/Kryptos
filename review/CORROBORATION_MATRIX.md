# Cross-reviewer corroboration matrix

*Built by the synthesizing editor, 2026-09-02, from the five filed reports. The five reviewers worked
independently and could not see each other's findings; the workflow audit was produced by the editor.
Where several seats reached the same conclusion by different routes, that is corroboration, not
duplication — and it is the strongest signal in this review.*

Seats: **EIC** · **MET** (methodology) · **DOM** (domain) · **PER** (perspective) · **WF** (workflow
audit) · **DA** (devil's advocate — pending).

---

## T1 · The "dark = novel" chain is inflated at every link — 4 seats, 4 independent routes

This is the review's central convergent finding. Nobody was told to look for it; four seats found it.

| seat | route taken | measured deflation |
|---|---|---|
| **DOM** F1 | `phmmer` of 500 sampled dark-only ORFs against the project's *own* 78,397 named proteins | 22.8% have a named homolog at E<1e-5 (best 2.6e-204); dark-only **85.9% → ~66%** |
| **PER** F1 | same Pfam-A, same 3,343 representatives, `--cut_ga` → `-E 1e-3`; plus FESNov's own clustering criterion | **948/3,343 (28.4%)** recover a name; **51/93 (54.8%)** of the core |
| **MET** C4 | audited `pfam_verdict.py`'s scope regex for recall (`RelE` missing while `RelB` present; `TrfA`, `RepC`, `TrwC`, `Rop` filed novel) | "47.3% inside PlasAnn's scope" → **≥54.8%** |
| **PER** F6 | representative-only search used to label all 378,552 ORFs; partial-domain rescues counted as full | unquantified, additive |

**What every route recovers is the same thing:** `repB/repA/repL/mobV/mobC/mobA/traG` (DOM),
`HTH_17/NikA-like/RHH_1/Rep_1/Relaxase/PhdYeFM` (PER). **The plasmid backbone, unlabelled.** This is
the project's own §4.3 finding — "the artifact half is worse than a coverage gap" — recurring one
level deeper than the project took it. The programme found this pattern once and did not re-apply it
to its own residue.

**Editor's note:** these deflations are *not* additive in any simple way (DOM's is on the mix30
dark-only fraction, PER's on the post-Pfam deliverable, MET's on the scope classification). They must
be measured jointly, not summed. But every one of them moves the same direction.

## T2 · "Count the labels, not the entities" — a systematic bug class, 3 seats, 4 columns

| column | seat | defect | magnitude |
|---|---|---|---|
| `pf_n_inc` | pack §7.2 · **EIC** F2 · **MET** C1 | counts allele names; `pf_inc_families` is computed and thrown away | multireplicon 24,583 vs 22,422 family-level; one row reads 25 for 7 families |
| `pf_n_inc` (threshold) | **MET** D1 | rule documented as "identical criteria to CGE"; CGE's default is 95% id, not 80% | multireplicon **24,583 → 8,574 (−65%)**; full span 8,169–24,583 |
| `card_multidrug` | **MET** C2 | unions CARD drug-class *label strings* | **4,402 (21.2%) of the 20,741 "multidrug" plasmids carry exactly one ARO** |
| `plsdb_n_inc` | **MET** C3 | same allele inflation, inherited from PLSDB | — |
| `mob_rep_types` | **DOM** F4 | the locked fallback has its own version: 23.0% of the 18,685 are IncF-only; **60.4%** contain ≥1 uncalibrated `rep_cluster_NNNN` token | — |

**This is the most consequential resource-layer finding in the review.** It sits directly under the
question the project was built to answer, it is documented backwards in
`master_table_data_dictionary.md:129` ("count of distinct Inc **families**… the paper's definition"),
and the standing "mob_typer only" rule does not escape it. Four columns, one root cause.

## T3 · Corrections are made, then not propagated — 3 seats

| seat | finding |
|---|---|
| **EIC** F1 | `amr_onehealth_methodology.md:248` still asserts "Mobility rivals provenance"; zero forward pointers to the reassessment that killed it |
| **EIC** F3 | `dark_orf_clustering.md`, `pfam_dark_validation.md` and `dark_plasmidome.ipynb` (touched Sep 2, a day *after* the sweep) all still say the threshold is "untested at 40% and 50%"; none mentions the 89.12% compliance failure |
| **EIC** F8 | the refuted manuscript carries no warning in its own file |
| **PER** F2 | the flagship §R2 table cites **five representatives that do not exist in the shipped clustering** — including `GenBank_CP048556.1\|10`, source of the "224 lineages" headline — in a report stating "every number above is post-fix" |
| **WF** W1/W2/W4 | no version control, so no claim about which document was written against which data state is checkable; front-door doc blind to half the programme |

**PER F2 is the sharpest instance**: the single most-quoted number in Arc 2 is attached to a
family representative that the current output does not contain.

## T4 · The recurrence/lineage-spread argument does not survive its null — 2 seats, converging

| seat | finding |
|---|---|
| **PER** F3 | built the sampling null nobody had built: expected distinct `mob_cluster`s for a random draw of the same carrier count. Core families median obs/exp **0.18**, deliverable 0.20, top family 0.10. **Not one core family exceeds random expectation.** Against Pfam-named families (0.312 vs 0.197, p=2e-72) the candidate-novel families travel **less** than the known backbone |
| **MET** D7 | the compliance audit samples families uniformly when the claim is about ORFs (real in mechanism: compliance falls to 75% in families ≥50 — but moves the headline only 89.0→88.6%, so that conclusion **stands**) |
| **DOM** F6 · **MET** D8 | the host-range inversion is measuring the predictor: neighbourhood-based call, outcome-correlated selection |

**This directly contradicts a sentence in the EIC's own proposed abstract** ("63.9% sit in families
spanning two or more independent MOB lineages, so the recurrence is not clonal redundancy"). The EIC
could not see PER's null. Adjudicated in the decision below.

## T5 · The dispersal paradox is circular — 2 seats

| seat | finding |
|---|---|
| **DOM** F3 | `mob_orit` is populated for **zero** of 79,388 non-mobilizable plasmids *by construction* — mob_typer calls them non-mobilizable *because* it found no oriT. PlasAnn independently flags an oriT on **3,070 (3.9%)**, a floor (PlasAnn oriT sensitivity 38%). `orit_db_folder/` indexed since 10 July, unused |
| **MET** D2 | the rarefaction counts `Unlabelled` and `Geography only (no habitat)` as habitats; conjugative arm carries 1.8× more. Ratio 0.917 → **0.853**, p 4.9e-4 → 6.8e-5 — *the bug was weakening the project's own reported difference* |
| **MET** D3 | the K-ladder the report itself publishes declines monotonically (0.929→0.800 habitats, 0.878→0.793 countries). "Near-parity" is contradicted by the project's own table |
| **MET** D4 | rank-biserial published with **inverted sign** (+0.134 / +0.246 where Kerby gives −0.134 / −0.241) |

Note the direction: **MET D2 and D3 make the claim weaker; both are corrections against the
project's interest that the project did not make.** DOM F3 removes the paradox's foundation entirely.

## T6 · The catalogue is more robust than claimed, and less usable than needed — 2 seats, opposite signs

| seat | finding | direction |
|---|---|---|
| **EIC** F4 | re-derived the 3,343 catalogue under `--cluster-reassign` — nobody had: **3,490 families / 112,413 ORFs**, Jaccard 0.928, 94.8% ORF retention; 2,878 at UniRef90 | **under-claimed** |
| **PER** F5 | no sequences, no HMMs, no MSAs, no accessions, no confidence, drops the length column; representatives keyed on ids PER F2 shows are unstable | unusable |
| **PER** F4 | coverage is the dominant axis and was not swept: 92,752 → 72,335 families (−22%) under FESNov's rule; the clusterings are **not nested** | over-claimed |

**The catalogue is stable to re-clustering and unstable to how novelty is defined.** Those are
different axes and the project has conflated them.

## T7 · Things reviewers tried to break and could not — reported because they cost budget

| seat | test | result |
|---|---|---|
| **MET** C11 | can the master build multiply rows? | **No** — all six joins strictly 1:1 |
| **MET** C9 | is the CARD NaN→0 fill a real-zero fill? | Benign — RGI covered 208,245/208,248 |
| **MET** C7 | is the dedup safe? | Byte-identical across all 147 duplicated plasmids |
| **MET** — | artefact leakage into the One Health set | Zero; `mob_mobility` missing for exactly one plasmid |
| **MET** C12/D6 | is the Pfam propagation denominator wrong (editor's suspicion)? | **Editor was wrong** — 3,678 is correct, κ=0.909 |
| **MET** — | is the "95,442" figure irreproducible (editor's pack §7.3)? | **Editor was wrong** — `pf\|mob\|plsdb` = 95,442 exactly; pack retracted |
| **PER** Q1–Q7 | are dark ORFs gene-calling artefacts? | **Largely cleared.** Composition retains selection signatures; dark ORFs overlap neighbours *less* than named CDS (4.15% vs 7.05%); against 66,013 uncalled ORFs from the same replicons, 83.0% of dark Pfam-negative score >0 on codon LLR vs 16.5%. Residual risk bounded at **14–22 pp** |
| **EIC** F4 | does the 93→5 collapse nest as asserted? | **Yes** — all five 90% survivors map into the 93 |

PER's Q1–Q7 is the most valuable negative result in the review: it removes the objection that would
otherwise have been fatal, and no reviewer was going to accept the claim without it.

## T8 · Unanimous on the strengths

Every seat, unprompted, credited the same four things:

1. **The Pfam positive control** (NAMED 76.7% vs DARKREP 26.3%) — EIC, DOM, PER, MET all noted that
   most published dark-matter surveys omit it and that it is what makes the 26.3% interpretable.
2. **Publishing the CDS-indexing bug** — a defect that *inflated the project's own novelty claim*,
   written up with the corrupted cross-tab printed beside the corrected one.
3. **`manuscript/reassessment_round1.md`** — killing your own title claim with the right statistic.
4. **`clustering_threshold_sensitivity.md`** — DOM calls it "the best document in the repository";
   PER credits the self-initiated 89.12% compliance disclosure as "rare and admirable".

DOM adds the dbAPIS adjudication (catching `AcrIIA21` firing on RepA_N replication initiators) and
the defence-carriage-falls-with-lineage-count result. PER adds that the project "found and published
its own previous biggest deflation, which is why I expect it to absorb this one."

---

## Score consolidation

| dimension | EIC | MET | DOM | PER | spread |
|---|---:|---:|---:|---:|---|
| internal / statistical / biological validity | 72 | 66 | 56 | — | 56–72 |
| external validity | 58 | — | — | — | — |
| reproducibility | 68 | 65 | — | — | 65–68 |
| claim–evidence calibration | 55 | 58 (prose fidelity) | 73 (interpretive) | — | **55–73, the widest spread** |
| contribution | 74 | — | 61 | — | 61–74 |
| documentation | 63 | — | 34 (literature) | 25 (catalogue usability) | 25–63 |
| novelty-claim validity | — | — | — | 35 | — |
| methods currency | — | 72 (code) | — | 40 | — |

**The spread on calibration is itself a finding.** DOM scores *interpretive* calibration 73 — the
project's best axis — while EIC scores *claim–evidence* calibration 55 and MET scores *methods-prose
fidelity* 58. These are not in conflict: the project reasons about its evidence exceptionally well in
the newest documents and describes it inaccurately in the older ones and in the code comments. EIC's
phrase for it: **"the programme doubts the wrong things."**
