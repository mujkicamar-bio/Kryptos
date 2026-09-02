# Reassessment of Round-1 Review Findings

**Manuscript under revision:** `manuscript/amr_onehealth_NAR.md`
**Review being answered:** `manuscript/review_round1.md`
**Produced entirely by:** `scripts/reassess_review_findings.py`
**Date:** 2026-07-27

Regenerate:

```bash
export MAMBA_ROOT_PREFIX=$PWD/envs/.micromamba_root
ENV=$MAMBA_ROOT_PREFIX/envs/stats_env
PYTHONNOUSERSITE=1 PATH="$ENV/bin:$PATH" $ENV/bin/python scripts/reassess_review_findings.py
```

Outputs: 6 tables `data/plasmidscope_primary/reassess_*.tsv`, 1 figure
`reports/figures/onehealth/reassess_compartment_or_shift.png`.

**Bottom line: the reviewers were right on the substantive points, and the manuscript's central
claim does not survive.** Four of six testable objections are upheld, one is upheld with a
qualification, and one is empirically refuted in the paper's favour. The paper needs reframing, not
patching. Details below, in order of consequence.

---

## 1 · UPHELD — the title claim is refuted by its own data (DA1 / M3)

The manuscript argued "mobility rivals provenance" from the raw Table 3 contrast: conjugative
plasmids from natural environments carry resistance at 27.4%, above non-mobilizable clinical plasmids
at 23.9%. The reviewers objected that this compares cells differing enormously in plasmid size.

They were right, and the direction of the error is worse than suspected. The two cells differ 3.4-fold
in median length — **99 516 bp for environmental conjugative versus 29 120 bp for clinical
non-mobilizable**. Computing the same contrast as adjusted predicted probabilities at matched size and
median GC:

| Plasmid size | Env-natural, conjugative | Clinical, non-mobilizable |
|---|---:|---:|
| 5 kb | 4.1% | **7.4%** |
| 10 kb | 6.5% | **11.5%** |
| 50 kb | 17.3% | **28.2%** |
| 100 kb | 25.3% | **38.8%** |

**At every matched size the clinical non-mobilizable plasmid is the more likely carrier, by roughly
1.5–1.8-fold.** The raw ordering reversed only because environmental conjugative plasmids happen to
be 3.4× larger. The sentence the paper is named after is an artefact of size composition.

Mobility remains a genuine and large effect (conjugative OR 4.77 adjusted), so the broader point that
mobility matters is intact. But "mobility rivals provenance" as a ranking claim is not supported:
matched for size, provenance wins at every point tested. The title, abstract framing and Discussion
¶2 must change.

## 2 · UPHELD — roughly half the compartment effect is host lineage (D1 / DA3)

Host taxonomy was joined from PlasmidScope `Host` via the provenance key. Adding host genus (top 20 +
"other") to the model, on the subset with host data (*n* = 55 288), attenuates every major
compartment effect by close to half:

| Term | OR, compartment only | OR + Enterobacteriaceae flag | OR + host genus | Change |
|---|---:|---:|---:|---:|
| Human — clinical | 7.66 | 5.27 | **3.30** | −56.9% |
| Human — community | 6.03 | 4.13 | **3.15** | −47.8% |
| Animal — livestock & poultry | 5.22 | 3.38 | **2.85** | −45.5% |
| Environment — engineered interface | 2.99 | 2.16 | **1.61** | −46.2% |
| Animal — companion & aquaculture | 2.26 | 1.92 | 2.26 | +0.1% |
| Plant & food | 1.08 | 1.16 | 1.35 | +25.5% |
| Animal — wildlife & other | 1.19 | 0.96 | 1.17 | −1.8% |
| conjugative | 4.22 | 3.51 | 3.38 | −20.0% |
| mobilizable | 2.30 | 2.30 | 1.81 | −21.1% |

Pseudo-*R*² rises from 0.243 to 0.299 when host genus enters, confirming it carries substantial
independent information the original model omitted.

The compartment effects **survive** — clinical remains 3.3× the environmental baseline and highly
significant — but at roughly half the magnitude the manuscript reports. The claim "plasmids from
human clinical sources remain roughly nine times likelier to carry resistance" becomes approximately
**three times**, and must be restated.

**Important qualification on this correction.** Host coverage is severely uneven across compartments:
99.2% for clinical, 98.6% livestock, but only 25.6% for natural environment and 20.1% for wildlife.
The host-adjusted model therefore runs on a subset that over-represents exactly the compartments with
the strongest effects. The before/after comparison above is internally valid — both models are fitted
to the same rows — but generalising the ~50% attenuation to the full collection assumes a
representativeness the data do not establish. The honest statement is that host lineage accounts for
a substantial fraction of the apparent compartment effect, with the precise fraction uncertain.

## 3 · UPHELD — reported precision was overstated 2–9-fold (M1)

Refitting with standard errors clustered on MOB cluster (6182 clusters in the modelling set) widens
every confidence interval substantially:

| Term | OR | Nominal 95% CI | Clustered 95% CI | CI width ratio | Still *P* < 0.05 |
|---|---:|---|---|---:|---|
| Human — clinical | 8.94 | 8.27–9.66 | 6.34–12.61 | 4.5× | yes |
| Animal — livestock & poultry | 5.92 | 5.38–6.53 | 4.01–8.74 | 4.1× | yes |
| Human — community | 4.99 | 4.62–5.40 | 3.72–6.71 | 3.8× | yes |
| Environment — engineered interface | 2.94 | 2.67–3.23 | 2.22–3.88 | 3.0× | yes |
| Animal — companion & aquaculture | 2.51 | 2.14–2.94 | 1.69–3.71 | 2.5× | yes |
| conjugative | 4.77 | 4.54–5.02 | 3.30–6.90 | 7.4× | yes |
| mobilizable | 2.77 | 2.63–2.91 | 1.84–4.17 | 8.3× | yes |
| log₁₀ size | 4.89 | 4.71–5.09 | 3.54–6.76 | 8.6× | yes |
| GC content | 0.99 | 0.990–0.995 | 0.978–1.007 | 6.5× | **no** (*P* = 0.31) |
| Plant & food | 1.09 | 0.98–1.21 | 0.80–1.50 | 3.0× | no (consistent) |
| Animal — wildlife & other | 0.85 | 0.71–1.02 | 0.59–1.22 | 2.0× | no (consistent) |

All substantive conclusions survive, but **the GC effect disappears entirely** — it was an artefact of
treating clustered observations as independent. Every *P* value in the manuscript's Table 4 must be
replaced with the clustered version, and the *P* < 1 × 10⁻³⁰⁰ figures removed.

## 4 · UPHELD, reframed — size confounding is larger than collinearity (M2)

Formal collinearity is not the problem: `log_size` has VIF 1.62 and GC 1.11. The compartment dummies
show VIF 4.4–15.5, but that is largely intrinsic to an eight-level categorical and not itself
diagnostic.

The real finding is confounding. Removing `log_size` from the model changes the conjugative
coefficient from **OR 4.77 to OR 12.23** — meaning **about 61% of the raw conjugative effect is
plasmid size**, not transfer capability. (Mobilizable moves the other way, 2.02 → 2.77.) The
manuscript's description of mobility and size as independent axes is wrong; size is a major mediator
or confounder of the mobility effect and must be described as such.

## 5 · UPHELD — the multireplicon section is not salvageable in this paper (D3)

PlasmidFinder typing coverage varies almost nine-fold across the compartments being compared:

| Compartment | Typing coverage |
|---|---:|
| Human — clinical | 38.1% |
| Animal — livestock & poultry | 34.0% |
| Human — community | 26.2% |
| Environment — engineered interface | 16.7% |
| Plant & food | 15.7% |
| Animal — companion & aquaculture | 15.8% |
| Animal — wildlife & other | 6.8% |
| Environment — natural | **4.3%** |

The "fewer than two replicons" comparator is thus composed almost entirely of untyped plasmids in the
environmental compartments (95.7% untyped in natural environment) and much less so in the clinical
one. Cross-compartment comparison of multireplicon odds ratios compares groups whose reference
categories mean different things. **Recommend cutting the section**, as the domain reviewer advised.

## 6 · REFUTED — the 50% exclusion is not an AMR-biased exclusion (P2)

This is the one objection the data clear. Comparing resistance prevalence in assigned versus excluded
plasmids:

| Group | *n* | AMR prevalence |
|---|---:|---:|
| Compartment-assigned | 104 169 | **18.81%** |
| Excluded — no habitat label | 39 331 | **17.86%** |
| Excluded — simulated communities | 64 658 | 4.29% |
| Excluded — other | 90 | 2.22% |

The 39 331 plasmids excluded for lacking a habitat label carry resistance at essentially the same rate
as the analysed set — a 0.95-point difference. Habitat-label recovery is therefore close to
non-informative with respect to resistance, which substantially defuses the manuscript's own
"single largest threat" framing. (The simulated-community bucket differs sharply, but those are
excluded as sampling artefacts by a locked prior decision, not on any AMR-related ground.)

Limitations should be **revised downward** here, reporting this test rather than the current
unquantified worst-case statement.

## 7 · QUANTIFIED — the clinical-flag circularity is modest (D4)

Of the 18 171 plasmids in the clinical compartment:

| Basis for the clinical flag | *n* | Fraction |
|---|---:|---:|
| Clinical body site (blood/urine/respiratory/skin-wound/other clinical) | 12 614 | 69.4% |
| No body site, PLSDB disease tag only | 4 128 | 22.7% |
| No body site, no disease tag — context only | 1 429 | 7.9% |

Roughly 69% of clinical assignments rest on an independent body-site annotation, so the circularity
concern touches about a quarter of the compartment. Worth a sensitivity analysis restricted to
body-site-defined clinical plasmids, but not a threat to the compartment's validity.

---

## What the revised paper should claim

The reassessment does not destroy the study; it redirects it. Three results are robust to every
correction applied above:

1. **A genuine compartment gradient exists**, but at roughly half the reported magnitude — clinical
   about 3.3× the environmental baseline after host-taxon adjustment, not 8.9×.
2. **Mobility is a large independent effect** (conjugative OR ≈ 3.4 after host adjustment, ≈ 4.8
   before), but it does not out-rank provenance at matched size, and about 61% of its raw magnitude
   is plasmid length.
3. **Plasmid size is the dominant single predictor** throughout, and both the original mobility claim
   and two of the original compartment effects were size artefacts.

The defensible thesis is no longer "mobility rivals provenance." It is closer to: *apparent One Health
structure in plasmid resistomes is substantially explained by host lineage and plasmid size, and
roughly halves when either is controlled — but a real compartment gradient and a real mobility effect
both survive.* That is a more conservative claim than the current draft makes and, because it
quantifies how much of a widely assumed pattern is compositional, arguably a more useful one.

## Recommended disposition

Reframe and re-run rather than patch. Concretely: retitle; rebuild Table 4 as a model sequence
(compartment only → + size/mobility → + host taxon), all with clustered standard errors; replace the
Table 3 headline contrast with the matched-size predicted probabilities in §1; cut the multireplicon
section; revise the Limitations exclusion paragraph downward per §6 and add the host-coverage
imbalance from §2 as a new limitation. Items already verified as unaffected — the reproducibility
apparatus, the compartment mapping, the drug-class signatures and the Jaccard repertoire structure —
carry over unchanged.
