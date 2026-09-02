# Peer Review — Round 1

**Manuscript:** *Mobility rivals provenance: plasmid-borne antimicrobial resistance across One
Health compartments in 104 169 complete plasmids* (`manuscript/amr_onehealth_NAR.md`)
**Target venue:** Nucleic Acids Research
**Review date:** 2026-07-27
**Mode:** `full` — 5-seat panel (EIC + 3 peer reviewers + Devil's Advocate) + editorial synthesis

*This document is review output only. The manuscript was not modified. Objections marked
**[verified]** were empirically tested against the source data during review; the test and its result
are stated inline so the authors can reproduce or contest them.*

---

## PANEL CONFIGURATION (Phase 0)

| Seat | Identity | Angle |
|---|---|---|
| EIC | NAR database/genomics section editor; comparative microbial genomics | Fit, novelty, significance |
| R1 | Statistical genomicist; large-N observational inference | Design, model specification, inference validity |
| R2 | AMR domain specialist; resistome databases and host ecology | Domain accuracy, detection bias, literature |
| R3 | One Health / AMR surveillance scientist | Cross-domain framing, practical implications |
| DA | Devil's Advocate | Core-claim challenge, logical structure |

Primary discipline: microbial comparative genomics. Paradigm: observational, database-derived,
cross-sectional. Maturity: complete analysis, well-documented, single-round draft.

---

## REVIEWER 1 — METHODOLOGY

**Recommendation: Major Revision.** The analysis is competent, reproducible and unusually well
documented. Two specification problems, both verifiable from the authors' own data, currently
undermine the inferential claims.

### M1 [CRITICAL] — Observations are treated as independent when they are severely clustered

The logistic regression (Table 4) and the χ² test are fitted to 104 168 plasmids as if independent.
They are not. **[verified]** MOB-suite resolves only 7054 clusters across 208 248 plasmids; the
single largest cluster contains **39 837 plasmids — 19% of the entire collection** — and the ten
largest clusters together hold 30.7%. Median cluster size is 6.

The effective sample size is therefore smaller than the nominal *n* by orders of magnitude in parts
of the distribution, and every standard error, confidence interval and *P* value in Tables 1 and 4 is
correspondingly overstated. Reporting *P* < 1 × 10⁻³⁰⁰ against this dependence structure is not
meaningful. The §"Clonal redundancy" analysis shows the authors are aware of the redundancy, but it
is applied only to descriptive prevalence — never to the inferential model that carries the paper's
central claims.

*Fix:* refit with cluster-robust standard errors grouped on `mob_cluster`, or as a mixed-effects
model with a cluster random intercept. Report both nominal and clustered intervals. The point
estimates will likely survive; the precision claims will not.

### M2 [MAJOR] — Mobility and plasmid size are strongly collinear, and this is not diagnosed

Table 4 presents mobility (conjugative OR 4.77) and log₁₀ size (OR 4.89) as separable contributions.
**[verified]** They are substantially entangled: median length is **94 568 bp for conjugative
plasmids versus 5968 bp (mobilizable) and 5744 bp (non-mobilizable)** — a ~16-fold difference — with
Spearman ρ = 0.484 between log₁₀ size and the conjugative indicator.

At ρ ≈ 0.48 the model is estimable and the coefficients are not meaningless, but no collinearity
diagnostic is reported, and the manuscript's repeated framing of the two as "independent" axes
(Abstract; §"Adjustment separates…"; Discussion ¶2) overstates what the design can separate. A reader
cannot currently tell how much of the mobility effect is conjugative plasmids simply being large.

*Fix:* report VIF or condition indices; show the mobility coefficient with and without `log_size` in
the model; and soften "independent" to "separable but correlated" wherever it appears.

### M3 [MAJOR] — The headline comparison is unadjusted

The sentence the title rests on — a conjugative environmental plasmid (27.4%) out-carrying a
non-mobilizable clinical plasmid (23.9%) — is read directly off Table 3, which controls for nothing.
Given M2, those two cells differ ~16-fold in median plasmid length. The comparison is real as a
description of the raw data but cannot bear the interpretive weight placed on it in the Abstract and
Discussion.

*Fix:* produce the same contrast as adjusted predicted probabilities from the Table 4 model, at a
fixed plasmid size. If the ordering survives, the claim is earned and much stronger; if it does not,
the framing must change.

### M4 [MINOR] — "No Simpson's paradox" is imprecise

Monotonicity of the mobility ordering within every compartment is not what Simpson's paradox means;
the paradox concerns reversal between aggregated and stratified associations. The observation in
§"The mobility ordering holds inside every compartment" is worth keeping — just describe it as
consistency of the stratified association, not absence of Simpson's paradox.

### M5 [MINOR] — Multiple-testing families are not declared

Benjamini–Hochberg is applied within the 7 pairwise contrasts and within the 8 multireplicon
contrasts separately. State this explicitly, or correct family-wide.

**Scores** — Design 6/10 · Statistical validity 4/10 · Reproducibility 9/10 · Transparency 9/10.

*Reproducibility is a genuine strength and should be said plainly: a single seeded script regenerates
every number, tool versions are pinned, and the compartment mapping is fully enumerated. This is
above the norm for the venue.*

---

## REVIEWER 2 — DOMAIN

**Recommendation: Major Revision.** The compartment resolution is a real contribution. Two
domain-specific confounders, either of which could generate the central gradient on its own, are
currently unaddressed.

### D1 [CRITICAL] — The compartment gradient is confounded with host taxonomy

Resistance carriage is a property of bacterial lineage at least as much as of sampling environment,
and the compartments here are taxonomically very different. **[verified]** Among plasmids with a
species assignment, the human clinical compartment is **32% *Klebsiella pneumoniae* + 20%
*Escherichia coli*** (>half from two species), and livestock is **48% *E. coli* + 19% *Salmonella
enterica*** (two-thirds from two species). The natural-environment compartment has no species above
**5%** — it is taxonomically diffuse.

The manuscript therefore cannot distinguish "plasmids from clinical settings carry more resistance"
from "plasmids from Enterobacteriaceae carry more resistance, and clinical sampling is
overwhelmingly Enterobacteriaceae." This is the most consequential gap in the paper, and it is not
mentioned in the Results, Discussion or Limitations. The master table already carries
`plsdb_species` and `mob_host_range`, so the analysis is feasible with data in hand.

*Fix:* add host taxon (at least an Enterobacteriaceae indicator, ideally family level) to the Table 4
model, and report the compartment coefficients before and after. If the compartment effect attenuates
substantially, that is itself the paper's most interesting result and should be reported as such
rather than avoided.

### D2 [MAJOR] — CARD's clinical curation bias is not acknowledged

CARD is curated predominantly from clinically observed resistance determinants; environmental and
soil resistome diversity is systematically under-represented in its reference models. Detection
sensitivity is therefore not constant across the compartments being compared — it is highest in
exactly the compartment reported to have the most resistance. This biases the paper's central
gradient in its own favour.

The Limitations section covers conservative *thresholds* (Perfect+Strict) thoroughly, but threshold
conservatism is a different problem from database composition bias, and only the former is discussed.

*Fix:* add this explicitly to Limitations. Consider a sensitivity analysis restricted to ARO terms
with documented environmental provenance, or at minimum quantify what fraction of detected AROs come
from clinically-derived models.

### D3 [MAJOR] — PlasmidFinder's Enterobacteriaceae bias compounds the multireplicon section

The manuscript notes that PlasmidFinder types only ~20% of the set and correctly flags the mixed
comparator. It does not note that the untyped fraction is *non-randomly* distributed: PlasmidFinder's
reference replicons are Enterobacteriaceae-centric, so typing coverage is far better in the clinical
and livestock compartments than in environmental ones. That makes the cross-compartment comparison of
multireplicon odds ratios — the section's actual claim — structurally unsound, not merely imprecise.

*Fix:* I would cut the multireplicon section from this paper. Reporting odds ratios the authors
themselves label upper bounds, from a comparator they know is biased differently across the very
groups being compared, weakens an otherwise careful manuscript. It is a good separate analysis once a
replicon-typed-only baseline exists.

### D4 [MINOR] — Possible circularity in the clinical flag

`is_clinical` derives partly from PLSDB disease tags, and PLSDB records are themselves enriched for
resistance-carrying plasmids. Some of the clinical compartment's elevation may be definitional.
Quantify how many clinical assignments rest on the disease tag alone.

### D5 [MINOR] — Missing literature

Smillie *et al.* (2010, *Microbiol Mol Biol Rev* 74:434–452) established the mobility taxonomy the
paper depends on and should be cited where mobility classes are introduced. Given the strong size
effect, the plasmid-size/ARG literature deserves engagement.

**Scores** — Domain accuracy 6/10 · Literature coverage 6/10 · Contribution 7/10.

---

## REVIEWER 3 — PERSPECTIVE

**Recommendation: Minor Revision** (conditional on R1/R2 issues being resolved).

### P1 [MAJOR] — The surveillance recommendation outruns the evidence

The Discussion's operational claim — that surveillance should weight plasmid mobility on par with
provenance — is the paper's most useful contribution and the reason it will be read outside genomics.
But it rests on M3's unadjusted comparison. If the adjusted analysis supports it, state it as a
quantified recommendation (e.g. expected yield per plasmid screened by mobility class); if it does
not, the recommendation should be withdrawn rather than softened. Right now it sits between the two.

### P2 [MAJOR] — The 50% compartment-assignment rate is under-weighted relative to its severity

Limitations names this as the largest threat, which is correct and commendably direct. But the
Abstract states "assigned 104 169 of them to eight One Health compartments" without disclosing that
this is half the collection, and the framing throughout implies broader coverage than exists. A
reader who stops at the Abstract will materially misjudge the evidence base.

*Fix:* state the fraction in the Abstract. Additionally, test whether assignment is associated with
resistance status — compare AMR prevalence in assigned versus unassigned plasmids. That is a
one-line check that would substantially bound the concern either way.

### P3 [MINOR] — The engineered-interface result is the most policy-relevant and is under-developed

The intermediate prevalence, biocide signature and repertoire position of the engineered interface
form the clearest actionable finding for wastewater surveillance, and the Pal *et al.* concordance is
well drawn. This deserves more than the two paragraphs it currently receives.

### P4 [MINOR] — Reporting standard

Consider stating conformity to a reporting checklist for observational genomic studies; NAR readers
increasingly expect it.

**Scores** — Significance 7/10 · Practical impact 7/10 · Clarity 8/10.

---

## DEVIL'S ADVOCATE

### Strongest counter-argument

*The paper's title claim may be an artefact of plasmid size, and the paper contains the evidence to
show it.*

"Mobility rivals provenance" rests on the Table 3 contrast between conjugative environmental plasmids
(27.4%) and non-mobilizable clinical ones (23.9%). But conjugative plasmids in this dataset have a
median length of 94 568 bp against 5744 bp for non-mobilizable ones **[verified]**, and the paper's
own model identifies plasmid length as its single strongest predictor (OR 4.89 per log₁₀ bp,
Table 4). A ~16-fold length difference across roughly 1.2 log₁₀ units predicts a large prevalence gap
on size alone. The comparison that names the paper is thus between two cells that differ enormously
in the very covariate the authors identify as dominant — and it is never adjusted for it.

There is a defensible version of this claim: mobility retains OR 4.77 in the adjusted model, so it is
not merely a size proxy. But that adjusted result is not the comparison the title, Abstract and
Discussion actually deploy. The paper argues its headline from the unadjusted table while its own
adjusted model sits one section away. Either the adjusted contrast is computed and the claim is made
properly, or the title overstates what was shown.

### Issue list

| ID | Severity | Dimension | Location | Issue |
|---|---|---|---|---|
| DA1 | **CRITICAL** | Logic | Title, Abstract, Table 3, Discussion ¶2 | Headline claim argued from unadjusted data while the paper's own model identifies a dominant confounder; adjusted version never computed |
| DA2 | **MAJOR** | Internal coherence | Abstract vs §"Adjustment separates…" | Abstract leads with "24-fold, from 44.6% to 1.8%", anchored on the wildlife compartment — which the Results then show is an artefact that vanishes under adjustment (OR 0.85, *P* = 0.084). The paper's headline range is anchored on a value it later disowns |
| DA3 | **MAJOR** | Alternative explanation | Throughout | Host taxonomy (R2/D1) is a complete alternative explanation for the gradient and is never entertained |
| DA4 | **MAJOR** | Selective presentation | §Multireplicon | Results the authors label upper bounds from a knowingly biased comparator are reported as extending prior work; a caveat is not a substitute for a valid comparator |
| DA5 | MINOR | Overgeneralization | Discussion ¶3 | "We would expect this artefact to affect any habitat-stratified plasmid resistome comparison that does not adjust for size" — generalizes from one dataset to all such studies without support |

### Alternative explanations not considered

1. **Host phylogeny** — the gradient may be a lineage effect (D1).
2. **Database detection bias** — CARD's clinical curation (D2) predicts this gradient with no
   compartment-level selection required.
3. **Assembly and deposition practice** — clinical isolates are more often closed with long reads,
   and completeness correlates with recovering large conjugative plasmids. Compartments may differ in
   the *kind* of plasmid recoverable, not the plasmids present.

### "So what?" test

**Passes, conditionally.** If the gradient survives host-taxon adjustment, the compartment-resolved
resistome map is a genuine contribution and the surveillance implication is actionable. If it does
not survive, the paper still has a publishable and arguably more interesting finding — that apparent
One Health structure in plasmid resistomes is substantially lineage structure. Both outcomes are
worth publishing. What is not currently defensible is the claim as stated.

### Observations (non-defects)

- The authors' own falsification of the plant/food and wildlife effects is intellectually honest and
  rare; it should be foregrounded, not buried mid-Results.
- One suspicion I tested did **not** hold: I hypothesised the two-block Jaccard structure was an
  artefact of differing ARO richness across compartments. **[verified — refuted]** Spearman between
  unique-ARO richness and mean off-diagonal Jaccard is −0.071 (*P* = 0.87). No richness effect. The
  repertoire structure is not explained by compartment size. Note the test has only 8 points and is
  underpowered, so I raise it as insufficient grounds for objection rather than as clearance.

---

## EDITORIAL DECISION

### Decision: **MAJOR REVISION**

### Consensus across the panel

All five seats independently reached the same structural conclusion: the analysis is careful,
transparent and reproducible, and the descriptive results are almost certainly correct — but the
causal-sounding claims are not yet supported by the inference as performed. Three findings were
corroborated by more than one seat:

| Issue | Raised by | Severity |
|---|---|---|
| Central claim argued from unadjusted data | R1 (M3), DA (DA1) | CRITICAL |
| Host-taxon confounding unaddressed | R2 (D1), DA (DA3) | CRITICAL |
| Clustered non-independence invalidates reported precision | R1 (M1) | CRITICAL |
| Mobility/size collinearity undiagnosed | R1 (M2), DA (DA1) | MAJOR |
| CARD clinical curation bias absent from Limitations | R2 (D2) | MAJOR |
| Multireplicon section rests on a biased comparator | R2 (D3), DA (DA4) | MAJOR |
| 50% assignment rate under-disclosed in Abstract | R3 (P2) | MAJOR |

### Adjudication of Devil's Advocate CRITICAL (Checkpoint Rule #4)

**DA1 is validated.** It is independently corroborated by R1's M3, and its factual basis — the
16-fold median size difference between mobility classes — was verified against the source data during
review. It therefore blocks acceptance in the present form and is carried into the roadmap as item 1.
It is not fatal: the adjusted model already exists and needs only to be queried for predicted
probabilities, so this is a tractable revision rather than a redesign.

### Dissent

R3 recommended Minor Revision against Major from the other three seats. R3's assessment was scoped to
framing and impact rather than inference, and R3 explicitly conditioned it on R1/R2 issues being
resolved. The condition is not met, so the decision resolves to Major Revision. R3's dissent is
recorded, not overridden.

### What the paper does well

Stated plainly because it is unusual and should survive revision intact: the analysis is regenerated
end-to-end by one seeded script with pinned tool versions; the compartment mapping is fully
enumerated rather than described; the conservative detection threshold is validated against an
independent tool; and the authors falsify two of their own eight compartment effects. The Limitations
section is more candid than is typical. None of the issues below concern honesty or care — they
concern inferential reach exceeding design.

### Revision roadmap (priority order)

1. **Compute the adjusted headline contrast** (DA1, M3). Predicted probabilities from the Table 4
   model for conjugative-environmental versus non-mobilizable-clinical at fixed plasmid size. Retitle
   and reframe if the ordering does not survive. *This single analysis determines the paper's thesis.*
2. **Add host taxon to the model** (D1, DA3). Report compartment coefficients before and after.
   Whichever way it resolves is reportable; not testing it is not an option.
3. **Refit with cluster-robust standard errors** on `mob_cluster` (M1). Report nominal and clustered
   intervals side by side.
4. **Diagnose the mobility/size collinearity** (M2): VIF, and the mobility coefficient with and
   without `log_size`. Replace "independent" with "separable but correlated" throughout.
5. **Cut or rebuild the multireplicon section** (D3, DA4). Recommend cutting for this paper.
6. **Add CARD curation bias to Limitations** (D2), distinct from the existing threshold discussion.
7. **Disclose the 50% assignment rate in the Abstract**, and test AMR prevalence in assigned versus
   unassigned plasmids (P2).
8. **Resolve the Abstract/Results tension** on the wildlife-anchored 24-fold span (DA2).
9. Minor: Simpson's paradox terminology (M4), multiple-testing families (M5), `is_clinical`
   circularity check (D4), Smillie 2010 citation (D5), soften the Discussion ¶3 generalization (DA5),
   expand the engineered-interface treatment (P3).

Items 1–3 are prerequisites for reassessment. Items 1, 2 and 4 are all computable from the existing
master table with no new annotation.

### Overall scores

| Dimension | Weight | Score |
|---|---:|---:|
| Originality | 20% | 6.5/10 |
| Methodological rigour | 25% | 4.5/10 |
| Evidence sufficiency | 25% | 5.5/10 |
| Argument coherence | 15% | 5.0/10 |
| Writing quality | 15% | 8.0/10 |
| **Weighted total** | | **5.7/10** |

Consistent with Major Revision for NAR. The gap between writing quality and methodological rigour is
the paper's defining characteristic: it is a well-written, transparent presentation of an analysis
whose inferential claims currently outrun its design. Items 1–3 are tractable with data already in
hand, and a revision that addresses them would be a materially stronger paper.
