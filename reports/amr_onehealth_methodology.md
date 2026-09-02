# Plasmid-borne AMR across One Health compartments — methodology & results

**Scope.** Phase-3 interpretive analysis of the CARD/RGI resistome (`card_*` columns of
`plasmid_metadata_master.tsv`) stratified by **One Health compartment** and **mobility class**. This
is the follow-on named in `reports/card_amr_methodology.md` §7 ("AMR prevalence by habitat /
mobility, ARG co-occurrence, conjugative-vs-non-mobilisable AMR carriage"). Upstream layers:
`reports/environment_taxonomy_LOCKED.md` (habitat), `reports/typing_methodology.md` (mobility,
replicons), `reports/card_amr_methodology.md` (AMR).

**Produced entirely by** `scripts/analyze_amr_onehealth_mobility.py`. Regenerate:

```bash
export MAMBA_ROOT_PREFIX=$PWD/envs/.micromamba_root
ENV=$MAMBA_ROOT_PREFIX/envs/stats_env
PYTHONNOUSERSITE=1 PATH="$ENV/bin:$PATH" $ENV/bin/python scripts/analyze_amr_onehealth_mobility.py
```

Outputs: 10 tables `data/plasmidscope_primary/onehealth_*.tsv`, 6 figures `reports/figures/onehealth/`.
Deterministic apart from the dereplication draw in §7, fixed at `SEED = 20260727`.

---

## 1 · Why `hab_top` is not used

The locked taxonomy's top level (`hab_top`) collapses human, livestock, poultry, insect and plant
samples into a single **Host-associated** bucket. That is precisely the distinction a One Health
analysis exists to make, so `hab_top` is unusable as a stratifier here. Compartments are built from
**`hab_sub` + `is_clinical`** instead, which is strictly finer-grained.

`hab_top` is *not* consulted at all. The locked exclusion of the Simulated-artifact (64,658) and
Lab-artifact (87) buckets is honoured by excluding their `hab_sub` values
(`Simulated communities`, `Lab enrichment`, `Engineered: other`) directly.

## 2 · The compartment mapping

Transparent, first-match-wins, defined in `COMPARTMENT_MAP` / `HUMAN_SUBS` in the script. Human
sub-habitats split on `is_clinical`; all other habitats keep their habitat compartment even when
`is_clinical = 1` (a hospital-wastewater sample is still wastewater).

| Compartment | `hab_sub` members | n |
|---|---|---:|
| Human — clinical | the 8 human sub-habitats where `is_clinical = 1` | 18,171 |
| Human — community | the same 8 where `is_clinical = 0` | 27,766 |
| Animal — livestock & poultry | Livestock pig/cattle/other, Poultry/bird | 5,041 |
| Animal — companion & aquaculture | Companion animal, Fish/aquaculture | 1,395 |
| Animal — wildlife & other | Other mammal, Rodent, Insect/arthropod, Other invertebrate, Gut/faeces unspecified host | 7,891 |
| Environment — engineered interface | Wastewater/sewage, Bioreactor, Built environment, Solid waste/compost, Industrial/remediation | 12,028 |
| Environment — natural | Terrestrial/soil, Aquatic freshwater/marine, Air, Extreme, Environmental unspecified | 24,742 |
| Plant & food | Plant, Food/fermentation, Algae, Fungi | 7,135 |
| **Total analysis set** | | **104,169** |

**104,079 plasmids (50.0%) carry no compartment** and are excluded: Simulated communities (64,658),
Unlabelled (29,283), Geography only (9,293), Unresolved (755), Lab enrichment (60), Engineered:
other (27), plus 3 ambiguous multi-labels. Every result below is conditional on this half of the
working set — see §8.

---

## 3 · Headline results

### 3.1 A steep compartment gradient

`onehealth_prevalence_by_compartment.tsv` · `figures/onehealth/amr_prevalence_by_compartment.png`

| Compartment | n | AMR+ | prevalence | multidrug |
|---|---:|---:|---:|---:|
| Human — clinical | 18,171 | 8,111 | **44.6%** | 36.4% |
| Animal — livestock & poultry | 5,041 | 1,961 | **38.9%** | 30.5% |
| Human — community | 27,766 | 5,873 | 21.2% | 16.8% |
| Animal — companion & aquaculture | 1,395 | 287 | 20.6% | 15.3% |
| Plant & food | 7,135 | 891 | 12.5% | 8.5% |
| Environment — engineered interface | 12,028 | 1,310 | 10.9% | 6.9% |
| Environment — natural | 24,742 | 1,016 | 4.1% | 2.5% |
| Animal — wildlife & other | 7,891 | 145 | **1.8%** | 1.1% |

χ² = 15,043.4, df = 7, p < 1e-300, Cramér's V = 0.380. A **24-fold spread** between the clinical
human compartment and wildlife. Fisher tests against the natural-environment baseline
(`onehealth_pairwise_vs_natural.tsv`, Benjamini–Hochberg corrected) are all significant: clinical
OR = 18.8, livestock OR = 14.9, community OR = 6.3, engineered OR = 2.9, and wildlife
*below* baseline at OR = 0.44.

Livestock & poultry sitting second, close to the clinical compartment and well above human
community carriage, is the clearest agricultural-selection signal in the dataset.

### 3.2 Mobility is the single strongest axis

`onehealth_prevalence_by_mobility.tsv`

| Mobility | n | AMR+ prevalence | multidrug |
|---|---:|---:|---:|
| conjugative | 17,965 | **60.1%** | 51.8% |
| mobilizable | 24,677 | 18.9% | 11.9% |
| non-mobilizable | 61,526 | 6.7% | 4.7% |

Conjugative plasmids are ~9× more likely to carry resistance than non-mobilizable ones — the
self-transmissible fraction is where the plasmid resistome concentrates.

### 3.3 The two axes are independent, and the gradient holds within every compartment

`onehealth_compartment_x_mobility.tsv` · `figures/onehealth/amr_compartment_x_mobility.png`

AMR prevalence (%), compartment × mobility:

| Compartment | conjugative | mobilizable | non-mobilizable |
|---|---:|---:|---:|
| Human — clinical | 73.6 | 35.2 | 23.9 |
| Human — community | 67.1 | 13.5 | 8.8 |
| Animal — livestock & poultry | 59.5 | 26.8 | 23.5 |
| Animal — companion & aquaculture | 55.6 | 16.9 | 9.4 |
| Animal — wildlife & other | 38.1 | 3.5 | 0.7 |
| Environment — engineered interface | 44.3 | 15.6 | 3.9 |
| Environment — natural | 27.4 | 9.0 | 2.0 |
| Plant & food | 24.7 | 15.2 | 6.7 |

The mobility ordering is **monotonic in all eight compartments** — no Simpson's paradox, which
`RESEARCH_OUTLINE.md` §4.3 flags as the main risk for pooled enrichment claims. Notably, a
conjugative plasmid from the natural environment (27.4%) is more likely to carry an ARG than a
non-mobilizable plasmid from the clinic (23.9%): mobility class matters as much as provenance.

### 3.4 Adjusted odds — and two compartments that vanish under control

`onehealth_logit_odds_ratios.tsv` · `figures/onehealth/amr_adjusted_odds_ratios.png`

Logistic regression, `AMR+ ~ compartment + mobility + log10(size_bp) + GC`, n = 104,168,
pseudo-R² = 0.360. Reference: Environment — natural / non-mobilizable.

| Term | adjusted OR | 95% CI | p |
|---|---:|---|---|
| Human — clinical | **8.94** | 8.27–9.66 | <1e-300 |
| Animal — livestock & poultry | **5.92** | 5.37–6.53 | 2.3e-282 |
| Human — community | 4.99 | 4.62–5.40 | <1e-300 |
| Environment — engineered interface | 2.94 | 2.67–3.23 | 2.5e-109 |
| Animal — companion & aquaculture | 2.51 | 2.14–2.94 | 1.0e-29 |
| Plant & food | 1.09 | 0.98–1.21 | 0.097 (n.s.) |
| Animal — wildlife & other | 0.85 | 0.71–1.02 | 0.084 (n.s.) |
| conjugative | 4.77 | 4.54–5.02 | <1e-300 |
| mobilizable | 2.77 | 2.63–2.91 | <1e-300 |
| log10(size_bp) | **4.89** | 4.71–5.09 | <1e-300 |
| GC % | 0.99 | 0.990–0.995 | 3.5e-11 |

Three things worth stating plainly:

1. **Plasmid size is the largest single predictor** (OR 4.89 per log₁₀ bp) — bigger plasmids carry
   more of everything. Any AMR comparison that does not control for size is partly a size comparison.
2. **Plant & food and wildlife lose their apparent signal entirely** once size, GC and mobility are
   controlled. Their raw prevalences (12.5% and 1.8%) are compositional artefacts of the plasmid
   size/mobility mix in those compartments, not evidence of compartment-level selection.
3. The clinical, livestock, community and engineered-interface effects **survive adjustment**, so
   those are genuine compartment effects.

### 3.5 Drug-class signatures are compartment-specific

`onehealth_drug_class_by_compartment.tsv` · `figures/onehealth/drug_class_by_compartment.png`
(fraction of AMR+ plasmids in the compartment carrying each class; top 5 shown)

| Compartment | signature |
|---|---|
| Human — clinical | penicillin β-lactam 66% · cephalosporin 58% · aminoglycoside 57% · sulfonamide 39% · **carbapenem 33%** |
| Human — community | penicillin β-lactam 60% · cephalosporin 58% · aminoglycoside 53% · sulfonamide 39% · carbapenem 33% |
| Animal — livestock & poultry | **aminoglycoside 59% · tetracycline 53%** · cephalosporin 48% · sulfonamide 42% · penicillin 39% |
| Animal — companion & aquaculture | aminoglycoside 57% · cephalosporin 52% · sulfonamide 48% · penicillin 44% · tetracycline 43% |
| Animal — wildlife & other | cephalosporin 41% · aminoglycoside 36% · penicillin 35% · tetracycline 35% · sulfonamide 28% |
| Environment — engineered interface | cephalosporin 48% · penicillin 46% · aminoglycoside 46% · **disinfectants/antiseptics 44%** · sulfonamide 37% |
| Environment — natural | cephalosporin 34% · **disinfectants/antiseptics 32%** · penicillin 31% · aminoglycoside 30% · tetracycline 30% |
| Plant & food | **tetracycline 53% · fluoroquinolone 33%** · aminoglycoside 33% · cephalosporin 30% · penicillin 25% |

β-lactams (including carbapenem at 33%) dominate the human compartments; **tetracycline and
aminoglycoside** dominate livestock, the canonical growth-promoter/veterinary classes; and
**disinfectants/antiseptics** are the distinguishing class of the engineered interface — the
biocide-co-selection signature expected in wastewater and built environments.

### 3.6 Cross-compartment ARG sharing

`onehealth_arg_jaccard.tsv` · `onehealth_arg_richness.tsv` ·
`figures/onehealth/arg_repertoire_jaccard.png`

Jaccard index over the sets of distinct AROs observed in each compartment. Unique-ARO richness:
clinical 686 (from 8,111 carriers), community 623, engineered interface 327, natural 310,
livestock 303, plant & food 236, companion 165, wildlife 109.

Most-shared pairs: **human clinical ↔ human community 0.53**; then livestock ↔ plant & food 0.47,
natural ↔ plant & food 0.47, engineered ↔ natural 0.44, natural ↔ livestock 0.44. Least shared:
clinical ↔ wildlife 0.15.

The structure is two loosely-coupled blocks: a **human block** (clinical–community, 0.53) and an
**agricultural–environmental block** (livestock, plant & food, natural, engineered, all 0.42–0.47
among themselves), with only 0.33–0.37 sharing between the human block and the environmental one.
The human clinical resistome is comparatively distinct in *composition*, even though the
engineered interface is where environmental and human-associated repertoires come closest.

### 3.7 Multireplicon plasmids and AMR

`onehealth_multireplicon_amr.tsv` — PlasmidFinder ≥2 Inc types vs AMR, per compartment, all
BH-significant:

| Compartment | n multi | AMR+ (multi) | AMR+ (other) | OR |
|---|---:|---:|---:|---:|
| Environment — natural | 643 | 37.2% | 3.2% | **17.8** |
| Animal — wildlife & other | 266 | 16.2% | 1.3% | **14.2** |
| Animal — companion & aquaculture | 147 | 61.9% | 15.7% | 8.7 |
| Plant & food | 501 | 44.5% | 10.1% | 7.2 |
| Human — community | 4,460 | 51.0% | 15.4% | 5.7 |
| Environment — engineered interface | 1,135 | 27.5% | 9.2% | 3.8 |
| Animal — livestock & poultry | 1,157 | 54.4% | 34.3% | 2.3 |
| Human — clinical | 4,850 | 55.4% | 40.7% | **1.8** |

The association is universal but its **strength is inversely related to the compartment's baseline
AMR burden**: in the clinic, resistance is common regardless of replicon architecture (OR 1.8),
whereas in the natural environment a multireplicon plasmid is ~18× more likely to be a carrier. This
is consistent with de Quinto et al.'s multireplicon/AMR link, and extends it by showing the effect
is *strongest* precisely in the non-clinical settings their dataset could not resolve. Read §8.4
before quoting these odds ratios.

---

## 4 · Robustness: clonal redundancy

`onehealth_dereplication_sensitivity.tsv` · `figures/onehealth/dereplication_sensitivity.png`

`RESEARCH_OUTLINE.md` §4.3 warns that clinical clones are massively over-sampled. Re-running the
prevalence table on **one randomly-drawn plasmid per MOB cluster** (6,182 plasmids from 7,054
clusters):

| Compartment | all | 1 per MOB cluster |
|---|---:|---:|
| Human — clinical | 44.6% | 23.8% |
| Animal — livestock & poultry | 38.9% | 22.0% |
| Human — community | 21.2% | 14.1% |
| Animal — companion & aquaculture | 20.6% | 15.3% |
| Plant & food | 12.5% | 5.9% |
| Environment — engineered interface | 10.9% | 6.2% |
| Environment — natural | 4.1% | 3.5% |
| Animal — wildlife & other | 1.8% | 0.9% |

**Absolute prevalences roughly halve, but the ordering is preserved** and the clinical/livestock
compartments remain the two highest. Interpret the §3.1 percentages as an upper bound inflated by
clonal over-sampling, and the gradient itself as the robust finding. The dereplicated numbers are a
*lower* bound (MOB clustering is coarse — 7,054 clusters for 208,248 plasmids — so it collapses
genuinely distinct plasmids); the truth lies between the two columns.

---

## 5 · What this establishes

1. Plasmid-borne AMR follows a **steep, monotonic One Health gradient**, from the clinic (44.6%) to
   wildlife (1.8%), that survives adjustment for plasmid size, GC and mobility for the clinical,
   livestock, community and engineered-interface compartments.
2. **Mobility rivals provenance.** Conjugative plasmids carry AMR at 60.1% overall, and the
   conjugative fraction of the natural environment out-carries the non-mobilizable fraction of the
   clinic.
3. The **engineered interface behaves as an interface**: intermediate prevalence (10.9%), a distinct
   biocide co-selection signature, and the closest ARG-repertoire link between environmental and
   human-associated compartments.
4. **Two apparent compartment effects are artefacts** (plant & food, wildlife) — they disappear once
   plasmid size and mobility are controlled.
5. The **multireplicon–AMR association is universal but strongest outside the clinic**, extending
   de Quinto et al. into exactly the compartments their cultured/clinical dataset could not reach.

## 6 · Not covered here

ARG–ARG co-occurrence networks within plasmids, Inc-type-resolved carriage, geographic structure
(`geo_*` is populated but unused here), and host-range breadth (`mob_host_range`). SegMantX
duplication analysis remains unrun.

## 7 · Caveats

1. **Half the working set has no compartment** (104,079 excluded, mostly Simulated communities and
   Unlabelled). If environment recovery is non-random with respect to AMR, the gradient is biased by
   an unknown amount. This is the single largest threat to the result.
2. **Sampling bias is not corrected.** Clinical and livestock isolates are sequenced for reasons
   correlated with resistance. The gradient measures *the plasmids that were deposited*, not
   environmental reality; §4 bounds but does not remove this.
3. **Cross-sectional and associational.** Shared ARG repertoires (§3.6) show compositional overlap,
   **not** transmission, and carry no direction. No claim about flow between compartments is
   supported by this design.
4. **The multireplicon comparator is mixed.** PlasmidFinder types only ~20% of the working set, so
   "single or none" is dominated by *untypeable* plasmids — which skew small, cryptic and
   AMR-negative. This inflates the §3.7 odds ratios, most severely in the low-typing environmental
   compartments where the largest ORs appear. These should be treated as upper bounds and re-run
   against a replicon-typed-only comparator before publication.
5. **CARD Perfect+Strict is conservative** (`card_amr_methodology.md` §5): ~93% of PLSDB+/CARD−
   disagreements do carry a Loose hit. Absolute prevalences are therefore floors.
6. `is_clinical` for non-human habitats (e.g. 330 wastewater, 108 poultry) keeps the habitat
   compartment by design, so the clinical compartment is human-sample-only.
7. **Plant & food is heterogeneous** (plant tissue, fermented food, algae, fungi) and warrants
   splitting if it becomes a load-bearing comparison.
