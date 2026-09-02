# Mobility rivals provenance: plasmid-borne antimicrobial resistance across One Health compartments in 104 169 complete plasmids

*Manuscript draft for Nucleic Acids Research. All quantitative results derive from
`reports/amr_onehealth_methodology.md`, produced by `scripts/analyze_amr_onehealth_mobility.py`.
Author list, affiliations and CRediT statement are deliberately omitted from this draft.*

---

## ABSTRACT

Plasmids carry much of the clinically important resistome, but the relative contribution of *where*
a plasmid was sampled and *how* it transfers has not been measured on a single, uniformly annotated
collection spanning the human, animal and environmental compartments of One Health. We assembled
208 248 complete, non-redundant bacterial plasmids from PlasmidScope, annotated them uniformly with
CARD/RGI, MOB-suite and PlasmidFinder, and assigned 104 169 of them to eight One Health compartments
derived from curated sampling metadata. Resistance prevalence spans 24-fold, from 44.6% of plasmids
in the human clinical compartment to 1.8% in wildlife, with livestock and poultry second at 38.9%.
Mobility is an axis of comparable magnitude: conjugative plasmids carry resistance at 60.1% versus
6.7% for non-mobilizable plasmids, and the ordering is monotonic within every compartment. The two
effects are independent in a logistic model adjusting for plasmid size and GC content, in which
plasmid length is the single largest predictor (OR 4.89 per log₁₀ bp). Adjustment eliminates the
apparent signal in two compartments — plant/food and wildlife — showing them to be compositional
artefacts. Resistance gene repertoires separate into a human block and an agricultural–environmental
block, with the engineered interface (wastewater, bioreactors, built environments) closest to both
and carrying a distinctive biocide co-selection signature. Conjugative plasmids from natural
environments carry resistance more often than non-mobilizable plasmids from the clinic, indicating
that transfer capability deserves parity with sampling origin in genomic AMR surveillance.

**Keywords:** plasmid; antimicrobial resistance; One Health; conjugation; mobility; CARD;
comparative genomics

---

## INTRODUCTION

Bacterial antimicrobial resistance was associated with an estimated 4.95 million deaths in 2019 and
directly attributable to 1.27 million (1), and the mobile genetic elements that disseminate
resistance determinants between lineages are central to that burden. Plasmids occupy a special place
among those elements: they replicate autonomously, frequently encode their own conjugative
machinery, and assemble cargo — resistance, virulence, metal and biocide tolerance — into
transmissible packages that cross species and genus boundaries (2). Mapping plasmid diversity at
scale has consequently become a priority, and the plasmid taxonomic unit framework showed that the
global plasmidome is discretely structured, with individual units differing systematically in host
range (3).

Two limitations recur in this literature. The first is provenance bias. Curated, non-redundant
plasmid collections are assembled predominantly from cultivable organisms of clinical or
agricultural interest, so any resistance enrichment they display is partly a statement about what
gets sequenced. The recent demonstration that more than 30% of 24 000 non-redundant plasmids encode
multiple replicons, and that such multireplicon plasmids are larger, more mobile, broader in host
range and enriched in resistance cargo (4), is a case in point: the underlying collection cannot
distinguish an intrinsic property of plasmid biology from a signature of clinical sampling. The
second limitation is that mobility and provenance are seldom analysed jointly. Studies typically
stratify by habitat *or* report mobility class, leaving open whether a plasmid's resistance burden
tracks the environment it came from, its capacity to transfer, or merely its size.

The One Health framework asserts that human, animal and environmental resistance reservoirs are
interconnected, and sewage metagenomics has demonstrated globally structured resistance abundance
that correlates with socio-economic and health indicators (5). Yet One Health claims about plasmids
are usually assembled from heterogeneous studies using different resistance databases, different
detection thresholds and incompatible habitat vocabularies. What has been missing is a single
collection, uniformly annotated end to end, in which compartment and mobility can be measured on the
same plasmids with the same thresholds.

We provide that here. From PlasmidScope (6) we assembled 208 248 complete, non-redundant plasmids,
annotated the resistome of every one of them with CARD/RGI under a fixed confidence threshold,
typed mobility with MOB-suite and replicons with PlasmidFinder, and reconciled sampling metadata
into eight One Health compartments. We ask three questions. How steep is the resistance gradient
across One Health compartments when detection is held constant? Does mobility act independently of
provenance, or does one explain the other? And which apparent compartment effects survive adjustment
for the compositional properties — plasmid size above all — that differ between compartments?

---

## MATERIALS AND METHODS

### Plasmid collection and working set

Plasmid sequences and metadata were obtained from PlasmidScope (6), which aggregates plasmids from
ten public repositories including IMG/PR (7), PLSDB (8), RefSeq, GenBank, COMPASS, mMGE, ENA and
DDBJ, and which is deduplicated at 100% identity and coverage, so no plasmid is counted twice across
constituent databases. From the 852 600 deduplicated entries we retained only sequences flagged
complete or closed (208 360) and then removed 112 laboratory-made or synthetic constructs, giving a
working set of **208 248 plasmids**. Three plasmids, each contributed by a single source database, had no retrievable nucleotide
sequence, so the sequence-based tools ran on the remaining 208 245.

### Resistome annotation

Antimicrobial resistance genes were called with RGI 6.0.8 against CARD v4.0.1 (9), using
`rgi main --input_type contig -a DIAMOND -n 2 --clean --local -d plasmid --include_loose`. RGI
predicts open reading frames with Prodigal (10) and aligns the inferred proteins to CARD with
DIAMOND 2.2.4 (11). The working set was split into 595 shards and processed as an HPC array job;
all 595 shards completed successfully with no plasmid dropped.

Only **Perfect and Strict** hits were counted, following CARD's standard convention; Loose hits were
retained in the raw output for audit but never contribute to any figure reported here. Under this
threshold 29 396 plasmids (14.1%) carry at least one resistance gene and 20 741 (9.96%) carry genes
against two or more drug classes. Presence/absence agreement with AMRFinderPlus calls curated in
PLSDB, over the 45 743 plasmids present in both, is 90.5%; 93% of the discordant plasmids do carry a
CARD Loose hit, confirming that the disagreement reflects our deliberately conservative threshold
rather than missed genes. Reported prevalences are therefore floors.

### Mobility and replicon typing

Mobility class (conjugative, mobilizable, non-mobilizable), relaxase and MPF content, *oriT*, MOB
cluster and predicted host range were assigned with MOB-suite `mob_typer` (12), which returned a
call for 100% of the working set and resolved 7054 distinct MOB clusters. Replicon and
incompatibility typing used PlasmidFinder (13) at the CGE thresholds of ≥80% identity and ≥60%
coverage; 41 555 plasmids (20%) carry at least one detectable replicon and 24 583 carry two or more.
Gene content and functional categories were annotated with PlasAnn (14).

### One Health compartment assignment

Sampling environment was recovered from GOLD ecosystem paths, NCBI BioSample and SRA records and
reconciled into a locked habitat taxonomy with a raw label available for 85.9% of the set. For this
analysis we deliberately did **not** use the coarse top level of that taxonomy, whose
"Host-associated" category merges human, livestock, poultry, companion animal, insect and plant
samples into one bucket and therefore cannot express the human/animal distinction that One Health
requires. Compartments were instead built from the 40-level sub-habitat field combined with a
clinical flag (clinical body site, curated disease tag, or hospital/patient context), using an
explicit first-match-wins mapping:

| Compartment | Sub-habitats | *n* |
|---|---|---:|
| Human — clinical | eight human sub-habitats where the clinical flag is set | 18 171 |
| Human — community | the same eight where it is not | 27 766 |
| Animal — livestock & poultry | pig, cattle, other livestock, poultry/bird | 5 041 |
| Animal — companion & aquaculture | companion animal, fish/aquaculture | 1 395 |
| Animal — wildlife & other | other mammal, rodent, insect/arthropod, other invertebrate, unspecified-host gut | 7 891 |
| Environment — engineered interface | wastewater/sewage, bioreactor, built environment, solid waste/compost, industrial/remediation | 12 028 |
| Environment — natural | soil, freshwater, marine, air, extreme, unspecified environmental | 24 742 |
| Plant & food | plant, food/fermentation, algae, fungi | 7 135 |
| **Total** | | **104 169** |

Human sub-habitats split on the clinical flag; non-human habitats retain their habitat even when the
flag is set, so that hospital wastewater remains wastewater. Simulated-community and laboratory-
enrichment records were excluded as sampling artefacts rather than habitats, together with unlabelled
and geography-only records: 104 079 plasmids (50.0% of the working set) carry no compartment and are
absent from all compartment analyses.

### Statistical analysis

Prevalences are reported with Wilson 95% confidence intervals. Association between compartment and
resistance status was tested by χ² with Cramér's *V* as effect size; pairwise contrasts against the
natural-environment baseline used Fisher's exact test with Benjamini–Hochberg control of the false
discovery rate (15). Independent contributions of compartment and mobility were estimated by logistic
regression of resistance status on compartment, mobility class, log₁₀ plasmid length and GC content,
with the natural environment and the non-mobilizable class as reference levels. Compartment
resistance-gene repertoires were compared by Jaccard index over the sets of distinct ARO accessions
observed. Robustness to clonal redundancy was assessed by redrawing one plasmid at random per MOB
cluster under a fixed seed. Analyses used Python with pandas, SciPy and statsmodels; the complete
analysis is a single script and reruns are deterministic.

---

## RESULTS

### A 24-fold resistance gradient across One Health compartments

Resistance prevalence differs sharply and systematically between compartments (Table 1; χ² = 15 043.4,
df = 7, *P* < 1 × 10⁻³⁰⁰, Cramér's *V* = 0.380). Nearly half of plasmids sampled from human clinical
sources carry at least one confident resistance gene (44.6%), against fewer than one in fifty from
wildlife (1.8%).

**Table 1.** Resistance prevalence by One Health compartment.

| Compartment | *n* | AMR⁺ | Prevalence | Multidrug | OR vs natural |
|---|---:|---:|---:|---:|---:|
| Human — clinical | 18 171 | 8 111 | 44.6% | 36.4% | 18.8 |
| Animal — livestock & poultry | 5 041 | 1 961 | 38.9% | 30.5% | 14.9 |
| Human — community | 27 766 | 5 873 | 21.2% | 16.8% | 6.3 |
| Animal — companion & aquaculture | 1 395 | 287 | 20.6% | 15.3% | 6.0 |
| Plant & food | 7 135 | 891 | 12.5% | 8.5% | 3.3 |
| Environment — engineered interface | 12 028 | 1 310 | 10.9% | 6.9% | 2.9 |
| Environment — natural | 24 742 | 1 016 | 4.1% | 2.5% | 1 (ref) |
| Animal — wildlife & other | 7 891 | 145 | 1.8% | 1.1% | 0.44 |

All contrasts against the natural-environment baseline remain significant after Benjamini–Hochberg
correction. Multidrug carriage tracks overall prevalence closely, so the gradient reflects the number
of plasmids carrying resistance rather than a change in the intensity of carriage among carriers.

Two features deserve emphasis. Livestock and poultry rank second, close behind the clinical
compartment and well above human community carriage — the clearest agricultural selection signal in
the dataset. And wildlife falls significantly *below* the natural-environment baseline (OR 0.44),
which we return to below, because adjustment changes its interpretation.

### Mobility is an axis of comparable magnitude

Stratifying the same plasmids by MOB-suite mobility class produces a spread of similar size
(Table 2). Conjugative plasmids carry resistance at 60.1%, roughly nine times the 6.7% of
non-mobilizable plasmids, with mobilizable plasmids intermediate at 18.9%.

**Table 2.** Resistance prevalence by mobility class.

| Mobility | *n* | Prevalence | Multidrug |
|---|---:|---:|---:|
| Conjugative | 17 965 | 60.1% | 51.8% |
| Mobilizable | 24 677 | 18.9% | 11.9% |
| Non-mobilizable | 61 526 | 6.7% | 4.7% |

The self-transmissible fraction is where the plasmid resistome concentrates, and the multidrug
figures separate even more strongly than overall prevalence: a conjugative plasmid is eleven times
likelier than a non-mobilizable one to carry genes against two or more drug classes.

### The mobility ordering holds inside every compartment

Because compartments differ in mobility composition, the two univariate gradients could in principle
be the same phenomenon viewed twice. They are not. Crossing the two axes (Table 3) shows the mobility
ordering conjugative > mobilizable > non-mobilizable holding monotonically in all eight compartments,
with no reversal — that is, no Simpson's paradox of the kind that pooled enrichment analyses of
plasmid cargo are vulnerable to.

**Table 3.** Resistance prevalence (%) by compartment and mobility class.

| Compartment | Conjugative | Mobilizable | Non-mobilizable |
|---|---:|---:|---:|
| Human — clinical | 73.6 | 35.2 | 23.9 |
| Human — community | 67.1 | 13.5 | 8.8 |
| Animal — livestock & poultry | 59.5 | 26.8 | 23.5 |
| Animal — companion & aquaculture | 55.6 | 16.9 | 9.4 |
| Animal — wildlife & other | 38.1 | 3.5 | 0.7 |
| Environment — engineered interface | 44.3 | 15.6 | 3.9 |
| Environment — natural | 27.4 | 9.0 | 2.0 |
| Plant & food | 24.7 | 15.2 | 6.7 |

The comparison that most directly motivates our title sits inside this table. A conjugative plasmid
sampled from a natural environment carries resistance at 27.4%, *higher* than a non-mobilizable
plasmid sampled from a human clinical source at 23.9%. Knowing that a plasmid can transfer itself is
at least as informative as knowing it came from a hospital.

### Adjustment separates genuine compartment effects from composition

Compartments differ not only in mobility mix but in plasmid size and base composition, both of which
independently predict gene content. We therefore fitted a logistic model of resistance status on
compartment, mobility, log₁₀ length and GC content (*n* = 104 168, pseudo-*R*² = 0.360; Table 4).

**Table 4.** Adjusted odds of carrying ≥1 resistance gene. Reference: natural environment,
non-mobilizable.

| Term | Adjusted OR | 95% CI | *P* |
|---|---:|---|---|
| Human — clinical | 8.94 | 8.27–9.66 | <1 × 10⁻³⁰⁰ |
| Animal — livestock & poultry | 5.92 | 5.37–6.53 | 2.3 × 10⁻²⁸² |
| Human — community | 4.99 | 4.62–5.40 | <1 × 10⁻³⁰⁰ |
| Environment — engineered interface | 2.94 | 2.67–3.23 | 2.5 × 10⁻¹⁰⁹ |
| Animal — companion & aquaculture | 2.51 | 2.14–2.94 | 1.0 × 10⁻²⁹ |
| Plant & food | 1.09 | 0.98–1.21 | 0.097 |
| Animal — wildlife & other | 0.85 | 0.71–1.02 | 0.084 |
| Conjugative | 4.77 | 4.54–5.02 | <1 × 10⁻³⁰⁰ |
| Mobilizable | 2.77 | 2.63–2.91 | <1 × 10⁻³⁰⁰ |
| log₁₀ length | 4.89 | 4.71–5.09 | <1 × 10⁻³⁰⁰ |
| GC content | 0.99 | 0.990–0.995 | 3.5 × 10⁻¹¹ |

Three results follow. First, **plasmid length is the largest single predictor** in the model, at
OR 4.89 per log₁₀ base pairs: any comparison of resistance carriage that does not control for size is
in part a comparison of size. Second, **the clinical, livestock, community and engineered-interface
effects survive adjustment**, remaining between three- and ninefold above the environmental baseline,
so these are genuine compartment effects and not compositional shadows. Third, **the plant/food and
wildlife effects disappear entirely** (OR 1.09, *P* = 0.097; OR 0.85, *P* = 0.084). Their striking
univariate prevalences of 12.5% and 1.8% are products of the size and mobility distributions of the
plasmids sampled from those settings, not evidence of compartment-level selection. The apparent
depletion of resistance in wildlife, in particular, is not a wildlife effect: it is what small,
non-mobilizable plasmids look like anywhere.

Mobility retains a large independent effect after adjustment (conjugative OR 4.77, mobilizable
OR 2.77), confirming that the two axes contribute separately.

### Resistance repertoires form a human block and an agricultural–environmental block

Prevalence describes how often resistance occurs; it does not say whether the same genes are
involved. Unique ARO richness ranges from 686 accessions among the 8111 clinical carriers to 109
among the 145 wildlife carriers. Pairwise Jaccard comparison of compartment repertoires reveals a
two-block structure. The strongest similarity is between the two human compartments (clinical ↔
community, 0.53). A second block comprises livestock, plant/food, natural and engineered
compartments, which resemble one another at 0.42–0.47 — livestock ↔ plant/food and natural ↔
plant/food both reach 0.47. Sharing *between* the human block and the environmental block is
consistently lower, at 0.33–0.37, and the human clinical ↔ wildlife comparison is the lowest overall
at 0.15.

The human clinical resistome is therefore comparatively distinct in composition, not merely larger.
The engineered interface is the environmental compartment that comes closest to the human
repertoires, consistent with its position as the physical interface between them.

### Compartment-specific drug-class signatures

Drug-class composition among carriers separates the compartments in a way that matches known
selection pressures. β-lactams dominate both human compartments, with penicillin β-lactam resistance
in 66% of clinical carriers, cephalosporin in 58% and carbapenem in 33%. Livestock and poultry
carriers are instead led by aminoglycoside (59%) and tetracycline (53%) resistance, the classes with
the longest veterinary and growth-promotion histories. Plant and food carriers are dominated by
tetracycline (53%) and fluoroquinolone (33%). Most distinctively, resistance to disinfecting agents
and antiseptics appears in 44% of engineered-interface carriers and 32% of natural-environment
carriers, versus a much smaller share of clinical ones — a biocide co-selection signature
concentrated exactly where quaternary ammonium compounds and related biocides are applied.

### Multireplicon architecture and resistance outside the clinic

Plasmids carrying two or more PlasmidFinder replicons are more likely to carry resistance in every
compartment, but the strength of that association is inversely related to the compartment's baseline
burden. The odds ratio is 1.8 in the human clinical compartment, where resistance is common
regardless of architecture, but 17.8 in the natural environment and 14.2 in wildlife. Multireplicon
architecture is thus most informative precisely in the settings where resistance is otherwise rare.
This extends the reported association between multireplicon structure and resistance (4) into
non-clinical compartments that the underlying 24 000-plasmid collection could not resolve. We stress
that these particular odds ratios are upper bounds; the comparator group is heterogeneous, as
detailed under Limitations.

### Clonal redundancy inflates absolute prevalence but preserves the gradient

Public plasmid collections over-represent clonal lineages from intensively sequenced clinical and
agricultural settings. Redrawing one plasmid per MOB cluster (6182 plasmids from 7054 clusters)
approximately halves absolute prevalence — the clinical compartment falls from 44.6% to 23.8%,
livestock from 38.9% to 22.0%, plant/food from 12.5% to 5.9% — while leaving the ordering intact and
the clinical and livestock compartments in the top two positions. Because MOB clustering is coarse,
this procedure collapses genuinely distinct plasmids and so under-counts diversity; the dereplicated
values are a lower bound as surely as the raw values are an upper bound. The gradient itself,
rather than any single percentage, is the robust finding.

---

## DISCUSSION

Analysing provenance and mobility together on one uniformly annotated collection changes what each
appears to explain. The One Health gradient is real and steep: after adjusting for plasmid size, GC
content and mobility, plasmids from human clinical sources remain roughly nine times likelier to
carry resistance than plasmids from natural environments, and plasmids from livestock and poultry
roughly six times. The position of livestock immediately behind the clinic, and clearly ahead of
human community carriage, is consistent with agricultural antimicrobial use acting as a selective
force of magnitude comparable to clinical use.

The finding we consider most consequential for surveillance is that mobility class carries
information of the same order as sampling origin, and does so independently. A conjugative plasmid
from soil or freshwater is likelier to carry a resistance gene than a non-mobilizable plasmid from a
clinical specimen. Genomic AMR surveillance that ranks isolates by provenance while treating plasmid
mobility as a secondary annotation is discarding a predictor of comparable weight. Because
conjugative plasmids are also the subset capable of moving resistance between hosts unaided, the
practical case for weighting them is stronger still than the prevalence figures alone imply.

Our results also caution against a common inferential shortcut. Two of the eight compartment effects
— plant/food and wildlife — do not survive adjustment. Both would have been reportable as findings
from univariate analysis, and the wildlife result in particular is the kind of number that invites a
narrative about pristine environments harbouring little resistance. Once plasmid size and mobility
are accounted for, neither compartment differs from the environmental baseline. Since plasmid length
is the strongest single predictor in our model, compartments that happen to yield small plasmids will
appear resistance-poor for reasons that have nothing to do with selection. We would expect this
artefact to affect any habitat-stratified plasmid resistome comparison that does not adjust for size.

The repertoire analysis adds a compositional dimension that prevalence alone conceals. The two-block
structure — a human block, and an agricultural–environmental block whose members resemble one another
more than either resembles the clinic — indicates that compartments differ in *which* resistance
genes they carry, not merely how many plasmids carry any. The engineered interface behaves as its
name suggests: intermediate in prevalence, closest of the environmental compartments to the human
repertoires, and marked by a distinctive biocide signature. This last observation aligns closely with
earlier work showing that co-occurrence of biocide/metal and antibiotic resistance genes on plasmids
is rare in external environments but substantially more common on plasmids of human and domestic
animal origin, and that quaternary ammonium resistance is among the few biocide determinants
routinely co-located with resistance genes on plasmids (16). Our engineered-interface result is
consistent with wastewater and built environments being where that co-selection is concentrated, and
with sewage-based surveillance capturing a genuinely mixed reservoir (5).

Finally, the multireplicon result illustrates the value of testing a clinically derived rule outside
the clinic. The association between multireplicon architecture and resistance holds in all eight
compartments, but it is weakest where it was originally characterised and strongest in the natural
environment and wildlife. A rule inferred from cultured, largely clinical plasmids is not merely
reproduced outside that setting; its diagnostic value is greater there. Confirming this will require
re-running the comparison against a replicon-typed-only baseline, for the reason given below.

### Limitations

Five constraints bound these conclusions.

**Half the collection lacks a compartment.** Of 208 248 plasmids, 104 079 (50.0%) could not be
assigned — predominantly simulated-community records and entries with no recoverable environment
label. If metadata recovery is non-random with respect to resistance, the gradient is biased by an
unknown amount. This is the single largest threat to our conclusions, and it cannot be resolved
without better source metadata.

**Sampling bias is bounded, not removed.** Clinical and livestock isolates are sequenced for reasons
correlated with resistance. The gradient describes the plasmids that were deposited in public
databases, not the environment as it is. Our dereplication analysis brackets the effect but does not
eliminate it.

**The design is cross-sectional and associational.** Shared repertoires between compartments
demonstrate compositional overlap. They carry no directionality and are not evidence of transmission
between compartments; nothing here supports a claim about flow from farm to clinic or clinic to
river.

**The multireplicon comparator is heterogeneous.** PlasmidFinder types only about 20% of the working
set, so the "fewer than two replicons" comparator is dominated by plasmids with *no detectable*
replicon, which skew small, cryptic and resistance-negative. This inflates the odds ratios reported
above, most severely in the environmental compartments where typing coverage is lowest and where the
largest ratios appear. These should be read as upper bounds pending a replicon-typed-only comparison.

**Detection is deliberately conservative.** Restricting to CARD Perfect and Strict hits means
approximately 93% of plasmids that AMRFinderPlus calls resistant but CARD does not still carry a CARD
Loose hit. All prevalences reported here are floors, and the compartment gradient would be expected
to shift upward, though not necessarily uniformly, under a permissive threshold.

### Conclusion

Across 104 169 uniformly annotated complete plasmids, resistance carriage spans 24-fold between One
Health compartments, and mobility class spans a comparable range independently of provenance. Plasmid
size is a stronger predictor than most compartment memberships, and controlling for it dissolves two
of eight apparent compartment effects. The compartments differ in which resistance genes they carry
as well as how often, with the engineered interface positioned between the human and environmental
blocks and marked by biocide co-selection. For surveillance, the operational implication is that a
plasmid's capacity to transfer itself deserves parity with the setting it was isolated from.

---

## DATA AVAILABILITY

All source data are public. Plasmid sequences and metadata are from PlasmidScope
(https://plasmid.deepomics.org/), with constituent records traceable to IMG/PR
(https://img.jgi.doe.gov/pr/), PLSDB (https://www.ccb.uni-saarland.de/plsdb2025), RefSeq, GenBank,
COMPASS, mMGE, ENA and DDBJ. Resistance annotation used CARD v4.0.1 (https://card.mcmaster.ca/).

The analysis is a single deterministic script, `scripts/analyze_amr_onehealth_mobility.py`, which
regenerates every number, table and figure reported here from the per-plasmid master table. Its ten
output tables and six figures, together with the full methodological write-up including the complete
compartment mapping, are provided as supplementary material. Tool versions are pinned (CARD v4.0.1,
RGI 6.0.8, DIAMOND 2.2.4, BLAST 2.16.0+); the dereplication draw is seeded, and reruns reproduce
identical values.

## SUPPLEMENTARY DATA

Supplementary tables S1–S10 (compartment prevalence, mobility prevalence, compartment × mobility
joint counts, logistic-regression coefficients, pairwise contrasts, drug class by compartment,
repertoire Jaccard matrix, ARO richness, multireplicon contrasts, dereplication sensitivity) and
figures S1–S6 accompany this manuscript.

## FUNDING

[To be completed.]

## CONFLICT OF INTEREST

None declared.

## ETHICS DECLARATION

This study analysed only publicly available, previously published sequence data and associated
metadata. No human subjects, animal subjects, or identifiable personal data were involved, and no
ethical approval was required. Sample-level clinical descriptors used for compartment assignment
(body site, disease tag, hospital context) derive from public database annotations and carry no
individual identifiers.

## DECLARATION ON THE USE OF AI TOOLS

A large language model (Claude, Anthropic) was used to assist with the implementation of the
analysis script, the statistical workflow, and the drafting and editing of this manuscript. All
quantitative results were generated by the cited script executed on the source data, and were
verified against its output; all references were checked against primary sources. The authors take
full responsibility for the content, the correctness of the analyses, and the conclusions drawn.

---

## REFERENCES

1. Murray,C.J.L., Ikuta,K.S., Sharara,F., Swetschinski,L., Robles Aguilar,G., Gray,A., Han,C.,
   Bisignano,C., Rao,P., Wool,E. *et al.* (2022) Global burden of bacterial antimicrobial resistance
   in 2019: a systematic analysis. *Lancet*, **399**, 629–655. doi:10.1016/S0140-6736(21)02724-0

2. Rodríguez-Beltrán,J., DelaFuente,J., León-Sampedro,R., MacLean,R.C. and San Millán,Á. (2021)
   Beyond horizontal gene transfer: the role of plasmids in bacterial evolution. *Nat. Rev.
   Microbiol.*, **19**, 347–359. doi:10.1038/s41579-020-00497-1

3. Redondo-Salvo,S., Fernández-López,R., Ruiz,R., Vielva,L., de Toro,M., Rocha,E.P.C.,
   Garcillán-Barcia,M.P. and de la Cruz,F. (2020) Pathways for horizontal gene transfer in bacteria
   revealed by a global map of their plasmids. *Nat. Commun.*, **11**, 3602.
   doi:10.1038/s41467-020-17278-2

4. de Quinto,I., Ramiro Martínez,P., Silva Rosa,R., Herencias,C., Lanza,V.F. and
   Rodríguez-Beltrán,J. (2026) Multireplicon plasmids emerge under predictable rules and drive the
   spread of antimicrobial resistance across bacterial hosts. *bioRxiv*, 2026.05.05.722934.
   doi:10.64898/2026.05.05.722934

5. Hendriksen,R.S., Munk,P., Njage,P., van Bunnik,B., McNally,L., Lukjancenko,O., Röder,T.,
   Nieuwenhuijse,D., Pedersen,S.K., Kjeldgaard,J. *et al.* (2019) Global monitoring of antimicrobial
   resistance based on metagenomics analyses of urban sewage. *Nat. Commun.*, **10**, 1124.
   doi:10.1038/s41467-019-08853-3

6. Li,Y., Feng,X., Chen,X., Yang,S., Zhao,Z., Chen,Y. and Li,S.C. (2025) PlasmidScope: a
   comprehensive plasmid database with rich annotations and online analytical tools. *Nucleic Acids
   Res.*, **53**, D179–D188. doi:10.1093/nar/gkae930

7. Camargo,A.P., Call,L., Roux,S., Nayfach,S., Huntemann,M., Palaniappan,K., Ratner,A., Chu,K.,
   Mukherjeep,S., Reddy,T.B.K. *et al.* (2024) IMG/PR: a database of plasmids from genomes and
   metagenomes with rich annotations and metadata. *Nucleic Acids Res.*, **52**, D164–D173.
   doi:10.1093/nar/gkad964

8. Molano,L.G., Hirsch,P., Hannig,M., Müller,R. and Keller,A. (2025) The PLSDB 2025 update: enhanced
   annotations and improved functionality for comprehensive plasmid research. *Nucleic Acids Res.*,
   **53**, D189–D196. doi:10.1093/nar/gkae1095

9. Alcock,B.P., Huynh,W., Chalil,R., Smith,K.W., Raphenya,A.R., Wlodarski,M.A., Edalatmand,A.,
   Petkau,A., Syed,S.A., Tsang,K.K. *et al.* (2023) CARD 2023: expanded curation, support for
   machine learning, and resistome prediction at the Comprehensive Antibiotic Resistance Database.
   *Nucleic Acids Res.*, **51**, D690–D699. doi:10.1093/nar/gkac920

10. Hyatt,D., Chen,G.-L., LoCascio,P.F., Land,M.L., Larimer,F.W. and Hauser,L.J. (2010) Prodigal:
    prokaryotic gene recognition and translation initiation site identification. *BMC
    Bioinformatics*, **11**, 119. doi:10.1186/1471-2105-11-119

11. Buchfink,B., Reuter,K. and Drost,H.-G. (2021) Sensitive protein alignments at tree-of-life scale
    using DIAMOND. *Nat. Methods*, **18**, 366–368. doi:10.1038/s41592-021-01101-x

12. Robertson,J. and Nash,J.H.E. (2018) MOB-suite: software tools for clustering, reconstruction and
    typing of plasmids from draft assemblies. *Microb. Genom.*, **4**, e000206.
    doi:10.1099/mgen.0.000206

13. Carattoli,A., Zankari,E., García-Fernández,A., Voldby Larsen,M., Lund,O., Villa,L.,
    Møller Aarestrup,F. and Hasman,H. (2014) In silico detection and typing of plasmids using
    PlasmidFinder and plasmid multilocus sequence typing. *Antimicrob. Agents Chemother.*, **58**,
    3895–3903. doi:10.1128/AAC.02412-14

14. Islam,H., Sharma,A., Blair,J. and Lopatkin,A.J. (2026) PlasAnn: a curated plasmid-specific
    database and annotation pipeline for standardized gene and function analysis. *Nucleic Acids
    Res.*, **54**, gkaf1507. doi:10.1093/nar/gkaf1507

15. Benjamini,Y. and Hochberg,Y. (1995) Controlling the false discovery rate: a practical and
    powerful approach to multiple testing. *J. R. Stat. Soc. B*, **57**, 289–300.

16. Pal,C., Bengtsson-Palme,J., Kristiansson,E. and Larsson,D.G.J. (2015) Co-occurrence of resistance
    genes to antibiotics, biocides and metals reveals novel insights into their co-selection
    potential. *BMC Genomics*, **16**, 964. doi:10.1186/s12864-015-2153-5

---

## FIGURE LEGENDS

**Figure 1.** Plasmid-borne resistance prevalence across eight One Health compartments, with Wilson
95% confidence intervals (`reports/figures/onehealth/amr_prevalence_by_compartment.png`).

**Figure 2.** Resistance prevalence by compartment and mobility class, with cell sample sizes; the
mobility ordering is monotonic in every compartment
(`reports/figures/onehealth/amr_compartment_x_mobility.png`).

**Figure 3.** Adjusted odds of resistance carriage from logistic regression controlling for plasmid
size and GC content, reference natural environment / non-mobilizable
(`reports/figures/onehealth/amr_adjusted_odds_ratios.png`).

**Figure 4.** Drug-class spectrum of the plasmid resistome by compartment, as the fraction of
carriers in each compartment carrying each class
(`reports/figures/onehealth/drug_class_by_compartment.png`).

**Figure 5.** Pairwise Jaccard similarity of compartment resistance-gene repertoires over distinct
ARO accessions (`reports/figures/onehealth/arg_repertoire_jaccard.png`).

**Figure 6.** Clonal-redundancy robustness: prevalence computed on all plasmids versus one plasmid
per MOB cluster (`reports/figures/onehealth/dereplication_sensitivity.png`).
