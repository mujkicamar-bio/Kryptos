# Domain review — Reviewer 2 (plasmid biology / MGE ecology)

**Date:** 2026-09-02 · **Scope:** biology, not arithmetic. `review/VERIFICATION_PACK.md` read first;
the numbers are taken as confirmed and are not re-litigated except where a *definition* rather than a
computation is at fault.

**Read:** `PROJECT_OVERVIEW.md`, `RESEARCH_OUTLINE.md`, `METHODOLOGY.txt`, `papers_of_interest.txt`,
`reports/{typing_methodology,environment_taxonomy_LOCKED,small_cryptic_methodology,dark_orf_clustering,
pfam_dark_validation,clustering_threshold_sensitivity,amr_onehealth_methodology}.md`,
`manuscript/reassessment_round1.md`, `reports/dark_plasmidome_deck.html`, and the markdown narrative of
all nine `notebooks_my/*.ipynb`.

**Read-only analyses I ran** (nothing in the repo was modified; all working files in scratch):
a permissive Pfam scan of `data/zoom_two_orfs/two_reps.faa`; `phmmer` of those two ORFs and of a
seeded random sample of 500 "dark-only" ORFs against the project's own 78,397 named plasmid proteins;
a category census of `data/plasann_run/annot/*shard_*.tsv.gz` restricted to
`data/dark_orf_run/cryptic_small_ids.txt`; and a decomposition of `mob_rep_types` and of the oriT
columns in `plasmid_metadata_master.tsv`.

---

## Summary

The engineering is strong and the self-criticism is unusually honest, but three of the four headline
biological claims rest on definitions that do not mean what the text says they mean. I tested them
directly. (1) "85.9% of dark ORFs have no named homolog anywhere in the annotated set" is a
**coverage-threshold artefact**: 22.8% of a random sample of those ORFs have a named homolog in the
project's *own* corpus at E < 1e-5, down to E = 2.6e-204, and they are `rep` and `mob` genes.
(2) The ColE1 cassette in `zoom_two_dark_orfs.ipynb` is almost certainly the **mbe/mob accessory
module**, and PlasAnn already names its orthologues elsewhere in the same dataset as `traD` and
`sugE`. (3) "Non-mobilizable" dispersal is circular — `mob_orit` is populated for **zero** of 79,388
non-mobilizable plasmids by construction, while PlasAnn independently flags an oriT on 3,070 of them.
"Payload-free" is a statement about four columns, not about a plasmid. Literature engagement is the
weakest axis: the entire de la Cruz/Rocha mobilization corpus and all of ColE1 molecular biology are
absent.

---

## Findings

### [CRITICAL] F1 — "No named homolog" is a length filter, not a homology test; the dark-only fraction is inflated by ~a fifth

**Evidence anchor.** `reports/dark_orf_clustering.md` R3: *"**85.9% of dark ORFs belong to families
with no named homolog anywhere in the annotated set**, and **5,028 dark-only families have ≥10
members**"*. Produced by `mmseqs easy-cluster … --min-seq-id 0.3 -c 0.8 --cov-mode 0`
(`scripts/cluster_dark_orfs.sh`).

I reconstructed the dark-only set from `data/dark_orf_run/mix30_cluster.tsv` (325,333 dark ORFs in
clusters with no `|N` member — matches the report exactly), drew a seeded random sample of 500, and
searched them with `phmmer` against the 78,397 named proteins extracted from the project's own
`all_cds.faa`:

| | n / 500 | % |
|---|---:|---:|
| named homolog at E < 1e-5 | **114** | **22.8** |
| named homolog at E < 1e-20 | 43 | 8.6 |
| best-hit E-value, min / median | 2.6e-204 / 4.1e-14 | |

The PlasAnn gene names of those homologs are the plasmid backbone, again:
`repB`(5) `repA`(5) `repL`(7) `reP`(10) `repM`(3) · `mobV`(10) `mobC`(8) `mobA`(4) `mobQ`(3)
`mobB`(2) `traG`(3) `rlxA`(2) · `dinJ` `higA` `parE` `copG` `immR` `exc1`(2) `roM` — by category,
35 replication/maintenance, 33 conjugation, 12 toxin–antitoxin, 6 metabolism.

**What is wrong biologically.** `-c 0.8 --cov-mode 0` demands 80% coverage *of both* sequences.
Plasmid `rep` and `mob` genes are the single most fusion-, truncation- and frameshift-prone gene class
in the mobilome, and MMseqs2 at 30% identity is a fast clustering heuristic, not a sensitive homology
search — profile-free HMM search routinely detects homology at 12–20% identity that MMseqs2 at 30%
cannot see. So "dark-only" measures *"MMseqs2 did not put it in the same box"*, and the report reads
it as *"no homolog exists in the corpus"*. The length test rules out the simplest story: only 25% of
the 114 are short fragments of their homolog (median length ratio 1.01, IQR 0.77–1.18); most are
full-length proteins that MMseqs2 simply missed. The project already noticed the mechanism —
`zoom_two_dark_orfs.ipynb` §6 shows a 208 aa representative and a 155 aa protein *"can never share a
family at any identity (0.75 < 0.80)"* — but drew only a ranking lesson from it, not the consequence
for R3 and for `dark_novel_families.tsv`.

**Concrete fix.** Replace the co-clustering test with a sensitive search: `phmmer` (or `jackhmmer`,
2 iterations) of every dark family representative against `named.faa`, at E < 1e-5. Report the
corrected dark-only percentage — my sample puts it near **66%, not 85.9%** — and re-filter
`data/pfam_run/dark_novel_families.tsv` and `core_dark_families.tsv` through the same screen before
either is offered as a deliverable. Cost: minutes.

---

### [CRITICAL] F2 — Even the five "strongest candidates" have not been screened this way; one is an entry-exclusion protein

**Evidence anchor.** `notebooks_my/dark_families_90pct.ipynb`: *"Short proteins held at ≥90%
amino-acid identity across 14–42 independent MOB lineages … with no Pfam-A match at `--cut_ga`.
That combination is hard to explain by sampling"*; the deck offers them as *"the sequences I would
like this room to look at … a plausible synthesis order"*.

`phmmer` of `data/dark_orf_run/recluster/core_dark_families_90.faa` against the named proteome:

| candidate | len | best named homolog | E |
|---|---:|---|---:|
| `IMGPR_plasmid_2551306362_000001\|6` (*A. pittii*) | 105 | `IMGPR_plasmid_2648501511_000003` CDS 1 = **`excL` "entry exclusion protein 1 (Exc1)"** | 2.5e-08 |
| `IMGPR_plasmid_2903362858_000001\|1` (*K. pneumoniae*) | 78 | `GenBank_CP039323.1` CDS 1 = `insB` (IS1 transposase B) | 3.5e-09 |
| other three | | none below E = 0.019 | |

**What is wrong biologically.** Two of five survive to the shortlist only because Pfam-A is silent and
MMseqs2 is silent; a five-minute search inside the project's own data names one of them as an
**entry-exclusion protein** — a well-defined plasmid function, and (see F5) exactly the function the
project hypothesised for a different protein without connecting the two. Note also that the Exc1
carrier plasmid carries `roM` (Rom/Rop) and `RNAI`: it is a ColE1-type replicon, so this is not a
coincidental hit.

**Concrete fix.** Run the F1 screen over the 90% survivors and the 30%/50% target sets *before* any
of them is proposed for folding or synthesis; annotate each with its nearest named in-corpus homolog
and E-value in the deliverable TSV. Candidates that survive that screen are worth a slot; the two
above are not.

---

### [CRITICAL] F3 — The dispersal paradox is circular, and the test that would break the circle has existed since July

**Evidence anchor.** `reports/small_cryptic_methodology.md` §3: *"The defensible statement is
**near-parity**: elements encoding no transfer machinery achieve essentially the same habitat and
geographic spread per lineage as conjugative plasmids"*; Limitation 3: *"'Non-mobilizable' is a
detection statement … the first thing the follow-up work should test (relaxed oriT search against
`orit_db_folder/`)"*. `orit_db_folder/` has been on disk, indexed, since **2026-07-10**.

Direct check of the master table:

| | n | mob_typer oriT | PlasAnn oriT feature |
|---|---:|---:|---:|
| non-mobilizable (analysis set) | 79,388 | **0 (0.0%)** | **3,070 (3.9%)** |
| mobilizable | 38,467 | 38,467 | 14,665 (38.1%) |
| small-cryptic, non-mobilizable | 39,711 | **0** | **1,592 (4.0%)** |

**What is wrong biologically.** mob_typer assigns *mobilizable* precisely when it finds a relaxase or
an oriT, so `mob_orit` is empty for the non-mobilizable class **by construction** — there is not one
byte of independent transfer evidence for the group whose transfer capability the paradox is about.
Meanwhile a second, entirely independent annotator already puts an oriT on 3,070 of them. And PlasAnn's
oriT sensitivity is only 38.1% on plasmids mob_typer calls mobilizable, so 3.9% is a **floor**;
scaling by that sensitivity gives a plausible true rate near 10%. Beyond oriTs, the project has ruled
out none of the standard alternatives for how a replicon with no *cis* transfer machinery travels:
**mobilization in trans** by a co-resident conjugative plasmid (the normal life of ColE1/Col-type
plasmids, which carry only an oriT plus a relaxase and borrow the T4SS), **natural transformation**
(the dominant route in *Acinetobacter*, *Neisseria*, *Streptococcus*, *Vibrio* — all present in this
dataset), **generalized and lateral transduction**, and **membrane-vesicle-mediated transfer**.
mob_typer's oriT database is small and relaxase-centric, so "non-mobilizable" is closest in meaning
to *"carries no relaxase MOB-suite recognises"*.

**Concrete fix.** Three tests, in order. (i) Run the blastn against `orit_db_folder/` that the report
already names, at a relaxed threshold, and add a fourth mobility class `oriT-only / mobilizable in
trans`; re-run §3 with that class separated out. (ii) Add `plasann_has_orit` as an independent
covariate to the dispersal model *now* — it is already in the master table, column 40. (iii) Test
trans-mobilization directly where the data allow it: for isolate-derived plasmids, ask whether
non-mobilizable small plasmids are over-represented in genomes that also carry a conjugative plasmid.
Until (i)–(iii) are done, the sentence must read *"elements in which MOB-suite detects no relaxase"*,
never *"elements encoding no transfer machinery"*.

---

### [MAJOR] F4 — Four "replicon" columns measure four different things, and the mob_typer-only rule inherits a version of the same defect

**Evidence anchor.** `PROJECT_OVERVIEW.md` §1: *"The central biological question driving the design is
the **multireplicon** question"*. VERIFICATION_PACK §7.2 correctly kills `pf_n_inc` (allele-level
inflation, `IncFII(29)`…`IncFII(pSFO)` counted as 11 replicons). The standing rule is now
`mob_rep_types` only, and the deck reports **24.7% of typed**.

The four columns are not commensurable, and this should be stated once, plainly:

- `pf_inc_types` — PlasmidFinder **allele** names (482 references, CGE 80/60 thresholds). A *locus*,
  matched many times.
- `pf_inc_families` — the same collapsed to Inc family. Computed by `scripts/harvest_typing.py` and
  then **not used for the count** (VERIFICATION_PACK §7.2).
- `mob_rep_types` — MOB-suite's own rep database: a **mixed nomenclature**.
- `plasann_replicons` — permissive 60% blastn, an annotation signal only.

Decomposing the 18,685 mob-based "multireplicon" plasmids:

| | n | % of multireplicon |
|---|---:|---:|
| every token is `IncF*` (IncFIA/FIB/FIC/FII) | **4,300** | **23.0** |
| ≥ 1 token is an unnamed `rep_cluster_NNNN` | **11,291** | **60.4** |
| *every* token is an unnamed `rep_cluster_NNNN` | 4,904 | 26.2 |

Commonest sets: `IncFIB+IncFII` (1,255), `IncFIB+IncFII+rep_cluster_2183` (1,002),
`IncFIA+IncFIC` (829), `IncFIA+IncFIB+IncFIC` (610). Distinct tokens overall: **1,759**.

**What is wrong biologically.** An IncF plasmid *is defined* by carrying RepFIA + RepFIB + RepFIC/FII
— that is what the FAB replicon-sequence-typing formula encodes (Villa et al. 2010). Counting
IncFIB+IncFII as "two distinct replicons" makes the canonical single IncF replication system the
archetype of multireplicon, which is circular: 23% of the multireplicon set is one well-known
replicon system, counted twice or three times. Symmetrically, `rep_cluster_NNNN` ids are MOB-suite
clusters of its rep database at an internal similarity threshold; their granularity is **not
calibrated to Inc families**, so for 26% of the "multireplicon" set nobody can say which two replicon
families are present, and two `rep_cluster` tokens may be two sub-variants of one replicon. The
ColE1 case is the mirror image and is visible in the project's own top families
(*"Col / ColRNAI / Col440I"*, `typing_methodology.md` §1): a single ColE1-type origin hits several Col
references simultaneously, and `GenBank_CP039323.1` in this dataset carries both `Col440II` and
`Col(pHAD28)` calls 800 bp apart — one replicon, two names. So the mob_typer-only rule is the right
*relative* call (it avoids the allele explosion) but it is not a clean one.

**Concrete fix.** Define "multireplicon" once, in replicon-family space, with the collapse rules
written down: (a) all `IncF*` → one `IncF` system unless a non-F replicon is also present; (b) all
`Col*` / ColE1-type → one `ColE1-type` system; (c) declare `rep_cluster_NNNN` **untyped**, not a
family, and report multireplicon prevalence separately for the named-only subset. Then report three
numbers side by side — PlasmidFinder families, mob_typer collapsed families, named-only subset — and
say which is the headline. Also: `pf_n_inc ≥ 2 = 24,583` must be removed or footnoted in
`PROJECT_OVERVIEW.md` §4 and `reports/typing_methodology.md` §§1,3, which still publish it without
qualification even though the deck (slide "2 · PlasmidFinder inflates multireplicon counts") already
knows better.

---

### [MAJOR] F5 — "Payload-free" is a property of four columns; the plasmids are not empty

**Evidence anchor.** `reports/small_cryptic_methodology.md`: *"**cryptic** — no AMR (`card_n_arg`),
no virulence, no metal/biocide, and no conjugation machinery"*; `defense_systems.ipynb`:
*"A plasmid carrying a complete Type II RM system is still 'payload-free' under the definition we
have been using."*

Census of what the 65,805 payload-free plasmids with gene calls actually encode
(`data/plasann_run/annot/*shard_*.tsv.gz`, CDS rows only, restricted to `cryptic_small_ids.txt`):

| PlasAnn category | CDS | plasmids | % |
|---|---:|---:|---:|
| Conjugation | 26,397 | 15,773 | **24.0** |
| Plasmid maintenance / replication / regulation | 18,007 | 14,742 | 22.4 |
| Other | 17,998 | 10,766 | 16.4 |
| Toxin–antitoxin | 5,790 | 3,819 | 5.8 |
| Non-conjugative DNA mobility | 5,657 | 3,573 | 5.4 |
| **Metabolism** | 2,431 | **1,480** | 2.2 |
| Stress response | 1,410 | 1,103 | 1.7 |
| **Antibiotic resistance** | 703 | **562** | 0.9 |

Plus **3,362 plasmids carrying 3,647 DefenseFinder systems** (`defense_systems.ipynb`).

Two things this shows that the reports do not say:

1. **The leak is inside PlasAnn's own categories, not only PlasAnn-vs-CARD.** The acknowledged
   "562-plasmid AMR definitional leak" is real (I reproduce **562** exactly from the gene shards —
   VERIFICATION_PACK §7.4's 694 is the master-column figure, so the shard number is the correct one
   for this set). But the `Metabolism` and `Other` buckets are never checked, and they contain
   genuine payload: `streptomycin adenylyltransferase Str` (89 CDS), `NAD(+)–rifampin
   ADP-ribosyltransferase` (Arr, 35), `chromate efflux transporter` + `chromate resistance protein`
   (70), `MccB`/`MccC` (microcin C7 biosynthesis, 94), a complete leucine-biosynthesis operon
   (isopropylmalate synthase/dehydratase/dehydrogenase, ~170), `catalase`, `glutathione reductase`,
   `TonB-dependent receptor`. None of these is counted by `plasann_n_amr`,
   `plasann_n_metal_biocide` or `plasann_n_virulence`.
2. **The 20 kb definition includes conjugation machinery on 24% of its members** — 15,773 plasmids.
   This is VERIFICATION_PACK §7.5 quantified, and it matters most for the dark proteome, every number
   of which is conditioned on this looser set.

**What is wrong biologically.** "Cryptic" in the plasmid literature means *no phenotype has been
demonstrated for this plasmid* — a statement about experiments, historically applied to a handful of
sequenced replicons. It has never meant "four annotation columns are zero". Beyond the categories
above, the payload list omits, in rough order of how often a referee will raise them:
restriction–modification and other defence (quantified: 3,647 systems), toxin–antitoxin (deliberately
excluded — defensible, but TA is the archetypal plasmid *addiction* payload and its exclusion must be
argued, not just declared), **partitioning and multimer resolution** (`par`, `cer`/Xer, `rom`/`rop`
copy-number control — `plasann_n_maintenance` exists as column 49 and is simply not used in the
filter), **entry/surface exclusion** (`excL`/Exc1 is filed under `Other`, so it is invisible),
bacteriocins outside the colicin/microcin naming (`pys2`, pyocin S2, 635 corpus-wide, sits in
`Other`), metal resistance outside the biocide list (`terF`, tellurite, 782, also `Other`), catabolic
and degradative pathways, nutrient acquisition (siderophores, TonB receptors), and
symbiosis/nodulation. One category error is also worth fixing: colicins and microcins (`ceA`, `cvaC`,
`cbA`, `ckA`, `mcmM`, `imM`, `cmI` …) *are* captured, but as **"Virulence and Defense Mechanism"**,
so a colicinogenic Col plasmid — the archetype of the small plasmid this project studies — is
excluded from the cryptic set as a *virulence* plasmid.

**Concrete fix.** Rename the compartment throughout to **"payload-undetected"** or
**"no annotated accessory cargo"**, and state the operational definition inline every time. Add
`plasann_n_maintenance`, a `Metabolism` count and a `Stress Response` count as columns and report
carriage rather than filtering on them. Run DefenseFinder once over the whole set and add a defence
column. And report, once, the fraction of the compartment that carries *something* under any of these
heads — on my census that is on the order of 5–6% before conjugation is even considered, which is a
finding in its own right and disarms the obvious referee objection.

---

### [MAJOR] F6 — The host-range inversion is measuring the predictor, not the host range

**Evidence anchor.** `reports/small_cryptic_methodology.md` §4: *"Multi-phylum calls are **9.3×** more
common in small cryptic than in conjugative plasmids"*, with the artefact control being the
mash-distance stratification.

**What is wrong biologically.** `mob_host_range` is MOB-suite's cluster-neighbourhood convergence
prediction: it reports the lowest common taxonomic rank of the *observed hosts of the plasmids in the
same MOB cluster*. Its resolution therefore scales with how densely that cluster is populated by
host-annotated members. Small cryptic plasmids are 73% IMG-PR-only metagenomic, are 71% untyped, and
sit in a sparse, poorly host-annotated cluster space; when a cluster contains three plasmids from
three unrelated metagenomes, "multi-phylum" is the *only* answer the algorithm can return. Conjugative
plasmids sit in dense, isolate-derived, *Enterobacterales*-saturated clusters, where the algorithm
converges to "order". A 9.3× ratio is close to what you would predict from that asymmetry alone. The
mash-distance control in §6d is not the right control: it holds neighbourhood *similarity* fixed but
not neighbourhood *size* or *host-annotation density*, which are the variables that drive the
prediction.

That said, the phenomenon is plausible and there is a real mechanism if it survives: small plasmids
are disproportionately **rolling-circle replicons** (Rep_1/Rep_2/Rep_trans/Mob_Pre families, ssDNA
intermediates) whose Rep proteins and MOBV/MOBQ relaxases are documented across Firmicutes,
Proteobacteria, Actinobacteria and Bacteroidota; broad-host-range is a known property of the pMV158
and pC194/pT181 families, and of ColE1-type mobilizable plasmids more generally.

**Concrete fix.** Replace the predictor with observation wherever possible. (i) Condition every
host-range comparison on cluster size and on the number of host-annotated members in the cluster,
and report the ratio within strata. (ii) Use the real host taxonomy already joined for the AMR work
(PlasmidScope `Host`, used in `reassessment_round1.md` §2) and compute *observed* host breadth per
MOB cluster; report the observed and predicted versions side by side. (iii) Split the small-cryptic
group by replication mode (RCR vs theta, from the Pfam rescue: `Rep_1`, `Rep_2`, `Rep_trans`,
`Mob_Pre` vs `Rep_3`/`RepA_N`) and test whether the multi-phylum excess is carried by the RCR half —
that is the mechanistic prediction, and it is testable with data already on disk.

---

### [MAJOR] F7 — The dark-plasmidome claim is a well-known gap being quantified for the first time, and is not framed that way

**Evidence anchor.** Deck: *"a discovery or a measurement of how bad plasmid annotation is? We think
it is both."* `reports/pfam_dark_validation.md` §5: *"a large, lineage-spanning protein space in the
small plasmidome is invisible to both plasmid annotation and Pfam."*

**What is wrong biologically.** To a plasmid biologist "most genes on small cryptic plasmids have no
assigned function" is not surprising — it is the field's standing complaint, stated in every
plasmid-genomics survey since the first Col plasmid sequences, and re-stated quantitatively for gut
metagenomes by PlasX/MobMess (Yu, Fogarty & Eren 2024). The genuinely valuable contributions here are
(a) the **magnitude at this scale**, (b) the **positive control** (`NAMED` at 76.7% vs `DARKREP` at
26.3%) which most such claims omit, and (c) the **decomposition** into ~29% annotation failure vs
~57% residue. Framed as novelty, the claim will draw the response "we knew that"; framed as the first
calibrated measurement of the gap, with the positive control front and centre, it is publishable.
The 57.4% residue also needs three deflations stated before a referee states them: F1 above (a fifth
of it has a named homolog after a proper search), **gene-calling artefacts** (PlasAnn's ORF calls are
inherited and unverified — a Prodigal/Pyrodigal cross-call is a one-afternoon control), and the
**phage-plasmid / ssDNA-element boundary** — small, ORF-dense, untypeable circular metagenomic
contigs are exactly where phage-plasmids and CRESS/Microviridae-like elements hide. The deck lists
this as future work; it belongs in the limitations of every current claim, because those elements
would contribute genuinely novel protein families that are *not* plasmid biology.

**Concrete fix.** Retitle the claim from novelty to calibration. Add the Prodigal cross-call and a
geNomad/phage-plasmid screen (geNomad is already cited in `papers_of_interest.txt` and produced the
IMG/PR calls). Report the residue after F1's screen.

---

### [MINOR] F8 — PlasAnn label noise runs in both directions, and one direction inflates the payload filter

`data/plasann_run/annot/` contains a 187 aa protein named `sugE` /
**"quaternary ammonium compound-resistance protein"**, categorised **Antibiotic Resistance**,
on 304 plasmids corpus-wide. It is 141/141 residues identical to the 177 aa dark ORF-A of §F5/§ColE1
below (see next section) — i.e. a mobilisation accessory protein being counted as a biocide-resistance
gene. Small in absolute terms, but it means the payload filter has false positives as well as false
negatives, and those *remove* plasmids from the cryptic compartment. **Fix:** spot-audit the 20 most
frequent products in the AMR and Metal/Biocide categories against CARD/AMRFinder before either count
is used as a filter; the concordance work in `card_amr_methodology.md` already provides the machinery.

---

### [MINOR] F9 — The surviving One Health claim is correct and biologically interesting, but "size" needs a mechanism

`manuscript/reassessment_round1.md` states it well: *"apparent One Health structure in plasmid
resistomes is substantially explained by host lineage and plasmid size, and roughly halves when either
is controlled — but a real compartment gradient and a real mobility effect both survive."* That is
defensible, correctly hedged (including the honest note that host coverage is 99.2% clinical vs 20.1%
wildlife), and genuinely useful — quantifying how much of an assumed pattern is compositional is a
contribution. Two refinements. First, "plasmid size" is not a mechanism; state the biology — a larger
plasmid has more accessory-gene slots, and large size is itself a consequence of carrying a conjugative
backbone (~35–60 kb of T4SS), so size is partly a *mediator* on the causal path from mobility to
cargo, not only a confounder. Adjusting for it therefore removes some real mobility effect, which is
why the conjugative OR moves 12.23 → 4.77. Say so explicitly rather than presenting the adjusted
number as simply the truer one. Second, decision §5 (cut the multireplicon section on 9-fold
differential typing coverage) is exactly right and should be applied to *every* cross-compartment
comparison conditioned on a typing call, including the host-range work in F6.

---

## The ColE1 cassette — what ORF-A and ORF-B are

This is the most actionable item in the review. `notebooks_my/zoom_two_dark_orfs.ipynb` describes the
object precisely and correctly:

> *"They flank the relaxase. The offset between them is exactly ±2 genes on 87% of the 336 carriers …
> and that gene is a mobilisation gene (`moB` or `mobA`) in 72% of cases. The modal architecture is
> `reP — ORF-B — moB — ORF-A — yibT` on a **Col156** replicon (333 of 336 carriers). These are
> textbook small ColE1-family mobilizable plasmids … Two of the five or six genes on the canonical
> ColE1-like backbone have no name."*

and then stops at composition: *"one membrane-associated protein and one basic soluble protein,
flanking a relaxase … would also fit entry exclusion or a small regulatory system."*

**This is a literature gap on the authors' side, not an annotation gap.** The canonical ColE1
mobilisation region has been defined since 1989 as **`mbeA`–`mbeB`–`mbeC`–`mbeD`** (Boyd, Archer &
Sherratt 1989), with reported products of **60, 19.5, 13 and 11 kDa**; MbeA is the MOBHEN/MOBP
relaxase (Varsaki et al. 2003) and MbeC is a **ribbon–helix–helix nicking accessory protein** that
binds oriT (Varsaki et al. 2009). The `mob`/`mbe` accessory nomenclature across plasmid families is
catalogued in Francia et al. 2004. Every observation in the notebook maps onto this:

| observation (notebook) | ColE1 expectation |
|---|---|
| ORF-B: **132 aa**, net charge **+6**, 20% K/R, no hydrophobic window, sits beside the relaxase | **MbeC / MobC** — RHH relaxosome accessory, ~13 kDa, basic, DNA-binding |
| ORF-A: **exactly 177 aa in 98% of 411 copies**, N-terminal charged coiled-coil + one C-terminal TM helix, immediately downstream of the relaxase | **MbeB / MobB** — reported at **19.5 kDa ≈ 177 aa**, membrane-associated |
| always ±2 genes apart, never adjacent, relaxase between them | the `mbe` operon |
| Col156 on 333/336, median 5,192 bp, all `mobilizable`, zero AMR, mostly *E. coli* | a ColE1-family mobilizable Col plasmid |
| *"almost invariant … under selection"* | core backbone, not cargo |

**Two independent lines of evidence I generated confirm the direction, and one of them is already
inside the project's own data.**

1. **Permissive Pfam scan** (`hmmscan -E 100`, no gathering thresholds) of
   `data/zoom_two_orfs/two_reps.faa`: ORF-B's best Pfam hit is **`TraD` (PF06412, "Conjugal transfer
   protein TraD")** at E = 0.047 — a relaxosome accessory family, and one already on the project's own
   list of top Pfam-rescued dark families (`pfam_dark_validation.md` §4.3). ORF-A's best is
   `DUF2730` (PF10805) at E = 0.037. Both far below significance, as expected for a family Pfam
   models only through divergent relatives — Pfam *does* carry `MobC` (PF05713), `MobC_2` (PF19514),
   `Mobilization_B` (PF17511), **`MbeD_MobD` (PF04899)**, `Mob_Pre` (PF01076), `Relaxase` (PF03432)
   and `Exclusion-determining_protein` (PF27689), so the right test is a targeted, threshold-free
   comparison against those models, not `--cut_ga` silence.

2. **`phmmer` against the project's own 78,397 named plasmid proteins** — decisive:

   | query | best named homolog | PlasAnn name on that plasmid | E |
   |---|---|---|---:|
   | ORF-A `GenBank_CP063716.1\|6` (155 aa) | `EMBL_OW967971.1\|1`, `GenBank_CP042891.1\|2`, `GenBank_CP048365.1\|4` | **`sugE` — "quaternary ammonium compound-resistance protein"** | **9.5e-76** |
   | ORF-B `IMGPR_plasmid_2609460067_000008\|1` (132 aa) | `IMGPR_plasmid_2734482010_000011\|4` | **`traD` — "conjugal transfer protein TraD (Potential coupling protein)"** | **1.7e-36** |

   ORF-B has **214 named hits** in total at E < 1e-3; ORF-A has 4. The ORF-A alignment is exact for
   **141 consecutive residues** (`MASTEEKLLTLLLSVEHLQHTAM…SWKWLGGVSVAFLVMFCVLGWGMKTM`), after which
   the "named" copy runs into a different C-terminus — i.e. a frameshifted or mis-called variant of
   the same gene that PlasAnn happened to label. The `sugE` label is itself wrong (SugE is a 105 aa
   four-TM SMR transporter; this is a 187 aa single-TM protein sitting immediately downstream of a
   real `sugE` in a guanidine-riboswitch/`gdx` module) — see F8.

**So the honest headline is stronger, not weaker, than the current one.** Not *"two conserved proteins
nobody can name"* but: *"two of the four genes of the canonical ColE1 mobilisation operon are
unnamed by a plasmid-specific annotator and by Pfam-A at its curated thresholds, on 336 plasmids
across 30 countries — and the same annotator names their orthologues elsewhere in the same dataset as
`traD` and as a biocide-resistance gene."* That is a sharper, more defensible and more publishable
claim, and it makes the project's annotation-failure thesis concrete on a named, textbook object.

**What to do, in order (about a day's work, no new compute):**

1. **Name them.** Pull ColE1 (GenBank `J01566`) MbeA/MbeB/MbeC/MbeD and the MOB reference relaxosome
   accessory proteins, and run `phmmer`/`jackhmmer` + HHpred/HHblits with ORF-A and ORF-B as queries.
   Also run both against Pfam PF05713/PF19514/PF17511/PF04899/PF06412 **without** `--cut_ga`. My
   prediction, stated so it can be falsified: **ORF-B = MbeC/MobC (RHH nicking accessory)** and
   **ORF-A = MbeB/MobB**, with `MbeD_MobD` the main alternative for ORF-A.
2. **Check the flanks.** `yibT` and `reP` on the same 336 carriers should resolve to the rest of the
   ColE1 backbone; look specifically for **`cer`/Xer recombination sites**, **`rom`/`rop`** (PlasAnn
   *does* name `roM` elsewhere in the corpus — see F2) and **RNAI/RNAII** copy-number control, all of
   which the pipeline detects as non-coding features and none of which appears in the cassette
   description. A Col156 plasmid that has `RNAI` but no `cer` and no `rom` is worth a sentence.
3. **Rewrite the deck slide and the notebook conclusion** around "canonical operon, unnamed by the
   tools", and keep the wet-lab proposal — the mobilisation-frequency knockout and the entry-exclusion
   assay are both exactly the right experiments for MbeB/MbeC/MbeD, and MbeD in particular has been
   linked to exclusion in ColE1-family plasmids. Naming the genes makes the experiment fundable
   instead of speculative.
4. **Keep the denominator discipline.** The 365-pairings-vs-336-carriers note is good practice and
   should survive the rewrite.

---

## Missing literature

`papers_of_interest.txt` holds 15 references, 11 marked `[?]` ("from memory"). The gaps are systematic:
the whole de la Cruz / Rocha mobilization corpus, all ColE1 molecular biology, and the structural
novel-family precedents. Named below with what each one changes.

**Mobilization and the dispersal paradox (F3) — the most consequential gap**
1. **Ares-Arroyo M, Coluzzi C, Rocha EPC. "Origins of transfer establish networks of functional
   dependencies for plasmid transfer by conjugation." *Nucleic Acids Research* 51:3001–3016 (2023).**
   The reference oriT-based method and database; shows a large fraction of plasmids without relaxases
   carry oriTs and are mobilizable in trans. *Changes:* makes the relaxed oriT search mandatory, and
   supplies the method and the comparison numbers for §3.
2. **Coluzzi C, Garcillán-Barcia MP, de la Cruz F, Rocha EPC. "Evolution of plasmid mobility: origin
   and fate of conjugative and non-conjugative plasmids." *Mol Biol Evol* 39:msac115 (2022).**
   *Changes:* "non-mobilizable" is an evolutionary state that plasmids enter and leave, not a class;
   directly reframes the paradox.
3. **Smillie C, Garcillán-Barcia MP, Francia MV, Rocha EPC, de la Cruz F. "Mobility of plasmids."
   *Microbiol Mol Biol Rev* 74:434–452 (2010).** Listed but marked `[?]` and unused. *Changes:* the
   canonical mobilizable/conjugative/non-mobilizable tripartition and its detection limits.
4. **Garcillán-Barcia MP, Francia MV, de la Cruz F. "The diversity of conjugative relaxases and its
   application in plasmid classification." *FEMS Microbiol Rev* 33:657–687 (2009).** *Changes:*
   MOB typing is the classification mob_typer implements; needed to say what `mob_relaxase` means.
5. **Ramsay JP, Firth N. "Diverse mobilization strategies facilitate transfer of non-conjugative
   mobile genetic elements." *Curr Opin Microbiol* 38:1–9 (2017).** *Changes:* enumerates the
   alternatives (trans-mobilization, transduction, transformation, vesicles) the project must rule out.
6. **Johnston C, Martin B, Fichant G, Polard P, Claverys J-P. "Bacterial transformation: distribution,
   shared mechanisms and divergent control." *Nat Rev Microbiol* 12:181–196 (2014).** *Changes:*
   supplies the second dispersal route for the *Acinetobacter*/*Streptococcus*/*Vibrio* fraction.
7. **Chen J, Quiles-Puchalt N, Chiang YN, et al. "Genome hypermobility by lateral transduction."
   *Science* 362:207–212 (2018).** *Changes:* the third route; relevant to gut/wastewater carriers.

**ColE1 and the cassette (dedicated section)**
8. **Boyd AC, Archer JAK, Sherratt DJ. "Characterization of the ColE1 mobilization region and its
   protein products." *Mol Gen Genet* 217:488–498 (1989).** *Changes:* names ORF-A and ORF-B.
9. **Varsaki A, Lucas M, Afendra AS, Drainas C, de la Cruz F. "Genetic and biochemical
   characterization of MbeA, the relaxase involved in plasmid ColE1 conjugative mobilization."
   *Mol Microbiol* 48:481–493 (2003).** *Changes:* identifies the middle gene of the cassette.
10. **Varsaki A, Moncalián G, Garcillán-Barcia MP, Drainas C, de la Cruz F. "Analysis of ColE1 MbeC
    unveils an extended ribbon-helix-helix family of nicking accessory proteins." *J Bacteriol*
    191:1446–1455 (2009).** *Changes:* the RHH-family assignment for ORF-B, and predicts its oriT
    binding — a testable structural claim.
11. **Francia MV, Varsaki A, Garcillán-Barcia MP, Latorre A, Drainas C, de la Cruz F. "A
    classification scheme for mobilization regions of bacterial plasmids." *FEMS Microbiol Rev*
    28:79–100 (2004).** *Changes:* the general MOB-region architecture the cassette instantiates.
12. **Summers DK, Sherratt DJ. "Multimerization of high copy number plasmids causes instability:
    ColE1 encodes a determinant essential for plasmid monomerization and stability." *Cell*
    36:1097–1103 (1984).** *Changes:* `cer`/Xer is a maintenance module the payload definition and
    the cassette description both miss.
13. **Cascales E, Buchanan SK, Duché D, et al. "Colicin biology." *Microbiol Mol Biol Rev*
    71:158–229 (2007).** *Changes:* Col plasmids are named for colicins; grounds the decision about
    where colicins/microcins belong in the payload taxonomy (currently "Virulence and Defense").

**Small-plasmid and replicon biology (F4, F6)**
14. **del Solar G, Giraldo R, Ruiz-Echevarría MJ, Espinosa M, Díaz-Orejas R. "Replication and control
    of circular bacterial plasmids." *Microbiol Mol Biol Rev* 62:434–464 (1998)** and **Khan SA.
    "Plasmid rolling-circle replication: highlights of two decades of research." *Plasmid*
    53:126–136 (2005).** *Changes:* the Rep families the Pfam rescue recovered (`Rep_1`, `Rep_3`,
    `RepL`, `Mob_Pre`) and the RCR/broad-host-range mechanism for F6.
15. **Villa L, García-Fernández A, Fortini D, Carattoli A. "Replicon sequence typing of IncF plasmids
    carrying virulence and resistance determinants." *J Antimicrob Chemother* 65:2518–2529 (2010)**
    and **Carattoli A, Bertini A, Villa L, et al. "Identification of plasmids by PCR-based replicon
    typing." *J Microbiol Methods* 63:219–228 (2005).** *Changes:* establishes that IncF plasmids
    carry FIA+FIB+FII **by definition** — this is the reference that makes F4's 23% figure a problem
    rather than a curiosity.
16. **Robertson J, Bessonov K, Schonfeld J, Nash JHE. "Universal whole-sequence-based plasmid typing
    and its utility to prediction of host range and epidemiological surveillance." *Microb Genom*
    6:mgen000435 (2020).** *Changes:* this — not the 2018 MOB-suite paper currently listed — is the
    source of `mob_cluster` and `mob_host_range`, including the authors' own stated accuracy and the
    sampling dependence that F6 turns on.
17. **Million-Weaver S, Camps M. "Mechanisms of plasmid segregation: have multicopy plasmids been
    overlooked?" *Plasmid* 75:27–36 (2014).** *Changes:* explains how a "module-free" 4 kb plasmid
    persists without par — directly answers §6's *"62.4% … no recognisable module of any kind"*.
18. **Garcillán-Barcia MP, Redondo-Salvo S, de la Cruz F. "Plasmid classifications." *Plasmid*
    126:102684 (2023).** *Changes:* the authoritative comparison of replicon/MOB/PTU/MOB-cluster
    schemes — exactly the "are these measuring the same thing?" question of F4.

**Dark matter, novelty and prior art (F1, F7)**
19. **Yu MK, Fogarty EC, Eren AM. "Diverse plasmid systems and their ecology across human gut
    metagenomes revealed by PlasX and MobMess." *Nat Microbiol* 9:830–847 (2024).** Cited once in the
    threshold report but absent from the reading list. *Changes:* the closest prior art for the dark
    plasmidome claim; the novelty framing of F7 must be written against it.
20. **Rodríguez del Río Á, Giner-Lamia J, Cantalapiedra CP, et al. "Functional and evolutionary
    significance of unknown genes from uncultivated taxa." *Nature* 626:377–384 (2024)** (FESNov).
    *Changes:* the standard for what evidence a "novel family" claim must carry.
21. **Durairaj J, Waterhouse AM, Mets T, et al. "Uncovering new families and folds in the natural
    protein universe." *Nature* 622:646–653 (2023)** and **Barrio-Hernandez I, Yeo J, Jänes J, et al.
    "Clustering predicted structures at the scale of the known protein universe." *Nature*
    622:637–645 (2023).** *Changes:* the structural-clustering precedent the project's own
    "foldseek next" plan needs, and the demonstration that sequence-family silence is a weak
    novelty criterion — the argument of F1 and F2.
22. **van Kempen M, Kim SS, Tumescheit C, et al. "Fast and accurate protein structure search with
    Foldseek." *Nat Biotechnol* 42:243–246 (2024).** *Changes:* the named tool for the step the
    project repeatedly says is next.
23. **Pfeifer E, Moura de Sousa JA, Touchon M, Rocha EPC. "Bacteria have numerous distinctive groups
    of phage–plasmids with conserved phage and variable plasmid gene repertoires." *Nucleic Acids
    Research* 49:2655–2673 (2021).** *Changes:* the phage-plasmid contamination limitation in F7
    goes from a deck footnote to a quantified screen.
24. **Tesson F, Hervé A, Mordret E, et al. "Systematic and quantitative view of the antiviral arsenal
    of prokaryotes." *Nat Commun* 13:2561 (2022).** DefenseFinder is used but uncited.
    *Changes:* required citation, and supplies the plasmid-vs-chromosome defence baseline that
    `defense_systems.ipynb` currently lacks.

**Also worth confirming rather than adding:** the 11 `[?]` entries in `papers_of_interest.txt` are
still unverified as of 2026-09-02, three months after the note was written. The de Quinto et al.
preprint that the whole design descends from should be re-checked for a peer-reviewed version before
`RESEARCH_OUTLINE.md` §1 is used in a manuscript.

---

## Scores

| axis | score | justification |
|---|---:|---|
| **Biological validity** | **56 / 100** | The measurement layer is sound and the arithmetic is verified. But four headline biological constructs — "multireplicon", "payload-free", "no named homolog", "non-mobilizable" — each mean something narrower than the text asserts, and I was able to break three of them empirically in under an hour using only data already on disk. The dispersal paradox (F3) is circular as stated; the host-range inversion (F6) is very likely an artefact of the predictor; the dark-only fraction (F1) is inflated by roughly a fifth. None of these is fatal — all are fixable, and the fixes make the claims sharper — but as currently written they would not survive a plasmid-literate referee. |
| **Literature engagement** | **34 / 100** | 15 references, 11 marked "from memory" and still unverified after three months. The de la Cruz/Rocha mobilization corpus, the ColE1 molecular genetics that would name the project's own flagship finding, the MOB-suite host-range paper that produced two of its columns, and the structural novel-family precedents are all absent. Two partial offsets: the reference table in `clustering_threshold_sensitivity.md` §2 is genuinely good comparative scholarship, and the AMR manuscript's 16 references were verified against primary sources. The problem is not carelessness, it is that the reading has not followed the project as it moved from multireplicon epidemiology into small-plasmid molecular biology. |
| **Domain contribution** | **61 / 100** | The resource is real: 208,248 complete plasmids with seven annotation layers, 100% mob_typer coverage, a locked habitat taxonomy, and a verified master table. The calibrated dark-matter decomposition (29% annotation failure / 57% residue, with a positive control) is a genuine and citable contribution. The 336-carrier ColE1 cassette is a real, well-characterised object. Deductions: the "novel family" deliverables shrink substantially under F1/F2, the multireplicon question the project was built around is not yet answered in a form a plasmid biologist would accept (F4), and the strongest single finding is currently mis-framed as an absence rather than a named gene. |
| **Interpretive calibration** | **73 / 100** | The project's best axis, and unusually good. `reassessment_round1.md` refutes its own title claim with matched-size predicted probabilities rather than defending it. The Pfam positive control, the dbAPIS adjudication that caught `acrIIA21` → RepA_N, the clean anti-defence negative, the threshold sweep that dissolves the project's own 93-family core down to 5, and the public documentation of two pipeline defects are all above the norm for this kind of work. Marked down because the discipline is applied unevenly: the same scepticism is not turned on the definitions themselves (F1, F3, F5), and specific numbers — 85.9%, 9.3×, "near-parity", "no name" — are stated flatly in the deck and reports after the caveat has been made elsewhere. |

---

## Genuine strengths

1. **The reassessment is a model of scientific honesty.** `manuscript/reassessment_round1.md` §1
   computes matched-size adjusted probabilities and concludes *"the sentence the paper is named after
   is an artefact of size composition."* Very few projects kill their own title claim, and fewer
   still do it with the right statistic.
2. **The Pfam design is correct.** Searching the 78,397 PlasAnn-**named** proteins as a positive
   control (76.7% vs 26.3%) is the one thing that makes a "dark matter" claim interpretable, and it is
   routinely omitted in published versions of this analysis.
3. **The negative results are properly negative.** Anti-defence is rejected at three stringency tiers
   and by two independent tools, *and* the adjudication that caught `acrIIA21` firing on RepA_N
   replication initiators is exactly the kind of check that separates a real screen from a hit list.
4. **The clonal-redundancy control is applied consistently and produces a non-obvious result** —
   defence carriage *falls* as MOB lineages grow (10.5% → 3.0%), the opposite of what clonal inflation
   would give. That is a real ecological observation.
5. **`clustering_threshold_sensitivity.md` is the best document in the repository.** It benchmarks
   against UniRef, FESNov, NMPFamsDB and PlasX, finds that the project's own clustering was only 89%
   compliant with its stated criterion, and then reports honestly that its 93-family core collapses to
   5. Publishing that took discipline.
6. **Denominator discipline.** The 365-pairings-vs-336-carriers warning in `zoom_two_dark_orfs.ipynb`,
   and the `*shard_*` globbing pitfall recorded for reuse, are the marks of someone who has been
   burned and wrote it down.
7. **Defects are documented, not silently patched** — the CDS-indexing bug that broke the dark↔all-CDS
   join for 10.7% of ORFs is written up with the corrupted and corrected cross-tabs side by side, and
   an `assert` added so it cannot recur.

---

## Three highest-value biological next steps, ranked

**1 · Name the ColE1 cassette (hours, no new compute, highest payoff).**
Run ORF-A and ORF-B against ColE1 `J01566` MbeA/B/C/D, against Pfam PF05713 / PF19514 / PF17511 /
PF04899 / PF06412 *without* `--cut_ga`, and against the project's own named proteome (which already
returns `traD` at E = 1.7e-36 and a mislabelled `sugE` at E = 9.5e-76). Then check the 336 carriers
for `cer`, `rom`/`rop` and RNAI. This converts the project's flagship object from "two anonymous
conserved proteins" into "the canonical ColE1 mobilisation operon is unnamed by PlasAnn and by Pfam-A
at its gathering thresholds on 336 plasmids in 30 countries" — a sharper claim, a better paper, and
a fundable knockout experiment (mobilisation frequency ± helper; entry exclusion; serial-passage
stability) because the genes now have names and priors.

**2 · Re-derive "novel" with a sensitive homology search, then structure (days).**
Replace MMseqs2 silence with `phmmer`/`jackhmmer` against `named.faa` plus Pfam without `--cut_ga`,
recompute the dark-only fraction (my sample says ~66%, not 85.9%), and re-filter
`dark_novel_families.tsv` and both core sets. Then run **Foldseek against the AlphaFold DB** on what
survives — the step the project names as next in four separate documents and has not taken. Add the
two cheap controls that would otherwise be a referee's first two questions: an independent
Prodigal/Pyrodigal gene call, and a geNomad phage-plasmid screen of the small untypeable metagenomic
fraction. Expect the residue to shrink again; a smaller number that survives all four screens is worth
far more than 57.4%.

**3 · Break the dispersal circularity (days).**
Run the relaxed oriT search against `orit_db_folder/` that has been ready since July, add relaxase
HMM profiles (MOBscan / MOBF-MOBP-MOBQ-MOBV-MOBC), and introduce a fourth mobility class —
*oriT-only / mobilizable in trans*. Immediately, and for free, add `plasann_has_orit` (column 40,
positive on 3,070 "non-mobilizable" plasmids) as a covariate to §3. Then test trans-mobilization
where the data allow: are non-mobilizable small plasmids over-represented in isolate genomes that also
carry a conjugative plasmid? If near-parity survives all of that, it is a genuinely interesting
result about transformation, transduction and vesicles; if it does not, the project has avoided
publishing a detection artefact as a biological paradox.
