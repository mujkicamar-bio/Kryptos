# Design — a controlled plasmid annotation pipeline for dark-ORF discovery

**Status:** design, awaiting approval
**Date:** 2026-09-02
**Supersedes:** the ad-hoc dark-ORF pipeline documented in `reports/dark_orf_clustering.md`
**Review that motivated it:** `review/EDITORIAL_DECISION.md` and the five reports under `review/reports/`

---

## 1 · Purpose

Build a reusable, controlled pipeline that annotates every ORF on every plasmid in the working set,
and that emits, as its primary product, **fully annotated plasmids**. The dark proteome is then
defined as what the annotation could not resolve — an output of the pipeline rather than an input
assumption — and is mined for lab-screenable targets.

### Goals

1. **Fully annotated plasmids.** Per-plasmid GFF3 / GenBank / TSV with every CDS carrying its
   coordinates, annotation, the evidence that produced it, and the tier at which it was resolved.
2. **A defensible dark set.** Proteins that survive a graded, pre-registered homology cascade.
3. **Class-resolved biology.** Which protein families are native to big, small, and cryptic
   plasmids — answered by stratifying one table, not by running separate studies.
4. **Screenable targets** in three tiers: pooled discovery, arrayed follow-up, cloning assays.
5. **Better annotation than existing tools**, measurably — the pipeline's own annotation rate against
   PlasAnn's on the same input is a reportable result.

### Non-goals

- Re-deriving the working set. The 208,248-plasmid assembly and the locked habitat taxonomy are
  reused unchanged (§3).
- Re-running geNomad. Explicitly excluded (§5.2).
- Rebuilding the AMR, PLSDB or environment layers. They verified clean.
- Writing any paper. Numbers will move; writing waits.

---

## 2 · Design principles

These are the five inversions that separate this pipeline from its predecessor. Each maps to a
finding in `review/EDITORIAL_DECISION.md`.

| # | Principle | Defect it removes |
|---|---|---|
| **P1** | **Darkness is an output, not an input.** No tool's silence defines a category. A protein is dark because a recorded cascade failed, at thresholds fixed before the run. | The whole "is PlasAnn's `ORF` label real?" problem, and the 28.5%-artifact controversy |
| **P2** | **Every sequence is annotated individually.** No label is ever inferred from a cluster representative. | PER F6 (representative-only search labelling 378,552 proteins); the 96.1%-not-100% propagation error |
| **P3** | **Homology depth and functional resolution are orthogonal.** Hitting a "hypothetical protein" in NR is a homolog, not an annotation. The two are recorded separately. | The core annotation-quality failure that no standard pipeline handles (§6) |
| **P4** | **Clonal structure is controlled once, centrally.** Not per notebook, not optionally. | DA F1 (the control existed and was absent from every dark notebook) |
| **P5** | **Thresholds are pre-registered in config.** A tier resolves a protein only against a threshold declared before the run and recorded in the output row. | MET D1 (an undocumented 80%-vs-95% choice moved the founding number by 65%) |

---

## 3 · Reused vs rebuilt

Reuse is deliberate and evidence-based: these layers were independently recomputed
(`review/VERIFICATION_PACK.md`, 47/52 exact) and survived adversarial testing.

| layer | decision | why |
|---|---|---|
| Working set (208,248) + FASTA | **reuse** | Recomputed exactly; dedup byte-identical across all 147 duplicated plasmids |
| Analysis set (143,503, artifacts excluded) | **reuse** | Locked decision, applied consistently; verified |
| Habitat taxonomy, geography, `is_clinical` | **reuse** | Locked before analysis; verified |
| CARD/RGI AMR layer | **reuse** | RGI covered 208,245/208,248; zero-fill benign (MET C9) |
| PLSDB enrichment, MOB-suite mobility/cluster | **reuse** | 100% coverage; joins 1:1 (MET C11) |
| Element provenance (`sources`, IMG-PR flag) | **reuse** | Needed for the provenance stratification that defeated the artifact attack |
| **Gene calls (PlasAnn)** | **rebuild** | The one inherited, unverified dependency under every dark claim; and coordinates were discarded (MET C6) |
| **All protein annotation** | **rebuild** | Single-database novelty; no cascade; representative propagation |
| **Dark-ORF clustering + shortlists** | **rebuild** | No clonal control; stale representatives (PER F2); threshold-dependent |
| `pf_n_inc`, `card_multidrug`, `plsdb_n_inc` | **fix in place** | Label-counting defect; fix at source, do not carry forward (§10) |

**Consequence, stated plainly:** Pyrodigal will not reproduce PlasAnn's ORF set. Protein counts,
boundaries and family sizes will all differ. The existing reports become **historical**, not
superseded — they describe a different object. This is a one-way door and is accepted.

---

## 4 · Stage graph

```
S0   INPUT              143,503 non-artifact plasmids + FASTA + reused metadata
      │
S1   ORF CALLING        Pyrodigal (meta mode), all plasmids
      │                 → proteins.faa + orfs.gff3 (coordinates, strand, partial flags)
      │                 → ~8–12M ORFs  [pilot-measured]
      │
S2   DEREPLICATION      exact-identity collapse (SHA-256 of sequence)   LOSSLESS
      │                 → unique_proteins.faa (~4–6M) + protein↔ORF map
      │
S3   ANNOTATION CASCADE per unique sequence; tier N input = tier N−1 unresolved
      │                 T1 Pfam-A 38.2 --cut_ga
      │                 T2 Pfam sub-GA + composition-shuffled null
      │                 T3 VOGDB + PHROGs                  (phage proteins)
      │                 T4 UniRef90 → UniRef50 (MMseqs2)
      │                 T5 eggNOG-mapper                   (orthology + function)
      │                 T6 NR + Swiss-Prot (DIAMOND)
      │                 T7 HHblits vs PDB70                (remote / structural profile)
      │                 T8 ESMFold + Foldseek vs AFDB      (survivors only)
      │
S4   ►► ANNOTATED PLASMIDS ◄◄   per-plasmid GFF3 + GenBank + TSV     PRIMARY DELIVERABLE
      │                 every CDS: coords, label, tier, evidence, confidence, functional class
      │
S5   DARK SET           CDS with functional_class ∈ {UNCHARACTERIZED_HOMOLOG, NONE}
      │                 + invariant check: co-cluster dark with annotated, assert no leakage
      │
S6   CLUSTERING         MMseqs2 --cluster-reassign, identity sweep, DARK ONLY
      │
S7   STRATIFICATION     size / payload / mobility / element_class / host / habitat / provenance
      │                 → per-class family enrichment: what is native to which plasmid class
      │
S8   SYSTEMS            operons, gene neighbourhoods, BGC (antiSMASH/BAGEL), defence, SignalP
      │                 uses S1 coordinates; prophage-boundary aware
      │
S9   TARGETS            clonal-controlled ranking → Tier 1 pooled / Tier 2 arrayed / Tier 3 cloning
```

---

## 5 · Stage specifications

### S1 · ORF calling

| | |
|---|---|
| Tool | Pyrodigal ≥3.x (container; `prodigal/2.6.3` module as cross-check) |
| Mode | `meta` — plasmids are short and taxonomically heterogeneous; single-genome training is unsafe below ~20 kb |
| Input | `working_set.fna.gz`, subset to the 143,503 analysis set |
| Output | `proteins.faa`, `orfs.gff3`, `orf_index.tsv` |
| Sharding | 600 shards, SLURM array, `pelle` |

**`orf_index.tsv` schema** — one row per ORF, and the join key for everything downstream:

```
orf_id  plasmid_id  contig  start  end  strand  length_aa  partial  rbs_motif
start_type  gc_cont  seq_sha256
```

`orf_id` = `{plasmid_id}|{ordinal}` where ordinal counts **every** ORF on the plasmid in coordinate
order, regardless of any later filtering. (This is the direct fix for the 2026-08-28 CDS-indexing
defect, where indices were assigned within a filtered subset and silently diverged between files.)

**Invariants (asserted, not assumed):**
- every `orf_id` unique; ordinals contiguous per plasmid
- `end > start`; `length_aa == (end-start+1)/3 - 1` for complete ORFs
- ORF count per plasmid > 0 for all plasmids ≥ 1 kb, else flagged for review
- protein FASTA record count == `orf_index` row count

**Cross-check gate:** run Prodigal 2.6.3 on a 1,000-plasmid sample and require ≥95% boundary
concordance with Pyrodigal. This is the control the previous pipeline never had.

### S2 · Dereplication

Exact-sequence collapse only — **not** clustering. Identical sequences receive identical annotations
by definition, so this is lossless and is purely a compute optimisation.

```
unique_proteins.faa           one record per distinct sequence, id = seq_sha256[:16]
protein_map.tsv               seq_sha256 → [orf_id, ...]
```

Expected reduction is large: the devil's advocate measured 29.9% of dark-ORF carriers sharing an
exact proteome with another plasmid. **Measure it in the pilot; do not assume it.**

**Invariant:** `sum(len(orf_ids) for each unique)` == total ORF count. Asserted after every run.

### S3 · Annotation cascade

Detailed in §6. Structural rules:

- Each tier reads the previous tier's **unresolved** set. No tier ever sees a protein another tier
  resolved, so cost falls monotonically.
- Every tier writes one row per query it resolves, including the threshold it applied.
- Tier order is fixed in config. Reordering tiers is a config change that invalidates the DAG
  downstream — by design.
- T8 (ESMFold) runs on the `gpu` partition (`l40s`), all other tiers on `pelle`.

### S4 · Annotated plasmids — the primary deliverable

Per plasmid, three files generated from one source of truth:

```
annotated/{plasmid_id}.gff3       CDS features + annotation attributes
annotated/{plasmid_id}.gbk        GenBank flat file, for tooling that expects it
annotated/plasmid_annotation.tsv  one row per CDS, all plasmids (the analysis table)
```

**`plasmid_annotation.tsv` schema:**

```
orf_id  plasmid_id  start  end  strand  length_aa  partial
annot_label  annot_source_db  annot_accession
annot_tier  annot_evidence  annot_confidence
functional_class  homology_depth
neighbourhood_id  operon_id
```

This table is the join target for S5–S9 and is the object other people will reuse.

### S5 · Dark set definition + leakage invariant

Dark = `functional_class ∈ {UNCHARACTERIZED_HOMOLOG, NONE}` (§6.2).

**`DOMAIN_ONLY` is deliberately *not* dark, but *is* target-eligible.** A protein carrying a bare
HTH or a domain covering < 50% of its length has a structural hint and no functional assignment. It
does not belong in the dark headline — that would re-import the over-claiming this design exists to
remove — but excluding it from screening would discard real candidates. It is therefore reported as
its own class, and enters S9 ranked below `UNCHARACTERIZED_HOMOLOG`. Every headline dark count is
reported with and without it, so the choice is visible rather than buried.

**Leakage invariant.** Cluster the dark set together with the annotated set once, at 30%/80%, purely
as a check. Assert that no dark protein lands in a family containing a protein annotated
`FUNCTIONAL` at confidence HIGH. Any violation means the cascade has a gap; the run fails and the
gap is investigated before targets are built. This converts a former *finding* (14.1% of dark ORFs
sat in families with a named member) into an automated gate.

### S6 · Clustering the dark proteome

```
mmseqs easy-cluster --min-seq-id {0.3,0.5,0.7,0.9} -c 0.8 --cov-mode 0 --cluster-reassign
```

`--cluster-reassign` is **mandatory and non-configurable**. Its absence in the previous run meant
only 89.12% of members met the criterion the methods claimed.

Coverage is swept as well as identity — `-c {0.5 --cov-mode 1, 0.8 --cov-mode 0}` — because the
perspective reviewer showed coverage is the dominant axis (−22% of families under FESNov's rule) and
that the clusterings are not nested. Two axes, not one.

**Compliance audit is a pipeline rule, not an afterthought:** sample 200 families per threshold,
align every member to its own representative, assert ≥99% compliance. The audit that found the
original defect becomes a gate that prevents it.

### S7 · Stratification and class-native families

Every plasmid carries its class columns; no analysis pre-filters on them.

```
size_class        small (<10kb) | mid (10–20kb) | large (≥20kb)     [thresholds in config]
payload_class     cryptic | cargo    (payload definition in config, versioned)
mobility_class    conjugative | mobilizable | non-mobilizable
element_class     plasmid | phage-plasmid-like | ambiguous          (§5.2)
provenance_class  IMG-PR-only | isolate-derived | both
```

**Class-native family analysis.** For each family and each class, compute observed vs expected
membership under a null that preserves family size and class sizes, with the **clonal control
applied first** (§8). Report enrichment with a multiple-testing correction across families. This is
what answers "how many proteins are native to each class", and it is a query on one table.

### S8 · Systems and context

Enabled by S1 coordinates, which the previous pipeline discarded.

- **Neighbourhoods:** ±5 genes around each dark ORF, with strand and intergenic distance.
- **Operon calling:** same-strand runs with intergenic gaps below a configured threshold.
- **BGC / small molecules:** antiSMASH and BAGEL (containers) for bacteriocins, microcins and RiPPs.
  Motivated directly — the domain reviewer found MccB/MccC (microcin biosynthesis) on "payload-free"
  plasmids, filed by PlasAnn under `Metabolism`/`Other`. Small molecules are a real seam here.
- **Defence systems:** DefenseFinder / PADLOC (containers).
- **Secretion:** SignalP 6.0 (module) — secreted small proteins are prime bacteriocin candidates and
  prime screening targets.
- **Prophage boundaries:** operons are never called across a putative provirus junction.

### S9 · Target selection

See §9.

---

## 5.2 · Phage handling (no geNomad)

Per-element phage classification is derived **from the annotation itself**, not from a separate
classifier:

```
element_class = phage-plasmid-like
    IF the plasmid carries annotated hallmark genes from ≥ N of:
       {large terminase, portal, major capsid, tail, holin/endolysin, phage integrase}
    N declared in config/classes.yaml (default 3), per principle P5 — not hardcoded
```

Rationale: (a) it avoids re-running the classifier that produced IMG/PR's calls in the first place —
62% of the set — which would largely be the tool agreeing with itself; (b) the rule is transparent
and auditable, unlike a model score; (c) it costs nothing, because the annotation already exists.

Per-protein phage annotation is tier T3 (VOGDB, ready as an MMseqs2 database; PHROGs by container).

**Phage-plasmid elements are flagged and kept, never excluded.** Anti-defence, host takeover and
superinfection exclusion live there, and the project's database-based anti-defence screen was
negative — which is exactly the case where a wet screen is worth more than another database.
`element_class` is reported on every S9 target so it is always known whether a target is plasmid
biology or phage biology.

---

## 6 · The annotation model

This is the intellectual core and the answer to "fixing the annotation problems other tools have".

### 6.1 · Two orthogonal axes

Existing pipelines record one bit: named or not. That bit conflates two different facts and is why
"hypothetical protein" pollutes every large annotation effort.

**Axis 1 — `homology_depth`:** how hard was it to find *any* homolog.

| tier | method | what a hit means |
|---|---|---|
| T1 | Pfam-A `--cut_ga` | curated family membership |
| T2 | Pfam sub-GA vs shuffled null | real but sub-threshold domain similarity |
| T3 | VOGDB / PHROGs | viral orthologous group |
| T4 | UniRef90 → UniRef50 | close then mid-range sequence homology |
| T5 | eggNOG-mapper | orthologous group with functional transfer |
| T6 | NR + Swiss-Prot (DIAMOND) | exhaustive sequence homology |
| T7 | HHblits vs PDB70 | remote / structural profile homology |
| T8 | ESMFold + Foldseek vs AFDB | structural homology |
| — | UNRESOLVED | nothing found anywhere |

**Axis 2 — `functional_class`:** does the best hit actually carry a function.

| value | meaning |
|---|---|
| `FUNCTIONAL` | best hit names a specific function |
| `DOMAIN_ONLY` | a domain matched but covers < 50% of the query, or is a bare structural motif (e.g. HTH) with no functional assignment |
| `UNCHARACTERIZED_HOMOLOG` | homologs exist, but every one is itself unnamed — "hypothetical protein", "DUF*", "uncharacterized", "putative uncharacterized" |
| `NONE` | no homolog at any tier |

### 6.2 · The class that matters most

The 2×2 gives four biologically distinct outcomes, and the previous pipeline could not express two
of them:

```
                    FUNCTIONAL          UNCHARACTERIZED_HOMOLOG      NONE
  shallow (T1-T3)   annotated           ← conserved & unknown        —
  deep    (T4-T8)   annotated late      ← conserved & unknown        —
  none                  —                        —                 orphan
```

**`UNCHARACTERIZED_HOMOLOG` is the prize class.** These proteins are *certainly real* (they have
homologs, often many, across taxa) and *certainly unknown* (nobody has ever characterised any of
them). For high-throughput lab screening that combination is worth more than an orphan with no
homolog at all, which is more likely to be lineage-specific or a calling artefact.

This distinction is the single most important thing the pipeline adds over PlasAnn, Pfam-alone, or a
standard Prokka/Bakta run — none of which separate "unknown" from "nobody has looked".

### 6.3 · Anti-gaming rules

Because `functional_class` decides what is dark, it must not be promotable by a weak hit.

1. **Pre-registered thresholds.** Every tier's threshold lives in `config/cascade.yaml`, is fixed
   before the run, and is written into every output row. Changing one invalidates downstream DAG
   outputs automatically.
2. **The uninformative-label list is versioned and recall-tested.** `hypothetical`, `DUF####`,
   `uncharacterized`, `putative uncharacterized`, `conserved protein`, `ORF`, `protein of unknown
   function`. Methodology finding C4 showed the previous regex had good precision and poor recall
   (`RelE` missing while `RelB` present; `TrfA`, `RepC`, `TrwC`, `Rop` filed as novel), which moved a
   published percentage by ≥7 points. This list therefore ships with a labelled test set and a
   measured recall figure, and CI fails if recall drops.
3. **Coverage floor.** A hit covering < 50% of the query is `DOMAIN_ONLY`, never `FUNCTIONAL`.
   Median Pfam alignment covered ~52% of the dark proteins in the previous run, so this matters.
4. **T2 requires its null.** Sub-GA Pfam hits count only if they beat a composition-shuffled control
   at the same threshold — otherwise low-complexity sequence promotes itself.
5. **Confidence is recorded, not inferred.** `annot_confidence ∈ {HIGH, MEDIUM, LOW}` from
   tier-specific evidence bands declared in config.
6. **No tier may overwrite a shallower tier's `FUNCTIONAL` call.** Cascade order is authority order.

---

## 7 · Workflow engine and control

**Snakemake 8.27** (`snakemake/8.27.0-foss-2024a`), SLURM executor plugin, apptainer for tools with
no module.

```
plasmid-annotation/
├── Snakefile
├── config/
│   ├── config.yaml            paths, partitions, shard counts
│   ├── cascade.yaml           tier order, thresholds, evidence bands   ← pre-registered
│   └── classes.yaml           size/payload/mobility/element definitions
├── workflow/rules/            s1_orfcall.smk … s9_targets.smk
├── workflow/envs/             one conda env or container per rule
├── workflow/scripts/          thin; logic lives in tested modules
├── src/plasmidann/            importable, unit-tested package
└── tests/                     unit + invariant + one 500-plasmid integration fixture
```

**How each historical defect becomes structurally impossible:**

| defect | structural prevention |
|---|---|
| Stale input silently consumed (the §R2 table with representatives that no longer exist) | Snakemake rebuilds any output whose input changed; a stale file cannot be read |
| Parameter change not propagated (`--cluster-reassign`) | Threshold lives in `cascade.yaml`; changing it invalidates every downstream output |
| Silent partial failure (one shard of 594) | Rule-level success criteria; the DAG does not complete |
| CDS index divergence between files | One `orf_id` assigned once in S1, asserted unique, never recomputed |
| Numbers existing only as hand-typed prose | Every reported number comes from a rule output file; notebooks read files, never retype |
| Environment not reproducible | Per-rule conda/container, pinned in the repo |
| No provenance | `snakemake --report` emits parameters, versions and runtimes for every rule |

**Version control.** The repository is currently not under git at all
(`review/reports/workflow_audit.md` W1). `git init` is a prerequisite for this design, not a
follow-up: without it, none of the above provenance claims are verifiable.

---

## 8 · Clonal control (applied once, centrally)

The devil's advocate established that the project's own dereplication control — one plasmid per MOB
cluster — was used in `cryptic_plasmids.ipynb` and `defense_systems.ipynb` but appeared in **no**
dark-family notebook, and that applying it moved recurrence from 84.7% to 39.2% and the substantial
family count from 3,343 to ~40.

**Standard, non-optional:** `cap_k = 10` plasmids per `mob_cluster`, applied at S7 before any
enrichment or ranking, seeded and repeated 20×.

Rationale for cap-10 rather than 1-per-cluster: full dereplication is too harsh because `mob_cluster`
AA379 alone absorbs ~26% of the small-plasmid subset and spans 236–19,992 bp — it is a junk drawer,
not a lineage. Cap-10 retained ~75.7% recurrence and ~1,040 substantial families in the DA's test.

**Every reported count is emitted at three levels**, always together, never one alone:

```
pooled          no control          (comparable to the old numbers)
cap-10          standard            ← the headline
1-per-cluster   maximal control     (lower bound)
```

A sampling null (observed vs expected lineage spread for a random draw of the same carrier count)
is computed for every family, following the perspective reviewer's method, and reported alongside
raw spread. Raw lineage spread alone is never reported — it was the statistic that inverted.

---

## 9 · Target selection — three tiers

All tiers draw from one ranked table; they differ only in filter profile. No assay-specific logic is
hardcoded, per the decision to keep assay options open.

**Shared ranking features** (each a column, weights in config):

```
functional_class == UNCHARACTERIZED_HOMOLOG     (the prize class, §6.2)
homology_depth                                  deeper = more thoroughly unknown
family size after cap-10 clonal control         real prevalence, not deposits
lineage spread vs sampling null (obs/exp)       conservation beyond chance
habitat / host breadth
length_aa                                       synthesis cost
partial flag                                    complete ORFs only for constructs
SignalP / TM prediction                         secreted or membrane
neighbourhood context                           near defence, BGC, or transfer loci
element_class                                   plasmid vs phage-plasmid biology
predicted structure + pLDDT (T8)                foldable and confident
Foldseek best hit                               structural analogy without sequence homology
```

| tier | n | filter emphasis | output |
|---|---|---|---|
| **1 — pooled discovery** | 1,000–10,000 | breadth across families; length ceiling for synthesis; complete ORFs; no lethality-predicting features unless lethality is the readout | synthesis-ready construct table + barcode assignments |
| **2 — arrayed follow-up** | 100–1,000 | one representative per high-ranked family; unambiguous per-well phenotype; host compatible with expression strain | per-well plate map |
| **3 — cloning assays** | 10–50 | highest confidence × highest novelty; structural prediction available; genomic context interpretable | full dossier per target |

**Every target ships with a dossier**: sequence, coordinates, source plasmid, host, habitat, family
membership at all four identity thresholds, cascade result at every tier, clonal-controlled
prevalence, obs/exp spread, structure and pLDDT, neighbours, and **the reason it was not annotated**.
That last field is what makes a negative screen interpretable.

---

## 10 · Resource-layer fixes folded in

These are not part of the pipeline proper but must land in the same pass, because S7 stratification
consumes them (`review/EDITORIAL_DECISION.md` CF2):

- `harvest_typing.py`: emit `pf_n_fam` (already computed, currently discarded); publish multireplicon
  as a range across the identity threshold (8,169–24,583), with the 95%-CGE value primary.
- `aggregate_card.py`: `card_multidrug` must count distinct AROs, not drug-class label strings
  (4,402 plasmids — 21.2% — are currently mislabelled).
- `master_table_data_dictionary.md:129`: corrected — it documents `pf_n_inc` as a family count when
  it is an allele count.
- `mob_rep_types`: exclude uncalibrated `rep_cluster_NNNN` tokens from any replicon count (60.4% of
  mob multireplicon calls contain one).

---

## 11 · Compute

All figures are **estimates to be replaced by pilot measurements**. Cluster: `pelle` (80×96 cores,
772 GB), `gpu` (l40s ×10, h100 ×2), `fat` (2.3 TB).

| stage | scale | estimate | partition |
|---|---|---|---|
| S1 Pyrodigal | 143,503 plasmids | 1–2 h, 600-way array | pelle |
| S2 dereplication | ~10M → ~5M | < 1 h | pelle |
| T1 Pfam `--cut_ga` | ~5M unique | 4–8 h, wide array | pelle |
| T2 sub-GA + null | unresolved | 6–12 h | pelle |
| T3 VOGDB/PHROGs | unresolved | 1–2 h | pelle |
| T4 UniRef90/50 | unresolved | 4–8 h | pelle / fat |
| T5 eggNOG-mapper | unresolved | 8–16 h | pelle |
| T6 NR (DIAMOND) | unresolved | 12–24 h | fat |
| T7 HHblits vs PDB70 | 10⁴–10⁵ | 12–24 h | pelle |
| T8 ESMFold + Foldseek | 10³–10⁴ | 6–24 GPU-h | gpu (l40s) |
| S6 clustering sweep | dark set only | 1–2 h | fat |

The cascade is self-narrowing, so the expensive tiers never see the full set. **T6 is the cost
driver and its query set size is the number to watch in the pilot.**

---

## 12 · Milestones and gates

Each milestone has a verification that must pass before the next begins.

| # | Milestone | Gate |
|---|---|---|
| **M0** | `git init`, repo skeleton, Snakemake scaffold, config schemas | `snakemake --lint` clean; CI runs unit tests |
| **M1** | S1+S2 on a 5,000-plasmid pilot | All S1/S2 invariants pass; Prodigal cross-check ≥95% boundary concordance; **measured** dereplication ratio and ORF count |
| **M2** | Full cascade on the pilot | Every tier emits rows with thresholds recorded; `annotation_depth` distribution sane; uninformative-label recall measured on the labelled test set |
| **M3** | Full-scale S1–S4: **annotated plasmids delivered** | Per-plasmid GFF3/GBK/TSV for all 143,503; annotation rate reported **against PlasAnn's on the same input** — the "we annotate better" result |
| **M4** | S5+S6: dark set + clustering | Leakage invariant passes; compliance audit ≥99% at every threshold |
| **M5** | S7+S8: stratification, class-native families, systems | Clonal control applied; all counts reported at three levels; sampling null computed |
| **M6** | S9: three target tiers with dossiers | Every target carries its full dossier including why it was not annotated |

**M3 is the point at which the project has a deliverable regardless of what the dark analysis
finds.** That is deliberate: it de-risks the programme, because "we annotated 143,503 plasmids
better than the existing tools, and here is the resource" stands on its own.

---

## 12.1 · Scope note — this spec is two implementation plans

M0–M3 and M4–M6 are separable projects with different risk profiles, and should get separate
implementation plans rather than one:

**Plan A — the annotation engine (M0–M3).** Infrastructure plus a standalone resource. Its success
does not depend on the dark analysis finding anything. Deliverable: 143,503 fully annotated plasmids
and a head-to-head annotation-rate comparison against PlasAnn.

**Plan B — dark analysis and targets (M4–M6).** Consumes Plan A's output. Its design may need
revision once M3 reveals how large the residue actually is — if the cascade resolves most of it,
Plan B's centre of gravity shifts from "catalogue the dark set" to "characterise the
`UNCHARACTERIZED_HOMOLOG` class", which is a different analysis.

Write Plan A now. Write Plan B after M3, with M3's numbers in hand.

---

## 13 · Risks and open decisions

| risk | mitigation |
|---|---|
| Pyrodigal ORF set differs from PlasAnn's, breaking comparability | Accepted and stated (§3). M1 quantifies the difference; the M3 head-to-head turns it into a result |
| T6 (NR) query set larger than expected → cost blowout | Measured at M2 on the pilot; if too large, T5/T6 order swaps or NR is restricted to `UNCHARACTERIZED_HOMOLOG` candidates |
| The dark residue shrinks to near-nothing after a full cascade | This would be a **finding**, not a failure — and it is the honest version of the current claim. §6.2's `UNCHARACTERIZED_HOMOLOG` class is the target set either way |
| Container availability (antiSMASH, PHROGs, DefenseFinder, ESMFold, Foldseek) | Verify all pulls at M0 before depending on them |
| Rewrite scope creep into a paper | Writing is explicitly out of scope until M5 |

**Open decisions, deferred deliberately:**

1. Lab assay format — decided after M6 shows how many high-confidence targets exist.
2. Whether to extend to the 64,658 excluded Simulated-artifact plasmids as a flagged class.
3. Whether the annotated-plasmid resource (M3) is published separately from the dark-ORF work.

---

## 14 · Available cluster resources (verified 2026-09-02)

**Modules:** Nextflow 25.10.2 · snakemake 8.27.0 · MMseqs2 18-8cc5c · DIAMOND 2.1.24 · HMMER 3.4 ·
PyHMMER 0.12.0 · HH-suite 3.3.0 · BLAST+ 2.17.0 · prodigal 2.6.3 · eggnog-mapper 2.1.13 ·
SignalP 6.0h · Infernal 1.1.5 · SeqKit 2.10.1 · AlphaFold 3.0.1 / 2.3.2-CUDA · apptainer 1.5.3

**Prebuilt MMseqs2 databases** (`/sw/data/MMseqs2_data/latest/`): UniRef50 · UniRef90 · UniRef100 ·
NR · UniProtKB (+Swiss-Prot, TrEMBL) · eggNOG · CDD · Pfam-A (full/seed) · Pfam-B · **PDB70** ·
**VOGDB** · GTDB · SILVA · Resfinder · dbCAN2 — all with taxonomy mapping.

**Other data:** `/sw/data/blast_databases/` (nr 248 volumes, swissprot, env_nr) ·
`/sw/data/Pfam/35.0` (project's local 38.2 is newer — keep it) ·
`/sw/data/alphafold_dataset` (9.5 TB, incl. 3.0.1) · `/sw/data/InterProScan_data/5.78-109.0`

**Needed as containers:** Pyrodigal · PHROGs · antiSMASH · BAGEL · DefenseFinder/PADLOC · ESMFold ·
Foldseek · (CheckV, optional)

**Storage:** `/gorilla/proj` 4.8 PB free.

**Note on FlashFold:** `data/flashfold_run/install.log` shows the failure is a hardcoded path into a
staged conda tree (`/sw/apps/conda/latest/rackham_stage/envs/flashfold/bin/python3`), not a quota
problem, and the AlphaFold2 weights downloaded successfully. It is fixable, but ESMFold via
container is the simpler route for T8 and is what this design assumes.
