# Scientific and feasibility review: `plans/2026-09-10-pipeline-v2-design.md`

Reviewed 2026-09-16 against `docs/annotation_statistics.md`, `docs/PARAMETER_PROVENANCE.md`,
`docs/PIPELINE_CODE.md` and `config/cascade.yaml`.

Scale used throughout, not re-derived: 143,503 plasmids in the analysis set, 9,317,050 ORFs,
3,497,616 unique protein sequences. Hardware: one node, 96 cores, 600 GB RAM, 7-day wall limit
= **16,128 core-hours per job**. Present on disk: Pfam-A 38.2, AntiFam, PHROGs, `nr.dmnd`
(350 GB), `swissprot.dmnd`, foldseek pdb + prostt5. Absent: InterProScan, eggNOG-mapper
databases, HH-suite.

---

## 1. Internal contradictions

### 1.1 The document specifies a hub; the measured design is a cascade
§10 runs seven sources over the full protein set and §11.2 says "allow adjudication after all
evidence is collected". The implemented and *measured* design is a narrowing cascade:
`narrow_at: 0.9`, `min_explained: 0.5`, a 2% sweep cohort, and tier order as authority order
(`config/cascade.yaml`). The spec never mentions narrowing, `narrow_at`, `min_explained`, the
sweep cohort or `-Z` — the four decisions `docs/annotation_statistics.md` exists to defend. This
is not a detail: under §10 every tier including nr sees all 3,497,616 proteins instead of the
residue, which multiplies the dominant cost by 1/(residual fraction). Either §10 supersedes the
cascade and the compute estimates in §3 below apply, or the spec is stale. It does not say which.

### 1.2 Tier identifiers collide with `config/cascade.yaml`
§10 numbers T1–T7 with T3 = InterPro, T5 = DIAMOND Swiss-Prot, T7 = DIAMOND nr. `cascade.yaml`
numbers T1–T4 with T3 = Swiss-Prot and T4 = nr, and reserves a T4 for PHROGs. `homology_depth`
is defined as the tier's position in that list, and tier ids are stamped into output rows. The
same string `T3` means InterPro in the design document and Swiss-Prot in the artefact it governs.
The spec also has no PHROGs tier at all, though the phage-plasmid rationale for one is written
into both `cascade.yaml` and `docs/PIPELINE_CODE.md` §15.

### 1.3 §10.1 specifies `hmmscan`, which breaks the `-Z` convention
`hmmscan` and `hmmsearch` interpret `-Z` differently: in `hmmsearch` it is the number of
*sequences* searched; in `hmmscan` it is the number of *profiles*. `cascade.yaml` pins
`hmmer_z: 3497616` precisely because it is the analysis-set sequence count
(`annotation_statistics.md` §3). Passed to `hmmscan` that number is a claim about a 3.5-million-
profile database against a real Pfam-A of ~21k families — every E-value inflated by two orders of
magnitude. The single measured throughput figure in the repo (65.8 s fixed + 0.0365 s/protein at
4 threads) is an `hmmsearch` slope and does not apply to `hmmscan` either.

### 1.4 20 aa versus 30 aa versus `min_orf_aa: 20`
§7.3 and §7.4 exclude proteins `< 20` aa from the discovery universe; §31's config block says
`min_dark_length_aa: 30`; `config/config.yaml` sets `min_orf_aa: 20` as a *gene-calling* floor.
Three names, two values, two different semantics. Worse, §7.4 mandates that sub-threshold
proteins be "labeled ... rather than deleted", but `min_orf_aa` is a Pyrodigal floor that deletes
them before dereplication, so under the current implementation `discovery_excluded_short` can
never be TRUE. At 30 aa the spec would also discard the length class that
`PARAMETER_PROVENANCE.md` cites FESNov's validated 36-residue peptide to protect.

### 1.5 Missing sections §2.6 and §13.4
Headings run 2.5 → 2.7 and 13.3 → 13.5. §13.5 is titled "Recommended *additional* field", so the
deleted §13.4 held the primary dark fields. The consequence is visible: §2.5 introduces the dark
sublabels `DARK_HYPOTHETICAL`, `DARK_DARK_DUF` (a doubled prefix, and a typo), and
`DARK_NO_SEQUENCE_FUNCTION`, and no column anywhere in §25 or §26 holds them — §25.4 carries only
`dark_status` and `dark_evidence_level`. The sublabels are defined and orphaned.

### 1.6 §19.3 collapses two orthogonal axes into one column
The four structural states are `NO_STRUCTURAL_MODEL` / `STRUCTURAL_MODEL_AVAILABLE` (was a model
produced) and `DARK_FOLD_KNOWN` / `DARK_NO_STRUCTURE` (did it match). §25.13 gives one column,
`structural_status`, to hold all four. A protein with a model *and* a hit is simultaneously
`STRUCTURAL_MODEL_AVAILABLE` and `DARK_FOLD_KNOWN`; a protein with no model is neither
`DARK_FOLD_KNOWN` nor `DARK_NO_STRUCTURE`, yet §2.7 says only those two states are used.

### 1.7 "No global score" versus `dark_evidence_level` and §35
§1.2 forbids "a final biological ranking score", §23 opens "must not collapse all evidence into
one numerical score", and Rule 8 (§30) repeats it. §13.5 then defines a 0–5 ordinal that is
exactly a collapse of recurrence, family conservation and genomic context into one number, and
§25.4 puts it in the primary deliverable. §35 completes the contradiction by listing "scores"
among what the annotation pipeline produces, in the section whose purpose is to draw that boundary.

### 1.8 §5.3's single `NA` destroys the distinction §5.4 requires
§5.3 mandates one missing-value token and forbids any alternative. `docs/PIPELINE_CODE.md` §13
records the opposite decision, taken for cause: `TOO_FEW_MEMBERS`, `NO_DIVERGENCE`, `SATURATED`,
`NO_OUTPUT` are distinct statuses "because a blank reads as a failed test rather than a test that
never ran". Under §5.3 a dN/dS that saturated, a family too small to test, and a crashed job are
the same cell. §32.3's "check missingness" then cannot distinguish them, and §5.4's own principle
— a label must coexist with the measurement behind it — is violated by the rule three paragraphs
above it.

### 1.9 §7.1 forbids the artefact detectors the implementation already runs
"Version 1 uses only Prodigal/Pyrodigal output properties and AntiFam. No additional artifact
detector should be introduced without a deliberate pipeline revision." The implementation runs a
tantan low-complexity screen (`max_low_complexity_fraction: 0.5`, `cascade.yaml`) and an
antisense-strand RNAcode check that `PIPELINE_CODE.md` §13 describes as "the artefact class
nothing else here catches". §20 and §20.1 then require `low_complexity_fraction` and a
`LOW_COMPLEXITY` liability label, i.e. the same screen under another name.

### 1.10 §7.2 specifies AntiFam with no threshold
§7.2 says only "Run AntiFam" and records `antifam_score` / `antifam_evalue`, and §31 says
`antifam: enabled: true`. The project has *measured* that a blanket `-E 1e-5` overrides the
curator on 274 of 278 profiles (98.6%), median curated cut E = 7.3e-4, and settled on `--cut_ga`
with a regression test (`annotation_statistics.md` §7; `PARAMETER_PROVENANCE.md`). The design
document that is supposed to be upstream of that decision does not record it, so an implementer
following the spec would reintroduce the exact defect the test exists to catch.

### 1.11 Deliverable format and row-unit disagree with the implementation
§25 and §38 name `.tsv`; `PIPELINE_CODE.md` §14 writes `.csv` for a stated reason (DIAMOND
`stitle` routinely contains commas). §38 lists both `protein_annotations_complete.tsv` and
`annotation_complete.tsv` without saying which is authoritative. §25 fixes one row per ORF
occurrence — 9.3M rows — while §25.10–§25.13 are per-family and per-unique-protein quantities;
replicating a family's dN/dS and Foldseek hit across every occurrence invites occurrence-level
counting of family-level evidence in any downstream `group by`.

### 1.12 Rule 1 (§30) is unsatisfiable given the rest of the document
"Do not invent biological thresholds that are not defined in this document or configuration."
Section 4 below lists ~25 places where the document defines no threshold, and §31's example
config writes the literal string `configurable` as the *value* of `relaxed_evalue` and of all
three MMseqs resolutions. An implementer obeying Rule 1 cannot execute §10.2, §15.5, §18.3,
§18.4, §19.4 or §13.

### 1.13 Plasmid count
§3.1 and §37 say 143,504. `config/config.yaml` and `annotation_statistics.md` §3 say 143,503
after `exclude_hab_top: [Simulated-artifact, Lab-artifact]`. One record's worth, but it is the
denominator of every prevalence in §21 and it is stamped into rows.

### 1.14 §32.6's inequalities are asserted and then withdrawn
"Verify: `host_count <= plasmid_count`, `genus_count <= host_count` ... Exact inequalities depend
on the unit definitions and must be encoded explicitly." A validation check that the same
paragraph declares undefined is not a check. `genus_count <= host_count` is true only if `host`
means strain; if `host` means species it is false whenever a genus is represented by one species.

---

## 2. Scientific validity

### 2.1 §13.2/§13.3 — the study's primary outcome variable is undefined and circular
`ANNOTATED` = "supports a sufficiently specific assignment according to the configured
adjudication rules"; `DARK` = "no sufficiently supported assignment ... according to the
configured adjudication rules". The rules are never given. The size, composition and every
downstream property of the dark set — the entire result — is set by a parameter the design does
not contain. The implementation does define it (`min_explained: 0.5` on informative coverage),
and `PARAMETER_PROVENANCE.md` records that value as **UNREFERENCED**. The spec therefore removes
the only concrete definition in the project and replaces it with a placeholder. No claim of the
form "X% of the plasmidome is dark" is supportable until this is a number with a defence.

### 2.2 §23 — `independent_evidence_categories` overstates independence
§23.2 correctly identifies Pfam / InterPro / HHsearch / DIAMOND / MMseqs as one homology signal.
§23.1 then lists `sequence_homology`, `domain_homology`, `orthology` and `structural_evidence` as
four separate categories — all four are homology detection at different sensitivities. The
structural category is the worst case: ProstT5 is a *sequence* language model that emits 3Di
directly, with no structure prediction step (`PARAMETER_PROVENANCE.md` records `min_plddt` being
removed for exactly this reason). A ProstT5+Foldseek hit is a learned sequence→sequence mapping,
not evidence independent of sequence. `evolutionary_conservation` and `distribution` are also not
independent: both are computed from the same MMseqs family and the same occurrence table. At most
three separable sources exist — homology (all tiers + foldseek + mmseqs), genomic context, and
coding-potential/selection — and the third is derived from the first's family assignment.

### 2.3 §14 — deposit counts measure submission sociology, not biology
`hypothetical_deposit_count`, `independent_database_count` and `independent_plasmid_count` are
proposed as evidence that a protein is real. All three are confounded, and in the same direction:
(i) RefSeq `WP_` records are propagated by PGAP annotation transfer and collapsed under
MULTISPECIES, so a "deposit" is a genome submission, not an independent gene call; (ii)
PlasmidScope, PLSDB, IMG/PR and NCBI share upstream records, so `independent_database_count`
counts aggregators of the same accessions; (iii) clonal over-sampling dominates — thousands of
near-identical ST131 / ST258 plasmids enter from single outbreak studies. The spec specifies no
dereplication by BioProject, submitter or ANI, though `workflow/scripts/clonal_registry.py`
already exists for it. §34.3 then elevates this to "high-confidence evidence of protein
existence", which is the strongest claim in the document resting on the weakest measurement.

### 2.4 §14 is maximally confounded with the artefact class §7.1 refuses to screen for
A shadow ORF on the reverse complement of a conserved gene is recurrent for precisely the reason
the real gene is recurrent: it is deposited in every genome carrying that gene, in every database,
across the whole taxonomic range. It therefore scores maximally on all five of §14's recurrence
fields. `cascade.yaml` states this explicitly — artefacts "pass the entire cascade cleanly"
because "our selection criterion is nothing named it", and they "survive multi-lineage and
selection tests". §14's recurrence evidence is not merely unable to discriminate real proteins
from shadow ORFs; it is *positively* enriched for them, while §7.1 forbids adding any detector
beyond AntiFam.

### 2.5 §13.5 — the 0–5 scale imposes an order the data do not have
Levels 2 (sequence recurrence), 4 (family conservation) and 5 (family conservation *plus*
genomic context) are three different axes presented as one monotone quantity. A singleton with a
conserved defence-adjacent context and a 10,000-member family with no context cannot be ordered,
yet the scale forces it. Level 3 ("independently recurrent hypothetical") inherits every
confound in §2.3 above. The scale is also a lossy recode of columns already in §25 — it adds no
information and one arbitrary ordering — and §13.5's disclaimer ("do not use as the final
experimental score") is unenforceable once §25.4 ships it in the primary table and §36's
candidate classes are the same axes.

### 2.6 §18.4 — dN/dS is not applicable to most plasmid gene families
Three independent failures. (i) *Saturation*: families are built at 30% amino-acid identity
(FESNov, `PARAMETER_PROVENANCE.md`), where synonymous sites are saturated and dS is not
estimable; the implementation's `SATURATED` status exists because this is the common case.
(ii) *Within-population comparison*: most family members are near-identical sequences from the
same species, and dN/dS computed on within-population polymorphism is biased toward 1 and is not
a test of selection (Kryazhimskiy & Plotkin 2008) — so both tails of the divergence distribution
are invalid and only a narrow middle band is interpretable. (iii) *Recombination*: plasmid genes
recombine, and pairwise Nei–Gojobori with no tree conflates recombination breakpoints with
selection. §18.4's four gates ("enough", "adequate", "defensible", "appropriate") are the places
these would have to be excluded, and all four are undefined. `PIPELINE_CODE.md` §15 adds that
Nei–Gojobori is now the only estimator, with nothing able to contradict it.

### 2.7 §18.3 — RNAcode, and an uncorrected p-value across ~10^5–10^6 families
RNAcode tests for a coding signature in an alignment of orthologous nucleotide sequence and was
calibrated on moderately divergent bacterial alignments. At 30% amino-acid identity the
underlying codon alignment is near saturation and RNAcode loses power; at the clonal end there is
no divergence and it has no signal. Those two regimes bracket most of the dark set. Separately,
`rnacode_max_p: 0.05` is cited to FESNov as a per-alignment threshold, and the spec specifies no
multiple-testing correction — applied across even 300,000 families that is ~15,000 families firing
a coding claim by chance. §18.3's gate, "where sufficient homologous nucleotide sequence is
available", never says how much.

### 2.8 §2.7 / §19 — `DARK_NO_STRUCTURE` is correlated with darkness by construction
ProstT5's 3Di prediction is learned from the same evolutionary signal that the sequence searches
use. For a protein with no detectable homologues — which is the definition of the dark set — the
predicted 3Di string is least reliable, so a Foldseek miss is expected whether or not the fold is
known. `DARK_NO_STRUCTURE` therefore does not partition the dark set on structure; it partly
re-measures the sequence darkness that put the protein there. §34.5 and §1.2 correctly forbid
reading it as "novel fold", but §36 then lists "dark proteins with no known structural match" as
a candidate class, which is the same claim under a different name.

### 2.9 §19.1 — `structure_confidence` cannot be produced by the specified workflow
§19.1 names ProstT5 and requires `structure_confidence` and `model_length`. ProstT5 does not
predict a structure and emits no pLDDT; `PARAMETER_PROVENANCE.md` records `structure.min_plddt`
being *removed* for exactly this reason, after it was found declared, schema-required and read by
nothing. Reintroducing a confidence field without naming its estimator reinstates the defect: the
`novel_fold` stratum becomes silently ungated again.

### 2.10 §12.1 — the evidence hierarchy inverts the project's own authority order
Level 3 ("high-confidence orthology with informative function", i.e. eggNOG) outranks level 4
("strong curated sequence homology", i.e. Swiss-Prot). `cascade.yaml` places Swiss-Prot at T3 and
describes it as "highest label quality available"; eggNOG OG assignment is an automated ortholog
transfer. No justification is given for the inversion, and §12.2's `primary_annotation_confidence`
has no defined value set, so the hierarchy cannot be applied as written.

### 2.11 §9.4 — the required fields cannot support the required decision
"An annotation is not considered functionally characterized merely because it appears in
Swiss-Prot ... Evidence provenance must be inspected." The provenance that matters is the
*target's*: UniProt ECO evidence codes, protein-existence (PE) level, annotation score. §10.5
requires none of these — only target id, description, bitscore, E-value, identity and coverages.
The `FUNCTIONALLY_CHARACTERIZED` versus `CURATED_UNKNOWN_FUNCTION` distinction in §9.2 is
therefore unassignable from the data the spec collects.

### 2.12 §21 — normalization is mandated but no model is specified, and the dominant confound is absent
§21 correctly names length, gene count, taxonomy, database and sampling bias, then defers "exact
statistical tests" to a later module — while §21.3 and §25 require a populated
`normalized_prevalence` column. You cannot emit a normalized prevalence without choosing the
model. The confound that dominates this particular dataset is not listed: clonal redundancy, which
makes `plasmid_count` not a count of independent observations. `habitat` is used as a stratifier
in §21.2 and §21.4 with no definition, against this project's standing rule that comparisons are
built from `hab_sub` + `is_clinical` and never from `hab_top`.

### 2.13 §16.1/§13 — darkness is per sequence, context is per occurrence, and the spec never joins them
Dark status is adjudicated on a unique protein (§13 operates on the annotation evidence set);
genomic context, synteny and every §17 conservation measure are properties of occurrences. §25
puts `dark_status` on an occurrence row. A protein present on 4,000 near-identical plasmids gets
4,000 "independent" neighbourhood observations that are one observation. §17.3 warns against
over-reading context recurrence but no field records how many *independent* backbones the
recurrence spans.

---

## 3. Feasibility at scale

Budget: **16,128 core-hours per 7-day job** (96 cores). All figures are for 3,497,616 unique
proteins. Where the repo has a measurement it is used; where it does not, the arithmetic is shown
and bracketed, because neither the spec nor `docs/` benchmarks any stage except T1
(`PIPELINE_CODE.md` §15 lists the nr tier as unbenchmarked, and structure not at all).

### 3.1 §10 annotation hub

| §10 source | core-hours | wall on 96 cores | verdict |
|---|---:|---:|---|
| T1 Pfam GA, `hmmsearch` (measured 0.146 core-s/protein) | ~142 | 1.5 h | fine |
| T1 Pfam GA **as §10.1 specifies, `hmmscan`** (10–20x slower) | 1,400–2,800 | 15–30 h | wasteful; breaks `-Z` (§1.3) |
| T2 Pfam relaxed, full set rather than residue | 140–2,800 | 1.5–30 h | fine |
| T3 InterProScan — **not installed** | 3,000–20,000 | 1.3–8.7 d | must be subset |
| T4 eggNOG-mapper — **databases not installed** | 1,500–4,000 | 16–42 h | feasible if installed |
| T5 DIAMOND Swiss-Prot | 20–60 | <1 h | fine |
| T6 HH-suite — **not installed** | 20,000–100,000 | 9–43 d | **impossible** |
| T7 DIAMOND nr `--very-sensitive`, all 3.5M | 10,000–50,000 | 4–22 d | **at/over the limit** |

**T6 HMM-HMM is the clearest impossibility.** `hhblits` needs 2–4 iterations against a
UniClust/BFD-scale database (100–270 GB, not on disk) to build each query MSA, at ~10–60 core-s
per query; 3.5M queries is 10^4–10^5 core-hours, i.e. 1.2–6 times the *entire* 7-day allocation
even before `hhsearch`. It can only run on family representatives: at ~600k dark representatives
it is still 3,300–10,000 core-hours, and at a shortlist of 10k it is ~60–170 core-hours. The spec
already concedes the tool is optional ("if configured", §10.6), so the honest fix is to strike it
from §10 and define it as a shortlist-only confirmatory stage.

**T7 nr is the second, and the spec makes it worse.** §10's hub removes the narrowing that made
nr affordable: under `cascade.yaml` T4 sees only the unresolved residue, under §10 it sees all
3,497,616. Add the fixed cost the repo flags as unmeasured: DIAMOND streams the whole 350 GB
database per invocation, so 64 shards pay it 64 times — 6.2 h of pure I/O at 1 GB/s, 31 h at a
realistic 200 MB/s shared filesystem, before a single alignment. `workflow/bench_nr.sbatch` exists
to measure this and has not been run; no schedule for §10 is credible until it has.

**T3 InterProScan has a second problem beyond not being installed.** Its cheap path is the
precalculated match lookup against UniParc, which novel plasmid ORFs miss by construction, so
essentially all 3.5M proteins take the full local calculation across Gene3D, SUPERFAMILY, PANTHER
and CDD. Those are also precisely the members that §10.3 itself declares non-independent of Pfam,
which is already T1/T2. The most expensive tool in the hub is being bought for the evidence the
same section says not to count twice.

### 3.2 Downstream stages

| stage | core-hours | wall on 96 cores | verdict |
|---|---:|---:|---|
| §15 MMseqs2, 3 resolutions (§15.5) | 150–600 | 2–6 h | fine; RAM is the constraint, needs `--split-memory-limit` |
| §16.5 DefenseFinder over 143,503 replicons | 400–1,500 | 4–16 h | fine **if batched**; 143,503 MacSyFinder process launches is not |
| §16.6 IntegronFinder `--local-max` over 143,503 | 1,200–4,000 | 13–42 h | fine |
| §18.3 RNAcode + codon alignments, 0.3–0.8M families | 300–1,500 | 3–16 h | fine (alignment dominates, not RNAcode) |
| §18.4 dN/dS, Nei–Gojobori, <=50 members (1.95 s worst family) | 100–400 | 1–4 h | fine |
| §19 ProstT5, **dark family representatives** (~0.6M) | ~12,000 | ~5 d | at the limit, CPU-only |
| §19 ProstT5, **all dark proteins** (~2.6M) as §19 demands | ~52,000 | ~22 d | **impossible** |
| §19.2 Foldseek vs PDB, 2.6M 3Di queries | 50–300 | <3 h | fine |
| §19.2 Foldseek vs AFDB, if that is the "appropriate" database | 5,000–20,000 | 2–9 d | over |
| §20 SignalP6 / DeepTMHMM / disorder over 3.5M | 300–1,200 | 3–13 h | fine; licence, not compute, is the blocker |

### 3.3 §19's "all dark proteins, not family representatives" — the cost quantified
ProstT5 is a 3-billion-parameter T5. A 300-residue protein costs roughly 2 x 3e9 FLOP/token over
~600 encoder+decoder tokens ~= 3.6 TFLOP. A core sustains ~50 GFLOPS, so ~70 core-s per protein;
**no GPU appears in the stated hardware**. Taking the dark set at ~2.6M unique proteins (60–90% of
3.5M) and dark family representatives at ~600k after 30%/50% clustering:

- all dark proteins: ~52,000 core-hours = **22 days wall on the full node, 3.2x the 7-day budget**
- representatives only: ~12,000 core-hours = ~5 days, inside one job
- **the §19 choice costs ~4.3x, and is the difference between a stage that runs and one that cannot**

§19 gives no reason for it. The reason it would need is that family members' 3Di strings differ
usefully — but members within a 30%-identity cluster are exactly the set ProstT5 will map to
near-identical 3Di, so the extra 40,000 core-hours buys a prediction that is redundant by
construction. On a GPU (~0.04 s/protein) the whole dark set is ~26 h and the question disappears;
acquiring one is the cheaper fix than restricting the stage.

### 3.4 Totals
§10 as written, at full scale on all 3.5M proteins: **~36,000–180,000 core-hours = 15–78 days of
the entire node**, i.e. 2–11 times a single 7-day allocation. Adding §19's all-dark structural
stage puts the floor near 90,000 core-hours. Three stages must be restricted to a subset or
struck — T6 HH-suite (strike or shortlist), T3 InterProScan (representatives only), T7 nr (restore
the cascade residue) — and §19 must run on representatives unless a GPU is provisioned. The rest
of the pipeline, roughly 2,500–9,000 core-hours, fits comfortably.

---

## 4. Missing definitions that block implementation

Each line is a biological threshold an implementer must invent, in violation of §30 Rule 1.

| § | phrase | what must be decided |
|---|---|---|
| 13.2, 13.3 | "sufficiently specific / supported ... configured adjudication rules" | **The adjudication rule itself** — the definition of the study's outcome variable. Coverage fraction, per-tier E-values, whether a DUF counts as annotated |
| 13.5 | "numerical boundaries should be treated as configuration" | five boundaries for levels 1–5, and which axis decides each |
| 12.4 | "configured annotation-transfer criteria" | identity + coverage + E-value permitting transfer (no twilight-zone cut is named) |
| 19.4, 2.7 | "convincing structural match" | Foldseek E-value, TM-score, and alignment-coverage cut |
| 19.2 | "an appropriate structural database" | PDB vs AFDB vs both — 100x cost difference, and it changes what `DARK_NO_STRUCTURE` means |
| 19.1 | `structure_confidence` | the estimator; ProstT5 supplies none (§2.9) |
| 18.3 | "where sufficient homologous nucleotide sequence is available" | minimum members, minimum aligned length, divergence window, multiple-testing correction |
| 18.4 | "enough", "adequate", "defensible", "appropriate" | four separate gates, plus the test behind `selection_significance` and its threshold |
| 18.2 | "measure where appropriate"; `conserved_motifs` | which motif finder, what significance, what counts as conserved |
| 16.7 | seven context labels | the adjacency rule and distance cut that makes each TRUE; §24's example implies "median 2 genes, 17 plasmids" and neither is stated as a rule |
| 17.2 | six conservation measures | the formula for each: over occurrences or family members, and neighbour identity by family / Pfam / function class |
| 16.5, 16.6 | "where configured" | whether DefenseFinder and IntegronFinder run at all |
| 22 | six rarity/conservation labels | the count or breadth cut for `rare`, `lineage_specific`, `widely_conserved`, `cross_MOB`, `cross_host`, `cross_taxon` |
| 15.5, 31 | `close` / `intermediate` / `broad` = `configurable` | three identity + coverage + coverage-mode triples |
| 10.2, 31 | `relaxed_evalue: configurable` | the T2 E-value — and `annotation_statistics.md` §2 shows any value above ~1e-3 is needed to open a sub-GA band at all |
| 21.3 | `normalized_prevalence` | the normalization model; §21.5 lists four candidates and picks none |
| 23.3 | `independent_evidence_categories` | the machine-checkable independence rule; §23.2 is prose |
| 24 | `functional_hypothesis` | the firing condition for each of the ten hypotheses |
| 20, 20.1 | `disorder_fraction`, `coiled_coil`, `secondary_structure_summary`, six liability labels | tool for each, and the cut for each label |
| 9.2, 9.3 | ten normalized classes, `normalization_rule` | the string→class mapping, and how `FUNCTIONALLY_CHARACTERIZED` is decided (§2.11) |
| 11 | search / evidence / reporting thresholds | three named threshold types; the document sets no value for any of them anywhere |
| 7.2 | AntiFam | no threshold given (§1.10) |
| 32.6 | "exact inequalities depend on the unit definitions" | the unit definitions, and therefore the checks |
| 14 | `hypothetical_recurrence_level` | its levels, and what makes a source "independent" |
| 5.2 | "input dataset hash when appropriate" | when |

§39's success criteria inherit the problem: "evolutionary measurements are available **where
valid**" and "structural status is available **where attempted**" are satisfied by a run that
produced zero dN/dS values and attempted no structures.

---

## 5. Thresholds without a source

Standing rule (`docs/PARAMETER_PROVENANCE.md`): every parameter carries a citation to a paper
that used the same value, or an explicit measured justification for deviating. New or implied by
this spec, with neither:

1. **§31 `min_dark_length_aa: 30`** — no source, and it contradicts §7.3's 20 and
   `config.yaml`'s `min_orf_aa: 20`. The 20 has a measured defence (FESNov's validated 36-residue
   peptide; the floor costs ~0.14% of ORFs); the 30 has nothing and would begin cutting into that
   length class.
2. **§13.5 levels 0–5** — five boundaries on the field that §25.4 ships in the primary table.
3. **§16.2 / §31 `neighborhood_genes: 10`** — the +/-3 summary window is FESNov-cited
   (`context.neighbourhood_window: 3`); the +/-10 extraction window is new and uncited.
4. **§24 "median distance = 2 genes"** and **"observed in 17 independent plasmids"** — written as
   an illustration, but they are the de facto firing conditions for `defence_associated` and no
   other numbers are offered.
5. **§10.2 relaxed Pfam E-value** — undeclared. Any value is a deviation from the *measured*
   `DECIDED` entry in `cascade.yaml` T2, which records that 1e-5 stays and explains what loosening
   costs at 3.5M queries.
6. **§15.5 the second and third clustering resolutions** — FESNov supplies one (30% id / 50% cov);
   `intermediate` and `broad` have no source and no measurement.
7. **§11.1** forbids "a universal 90% explained threshold" while the implementation runs
   `narrow_at: 0.9`, already listed UNREFERENCED. The spec neither justifies nor replaces it, so
   the number is simultaneously in use and prohibited.
8. **§19.4 the Foldseek significance cut** — `structure.max_evalue: 1e-3` is already flagged
   UNREFERENCED; the spec adds `foldseek_tm_score_if_available` and `foldseek_alignment_coverage`
   as further ungated axes.
9. **§14 recurrence-level boundaries** — every cut on `hypothetical_deposit_count`,
   `independent_source_count`, `independent_database_count`, `independent_plasmid_count`.
10. **§22 the six rarity/conservation cuts.**
11. **§18.3 multiple-testing correction** — `rnacode_max_p: 0.05` is FESNov-cited per alignment;
    applying it across ~10^5–10^6 families without correction is a deviation from the cited use
    and needs its own justification (§2.7).
12. **§12.1's eight-level ordering** — a claim, not a measurement, and it inverts `cascade.yaml`'s
    authority order at levels 3/4 (§2.10). It needs the same explicit methods defence that
    `targets.RANK_PRIORITY` is given in `PARAMETER_PROVENANCE.md`.
13. **§21.5 the normalization method** — "per-gene", "per-plasmid", "size-matched", "stratified"
    are offered as alternatives with no basis for choosing, and the choice determines every
    enrichment claim in §21.4.

Not resolved by this spec and still outstanding from `PARAMETER_PROVENANCE.md`: the five
`peptide.py` constants outside config entirely, which gate 325 of the 1,000 constructs.
