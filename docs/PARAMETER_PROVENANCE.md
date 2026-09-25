# Parameter provenance

Standing rule: **every parameter carries a citation to a paper that used the same value, or
an explicit measured justification for deviating.** A value with neither is a known weakness,
recorded here as such rather than quietly kept.

An unreferenced parameter is a researcher degree of freedom. There are ~25 here; even at
three defensible values each that is 10^12 possible pipelines, and a result can drift toward
the desired answer with nobody intending it. This project has already been damaged twice this
way — the unswept cascade thresholds in v1, and the invented scoring coefficients in v2.

**A citation is not automatic justification.** Pavlopoulos' >=100-member cluster minimum is
wrong for this dataset (they had 1.17 billion sequences; we have 3.5 million), so deviating
is correct — but the deviation is recorded with the measurement that motivated it.

## Sources

| key | paper |
|---|---|
| **FESNov** | Rodríguez del Río Á. *et al.* Functional and evolutionary significance of unknown genes from uncultivated taxa. *Nature* **626**, 377–384 (2024) |
| **Pavlopoulos** | Pavlopoulos G.A. *et al.* Unraveling the functional dark matter through global metagenomics. *Nature* **622**, 594–602 (2023) |
| **Durairaj** | Durairaj J. *et al.* Uncovering new families and folds in the natural protein universe. *Nature* **622**, 646–653 (2023) |
| **ECLIPSE** | Lata S. & Heinz D.W. ECLIPSE: exploring the dark proteome of ESKAPE pathogens. *Bioinformatics* **42**(8) (2026) |
| **Pfam** | Mistry J. *et al.* Pfam: The protein families database in 2021. *Nucleic Acids Res.* **49**, D412–D419 |
| **AntiFam** | Eberhardt R.Y. *et al.* AntiFam: a tool to help identify spurious ORFs in protein annotation. *Database* **2012**, bas003 |
| **DefenseFinder** | Tesson F. *et al.* Systematic and quantitative view of the antiviral arsenal of prokaryotes. *Nat. Commun.* **13**, 2561 (2022) |
| **HMMER** | Eddy S.R. Accelerated profile HMM searches. *PLoS Comput. Biol.* **7**, e1002195 (2011) |
| **PlasmidScope** | Li Y. *et al.* PlasmidScope: a comprehensive plasmid database with rich annotations and online analytical tools. *Nucleic Acids Res.* **53**, D179–D188 (2025), doi:10.1093/nar/gkae930 |
| **ISEScan** | Xie Z. & Tang H. ISEScan: automated identification of insertion sequence elements in prokaryotic genomes. *Bioinformatics* **33**, 3340–3347 (2017) |
| **eggNOG-mapper** | Cantalapiedra C.P. *et al.* eggNOG-mapper v2. *Mol. Biol. Evol.* **38**, 5825–5829 (2021) |
| **UniRef** | Suzek B.E. *et al.* UniRef clusters: a comprehensive and scalable alternative for improving sequence similarity searches. *Bioinformatics* **31**, 926–932 (2015) |
| **Smillie** | Smillie C. *et al.* Mobility of plasmids. *Microbiol. Mol. Biol. Rev.* **74**, 434–452 (2010) |
| **CheckV** | Nayfach S. *et al.* CheckV assesses the quality and completeness of metagenome-assembled viral genomes. *Nat. Biotechnol.* **39**, 578–585 (2021) |
| **NCBI PGAP** | Li W. *et al.* RefSeq: expanding the Prokaryotic Genome Annotation Pipeline reach with protein family model curation. *Nucleic Acids Res.* **49**, D1020–D1028 (2021) |

---

## Referenced — value matches a published use

| parameter | value | source | what they used it for |
|---|---|---|---|
| `clustering.min_seq_id` | 0.30 | FESNov | MMseqs2 deep-homology clustering of 400M genes |
| `clustering.coverage` | 0.50 | FESNov | same |
| `clustering.cov_mode` / `cluster_mode` | 1 / 2 | FESNov | verbatim from its Methods: "--min-seq-id 0.3 -c 0.5 --cov-mode 1 --cluster-mode 2 -e 0.001"; -e 0.001 is the MMseqs2 default. Applied at all three resolutions (was 0 / 0, which required 50% coverage of both proteins and had no source) |
| cascade selection (S2s) | family holds an unexplained small-plasmid protein | UniRef50 (Suzek et al. 2015); Durairaj et al. 2023 | the selecting clusters ARE the families: `clustering.primary`, intermediate since 2026-09-25 (50% identity, 80% coverage of the member, cov-mode 1, cluster-mode 2; the UniRef50 analogue Durairaj et al. built their network on). Made before the cascade (`protein_clustering.py`) and applied by `plasmidann.selection.select`. No annotation is transferred at this level |
| `search_clustering.min_seq_id` | 0.90 | UniRef | UniRef90: members of a 90% cluster are represented by one sequence for similarity search. The representative's cascade result is copied to its members (`cascade_resolve.py`, `annot_source = representative`) |
| `evolution.min_members_for_dnds` | 3 | FESNov | ">=3 complete genes" required per family |
| `evolution.dnds_purifying_max` | 0.5 | FESNov | purifying-selection filter on novel families |
| `evolution.rnacode_max_p` | 0.05 | FESNov | RNAcode coding-potential test. Now actually applied, on both strands. |
| `context.max_operon_gap` | 100 nt | FESNov | "absence of intergenic regions >100 nucleotides" |
| `context.neighbourhood_window` | 3 | FESNov | neighbouring genes at +/-3 positions |
| `context.min_context_conservation` | 0.50 | FESNov | 52,793 families annotated at >=50% confidence |
| `context.high_confidence_conservation` | 0.90 | FESNov | 4,349 families at >=90% confidence |
| `quality_gate.min_control_recall` | 0.99 | ECLIPSE | 99.2–100% of 246 virulence / 42 AMR / 75 essential genes recovered as annotated |
| T2, T3, T4 `max_evalue` | 1e-5 | FESNov | Pfam, AntiFam and pVOG searches (pVOG with 50% coverage) |
| `orf.min_terminal_repeat_bp` | 20 | CheckV | direct-terminal-repeat criterion; S0 removes one copy of the repeat from a circular record (measured: 400 of 400 sampled 'direct terminal repeat' records carry one, 76% of a length not divisible by 3) |
| `DOMAIN_NAMED` ("X domain-containing protein", "X family protein") | DOMAIN_ONLY | NCBI PGAP | PGAP gives these names to domain-level or family-level HMM assignments, not full-length functions; such a name alone never makes a protein FUNCTIONAL (`cascade.classify`) |
| `network.min_cov` | 0.50 of either protein | Durairaj | edge rule of the UniRef50 sequence similarity network |
| `network.max_evalue` | 1e-4 | Durairaj | same |
| `network.max_out_edges` | 4 | Durairaj | "a maximum of four outbound edges were considered per node" |
| `network.dark_brightness` | 0.05 | Durairaj | dark = functional brightness not above 5% |
| T5 database | NCBI ClusteredNR 2026-08-30 | UniRef | nr clustered at 90% identity / 90% length, searched by representative - the practice UniRef90 established for similarity search. Replaces the 2025-03-03 nr: newer, and 0.64x the residues (`tools/download_clustered_nr.sh`) |
| T1 `--cut_ga` | — | Pfam | curated per-family gathering thresholds define family membership |
| DefenseFinder `inter_gene_max_space` etc. | per model | DefenseFinder | shipped in the 711 system definition XMLs; referenced by construction |

## Deviating from a paper, with a measurement

| parameter | our value | paper value | why we deviate |
|---|---|---|---|
| `artefact_screen.antifam_args` | `--cut_ga` | FESNov used E <= 1e-5 | AntiFam ships a curated GA line for **all 278 profiles**. Measured directly from `data/refs/antifam/AntiFam.hmm` by converting each profile's GA to an E-value with its own `STATS LOCAL FORWARD` parameters at Z=3,497,616: **274 of 278 (98.6%) curated thresholds are LOOSER than E=1e-5**, the median curated cut corresponds to E=7.3e-4 (73x looser) and the loosest, `Spurious_ORF_67`, to E=0.195. A blanket 1e-5 floor therefore overrides the curator across essentially the whole artefact database, in the one stage that stops shadow ORFs reaching the plate. This is the same argument that already sets `--cut_ga` on T1, applied where the consequence is worse: a missed artefact is not a missed annotation, it is a non-protein sent to the bench. Regression test: `tests/test_scripts_smoke.py::test_the_artefact_screen_uses_antifams_curated_thresholds`. |
| `search_clustering.coverage` / `cov_mode` | 0.8 of BOTH sequences (cov-mode 0) | UniRef: 80% overlap with the seed (the longest member) | Stricter on purpose. With coverage of the member only (cov-mode 1, as in the families) a 120-aa member can sit inside a 300-aa representative, and the representative's Pfam hit - copied to the member - may lie in the part the member lacks. Requiring 80% of both keeps lengths within ~20%. Cost, measured on PlasmidScope's proteins (2026-09-24): 216,546 representatives instead of 196,973 of 356,959 selected proteins, about 8 h more nr |
| T5 `max_evalue` | 1e-5, no coverage floor | FESNov: 1e-3 with >50% query coverage for eggNOG and RefSeq | Kept at the cascade's 1e-5, stricter than FESNov's E-value, because a false positive at the deepest tier removes a genuine dark protein for good; and no coverage floor, because coverage decides FUNCTIONAL versus DOMAIN_ONLY afterwards (`min_coverage`) rather than whether a hit exists. DIAMOND scales E-values by database size, and ClusteredNR is 0.64x the letters of nr, so the same 1e-5 is ~0.65 bits more permissive than on full nr |
| T5 `skip_if_named_by` | pfam, swissprot | none | Decided 2026-09-24: a protein Pfam or Swiss-Prot named keeps that curated name and is not searched against nr, so nr cannot replace a curated label with free text (spec section 20). Consequence stated: the unexplained part of a DOMAIN_ONLY protein Pfam named is not searched in nr |
| `input.max_plasmid_size_bp` | 20,000 | Smillie describe the bimodal plasmid size distribution; no paper fixes a cut-off | The antimode of this collection, measured: Gaussian KDE of log10(size_bp) over the 143,503 plasmids (SciPy `gaussian_kde`, Scott bandwidth) has modes at 4,831 and 96,828 bp and its antimode at 19,011 bp; 20 kb is that value rounded. The earlier small-cryptic analysis used 10 kb, which remains a subset |
| cluster size minimum | none (ORPHAN label) | Pavlopoulos required >=100 members | They clustered 1.17 billion sequences; we have 3.5 million, 300x smaller. A 100-member floor would fragment real plasmid families into singletons. FESNov's >=3 is the right order of magnitude, and even that is reported rather than gated. |
| `library.min_length` | none | ECLIPSE used >=300 aa | Their floor was calibrated to virulence-gene lengths in one pathogen. We are hunting small ORFs: FESNov's validated antimicrobial peptide was **36 residues**. |

## Removed, because the method cannot honour it

| parameter | what happened |
|---|---|
| `structure.min_plddt` (was 70, cited to FESNov) | Declared, schema-required, and **read by nothing** — so the `novel_fold` stratum was silently ungated and a family whose structure was too poor to trust looked identical to a confident novel fold. It cannot be honoured here in principle: pLDDT is emitted by a structure *predictor*, and Foldseek with ProstT5 does not predict a structure — it translates sequence straight into the 3Di alphabet, which is what makes it cheap enough to run over the whole dark set. Removed from config and schema rather than left as a citation for something that was not happening. The confidence axis is the Foldseek E-value; restoring a pLDDT gate needs ColabFold or ESMFold as a confirmatory stage on the shortlist. |

## Thresholds that live in Python, not config

P4 says every threshold is declared in config, schema-validated, and stamped into the row it
governs. These are the exceptions, found by review. The reality-test thresholds have since
been moved; the peptide constants have not.

| parameter | value | status |
|---|---|---|
| `targets.REALITY_TESTS` dN/dS cutoff | 0.5 | **Fixed.** Was a literal in a lambda that duplicated `evolution.dnds_purifying_max`. Now read from that config value, so the two cannot diverge. |
| `targets.REALITY_TESTS` member minimum | 3 | **Fixed.** Was a literal duplicating `evolution.min_members_for_dnds`. Now read from it, and `targets.check_reality_config` refuses a config where they differ — because the `purifying_selection ⇒ is_family` nesting the eligibility count relies on is true only while they are equal. |
| `evidence.reality_thresholds` lineage minimum (`min_lineages`) | 2 | Definitional: "more than one". Counted over Stage 6 Mash lineages since 2026-09-25; it was a MOB-suite cluster minimum, and MOB-suite assigns the nearest reference's cluster however distant. |
| `rarity.cross_min_hosts` / `rarity.cross_min_genera` | 2 / 2 | Definitional: "more than one" observed species / genus (were 5 and 3, no source). `SINGLE_MOB`/`CROSS_MOB` are fixed at 1 and >= 2 MOB-suite clusters; `cross_min_mob` was removed. |
| `peptide.AMP_MAX_LENGTH` | 100 | **Outside config. No source.** |
| `peptide.AMP_MIN_CHARGE` | 2.0 | **Outside config. No source.** |
| `peptide.AMP_MIN_HYDROPHOBIC_FRACTION` | 0.3 | **Outside config. No source.** |
| `peptide.TM_WINDOW` | 19 | **Outside config. No source.** Approximates a membrane-spanning helix length, which is conventional, but no cited calibration. |
| `peptide.TM_THRESHOLD` | 1.6 | **Outside config. No source.** |

The five `peptide.py` constants gate two strata — `small_cationic_peptide` (125) and
`membrane_or_secreted` (200) — so **325 of the 1,000 constructs, 32.5% of the deliverable,
are assigned by numbers with no citation and no schema.** None of the seven source papers
covers antimicrobial-peptide or transmembrane biophysical cutoffs, so a citation will have to
come from a different literature. Until then this is the largest single provenance gap in the
pipeline, and it is larger than the count below suggests, because these do not appear in it.

## UNREFERENCED — known weakness, must be resolved

These have no published source. Each is a researcher degree of freedom until it gets one.

| parameter | value | status |
|---|---|---|
| `narrow_at` | 0.9 | **No source.** Chosen to be permissive enough to preserve the counterfactual. The 2% sweep cohort measures its cost, which is a partial substitute for a citation but not a replacement. |
| `min_explained` | 0.5 | **No source.** Applied post hoc and sweepable, so its cost is measurable — but the value itself is unjustified. |
| `min_coverage` | 0.5 | Partial: 50% coverage appears in FESNov and ECLIPSE for *clustering*, not for FUNCTIONAL/DOMAIN_ONLY classification. Not the same use. |
| `full_at` / `partial_at` | 0.8 / 0.5 | **No source.** |
| `max_target_seqs` | 25 | DIAMOND's default. Raised from 5 (no source) after 44% of nr benchmark queries filled all five slots with uninformative hits; the limit also changes DIAMOND's search, not only its reporting (Shah et al., *Bioinformatics* **35**, 1613 (2019)). |
| `sweep_cohort_fraction` | 0.02 | **No source.** A cost/benefit choice, not a scientific one — but still undeclared in the literature. |
| `max_low_complexity_fraction` | 0.5 | **No source.** |
| `evolution.min_codons` | 20 | **No source.** Now actually applied: `dnds_detail(..., min_codons=)` is passed the configured value per pair by S7b, and the arithmetic floor of 3 remains underneath it. Until that wiring existed the declared 20 was inert and an 8-codon fragment could fire `purifying_selection`, the strongest of the four reality tests. |
| `evolution.max_members_aligned` | 50 | **No source.** A compute cap; measured worst-case family 1.95 s. |
| `structure.max_evalue` | 1e-3 | **No source.** |
| `prioritisation.min_reality_lines` | 2 | **No source.** Sensitivity across 1–4 is reported, which shows its cost but does not justify the choice. |
| `portfolio.strata` quotas | 200/175/175/175/125/75/75 | **No source.** A design judgement about experimental portfolio balance; arguably not the kind of number a paper can supply, but it must be defended explicitly in the methods. |
| `library.length_liability_above_aa` | 400 | **No source.** |
| `hmmer_z` | 3,498,616 | Not a threshold — it is every unique protein of the analysis set (3,497,616) plus the 1,000 controls, searched or not, and pinning it is what makes E-values comparable across tiers (HMMER). Proteins resolved at Tier 0 stay in it, so E-values do not depend on PlasmidScope's coverage. The *practice* is standard; the specific value is simply our data. |

## Declared procedures that are not parameters

Two things in S9 look like tuning and are not. They are recorded here so a reader does not
go hunting for a citation that cannot exist.

| item | what it is |
|---|---|
| S8e, ISEScan | Not a threshold of ours. ISEScan 1.7.3 is run with its published defaults (ISEScan), so partial elements are kept; `--removeShortIS` is not used, because a partial IS on a plasmid is still an IS-derived region. |
| Tier 0, PlasmidScope transfer | Not a threshold of ours. A protein identical in sequence to a PlasmidScope protein takes PlasmidScope's published eggNOG-mapper 2.1.12 result (PlasmidScope; eggNOG-mapper, default settings as published). It counts as annotated when that result has a KEGG KO, an EC number, or a Pfam family that passes `cascade.is_informative` (so not a DUF/UPF family) — the eggNOG fields as released, with no score cut-off applied by us. A COG/OG category letter alone does not count. |
| `targets.RANK_PRIORITY` | The lexicographic ranking order within a stratum, six keys deep. Each key is a **claim**, not a coefficient: "more independent evidence beats a better hypothesis", "a family on more independent plasmid backbones is better supported than one on fewer". A reader can agree or disagree with each in turn, which is exactly what a weighted blend prevents. The last key is the family id, present so that ties are reproducible and *stated* rather than decided by hash order. |
| `targets.IMPLIED_BY` | Which reality tests entail which others. `purifying_selection` requires a measured dN/dS, which requires >=3 members, which is `is_family` — so those two are one line of evidence, not two, and `min_reality_lines` counts only the independent ones. |

**Count: 13 referenced, 3 justified deviations, 16 unreferenced in config, plus 5 more
that live in `peptide.py` outside config entirely — 21 unreferenced in total.**

Two parameters left this table on 2026-09-14 when S7d, the HyPhy BUSTED codon model on a
per-family FastTree tree, was removed from the analysis: `evolution.busted_max_p` (which
was referenced, to Murrell et al. 2015) and `evolution.min_members_for_busted` (which was
not). Removing a stage is the only way an unreferenced parameter leaves this table without
someone measuring something.

That ratio is the honest state of the pipeline. The unreferenced group is concentrated in the
cascade thresholds and the portfolio design — which is unsurprising, because those are the
parts with no direct precedent in the three source papers. They need either a supporting
citation from a different literature, or an explicit methods paragraph defending them as
design choices with measured sensitivity.

## Context rates (S8c)

Produced by `workflow/scripts/context_features.py`: per dark family, the fraction of its
plasmids on which a member has each context feature (`cons_*`). The enrichment test that
used to sit here - Fisher's exact test, Benjamini-Hochberg correction, a stratified
background (`background.covariates`, `background.min_stratum_size`) and the label-category
layer - was removed on 2026-09-25, and its parameters with it.

| parameter | value | source |
|---|---|---|
| unit of observation | plasmid | A family's members are homologs and its plasmids are frequently near-identical. Counting members makes sequencing effort look like evidence: one clinical plasmid sequenced forty times is forty members and one observation. |

**Outstanding.** The plasmid unit removes copy-number inflation *within* a plasmid but not
clonal redundancy *between* plasmids: forty independent depositions of the same clinical
plasmid remain forty units. `workflow/scripts/clonal_registry.py` holds the registry for
that correction; it is not applied to these rates.
