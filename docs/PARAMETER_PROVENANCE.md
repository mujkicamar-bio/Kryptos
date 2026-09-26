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
| **PlasAnn** | Islam H., Sharma A., Blair J. & Lopatkin A.J. PlasAnn: a curated plasmid-specific database and annotation pipeline for standardized gene and function analysis. *Nucleic Acids Res.* **54**(3), gkaf1507 (2026), doi:10.1093/nar/gkaf1507. Cited for its identity and coverage tiers only; its database, labels and tool are not used |
| **CARD** | Alcock B.P. *et al.* CARD 2023: expanded curation, support for machine learning, and resistome prediction at the Comprehensive Antibiotic Resistance Database. *Nucleic Acids Res.* **51**, D690–D699 (2023), doi:10.1093/nar/gkac920 |
| **RGI** | the Resistance Gene Identifier source code, github.com/arpcard/rgi: `app/Diamond.py` (DIAMOND `--more-sensitive`) and `app/HomologModel.py` (Strict at bit score `>=` the model cut-off) |
| **AMRFinderPlus** | Feldgarden M. *et al.* AMRFinderPlus and the Reference Gene Catalog facilitate examination of the genomic links among antimicrobial resistance, stress response, and virulence. *Sci. Rep.* **11**, 12728 (2021), doi:10.1038/s41598-021-91456-0 |
| **TADB** | Guan J. *et al.* TADB 3.0: an updated database of bacterial toxin–antitoxin loci and associated mobile genetic elements. *Nucleic Acids Res.* **52**, D784–D790 (2024), doi:10.1093/nar/gkad962 |
| **BacMet** | Pal C. *et al.* BacMet: antibacterial biocide and metal resistance genes database. *Nucleic Acids Res.* **42**, D737–D743 (2014), doi:10.1093/nar/gkt1252 (no separate paper for version 2.0) |
| **oriTDB** | Liu G. *et al.* oriTDB: a database of the origin-of-transfer regions of bacterial mobile genetic elements. *Nucleic Acids Res.* **53**, D163–D168 (2025), doi:10.1093/nar/gkae869 |
| **mobileOG-db** | Brown C.L. *et al.* mobileOG-db: a manually curated database of protein families mediating the life cycle of bacterial mobile genetic elements. *Appl. Environ. Microbiol.* **88**(18), e00991-22 (2022), doi:10.1128/aem.00991-22 |
| **dbAPIS** | Yan Y., Zheng J., Zhang X. & Yin Y. dbAPIS: a database of anti-prokaryotic immune system genes. *Nucleic Acids Res.* **52**, D419–D425 (2024), doi:10.1093/nar/gkad932 |
| **Anti-CRISPRdb** | Dong C. *et al.* Anti-CRISPRdb v2.2: an online repository of anti-CRISPR proteins including information on inhibitory mechanisms, activities and neighbors of curated anti-CRISPR proteins. *Database* **2022**, baac010 (2022), doi:10.1093/database/baac010 |
| **CONJScan** | Cury J. *et al.* Identifying conjugative plasmids and integrative conjugative elements with CONJscan. *Methods Mol. Biol.* **2075**, 265–283 (2020), doi:10.1007/978-1-4939-9877-7_19 |
| **Benjamini-Hochberg** | Benjamini Y. & Hochberg Y. Controlling the false discovery rate: a practical and powerful approach to multiple testing. *J. R. Stat. Soc. B* **57**, 289–300 (1995), doi:10.1111/j.2517-6161.1995.tb02031.x |
| **Coluzzi** | Coluzzi C., Garcillán-Barcia M.P., de la Cruz F. & Rocha E.P.C. Evolution of plasmid mobility: origin and fate of conjugative and nonconjugative plasmids. *Mol. Biol. Evol.* **39**(6), msac115 (2022), doi:10.1093/molbev/msac115 |
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
| `context.min_context_conservation` (was 0.50) and `context.high_confidence_conservation` (was 0.90), cited to FESNov | Removed from config on 2026-09-21 (commit cfffef3) because nothing read them; their rows stayed in this document until 2026-09-25. FESNov's two confidence levels, 50% and 90%, are used again only as the precision targets of `tools/calibrate_context.py` (see "Context terms and their calibration" below), which runs after the pipeline and sets no value inside it. |
| `peptide.AMP_MAX_LENGTH` (100), `AMP_MIN_CHARGE` (2.0), `AMP_MIN_HYDROPHOBIC_FRACTION` (0.3), `TM_WINDOW` (19), `TM_THRESHOLD` (1.6); `prioritisation.min_reality_lines` (2); `portfolio.strata` quotas; `library.length_liability_above_aa` (400); `targets.RANK_PRIORITY` | Removed with Layer C (experimental prioritisation, the portfolio and library design) on 2026-09-17 (commit efca5c3; spec section 76). `peptide.py`, `targets.py`, `prioritise.py` and `library_design.py` no longer exist, so none of these values reaches any output. Their rows remained below, as live and unreferenced, until 2026-09-25. |

## Thresholds that live in Python, not config

P4 says every threshold is declared in config, schema-validated, and stamped into the row it
governs. These are the exceptions, found by review. The reality-test thresholds have since
been moved into config or derived from it; the peptide constants that used to be listed here
went with Layer C (table above).

| parameter | value | status |
|---|---|---|
| `evidence.REALITY_TESTS` dN/dS cutoff | 0.5 | **Fixed.** Was a literal in a lambda that duplicated `evolution.dnds_purifying_max`. Now read from that config value by `evidence.reality_thresholds`, so the two cannot diverge. |
| `evidence.REALITY_TESTS` member minimum | 3 | **Fixed.** Was a literal duplicating `evolution.min_members_for_dnds`. `evidence.reality_thresholds` now takes it from that config value, so the `purifying_selection ⇒ is_family` nesting that `IMPLIED_BY` relies on cannot drift apart (`check_reality_config`, which used to refuse a config where they differed, went with Layer C). |
| `evidence.reality_thresholds` lineage minimum (`min_lineages`) | 2 | Definitional: "more than one". Counted over Stage 6 Mash lineages since 2026-09-25; it was a MOB-suite cluster minimum, and MOB-suite assigns the nearest reference's cluster however distant. |
| `rarity.cross_min_hosts` / `rarity.cross_min_genera` | 2 / 2 | Definitional: "more than one" observed species / genus (were 5 and 3, no source). `SINGLE_MOB`/`CROSS_MOB` are fixed at 1 and >= 2 MOB-suite clusters; `cross_min_mob` was removed. |

## UNREFERENCED — known weakness, must be resolved

These have no published source. Each is a researcher degree of freedom until it gets one.

| parameter | value | status |
|---|---|---|
| `narrow_at` | 0.7 | **No source. User decision, 2026-09-25** (was 0.9). A protein stops at the tier where it reaches 0.7 explained. Measured sensitivity: in v1, moving the single threshold from 0.3 to 0.8 changed the deep-tier set by 76%, so the value matters. Consequences: `min_explained` can be swept up to 0.7 without a re-run, not beyond; a protein 0.7-1 explained after an early tier gets no deeper-tier result and can stay PARTIAL where a deeper tier would have made it FULL (`full_at` 0.8). The 2% sweep cohort bypasses narrowing and measures what stopping at 0.7 costs, which is a partial substitute for a citation but not a replacement. `check_thresholds` requires `narrow_at >= min_explained`. Set in `config/cascade.yaml` and its test and benchmark counterparts, the only files the workflow reads it from. |
| `min_explained` | 0.5 | **No source.** Applied post hoc and sweepable, so its cost is measurable — but the value itself is unjustified. |
| `min_coverage` | 0.5 | Partial: 50% coverage appears in FESNov and ECLIPSE for *clustering*, not for FUNCTIONAL/DOMAIN_ONLY classification. Not the same use. |
| `full_at` / `partial_at` | 0.8 / 0.5 | **No source.** |
| `max_target_seqs` | 25 | DIAMOND's default. Raised from 5 (no source) after 44% of nr benchmark queries filled all five slots with uninformative hits; the limit also changes DIAMOND's search, not only its reporting (Shah et al., *Bioinformatics* **35**, 1613 (2019)). |
| `sweep_cohort_fraction` | 0.02 | **No source.** A cost/benefit choice, not a scientific one — but still undeclared in the literature. |
| `max_low_complexity_fraction` | 0.5 | **No source.** |
| `evolution.min_codons` | 20 | **No source.** Now actually applied: `dnds_detail(..., min_codons=)` is passed the configured value per pair by S7b, and the arithmetic floor of 3 remains underneath it. Until that wiring existed the declared 20 was inert and an 8-codon fragment could fire `purifying_selection`, the strongest of the four reality tests. |
| `evolution.max_members_aligned` | 50 | **No source.** A compute cap; measured worst-case family 1.95 s. |
| `structure.max_evalue` | 1e-3 | **No source.** |
| `hmmer_z` | 3,498,616 | Not a threshold — it is every unique protein of the analysis set (3,497,616) plus the 1,000 controls, searched or not, and pinning it is what makes E-values comparable across tiers (HMMER). Proteins resolved at Tier 0 stay in it, so E-values do not depend on PlasmidScope's coverage. The *practice* is standard; the specific value is simply our data. |

## Declared procedures that are not parameters

The items below look like tuning and are not. They are recorded here so a reader does not
go hunting for a citation that cannot exist.

| item | what it is |
|---|---|
| S8e, ISEScan | Not a threshold of ours. ISEScan 1.7.3 is run with its published defaults (ISEScan), so partial elements are kept; `--removeShortIS` is not used, because a partial IS on a plasmid is still an IS-derived region. |
| Tier 0, PlasmidScope transfer | Not a threshold of ours. A protein identical in sequence to a PlasmidScope protein takes PlasmidScope's published eggNOG-mapper 2.1.12 result (PlasmidScope; eggNOG-mapper, default settings as published). It counts as annotated when that result has a KEGG KO, an EC number, or a Pfam family that passes `cascade.is_informative` (so not a DUF/UPF family) — the eggNOG fields as released, with no score cut-off applied by us. A COG/OG category letter alone does not count. |
| `evidence.IMPLIED_BY` | Which reality tests entail which others. `purifying_selection` requires a measured dN/dS, which requires >=3 members, which is `is_family` — so those two are one line of evidence, not two, and `reality_n` counts only the independent ones. |

**Count (rows of the tables above, recounted 2026-09-25): 21 referenced, 7 justified
deviations, 4 definitional or config-derived thresholds in Python, and 11 rows in the
unreferenced table - 8 with no source at all, `min_coverage` with a partial one,
`max_target_seqs` at DIAMOND's default, and `hmmer_z`, which is not a threshold. The count
recorded before (13 / 3 / 21) had not been kept up to date; the 5 `peptide.py` constants
and 3 Layer C settings it included were removed on 2026-09-17. The parameters added on
2026-09-25 are counted in their own sections below.**

Two parameters left this table on 2026-09-14 when S7d, the HyPhy BUSTED codon model on a
per-family FastTree tree, was removed from the analysis: `evolution.busted_max_p` (which
was referenced, to Murrell et al. 2015) and `evolution.min_members_for_busted` (which was
not). Removing a stage is the only way an unreferenced parameter leaves this table without
someone measuring something.

That ratio is the honest state of the pipeline. The unreferenced group is concentrated in the
cascade thresholds - which is unsurprising, because they are the part with no direct
precedent in the three source papers. They need either a supporting
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
that correction; it is not applied to these rates. The context terms below, and synteny,
are counted over Stage 6 lineages and so do correct it; the `cons_*` rates, `cons_conj`
included, do not.

## Plasmid label databases (S4d, added 2026-09-25)

Produced by `workflow/scripts/label_databases.py`; the rules are constants of
`src/plasmidann/labeldb.py`, cited there, not config settings. PlasAnn's database, labels
and tool are not used (spec section 24b); only its published tiers are cited.

| parameter | value | source |
|---|---|---|
| tier 1 (`TIER1_IDENTITY`, `TIER1_COVERAGE`) | identity >= 80% and coverage >= 90% | PlasAnn, Methods ("Annotation pipeline"): "≥80% identity and ≥90% coverage in the primary tier". Applied to TADB, BacMet, oriTDB, mobileOG-db, dbAPIS and Anti-CRISPRdb. PlasAnn applied the tiers to its own database, which it had corrected against TADB 3.0 and BacMet 2.0; here they are applied to each database directly |
| tier 2 (`TIER2_IDENTITY`, `TIER2_COVERAGE`) | identity > 60% and coverage > 70% | PlasAnn: "a secondary tier of >60% identity and >70% coverage for more divergent homologs" |
| coverage measured on | the query AND the subject (the smaller of the two) | **Deviation, stricter.** The paper does not say which sequence the coverage is of, and its released code (plasann 1.1.6) applies no coverage filter at all. Query coverage alone lets a fragment carry the label of a full-length reference; subject coverage alone lets a multidomain protein carry the label of one domain. Measured on the 100-plasmid test set against the PlasAnn database: 1,188 proteins labelled with both coverages, 1,267 with query coverage alone |
| DIAMOND E-value | 1e-5 | PlasAnn: its database curation used "e-value < 1 × 10–5", and its released code searches with `-evalue 1e-5`. The same value as the cascade's T2–T4 (FESNov) |
| DIAMOND sensitivity | `--more-sensitive`, for every database | RGI (`app/Diamond.py` passes `--more-sensitive`); one mode for every database so that the sources are searched with the same sensitivity |
| target limit | `--max-target-seqs 0`, pre-filtered at `--id 60 --query-cover 70 --subject-cover 70` | **Measured.** A target limit keeps the top hits by bit score, not by tier, so a partial high-scoring hit can push the tier-1 hit out: on the 100-plasmid test set, 50 targets lost 1 oriTDB and 3 mobileOG-db labels and gave 2 + 2 more the wrong tier, and even 1,000 targets truncated oriTDB queries. The pre-filter is the tier-2 minimum, so it removes no hit that could qualify |
| best hit per protein and database | best tier, then highest bit score, then reference order | a declared order, not a threshold: the tier is the statement a label is cited with, and the last key makes ties independent of hit order |
| CARD Perfect | 100% identity over the full length of the reference | CARD: the "'Perfect' algorithm detects perfect matches to the curated reference sequences"; RGI |
| CARD Strict | bit score >= the model's curated cut-off (`model_param.blastp_bit_score` in card.json) | CARD: the "'Strict' algorithm predicts variants of known ARGs ... using CARD's curated bit-score cut-offs"; RGI (`app/HomologModel.py` tests `>=`). Loose hits are discarded, as CARD's own resistome predictions keep only Perfect and Strict |
| CARD models used | protein homolog models only (6,059 in CARD 4.0.2) | variant, rRNA, overexpression and knockout models detect resistance from mutations or absence, which the presence of a similar protein cannot show |
| CARD DIAMOND settings | `--more-sensitive`, DIAMOND's defaults otherwise | RGI (`app/Diamond.py`); the curated cut-off, not the E-value, decides a Strict call |
| AMRFinderPlus | `amrfinder -p <proteins> --plus`, its own rules and curated cut-offs; 4.2.7, database 2026-08-07.1 | AMRFinderPlus: a Reference Gene Catalog with "manually curated cutoffs"; `--plus` adds the stress-response and virulence genes. Nothing is overridden |
| TADB entries | experimentally validated protein entries (`*_exp` files), 963 | TADB distributes validated and predicted entries separately; the predicted set is not used (plan decision 2026-09-25). The 114 RNA toxins and antitoxins are nucleotide entries and cannot be searched with proteins |
| BacMet entries | BacMet2_EXP, 753 | BacMet: genes "with experimentally confirmed function described in the scientific literature"; the predicted set is not used |
| oriTDB entries | relaxase, auxiliary protein and T4CP, validated plus predicted (`_all` files), 16,834 | oriTDB holds experimentally validated and predicted entries; both are used. oriTDB offers no T4SS protein download |
| mobileOG-db entries | every entry: Manual, Homology and Keyword Search, 775,257 | user decision 2026-09-25; the evidence class of the reference entry is written into each label's sub_label |
| dbAPIS entries | verified APIS proteins and their sequence homologues, 17,414; the anti-CRISPR entries of its FASTA are not installed | dbAPIS: "experimentally verified APIS genes ... sequence and structural homologs", and it excludes anti-CRISPRs, whose source is Anti-CRISPRdb. Evidence class in sub_label |
| Anti-CRISPRdb entries | every entry of the core dataset: Verified, PLiterature and Putative, 3,692 | Anti-CRISPRdb: Verified are "experimentally validated Acrs in literatures", PLiterature are reported without an experiment, and putative entries were "retrieved from prokaryotes via PSI-BLAST alignment". Evidence class in sub_label |
| AMRFinderPlus element types given a context term | AMR -> amr:; STRESS METAL and BIOCIDE -> metal:; STRESS ACID, HEAT and VIRULENCE -> none | a mapping onto the term vocabulary, not a threshold. An element type the mapping does not know stops the stage rather than being guessed |

`label_disagreements.tsv` compares gene names after removing case and punctuation, and
counts one name that is a prefix of the other of at least three characters as the same
gene (`labeldb.same_gene`). Three characters is the length of a bacterial gene-symbol stem
(tet, sul, mer); the code cites Demerec et al. 1966 (Genetics 54:61) for it, whose full text
could not be retrieved to re-verify for this document. The rule changes no label: it only
decides which pairs are listed for review.

## CONJScan and DefenseFinder (S8f, S8a, added 2026-09-25)

| parameter | value | source |
|---|---|---|
| CONJScan models and their quorums | CONJScan 2.1.0, `CONJScan/Plasmids all`, as shipped | CONJScan; Coluzzi. Referenced by construction: the definitions are the published models. The package recommends running all models of one set together |
| plasmid mobility class | pCONJ (a T4SS_type model), else pdCONJ (dCONJ_type), else pMOB (MOB), else pMOBless | Coluzzi: conjugative when a plasmid encodes "a presumably complete machinery for conjugation (relaxase and MPF)", decayed conjugative when "a relaxase but an incomplete MPF machinery", mobilisable when "a relaxase gene but no or very few MPF genes", and pMOBless when it lacks a relaxase |
| MacSyFinder version for CONJScan | >= 2.1.6 (envs/conjscan pins 2.1.6 with MacSyLib 1.0.4) | CONJScan README, changelog entry 2.0.2: the definitions use model grammar 2.1 and need "MacSyFinder >= 2.1.6 and MacSylib >= 1.0.4". MacSyFinder 2.1.4, which DefenseFinder pins, refuses them (measured) |
| hit selection (`--i-evalue-sel`) | 0.001, MacSyFinder's default, for CONJScan and DefenseFinder | MacSyFinder's default, not ours. HMMER's i-evalue scales with the number of sequences searched, so both stages now search ONE database: measured on the 100-plasmid test set, 8 per-core chunks called 214 ORFs in 56 CONJScan systems and one database 212 in 55. The E-values still depend on the size of the collection searched, as with any single MacSyFinder database |
| replicon topology | circular, as for DefenseFinder | measured: a per-replicon topology file gave identical CONJScan calls on the test set (212 of 212 ORFs), because the T4SS and dCONJ plasmid models allow up to 500 intervening genes |

## Synteny (Stage 9, changed 2026-09-25)

The occurrence statistic and its minimum of two occurrences are replaced; see spec section 42.

| parameter | value | source |
|---|---|---|
| `synteny.levels` | close (90% identity) and intermediate (`clustering.primary`) | the two clusterings already cited above: UniRef90 for close, the UniRef50 analogue of Durairaj et al. for the family. The primary level must be listed |
| unit of observation | Stage 6 lineage (Mash distance <= 0.05, single linkage), fractional vote: each voting lineage has weight 1, split equally over its values | a counting rule, not a threshold. Forty copies of one redeposited plasmid are one lineage and one vote |
| `synteny.min_lineages` | 2 | arithmetic minimum: conservation across one lineage is a single observation (spec section 2.9). Below it the status is TOO_FEW_LINEAGES |

## Dark family co-occurrence (S8g, added 2026-09-25)

Produced by `workflow/scripts/dark_cooccurrence.py` with the rules of
`src/plasmidann/cooccurrence.py`; configured in `config/targets.yaml` `cooccurrence`.

| parameter | value | source |
|---|---|---|
| together | a member ORF of each family on the same plasmid | a definition, not a threshold. Two families in one lineage on different plasmids are not together |
| unit of observation | Stage 6 lineage; N = the distinct lineages in `plasmid_lineage.tsv` | as for synteny: clonal copies of one plasmid are one observation |
| test | hypergeometric upper tail P(X >= k) over lineages (K lineages carry A, n carry B, k share a plasmid), in log space with `math.lgamma` | the standard test for the overlap of two sets drawn from one population; exact, so no approximation is tuned. Conservative relative to lineage-level overlap, because k counts shared plasmids. Checked against exact enumeration in `tests/test_cooccurrence.py` |
| `cooccurrence.min_lineages_together` | 2 | **Measured justification, no published source.** One shared lineage is one observation, and for two singleton families P(X >= 1) = 1/N, so every pair of singleton families on one plasmid would appear significant. The same minimum as `synteny.min_lineages` (spec section 2.9). Because the minimum depends on the outcome, Benjamini-Hochberg is over the tested set only |
| multiple-testing correction | Benjamini-Hochberg across the tested pairs | Benjamini-Hochberg |
| `cooccurrence.fdr` | 0.05 | the conventional false discovery rate, used only to count partners in the report; every q-value is written to `dark_cooccurrence.tsv`, so another level can be applied |

Measured on the 100-plasmid test set: 93 lineages, 6 pairs tested, all at q <= 0.05, each
with K = n = k = 2 and p = 1/C(93, 2) = 2.34e-4. At full scale, with N of the order of 1e5
lineages, two rare families together in two lineages get p = 1/C(N, 2); a significant pair
can therefore rest on two observations, and significance means genetic linkage or
co-transfer, not a shared function.

## Context terms and their calibration (S8c and a post-run tool, added 2026-09-25)

`family_context_terms.tsv` is written by `workflow/scripts/context_features.py` with the
rules of `src/plasmidann/context_terms.py`; `tools/calibrate_context.py` runs after the
pipeline and is not a rule.

| parameter | value | source |
|---|---|---|
| gene-label neighbour | within `context.neighbourhood_window` (3) genes AND in the ORF's directon: same strand, intergenic gaps <= `context.max_operon_gap` (100 nt) | FESNov: neighbours "on the same orientation as FESNov family members" and "with no intergenic regions longer than 100 nucleotides", at +/-3 positions |
| system neighbour | the ORF is a component, or a component lies within 3 genes on either strand | **Deviation from the FESNov strand rule, by definition of a system.** A DefenseFinder or CONJScan system is a co-localised multi-gene call whose components sit on both strands; requiring the ORF's strand would split one system into two |
| tandem paralogues | a neighbour in the focal family is excluded | a counting rule: such a neighbour's label says what the family is, not what surrounds it |
| `context_terms.MIN_LINEAGES` | 2 | arithmetic minimum, as for synteny. Below it the status is TOO_FEW_LINEAGES |
| unit of observation | Stage 6 lineage | as for synteny |
| calibration precision targets (`LEVELS`) | 0.5 and 0.9 | FESNov: thresholds were set on 108,823 families with known KEGG pathways at the levels that recover the function of "at least 50 or 90% of all families" (their confidence scores). The levels are taken; the measure differs - here it is precision, the fraction of known families at or above a conservation value that carry the term themselves |
| benchmark membership | a known family carries a term when more than half of its member proteins carry it | a majority, definitional. Families where a minority carry it are excluded from both sides |
| calibration `MIN_FAMILIES` | 10 benchmark families and 10 negatives per term, and 10 families at or above a threshold | arithmetic: 90% precision cannot be observed on fewer than 10 families, since one false positive in 9 is 8/9 = 88.9%. With no negatives every level is 100% precise by construction. A term below the minimum is UNCALIBRATED and gets no threshold |
| threshold choice | the lowest conservation whose precision reaches the level and stays there at every higher point holding at least 10 families | a declared rule, not a tuned value |

**Count for the sections added on 2026-09-25:** 19 label-database rows, 5 CONJScan and
DefenseFinder rows, 3 synteny rows, 6 co-occurrence rows and 9 context-term rows. None is
unreferenced: each has a published source (a paper, or the source code of the tool a paper
describes), a measurement, or is an arithmetic minimum (2 lineages, 10 families), a
counting rule or a mapping declared as such; `cooccurrence.min_lineages_together` is
justified by a measurement rather than a paper, and `cooccurrence.fdr` is the conventional
level. Two are deviations, both stricter or definitional and stated:
coverage on both sequences, and system terms on either strand.
