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
| **ECLIPSE** | Lata S. & Heinz D.W. ECLIPSE: exploring the dark proteome of ESKAPE pathogens. *Bioinformatics* **42**(8) (2026) |
| **Pfam** | Mistry J. *et al.* Pfam: The protein families database in 2021. *Nucleic Acids Res.* **49**, D412–D419 |
| **AntiFam** | Eberhardt R.Y. *et al.* AntiFam: a tool to help identify spurious ORFs in protein annotation. *Database* **2012**, bas003 |
| **DefenseFinder** | Tesson F. *et al.* Systematic and quantitative view of the antiviral arsenal of prokaryotes. *Nat. Commun.* **13**, 2561 (2022) |
| **HMMER** | Eddy S.R. Accelerated profile HMM searches. *PLoS Comput. Biol.* **7**, e1002195 (2011) |

---

## Referenced — value matches a published use

| parameter | value | source | what they used it for |
|---|---|---|---|
| `clustering.min_seq_id` | 0.30 | FESNov | MMseqs2 deep-homology clustering of 400M genes |
| `clustering.coverage` | 0.50 | FESNov | same |
| `evolution.min_members_for_dnds` | 3 | FESNov | ">=3 complete genes" required per family |
| `evolution.dnds_purifying_max` | 0.5 | FESNov | purifying-selection filter on novel families |
| `evolution.rnacode_max_p` | 0.05 | FESNov | RNAcode coding-potential test. Now actually applied, on both strands. |
| `context.max_operon_gap` | 100 nt | FESNov | "absence of intergenic regions >100 nucleotides" |
| `context.neighbourhood_window` | 3 | FESNov | neighbouring genes at +/-3 positions |
| `context.min_context_conservation` | 0.50 | FESNov | 52,793 families annotated at >=50% confidence |
| `context.high_confidence_conservation` | 0.90 | FESNov | 4,349 families at >=90% confidence |
| `quality_gate.min_control_recall` | 0.99 | ECLIPSE | 99.2–100% of 246 virulence / 42 AMR / 75 essential genes recovered as annotated |
| T2, T3, T4 `max_evalue` | 1e-5 | FESNov | AntiFam and pVOG screening threshold |
| T1 `--cut_ga` | — | Pfam | curated per-family gathering thresholds define family membership |
| DefenseFinder `inter_gene_max_space` etc. | per model | DefenseFinder | shipped in the 711 system definition XMLs; referenced by construction |

## Deviating from a paper, with a measurement

| parameter | our value | paper value | why we deviate |
|---|---|---|---|
| `artefact_screen.antifam_args` | `--cut_ga` | FESNov used E <= 1e-5 | AntiFam ships a curated GA line for **all 278 profiles**. Measured directly from `data/refs/antifam/AntiFam.hmm` by converting each profile's GA to an E-value with its own `STATS LOCAL FORWARD` parameters at Z=3,497,616: **274 of 278 (98.6%) curated thresholds are LOOSER than E=1e-5**, the median curated cut corresponds to E=7.3e-4 (73x looser) and the loosest, `Spurious_ORF_67`, to E=0.195. A blanket 1e-5 floor therefore overrides the curator across essentially the whole artefact database, in the one stage that stops shadow ORFs reaching the plate. This is the same argument that already sets `--cut_ga` on T1, applied where the consequence is worse: a missed artefact is not a missed annotation, it is a non-protein sent to the bench. Regression test: `tests/test_scripts_smoke.py::test_the_artefact_screen_uses_antifams_curated_thresholds`. |
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
| `targets.REALITY_TESTS` MOB-cluster minimum | 2 | **Fixed.** Had no config counterpart at all. Now `prioritisation.min_mob_clusters`. |
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
| `max_target_seqs` | 5 | **No source.** |
| `sweep_cohort_fraction` | 0.02 | **No source.** A cost/benefit choice, not a scientific one — but still undeclared in the literature. |
| `max_low_complexity_fraction` | 0.5 | **No source.** |
| `evolution.min_codons` | 20 | **No source.** Now actually applied: `dnds_detail(..., min_codons=)` is passed the configured value per pair by S7b, and the arithmetic floor of 3 remains underneath it. Until that wiring existed the declared 20 was inert and an 8-codon fragment could fire `purifying_selection`, the strongest of the four reality tests. |
| `evolution.max_members_aligned` | 50 | **No source.** A compute cap; measured worst-case family 1.95 s. |
| `context.min_enrichment` | 2.0 | **No source.** |
| `structure.max_evalue` | 1e-3 | **No source.** |
| `prioritisation.min_reality_lines` | 2 | **No source.** Sensitivity across 1–4 is reported, which shows its cost but does not justify the choice. |
| `portfolio.strata` quotas | 200/175/175/175/125/75/75 | **No source.** A design judgement about experimental portfolio balance; arguably not the kind of number a paper can supply, but it must be defended explicitly in the methods. |
| `library.length_liability_above_aa` | 400 | **No source.** |
| `n_cascade_shards` | 64 | **No source, and not the kind of number a paper supplies** — it is an engineering trade-off. Finer shards reduce the cost of a failure and widen parallelism; they also multiply the fixed cost of reading the search database, because DIAMOND streams the whole of nr per invocation. Measured side: T1 is 65.8 s fixed + 0.0365 s/protein at 4 threads, so 64 shards is ~34 min per shard. **Unmeasured side: the nr fixed cost.** `workflow/bench_nr.sbatch` exists to measure it and needs a project allocation. |
| `hmmer_z` | 3,497,616 | Not a threshold — it is the analysis-set size, and pinning it is what makes E-values comparable across shards (HMMER). The *practice* is standard; the specific value is simply our data. |

## Declared procedures that are not parameters

Two things in S9 look like tuning and are not. They are recorded here so a reader does not
go hunting for a citation that cannot exist.

| item | what it is |
|---|---|
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

## Context enrichment significance (S8c)

Produced by `src/plasmidann/enrich.py`, consumed by `workflow/scripts/context_features.py`.

| parameter | value | source |
|---|---|---|
| test | Fisher's exact, one-sided (`alternative="greater"`) | Exact rather than chi-squared because many categories are rare and expected cell counts fall well below the 5 the chi-squared approximation requires. One-sided because the hypothesis is over-representation; a category a dark family avoids is not a screening hypothesis, and a two-sided test would spend half its power looking for one. |
| multiple-testing correction | Benjamini-Hochberg FDR, across the categories tested per family | Benjamini and Hochberg 1995, *J R Stat Soc B* 57:289. FDR rather than family-wise error because these values rank candidates for a 1,000-construct screen: a false positive costs a well, a false negative costs a discovery. Correction is needed at all only because the label vocabulary is open - the six hand-picked features it replaces did not need it. |
| unit of observation | plasmid | A family's members are homologs and its plasmids are frequently near-identical. Counting members makes sequencing effort look like evidence: one clinical plasmid sequenced forty times is forty members and one observation. |
| `min_units` | 2 | One plasmid is an anecdote. Below two units no test is performed and the status is `TOO_FEW_MEMBERS`, so an untested family cannot rank beside a tested one. |
| background | the rest of the corpus, not the whole corpus | Testing a family against a background it contributes to shrinks any real effect, and for a family covering most of the corpus it shrinks it to nothing. |

**Outstanding.** The plasmid unit removes copy-number inflation *within* a plasmid but not
clonal redundancy *between* plasmids: forty independent depositions of the same clinical
plasmid remain forty units. `workflow/scripts/clonal_registry.py` holds the registry for
that correction; it is not applied here, and it should be before the enrichment values are
quoted in a manuscript.
