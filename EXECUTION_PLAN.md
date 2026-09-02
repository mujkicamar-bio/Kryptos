# Execution Plan — multireplicon assembly rules, ARG spread routes, environmental distribution

*Derived from `RESEARCH_OUTLINE.md` (themes A-D) and a 2026-06-30 audit of the actual downloaded IMG/PR data. This plan re-prioritizes the outline around three explicit asks: (1) complete plasmids only, (2) routes of ARG spread via transposons/IS and plasmid backbone, (3) the assembly rules and gene-content cost of multireplicon fusion (analyzed with SegMantX, per Hanke & Dagan 2025), plus (4) environmental enrichment of multireplicon plasmids. Other outline themes (B1/B3/B4, C1-C3, D) are kept as secondary/stretch work, not dropped.*

## 0. Confirmed data state (2026-06-30 audit)

- `data/Cus_PR/IMG_VR_2023-08-08_1/IMGPR_nucl.fna` — 699,973 sequences (`grep -c '^>'` confirmed). 17 GB uncompressed.
- `data/Cus_PR/IMG_VR_2023-08-08_1/IMGPR_plasmid_data.tsv` — 699,973 data rows, confirmed reconciling with FASTA count. Columns: `plasmid_id, ptu, taxon_oid, scaffold_oid, source_type, ecosystem, length, gene_count, genomad_score, putatively_complete, topology, mob_genes, t4cp_genes, t4ss_atpase_genes, other_conjugation_genes, complete_mpf_family, origin_of_transfer, arg_genes, putative_phage_plasmid, host_prediction_method, host_taxonomy, closest_reference, closest_reference_ani_percent, closest_reference_af_percent`.
- **No protein FASTA downloaded** (`IMGPR_prot.faa.gz` absent). **No lat/long coordinates** in metadata — only the categorical, semicolon-delimited GOLD `ecosystem` field (e.g. `Engineered;Bioreactor;Anaerobic;Sludge`). **No replicon-type column** — must be derived.
- Useful pre-computed columns we do NOT need to re-derive: `arg_genes` (Resfams), `mob_genes`/`t4cp_genes`/`t4ss_atpase_genes`/`complete_mpf_family` (CONJscan mobility class), `putatively_complete`, `topology`, `ptu`.
- Compute: SLURM available (`sbatch`/`squeue`/`sinfo`), 48 cores / 491 GB RAM confirmed on this node.

## Decisions locked in (2026-06-30)

1. **Geography:** ecosystem-only. No coordinate recovery from IMG/M/GOLD. Outline's B3 (spatial autocorrelation/Moran's I) is dropped; "environmental spread" is operationalized entirely via the categorical `ecosystem` field.
2. **Proteins:** do not block on the JGI protein FASTA. Call our own ORFs with **Prodigal** on the complete-plasmid nucleotide subset where a tool below doesn't already supply gene calls, decoupling protein-dependent analyses from needing a fresh JGI token. The JGI `IMGPR_prot.faa.gz` remains an optional later cross-validation step if a token becomes available.
3. **Tooling (user-specified, 2026-06-30):**
   - **MOB-suite / `mob_typer`** (Robertson & Nash 2018, already in `papers_of_interest.txt`) is the primary replicon + relaxase (MOB family) + predicted mobility typing tool — **not** PlasmidFinder. MOB-suite's database is broader than PlasmidFinder's Enterobacteriaceae-centric set but is still reference-based, so the "report untypeable %" caveat still applies.
   - **PlasAnn** is the annotation tool for the functional/gene-level annotation needed across phases (replicon/AMR/mobility/feature annotation). Exact install path and capability surface are unconfirmed — first use is a Phase 0 verification task, not an assumption.
   - **SegMantX** (Hanke & Dagan 2025, *Mol Biol Evol* 42(10):msaf242, "a novel tool for detecting DNA duplications uncovers prevalent duplications in plasmids") is the core method for the **merger/fusion analysis of multireplicon plasmids** (Phase 2), used per the approach in that paper rather than the ad-hoc IS-junction heuristic originally sketched in the outline. IS/transposase density scanning is retained, but downgraded from "primary fusion-detection method" to "complementary signature feeding the ARG/Tn spread-route question in Phase 3."
4. **This document is a plan, not yet execution.** No code/data was changed by this audit beyond read-only inspection.

## Decisions locked in (2026-06-30, post-Phase-0)

5. **Exclude "Simulated communities" from all ecosystem-based analyses.** Phase 0.4 QC found that
   `Engineered;Modeled;Simulated communities (contig mixture)` — an IMG/PR in-silico contig-mixture
   category, not a sampled environment — accounts for **51.44%** of the complete-plasmid subset
   (79,563 / 154,680). Combined with blank/"Unclassified" (16.42%), only **32.14%** of complete
   plasmids carry a genuine, usable ecosystem label. Every Phase 4 analysis (niche-breadth,
   ecosystem-stratified ARG enrichment in Phase 3.4/3.5, etc.) must filter to
   `is_blank_or_unclassified == 0 AND is_simulated_community == 0` in
   `data/processed/ecosystem_levels.tsv` before computing any ecosystem-based statistic — pooling
   the simulated-community rows in with real environmental samples would skew enrichment/niche
   results toward an artificial category. This shrinks the Phase 4 working set from 154,680 to
   ~49,718 complete plasmids; report this denominator explicitly in every Phase 4 result.
6. **Every report must cite the script(s)/pipeline that produced it.** Each report in `reports/`
   states, near the top, the exact script path(s) and tool invocation(s) used, not just prose
   findings — so any number can be regenerated or audited without reverse-engineering which script
   wrote it.

## Scope ordering (why)

The outline's own recommendation was "lead with B + C1/C2." This plan instead leads with **Theme A (assembly rules + fusion cost, via SegMantX) and the ARG/Tn/backbone question**, because that is what was explicitly asked for; ecosystem-only B2 and the small-plasmid-as-monomer angle (C4) are pulled in as direct support for that question rather than treated as a separate track. B1/B3/B4 (generalist/specialist mapping, AMR reservoirs) and C2/C3 (cargo, hitchhiking) and Theme D (D1/D2/D5) are kept as later phases / explicit stretch goals.

---

## Phase 0 — Foundation (no external blockers)

| # | Task | Output | Verify |
|---|---|---|---|
| 0.1 | Build analysis env (conda/mamba): `mob_suite` (mob_typer), `seqkit`, `csvtk`, `prodigal`, `hmmer`, `mmseqs2` or `blast+`, Pfam-A HMM DB, Python (`pandas`, `scipy`, `statsmodels`, `networkx`, `seaborn`); locate and install **PlasAnn** and **SegMantX** (confirm source repo, dependencies, license, install method — do not assume) | `environment.yml` + `tools/README.md` documenting where PlasAnn/SegMantX came from | `conda env export` reproducible; PlasAnn and SegMantX run successfully on a 5-plasmid test set |
| 0.2 | Extract `putatively_complete == "Yes"` subset from nucleotide FASTA + metadata | `data/processed/complete_plasmids.{fna,tsv}` | row count matches outline's stated 154,680 (or current true count — re-verify, don't assume the paper's number) and FASTA headers == TSV `plasmid_id` 1:1 |
| 0.3 | QC metadata completeness: % non-blank/non-"Unclassified" `ecosystem`, % with `arg_genes`, % with `mob_genes`/`t4cp_genes`/`t4ss_atpase_genes`, % with `host_taxonomy` | `reports/metadata_coverage.md` | numbers reported before any downstream claim (outline caveat §4) |
| 0.4 | Parse hierarchical `ecosystem` string into levels (Domain;Type;Subtype;...); build a clean categorical table | `data/processed/ecosystem_levels.tsv` | spot-check 20 random rows by hand |
| 0.5 | Call ORFs on the complete-plasmid subset with Prodigal (for any analysis PlasAnn doesn't already cover) | `data/processed/complete_plasmids.{faa,gff}` | per-plasmid ORF count correlates with IMG/PR's own `gene_count` column (sanity check, not exact match expected) |

**Dependencies:** none. **Risk:** PlasAnn's exact capabilities/install path are unconfirmed at planning time — 0.1 must resolve this before Phase 1 can rely on it; if PlasAnn turns out not to cover something we need, fall back to the Prodigal+HMM path already planned for 0.5/1.2.

---

## Phase 1 — Replicon typing (everything downstream depends on this)

| # | Task | Output | Verify |
|---|---|---|---|
| 1.1 | Primary replicon + relaxase + mobility typing via **MOB-suite `mob_typer`** on complete-plasmid nucleotide subset | `replicon_calls_mobtyper.tsv` | spot-check known Inc-type plasmids resolve correctly; cross-check predicted mobility against IMG/PR's own `mob_genes`/`t4cp_genes`/`t4ss_atpase_genes`/`complete_mpf_family` columns for agreement |
| 1.2 | De-novo Rep-protein HMM layer (Pfam Rep_1/Rep_2/Rep_3/RepA_N/Rep_trans/RepL etc.) against Phase 0.5 proteins, to catch replicons outside MOB-suite's reference database | `replicon_calls_hmm.tsv` | manually inspect a sample of HMM-only (MOB-typer-negative) hits for plausibility |
| 1.3 | **PlasAnn** annotation pass on complete-plasmid subset (replicon/AMR/mobility/feature calls, scope depends on 0.1 findings) — cross-validate against 1.1/1.2 | `plasann_calls.tsv` | agreement rate with MOB-typer reported; discrepancies reviewed, not silently dropped |
| 1.4 | Merge 1.1+1.2+1.3 per plasmid → `replicon_count`, `replicon_types[]`, `typed`/`untyped` flag | `data/processed/replicon_typed.tsv` | **report untypeable % explicitly** (outline caveat §2) — this number is itself a result, not just a caveat |
| 1.5 | Classify single- vs multi-replicon; cross-check prevalence against de Quinto's reported >30% | `reports/multireplicon_prevalence.md` | = outline **Milestone 2**: result falls in a defensible, explainable range |

**Dependencies:** Phase 0 complete (incl. tool verification). **Risk:** the untypeable fraction may be large for environmental/uncultured plasmids, shrinking the classifiable set — treat this explicitly as a finding ("how much of the environmental plasmidome is dark to current typing tools"), not a failure. **Risk:** MOB-typer, HMM layer, and PlasAnn may disagree on some calls — define an explicit reconciliation rule (e.g. majority vote, or MOB-typer + HMM as primary with PlasAnn as cross-validation) before merging.

---

## Phase 2 — Multireplicon assembly rules & cost of fusion (lead theme — user's primary ask #1)

| # | Task | Maps to outline | Output |
|---|---|---|---|
| 2.1 | Replicon co-occurrence/fusion network on complete multireplicon plasmids (edges = co-occurring replicon pairs, weighted by frequency) | A2 | `network_replicon_fusion.gml` + figure |
| 2.2 | Attempt to obtain de Quinto et al.'s reported "repeatedly-fusing" pairs (preprint supplementary; flagged in memory as Cloudflare-blocked previously) for a quantitative Jaccard/network comparison; if unobtainable, fall back to qualitative comparison against the abstract-reported pattern | A2 | comparison note in `reports/` |
| 2.3 | **SegMantX duplication/merger analysis** (Hanke & Dagan 2025 methodology) on complete multireplicon plasmids: detect within-plasmid DNA duplications consistent with replicon-fusion events; localize duplicated segments relative to replicon boundaries | new, primary mechanism for "rules of assembly" | `segmantx_duplications.tsv` + figure |
| 2.4 | IS/transposase density at replicon-junctions (Pfam-scan Phase 0.5/1.3 proteins for transposase domains, e.g. DDE_Tnp_*), as a **complementary** signature alongside SegMantX duplications | A4 | `is_junction_density.tsv` — **feeds Phase 3 (ARG/Tn spread)** |
| 2.5 | **"Cost of fusion" gene-content analysis.** For each replicon type, compare gene content (Pfam categories: partition genes `parA`/`parB`, toxin–antitoxin domains, replication-initiator/regulator families, redundant mobility genes) between (a) that replicon as the *sole* replicon on a standalone small plasmid and (b) the same replicon as *one component* of a multireplicon plasmid. Use SegMantX output (2.3) to flag degraded/partial duplicate copies as candidate "remnants of what was lost" during fusion resolution. Hypothesis: redundant maintenance machinery is systematically depleted in the fused state. | directly answers user's "what is lost when two single-replicon plasmids merge" | `reports/fusion_gene_loss.md` — table of Pfam categories with depletion/enrichment odds ratios + significance |
| 2.6 | Fusion-intermediate detection: combine completeness/topology + IS-junction (2.4) + SegMantX duplication (2.3) signatures to flag putative fusion/deletion intermediates at population scale | D4 | candidate list + manual review of top examples |
| 2.7 | Host-range expansion test within PTUs: do multireplicon PTU members reach broader `host_taxonomy` breadth than single-replicon relatives of the same PTU? | D3 | `reports/host_range_by_replicon_count.md` |

**Dependencies:** Phase 1 (+ Phase 0.1 SegMantX install confirmed). **Risk:** 2.5 needs a curated TA-system/partition-gene Pfam list — start with freely available Pfam-A domains (known TA and par Pfams) rather than a licensed TA database; scope can expand later. **Risk:** 2.2's de Quinto pair list may not be recoverable — don't block the phase on it.

---

## Phase 3 — ARG spread routes via Tn/IS and plasmid backbone (lead theme — user's primary ask #2)

| # | Task | Maps to outline | Output |
|---|---|---|---|
| 3.1 | Use existing `arg_genes` column directly (Resfams) — no re-calling needed; cross-check coverage against PlasAnn's AMR calls if available | — | filtered ARG-positive complete-plasmid table |
| 3.2 | Define "backbone" per plasmid = replicon type(s) (Phase 1, MOB-typer-led) + mobility class from `mob_genes`/`t4cp_genes`/`t4ss_atpase_genes`/`complete_mpf_family` (conjugative / mobilizable / non-mobilizable), cross-validated against MOB-typer's own mobility prediction | — | `data/processed/plasmid_backbone.tsv` |
| 3.3 | Test whether ARG genes are preferentially flanked by IS/transposase density (from 2.4) compared to non-ARG genes on the same plasmids — composite-transposon-style mobilization signature | A4, extended | `reports/arg_is_flanking.md` |
| 3.4 | Map ARG×backbone combinations across `ecosystem` categories (ecosystem-only, per locked decision) — identify backbone/replicon classes that act as cross-ecosystem ARG "vehicles" | B4 (simplified) | `reports/arg_backbone_ecosystem_map.md` + figure |
| 3.5 | Within-ecosystem stratified test: is ARG density still enriched in multireplicon vs single-replicon plasmids outside host-associated/clinical ecosystems? (avoids Simpson's paradox per outline caveat) | A3 | `reports/arg_enrichment_by_ecosystem.md` |

**Dependencies:** Phase 1 (backbone needs replicon typing) + Phase 2.4 (IS density). **Risk:** none major; this phase is mostly statistics over already-available + Phase 1/2 columns.

---

## Phase 4 — Environmental distribution of multireplicon plasmids (user's primary ask #3, ecosystem-only)

| # | Task | Maps to outline | Output |
|---|---|---|---|
| 4.1 | Replicon-type × ecosystem occurrence matrix; niche-breadth index (Shannon/Levins' B) per replicon type | B1 | first descriptive figure |
| 4.2 | Ecosystem-breadth of multireplicon vs single-replicon plasmids, controlling for plasmid size and per-ecosystem sampling effort (rarefaction) | B2 | `reports/multireplicon_ecosystem_breadth.md` |
| 4.3 | Statistical test: multireplicon prevalence by ecosystem category (logistic regression / chi-sq, sampling-effort corrected) | B2 | effect sizes + p-values table — directly answers "are multireplicon plasmids more present in certain environments" |
| 4.4 | **Fragments-included sensitivity check (secondary, explicitly labeled).** Re-run 4.1-4.3 on the broader set of *all* plasmids (complete + fragments) with a genuine ecosystem label (`is_blank_or_unclassified==0 AND is_simulated_community==0` in `data/processed/ecosystem_levels.tsv`), to test whether 4.1-4.3's complete-only conclusions are an artifact of the small (~49,718-plasmid) complete+genuine-ecosystem working set. **Replicon typing on fragments is for ecosystem/prevalence breadth only — never for fusion-mechanism, gene-loss, or architecture claims** (a fragment can look single-replicon when it's actually part of a multireplicon plasmid split across contigs; RESEARCH_OUTLINE.md caveat §4.1). Report 4.1-4.3 (complete-only) as the primary result in all cases; report 4.4 as a clearly-labeled robustness/sensitivity appendix, never merged into or substituted for the primary numbers. | B2 (extended) | `reports/ecosystem_sensitivity_fragments.md` |

**Dependencies:** Phase 1 (and Phase 1's replicon-typing pipeline applied to the fragment set for 4.4). Can run in parallel with Phase 2/3. **Decision locked 2026-06-30** (user-confirmed): complete-only stays primary for every Phase 4 result; the fragments-included run (4.4) is an explicitly-labeled secondary sensitivity check only.

---

## Phase 5 — Small plasmids as fusion building blocks (supports Phase 2's "what merges")

| # | Task | Maps to outline | Output |
|---|---|---|---|
| 5.1 | Size spectrum of complete plasmids by `source_type`; define small-plasmid threshold (<10 kb, per outline) | C1 | size-distribution figure |
| 5.2 | **C4 — direct mechanistic test:** do replicon types found standalone in small single-replicon plasmids recur, enriched, as components of multireplicon plasmids? (Fisher/hypergeometric enrichment test) | C4 | `reports/replicon_building_blocks.md` — direct evidence for/against "small plasmids are the monomers" |
| 5.3 *(stretch)* | Functional/cargo content of small plasmids (TA systems, defense systems) vs "truly cryptic" fraction | C2 | only if time remains |
| 5.4 *(stretch)* | Same-sample/same-host co-occurrence of mobilizable small plasmids with conjugative plasmids ("hitchhiking") | C3 | only if time remains; needs `taxon_oid`/`scaffold_oid` linkage logic, nontrivial |

**Dependencies:** Phase 1 (+ Phase 0.5/1.3 proteins for 5.3).

---

## Phase 6 — Synthesis

| # | Task |
|---|---|
| 6.1 | Cross-theme figures: fusion network with SegMantX duplication + IS/Tn annotation; ecosystem-stratified ARG/multireplicon enrichment; gene-loss-on-fusion summary; small-plasmid building-block enrichment |
| 6.2 | Caveats write-up: untypeable %, ecosystem metadata coverage %, sampling-effort corrections applied, ecosystem-only geography limitation explicitly stated, MOB-typer/HMM/PlasAnn reconciliation method stated |
| 6.3 *(optional)* | If a JGI token becomes available later: download `IMGPR_prot.faa.gz`, cross-validate Prodigal-based protein calls against IMG/PR's official gene calls as a robustness check |

---

## Stretch / explicitly deprioritized (outline Theme D minus D3/D4 already folded into Phase 2)

- **D1** Plasmid "guild" ecology (bipartite replicon×ecosystem community detection) — natural extension of Phase 4.1, pursue only after Phases 1-5 land.
- **D2** Defense-system cargo as an alternative spread driver — needs defense-system HMMs (PADLOC/DefenseFinder-style); follow-up to Phase 5.3.
- **D5** Metatranscriptome in-situ activity — IMG/PR has metatranscriptome-derived plasmids (`source_type`); genuinely novel but out of scope for this pass.

---

## Risks summary

1. **PlasAnn and SegMantX are unconfirmed installs at planning time** — Phase 0.1 must resolve source/install/capability before Phases 1-2 can rely on them. If either tool is unavailable or doesn't do what's expected, fall back to the Prodigal+Pfam-HMM path already planned as a parallel track.
2. Untypeable replicon fraction may be large for environmental plasmids (Phase 1) — report as a finding, not just a limitation.
3. De Quinto's exact fusion-pair list may be unobtainable (Phase 2.2) — don't block on it.
4. TA-system/partition-gene Pfam coverage for the fusion-cost analysis (Phase 2.5) may be incomplete — start with freely available Pfam-A, expand if needed.
5. `ecosystem` field is hierarchical free text with likely blank/"Unclassified" entries — must be parsed and coverage-reported before use (Phase 0.3/0.4).
6. Self-called Prodigal ORFs vs IMG/PR's own `gene_count` may disagree somewhat — cross-check as QC (Phase 0.5), not a blocker.
7. MOB-typer, HMM layer, and PlasAnn replicon/mobility calls may disagree — define an explicit reconciliation rule before merging (Phase 1.4).
8. All HMM/network/SegMantX jobs on ~150k+ sequences should be SLURM-batched (48 cores/491GB available on this node), not run interactively in a single session.

## Success criteria (verifiable milestones)

- **M1** (done): nucleotide FASTA count == metadata row count == 699,973.
- **M1b**: complete-plasmid subset count reconciles 1:1 between FASTA and TSV.
- **M2**: multireplicon prevalence on complete plasmids reported with untypeable % stated; falls in a defensible, explainable range vs. de Quinto's 30%.
- **M3** (Phase 2.5): quantified table of which gene categories are depleted/enriched in fused vs. standalone replicon contexts, with effect sizes — the direct answer to "what is lost," grounded in SegMantX duplication evidence.
- **M4** (Phase 4): replicon×ecosystem matrix + niche-breadth ranking; multireplicon ecosystem-breadth test with effect size + p-value.
- **M5** (Phase 3.4): table of ARG×backbone combinations recurring across ≥2 distinct ecosystem categories — candidate spread routes.
- **M6** (Phase 5.2): enrichment statistic (odds ratio + CI) for small-plasmid standalone replicons recurring as multireplicon components.
