# Plasmid analysis — project overview & status

*Last updated 2026-07-27. This is the front-door document: what the project is, what has been
built so far, where everything lives, and how to reproduce it. Every table and column named here is
produced by a specific script, cited inline.*

---

## 1. What this project is

A large-scale characterization of **complete, non-redundant bacterial plasmids** drawn from public
databases, built up into a single per-plasmid **master table** carrying physical properties,
sampling environment, geography, replicon/Inc typing, mobility, annotated gene content, and
antibiotic-resistance (AMR) content. The central biological question driving the design is the
**multireplicon** question (how many plasmids carry ≥2 distinct replicons, and what distinguishes
them), with AMR-by-environment and AMR-by-mobility as the main downstream analyses.

The planning and rationale behind each decision live in the root planning documents
(`RESEARCH_OUTLINE.md`, `EXECUTION_PLAN.md`, `METHODOLOGY.txt`, and the per-phase `*_PLAN.md` files).
This overview summarises what was actually executed.

---

## 2. Current status (done vs. pending)

| Stage | Status | Central output |
|---|---|---|
| Phase 0 — tooling, environments, IMG/PR pilot | ✅ done | `tools/README.md`, envs |
| Source characterization (PlasmidScope `ALL`) | ✅ done | `reports/plasmidscope_primary_characterization.md` |
| Working-set assembly (completeness + lab-made gates) | ✅ done | `data/plasmidscope_primary/working_set.tsv` (208,248) |
| Environment recovery + locked habitat taxonomy | ✅ done | `reports/environment_taxonomy_LOCKED.md` |
| PLSDB curated metadata enrichment | ✅ done | `data/plasmidscope_primary/plsdb_enrichment.tsv` |
| PlasAnn functional annotation | ✅ done (99.97%) | `data/plasmidscope_primary/plasann_features.tsv` |
| Authoritative replicon/Inc + MOB-suite typing | ✅ done | `reports/typing_methodology.md` |
| CARD/RGI AMR annotation | ✅ done | `reports/card_amr_methodology.md` |
| **Master table assembled** | ✅ **67 columns × 208,248 plasmids** | `data/plasmidscope_primary/plasmid_metadata_master.tsv` |
| Phase 3 — AMR × One Health compartment × mobility | ✅ done | `reports/amr_onehealth_methodology.md` |
| Manuscript draft (Nucleic Acids Research) | ⚠️ draft, needs reframing | `manuscript/amr_onehealth_NAR.md` |
| Peer review + reassessment | ✅ done — thesis refuted, see below | `manuscript/review_round1.md`, `manuscript/reassessment_round1.md` |
| Phase 3 — ARG co-occurrence, Inc-resolved carriage, geography | ⏳ next | — |
| SegMantX duplication/segment analysis | ⏳ tool installed, not yet run | `tools/SegMantX/` |

The data layers are complete; what remains is the interpretive analysis on top of the master table.

---

## 3. The working set (inclusion funnel)

Built by `scripts/characterize_plasmidscope_all.py` (source characterization) and
`scripts/assemble_working_set.py` (assembly). PlasmidScope's `ALL` table is already deduplicated
(MMseqs2, 100% identity & coverage), so there are **no duplicate plasmids across source databases**.

| stage | n | note |
|---|---:|---|
| PlasmidScope `ALL` (deduplicated) | 852,600 | source universe |
| complete / closed only | 208,360 | −644,240 incomplete/unknown-completeness |
| − lab-made / synthetic | −112 | logged in `discarded_labmade.tsv` |
| **= WORKING SET** | **208,248** | the analysed set |

Source-DB membership (a plasmid may belong to several): IMG-PR 136,318 · RefSeq 56,042 ·
GenBank 55,428 · PLSDB 47,215 · COMPASS 12,084 · mMGE 7,207 · DDBJ 4,331 · ENA 4,049 · others.
3 single-source plasmids are absent from the FASTA snapshot, so sequence-based tools (PlasAnn,
PlasmidFinder, mob_typer, RGI) run on **208,245**.

---

## 4. Pipeline stages and their outputs

Each stage writes a TSV that `scripts/build_metadata_master.py` folds into the master table.

**Environment & geography** — `scripts/recover_env_local.py` (IMG/PR GOLD),
`scripts/fetch_biosample_env.py` (NCBI BioSample), `scripts/fetch_sra_env.py` (mMGE SRA),
`scripts/scrape_img_geo.py` + `scripts/geo_utils.py` (geography), reconciled by
`scripts/reconcile_environment.py`. A raw label was recovered for **85.9%** of the set. The locked
habitat taxonomy (see `reports/environment_taxonomy_LOCKED.md`) gives:
Host-associated 65,108 · Environmental 24,693 · Engineered 14,371 · Unknown 39,331, plus two
**excluded** buckets — Simulated-artifact 64,658 (GOLD `Simulated communities`) and Lab-artifact 87.
→ **Analysis set = 143,503** for environmental/geographic views. `is_clinical` = 19,446 (9.3%).

**PLSDB curated enrichment** — `scripts/enrich_plsdb_metadata.py` → `plsdb_enrichment.tsv`, joined
for the **45,743** plasmids with a PLSDB record (Inc typing, relaxase/mobility, MOB cluster, geo,
ecosystem/disease tags, species, AMRFinderPlus genes). See `METADATA_ENRICHMENT_PLAN.md`.

**PlasAnn functional annotation** — sharded SLURM run: `scripts/plasann_shard.py` →
`scripts/plasann_array.sbatch`/`plasann_array2.sbatch` → `scripts/plasann_task.sh`, harvested by
`scripts/harvest_plasann.py` and `scripts/aggregate_plasann_genes.py`. **99.97%** annotated
(68 failures in `data/plasann_run/final_missing.txt`). Adds `plasann_*` columns (CDS count,
replicons, oriV/oriT, AMR/metal-biocide/conjugation/mobile-element/toxin-antitoxin/virulence counts).

**Authoritative replicon/Inc + MOB-suite typing** — see `reports/typing_methodology.md`.
PlasmidFinder (`scripts/plasmidfinder_run.sbatch`, ≥80% id / ≥60% cov, CGE thresholds): **41,555
plasmids (20%)** carry ≥1 replicon, **24,583 are multireplicon**. mob_typer
(`scripts/mobtyper_array.sbatch` → `scripts/mobtyper_task.sh`): **100%** get mobility + MOB cluster
(7,054 distinct clusters); host range for 58%. Harvested by `scripts/harvest_typing.py`
(+`scripts/aggregate_mob_full.py`). Any-standard-threshold replicon coverage: **95,442 (46%)**,
up from the 12% PLSDB-only baseline.

**CARD / RGI AMR annotation** — see `reports/card_amr_methodology.md`. CARD **v4.0.1**, RGI **6.0.8**
(`data/card/CARD_VERSION.txt`), run over the 595 shards via `scripts/card_array.sbatch` →
`scripts/card_task.sh`, harvested by `scripts/harvest_card.py` + `scripts/aggregate_card.py`.
Confident (Perfect+Strict) calls: **AMR+ = 29,396 (14.1%)**, multidrug 20,741 (9.96%). Validated
against PLSDB/AMRFinderPlus at **90.5% presence/absence concordance**. Adds `card_*` columns.

---

## 5. The master table (central output)

`data/plasmidscope_primary/plasmid_metadata_master.tsv` — **208,248 plasmids × 67 columns**, built
by `scripts/build_metadata_master.py`. `scripts/build_analysis_table.py` derives the analysis-ready
subset. Columns by group:

- **Identity/physical** — `plasmid_id`, `sources`, `topology`, `size_bp`, `gc_percent`, `n_source_dbs`
- **Mobility (PlasmidScope MOB-suite)** — `predicted_mobility`, `mob_families`
- **Environment/geo** — `env_status`, `hab_top`, `hab_sub`, `hab_channel`, `is_clinical`,
  `geo_country`, `geo_admin1`, `geo_lat`, `geo_lng`, `geo_precision`, `geo_source`, `collection_date`
- **PLSDB curated (`plsdb_*`)** — Inc/rep types, relaxase, mobility, MOB cluster, AMR genes, drug
  classes, ecosystem/disease tags, species, accession (columns 14–26)
- **PlasAnn (`plasann_*`)** — annotation flag, feature/CDS counts, replicons, oriV/oriT, and
  functional-category counts (columns 34–49)
- **Authoritative typing** — PlasmidFinder `pf_*` (columns 50–52) and mob_typer `mob_*`
  (columns 53–59)
- **CARD AMR (`card_*`)** — ARG counts, ARO list, drug classes, mechanisms, gene families,
  `card_multidrug` (columns 60–67)

---

## 6. Environments

Built with **micromamba** rooted project-locally at `envs/.micromamba_root/`
(`$HOME` quota is too small — see `tools/README.md`). Activate any env with:

```bash
export MAMBA_ROOT_PREFIX=/gorilla/proj/h-mel-phylo/h-mel-phylo/private/amar/plasmid_analyis/envs/.micromamba_root
ENV=$MAMBA_ROOT_PREFIX/envs/<env_name>
PYTHONNOUSERSITE=1 PATH="$ENV/bin:$PATH" $ENV/bin/<tool> ...
```

**Live environments** (present now; fresh specs in `envs/live_env_specs/`):

| env | purpose | key packages |
|---|---|---|
| `rgi_env` | CARD/RGI AMR annotation | rgi 6.0.8, diamond 2.2.4, blast 2.16 |
| `segmantx_env` | SegMantX duplication/segment analysis | pandas<2.2, blast, biopython, plotly |
| `stats_env` | notebook analysis & plotting | numpy, pandas, scipy, statsmodels, matplotlib-base (no seaborn) |

**Decommissioned phase environments** — used earlier, then removed to reclaim space once their stage
finished. Only their exported specs remain, documented in `tools/README.md`:
`envs/plasmid_phase0.environment.yml` (Phase-0 general stack), `envs/plasann_env.environment.yml`
(PlasAnn), `envs/mob_suite_env.environment.yml` (mob_typer — needs pinned old pandas).
`envs/segmantx_env.environment.yml` is the original SegMantX export (superseded by the live spec).

---

## 7. Repository map

```
plasmid_analyis/
├── PROJECT_OVERVIEW.md        ← this file
├── CLAUDE.md                  working conventions for this repo
├── RESEARCH_OUTLINE.md        research question & design
├── EXECUTION_PLAN.md          phased execution plan
├── METHODOLOGY.txt            master methods narrative (cited by the report files)
├── METADATA_ENRICHMENT_PLAN.md · PLASMIDSCOPE_PRIMARY_PLAN.md · CARD_AMR_PLAN.md
├── papers_of_interest.txt · environment.yml · .env_root
│
├── data/            (15G) all inputs, intermediates, and final tables
│   ├── plasmidscope_primary/   ← THE working set + every final table (master, typing, card, plasann…)
│   ├── PlasmidScope/ PLSDB/ GOLD/ Cus_PR/   source downloads & metadata
│   ├── plasann_run/  (9.4G) PlasAnn intermediates: shards/ + gbk/ + annot/  ← raw byproducts
│   ├── card_run/     (1.4G) RGI raw per-shard output + SLURM logs           ← raw byproducts
│   ├── typing/       PlasmidFinder + mob_typer intermediates + DB
│   └── processed/ refs/
│
├── envs/            (3.7G) micromamba root (.micromamba_root/) + env specs
├── scripts/         all pipeline scripts (see §4; grouped by stage)
├── notebooks/       5 investigation notebooks (built by scripts/build_*_notebook.py)
├── manuscript/      NAR draft + round-1 review + reassessment (thesis needs reframing)
├── reports/         6 methodology write-ups + figures/ (48 figures, 5 stages)
├── tools/           SegMantX (vendored) + tool/env provenance README
├── orit_db_folder/ transposon_db_folder/   auxiliary reference databases
```

`data/plasmidscope_primary/` is the one directory that matters for analysis: it holds the working
set and every final table. `data/plasann_run/{shards,gbk}` and `data/card_run/rgi` are **raw
per-shard byproducts** already rolled up into those final tables (see §9).

---

## 8. Reports, notebooks, figures

**Reports** (`reports/`, each names the exact scripts that produced it):
`plasmidscope_primary_characterization.md`, `plasmidscope_primary_working_set.md`,
`environment_taxonomy_LOCKED.md`, `typing_methodology.md`, `card_amr_methodology.md`,
`amr_onehealth_methodology.md` (Phase-3 AMR × One Health compartment × mobility).

**Notebooks** (`notebooks/`, regenerated by the matching `scripts/build_*_notebook.py`):
`plasmid_data_investigation.ipynb`, `plasann_annotation_investigation.ipynb`,
`typing_investigation.ipynb`, `card_amr_investigation.ipynb`,
`onehealth_amr_investigation.ipynb` (Phase-3 analysis + round-1 reassessment; run with the
`genesis_nb` kernel).

**Figures** (`reports/figures/`): investigation 17 · card 10 · typing 8 · plasann 7 · onehealth 7
(49 total).

---

## 9. Reproduce / regenerate

The master table is regenerated by re-running the folding step; individual layers regenerate from
their SLURM jobs (see each report's "Reproduce" section). Quick pointers:

```bash
# rebuild master + analysis table from the layer TSVs
python3 scripts/build_metadata_master.py
python3 scripts/build_analysis_table.py

# per-layer (each also documented in its report):
python3 scripts/reconcile_environment.py          # environment taxonomy
python3 scripts/harvest_typing.py                 # PlasmidFinder + mob_typer
python3 scripts/harvest_card.py && python3 scripts/aggregate_card.py   # CARD AMR
```

**Raw intermediates.** `data/plasann_run/{shards,gbk}` (9.2G) and `data/card_run/rgi` (1.4G) are the
raw per-shard tool outputs. They are already harvested into the final tables in
`data/plasmidscope_primary/`, so they are **not needed for any downstream analysis** — only for a
from-scratch re-annotation (an HPC array job). They are candidates for deletion when space is needed.

---

## 10. Housekeeping (2026-07-27)

- Reclaimed **1.7 GB** by clearing the micromamba package cache (`micromamba clean --all`); the three
  live environments are untouched.
- Removed `scripts/__pycache__` (regenerated automatically).
- Added `envs/live_env_specs/` with fresh reproducible specs for the live environments (the older
  `envs/*.environment.yml` files describe the now-decommissioned phase environments).
- No data files were deleted. The ~11 GB of raw intermediates in §9 are kept in place by choice;
  see §9 for what is safely reclaimable later if space is needed.
