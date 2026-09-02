# The small cryptic plasmidome — methodology and results

*Written 2026-08-26. Every number in this report is produced by a named script; nothing here is
asserted without a pipeline behind it.*

## Provenance

| what | produced by |
|---|---|
| Notebook (all analysis) | `scripts/build_small_cryptic_notebook.py` → `notebooks/small_cryptic_investigation.ipynb` (execute with the **genesis_nb** kernel) |
| Result tables | `data/plasmidscope_primary/smallcryptic_*.tsv` (13 tables, written by the notebook) |
| Figures | `reports/figures/small_cryptic/fig1..fig5*.png` |

Inputs, all pre-existing:

| input | produced by | documented in |
|---|---|---|
| `data/plasmidscope_primary/plasmid_metadata_master.tsv` (67 cols × 208,248) | `scripts/build_metadata_master.py` | `PROJECT_OVERVIEW.md` §5 |
| `data/plasmidscope_primary/mob_full.tsv.gz` | `scripts/aggregate_mob_full.py` | `reports/typing_methodology.md` |
| `data/plasmidscope_primary/working_set.tsv` | `scripts/assemble_working_set.py` | `reports/plasmidscope_primary_working_set.md` |
| `data/plasann_run/annot/{shard,mshard}_*.tsv.gz` (10,094,230 CDS rows) | `scripts/plasann_harvest_shard.py` | `PROJECT_OVERVIEW.md` §4 |
| CARD/RGI `card_*` columns (Perfect+Strict) | `scripts/harvest_card.py`, `scripts/aggregate_card.py` | `reports/card_amr_methodology.md` |
| PlasmidFinder `pf_*` / mob_typer `mob_*` columns | `scripts/harvest_typing.py` | `reports/typing_methodology.md` |
| habitat taxonomy (`hab_top`, `hab_sub`, `is_clinical`) | `scripts/reconcile_environment.py` | `reports/environment_taxonomy_LOCKED.md` |

**Analysis set.** The locked exclusion of `hab_top ∈ {Simulated-artifact, Lab-artifact}` is applied
once, giving **143,503** plasmids (asserted in the notebook). Ecological strata use `hab_sub`.

> **Pitfall recorded for reuse:** the PlasAnn per-gene harvest wrote **two** shard series in
> `data/plasann_run/annot/` — `shard_*.tsv.gz` (595, primary run) and `mshard_*.tsv.gz` (477,
> reshard of initially-missing plasmids). Globbing only one pattern silently analyses ~27% of the
> data. The notebook globs `*shard_*.tsv.gz` and reports the resulting coverage (99.9%).

## Definitions

- **small** — `size_bp < 10,000`. A Gaussian KDE on log₁₀(size) puts the two modes at **4,748 bp**
  and **91,983 bp** and the antimode at **18,757 bp**. The 10 kb cut is *not* the antimode: it is a
  deliberately conservative line inside the small mode's shoulder, chosen so the group contains
  unambiguously small-mode material rather than the sparse valley between modes. Because that is a
  judgement call, every headline quantity is re-reported at 5/8/10/15/20 kb — a range bracketing the
  antimode — in `smallcryptic_threshold_sensitivity.tsv`; none of the conclusions move.
- **cryptic** — no AMR (`card_n_arg`, Perfect+Strict), no virulence, no metal/biocide, and no
  conjugation machinery (PlasAnn counts). Toxin–antitoxin and mobilization genes are *not*
  disqualifying; their prevalence is a result, not a filter.
- **groups** — `small-cryptic`, `small-cargo`, `large-conjugative` (≥10 kb and mob_typer
  `conjugative`), `large-other`.

## Results

### 1 — Scale (`smallcryptic_funnel.tsv`, `smallcryptic_group_profile.tsv`)

Plasmid size is bimodal with modes at 4,748 bp and 91,983 bp (antimode 18,757 bp). Small plasmids are **70,243 (48.9%)** of the
analysis set; **47,031 (32.8% of all complete plasmids)** are small *and* payload-free. Median small
cryptic plasmid: 4,148 bp, 43.6% GC, **5 CDS**.

### 2 — The typing blind spot (`smallcryptic_typing_blindspot.tsv`, `smallcryptic_typing_by_size.tsv`)

Three independent replicon-typing methods, % of each group receiving a call:

| group | n | PlasmidFinder | mob_typer | PlasAnn | any of the three | host-range call |
|---|---:|---:|---:|---:|---:|---:|
| small-cryptic | 47,031 | 10.5 | 27.5 | 17.6 | **28.9** | 35.3 |
| small-cargo | 23,212 | 34.6 | 67.5 | 49.0 | 70.2 | 84.9 |
| large-conjugative | 25,631 | 45.5 | 91.3 | 77.2 | 92.2 | 99.4 |
| large-other | 47,629 | 17.0 | 49.9 | 30.0 | 50.8 | 65.1 |

**71.1% of small cryptic plasmids (33,427 molecules, 23.3% of the entire analysis set) receive no
replicon call from any method.** Stable at 68.4–72.7% across all five size thresholds.

### 3 — Dispersal at matched sampling depth (`smallcryptic_dispersal_*.tsv`)

MOB-suite primary clusters are lineage proxies; breadth is measured by rarefaction (each cluster
subsampled to exactly *K* members, 200 replicates) because raw breadth scales with sampling effort.

Pooled (K = 10): cryptic lineages span **4.58** habitats vs **5.0** for conjugative (ratio 0.92,
p = 4.9 × 10⁻⁴) and **4.50** countries vs **5.40** (ratio 0.83, p = 8.2 × 10⁻⁹). Effect sizes are small
(rank-biserial 0.13–0.25) and the ordering is stable across K ∈ {5, 10, 20, 40}.

**Restricted to isolate-derived plasmids** (metagenome-only records excluded), the difference is not
significant: habitats 0.94×, p = 0.06; countries 1.05×, p = 0.11. The defensible statement is
**near-parity**: elements encoding no transfer machinery achieve essentially the same habitat and
geographic spread per lineage as conjugative plasmids.

### 4 — Host-range breadth is inverted (`smallcryptic_hostrange_rank.tsv`)

Among plasmids receiving a mob_typer host-range call:

| group | genus | family | order | class | phylum | multi-phyla |
|---|---:|---:|---:|---:|---:|---:|
| small-cryptic | 48.7 | 6.9 | 11.0 | 4.5 | 2.9 | **25.9** |
| small-cargo | 27.2 | 22.2 | 10.3 | 8.3 | 2.6 | **29.3** |
| large-conjugative | 15.9 | 10.4 | **51.0** | 8.0 | 11.9 | 2.8 |

Small plasmids are **bimodal** — narrow (genus) or cross-phylum — while conjugative plasmids sit at a
single intermediate rank (order, largely *Enterobacterales*). Multi-phylum calls are **9.3× more
common** in small cryptic than in conjugative plasmids.

*Artefact control* (`smallcryptic_hostrange_vs_mashdist.tsv`): the call is neighbourhood-based, so a
distant/noisy neighbourhood could inflate it. It does not — in the *identical-neighbour* (mash
distance 0) bin the multi-phylum rate is 25.5% for small-cryptic vs 2.5% for large-conjugative, and
for small-cryptic the rate *falls* with increasing distance (17.4% in the >0.1 bin).

### 5 — Coding dark matter (`smallcryptic_dark_matter.tsv`, `smallcryptic_category_mix.tsv`)

Over all 10,094,230 PlasAnn CDS (99.9% working-set coverage):

| group | total CDS | pooled % unannotatable | median per-plasmid dark fraction | % of plasmids with **every** CDS unannotatable |
|---|---:|---:|---:|---:|
| small-cryptic | 234,112 | **89.9** | 1.00 | **73.1** |
| small-cargo | 130,740 | 51.5 | 0.50 | 1.2 |
| large-conjugative | 3,753,993 | 43.0 | 0.33 | 0.1 |
| large-other | 5,000,639 | 73.3 | 0.80 | 14.2 |

### 6 — Module architecture (`smallcryptic_module_architecture.tsv`)

Presence of a detectable replication (`rep`), mobilization (`mob`) or toxin–antitoxin (`TA`) module:

| architecture | n | % |
|---|---:|---:|
| **none detected** | 29,333 | **62.4** |
| rep | 6,774 | 14.4 |
| rep + mob | 6,084 | 12.9 |
| mob | 2,490 | 5.3 |
| TA / mob+TA / rep+TA / rep+mob+TA | 2,350 | 5.0 |

**62.4% of small cryptic plasmids are complete, closed replicons with no recognisable replication,
mobilization or maintenance module of any kind.**

### 7 — Bias audit (`smallcryptic_provenance.tsv`, `smallcryptic_within_provenance.tsv`)

Small cryptic plasmids are **73% IMG-PR-only** (metagenome-derived) vs 11.6% for large-conjugative, so
provenance is a genuine confounder and is handled by stratification rather than adjustment. The blind
spot survives inside both strata:

| stratum | group | n | % untyped by all three |
|---|---|---:|---:|
| IMG-PR only (metagenomic) | small-cryptic | 34,330 | 82.7 |
| IMG-PR only (metagenomic) | large-conjugative | 2,961 | 15.4 |
| isolate-derived only | small-cryptic | 8,749 | **40.8** |
| isolate-derived only | large-conjugative | 20,256 | **6.0** |

A 6.8× gap remains in isolate-derived data alone. The effect is attenuated by provenance but not
explained by it.

## Limitations

1. `hab_sub` and `geo_country` are present for a minority of small cryptic plasmids; the dispersal
   analysis is restricted to clusters with ≥ K annotated members and is not representative of the
   whole compartment.
2. mob_typer host-range prediction is neighbourhood-based and therefore partly circular with cluster
   assignment. §6d bounds but does not eliminate this.
3. "Non-mobilizable" is a *detection* statement about the mob_typer relaxase/oriT databases, not a
   demonstration of immobility — this is the leading alternative explanation for §3 and is the first
   thing the follow-up work should test (relaxed oriT search against `orit_db_folder/`).
4. Metagenome-derived complete plasmids may include mis-assembled or non-plasmid circular elements;
   the isolate-only stratification in §7 is the control for this, not a fix.

## Reproduce

```bash
cd /gorilla/proj/h-mel-phylo/h-mel-phylo/private/amar/plasmid_analyis
python3 scripts/build_small_cryptic_notebook.py
jupyter nbconvert --to notebook --execute --inplace \
    --ExecutePreprocessor.timeout=3000 notebooks/small_cryptic_investigation.ipynb
```
