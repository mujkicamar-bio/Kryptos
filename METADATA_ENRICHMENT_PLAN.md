# Metadata enrichment plan — maximal metadata for every working-set plasmid

*Goal (2026-07-09): for each of the **208,248** working-set plasmids, recover **every metadata
dimension we can** — geography, ecosystem/habitat, host taxonomy, disease, replicon/Inc & MOB typing,
AMR/virulence, collection date, functional annotation — each tagged with the source it came from and
kept raw. Supersedes nothing; extends Phase A. Every field's provenance is recorded (no silent
merging). Reconciliation of overlapping labels stays deferred (see notebook §6).*

## 0. Metadata dimensions we want (the target schema)

| dimension | field(s) | why it matters |
|---|---|---|
| Provenance | source DB set, lifestyle | stratifier (done, Phase A) |
| Sequence type | size, topology, GC, completeness | intrinsic (done) |
| **Replicon typing** | Inc/rep type(s) + **count** | **the multireplicon question** — currently missing for most |
| Mobility | MOB relaxase, mpf, oriT, predicted mobility, cluster | conjugation/spread |
| **AMR / virulence** | ARG genes + drug class, VFGs | ARG-spread analysis (Phase 3) |
| **Geography** | country/name + **lat/lng** | spatial distribution |
| **Ecosystem / habitat** | hierarchical ecosystem tags | environmental distribution (Phase 4) |
| Host | host organism + taxonomy lineage | host range |
| Disease | disease/symptom ontology | clinical stratum |
| Collection date | year | temporal trends |
| Annotation | protein-coding genes, IS, integron, oriT | architecture |

## 1. Source-by-source strategy (what each layer fills, for which subset)

Nearly everything hangs off the **NCBI/INSDC accession backbone**; one identifier unlocks several
layers. Coverage numbers are of the 208,248 working set.

| Layer | Subset it covers | Fields it adds | Status |
|---|---|---|---|
| **PlasmidScope base** | all 208,248 | size, topology, GC, host, MOB_type, mobility, cluster | ✅ done (`analysis_table.tsv`) |
| **IMG/PR GOLD** | 136,316 IMG-PR | ecosystem hierarchy (curated) | ✅ done (`env_local.tsv`) |
| **NCBI BioSample** | ~64k INSDC isolates | isolation_source, host, geo_loc_name | ✅ done (`env_raw.tsv`) |
| **PLSDB curated 2024_05_31_v2** | **45,743** PLSDB-sourced | **Inc typing, MOB, AMR, lat/lng, ecosystem tags, disease, taxonomy** | ✅ **done now** (`plsdb_enrichment.tsv`) |
| **mMGE SRA** | 3,657 mMGE w/ run | metagenome sample context | ✅ done (`sra_env_raw.tsv`) |
| NCBI source-qualifier top-up | non-PLSDB INSDC gap | `collection_date`, `lat_lon`, `country` from GenBank source feature | ⬜ TODO (cheap efetch) |
| **Replicon + MOB typing (compute)** | the ~162k **not** in PLSDB | Inc/rep type & **count**, relaxase, mobility | ⬜ **TODO — needs tools**, the real gap |
| AMRFinderPlus (compute) | all / any subset | ARG + VFG uniformly | ⬜ optional (Phase 3) |
| PIPdb (pathogen DB) | pathogen subset only | VFG, HMRG, integron, IS, risk score | ⬜ optional complement — see §4 |

### What PLSDB just added (executed — `scripts/enrich_plsdb_metadata.py`)
Of our 45,743 PLSDB-matched plasmids (96.9% of the PLSDB slice): MOB typing **100%**, geo-coordinates
**89%**, ecosystem tags **86%**, PlasmidFinder Inc typing **54%**, AMR genes **38%**, disease **22%**.
PLSDB's `ECOSYSTEM_tags` are a **clean controlled vocabulary** (`host_associated`,
`circulatory_system/blood`, `gastrointestinal_system/fecal`, `urinary_system`, `respiratory_system`,
`soil`, `wastewater`, `food`, `anthropogenic/hospital`, `disease`) — a ready-made reconciliation
target for notebook §6.

## 2. The one dimension that cannot be fetched — replicon typing

PLSDB gives Inc/MOB typing only for its own 45,743 (and PlasmidFinder fires on ~54% of those). The
remaining **~162,000 plasmids have no replicon type or count** — and *multireplicon architecture is the
project's central question*. This is **compute, not a download**: run MOB-suite (`rep_type(s)` +
count) and/or PlasmidFinder + Rep-HMM on the working-set FASTA, reusing PLSDB's calls where they exist
so we don't recompute 45k. This is the gate to Phases 2–5 and the natural next step after this plan.

## 3. Deliverable — the master metadata table ✅ DONE

**`data/plasmidscope_primary/plasmid_metadata_master.tsv`** (`scripts/build_metadata_master.py`):
208,248 plasmids × 34 columns, merging base (`analysis_table.tsv`) + reconciliation
(`environment_reconciled_PROPOSED.tsv`) + PLSDB (`plsdb_enrichment.tsv`) + a coalesced geography block.
Every field keeps its source (column prefix / `geo_source`); blanks stay blank. **PIPdb is excluded
from the master and from all analysis (project decision) — see §4.** Replicon typing for the untyped
~162k is still to be folded in once computed.

### Geography — maximised (§3a)
Geography is coalesced per plasmid, priority: explicit **point** coords (BioSample `lat_lon` → else
PLSDB lat/lng) → **country centroid** from a place name → country-name-only. Centroids are **empirical**
(median of the 34,703 real coordinate points we already hold, giving 129 countries) + a 44-country
hand fallback for the tail (`scripts/geo_utils.py`); `lat_lon` parses 100%.

**Result: 63,441 plasmids (30.5%) now carry a country; all of them get coordinates** (46,516 exact
points + 16,926 country-centroids; 57 name-only). By lifestyle: **isolate 87.9%**, mixed 59.3%,
metagenomic **1.9%** — i.e. geography is essentially complete for the INSDC/isolate portion (63,441 =
88.9% of the 71,361 INSDC-reachable plasmids) and near-absent for the metagenomic majority. Top
countries: China 14.1k, USA 11.4k, South Korea 3.8k, Japan 3.6k, UK 3.1k.

### IMG/M sample geography — ✅ DONE (2026-07-09, via user's JGI session token)
The metagenomic majority (IMG/PR) had **no geography in any local source**. GOLD's bulk export stayed
403, but the user's **JGI session token** unlocked IMG's `MetaDetail` pages, which expose
Latitude/Longitude/Geographic Location/Country/Collection Date per sample. The 136k IMG/PR plasmids
collapse to **26,149 samples**, of which **10,582 are metagenomes** (isolate taxa carry no coordinate on
IMG and are already INSDC-geolocated, so skipped).

`scripts/scrape_img_geo.py` (12 concurrent workers — IMG pages render ~6 s each; session as
`jgi_session` cookie via `IMG_SESSION` env var; resumable) scraped all 10,582 → `data/GOLD/img_taxon_geo.tsv`:
**10,466/10,582 (98.9%) have geography, 94.3% exact coordinates.** `build_metadata_master.py` fans these
onto plasmids via `taxon_oid` (priority: BioSample lat_lon → PLSDB → **IMG/M** → country centroid).

**⚠ Artifact caveat:** the Simulated-communities samples are geolocated to JGI's compute site
(Berkeley → USA); those 64,658 plasmids are flagged `geo_precision=artifact` / `geo_source=img_gold_simulated`
and **excluded from real geography**.

**Result: real geography jumped 44% → 76.5%** — of the 143,590 non-artifact plasmids, **109,786 (76.5%)
now have a country and 109,442 (76.2%) have coordinates** (real metagenomic 1.9% → **68%**; isolate
87.9%; the artifact-inclusive metagenomic figure is 83% but those coords are JGI's compute site). Top
real countries: USA 37k, China 16k, UK 5.3k, S.Korea 3.9k, Japan 3.8k. Remaining ~34k with no geography
= isolates/samples lacking any source coordinate. **Notebook §7 visualises this** (world map, coverage,
top countries, collection year — `scripts/build_investigation_notebook.py`).

**Order of execution:** (1) ✅ master built from what's in hand incl. maximised INSDC geography;
(2) ⬜ optional GOLD/IMG `taxon_oid` geography for metagenomic plasmids; (3) ⬜ replicon/MOB typing on
the untyped ~162k [compute, the gate]; (4) ⬜ fold typing into the master.

## 4. PIPdb (NAR 2025, 53(D1):D169) — assessed; **complement, not a replacement**

**What it is.** *Plasmids in Pathogens Database* — 792,964 **plasmid segment clusters (PSCs)** predicted
(Plasmer/Platon/geNomad) from 1,009,571 **assembled genomes** of **450 pathogenic species**, with rich
annotation (ARG, VFG, heavy-metal resistance, integron, IS, oriT, relaxase, T4CP) and a risk score, at
`nmdc.cn/pipdb`.

**Would it be better than PlasmidScope? No — for three structural reasons.**
1. **Pathogen-only.** 450 pathogenic species reintroduces exactly the clinical/host-associated bias this
   project fights; it discards the environmental & metagenomic plasmidverse (66% of our set) that is the
   whole point of the environmental-distribution aim.
2. **Not complete/closed plasmids.** PIPdb's units are *predicted plasmid contigs/segment clusters* from
   genome assemblies — fragments, not closed replicons. Our locked inclusion rule is complete/closed
   only; PIPdb's content type fails that gate.
3. **Not accession-joinable.** It is organised by sequence clusters, not INSDC accessions, so enriching
   our plasmids would need sequence-level matching (ANI/MMseqs), and much of its INSDC-derived content
   already sits inside PlasmidScope's sources — high integration cost, largely redundant sequence-wise.

**Where it *is* useful (optional).** As a **complementary annotation reference for the pathogen subset**
in Phase 3 (ARG/VFG spread): its VFG / heavy-metal / integron / risk-score layers and its host-range
metadata could cross-validate and enrich our clinical plasmids. Verdict: **keep PlasmidScope primary;
consider PIPdb only as a targeted annotation cross-reference later, not a data-source swap.** (Download
mechanism at `nmdc.cn/pipdb` is a JS interface; a bulk/accession export would need confirming before any
integration.)

## 5. Reproducibility

Raw PLSDB tables under `data/PLSDB/plsdb2025_meta/` (from the official META archive, v2024_05_31_v2;
3.65 GB archive discarded after extracting the 199 MB of curated tables). Scripts:
`enrich_plsdb_metadata.py` (done), `build_metadata_master.py` (next), typing pipeline (gate). Every
report/table names its script; every field records its source.
