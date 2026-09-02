# Master table — data dictionary

Column reference for **`data/plasmidscope_primary/plasmid_metadata_master.tsv`**: one row per
working-set plasmid, **208,248 rows × 67 columns**. Assembled by
[`scripts/build_metadata_master.py`](../scripts/build_metadata_master.py), which folds together the
per-layer tables (each layer's producing script is named in its section below). Coverage % is the
share of the 208,248 rows that are non-null.

**How to read coverage:** count/flag columns (`*_n_*`, `*_has_*`, `plasann_annotated`, `is_clinical`,
`card_multidrug`) are ~100% populated because non-carriers are recorded as `0`; the list/string
columns are populated **only where the feature exists**, so their coverage ≈ the prevalence of that
feature. The `plsdb_*` block is gated to the ~22% of plasmids with a PLSDB record.

---

## A · Identity, provenance & physical properties (cols 1–8)
Source: PlasmidScope `all_metadata` + [`assemble_working_set.py`](../scripts/assemble_working_set.py).

| # | column | type | cov. | meaning |
|--:|---|---|--:|---|
| 1 | `plasmid_id` | str | 100% | Unique id = the **representative** (first source) id of the deduplicated plasmid. Primary key / join key. |
| 2 | `sources` | str | 100% | `;`-joined set of source databases the plasmid was consolidated from (e.g. `GenBank;PLSDB;RefSeq`). 102 distinct combinations. |
| 3 | `topology` | str | 100% | Reported topology: `circular` / `linear` / … (4 values). |
| 4 | `size_bp` | int | 100% | Plasmid length in base pairs. |
| 5 | `gc_percent` | float | 100% | GC content (%). |
| 6 | `predicted_mobility` | str | 100% | **PlasmidScope's own** MOB-suite mobility: `conjugative` / `mobilizable` / `non-mobilizable`. (Mobility signal 1 of 3 — see §H.) |
| 7 | `mob_families` | str | 37.5% | PlasmidScope MOB relaxase families (e.g. `MOBP`); populated only for mobilizable/conjugative plasmids. |
| 8 | `n_source_dbs` | int | 100% | Number of distinct source DBs the plasmid came from (1–7). |

---

## B · Environment & habitat taxonomy (cols 9–13)
Source: [`reconcile_environment.py`](../scripts/reconcile_environment.py); rules locked in
[`environment_taxonomy_LOCKED.md`](environment_taxonomy_LOCKED.md).

| # | column | type | cov. | meaning |
|--:|---|---|--:|---|
| 9 | `env_status` | str | 100% | Was any label recovered? `labelled` (178,965) / `unknown` (29,283). |
| 10 | `hab_top` | str | 100% | Top-level domain: `Host-associated` / `Environmental` / `Engineered` (real) + `Simulated-artifact` / `Lab-artifact` (excluded) + `Unknown`. **Coarse — do not stratify One Health by this** (it lumps human, animal & plant together). |
| 11 | `hab_sub` | str | 100% | Fine-grained habitat (42 values, e.g. `Human: blood`, `Livestock: pig`, `Terrestrial/soil`, `Wastewater/sewage`). **Build habitat compartments from this** + `is_clinical`. |
| 12 | `hab_channel` | str | 85.6% | Which source supplied the environment label: `GOLD` / `BioSample` / `SRA` (provenance of the label, not the habitat). |
| 13 | `is_clinical` | int | 100% | Orthogonal 0/1 flag — clinical body site **or** PLSDB disease tag **or** hospital/patient/nosocomial context. `1` = 19,446 (9.3%). Splits `Human` into clinical vs community. |

> **Exclusions:** drop `hab_top ∈ {Simulated-artifact (64,658), Lab-artifact (87)}` from every
> environment/geography analysis — they are not real habitats. Analysis set = **143,503**.

---

## C · PLSDB curated enrichment (cols 14–26)
Source: [`enrich_plsdb_metadata.py`](../scripts/enrich_plsdb_metadata.py). All gated by `plsdb_acc` —
present only for the **45,743 (22%)** plasmids with a PLSDB record. AMR here is **AMRFinderPlus**
(a different tool/DB from CARD; note it also reports some virulence genes).

| # | column | type | cov. | meaning |
|--:|---|---|--:|---|
| 14 | `plsdb_acc` | str | 22.0% | PLSDB accession — the gate for this whole block. |
| 15 | `plsdb_inc_types` | str | 11.9% | PLSDB PlasmidFinder Inc types (allele-level). |
| 16 | `plsdb_n_inc` | float | 22.0% | Count of Inc types (`0` for a PLSDB plasmid with none; `NaN` if non-PLSDB). |
| 17 | `plsdb_rep_types` | str | 16.3% | PLSDB rep types. |
| 18 | `plsdb_relaxase` | str | 11.9% | PLSDB relaxase (MOB) type. |
| 19 | `plsdb_mobility` | str | 22.0% | PLSDB mobility call (mobility signal 2 of 3). |
| 20 | `plsdb_mob_cluster` | str | 22.0% | PLSDB MOB cluster id. |
| 21 | `plsdb_amr_genes` | str | 8.3% | `|`-joined AMRFinderPlus gene calls (may include virulence genes). |
| 22 | `plsdb_n_amr` | float | 22.0% | Count of AMRFinderPlus genes. |
| 23 | `plsdb_drug_classes` | str | 7.3% | Drug classes (AMRFinderPlus). |
| 24 | `plsdb_ecosystem_tags` | str | 18.9% | PLSDB ecosystem tags (e.g. `host_associated`). |
| 25 | `plsdb_disease_tags` | str | 4.9% | PLSDB disease tags (feeds `is_clinical`). |
| 26 | `plsdb_species` | str | 22.0% | PLSDB host species — the only curated host taxonomy in the master. |

---

## D · Geography & collection date (cols 27–33)
Source: [`build_metadata_master.py`](../scripts/build_metadata_master.py) +
[`geo_utils.py`](../scripts/geo_utils.py) + [`scrape_img_geo.py`](../scripts/scrape_img_geo.py).
Country is parsed from a **place-name** field (not reverse-geocoded); coordinates follow a
point → country-centroid → country-only ladder.

| # | column | type | cov. | meaning |
|--:|---|---|--:|---|
| 27 | `geo_country` | str | 83.8% | Country label (normalized). **Inflated by artifacts — see caveat.** |
| 28 | `geo_admin1` | str | 71.3% | Sub-national region / free-text locality. |
| 29 | `geo_lat` | float | 83.6% | Latitude (decimal). |
| 30 | `geo_lng` | float | 83.6% | Longitude (decimal). |
| 31 | `geo_precision` | str | 83.8% | `point` (real coord) / `country_centroid` / `country_only` / **`artifact`** (simulated, fake coordinate). |
| 32 | `geo_source` | str | 83.8% | Which route produced it: `img_gold`, `plsdb`, `biosample_latlon`, `geoname_centroid`, `geoname`, `img_gold_simulated`. |
| 33 | `collection_date` | str | 74.4% | Reported sampling date (string). |

> **Caveat:** the 64,543 `Simulated-artifact` plasmids carry the JGI/Berkeley compute-site coordinate
> and are labelled `geo_country = "USA"`, `geo_precision = "artifact"`. Raw `USA` = 101,806 but **real
> USA ≈ 37,263**. Filter `geo_precision != "artifact"` before any geographic analysis (leaves 109,795
> real country labels).

---

## E · PlasAnn functional annotation (cols 34–49)
Source: [`harvest_plasann.py`](../scripts/harvest_plasann.py) +
[`aggregate_plasann_genes.py`](../scripts/aggregate_plasann_genes.py). Count columns are `0` for
plasmids lacking the feature (hence ~100% populated).

| # | column | type | cov. | meaning |
|--:|---|---|--:|---|
| 34 | `plasann_annotated` | int | 100% | 1 if PlasAnn annotation succeeded (99.97% of the set). |
| 35 | `plasann_n_features` | float | 100% | Total annotated features. |
| 36 | `plasann_n_cds` | float | 100% | Coding sequences (CDS). |
| 37 | `plasann_replicons` | str | 30.7% | PlasAnn replicon calls — **loose 60% BLAST, an annotation signal, NOT authoritative typing** (use `pf_*`). |
| 38 | `plasann_n_replicons` | float | 100% | Count of PlasAnn replicons (0–7). |
| 39 | `plasann_has_oriv` | float | 100% | 1 if an origin of replication (oriV) was annotated. |
| 40 | `plasann_has_orit` | float | 100% | 1 if an origin of transfer (oriT) was annotated. |
| 41 | `plasann_n_amr` | float | 100% | PlasAnn AMR feature count — **loose signal; use `card_*` for AMR.** |
| 42 | `plasann_n_metal_biocide` | float | 100% | Metal/biocide-resistance feature count. |
| 43 | `plasann_n_conjugation` | float | 100% | Conjugation-machinery feature count. |
| 44 | `plasann_n_mob_dna` | float | 100% | Mobile-DNA feature count. |
| 45 | `plasann_n_mobile_element` | float | 100% | Mobile-element feature count. |
| 46 | `plasann_n_toxin_antitoxin` | float | 100% | Toxin–antitoxin system feature count. |
| 47 | `plasann_n_virulence` | float | 100% | Virulence feature count. |
| 48 | `plasann_n_ncrna` | float | 100% | Non-coding RNA feature count. |
| 49 | `plasann_n_maintenance` | float | 100% | Plasmid-maintenance feature count. |

---

## F · Replicon / Inc typing — PlasmidFinder (cols 50–52)
Source: [`harvest_typing.py`](../scripts/harvest_typing.py) (PlasmidFinder ≥80% id / ≥60% cov).
**This is the authoritative Inc typing.**

| # | column | type | cov. | meaning |
|--:|---|---|--:|---|
| 50 | `pf_inc_types` | str | 20.0% | PlasmidFinder Inc types, allele-level. |
| 51 | `pf_inc_families` | str | 20.0% | Inc families (allelic variants collapsed). |
| 52 | `pf_n_inc` | float | 100% | Count of distinct Inc families (`0` = untyped/none). **`pf_n_inc ≥ 2` = multireplicon** — the paper's definition. |

---

## G · Mobility typing — MOB-suite `mob_typer` (cols 53–59)
Source: [`harvest_typing.py`](../scripts/harvest_typing.py). **`mob_mobility` and `mob_cluster` are the
authoritative mobility call and plasmid-taxonomy id used in the analysis.**

| # | column | type | cov. | meaning |
|--:|---|---|--:|---|
| 53 | `mob_rep_types` | str | 45.1% | MOB-suite replicon calls. |
| 54 | `mob_relaxase` | str | 37.5% | Relaxase (MOB) type. |
| 55 | `mob_mpf` | str | 14.4% | Mating-pair-formation (MPF) type — conjugation apparatus. |
| 56 | `mob_orit` | str | 12.5% | oriT type. |
| 57 | `mob_mobility` | str | 100% | **Authoritative** mobility: `conjugative` / `mobilizable` / `non-mobilizable` (mobility signal 3 of 3). |
| 58 | `mob_cluster` | str | 100% | MOB cluster id — plasmid taxonomy (7,054 clusters). Used for clonal-redundancy dereplication. |
| 59 | `mob_host_range` | str | 57.7% | Predicted host range (coarse taxonomic-convergence prediction, not observed). |

---

## H · Antibiotic resistance — CARD / RGI (cols 60–67)
Source: [`harvest_card.py`](../scripts/harvest_card.py) +
[`aggregate_card.py`](../scripts/aggregate_card.py). **Confident set = Perfect + Strict** (Loose
excluded). CARD v4.0.1 / RGI 6.0.8. **This is the authoritative AMR layer.**

| # | column | type | cov. | meaning |
|--:|---|---|--:|---|
| 60 | `card_n_arg` | int | 100% | Number of confident ARG hits. **`card_n_arg > 0` = AMR-positive** (14.1% of the set). |
| 61 | `card_n_arg_unique` | int | 100% | Number of distinct ARGs. |
| 62 | `card_aro_list` | str | 14.1% | `;`-joined ARO gene names (carriers only). |
| 63 | `card_drug_classes` | str | 14.1% | `;`-joined drug classes conferred. |
| 64 | `card_n_drug_classes` | int | 100% | Count of distinct drug classes. |
| 65 | `card_resistance_mechanisms` | str | 14.1% | Resistance mechanisms (e.g. antibiotic inactivation, efflux). |
| 66 | `card_amr_gene_families` | str | 14.1% | AMR gene families. |
| 67 | `card_multidrug` | int | 100% | 1 if ≥2 distinct drug classes (multidrug). |

---

## Cross-cutting: which signal is authoritative?
Several biological properties appear more than once, from different tools. Use the **authoritative** one:

| property | signals present | **use** |
|---|---|---|
| **Mobility** | `predicted_mobility` (PlasmidScope), `plsdb_mobility` (PLSDB, 22%), `mob_mobility` (mob_typer) | **`mob_mobility`** (100%, our run) |
| **Replicon / Inc** | `plsdb_inc_types` (22%), `plasann_replicons` (loose 60%), `mob_rep_types`, `pf_inc_types`/`pf_inc_families` | **`pf_*`** (PlasmidFinder ≥80/≥60); `pf_n_inc ≥ 2` = multireplicon |
| **AMR** | `plsdb_amr_genes`/`plsdb_n_amr` (AMRFinderPlus, 22%), `plasann_n_amr` (loose) | **`card_*`** (RGI/CARD Perfect+Strict) |
| **Host taxonomy** | `plsdb_species` (22%), `mob_host_range` (predicted) | `plsdb_species` for observed; join PlasmidScope `Host` for wider (sparse) coverage |

*Generated from `plasmid_metadata_master.tsv`; column meanings trace to
[`build_metadata_master.py`](../scripts/build_metadata_master.py) and the per-layer scripts cited above.*
