# Environment taxonomy — LOCKED (2026-07-10)

**Pipeline:** `scripts/reconcile_environment.py` → `data/plasmidscope_primary/environment_reconciled.tsv`,
folded into `plasmid_metadata_master.tsv` by `scripts/build_metadata_master.py`. Regenerate:
```
python3 scripts/reconcile_environment.py && python3 scripts/build_metadata_master.py
```

## What each plasmid is labelled with

| axis | column(s) | meaning |
|---|---|---|
| **Database provenance** | `sources` | the set of source databases the plasmid comes from (IMG-PR, RefSeq, GenBank, PLSDB, COMPASS, ENA, DDBJ, mMGE, Kraken2, TPA) |
| **Environment (reconciled)** | `hab_top`, `hab_sub` | the habitat the plasmid was sampled from — *is it wastewater, human faeces, soil, blood…* |
| **Clinical tag** | `is_clinical` | additional flag; environment-based (clinical body site OR disease OR hospital), nested under Host-associated;Human |
| **Environment channel** | `hab_channel` | which source supplied the env label (GOLD / BioSample / SRA) |

**Dropped on purpose:** the `lifestyle` / `sample_nature` (metagenomic-vs-isolate) axis. The source-DB tag
does **not** reflect the true sample origin — e.g. clinical datasets contain metagenomic samples — so it
was confusing and is not used. Provenance is captured by `sources`; origin is captured by the environment.

## The locked taxonomy

**Top level (3 GOLD-aligned domains) + 2 non-habitat buckets:**

| hab_top | n | note |
|---|---:|---|
| Host-associated | 65,108 | human (with body site), livestock, poultry, insect, plant, … |
| Environmental | 24,693 | aquatic, terrestrial/soil, air, extreme |
| Engineered | 14,371 | wastewater/sewage, bioreactor, built-env, food, solid-waste, industrial — **real field/engineering sites only** |
| **Simulated-artifact** | 64,658 | GOLD `Modeled;Simulated communities` — **excluded from analysis** (not a real habitat) |
| **Lab-artifact** | 87 | GOLD `Lab enrichment` (60, lab cultures on defined media) + content-free bare `Engineered` (27) — **excluded from analysis** (no real field site) |
| **Unknown** | 39,331 | unresolved + unlabelled + geography-only, merged |

**Analysis set = 143,503** (208,248 − 64,658 simulated − 87 lab-artifact). Excluded buckets keep their
real sequences (counted in provenance/physical/typing) but are dropped from every environmental and
geographic view.

**Sub-habitats (`hab_sub`):** Human: blood / urine / respiratory / skin-wound / oral / gut-faeces /
other-clinical / unspecified · Livestock: pig / cattle / other · Poultry · Rodent · Companion · Fish ·
Insect · Plant · Algae · Fungi · Other invertebrate · Gut/faeces (unspecified host) · Aquatic: marine /
freshwater · Terrestrial/soil · Air · Extreme · Wastewater/sewage · Bioreactor · Built environment ·
Food/fermentation · Solid waste/compost · Industrial/remediation.

**`is_clinical` = 19,446 (9.3%)** — clinical body site (blood/urine/respiratory/skin-wound/other-clinical),
or a PLSDB disease tag, or a hospital/patient/nosocomial context.

## Locked resolutions to the four prior open questions

1. **Clinical** → stays `Host-associated;Human;<body site>` with an additional `is_clinical` tag (not a 4th
   top level). Environment-based, so it is unaffected by the source-DB confusion.
2. **Metagenomic-vs-isolate** → **dropped** (see above). Not a stratifier.
3. **Ambiguous bare terms** → `feces`/`stool` → `Human: gut/faeces`; bare `gut` → `Gut/faeces: unspecified host`.
4. **Geography-only** → folded into `Unknown` habitat.

Rules are transparent keyword + GOLD mappings (first-match-wins) in `scripts/reconcile_environment.py`;
adjust a rule and re-run to change a mapping.
