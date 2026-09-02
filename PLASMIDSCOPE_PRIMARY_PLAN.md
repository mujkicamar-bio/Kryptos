# PlasmidScope-primary plan — project restart

*Authoritative plan as of 2026-07-09. Supersedes the IMG/PR-centric foundation in `EXECUTION_PLAN.md`
(Phases 0–1). The **biology questions** (`RESEARCH_OUTLINE.md`) and the **analysis methodology** for
Phases 2–5 of `EXECUTION_PLAN.md` (fusion/SegMantX, ARG routes, small-plasmid building blocks) are
unchanged and reused — only the **dataset and its foundation** change.*

## 1. Reframe

**PlasmidScope is now THE dataset.** Treat the project as if started fresh on PlasmidScope (Tang et al.,
NAR 2025); IMG/PR is demoted to *one source repository within* PlasmidScope, not a separate study. The
prior IMG/PR-only results (Phase 0/1/M2, Phase 4/M4) become a **reusable component** covering
PlasmidScope's IMG-PR slice, not the endpoint.

## 2. Locked decisions (from the 2026-07-09 discussion)

1. **Universe = PlasmidScope's own non-redundant set** (the `ALL` master, 852,600), not our
   per-source union. PlasmidScope deduplicated with **MMseqs2 v15.6f452 at 100% identity & 100%
   coverage** (`--cov-mode 0 -c 1.0 --min-seq-id 1.0`), collapsing cross-source *and* IMG/PR-internal
   duplicates; each surviving row is tagged with the **set** of source DBs (and their original IDs)
   it was consolidated from. So there are **no duplicate plasmids** anywhere in this set.
   *Confirmed (Phase A.1-A.3, `reports/plasmidscope_primary_characterization.md`): filtering `ALL`
   to `Completeness=="complete"` gives **208,360 complete plasmids** (348,617 raw source-entries
   collapsed). Composition: 65.7% metagenomic-only, 31.1% isolate-only (64,834), 3.2% both; IMG-PR
   membership 136,318 — below our old 154,680 because IMG/PR internal dups are now collapsed.*
2. **Two hard discards only** — everything else is kept:
   - **not complete/closed** (incomplete, fragments, contigs) — trust PlasmidScope `Completeness`.
   - **lab-made / synthetic** (cloning/expression/shuttle vectors, `synthetic construct`) —
     conservative detection, every hit logged with the signal that caught it.
3. **Environment is NOT a filter.** Plasmids with blank/unrecoverable environment are **kept** and
   flagged `environment: unknown` — decided on later, never dropped for lacking a habitat.
4. **Environment labels stay raw / as-provided.** Recover from the source DB and attach verbatim; **no
   reconciliation/harmonisation now** — we inspect them once collected and decide how to reconcile.
5. **Provenance is a set.** A sequence appearing in several sources keeps **all** source labels (and
   all environments those sources provide).
6. **Trust PlasmidScope's annotations** — do not re-run their tools (completeness, mobility, host,
   clustering are taken as given).
7. **Data-only for now — no typing tools yet** ("not yet").
8. **Meticulous accounting / full reproducibility** — every discarded plasmid is counted with a
   reason; every report cites the exact script(s) that produced it (standing project rule).

## 3. The inclusion funnel (every plasmid accounted for)

```
PlasmidScope non-redundant set (ALL, ~852,600)
  └─ keep COMPLETE / closed only        → discard incomplete + fragments/contigs   (counted, reason)
      └─ remove LAB-MADE / synthetic      → discard cloning/expression vectors etc. (counted, signal logged)
          = WORKING SET  — everything else kept
```

Only two gates. `Completeness == "complete"` is trusted from PlasmidScope; the topology breakdown
(circular / linear / DTR / concatemer) is **reported** for transparency but not used as an extra gate.

**Lab-made detection (conservative, logged).** Flag as lab-made only on unambiguous signals:
host/organism ∈ {`synthetic construct` (NCBI taxon 32630), `Cloning vector`, `Expression vector`}, or
the record title explicitly names a cloning/expression/shuttle vector. Each flagged plasmid is logged
with the matching signal so the call is auditable and the threshold can be dialed later.

## 4. Per-plasmid tagging (two independent axes)

Every kept plasmid carries:
- **Provenance** — `{ all source DBs it appears in }` + isolate-vs-metagenome, reconstructed by
  cross-referencing the per-source tables against the non-redundant master.
- **Environment (raw)** — the verbatim label(s) recovered from the source:
  - IMG/PR-sourced → GOLD `ecosystem` string from `IMGPR_plasmid_data.tsv` (already in hand).
  - Isolate-sourced (RefSeq/GenBank/PLSDB/COMPASS/ENA/DDBJ/Kraken2) → NCBI BioSample
    `isolation_source` / `host` / `geo_loc_name` (Entrez, machinery already built).
  - mMGE → its metagenome study context.
  - none recoverable → `environment: unknown` (kept).

## 5. Phase A — data characterization (immediate next step; NO tools)

Produces the funnel numbers and the working set. All Python/stdlib + the existing download/fetch
scripts; zero SLURM, zero typing tools.

| # | Task | Output |
|---|---|---|
| A.1 | Acquire the PlasmidScope non-redundant master metadata (`ALL.plasmid_list`); reconcile its complete count against our per-source downloads | `data/plasmidscope_primary/all_metadata.tsv` + reconciliation note |
| A.2 | Apply gate 1 (complete/closed); log discards with reasons + topology breakdown | ledger rows |
| A.3 | Reconstruct multi-source provenance per non-redundant plasmid | provenance column (set) |
| A.4 | Flag & remove lab-made/synthetic (conservative, per §3); log each with its signal | ledger rows + `discarded_labmade.tsv` |
| A.5 | Recover **raw** environment per kept plasmid from source DBs (IMG/PR GOLD + BioSample + mMGE); blanks → `unknown` | `working_set_environment_raw.tsv` |
| A.6 | Assemble the working set; write the characterization report: the discard ledger + kept-set breakdown by provenance and by raw environment (incl. an explicit `unknown` bucket) | `data/plasmidscope_primary/working_set.tsv` + `reports/plasmidscope_primary_characterization.md` |

**Gate to later phases:** review the characterization report (working-set size, habitat composition,
unknown fraction) before any tool is run.

### Phase A RESULTS (executed 2026-07-09 → `reports/plasmidscope_primary_working_set.md`)

Funnel: `ALL` 852,600 → complete 208,360 → **− 112 lab-made** (cloning/expression/shuttle vectors +
synthetic constructs; `discarded_labmade.tsv`, each logged with its signal) → **WORKING SET =
208,248**. Deliverables: `data/plasmidscope_primary/{complete_provenance,env_local,working_set}.tsv`
+ `discarded_labmade.tsv`; scripts `characterize_plasmidscope_all.py`, `recover_env_local.py`,
`fetch_biosample_env.py`, `fetch_sra_env.py`, `assemble_working_set.py`.

- **Raw environment recovered for 178,972 / 208,248 (85.9%)**; the rest kept as `unknown`.
  Channels: IMG/PR GOLD 111,708 · NCBI BioSample (all INSDC sources) 63,762 · mMGE SRA 3,517.
- Lifestyle: metagenomic 65.7% · isolate 31.1% · mixed 3.2%.
- **mMGE caveat:** only 3,657/7,207 mMGE plasmids encode a real SRA run (env recovered for 96% of
  those); the other 3,550 carry assembly-contig IDs (`mMGEs_k119_…`) with no run — env `unknown`
  from the id; recovering them would need mMGE's own contig→sample table (optional follow-up).
- Optional enhancement not yet done: PLSDB's curated metadata layer for PLSDB-sourced plasmids.
- **Environment reconciliation (raw → habitat taxonomy) is deferred**, per plan.

## 6. Later phases (reuse existing methodology, on the new working set)

- **Replicon typing** — the one thing PlasmidScope does *not* provide (no Inc/rep type/count in its
  download). Required before any multireplicon/architecture claim. Approach & compute to be decided at
  the gate (likely MOB-suite + Rep-HMM as before, now on the working set).
- **Phase 2** (fusion rules & cost of fusion, SegMantX), **Phase 3** (ARG spread routes), **Phase 5**
  (small plasmids as building blocks), **Phase 4** (environmental distribution) — methodology per
  `EXECUTION_PLAN.md`, re-run on the working set, always **stratified by provenance/habitat** so the
  metagenomic-vs-isolate composition shift is visible, never pooled silently.

## 7. Deferred decisions (explicitly parked)

- **Environment reconciliation** — how to harmonise raw labels into a clean habitat taxonomy
  (and whether "clinical" is a top-level stratum or nested under Host-associated;Human).
- **Replicon typing approach** and whether to reuse the IMG-PR slice's existing calls.
- **How the `unknown`-environment plasmids are used** in each analysis.

## 8. Reproducibility

Raw under `data/plasmidscope_primary/`; scripts in `scripts/`; every report names its scripts. The
funnel is a deterministic, re-runnable pipeline: from the non-redundant master to the working set,
each plasmid has a recorded disposition (kept, or discarded + reason + signal).
