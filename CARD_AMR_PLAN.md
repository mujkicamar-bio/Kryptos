# CARD / RGI antibiotic-resistance annotation plan

**Goal.** Add an authoritative, ARO-normalised, whole-set antibiotic-resistance layer to the 208,248
working-set plasmids by running **RGI (Resistance Gene Identifier) `main`** against the **CARD**
database. This supersedes the two weak AMR signals already in the master — PLSDB's calls (present
only for the ~22% PLSDB subset) and PlasAnn's loose `plasann_n_amr` count (no drug class, no
resistance mechanism, no confidence tier). Feeds the project's ARG-spread question (EXECUTION_PLAN
Phase 2/3).

## Tool & database versions (pinned)
Recorded in `data/card/CARD_VERSION.txt`:
- CARD data (`card.json` `_version`): **4.0.1**, release 2025-05-29, downloaded 2026-07-23 from
  `https://card.mcmaster.ca/latest/data`.
- RGI **6.0.8**, DIAMOND **2.2.4**, BLAST **2.16.0+**.
- Env: isolated micromamba env `envs/.micromamba_root/envs/rgi_env` (bioconda `rgi`), run with
  `PYTHONNOUSERSITE=1` and env `bin/` prepended to `PATH` (project convention — see `tools/README.md`).
- Models: canonical CARD **protein homolog + SNP/variant** models (WILDCARD/prevalence not loaded —
  not needed for per-plasmid gene detection).

## Input
`data/plasmidscope_primary/working_set.fna.gz` (208,245 sequences; 3 single-source plasmids absent
from the PlasmidScope FASTA snapshot, listed in `working_set.fna.missing.txt`). FASTA header = the
working-set `plasmid_id` (e.g. `COMPASS_AY362554.1`); RGI's `Contig` output column carries it
verbatim → direct join key to the master.

## Pipeline (scripts, in order)
1. **Env + DB (login node, one-time).** Build `rgi_env`; download+extract CARD to `data/card/`;
   `rgi load --card_json card.json --local` → `data/card/localDB/`. A throwaway `rgi main` run
   pre-builds the DIAMOND index (`protein.db.dmnd`) inside `localDB/`, making it read-only at scale.
2. **Shard.** `scripts/plasann_shard.py 595` → `data/plasann_run/shards/shard_XXXX.fna.gz`
   (~350 plasmids/shard, round-robin so the rare megaplasmids spread evenly). Shards are reused
   across PlasAnn / mob_typer / CARD.
3. **RGI array.** `scripts/card_array.sbatch` → `scripts/card_task.sh` per shard, on
   `-A uppmax2025-2-42 -p pelle`, `--array=0-594%300`, `-c 2 --mem=8G -t 1:00:00`. Each task:
   copies `localDB` to node-local `$SNIC_TMP` once per node (flock — same fix that removed the
   mob_typer network-FS contention), runs
   `rgi main --input_type contig -a DIAMOND -n 2 --clean --local -d plasmid --include_loose`
   with output + all intermediates on node-local disk, and copies back only the final
   `data/card_run/rgi/shard_XXXX.txt` (the redundant `.json` is discarded). `--include_loose` keeps
   every cut-off in the raw output for audit. Steady-state ≈ 87 s/shard.
4. **Harvest.** `scripts/harvest_card.py` → `data/plasmidscope_primary/card_hits.tsv`: one row per
   ORF-level ARG hit (all cut-offs), analysis-relevant columns only, `Contig`→`plasmid_id`.
5. **Aggregate.** `scripts/aggregate_card.py` → `data/plasmidscope_primary/card_per_plasmid.tsv`:
   one row per plasmid, **Perfect+Strict only (Loose excluded)**. Multi-value CARD fields
   (Drug Class / Resistance Mechanism / AMR Gene Family) are `;`-split, stripped, unioned.
6. **Fold into master.** `scripts/build_metadata_master.py` (extended) joins `card_per_plasmid.tsv`
   on `plasmid_id`, filling non-carriers with 0/empty. New columns:
   `card_n_arg`, `card_n_arg_unique`, `card_aro_list`, `card_drug_classes`, `card_n_drug_classes`,
   `card_resistance_mechanisms`, `card_amr_gene_families`, `card_multidrug`.
7. **Report.** `reports/card_amr_methodology.md` — versions, thresholds, filtering, script index, and
   a PLSDB-vs-CARD concordance check on the overlap as a validation gate.

## Decisions
- **Confidence:** Perfect+Strict is the confident set folded into the master; Loose is retained in
  `card_hits.tsv` only (auditable, never counted).
- **Aligner:** DIAMOND (fast; adequate for CARD's small protein set).
- **Input type:** `contig` (RGI calls ORFs, then aligns) — correct for assembled complete plasmids;
  `rgi bwt` (reads) is not applicable.

## Caveats (baked in)
- **Circular wrap:** ORF calling linearises the sequence, so an ARG spanning the arbitrary
  breakpoint could be split — minor for short ARGs.
- **Model type:** SNP/mutational-resistance hits are flagged via `model_type` in `card_hits.tsv` so
  mutational calls (e.g. gyrA) are not conflated with acquired ARGs; the plasmid signal is
  overwhelmingly acquired homolog hits.
- **Determinism:** no RNG; versions pinned; the loaded `localDB` is captured on disk.

## Not in scope here (follow-on = Phase 3 analysis)
AMR prevalence, top ARGs, drug-class spectrum, and AMR × habitat / mobility / Inc-type — a separate
analysis/notebook step once the layer lands.
