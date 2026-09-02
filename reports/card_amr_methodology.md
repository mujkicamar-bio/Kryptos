# CARD / RGI antibiotic-resistance annotation — methodology & results

**Scope.** Authoritative, ARO-normalised antibiotic-resistance annotation of the full 208,248-plasmid
PlasmidScope working set with **RGI `main`** against **CARD**, replacing the two weak AMR signals the
master already held (PLSDB's AMRFinderPlus calls, present only for the ~22% PLSDB subset; PlasAnn's
loose `plasann_n_amr` count, no drug class / mechanism / confidence tier). Companion to
`CARD_AMR_PLAN.md`, `reports/typing_methodology.md`, and the top-level `METHODOLOGY.txt`.

Pipeline scripts (in order): `scripts/plasann_shard.py` (shared shards), `scripts/card_array.sbatch`
→ `scripts/card_task.sh`, `scripts/harvest_card.py`, `scripts/aggregate_card.py`; folded into the
master by `scripts/build_metadata_master.py`. Version record: `data/card/CARD_VERSION.txt`.

---

## 1 · Databases & tools (pinned)

- **CARD** `card.json` `_version` **4.0.1** (release 2025-05-29), downloaded 2026-07-23 from
  `https://card.mcmaster.ca/latest/data`; loaded with `rgi load --card_json card.json --local` →
  `data/card/localDB/`. Canonical **protein homolog + SNP/variant** models (WILDCARD/prevalence not
  loaded — not needed for per-plasmid gene detection).
- **RGI 6.0.8**, **DIAMOND 2.2.4**, **BLAST 2.16.0+**, in the isolated micromamba env
  `envs/.micromamba_root/envs/rgi_env`, always run with `PYTHONNOUSERSITE=1` and the env `bin/`
  prepended to `PATH` (project convention).

## 2 · Method

Input `data/plasmidscope_primary/working_set.fna.gz` (208,245 sequences; 3 single-source plasmids
absent from the FASTA snapshot, `working_set.fna.missing.txt`) was split into 595 gzipped multi-FASTA
shards (~350 plasmids each). Each shard was annotated with:

```
rgi main --input_type contig -a DIAMOND -n 2 --clean --local -d plasmid --include_loose
```

RGI calls ORFs (Prodigal) and aligns predicted proteins to CARD with DIAMOND, reporting one row per
ORF hit with a **`Cut_Off`** of Perfect / Strict / Loose. The RGI `Contig` column equals the
working-set `plasmid_id` (direct join key). `--include_loose` retains every cut-off in the raw output
for audit.

**Run infrastructure** (`card_task.sh`, mirroring the proven `mobtyper_task.sh`): the DIAMOND index
was pre-built once into `localDB/`, then copied to node-local `$SNIC_TMP` once per node (`flock`) so
tasks only read it; RGI output and all intermediates were written to node-local disk and only the
final `.txt` copied back — the quota'd project FS (500k-inode limit) never holds RGI's transient
files. SLURM array `0-594%300` on `-A uppmax2025-2-42 -p pelle`, `-c 2 --mem=8G -t 1:00:00`; ~87 s
per shard. **All 595 shards completed rc=0; all 208,245 plasmids processed (no drops).**

**Confidence filtering.** The confident set folded into the master is **Perfect + Strict only**;
Loose is retained in `card_hits.tsv` but never counted (CARD's standard convention — Loose hits fall
below the curated bitscore and are dominated by false positives).

## 3 · Outputs

- `data/plasmidscope_primary/card_hits.tsv` — one row per ORF hit, **all cut-offs**
  (634,948 rows: Perfect 77,250 + Strict 48,101 + Loose 509,597). Columns include `plasmid_id`,
  `orf_id`, `cut_off`, `best_hit_aro`, `aro`, `model_type`, `drug_class`, `resistance_mechanism`,
  `amr_gene_family`, `best_identities`, `best_hit_bitscore`, `pct_length_ref`, `nudged`.
- `data/plasmidscope_primary/card_per_plasmid.tsv` — one row per **carrier** plasmid (Perfect+Strict).
- **Master** `plasmid_metadata_master.tsv` (now **67 columns**), new fields:
  `card_n_arg`, `card_n_arg_unique`, `card_aro_list`, `card_drug_classes`, `card_n_drug_classes`,
  `card_resistance_mechanisms`, `card_amr_gene_families`, `card_multidrug`. Non-carriers = 0 / empty.

## 4 · Headline results

- **AMR+ (≥1 Perfect/Strict ARG): 29,396 / 208,248 = 14.1%.**
- **Multidrug (≥2 distinct drug classes): 20,741** (9.96% of all plasmids; 70.6% of carriers).
- Median ARGs among carriers = 2; maximum on a single plasmid = 74.
- **Top ARGs:** `sul1` (5,832), `TEM-1` (5,501), `qacEdelta1` (5,465), `sul2` (4,622), `APH(6)-Id`
  (4,418), `tet(A)` (4,119), `APH(3'')-Ib` (4,091), `mphA`, `AAC(6')-Ib10`, `aadA2/aadA`, `QnrS1`,
  `dfrA14` — the canonical class-1-integron / mobile clinical resistome.
- **Top drug classes:** penicillin β-lactam (14,613), aminoglycoside (14,020), cephalosporin
  (13,415), sulfonamide (9,478), tetracycline (8,671), disinfectants/antiseptics (7,919), monobactam,
  diaminopyrimidine, carbapenem (6,637), macrolide, fluoroquinolone, phenicol.
- **Mechanisms:** antibiotic inactivation (22,577) > efflux (12,803) > target replacement (10,409) >
  target alteration (6,596) > target protection (5,773).

## 5 · Validation — PLSDB vs CARD concordance

Independent check on the 45,743 plasmids with a PLSDB record (PLSDB AMR = AMRFinderPlus, a different
tool and database from CARD/RGI):

| | count |
|---|---|
| both AMR+ | 12,951 |
| both AMR− | 28,450 |
| CARD+ / PLSDB− | 22 |
| PLSDB+ / CARD-strict− | 4,320 |
| **presence/absence agreement** | **90.5%** |

The disagreement is almost entirely one-directional and explained by the **conservative threshold**,
not by CARD missing genes: **93% (4,038/4,320)** of the PLSDB+/CARD-strict− plasmids *do* carry a
CARD **Loose** hit (deliberately excluded); only **282 (0.6% of the overlap)** have no CARD hit at
any cut-off, and CARD calls AMR where PLSDB does not in just 22 cases. Counting Loose would raise
agreement to ~99.4%. This validates both the pipeline and the Perfect+Strict choice.

## 6 · Caveats

- **Circular wrap:** ORF calling linearises each sequence, so an ARG spanning the arbitrary
  breakpoint may be split — minor for short ARGs.
- **Model type:** SNP/mutational-resistance hits are flagged via `model_type` in `card_hits.tsv` so
  mutational calls (e.g. gyrA) are not conflated with acquired ARGs; the plasmid signal is
  overwhelmingly acquired homolog hits.
- **Determinism:** no RNG; tool + database versions pinned; the loaded `localDB` is captured on disk.

## 7 · Not covered here (follow-on = Phase 3 analysis)

AMR prevalence by habitat / mobility / Inc-type, ARG co-occurrence, and conjugative-vs-non-mobilisable
AMR carriage — a separate analysis/notebook step now that the layer exists.
