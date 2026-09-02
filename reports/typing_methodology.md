# Replicon (Inc) typing & MOB-suite typing — methodology

**Scope.** Authoritative, standard-threshold replicon/Inc typing and full MOB-suite typing of the
208,248-plasmid working set, to answer the project's central *multireplicon* question and to add plasmid
mobility, MOB cluster (taxonomy) and host range. Companion to the top-level `METHODOLOGY.txt` and
`reports/environment_taxonomy_LOCKED.md`.

Pipeline scripts: `scripts/plasmidfinder_run.sbatch`, `scripts/mobtyper_task.sh`,
`scripts/mobtyper_array.sbatch`, `scripts/harvest_typing.py`; folded into the master by
`scripts/build_metadata_master.py`. Output table: `data/plasmidscope_primary/typing_authoritative.tsv`.

---

## 0 · Why this run (what we already had)

Before this run the master already held **mobility for all 208,248 plasmids** (`predicted_mobility`,
from PlasmidScope's own MOB-suite) and PlasmidFinder Inc types for only the **24,872 PLSDB-covered**
plasmids (12%). PlasAnn's replicon calls (§notebook 4) exist but use a permissive **60%** BLAST and are
raw feature rows — an annotation signal, not a normalized Inc call. So the gap was **standard-threshold
Inc typing across the whole set** (needed for a reliable multireplicon count) plus the **MOB cluster /
host range**, which nothing else provided. The user chose to run **both PlasmidFinder and mob_typer**.

---

## 1 · PlasmidFinder (replicon / Inc typing)

**Method — identical criteria to CGE PlasmidFinder / PLSDB:** `blastn` of every working-set plasmid
against the PlasmidFinder replicon database, keeping a hit at **≥80% nucleotide identity** AND **≥60%
coverage of the replicon reference** (alignment length ÷ subject length).

- **Database.** The PlasmidFinder DB bundled with PlasAnn (`~/.plasann/Database/plasmidfinder.fasta`)
  was confirmed to be the genuine CGE database (482 replicon sequences, standard `IncFII(S)`,
  `IncI1-I(Alpha)`, `IncN`, `Col…`, `rep…` nomenclature). Its `>>` double-caret headers were cleaned to
  `>` and a nucleotide BLAST DB built once → `data/typing/plasmidfinder.*`.
- **Run.** One multithreaded SLURM job on `pelle` (`plasmidfinder_run.sbatch`): decompress the working
  set to node-local `$SNIC_TMP`, `blastn -task blastn -perc_identity 80 -evalue 1e-5 -max_target_seqs
  200 -num_threads 6 -mt_mode 1`, output `qseqid sseqid pident length slen` gzipped to
  `data/typing/plasmidfinder_blast.tsv.gz`. *(20 threads triggered a BLAST thread-creation error; 6
  threads + `-mt_mode 1` runs cleanly. Two independent runs agreed to <0.1%, confirming completeness.)*
- **Call rule at parse time** (`harvest_typing.py`): per plasmid, collect distinct replicon names
  passing ≥80/≥60; the Inc *family* is the name with allelic variants collapsed (e.g.
  `IncFII(pCoo)`→`IncFII`).

**Result:** **41,555 plasmids (20%)** carry ≥1 PlasmidFinder replicon; **24,583 are multi-replicon
(≥2 distinct)** = 59% of typed. Leading families: Col / ColRNAI / Col440I, then the IncF complex
(IncFII, IncFIB, IncFIA, IncFIC), IncN, IncR. *(20% overall is expected — the metagenomic majority
carries no PlasmidFinder-typeable replicon, and the DB is Enterobacterial / Gram-positive focused.)*

---

## 2 · MOB-suite `mob_typer` (replicon + relaxase + MPF + oriT + mobility + cluster + host range)

**Tool & environment.** `mob_typer` v3.1.9 in a purpose-built micromamba env `envs/mobtyper_env`
(mob_suite + **pandas 1.5.3** — mob_suite this version is incompatible with pandas ≥2; the env was built
with its package cache redirected to scratch and `--always-copy` so it is self-contained: its own BLAST,
mash and the 2.9 GB MOB database, 88 files).

**Run.** SLURM array over the same 595 gzipped shards used for PlasAnn (`mobtyper_array.sbatch` →
`mobtyper_task.sh`): each task decompresses its shard to `$SNIC_TMP`, runs
`mob_typer --multi --num_threads 1`, writes one TSV per shard to `data/typing/mob/shard_*.txt`.
`sample_id` is the raw FASTA header, so it joins to the master `plasmid_id` directly (no sanitization).

**The key performance fix.** mob_typer was initially **~81 s/plasmid** and threw
`sqlite3.OperationalError: disk I/O error`. Root cause: the host-range step reads `taxa.sqlite` (ete3
NCBI taxonomy) and the MOB blast/mash databases **from the shared network filesystem**, which cannot
serve 300 concurrent SQLite readers and is pathologically slow for SQLite. **Fix:** each node copies the
whole 2.9 GB MOB DB to node-local `$SNIC_TMP` once (flock-guarded, reused by that node's tasks) and
mob_typer reads the local copy → **~2 s/plasmid (≈40× faster)** and no I/O errors. The full run then
finished in ~30 min (`-t 2h`, `%300`). *(Two earlier attempts at 350 plasmids/shard timed out at 6 h and
12 h before this root cause was found — wasted cluster time, since corrected.)*

**Result (all 208,245 sequenced plasmids typed, 100%):**

| field | coverage | notes |
|---|---:|---|
| `mob_mobility` | 208,245 (100%) | non-mobilizable 59% / mobilizable 28% / conjugative 13% (matches PlasmidScope's independent MOB-suite — cross-validated) |
| `mob_cluster` (primary cluster id) | 208,245 (100%) | 7,054 distinct MOB clusters (plasmid taxonomy) |
| `mob_rep_types` | 93,925 (45%) | MOB-suite's own replicon calls |
| `mob_relaxase` | — | MOBP / MOBF / MOBQ / MOBV lead |
| `mob_mpf` | — | mating-pair-formation MPF_T / MPF_F / MPF_I |
| `mob_host_range` | 120,210 (58%) | Enterobacterales, Bacteroides, Pseudomonadota, Acinetobacter, … |

---

## 3 · Combined coverage & the master

`harvest_typing.py` → `typing_authoritative.tsv`, folded into `plasmid_metadata_master.tsv`
(now **59 columns**). New columns:

`pf_inc_types`, `pf_inc_families`, `pf_n_inc` · `mob_rep_types`, `mob_relaxase`, `mob_mpf`, `mob_orit`,
`mob_mobility`, `mob_cluster`, `mob_host_range`.

**Any standard-threshold replicon (PlasmidFinder ∪ MOB-suite ∪ PLSDB): 95,442 (46%)** — up from the
12% PLSDB-only baseline. **Multireplicon plasmids: 24,583** (PlasmidFinder ≥2 distinct).

---

## 4 · Caveats

- **20% PlasmidFinder / 45% MOB replicon coverage is biology, not failure**: metagenomic plasmids
  (the majority) frequently lack any typeable replicon, and both databases are skewed to
  Enterobacterial / Gram-positive isolates.
- PlasmidFinder replicon **families** collapse allelic variants; the full allele strings are kept in
  `pf_inc_types`.
- `mob_host_range` is a coarse taxonomic-convergence prediction, not an experimentally observed range.
- PlasAnn's loose 60% replicon calls remain in `plasann_replicons` as an annotation signal; the
  `pf_*` / `mob_*` columns are the authoritative typing.

---

## 5 · Reproduce

```
# PlasmidFinder DB (once)
sed 's/^>>/>/' ~/.plasann/Database/plasmidfinder.fasta > data/typing/plasmidfinder.fasta
makeblastdb -in data/typing/plasmidfinder.fasta -dbtype nucl -out data/typing/plasmidfinder
sbatch scripts/plasmidfinder_run.sbatch

# mob_typer (reuses the 595 PlasAnn shards; needs envs/mobtyper_env)
sbatch scripts/mobtyper_array.sbatch

# harvest + fold into the master
python scripts/harvest_typing.py
python scripts/build_metadata_master.py
```
