# Is the dark plasmidome a PlasAnn artifact? — Pfam-A validation

*Written 2026-08-28. Every number is produced by a named script; nothing is asserted without a
pipeline behind it. Two defects were found and fixed during this run — both are documented in
§6 rather than quietly corrected, because one of them silently corrupted an earlier cross-tab.*

## 1. The question

`reports/dark_orf_clustering.md` showed that the payload-free small plasmidome's unannotatable
proteome clusters into ~92,750 recurrent families, and that 86% of dark ORFs have no homolog among
proteins PlasAnn named. That comparison is **PlasAnn against itself**. It cannot distinguish:

- **artifact** — the proteins are well known, PlasAnn simply failed to name them; versus
- **real** — the proteins are genuinely uncharacterized.

This report settles it against an external, curated database.

## 2. Provenance

| what | produced by |
|---|---|
| Query set (3 sets, chunked) | `scripts/build_pfam_query.py` → `data/pfam_run/chunks/`, `query_manifest.txt` |
| The search | `scripts/pfam_array.sbatch` → `scripts/pfam_task.sh` → `data/pfam_run/domtbl/` |
| Hit-rate tables | `scripts/summarize_pfam.py` → `pfam_per_seq.tsv`, `darkfam_pfam.tsv`, `pfam_summary.txt` |
| The verdict | `scripts/pfam_verdict.py` → `dark_novel_families.tsv`, `pfam_verdict.txt` |
| Notebook | `notebooks_my/dark_plasmidome.ipynb` §4 (recomputes everything; nothing pasted) |

Inputs: `data/dark_orf_run/` (see `reports/dark_orf_clustering.md`) and the Pfam-A database below.

**Database.** Pfam-A **38.2** (2026-01, based on UniProtKB 2025_03), **30,134 families**, downloaded
from `https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/` into `data/refs/pfam/` and
`hmmpress`ed. **Software:** HMMER **3.4** (Aug 2023), `panaroo` conda env.

**Compute.** 88-task SLURM array on `pelle` (`-A uppmax2025-2-42`, `-c 4`, `--mem=12G`,
`--array=0-87%44`). Wall time ~6 minutes; all 264 task records `COMPLETED`.

## 3. Method

### 3.1 Three query sets, not one

A bare "X% of dark families hit Pfam" is uninterpretable. The run therefore searches three sets in
one pass, tagged in the FASTA header as `{SET}::{id}`:

| set | n | purpose |
|---|---:|---|
| `DARKREP` | 92,752 | the test set — one representative per dark family |
| `NAMED` | 78,397 | **positive control** — every protein PlasAnn *did* name |
| `VALID` | 3,678 | **propagation check** — all members of 200 sampled families (10–50 members, seed 20260828) |

`NAMED` is load-bearing. If Pfam also missed most PlasAnn-named plasmid proteins, a low `DARKREP`
rate would say nothing about novelty — only that Pfam covers plasmid proteins badly. `VALID` tests
the other assumption: that a representative's Pfam status may stand for its whole family.

### 3.2 Threshold

```bash
hmmsearch --cut_ga --noali --cpu 4 --domtblout <out> data/refs/pfam/Pfam-A.hmm <chunk>
```

`--cut_ga` uses Pfam's **curated per-family gathering thresholds** — the criterion Pfam itself uses
for family membership. No e-value was chosen by us, so the result cannot be tuned to taste.
`hmmsearch` (profiles vs sequence DB) is used rather than `hmmscan`; identical results, far faster.
Note that in `domtblout` the **target is our sequence** and the **query is the Pfam profile**.

### 3.3 Reproducing

```bash
cd /gorilla/proj/h-mel-phylo/h-mel-phylo/private/amar/plasmid_analyis
PY=/gorilla/home/amujkic/.conda/envs/genesis/bin/python
$PY scripts/build_pfam_query.py          # -> 88 chunks
sbatch scripts/pfam_array.sbatch         # ~6 min wall
$PY scripts/summarize_pfam.py
$PY scripts/pfam_verdict.py
```

## 4. Results

### 4.1 The control works — the dark set is genuinely depleted

| query set | queries | with a Pfam hit | % |
|---|---:|---:|---:|
| `NAMED` (positive control) | 78,397 | 60,124 | **76.7** |
| `DARKREP` | 92,752 | 24,357 | **26.3** |

Pfam recognises three-quarters of the plasmid proteins PlasAnn named, so the assay is sound. Dark
representatives hit at **0.34×** that rate. Weighted by family size, **149,160 of 378,552 dark ORFs
(39.4%)** sit in a family whose representative matches Pfam; **60.6% do not**. Among the 5,623
families with ≥10 members, **40.5%** hit.

### 4.2 The verdict: measurably both

Every dark ORF, cross-tabbed by what can name it (% of the 378,552):

| | no Pfam match | Pfam match |
|---|---:|---:|
| **dark-only** (no PlasAnn-named homolog) | **57.4** | **28.5** |
| has a PlasAnn-named homolog | 3.2 | 10.9 |

- **28.5% — PlasAnn annotation failure.** No named homolog in PlasAnn's own output, yet Pfam names
  it outright.
- **57.4% — the genuine residue.** Neither PlasAnn nor Pfam can name it.

### 4.3 The artifact half is worse than a coverage gap

Classifying what Pfam rescued by whether the domain falls inside a category PlasAnn *already has*
(pattern list in `scripts/pfam_verdict.py`, reproduced in the notebook):

| scope | % of rescued dark ORFs |
|---|---:|
| outside PlasAnn's categories | 52.7 |
| replication | 17.8 |
| mobilization / conjugation | 14.2 |
| toxin–antitoxin | 8.9 |
| transposition / recombination | 5.4 |
| partition / maintenance | 1.0 |

**47.3% of rescued dark ORFs — 18.6% of the entire dark proteome — carry a domain inside a category
PlasAnn already annotates.** The top rescued families are `Rep_1`, `Mob_Pre`, `Rep3_N`, `NikA-like`,
`HTH_3`, `RepL`, `Relaxase`, `MobA_MobL`, `Phage_integrase`, `Resolvase`, `ParE_toxin`, `RelB`,
`YoeB_toxin`, `TraD`, `T4SS-DNA_transf`.

These are replication, mobilization and toxin–antitoxin proteins that PlasAnn has slots for and
still labelled "Open reading frame". This also re-explains an earlier observation: the dark proteome
*looked* like replicate-and-mobilize machinery in `reports/dark_orf_clustering.md` because a
sizeable slice of it literally **is** that machinery, unlabelled.

### 4.4 What survives

**3,343 substantial families (≥10 members) have no Pfam domain at all**, holding **115,811 dark ORFs
(30.6% of the dark proteome)**. Of those families, **72.9% span ≥2 MOB lineages**, median **3
lineages / 4 habitats**. Widest-ranging:

| representative | proteins | plasmids | lineages | habitats |
|---|---:|---:|---:|---:|
| `IMGPR_plasmid_3300022589_000018\|2` | 2,487 | 2,468 | 88 | 17 |
| `IMGPR_plasmid_3300042343_000004\|14` | 747 | 743 | 84 | 18 |
| `GenBank_CP027661.1\|11` | 637 | 621 | 76 | 28 |
| `IMGPR_plasmid_2654587591_000001\|2` | 286 | 280 | 71 | 18 |
| `GenBank_CP096868.1\|2` | 246 | 236 | 65 | 14 |

Full list: `data/pfam_run/dark_novel_families.tsv`.

### 4.5 Propagation is good, not perfect

Across 200 families and 3,678 members, a member agrees with its representative's Pfam status
**96.1%** of the time; **23 of 200** families are internally mixed (5–95% of members hitting).
Family-level rates therefore carry a few percent of slop.

## 5. Interpretation

The dark plasmidome is **not** one thing. Roughly:

- **~29%** is a PlasAnn annotation failure that Pfam corrects, and most of that is core plasmid
  machinery PlasAnn was built to recognise.
- **~57%** resists both a plasmid-specific annotator and a curated 30,134-family domain database.

The novelty claim survives, but narrower than before, and now on external evidence rather than one
tool's silence: *a large, lineage-spanning protein space in the small plasmidome is invisible to
both plasmid annotation and Pfam.* The claim in `reports/dark_orf_clustering.md` §R3 that 86% of
dark ORFs lack a named homolog is arithmetically unchanged but **must no longer be read as evidence
of novelty** — a third of it is recognisable to Pfam.

## 6. Two defects found and fixed during this run

Both were caught by consistency checks, not by the pipeline failing.

### 6.1 Duplicate protein records (0.22%)

`data/plasann_run/gbk/` carries **two shard series** (`shard_*` and `mshard_*`), and **144 plasmids
appear in both**. Globbing both series is mandatory — globbing one drops ~27% of the data, a pitfall
already recorded in `reports/small_cryptic_methodology.md` — but it emits the overlap twice.

Detected as 845 duplicate FASTA headers (max multiplicity 2). Effect was small: dark ORFs
379,397 → **378,552**; recurrence rate 84.75% → 84.70%; 91 families were "recurrent" only because a
protein was counted against itself. Fixed in `scripts/cluster_dark_orfs.sh` with a
first-occurrence-wins dedup on the header (records are exactly two lines, so this is safe and
deterministic given the sorted glob). All numbers in this report are post-dedup, and
`reports/dark_orf_clustering.md` has been updated.

### 6.2 CDS indexing broke the dark↔all-CDS join (10.7%) — the serious one

`scripts/extract_plasann_proteins.py` originally numbered CDS **within the emitted subset**. In
`--mode dark` the counter advanced only on dark CDS; in `--mode all` it advanced on every CDS. The
same physical protein therefore received **different ids in the two output files** whenever its
plasmid also carried a named CDS.

Nothing failed. Ids stayed well-formed, the files looked correct, and the join silently dropped
**40,572 of 378,552 dark ORFs (10.7%)**, all of which were then miscounted as "no Pfam match".

It surfaced only because two independently computed totals disagreed — the family-size-weighted
Pfam-match figure (39.4%) against the cross-tab's row sums (35.5%). The corrupted cross-tab read
59.2 / 26.8 / 8.7 / 5.3; the correct one is **57.4 / 28.5 / 10.9 / 3.2**.

Fixed by advancing the index on **every** translated CDS regardless of mode, so dark ids are a
strict subset of all-CDS ids. Validated explicitly — 378,552 ids on each side, zero mismatches in
either direction — and the notebook now carries an `assert` on that join so the failure can never
again be silent.

## 7. Limitations

1. **No Pfam hit is not proof of novelty.** Profile HMMs miss remote homology that structure
   comparison finds. The honest next step is `foldseek` against the AlphaFold DB; until then the
   claim is bounded by "invisible to Pfam-A 38.2 at its gathering thresholds".
2. **Median Pfam alignment covers ~52% of the dark protein**, so many rescues explain only part of
   the sequence.
3. **Representative propagation is 96.1% concordant**, not 100% (§4.5).
4. **Scope classification (§4.3) is a regex over Pfam family names**, written by us. It is a
   reasonable mapping onto PlasAnn's categories, not an authoritative ontology; the split between
   "inside" and "outside" would shift by a few points under a different pattern list.
5. **One clustering threshold** upstream (30% identity / 80% coverage); untested at 40% / 50%.
6. **`mob_cluster` is a mash-based proxy, not a phylogeny**, and the subset is lineage-skewed
   (`AA379` holds 19,427 of 71,414 plasmids).
7. **PlasAnn's ORF calls are inherited, not verified** — no independent gene calling was run.
8. Per-gene coverage is 92.2%; the 562-plasmid CARD/PlasAnn AMR definitional leak carries through
   unchanged (0.24%, affects nothing here).

## 8. Outputs

Under `data/pfam_run/` and `data/refs/pfam/`:

| file | what |
|---|---|
| `data/refs/pfam/Pfam-A.hmm{,.h3f,.h3i,.h3m,.h3p}` | Pfam-A 38.2, pressed (~5.0 GB) |
| `chunks/q*.faa` | 88 query chunks, 174,827 sequences |
| `query_manifest.txt` | set sizes, chunk count, seed |
| `domtbl/q*.domtbl` | raw `hmmsearch --domtblout`, 88 files |
| `pfam_per_seq.tsv` | best Pfam hit per query sequence |
| `darkfam_pfam.tsv` | per dark family: size + Pfam status + family hit |
| `dark_novel_families.tsv` | **the deliverable** — 3,343 substantial Pfam-negative families |
| `pfam_summary.txt`, `pfam_verdict.txt` | the printed summaries |
