# How much of the dark-plasmidome result depends on the clustering threshold?

**Date:** 2026-09-01
**Plan:** `plans/conservative_reclustering.md`

## Scripts that produced this report

| Step | Script |
|---|---|
| Re-clustering (one threshold per call) | `scripts/recluster_dark_orfs.sh` |
| Sweep submission | `scripts/recluster_sweep.sbatch` (job 6595632) |
| Criterion-compliance audit | `scripts/audit_cluster_compliance.py` |
| Per-family stats | `scripts/recluster_family_stats.py` |
| Pfam search over the new representatives | `scripts/recluster_pfam_array.sbatch` (job 6595640) |
| Pfam roll-up per threshold | `scripts/recluster_pfam_rollup.py` |

Outputs: `data/dark_orf_run/recluster/dark{30,50,70,90}_{cluster.tsv,rep_seq.fasta,famstats.tsv,fampfam.tsv}`,
`data/pfam_run/recluster/`. **Nothing under `data/dark_orf_run/dark30_*` or `mix30_*` was
modified**, so `reports/dark_orf_clustering.md`, `reports/pfam_dark_validation.md` and the three
notebooks that read them remain reproducible.

## 1. Why this was done

The dark proteome had been clustered once, at `--min-seq-id 0.3 -c 0.8 --cov-mode 0`
(`scripts/cluster_dark_orfs.sh`). Two defects:

1. **Single threshold** — limitation #2 of `reports/dark_orf_clustering.md`.
2. **Unenforced criterion** — MMseqs2's help states cascaded clustering "can cluster sequence
   that do not fulfill the clustering criteria", corrected by `--cluster-reassign`, which was
   not passed.

## 2. What the literature does

| Study | Tool | Identity | Coverage | Min size |
|---|---|---|---|---|
| UniRef50 (Suzek et al., *Bioinformatics*) | CD-HIT | 50% | 80% overlap with seed | — |
| UniRef90 | CD-HIT | 90% | 80% | — |
| FESNov (Rodríguez del Río et al., *Nature* 2024) | MMseqs2 | 30% | **50%**, `--cov-mode 1` | ≥3 |
| NMPFamsDB (*NAR* 2024) | HipMCL on SSN | 70% | 70% aln length | ≥100 |
| PlasX / MobMess (*Nat Microbiol* 2024) | MMseqs2 | **sweep 0.05–0.9** | defaults | — |
| **this project (original)** | MMseqs2 | 30% | **80%, `--cov-mode 0`** | none |

Two conclusions. The coverage setting was **already stricter than any of them** — 80% of *both*
sequences, matching UniRef's rule (introduced 2013 to stop partial-sequence merging) and
stricter than the *Nature* 2024 novel-family catalogue. And no comparable study picks a single
identity threshold; PlasX ran twelve. So the sweep is on identity, with coverage held at 80%.

## 3. The criterion was genuinely not being enforced

`scripts/audit_cluster_compliance.py` samples 200 multi-member families per run and aligns every
member back to **its own representative** — no transitive links.

| clustering | member-rep pairs | compliant | fails identity | fails coverage |
|---|---:|---:|---:|---:|
| **30%, original (no reassign)** | 1,066 | **89.12%** | 5 | 111 |
| 30% + `--cluster-reassign` | 1,619 | **99.26%** | 1 | 11 |
| 50% + reassign | 944 | 99.58% | 3 | 1 |
| 70% + reassign | 1,286 | 99.92% | 1 | 0 |
| 90% + reassign | 719 | **100.00%** | 0 | 0 |

Roughly **one member in nine** of the published clustering did not meet the criterion the
methods claim, almost always on coverage. Adding the flag at the same 30% identity yields
**+6,610 families (+7.1%)** and moves **12.00% of members** to a different family.

This does not put homology in doubt — an audit of the largest family found its worst
member-to-representative E-value was 1.2e-43 — but the enforced criterion was looser than the
documented one.

## 4. Sensitivity of every headline number

378,552 dark ORFs on 63,993 carrier plasmids throughout. "Core" = Pfam-negative, on ≥100
plasmids, spanning ≥10 MOB lineages (the target set defined in `dark_families.ipynb`).

| | 30% orig | 30%+RA | 50% | 70% | 90% |
|---|---:|---:|---:|---:|---:|
| families | 92,752 | 99,362 | 117,393 | 137,149 | 160,518 |
| families on 1 plasmid | 63.0% | 62.9% | 64.5% | 66.7% | 69.6% |
| **dark ORFs in recurrent families** | **84.4%** | **83.3%** | **79.8%** | **75.7%** | **70.3%** |
| families holding 50% of dark ORFs | 3,009 | 4,177 | 6,612 | 10,186 | 15,985 |
| **dark ORFs in Pfam-named families** | 39.4% | 39.1% | 38.6% | 38.9% | 38.9% |
| **core families** | **93** | **61** | **39** | **19** | **5** |
| core, % of dark ORFs | 5.9% | 4.0% | 2.2% | 1.2% | 0.2% |
| most prevalent family: plasmids | 2,468 | 2,325 | 998 | 878 | 866 |
| most prevalent family: lineages | 88 | 83 | 8 | 3 | 3 |
| most prevalent family: Pfam | none | none | named | named | named |

### What survives

- **The dark proteome recurs.** Even at UniRef90 — 90% identity, 80% mutual coverage — **70.3%**
  of dark ORFs sit in families found on two or more plasmids. The central claim of
  `dark_plasmidome.ipynb` §1 is threshold-independent.
- **The Pfam verdict is flat.** The share of the dark proteome in Pfam-named families is
  **38.6–39.1% at every threshold** (39.4% originally). The artifact-versus-real conclusion of
  `reports/pfam_dark_validation.md` does not depend on the clustering.
- **The one-off majority.** 63–70% of families sit on a single plasmid at every threshold.

### What does not survive

- **The core target set.** 93 families at 30% → **61** once the criterion is actually enforced
  → **39** at UniRef50 → **5** at UniRef90. A third of the loss comes from `--cluster-reassign`
  alone, before any threshold change.
- **The most prevalent family.** `IMGPR_plasmid_3300022589_000018|2` is a 30% object. Its ORF
  falls into a family of 2,325 proteins / 83 lineages at 30%+RA, **89 proteins / 5 lineages** at
  50%, and **5 proteins / 1 lineage** at 70% and 90%. The "88 MOB lineages" figure is a property
  of the threshold, not of a protein.
- **The identity of the top family.** From 50% upward the most prevalent dark family is one Pfam
  *can* name — i.e. an annotation failure, not novel biology.

## 5. What this means for the existing notebooks

`dark_plasmidome.ipynb` §1–5 and the Pfam verdict stand as written. `dark_families.ipynb` §3–5
are threshold-conditional and should be read with this table beside them: the 93-family core is
an upper bound, and 61 (30% + enforced criterion) or 39 (UniRef50) are the defensible numbers.
The `widespread_dark_orfs.ipynb` results are unaffected — they use exact sequence identity, not
clustering.

## 6. Limits of this analysis

- One coverage setting (80%, `--cov-mode 0`). Coverage was not swept, because nothing in the
  literature reviewed uses a stricter value.
- The Pfam roll-up pools two hmmsearch runs (`data/pfam_run/domtbl/` and
  `data/pfam_run/recluster/domtbl/`). Same tool, same Pfam-A 38.2, same `--cut_ga`, so they are
  poolable, but they were run on different days.
- The compliance audit samples 200 families per run, not all of them.
- The big-plasmid dark proteome (5,135,708 ORFs, `data/fam_3300022589/big_dark_orfs.faa`) was
  **not** included. Task T7 of the plan.
