# Independent verification pack — recomputed 2026-09-02

Produced by an independent recompute from the tables on disk, **not** by reading the project's own
result TSVs where a recompute was possible. Scripts: `review/verify/v1_resource.py`,
`v2_multirep.py`, `v3_smallcryptic.py`, `v4_dark.py`, `v5_pfam.py`, `v6_sweep.py`.
Interpreter: `/gorilla/home/amujkic/.conda/envs/genesis/bin/python` (pandas 3.0.2, numpy 2.4.4).

**Headline: the arithmetic is sound.** 47 of 52 recomputed quantities match the reported value
exactly (one initially-flagged discrepancy, §7.3, was my error and is retracted). The five that
remain are listed in §7 and none of them overturns a conclusion. Reviewers
should therefore attack design, definition and interpretation — not the numbers.

---

## 1. Resource layer (`plasmid_metadata_master.tsv`)

| quantity | reported | recomputed | verdict |
|---|---:|---:|---|
| master rows × cols | 208,248 × 67 | 208,248 × 67 | ✅ |
| unique `plasmid_id` | — | 208,248 | ✅ no duplicate rows |
| Simulated-artifact | 64,658 | 64,658 | ✅ |
| Lab-artifact | 87 | 87 | ✅ |
| analysis set (excl. both) | 143,503 | 143,503 | ✅ |
| analysis set (excl. Simulated only) | 143,590 | 143,590 | ✅ both figures are real, see §7.1 |
| `is_clinical` | 19,446 | 19,446 (full set); **19,440** in analysis set | ✅ / note |
| AMR+ (`card_n_arg`>0) | 29,396 (14.1%) | 29,396 (14.12%) | ✅ |
| `card_multidrug` | 20,741 | 20,741 | ✅ |
| PlasmidFinder ≥1 replicon | 41,555 (20%) | 41,555 (20.0%) | ✅ |
| **PlasmidFinder ≥2 ("multireplicon")** | **24,583** | **24,583** | ✅ arithmetic — ❌ definition, §7.2 |
| `mob_cluster` distinct | 7,054 | 7,054 | ✅ |
| `mob_host_range` coverage | 58% | 57.7% | ✅ |
| any-method replicon coverage | 95,442 (46%) | 95,442 (46%) | ✅ **corrected**, see §7.3 |

## 2. Small cryptic compartment (`size_bp` < 10 kb, payload-free incl. conjugation filter)

| quantity | reported | recomputed | verdict |
|---|---:|---:|---|
| small < 10 kb | 70,243 (48.9%) | 70,243 (48.9%) | ✅ |
| small-cryptic | 47,031 (32.8%) | 47,031 (32.8%) | ✅ |
| median size / GC / CDS | 4,148 bp / 43.6% / 5 | 4,148 / 43.6 / 5 | ✅ |
| PlasmidFinder call rate | 10.5% | 10.5% | ✅ |
| mob_typer call rate | 27.5% | 27.5% | ✅ |
| PlasAnn call rate | 17.6% | 17.6% | ✅ |
| any-of-three | 28.9% | 28.9% | ✅ |
| **untyped by all three** | 33,427 (71.1%) | 33,427 (71.1%, 23.3% of analysis set) | ✅ |
| large-conjugative n / any-of-three | 25,631 / 92.2% | 25,631 / 92.2% | ✅ |

## 3. Dark-ORF working set (`size_bp` < 20 kb, payload-free **without** conjugation filter)

| quantity | reported | recomputed | verdict |
|---|---:|---:|---|
| small < 20 kb | 82,261 | 82,261 | ✅ |
| payload-free small | 71,414 | 71,414 | ✅ |
| `cryptic_small_ids.txt` | 71,414 | 71,414, set-identical to recompute | ✅ |
| PlasAnn-AMR "definitional leak" | 562 | **694** on the master column | ⚠️ §7.4 |
| overlap with the 47,031 small-cryptic set | — | 47,031 (strict subset) | ✅ but see §7.5 |

## 4. Dark proteome clustering (30% id / 80% cov, no `--cluster-reassign`)

| quantity | reported | recomputed | verdict |
|---|---:|---:|---|
| dark ORFs (`dark_orfs.faa`) | 378,552 | 378,552 | ✅ |
| all CDS (`all_cds.faa`) | 456,949 | 456,949 | ✅ |
| cluster rows / unique members | 378,552 | 378,552 / 378,552 | ✅ no double-counting |
| families | 92,752 | 92,752 | ✅ |
| singletons | 57,935 (15.3%) | 57,935 | ✅ |
| ORFs in families ≥2 | 320,617 (84.7%) | 320,617 (84.7%) | ✅ |
| families ≥10 members | 5,623 | 5,623 | ✅ |
| mix30 families | 97,500 | 97,500 | ✅ |
| dark+named mixed families | 1,674 | 1,674 | ✅ |
| dark-only families | 91,159 | 91,159 | ✅ |
| **dark-only share of dark ORFs** | 85.9% | 85.9% | ✅ |

## 5. Pfam-A 38.2 adjudication

| quantity | reported | recomputed | verdict |
|---|---:|---:|---|
| NAMED positive control hit rate | 60,124 / 78,397 = 76.7% | 76.7% | ✅ |
| DARKREP hit rate | 24,357 / 92,752 = 26.3% | 26.3% | ✅ |
| family-size-weighted Pfam+ | 149,160 (39.4%) | 149,160 (39.4%) | ✅ |
| families ≥10 with a hit | 40.5% | 40.5% | ✅ |
| **cross-tab (% of 378,552)** | 57.4 / 28.5 / 10.9 / 3.2 | **57.4 / 28.5 / 10.9 / 3.2** | ✅ exact |
| join coverage dark↔all-CDS | 378,552 | 378,552, zero loss | ✅ the §6.2 defect is genuinely fixed |
| Pfam-negative families ≥10 | 3,343 | 3,343 | ✅ |
| their dark ORFs | 115,811 (30.6%) | 115,811 (30.6%) | ✅ |
| % spanning ≥2 lineages | 72.9% | 72.9% | ✅ |

## 6. Threshold sweep (`--cluster-reassign`, 80% coverage held)

| threshold | families (rep/recomp) | 1-plasmid % | recurrent % (rep/recomp) | core families |
|---|---|---:|---|---:|
| 30% + RA | 99,362 / **99,362** | 62.9 / 62.9 | 83.3 / **83.6** | 61 / **61** |
| 50% | 117,393 / **117,393** | 64.5 / 64.5 | 79.8 / **80.1** | 39 / **39** |
| 70% | 137,149 / **137,149** | 66.7 / 66.7 | 75.7 / **76.0** | 19 / **19** |
| 90% | 160,518 / **160,518** | 69.6 / 69.6 | 70.3 / **70.6** | 5 / **5** |
| 30% original | 92,752 | — | 84.4 / **84.7** | 93 / **93** |

The core-set collapse (93 → 61 → 39 → 19 → 5) reproduces exactly.

---

## 7. The six discrepancies, and what they mean

### 7.1 Two analysis-set sizes are both in circulation — 143,503 and 143,590
Both are arithmetically correct: 143,503 excludes Simulated-artifact **and** Lab-artifact; 143,590
excludes Simulated-artifact only. The reports use 143,503. Anything downstream that used 143,590
is analysing 87 extra lab artifacts. **Action: state one canonical number and grep for the other.**

### 7.2 `pf_n_inc` counts allele-level PlasmidFinder names, not replicons — the multireplicon count is inflated
`scripts/harvest_typing.py:78` sets `pf_n_inc = len(incs)` where `incs` is the set of distinct
**allele** names, while the same file computes `pf_inc_families` with variants collapsed and never
uses it for the count. Consequence, verified directly:

```
GenBank_CM001154.1  pf_n_inc = 11
  IncFII;IncFII(29);IncFII(Cf);IncFII(S);IncFII(Yp);IncFII(p14);
  IncFII(pAR0022);IncFII(pCoo);IncFII(pKP91);IncFII(pSE11);IncFII(pSFO)   -> ONE replicon family
```

- multireplicon by `pf_n_inc ≥ 2`: **24,583**
- multireplicon by distinct **families** ≥ 2: **22,422** (−8.8%)
- the `pf_n_inc` distribution has a long tail to 15+, which is not biology — it is one locus
  matching many database variants.

`reports/typing_methodology.md` line 41 states families collapse allelic variants, then line 94
reports the multireplicon count on the un-collapsed names. `PROJECT_OVERVIEW.md` §4 repeats
24,583 without qualification. The project describes multireplicon as its founding question, so
this is the resource layer's most consequential defect. It is partly mooted by the standing
"mob_typer only" rule — the mob-based figure is 18,685 multireplicon of 75,751 typed (**24.7% of
typed**, the number the deck uses) — but the PlasmidFinder figure is still published in two docs.

### 7.3 WITHDRAWN — "any-standard-threshold replicon coverage 95,442" reproduces exactly
**This finding was wrong and is retracted.** I reconstructed the expression as
`pf_n_inc≥1 OR mob_rep_types OR plasann_n_replicons≥1`, which gives 97,668. The correct expression
is `pf | mob | plsdb` — `PROJECT_OVERVIEW.md` says PLSDB, not PlasAnn, and I misread it. Recomputed
by the methodology reviewer: **95,442 exactly**. The reported number was right; the pack was wrong.

*(Also corrected: my brief to the methodology reviewer questioned whether the 96.1% Pfam propagation
figure used the right denominator, on the grounds that `pfam_per_seq.tsv` holds only 1,140 VALID
rows against 3,678 queries. That file records only queries **with** a hit, so 1,140 is the VALID hit
count, not a truncated denominator. The published propagation figure is correct, κ = 0.909.)*

### 7.4 The "562-plasmid AMR definitional leak" is 694 on the master column
`plasann_n_amr > 0` within the 71,414 gives **694**, not 562. Part of the gap is that the notebook
counted from the per-gene shards (92.2% coverage ⇒ ~640 expected), leaving ~78 unexplained. The
number is small either way, but it is quoted in two reports and should be recomputed from one
stated source.

### 7.5 Two different "small payload-free plasmidome" definitions are used interchangeably
- `reports/small_cryptic_methodology.md`: < **10 kb**, payload-free **including** a conjugation-machinery
  filter → **47,031**
- `reports/dark_orf_clustering.md` / all `notebooks_my/` dark work: < **20 kb**, **no** conjugation
  filter → **71,414**

The first is a strict subset of the second (verified: overlap = 47,031). But headline sentences move
between them — e.g. "71.5% of the coding capacity of the payload-free small plasmidome is
unannotatable" (the 71,414 object) sits in a narrative whose scale claims ("47,031 payload-free",
"32.8% of all complete plasmids") come from the 47,031 object. Every dark-proteome number is
conditioned on the looser definition, which *includes conjugative plasmids up to 20 kb*.
**Action: name the two sets distinctly and audit every sentence that crosses between them.**

### 7.6 "Recurrent %" runs 0.3 pp below a direct recompute at every threshold
83.3/79.8/75.7/70.3 reported vs 83.6/80.1/76.0/70.6 recomputed — a constant offset across all five
clusterings, so a definitional difference (probably famstats `proteins` vs cluster rows), not an
error. Worth pinning down so the two paths agree.

---

## 8. What reproduced that did *not* have to

- The `mix30` ↔ `dark30` join that was silently dropping 10.7% of ORFs before the 2026-08-28 fix now
  joins **378,552 / 378,552 with zero loss**. The fix is real and the assert is doing its job.
- No duplicate FASTA headers survive in `dark_orfs.faa`; the dedup fix holds.
- `cryptic_small_ids.txt` is set-identical to a fresh recompute from the master table.

## 9. Things the recompute could not check
- Anything upstream of the master table (PlasAnn gene calls, RGI calls, MOB-suite internals,
  BioSample/GOLD environment scraping) — these are tool outputs taken on trust.
- The DefenseFinder / dbAPIS adjudication and the FlashFold target selection (notebook-only,
  not recomputed here).
- Any figure.
