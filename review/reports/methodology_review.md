# Internal rigour audit — Reviewer 1 (research design, statistics, code correctness, reproducibility)

**Date:** 2026-09-02 · **Scope:** methods, not arithmetic (per `review/VERIFICATION_PACK.md`)
**Interpreter used for every recompute below:** `/gorilla/home/amujkic/.conda/envs/genesis/bin/python`
(3.14.4, pandas 3.0.2, numpy 2.4.4); MMseqs2 18.8cc5c from
`/gorilla/home/amujkic/.conda/envs/panaroo/bin/mmseqs`.
**Iron rule observed:** no project file was modified. Scratch scripts live under `/tmp` and the
session scratchpad.

## Summary

The pipeline is unusually well-built for an academic codebase: every join is strictly 1:1, the
CDS-indexing fix is real, the dedup is provably safe, RGI covered 208,245/208,248 inputs, and the
propagation check the brief flagged as suspect is in fact correctly denominated (κ = 0.909). The
defects that matter are not arithmetic but *definitional*, and they share one signature: **counting
annotation labels instead of biological entities**, then describing the count in prose the code does
not implement. `pf_n_inc` counts PlasmidFinder alleles; `card_multidrug` calls a plasmid multidrug on
one gene; `pfam_verdict.py`'s regex misfiles ≥7.5% of rescued ORFs as "novel". Three statistical
claims do not survive recomputation: "near-parity" in dispersal, "61% of the conjugative effect is
size", and "identical criteria to CGE PlasmidFinder". Each is fixable without new data.

---

## Findings

Severity reflects effect on a published conclusion, not code aesthetics.

### CODE DEFECTS FOUND

---

**[MAJOR] C1 — `pf_n_inc` counts PlasmidFinder alleles; I traced the inflation to its source**

*Anchor:* `scripts/harvest_typing.py:80` (`"pf_n_inc": len(incs)`) vs `:77` (`fams` computed and
never used). Confirms pack §7.2 and adds the mechanism:

```
$ python  # parse data/typing/plasmidfinder_blast.tsv.gz sseqid column
distinct sseqid: 481 -> distinct parsed replicon names: 274
distinct Inc FAMILIES after collapsing variants: 206
families with the most allele variants: [('IncFII', 28), ('IncFIB', 17), ('Col', 14), ('IncFIA', 4)]
```

*Wrong:* the CGE database carries **28 IncFII allele variants**, so one IncFII locus can contribute
up to 28 to `pf_n_inc`. `reports/typing_methodology.md:41` states families collapse variants; `:94`
then reports 24,583 multireplicon on the un-collapsed names.

*Fix:* set `pf_n_inc = len(fams)`; keep the allele count as `pf_n_alleles`. Restate 24,583 → 22,422.

---

**[MAJOR] C2 — `card_multidrug` is the same bug in the AMR layer, and nobody has flagged it**

*Anchor:* `scripts/aggregate_card.py:47`, `"card_multidrug": int(len(drug) >= 2)` where `drug` is a
union over CARD's `Drug Class` **label strings**, not over genes.

```
AMR+ plasmids: 29396 ; card_multidrug=1: 20741 (70.6% of AMR+)
  AMR+ with exactly ONE distinct ARO: 12602 ; of those card_multidrug=1: 4402 (34.9%)
  -> 4402 plasmids (21.2% of all 'multidrug' plasmids) are 'multidrug' on a SINGLE gene
```

*Wrong:* a single promiscuous efflux ARO whose CARD record lists ≥2 drug classes makes the plasmid
"multidrug". **21.2% of the 20,741 multidrug plasmids carry exactly one ARG.** Any sentence reading
"multidrug plasmids" as multi-determinant carriage is false for a fifth of them.

*Fix:* keep `card_multidrug` but rename it `card_multiclass_annotation`, and add
`card_multidrug_genes = (card_n_arg_unique >= 2) & (n_drug_classes >= 2)`. Report both.

---

**[MINOR] C3 — `plsdb_n_inc` carries the same allele inflation, inherited from PLSDB**

*Anchor:* master column `plsdb_n_inc` (delimiter is `|`, not `;`).
`plsdb typed 24872 ; multireplicon by plsdb_n_inc>=2: 9155 ; by collapsed families>=2: 8923 (-2.5%)`
Small (2.5%) but it is the same class and should be documented in
`reports/master_table_data_dictionary.md`.

---

**[MAJOR] C4 — `pfam_verdict.py`'s scope regex is high-precision / low-recall, and the recall failure
biases the headline toward novelty**

*Anchor:* `scripts/pfam_verdict.py:23-38` (`PLASANN_SCOPE`, `scope_of`).

Precision is fine — I checked: the unanchored `re.I` toxin pattern picks up only 19 ORFs of genuine
false positives (`Endotoxin_N/C`, `Toxin_10`, `Toxin_YhaV`), and **no** Pfam name matches more than
one scope pattern, so first-match-wins never silently arbitrates. Recall is not fine:

```
AS PUBLISHED: inside = 70,575 / 149,160 = 47.3% of rescued ORFs
unambiguous plasmid-core Pfam names the regex misclassifies as 'outside': 1,136 families, 11,126 ORFs
  AbiEii Antitox_RHH DnaB_C Fic FtsK_SpoIIIE HicB HicB_lk_antitox KfrA_N PIN PIN_3 ParD_antitoxin
  ParG PhdYeFM_antitox RHH_1 RelE RepC Rop TacA1 Toprim_2 TraG-D_C TraX TrbL TrfA TrwC VbhA
  Zn_ribbon_DnaG wHTH_Par
=> corrected 'inside PlasAnn's scope' LOWER BOUND: 54.8%  (published 47.3%)
```

*Wrong:* `RelE` is absent while `RelB` is present — the toxin half of the canonical *relBE* system is
filed as novel. So are `TrfA` and `RepC` (replication initiators), `TrwC` (relaxase), `Rop` (ColE1
copy control), `KfrA_N`, `ParG`, `TraX`, `TrbL`. `reports/pfam_dark_validation.md:110` publishes
"47.3% of rescued dark ORFs … inside a category PlasAnn already annotates" — the true figure is
**≥54.8%**, and the "outside PlasAnn's categories" residue shrinks correspondingly. The error runs in
the direction that flatters the novelty claim.

*Fix:* replace string patterns with a Pfam **clan** map (`Pfam-A.clans.tsv`) plus a curated
accession list; regex on names cannot be made complete. At minimum add the 27 names above and
re-publish §4.3.

---

**[MINOR] C5 — `re.split` positional `maxsplit` is a hard deprecation on the project's own interpreter**

*Anchor:* `scripts/harvest_typing.py:47`, `re.split(r"_\d", s, 1)[0]` (the only occurrence in
`scripts/*.py`).

```
$ genesis/bin/python -W error::DeprecationWarning -c "import re; re.split(r'_\d','X_1',1)"
3.14.4
RAISES: DeprecationWarning 'maxsplit' is passed as positional argument
```
It becomes a `TypeError` in 3.15. *Fix:* `maxsplit=1`. (The parse itself is clean — all 481 database
sseqids resolve to 274 well-formed replicon names, zero truncated or unbalanced.)

---

**[MINOR] C6 — the protein extractor discards CDS coordinates, which makes the one artefact test that
matters impossible**

*Anchor:* `scripts/extract_plasann_proteins.py:70` detects `"     CDS "` purely as a state reset and
never parses the location string; only `{plasmid_id}|{cds_index}` survives.

*Wrong:* with no start/end/strand, nothing downstream can ask whether a "dark ORF" overlaps another
CDS, sits antisense to a named replication gene, or falls inside an oriV/RNAI region — the leading
artefact hypotheses (see §D15). *Fix:* emit a sidecar
`data/dark_orf_run/cds_coords.tsv (plasmid_id, cds_index, start, end, strand)`. The parser already
sees the line; capturing it is three lines of code and no re-run of PlasAnn.

*Note on the index fix itself:* the documented fix is real and internally consistent. The docstring
claims `cds_index` "counts EVERY translated CDS"; in fact it counts CDS having **both** a
`/category=` and a `/translation=` (`:64`, `:79`). Because both modes apply the identical rule, the
`dark`↔`all` join is still exact (pack §5 confirms 378,552/378,552, zero loss). The docstring is
wrong, not the code.

---

**[MINOR] C7 — the dedup is safe *in this data*, but nothing in the code makes it safe**

*Anchor:* `scripts/cluster_dark_orfs.sh:16-22`, `dedup() { awk '/^>/{keep = !seen[$0]++} keep' }`.

I tested the claim directly on every plasmid present in both PlasAnn shard series:

```
mshard plasmids: 16118  shard plasmids: 49834  IN BOTH: 147
records for the 147 duplicated plasmids: shard 993, mshard 993
headers in both: 993 ; only in shard_*: 0 ; only in mshard_*: 0
SAME header but DIFFERENT sequence: 0
plasmids where the two runs called a DIFFERENT number of CDS: 0 of 147
same plasmid|idx but DIFFERENT D/N category between runs: 0
```

*Verdict:* first-occurrence-wins is **empirically safe** — the two runs are byte-identical. But the
safety is a property of the data, not the algorithm: if a reshard ever called a different number of
CDS, the surplus higher-index records from the losing run would survive and the plasmid would carry a
chimera of two annotation runs, invisibly. The inline comment also says "144 plasmids"; the true
overlap is **147** in `faa_all`.

*Fix:* after dedup, assert `sort -u` on `(header, sequence)` yields the same line count as `sort -u`
on `header`; fail loudly otherwise. Correct 144 → 147.

---

**[MINOR] C8 — `audit_cluster_compliance.py` merges unguarded**

*Anchor:* `:79`, `own = sub.merge(h, on=["member","rep"], how="left")`. If MMseqs2 returns more than
one row per (member, rep) pair the denominator silently multiplies. *Fix:*
`assert len(own) == len(sub)`, or `.drop_duplicates(["member","rep"])` on `h` first (I had to add
exactly this to reproduce the audit).

*Second-order:* `(~aligned)` conflates "no homology" with "MMseqs2 prefilter miss at 30% identity".
That is conservative for the compliance estimate, but should be stated.

---

**[MINOR] C9 — `card_*` NaN→"0" fill in the master is a real-zero fill; verified safe today, unasserted**

*Anchor:* `scripts/build_metadata_master.py:235`,
`ca.get(k, "0" if k in card_num else "")` — a plasmid RGI never processed would silently become
AMR-negative. I checked whether that ever happens:

```
$ grep -h "done=" data/card_run/logs/*.out | awk ...
shards with a done line: 595  total inputs: 208245  nonzero rc: 0
```

208,245 of 208,248 inputs, every shard rc=0. The fill is safe. *Fix:* make it explicit — assert
`len(card_input_ids) == len(master)` in the build, so a future partial rerun cannot manufacture
14.1%-denominator zeros. The parallel `plasann_*` block at `:231` correctly leaves counts blank; only
the CARD block coerces.

---

**[MINOR] C10 — `num()` in the small-cryptic notebook makes unannotated plasmids automatically "cryptic"**

*Anchor:* `scripts/build_small_cryptic_notebook.py`, `def num(s): return
pd.to_numeric(s, errors="coerce").fillna(0)` feeding `real["cryptic"]`.

```
plasann_annotated==0 inside the analysis set : 71
small & PlasAnn-unannotated CALLED cryptic purely by fillna(0): 67
```

67 of 47,031 (0.14%) — negligible in magnitude, but those same 67 also enter §5's "73.1% of plasmids
with **every** CDS unannotatable" with zero CDS. *Fix:* drop rows with `plasann_annotated == 0` from
the cryptic call rather than zero-filling them, and say so.

---

**[MINOR] C11 — every master join is dict-keyed, so a duplicate key would be a silent last-wins overwrite**

*Anchor:* `scripts/build_metadata_master.py:89-98`. I verified all six inputs are currently 1:1:

```
analysis_table   208248 rows / 208248 keys   1:1
typing_authoritative 208245 / 208245        1:1
card_per_plasmid  29396 /  29396            1:1
plasann_features 208177 / 208177            1:1
plsdb_enrichment  45743 /  45743            1:1
environment_reconciled 208248 / 208248      1:1
```

**No row multiplication, no NaN-fill pathology, no fan-out anywhere in the master build.** The design
(dict lookup, `extrasaction="ignore"`, one output row per `analysis_table` row) is structurally
incapable of multiplying rows — this is better than a pandas `merge` chain would have been. *Fix:*
add a one-line duplicate-key assert per input so the guarantee is enforced rather than observed.

---

**[MINOR] C12 — the Pfam propagation check fails open**

*Anchor:* `scripts/summarize_pfam.py:110-119`. `_valid_targets` re-reads `data/pfam_run/chunks/*.faa`;
if those are ever cleaned up, `allv` is empty and `if len(allv):` skips the entire validation
**silently**, leaving `pfam_summary.txt` with no propagation section and no warning. *Fix:* `raise`
instead of skipping, or persist the VALID id list next to `query_manifest.txt`.

---

### DESIGN / STATISTICAL OBJECTIONS

---

**[MAJOR] D1 — the PlasmidFinder call rule is not "identical criteria to CGE", and the founding
number moves 65% at the tool's own default**

*Anchor:* `scripts/harvest_typing.py:45` (`pid >= 80 and length/slen >= 0.60`);
`reports/typing_methodology.md:27-28` "**Method — identical criteria to CGE PlasmidFinder / PLSDB**".

CGE PlasmidFinder's default is **95% identity / 60% coverage**; 80% is the lowest selectable setting
(and the ABRicate/PLSDB convention). I re-derived the whole layer from the raw BLAST table at both
thresholds:

```
PlasmidFinder call rule sensitivity (coverage held at 60%):
  pident >= 80 (as run)      >=1 call  41555   multirep by ALLELE 24583   by FAMILY 22422
  pident >= 95 (CGE default) >=1 call  25534   multirep by ALLELE  8574   by FAMILY  8169
```

*Wrong:* the headline "PlasmidFinder typed 41,555 (20%)" becomes 25,534 (12.3%) at the tool's
default, and multireplicon collapses **24,583 → 8,574 (−65%)**. Combined with C1, the defensible
multireplicon count spans 8,169–24,583 depending on two undocumented choices, for a quantity the
project calls its founding question. No sensitivity is reported anywhere.

*Fix:* attribute the 80% threshold to PLSDB (with citation), drop "CGE", and publish the 95%
sensitivity beside every PlasmidFinder headline.

---

**[MAJOR] D2 — the rarefaction counts "Unlabelled" and "Geography only (no habitat)" as habitats,
differentially between the two arms**

*Anchor:* `scripts/build_small_cryptic_notebook.py:364-375` — `sub = df[... df[field].notna()]` then
`len(np.unique(...))`. `hab_sub` has **zero NaN**; missingness is encoded as the *values*
`Unlabelled` (29,283) and `Geography only (no habitat)` (9,293), which therefore count as distinct
habitats.

```
pseudo-habitat share:  small-cryptic 15.8%   large-conjugative 29.1%

=== K=10 habitats ===
AS PUBLISHED (pseudo-habitats counted)  nA=623 nB=359 medA=4.58 medB=5.00 ratio=0.917 p=4.85e-04
pseudo-habitats treated as MISSING      nA=505 nB=257 medA=4.16 medB=4.88 ratio=0.853 p=6.84e-05
```

*Wrong:* the conjugative arm carries **1.8× more** pseudo-habitat records, so the artefact inflates
the comparator's breadth. Correcting it moves the ratio 0.917 → **0.853** and strengthens the
p-value by an order of magnitude. The published analysis therefore *understates* the gap — the error
runs against the paper's own "near-parity" conclusion, which is the honest thing to say and the
reason it must be fixed.

*Fix:* map the two pseudo-categories to NaN before rarefaction and re-run §3 end to end.

---

**[MAJOR] D3 — "near-parity" is contradicted by the K-ladder the report itself publishes**

*Anchor:* `reports/small_cryptic_methodology.md` §3; `smallcryptic_dispersal_Kladder.tsv`.

```
K   metric      ratio
5   habitats    0.929      5   countries   0.878
10  habitats    0.913     10   countries   0.845
20  habitats    0.896     20   countries   0.812
40  habitats    0.800     40   countries   0.793
```

*Wrong:* the report says "the ordering is stable across K ∈ {5,10,20,40}" — true but incomplete. The
**ratio declines monotonically in K on both metrics**, reaching 0.80. That is the signature of a real
breadth deficit that shallow rarefaction cannot resolve, not of parity. Adding D2 on top, the K=10
habitat ratio is 0.853, not 0.917.

Two further framing problems in the same section:
- §7 establishes provenance as "a genuine confounder handled by stratification"; the isolate-only
  stratum **is** therefore the confounder-controlled analysis, and it is null (p = 0.06, 0.11) on
  ~40% the sample size. Reporting the pooled (confounded) result as the headline and the stratified
  (adjusted, underpowered) result as a caveat inverts the evidential hierarchy.
- No confidence interval is attached to the ratio of medians anywhere, and 12 Mann–Whitney tests
  (pooled ×2, ladder ×8, isolate ×2) are reported without a multiplicity statement.

*Fix:* lead with the isolate-only stratum; bootstrap a CI on the ratio (cluster-resampled); state the
K-monotonicity explicitly; replace "near-parity" with "a modest but consistent 8–20% breadth deficit
that widens with sampling depth and is not resolvable in isolate-only data at current n".

---

**[MINOR] D4 — the rank-biserial effect size is reported with the sign inverted**

*Anchor:* `scripts/build_small_cryptic_notebook.py:404`, `rb = 1 - 2*u/(len(a)*len(b))` where `u` is
scipy's U₁ for sample `a` (small-cryptic). Kerby (2014) defines r = 2U₁/(n₁n₂) − 1.

```
habitats : rb(script)=+0.134  rb(Kerby)=-0.134
countries: rb(script)=+0.241  rb(Kerby)=-0.241
```

*Wrong:* `smallcryptic_dispersal_test.tsv` and the report publish **positive** 0.134 / 0.246 while
cryptic breadth is *lower* than conjugative. A reader applying the standard convention reads the
effect backwards. Magnitude is correct. *Fix:* `rb = 2*u/(len(a)*len(b)) - 1`, or state the
convention in the column name (`rank_biserial_conj_minus_cryptic`).

---

**[MAJOR] D5 — "about 61% of the raw conjugative effect is plasmid size" is computed on the wrong
scale, and the quantity is not identified in the first place**

*Anchor:* `manuscript/reassessment_round1.md:117` and `:187`;
`scripts/reassess_review_findings.py:212-218`. I refit both models:

```
conjugative: OR with size 4.77 (beta 1.5632) ; OR without size 12.23 (beta 2.5043)
  reduction on the OR scale   (as published): 61.0%
  reduction on the log-OR scale (correct)   : 37.6%
```

Three separate problems, in increasing order of severity:

1. **Scale.** The 61% is (12.23 − 4.77)/12.23, a percentage change of an odds *ratio*. Effect
   decomposition is linear on the log-odds scale, where the change is **37.6%**. The published figure
   overstates by 1.6×.
2. **Non-collapsibility.** Even with zero confounding and zero mediation, adding a strong outcome
   predictor to a logistic model moves the conditional OR *away* from 1. Here it moves *toward* 1, so
   the confounding/mediation signal genuinely dominates — but the change-in-estimate statistic is
   contaminated and cannot be read as "% of effect explained".
3. **Identification.** Size is partly *downstream* of conjugativity — a T4SS is ~25 kb of DNA — so
   conditioning on it is over-adjustment for a mediator, not confounder control. The reassessment
   acknowledges the ambiguity ("major mediator or confounder") at `:120` and then promotes the number
   to a headline claim at `:187` regardless.

*Fix:* delete the percentage. Report the two ORs as a model sequence, and if a decomposition is
wanted, do it properly — g-computation for a marginal risk-difference contrast, or a natural
direct/indirect effect decomposition with the causal ordering stated.

**Related [MINOR] D5b:** `scripts/reassess_review_findings.py:198-210` uses VIF to adjudicate the
size question. VIF measures variance inflation, not confounding; and VIFs on dummy columns from one
8-level factor are mechanically inflated regardless of any pathology (the reassessment half-notices
this at `:104`). VIF is the wrong instrument for the question being asked.

**Related [MINOR] D5c:** `:222` clusters SEs on `mob_cluster`, but one cluster holds **19,885 of
104,168 observations (19.1%)** across 6,182 clusters. Cluster-robust asymptotics assume many clusters
of comparable size; with this leverage the SEs are likely still anti-conservative. Report a
wild-cluster bootstrap, or a sensitivity dropping the mega-cluster.

**Related [MINOR] D5d:** `:247-252` sets the host-genus reference to `"other"`, a heterogeneous
grab-bag of all non-top-20 genera, so every genus OR is relative to a mixture; and `:128`
`host_raw.str.split().str[0]` maps "uncultured Klebsiella sp." to genus `uncultured`.

---

**[MAJOR] D6 — the Pfam positive control is not like-for-like (though the conclusion survives, which
the report should show rather than assume)**

*Anchor:* `reports/pfam_dark_validation.md:82-83` — 78,397 NAMED at 76.7% vs 92,752 DARKREP at 26.3%,
presented as one table.

Two asymmetries. **Unit:** NAMED is *every* named protein; DARKREP is one *representative per
family*, and 57,935/92,752 (62.5%) of those representatives are singletons — the sequences least
likely to have homologs anywhere. **Length:** Pfam `--cut_ga` hit rate is strongly length-dependent
and the two sets differ (DARKREP median 121 aa, NAMED 185 aa). I ran the missing check:

```
=== Pfam hit rate by protein length ===
             DARK n  DARK %hit  NAMED n  NAMED %hit  gap pp
<50            5913        3.2      825        21.0    17.8
50-74         16006        8.5     7479        74.9    66.3
75-99         15060       13.8    11291        65.2    51.3
100-149       18591       21.4    13956        75.2    53.7
150-199       11697       27.1     8296        55.6    28.5
200-299       12091       40.8    12715        83.7    42.9
300-499        9478       60.9    17766        92.1    31.2
>=500          3916       73.1     6069        80.5     7.4

NAMED crude 76.7% ; NAMED standardised to the DARKREP length distribution 70.6% ; DARKREP 26.3%
```

*Verdict:* **the control holds.** Length standardisation moves NAMED only 76.7 → 70.6, leaving a
44-point gap, and the gap is present in every length stratum. But this is a robustness check the
report needed and did not run, and the like-for-like comparison the report *should* headline is the
family-size-weighted 39.4% vs 76.7%, not 26.3% vs 76.7%.

*Fix:* add the length-stratified table above to §4.1; state the unit asymmetry explicitly; or cluster
NAMED at 30/80 too and compare representative-to-representative.

---

**[MAJOR] D7 — the compliance audit samples families uniformly when the claim is about ORFs**

*Anchor:* `scripts/audit_cluster_compliance.py:52-54`,
`multi = sizes[sizes >= 2].index.tolist(); picked = random.sample(multi, 200)`.

```
multi-member families: 34817 holding 320617 ORFs
family-uniform sample: E[size] =   9.2   (what the script samples)
ORF-weighted          : E[size] = 153.9   (the population the claim is about)
replicating the script's own draw (seed 0, n=200):
  mean size 6.3, max 83 ; represents 1266 of 320617 ORFs = 0.39%
  families >=50 in the draw: 3   (population has 937, holding 40.9% of ORFs)
```

I re-ran the audit with pairs sampled uniformly (probability ∝ family size), same criterion, same
MMseqs2 settings:

```
dark30 ORIGINAL (no --cluster-reassign)          dark30 + --cluster-reassign
 [family] 89.03% compliant (reported 89.12%)      [family] 99.08% (reported 99.26%)
 [orf   ] 88.55% compliant                        [orf   ] 97.62%
   by family size:  2-4  92.4% | 5-9  94.7%         by family size: 2-4 100% | 5-9  95.0%
                   10-49 90.5% | 50+  75.0%                        10-49 93.4% | 50+ 99.1%
```

*Verdict:* the audit **reproduces** (89.03 vs 89.12), and the mechanism I predicted is real —
compliance falls to 75% in families ≥50 in the original clustering — but the ORF-weighted headline
moves only 89.0 → 88.6 and 99.1 → 97.6. So the conclusion ("the criterion was not enforced;
`--cluster-reassign` fixes it") stands. Still: 41% of dark ORFs live in families the draw sees three
times, and the number is used to justify a pipeline-wide switch.

*Fix:* add `--weight orf` to sample member-rep pairs uniformly, and report the by-size-stratum
breakdown. Cost: one MMseqs2 run of a few minutes.

---

**[MAJOR] D8 — the host-range comparison conditions on an outcome-correlated selection, and the
mash-distance control does not address the circularity that matters**

*Anchor:* `reports/small_cryptic_methodology.md` §4 and its *Artefact control*; Limitation 2.

**(a) Circularity, quantified.** `mob_host_range` is essentially a deterministic function of
`mob_cluster`:

```
n with both cluster + host_range: 120210 over 6930 clusters
clusters where host_range is CONSTANT across all members: 5219/6930 (75.3%)
plasmids matching their cluster's modal host range: 91.3%
marginal accuracy of always predicting the GLOBAL modal host range: 18.5%
```

Knowing the cluster predicts the host range 91.3% of the time against an 18.5% base rate. The
mash-distance control conditions on the *focal plasmid's* distance to its neighbourhood; the
mechanism of concern is the *neighbourhood's own taxonomic breadth*, which that control leaves
untouched. Limitation 2 says the control "bounds but does not eliminate" the problem — it does not
bound it either.

**(b) Differential selection, unflagged.** §4 is computed "among plasmids receiving a host-range
call": **35.3% of small-cryptic vs 99.4% of large-conjugative** (the report's own §2 table). The
cryptic plasmids that get a call are precisely the ones with a recognisable mash neighbour — the
least novel third of the compartment — while the conjugative arm is essentially complete. The two
rows of the §4 table are not samples of the same construct.

*Fix:* (i) re-run §4 within mob_clusters containing both groups, which removes the cluster-assignment
channel entirely; (ii) report the host-range table with the 35.3%/99.4% call rates printed in the
same table so the selection is visible; (iii) get an assignment-independent host-range estimate —
e.g. the taxonomy of the source isolate/metagenome — for the subset where it exists.

---

**[MINOR] D9 — the threshold sweep varies one knob and calls it a sensitivity analysis**

*Anchor:* `scripts/recluster_dark_orfs.sh:12-14`; `reports/clustering_threshold_sensitivity.md` §2.

Identity is swept 30→90 with `-c 0.8 --cov-mode 0` and the default `-e` held fixed. The justification
(80% bidirectional coverage is stricter than UniRef and than Rodríguez del Río 2024) is correct and
well-cited. But 80% *bidirectional* coverage is not scale-free: for a 49-aa representative it demands
a ~39-residue alignment, which at `--min-seq-id 0.3` is ~12 identical residues — inside the noise
floor for short proteins, where MMseqs2's cascaded prefilter is also least reliable. Holding coverage
fixed therefore does not test the setting that binds hardest on exactly the short ORFs the project's
core families consist of.

*Fix:* add a two-point coverage arm (`-c 0.5 --cov-mode 1`, `-c 0.9`) at the 30% and 50% identity
levels, and a minimum-length arm (≥50 aa, ≥100 aa as NMPFamsDB uses) so short-ORF sensitivity is
separable from identity sensitivity.

---

**[MINOR] D10 — `geo_country` is not a country field**

*Anchor:* master column `geo_country`, used as the geographic dispersal unit in §3.

```
geo_country: distinct=430 ; entries appearing <5 times: 195
tail sample: TARA_112 TARA_068 TARA_031 ... 'Colorado' 'USA and various oceans'
             'Siachen glacier' 'Antarctic' 'Wellcamp field site' 'Andaman Sea' 'Puerto RIco'
```

430 values for ~195 real countries: TARA station IDs, sub-national units, sea names, field sites, and
at least one typo. *Effect is bounded* — only 115 records carry a `TARA_*` value (67 cryptic, 5
conjugative) — so this does not overturn §3. But "distinct countries per 10 members" is not counting
countries. *Fix:* pass `geo_country` through `geo_utils.ALIASES` + an ISO-3166 whitelist at master
build; route unresolved strings to `geo_locality` and exclude them from breadth.

---

**[MINOR] D11 — `hab_sub`/`geo_country` coverage is misstated in the limitations**

*Anchor:* `reports/small_cryptic_methodology.md` Limitation 1: "`hab_sub` and `geo_country` are
present for a minority of small cryptic plasmids".

```
small-cryptic      geo_country present 83.9% ; hab_sub present 100.0%
large-conjugative  geo_country present 75.9% ; hab_sub present 100.0%
```

Both are majority-present. The real limitation is D2 (the presence is partly pseudo-categorical), and
the ≥K filter — which I checked and which is *not* a problem: it retains 86.5%/87.6% of the two arms
for habitats and 72.0%/65.6% for countries, i.e. comparably. *Fix:* replace Limitation 1 with the
pseudo-category statement.

---

**[MINOR] D12 — item D: how much of the "dark plasmidome" is a gene-calling artefact?**

*Anchor:* `data/dark_orf_run/recluster/core_dark_families_90.tsv` (representative lengths 70, 49,
121, 105, 78 aa); PlasAnn's ORF calls are nowhere independently verified.

**Risk assessment — moderate, and currently untestable from the artefacts on disk.**

Evidence *against* wholesale artefact:

```
dark  n=378552 median=138 mean=190.8 sd=159.9 CV=0.84 P(len<100)=0.35
named n= 78397 median=185 mean=230.6 sd=171.5 CV=0.74 P(len<100)=0.25
   (a pure random-ORF population is geometric: CV -> 1.0)
   dark ORFs < 50 aa: 14,905 (3.9%)
```

The dark length distribution is protein-like, not geometric, and only 3.9% fall below 50 aa. And the
five core families cluster at **≥90% amino-acid identity with 80% bidirectional coverage across
14–42 distinct MOB lineages** — random ORFs on independently evolving DNA cannot do that.

Evidence that leaves it open:
- Conservation of the *protein* is not established, only conservation of a translated frame. A short
  ORF antisense to a conserved *rep* gene, or overlapping an iteron array or a ColE1 RNAI/RNAII
  countertranscript, is conserved at the amino-acid level as a by-product. `Rop` — the ColE1 copy
  control protein — appears in the Pfam rescue, so this neighbourhood is demonstrably present.
- The representatives skew short precisely because MMseqs2 representatives are singleton-enriched
  (62.5%); 5,913 of 92,752 DARKREP are <50 aa and hit Pfam at 3.2%.
- I could not test overlap or strand at all, because C6 discarded the coordinates.

Neighbourhood context I *could* extract (from `all_cds.faa` index adjacency):

```
family (rep)                        len   n  lin | neighbours: named  dark | sole CDS
COMPASS_KP718939.1|8                 70 140   42 |             107   121 |   0/140
IMGPR_plasmid_2563366806_000011|1    49 230   30 |             140   199 |   0/230
IMGPR_plasmid_2700989144_000001|10  121 105   30 |             108    94 |   0/105
IMGPR_plasmid_2551306362_000001|6   105 126   15 |              63   171 |   0/126
IMGPR_plasmid_2903362858_000001|1    78 187   14 |             107    88 |  12/187
```

All five sit beside named CDS as often as beside dark ones — consistent with real genes in real
operonic context, and inconsistent with isolated spurious calls in gene deserts, but adjacency is not
overlap.

**The minimal test that would settle it** (design only; ~800 sequences, hours not days):

1. **Recover coordinates.** Re-run `extract_plasann_proteins.py` with the C6 patch to emit
   `(plasmid_id, cds_index, start, end, strand)`, and pull the CDS nucleotide sequences for the five
   families' ~788 members.
2. **Overlap/antisense screen (the cheapest discriminator).** For each member, compute overlap with
   every other CDS on the same plasmid and the relative strand. An artefact family is *systematically*
   antisense to, or nested within, a named gene. A pass here removes the leading hypothesis in one
   pass over a table.
3. **Selection test (the decisive one).** Codon-aware alignment per family (MACSE or PRANK-codon),
   then dN/dS by counting or PAML M0. A protein-coding family under purifying selection gives
   dN/dS ≲ 0.3 with high confidence at n = 100+; a conserved *non-coding* element carrying an
   accidental ORF gives dN/dS ≈ 1 with the DNA still conserved. This distinguishes the two hypotheses
   directly and needs only ~20 members per family, so all five fit comfortably.
4. **Independent re-calling.** Re-call the ~150 carrier plasmids per family with Prodigal (`-p
   single`, plasmid-trained) and Bakta, and score whether the same frame is called. Cross-caller
   concordance plus a Shine–Dalgarno/RBS motif score is standard evidence.
5. **Structure (confirmatory only).** ESMFold the five representatives; pLDDT > 70 with a compact
   fold is strong positive evidence. The verification pack §9 notes FlashFold target selection already
   exists — fold these five first.

Steps 2 and 3 alone settle the question. Step 3 is the one a reviewer will demand.

---

### REPRODUCIBILITY GAPS

I checked two "Reproduce" blocks end to end for missing prerequisites. **Both would run** — every
input file, script, database and binary they name exists on disk (585/595 + 477 PlasAnn annot shards,
`mob_full.tsv.gz`, `working_set.tsv`, `dark_orfs.faa`, `Pfam-A.hmm`, `mmseqs`, `/scratch`). The gaps
are in the environment layer, not the data layer.

**[MAJOR] R1 — the two interpreters that produced the entire dark-ORF, Pfam and statistics layer are
undocumented personal envs outside the repo**

*Anchor:* `PROJECT_OVERVIEW.md` §6 documents exactly three environments — `rgi_env`, `segmantx_env`,
`stats_env` — all under `envs/.micromamba_root/`, all present, all with specs in
`envs/live_env_specs/`. But:

```
$ grep -rn "amujkic/.conda/envs" scripts/
 -> panaroo   in recluster_dark_orfs.sh, cluster_dark_orfs.sh, pfam_task.sh,
               recluster_pfam_array.sbatch, recluster_union.sbatch, fam_bigplasmid_search.sbatch,
               audit_cluster_compliance.py:24 (MMSEQS hardcoded)
 -> genesis   in plasann_task.sh:9, fam_bigplasmid_search.sbatch, and the genesis_nb kernel
```

Neither `genesis` nor `panaroo` has an exported spec anywhere in the repository, and both live in
`/gorilla/home/amujkic/.conda/` — a personal home directory §6 explicitly says is quota-constrained
and was the reason for going project-local in the first place. Everything in
`reports/dark_orf_clustering.md`, `reports/pfam_dark_validation.md`,
`reports/clustering_threshold_sensitivity.md` and `reports/small_cryptic_methodology.md` was produced
by these two.

*Fix:* `micromamba env export` both into `envs/live_env_specs/{genesis,panaroo}.yml` and add them to
§6's live table. Highest-value reproducibility action in the project.

**[MAJOR] R2 — `plasann_env` is referenced by live scripts but no longer exists**

*Anchor:* `scripts/plasann_task.sh:8` (`ENVBIN=/gorilla/home/amujkic/.conda/envs/plasann_env/bin`) and
`scripts/plasmidfinder_run.sbatch`. §6 lists it as decommissioned with only
`envs/plasann_env.environment.yml` retained — but the scripts still hardcode the dead absolute path,
so they fail with a confusing `No such file` rather than a clear "rebuild this env first".
*Fix:* add a `[ -x "$ENVBIN/python" ] || { echo "rebuild from envs/plasann_env.environment.yml"; exit 1; }`
guard at the top of each.

**[MINOR] R3 — the small-cryptic reproduce block does not pin its kernel, and `jupyter` resolves elsewhere**

*Anchor:* `reports/small_cryptic_methodology.md` Reproduce block. The prose says "execute with the
**genesis_nb** kernel"; the command is `jupyter nbconvert --to notebook --execute --inplace ...` with
no `--ExecutePreprocessor.kernel_name`. The notebook metadata does name `genesis_nb`
(`{'display_name': 'genesis_nb', 'name': 'genesis_nb'}`), so it works today — but `which jupyter` in a
bare shell gives `/sw/apps/conda/latest/rackham_stage/bin/jupyter`, the site-wide install, not
anything the project controls. *Fix:* add `--ExecutePreprocessor.kernel_name=genesis_nb` and give the
absolute `jupyter` path.

**[MINOR] R4 — `audit_cluster_compliance.py` hardcodes two machine-specific paths**

*Anchor:* `:24` `MMSEQS = "/gorilla/home/amujkic/.conda/envs/panaroo/bin/mmseqs"` and `:60`
`tempfile.TemporaryDirectory(dir="/scratch")`. Neither is a CLI flag or an env var. *Fix:*
`os.environ.get("MMSEQS", "mmseqs")` and `dir=os.environ.get("TMPDIR", "/scratch")`.

**[MINOR] R5 — `clustering_threshold_sensitivity.md` names no interpreter**

It is the only recent report with a full "Scripts that produced this report" table and no environment
statement — a small gap in an otherwise exemplary provenance habit. (Coverage across `reports/`:
`amr_onehealth` → stats_env; `card_amr` → rgi_env; `dark_orf_clustering`, `pfam_dark_validation` →
genesis + panaroo; `small_cryptic` → genesis_nb; `typing` → micromamba; the rest silent.)

**[MINOR] R6 — RNG hygiene in the notebook builder**

`scripts/build_small_cryptic_notebook.py:363` creates one module-level
`np.random.default_rng(20260826)` consumed sequentially by every `rarefied_breadth` call, so the
K-ladder values depend on top-to-bottom execution order — reproducible only for a full clean run, not
for re-executing a single cell. The boxplot jitter at `:417-419` uses the unseeded legacy global
`np.random` (cosmetic only). Also `reps` differs across analyses (200 pooled / 60 ladder / 100
isolate) without comment. *Fix:* instantiate a fresh seeded generator inside `rarefied_breadth`.

---

## Scores

| dimension | score | justification |
|---|---:|---|
| **Statistical rigour** | **66** | Genuine design thinking: rarefaction against sampling depth, a real positive control, a propagation check, cluster-robust SEs, provenance stratification, a threshold sweep, an artefact control. Losses: pseudo-categories counted as data (D2), a headline that contradicts the project's own K-ladder (D3), an effect decomposition on the wrong scale for an unidentified estimand (D5), a sign-inverted effect size (D4), no interval estimates anywhere in §3, and no sensitivity on the single most consequential threshold in the resource layer (D1). |
| **Code correctness** | **72** | Structurally sound where it counts: strictly 1:1 joins with no fan-out (C11), the CDS-index fix real and verified, dedup provably safe on the actual data (C7), RGI coverage complete so the zero-fill is benign (C9), clean regex parsing of all 481 database ids. Losses: a recurring "count labels, not entities" family across three columns (C1, C2, C3), a low-recall classifier driving a published percentage (C4), and safety properties that hold by luck rather than by assertion (C7–C9, C11–C12). |
| **Reproducibility** | **65** | Well above field norm on data: every report names the scripts that made it, seeds are set, all intermediates are on disk, both Reproduce blocks have their prerequisites, and the sweep was deliberately written not to overwrite what three notebooks depend on. Losses are concentrated in one place: the two interpreters behind the entire dark-ORF/Pfam/statistics layer are unexported personal envs (R1), one referenced env is gone (R2), and two machine-specific paths are hardcoded (R4). |
| **Methods-prose fidelity** | **58** | The weakest axis and the one most likely to be caught in review. "Identical criteria to CGE PlasmidFinder" when CGE's default is 95% (D1). "Families collapse allelic variants" one paragraph before counting alleles (C1). "47.3% inside PlasAnn's categories" when a name-matching classifier put ≥54.8% there (C4). "Near-parity" against a monotone ratio decline (D3). NAMED and DARKREP tabled as one comparison when they are different units (D6). "hab_sub present for a minority" when it is 100% (D11). Each is a sentence the code does not support. |

---

## Genuine methodological strengths

These are real, and several are things I tried to break and could not.

1. **The master build cannot multiply rows.** Six dict-keyed lookups over verified-unique keys, one
   output row per `analysis_table` row, `extrasaction="ignore"`. A `merge` chain would have been
   fragile; this is not. I checked all six inputs: 1:1 across the board.
2. **The Pfam propagation check is correctly designed and correctly denominated** — contrary to the
   suspicion I was asked to test. `_valid_targets` reads all 3,678 submitted VALID queries from the
   chunk FASTAs, not the 1,140 hit rows in `pfam_per_seq.tsv`. And the result is strong on the
   statistics that actually matter, not just on raw concordance:
   ```
   2x2 rep_hit x member_hit:  TT 1070  TF 73  FT 70  FF 2465
   overall concordance 96.11%   (all-negative baseline would be 69.0%)
   P(member hits | rep hits)     = 93.6%
   P(member no-hit | rep no-hit) = 97.2%
   Cohen kappa = 0.909
   ```
   Both conditionals exceed 93%; κ = 0.909. Family-level propagation is justified. (One real gap: the
   sample is drawn only from families of 10–50 members, and families >50 hold **34.5%** of all dark
   ORFs and are untested — worth a second draw.)
3. **The mash-neighbourhood lineage proxy is not circular in the way I expected.** `mob_cluster` has
   7,054 levels, median size 6, only 0.4% of plasmids in singleton clusters; among dark families with
   ≥2 members the lineages/plasmids ratio is 0.48 and Spearman ρ(plasmids, lineages) = 0.51.
   "Spanning ≥2 lineages" is not a relabelled plasmid count. (The `mob_host_range` circularity in D8
   is a different and real problem.)
4. **The rarefaction ≥K filter is not the selection trap it looks like.** It retains 86.5% / 87.6% of
   the two arms for habitats and 72.0% / 65.6% for countries — comparably, not differentially. And
   only 3% of conjugative clusters also appear in the cryptic arm, so the Mann–Whitney independence
   assumption is essentially satisfied.
5. **The One Health regression set is clean.** Zero Simulated-artifact or Lab-artifact plasmids leak
   into the 104,169 (the locked exclusion is honoured implicitly via `hab_sub`), and `mob_mobility` is
   missing for exactly **one** plasmid — so the conditioning at
   `reassess_review_findings.py:189` introduces no selection bias at all. I expected this to be a
   finding; it is not.
6. **`--cluster-reassign` was a correct and well-evidenced decision**, and I confirmed it independently
   under ORF weighting: 88.6% → 97.6% compliance. The literature table in
   `clustering_threshold_sensitivity.md` §2 is honest about where the project is stricter than
   published practice and where it is not, and the sweep was deliberately routed to
   `recluster/` so the three dependent notebooks stay reproducible. That is unusually disciplined.
7. **Verification pack §7.3 resolves in the project's favour** — I found the missing expression. The
   pack could not reconstruct "95,442 (46%)" because it used PlasAnn; `typing_methodology.md:93`
   states the set as PlasmidFinder ∪ MOB-suite ∪ **PLSDB**:
   ```
   pf | mob                   95331 (45.8%)
   pf | mob | plsdb           95442 (45.8%)   <- exact match to the reported figure
   pf | mob | plasann         97668 (46.9%)   <- what the pack tried
   ```
   The number is correct as published; the ⚠️ should be withdrawn. (The `(46%)` rounds 45.8%, which
   is worth tightening to 45.8%.)
8. **The self-critical documentation culture is the project's real asset.** Pitfalls are recorded
   inline for reuse (the two-shard-series trap), limitations sections name the leading alternative
   explanation rather than hedging, and `reassessment_round1.md` marks its own findings
   UPHELD/REFUTED/QUANTIFIED including one that clears an objection against interest. Most of what I
   found, this project is already equipped to fix.

---

## The three highest-value methodological fixes, ranked

**1. Re-derive every PlasmidFinder headline at the tool's own default, and define "multireplicon" on
collapsed families.** *(D1 + C1 — one line of code, two reports.)*
This is first because it is the project's founding question and the number is unstable by 65%:
multireplicon is 24,583 at 80%/alleles, 22,422 at 80%/families, and **8,574 at 95%/alleles**. Change
`pf_n_inc = len(fams)` in `harvest_typing.py:80`, publish the 95% column beside the 80% one in
`typing_methodology.md` §4, replace "identical criteria to CGE PlasmidFinder" with a PLSDB citation,
and grep `PROJECT_OVERVIEW.md` §4 and the deck for the bare 24,583.

**2. Rebuild small-cryptic §3 and retire "near-parity".** *(D2 + D3 + D4 + D11 — one notebook.)*
Map `Unlabelled` and `Geography only (no habitat)` to NaN before rarefaction (this alone moves the
habitat ratio 0.917 → 0.853 and p 4.9e-4 → 6.8e-5); fix the rank-biserial sign; bootstrap a
cluster-resampled CI on the ratio; lead with the isolate-only stratum since §7 already establishes
provenance as the confounder; and replace "near-parity" with the K-monotone statement the ladder
actually supports. The claim gets weaker and much harder to attack.

**3. Fix the two places where a classifier's recall silently became a scientific percentage, and the
one place a coefficient change became a causal share.** *(C4 + D5.)*
Swap `pfam_verdict.py`'s name-regex for a Pfam clan/accession map and re-publish §4.3 (47.3% →
≥54.8% inside PlasAnn's scope — `RelE`, `TrfA`, `RepC`, `TrwC`, `Rop` are not novel proteins). In the
same pass, delete "about 61% of the raw conjugative effect is plasmid size" from
`reassessment_round1.md:117` and `:187`: it is 37.6% on the log-odds scale, it is confounded by
non-collapsibility, and size is partly downstream of conjugativity. Report the model sequence, or do
g-computation properly.

*Honourable mention (cheapest of all):* export `genesis` and `panaroo` to `envs/live_env_specs/`
(R1). Ten minutes of work standing between this project and an outsider reproducing its central
result.
