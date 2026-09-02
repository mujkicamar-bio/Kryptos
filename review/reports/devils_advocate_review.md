# Devil's Advocate review — the dark plasmidome

*Reviewer 5 of 5, working independently. Written 2026-09-02. Read-only: nothing in the project was
modified. All recomputes were run from the tables on disk with
`/gorilla/home/amujkic/.conda/envs/genesis/bin/python` (pandas 3.0.2) and
`/gorilla/home/amujkic/.conda/envs/panaroo/bin/hmmsearch` (HMMER 3.4); scratch scripts live under
the session scratchpad, not in the repo.*

**Scope note.** I took `review/VERIFICATION_PACK.md` at its word that the arithmetic reproduces, and
attacked premises instead. I independently re-derived the headline cross-tab (57.4 / 28.5 / 10.9 /
3.2), the 3,343 Pfam-negative families, the 71.1% typing blind spot, and the 25.9% vs 2.8%
host-range inversion as a side effect of attacking them; all four reproduce exactly. Pack finding
§7.3 was retracted by the coordinator before I filed and is not used here.

---

## STRONGEST COUNTER-ARGUMENT

The dark plasmidome result is a measurement of sampling depth and annotation-tool coverage that has
been dressed as a measurement of biology, and the project's own designated control is the one
control it never applied to it.

Start with C4, because C1 rests on it. "The recurrence is not clonal redundancy" is established by a
criterion — *does a family contain members from ≥2 MOB clusters?* — that is close to the weakest
test that could have been written. It is satisfied by a family of 2,487 proteins in which 2,486 sit
in one clone and one does not. Apply instead the control this project uses everywhere else, one
plasmid per MOB cluster, and recurrence falls from 84.7% to 39.2%, the substantial Pfam-negative
family set falls from 3,343 to ~40, and the 93-family shortlist falls to ~30 families that still
reach ten carriers. That is not a sample-size artefact: a random subsample matched at the same ORF
count retains 61.3% recurrence, twenty-two points higher. The project applied MOB dereplication in
`cryptic_plasmids.ipynb` and `defense_systems.ipynb`, where it cut results it was willing to lose.
It appears in no dark-family notebook.

Then note that recurrence is not a constant of the object. Rarefying carriers gives 39.2% at 1,000
plasmids, 68.6% at 10,000, 84.7% at 63,993. "84.7% recurrent, therefore World B not World A" is
therefore not a discriminating test — it is a reading of where the sampling happened to stop, and no
null model for World A was ever computed.

Finally, C1's novelty is bounded by one database. Pfam-A at `--cut_ga` is a *family membership*
criterion, not a homology detector, and no DIAMOND or BLAST search against UniRef/NR was ever run
anywhere in this project. Re-searching 2,000 of the "genuine residue" families with no GA threshold,
against a composition-shuffled null, recovers real Pfam similarity for 27 percentage points of them
at E ≤ 1e-3 — and the domains recovered are HTH, Rep_1, Relaxase, Mob_Pre, NikA-like. The same
backbone the project already conceded once.

---

## FINDINGS

Ordered by severity. Every attack I tested is reported, including the ones that failed — three of my
most promising lines failed outright and I say so.

---

### [CRITICAL] F1 — The project's own clonal-redundancy control was never applied to the recurrence claim, and it takes most of the claim with it
**Attacks:** C4 (primary), C1 (consequentially)
**Verdict: SURVIVES.**

C4 is certified by one statistic: 43.8% of ≥2-member families, and 79.2% of ≥10-member families,
span ≥2 MOB clusters (`reports/dark_orf_clustering.md` R2). The project elsewhere uses a much
stronger control for exactly this failure mode — dereplication to one plasmid per MOB cluster, which
halved every AMR prevalence in Act I and is described in the deck as "the same MOB-cluster control
used throughout Act II". It is not used in Act II.

I applied it. Twenty random draws of one plasmid per MOB cluster, families as published:

```
scenario                              ORFs    fams   rec%  f>=10  pfneg   pfneg%
pooled (published)                  378552   92752   84.7   5623   3343     30.6
1 plasmid / MOB cluster (mean 20)    32763   23514   39.2    159     40      2.1
random 4004 plasmids (mean 20)       23680   13074   57.7    235    114      9.1
cap 5 per MOB cluster (mean 5)       91092   40989   69.8   1096    481     10.2
cap 10 per MOB cluster (mean 5)     125203   46938   75.7   2014   1040     16.8
cap 50 per MOB cluster (mean 5)     206608   56043   83.2   3537   2046     28.7
```

The matched control is the load-bearing row. Dereplication is not merely subsampling:

```
1-per-MOB-cluster: 32770 ORFs, recurrence 39.1%
   random  5500 plasmids -> 32553 ORFs, recurrence 61.3%
```

At the same ORF count, removing clonal structure costs 22 percentage points of recurrence. The
93-family core does not survive either:

```
of the 93, families still on >=10 distinct plasmids after dereplication: mean 30.0 (range 27-33)
median dereplicated plasmid count per core family: 7 (they were >=100 before)
```

A family selected for being "on ≥100 distinct plasmids" is, after dereplication, on a median of
seven. The ≥100-plasmid criterion was measuring deposits.

Confirmed that the control is absent from the dark work:

```
$ grep -rn "drop_duplicates.*mob_cluster" notebooks_my/*.ipynb
notebooks_my/defense_systems.ipynb:345:  DER = m[m.carrier].sample(...).drop_duplicates("mob_cluster")
notebooks_my/cryptic_plasmids.ipynb:85:   derep = conj.sample(...).drop_duplicates("mob_cluster")
```

Nothing in `dark_plasmidome`, `dark_families`, `dark_families_90pct`, `widespread_dark_orfs` or
`zoom_two_dark_orfs`.

**Fair counter, which I accept in part.** One-per-cluster is harsh: 4,174 clusters for 71,414
plasmids, and AA379 alone absorbs 26% of the subset (see F5), so it collapses genuine diversity.
The cap-K rows are the honest middle: at a cap of 10 per cluster, recurrence is 75.7% and 1,040
Pfam-negative substantial families remain. **The defensible statement is not "84.7%" but "70–76%
after capping clonal over-representation, and 39% if you dereplicate completely" — and the shortlist
is one order of magnitude smaller than published under any of them.**

---

### [CRITICAL] F2 — Recurrence is a function of sampling depth, and no null model for the rejected hypothesis was ever computed
**Attacks:** C1, C4
**Verdict: SURVIVES.**

The framing device of the whole project — World A (idiosyncratic junk) vs World B (annotation gap),
"these make opposite predictions about the *shape* of the dark proteome" — is adjudicated by a
statistic that is not a property of the dark proteome:

```
Full recurrence rarefaction (random plasmids):
     1000 plasmids ->    5836 ORFs, recurrence 39.2%
     4000 plasmids ->   23405 ORFs, recurrence 57.4%
    10000 plasmids ->   59433 ORFs, recurrence 68.6%
    20000 plasmids ->  118433 ORFs, recurrence 75.3%
    40000 plasmids ->  236863 ORFs, recurrence 81.3%
    63993 plasmids ->  378552 ORFs, recurrence 84.7%
```

Recurrence moves 45 points across the sampling range with the biology held fixed. A study of 10,000
carriers would have reported 68.6% and drawn the same qualitative conclusion; a study of 1,000 would
have reported 39.2% and drawn the opposite one. Nothing anywhere in the project states what World A
predicts at n = 63,993. Junk inherited by descent through a clonally over-sampled corpus produces
high recurrence too — that is precisely what F1 measures. **The two "opposite predictions" are not
opposite, and the test that is claimed to separate them does not.**

This is the deepest problem in Act II, and it is entirely fixable: the sensitivity analysis that was
done for the clustering threshold (thoroughly, to the project's credit) was never done for sampling
depth or for a null.

---

### [MAJOR] F3 — "Spans ≥2 lineages" is a near-vacuous criterion; under any median-standard version the majority of the recurrence is lineage-dominated
**Attacks:** C4
**Verdict: SURVIVES.**

For the 3,343 substantial Pfam-negative families — the project's headline deliverable — I computed
the share of members sitting in the single dominant MOB cluster:

```
dominant-lineage share: median 0.73  q25 0.50  q75 1.00
   families with >=50% of members in ONE mob_cluster: 2519 (75.4%) holding 65.6% of the 115,811 ORFs
   families with >=80% of members in ONE mob_cluster: 1497 (44.8%) holding 35.3%
   families with >=90% of members in ONE mob_cluster: 1219 (36.5%) holding 26.6%
   families with 100% of members in ONE mob_cluster:   905 (27.1%) holding 16.5%

the criterion actually used ("spans >=2 lineages") is met by 72.9%
a stricter one ("no lineage holds >50% of members") is met by 24.6%
```

And for the whole dark proteome, against R2's headline 63.9%:

```
ORFs in families spanning >=2 lineages:                63.9%
ORFs in families where no single lineage holds >50%:   26.8%
ORFs in families >=90% confined to one mob_cluster:    40.9%
```

Twenty-seven per cent of the "lineage-spanning" deliverable is entirely confined to one cluster —
those families are in the 3,343 because the ≥2-lineage filter is applied to family membership, not
to the 3,343 themselves, but they are nonetheless what the reader is invited to see as
lineage-crossing. And 40.9% of the dark proteome sits in families ≥90% confined to a single cluster.

The project *knew* this pattern. `dark_families.ipynb` records that the top family's "88 MOB
clusters are top-heavy — the largest three hold 70% of carriers". It flagged the concentration for
one family and did not compute it for the other 3,342. **The right summary statistic was one line of
pandas away and it changes 72.9% into 24.6%.**

---

### [MAJOR] F4 — The "genuine residue" is partly a gathering-threshold artefact, and Pfam-A is the only external database ever searched
**Attacks:** C1
**Verdict: SURVIVES, with an important correction to my own first attempt.**

`--cut_ga` is well chosen for what it is — Pfam's curated *family membership* criterion, with no knob
to tune — and the report is right that it makes the result untunable. But "is a curated member of a
Pfam family" and "has no detectable homology to anything in Pfam" are different claims, and the
project reports the first while asserting the second ("invisible to both plasmid annotation and
Pfam").

I re-searched 2,000 randomly sampled representatives of the 3,343 Pfam-negative substantial families
against the same Pfam-A 38.2 with no GA threshold, plus a composition-shuffled negative control of
the same 2,000 sequences (shuffling preserves amino-acid composition and length, destroys homology):

```
$ hmmsearch --noali --cpu 24 -E 1 --domtblout ... data/refs/pfam/Pfam-A.hmm novel_sample.faa
$ hmmsearch --noali --cpu 24 -E 1 --domtblout ... data/refs/pfam/Pfam-A.hmm novel_shuf.faa

E cutoff      real hit rate  shuffled(FPR)       excess
1.7e-06                4.0%           0.0%         4.0 pp     <- Bonferroni over 30,134 models
1e-05                  8.1%           0.0%         8.1 pp
0.0001                15.2%           0.4%        14.8 pp
0.001                 30.6%           3.4%        27.2 pp
0.01                  56.0%          22.3%        33.8 pp
```

**My first pass reported the 30.6% figure without the control and that would have been wrong** —
hmmsearch's sequence E-value does not carry the multiple-testing burden of 30,134 profiles, so a raw
E ≤ 1e-3 is badly contaminated. The shuffled control is what makes the number interpretable, and it
says the contamination at E ≤ 1e-3 is 3.4 points, leaving 27.2 points of genuine sub-GA homology.
Under strict Bonferroni it is 4 points. **The honest range is that between 4% and 27% of the
"genuine residue" families do have real, detectable Pfam-A homology that the gathering threshold
suppresses** — and the recovered domains are not novelty:

```
best-hit families among those rescued at E<=1e-3:
HTH_17 18 · NikA-like 14 · HTH_23 13 · RepB-RCR_reg 11 · RHH_1 10 · HTH_36 8 · Rep_1 8 ·
Relaxase 8 · DUF3847 7 · PhdYeFM_antitox 7 · Mob_Pre 7 · HTH_11 7
```

Exactly the backbone the project already identified as the artefact half. The artefact half is
bigger than 28.5%.

Worse, and this is the part I would press hardest at review: **Pfam-A is the only external database
this project has ever searched.**

```
$ grep -ril "uniref\|diamond\|blastp\|foldseek" scripts/ notebooks_my/ reports/*.md
# hits are: UniRef cited as a *clustering threshold* precedent; diamond inside the RGI env;
# foldseek named as future work. No sequence search against UniRef/NR/UniProt exists.
```

A protein can be absent from Pfam-A's 30,134 families and still have a thousand near-identical
homologs in UniProtKB. A DIAMOND run against UniRef90 is hours of compute and it is the single
cheapest test that could falsify C1. It has not been run, and until it is, "no tool can name them"
means "two tools, one of them a plasmid-specific annotator and the other a domain database, cannot
name them".

---

### [MAJOR] F5 — `mob_cluster` under-splits catastrophically at the small end; AA379 is a junk drawer, not a lineage
**Attacks:** C4, C5, and every "lineage" number in the project
**Verdict: SURVIVES as a semantics problem — but see F6, my over-splitting attack FAILED.**

```
== AA379 (the largest cluster in the subset)  n=18749 of 71,414 (26%)
size: median 4327 IQR 3032-7345  range 236-19992
GC:   median 45.0 IQR 38.7-54.5
mobility: non-mobilizable 18602 · mobilizable 147
AA379: distinct exact proteomes 15954 over 18749 plasmids (85.1% unique)
```

A single "primary cluster" containing 18,749 plasmids spanning 236 bp to 19,992 bp, GC from 38.7% to
54.5%, and 15,954 distinct protein complements is not a lineage. It is the bin into which
mash-distance clustering puts small plasmids it cannot resolve. Granularity confirms the pattern —
the smallest plasmids get the *fewest* clusters per capita:

```
== mob_cluster granularity by plasmid size (analysis set)
              n  clusters  plasmids_per_cluster  clusters_per_1000_plasmids
<5kb      40776      2286                 17.84                       56.06
5-10kb    29467      2466                 11.95                       83.69
10-20kb   12018      2416                  4.97                      201.03
20-50kb   18322      3031                  6.04                      165.43
>100kb    24604      2463                  9.99                      100.11
```

Consequences run in **both** directions and the project acknowledges neither. Under-splitting makes
"confined to one lineage" over-counted, which flatters nothing — it means the 56.2% of families
called clonal may not be. But it equally means 28.7% of the 3,343 deliverable families draw ≥50% of
their members from AA379 alone (15.0% draw ≥90%), and that "spans N lineages" has no stable
biological referent when one of the N is a bin of 18,749 unrelated replicons.

The limitation is stated in both reports ("`mob_cluster` is a proxy, not a phylogeny... AA379 alone
holds 19,427") and the deck asks the room about it directly (question 2). That is honest. But a
limitation acknowledged is not a limitation controlled, and every lineage-based control in Act II
runs on this column.

---

### [MAJOR] F6 — The typing blind spot is substantially a gene-count effect and partly definitional
**Attacks:** C2
**Verdict: PARTIALLY SURVIVES — the 9.1× gap is roughly 3.8×, not zero.**

Reproduced first: small-cryptic 71.1% untyped (n=47,031) vs large-conjugative 7.8% (n=25,631).

A 4 kb plasmid with 5 CDS has fewer chances to carry a detectable replicon than a 100 kb plasmid
with 120. Stratifying by PlasAnn CDS count:

```
grp     small-cryptic_untyped%  large-conjugative_untyped%   sc_n    lc_n
10-14                     67.9                        19.4   5302      36
15-24                     53.6                        27.0    735     274
25-49                     13.2                        16.2     38    2448
```

At 25–49 CDS the gap **reverses**. Restricting to the region where the two groups genuinely overlap:

```
== overlap-region-only comparison (10-49 CDS)
                      n  untyped
large-conjugative  2758     17.3
small-cryptic      6075     65.8
```

**3.8×, not 9.1×.** I also computed a full CDS-standardised reweighting (large-conjugative rises to
69.8%), but I am **not** relying on it and neither should anyone else: the two CDS distributions
barely overlap, so that number is extrapolated off 2–274 observations per bin. The overlap-region
comparison is the honest one, and it says roughly 60% of the headline gap is gene-count composition.

Separately, the definition does some of the work. "Cryptic" excludes plasmids carrying PlasAnn
conjugation genes (`mobA`, `mobC` are in that category), and those are exactly the plasmids most
likely to carry a detectable replicon:

```
small <10kb, payload-free WITHOUT conjugation filter: n=62416 untyped 61.7%
small <10kb, payload-free WITH    conjugation filter: n=47031 untyped 71.1%
the 15,385 removed by the conjugation filter alone:            untyped 32.9%
```

Nine points of the 71.1% headline are manufactured by filtering out well-annotated plasmids and then
reporting that the remainder is poorly annotated. **C2 is real but the magnitude is inflated roughly
2.4-fold by gene count and a further ~9 points by the payload definition.** The blind spot survives;
"6.8× gap in isolate-derived data alone" does not.

---

### [MAJOR] F7 — "Near-parity" (C3) is a non-significant underpowered subgroup reported as equivalence, and the pooled evidence points the other way
**Attacks:** C3
**Verdict: SURVIVES.**

From `smallcryptic_dispersal_Kladder.tsv` and `smallcryptic_dispersal_isolate_only.tsv`:

```
K   metric      median_cryptic  median_conj  ratio   p
5   habitats              3.07         3.30  0.929   1.7e-05
10  habitats              4.57         5.00  0.913   0.00037
20  habitats              6.58         7.35  0.896   0.00053
40  habitats              7.88         9.86  0.800   0.00017
5   countries             3.00         3.42  0.878   1e-09
40  countries             9.68        12.22  0.793   8.6e-06

isolate-only:  habitats  4.96 vs 5.27  ratio 0.941  p=0.0597  (235 vs 305 clusters)
               countries 5.66 vs 5.38  ratio 1.052  p=0.113   (181 vs 269 clusters)
```

Three problems, in ascending order.

1. The pooled result is not parity. It is a **consistent deficit**, significant at every K, in which
   non-mobilizable small plasmids reach fewer habitats and fewer countries per lineage.
2. The deficit **grows with K** — 0.93 → 0.80 for habitats, 0.88 → 0.79 for countries. The report
   quotes the K = 10 numbers (0.92, 0.83) as the headline and says only that "the ordering is stable
   across K". The ordering is stable; the magnitude doubles. At K = 40 it is a 20% deficit.
3. "Near-parity" is then asserted from the isolate-only subgroup, where n drops ~2.6× and p rises to
   0.06 and 0.11. **That is a failure to reject, in an underpowered subgroup, converted into a claim
   of equivalence.** No equivalence test (TOST or otherwise) was performed, and a ratio of 0.941 with
   p = 0.06 is fully consistent with a real 6% deficit the subgroup cannot see.

To the project's credit the isolate-only table is computed, reported, and the p-values are printed
in the report. The error is purely interpretive — but the interpretation is the claim, and the claim
is called a "paradox".

---

### [MODERATE] F8 — The host-range artefact control is non-monotone, and the report quotes only its endpoints
**Attacks:** C5
**Verdict: PARTIALLY SURVIVES — a framing failure; the core inversion is real.**

C5 reproduces exactly: multi-phyla 25.9% (small-cryptic) vs 2.8% (large-conjugative), 9.3×. The
project's artefact control conditions on mash neighbourhood distance and the report states that "for
small-cryptic the rate *falls* with increasing distance (17.4% in the >0.1 bin)". The project's own
stored table says otherwise:

```
$ cat data/plasmidscope_primary/smallcryptic_hostrange_vs_mashdist.tsv
small-cryptic  identical (0)  5143  25.5
small-cryptic  <=0.01         2702  42.3     <-- not quoted anywhere
small-cryptic  0.01-0.05      3347  24.8
small-cryptic  0.05-0.1       2979  20.1
small-cryptic  >0.1           2421  17.4
```

The rate rises to 42.3% before it falls. Quoting 25.5 and 17.4 and describing a fall is selecting
the two endpoints of a peaked curve. I also note the "identical (0)" bin is `mnd ≤ 0.001`, not zero;
at literally zero I get 16.1% (n = 1,755) vs 2.6% for large-conjugative.

**Why this still fails as a kill:** even at exact mash distance 0 the gap is 6.2×, and the
multi-phylum calls come from 226 distinct clusters with the top five holding only 37.2%. The
inversion is not an artefact of noisy neighbourhoods. But the control as written is stronger in the
report than it is in the data.

---

### [MODERATE] F9 — Redundancy the 100%-identity dedup did not catch
**Attacks:** C1, C4
**Verdict: PARTIALLY SURVIVES — real, but smaller than I expected.**

`PROJECT_OVERVIEW.md` §3 states flatly that PlasmidScope `ALL` is deduplicated at MMseqs2 100%
identity and coverage, so "there are no duplicate plasmids across source databases". Nucleotide
identity at 100% does not catch circular rotation, single-base differences, or independent deposits
of the same replicon. Hashing the sorted multiset of every CDS translation per plasmid:

```
distinct proteome hashes: 51872 over 65805 plasmids
plasmids sharing their exact proteome with >=1 other: 19645 (29.9%)
largest identical-proteome groups: 180, 166, 160, 86, 82, 81, 76, 75, 75, 71
collapsing to one plasmid per identical proteome leaves 51872 (78.8%)
```

Nearly a third of the corpus is a protein-level duplicate of something else in it. Recomputing the
headlines on the dereplicated set:

```
pooled (published)           ORFs 378552 | recurrent 84.7% | Pfam-neg>=10 3343 holding 30.6%
exact-proteome dereplicated  ORFs 318124 | recurrent 81.0% | Pfam-neg>=10 2654 holding 26.2%
```

**The attack lands but does not kill: 84.7 → 81.0 and 3,343 → 2,654.** The claim in
`PROJECT_OVERVIEW.md` that there are no duplicates should be corrected to "no *sequence-identical*
duplicates", but redundancy at this level is not what generates the recurrence result. F1 is where
the damage is.

---

### [MODERATE] F10 — The "median 138 aa" defence is about the pooled set; the residue being claimed as novel is much shorter
**Attacks:** C1
**Verdict: SURVIVES as framing — see F11, the spuriousness attack itself FAILED.**

The deck kills the short-ORF objection with one line: "the length distribution kills the cheapest
objection — these are not spurious short ORF calls; the median is 138 aa." That median is of the
whole dark proteome, including the 39.4% Pfam can name. Split it:

```
                       n       median   q25   q75   %<50aa  %<100aa
ALL dark ORFs      378552        138      86   249     3.9     34.6
Pfam-POSITIVE      149160        214     115   339     0.3     18.1
Pfam-NEGATIVE      229392        108      74   178     6.3     45.3
NAMED (PlasAnn)     78397        185     100   312     1.1     25.0
Mann-Whitney Pfam+ vs Pfam- ORF length: p = 0
```

The residue is median 108 aa with 45.3% under 100 aa. Small proteins are simultaneously (a) the class
where *ab initio* gene calling is least reliable, (b) the class Pfam covers worst, and (c) the class
where an HMM has least signal to work with. So the headline number is defended with a statistic
drawn from the half of the data that is not in dispute. **The defence should be restated on the
Pfam-negative subset, where it is weaker but still non-trivial.**

---

### [MINOR — attack FAILED] F11 — Spurious ORF calls on non-coding DNA
**Attacks:** C1
**Verdict: FAILS. I tested this hard and the data say the dark ORFs are real proteins.**

I expected this to be my best line and it collapsed under three independent tests.

*Amino-acid composition.* If the Pfam-negative residue were ORF-calling noise on non-coding DNA at
~44% GC, its composition would approach the random-codon expectation. It does not — it is nearly
indistinguishable from PlasAnn-named plasmid proteins:

```
L1 distance  NAMED vs darkPfam+ = 0.083
L1 distance  NAMED vs darkPfam- = 0.099
L1 distance  darkPfam+ vs darkPfam- = 0.028
L1 distance  NAMED vs random-GC44 = 0.331
L1 distance  darkPfam- vs random-GC44 = 0.303
```

Random ORFs are Leu/Ser/Arg-rich and Glu/Lys-poor. The dark set shows the opposite: E at 7.5% and K
at 7.4% against random expectations of 3.3% and 4.2%; Cys at 1.15% against 3.26%; Ser at 6.3%
against 9.8%. That is a translated-and-selected proteome.

*Coding density.* Over-calling inflates density. It is **low**, not high:

```
coding density (sum CDS nt / plasmid bp): median 0.687  q25 0.576  q75 0.766 ; >1.0: 0.0%
(bacterial replicons typically 0.85-0.90)
```

PlasAnn is under-calling on these plasmids, not over-calling. If anything the dark proteome is
larger than reported.

*Overlap.* Spurious calling produces mutually overlapping ORFs on opposite strands:

```
dark ORFs with >60bp overlap of an opposite-strand CDS: 130 / 17353 = 0.7%
```

**Attack 2 fails.** The dark ORFs are real coding sequence. The only survivor from this line is the
framing point in F10.

---

### [MINOR — attack FAILED] F12 — "Most of this is metagenomic assembly artefact"
**Attacks:** C1 (this was the coordinator's priority-1 line)
**Verdict: FAILS, decisively. This is the most important negative result in my review.**

73% of the small-cryptic set (68.0% of the 71,414 dark-ORF set) is IMG-PR-only, i.e. circular
contigs from metagenome assemblies with no isolate behind them. The project stratified typing and
dispersal by provenance but **never stratified the dark proteome**. I did:

```
                  n         darkonly/noPfam  darkonly/Pfam  named/noPfam  named/Pfam
ALL           378552              57.4           28.5           3.2         10.9
IMGPR_only    295702              56.8           28.4           3.1         11.6
isolate_only   74802              59.5           29.2           3.5          7.8
mixed           8048              59.6           27.8           2.4         10.2
```

The headline cross-tab is **invariant to provenance**, and the isolate-derived stratum is
*marginally darker* (59.5% vs 56.8%). Rebuilding family statistics from isolate-derived members only:

```
--- ISOLATE-ONLY members: 74802 ORFs, 26664 families
    ORFs in fams>=2: 74.9%
    families >=10 members: 1125; Pfam-negative of those: 663 holding 20234 ORFs (27.1%)
    of those Pfam-neg fams, % spanning >=2 lineages: 80.5; median lineages 4.0
```

Recurrence 74.9% (vs 84.7%), 663 substantial Pfam-negative families holding 27.1% of the stratum's
dark ORFs (vs 30.6%), and lineage spread is **higher**, not lower (80.5% vs 72.9%). The shortlist
shrinks 5× because the stratum is 5× smaller, but every proportion holds.

*Viral elements.* The deck concedes "small, ORF-dense, poorly typed replicons are exactly where
phage-plasmids hide. We have not screened for them." I ran the screen, using phage structural
hallmarks (terminase, capsid, portal, tail, baseplate, holin, lysozyme — deliberately excluding
integrase and resolvase, which ICEs and transposons share) propagated from family representatives:

```
dark family reps with a phage structural hallmark: 363 of 24,357 Pfam-positive reps
dark ORFs in hallmark families: 851 = 0.225% of the dark proteome
plasmids carrying >=1 hallmark dark ORF: 676 = 0.95% of 71,414
median size of hallmark carriers: 8042 bp (IQR 6271-11524)
top families: Zot 75 · Phage_capsid 27 · Por_Secre_tail 21 · PhageMin_Tail 20 · Phage_portal 15
              Terminase_4 14 · Phage_Coat_B 13 · Phage_connect_1 12 · Terminase_2 9
filamentous-phage-like (Zot / coat / encapsidation): 101 families, 333 ORFs, 298 plasmids
```

Under 1% of the set carries a detectable phage structural gene. Terminases and major capsid proteins
are among the best-modelled families in Pfam; if these were phage genomes, they would light up. They
do not.

**The bounded residual of this attack, stated fairly:** the screen can only see what Pfam names, and
Inoviridae (filamentous phages, ~6–8 kb, circular, few detectable genes, `Zot` the commonest hit
here) are exactly the class that would evade it. 298 plasmids is a floor, not a ceiling. A geNomad
or CheckV run would close this properly and has not been done. But "a large fraction of the small
payload-free plasmidome is assembly artefact or phage" is **not supported**, and the project's
decision to stratify rather than adjust is vindicated on the one analysis where it never bothered to
apply it.

---

### [MINOR — attack FAILED] F13 — "mash over-splits tiny plasmids, so '224 lineages' is inflated"
**Attacks:** C4, C5
**Verdict: FAILS.**

My hypothesis was that a 4 kb plasmid's mash sketch is nearly empty, so near-identical small plasmids
would scatter across many primary clusters and inflate every lineage count. Using the exact-proteome
duplicate groups from F9 as ground-truth near-identical plasmids:

```
duplicate groups: 5712 ; median members 2
groups whose members span >1 mob_cluster: 34 (0.6%)
distribution of lineages per duplicate group:  1 -> 5678 ; 2 -> 34
max lineages inside a single identical-proteome group: 2
weighted: 0.9% of plasmids in duplicate groups sit in a group spanning >1 cluster
```

mob_cluster essentially never splits identical plasmids. The error runs entirely the other way
(F5, under-splitting). **"Spans 224 distinct primary clusters" is not an over-splitting artefact** —
whatever else those 224 are, they are 224 genuinely distinct plasmid groups. C4's lineage counts are
not inflated by mash behaving badly on small sketches; they are undermined by the *criterion* (F3)
and by clonal over-sampling (F1), which is a different and more repairable problem.

---

### [MINOR] F14 — The 5 UniRef90 survivors withstood everything I threw at them
**Attacks:** C1 (the shortlist)
**Verdict: attack FAILS — recorded because it is the project's strongest asset.**

```
family                                     n    lin top-lin% spp known country kn IMGPRonly%
COMPASS_KP718939.1|8                     138     42      10%       43%        69%        36%
IMGPR_plasmid_2563366806_000011|1        227     30      19%       60%        81%        19%
IMGPR_plasmid_2700989144_000001|10       101     30      18%       42%        69%        40%
IMGPR_plasmid_2551306362_000001|6        121     15      25%       31%        72%        51%
IMGPR_plasmid_2903362858_000001|1        186     14      30%       58%        88%        11%

== the 5 under one-plasmid-per-MOB-cluster dereplication
   COMPASS_KP718939.1|8                 138 plasmids -> 42 after dereplication
   IMGPR_plasmid_2563366806_000011|1    227 plasmids -> 30
   IMGPR_plasmid_2700989144_000001|10   101 plasmids -> 30
   IMGPR_plasmid_2551306362_000001|6    121 plasmids -> 15
   IMGPR_plasmid_2903362858_000001|1    186 plasmids -> 14
```

Dominant lineage 10–30%; every one still on 14–42 carriers after full clonal dereplication; genuinely
multi-country (China/USA/UK/Korea/Spain/Germany); a mix of IMG-PR and isolate provenance. These are
the only objects in Act II that survive F1, F3, F9 and F12 simultaneously. Their weakness is only
that the host labels rest on 31–60% species coverage.

**This is where the project should plant its flag.** Not 3,343, not 115,811, not 57.4%.

---

## IGNORED ALTERNATIVE EXPLANATIONS

Ranked by how much damage each would do if true.

1. **Clonal over-sampling of the carrier population** (F1, F2). Not ignored in the project generally
   — it is the project's signature control — but ignored *here*, on the one result it most threatens.
   Damage: takes recurrence from 84.7% to 39–76% depending on severity, and the shortlist from 3,343
   to 40–1,040.
2. **The residue has homologs in UniProt/NR that Pfam-A simply has no model for** (F4). Never tested;
   the cheapest available test; would directly reduce the 57.4%.
3. **`--cut_ga` suppresses real homology** (F4). Tested by me: 4–27% of the residue rescued above a
   composition-matched null, and the rescues are backbone.
4. **The blind spot is gene count, not crypticity** (F6). Reduces C2's effect size ~2.4×.
5. **`mob_cluster` has no stable biological referent at this size** (F5). Does not flip any sign but
   makes every "N lineages" number uninterpretable as stated.
6. **Small proteins are a known, generic annotation blind spot** (F10). If the residue is "small
   proteins are badly annotated", that is a real finding but a much older one, and not specific to
   plasmids.
7. **The dark fraction is a whole-plasmid phenomenon, not a small-plasmid one** (see "So what?").
8. **Non-mobilizable ≠ immobile.** The project names this itself as the leading alternative for C3
   and proposes the relaxed oriT search against `orit_db_folder/`. Still not run.
9. **Filamentous phages / satellite elements below the Pfam detection floor** (F12 residual).
   Bounded at ≥298 plasmids; true ceiling unknown without geNomad/CheckV.
10. **Expression.** Nothing shows any of these ORFs is transcribed. The deck says this plainly and
    puts it last on its own "Not established" list, which is the right place for it.

---

## MISSING STAKEHOLDER PERSPECTIVES

**Geography.** The dark-family "spans 26 habitats and 30 countries" language rests on this:

```
== GEOGRAPHY of the 71,414 dark-proteome carriers
country known for 79.6% ; 326 countries represented
USA 25555 · China 4838 · United Kingdom 3023 · Canada 2303 · Denmark 1793 · Antarctica 1382 ·
Japan 1258 · Sweden 1247 · Germany 1217 · Australia 882
top 3 countries = 58.8% of located; top 10 = 76.6%
Nigeria 38 · Kenya 151 · Ethiopia 2 · Ghana 78 · South Africa 87 · Tanzania 71 · Uganda 38 · DRC 41
```

The USA alone is ~45% of located carriers. **Antarctica (1,382) contributes more plasmids to this
study than those eight African countries combined (506).** Every geographic-breadth statistic in Act
II — and the "30 countries" for the ColE1 cassette — measures where metagenome sequencing is funded.
The deck says this well for Act I ("it is not showing you where plasmids are, it is showing you which
countries run Enterobacterales surveillance") and then does not carry the caveat into Act II, where
country counts are used as evidence of dispersal.

**Taxonomy.** `plsdb_species` is populated for **14.9%** of the 71,414, and is Enterobacterales- and
clinically-skewed (E. coli 2,231, K. pneumoniae 1,551, S. enterica 568). The "commonest host" column
in the five-survivor table — the deck's single most concrete deliverable — is computed on 31–60%
coverage per family. Anaerobes, Bacteroidota, environmental Bacillota, archaeal plasmids and anything
from a low-income-country sampling programme are structurally absent.

**The 39,331 no-habitat plasmids.** 14,945 of them are inside the 71,414 dark set, and 12,688 of the
71,414 carry `hab_sub = "Unlabelled"`. The habitat-spread statistics (median 4–5 habitats per family)
are computed over a population where roughly a fifth of carriers have no habitat at all. The project
tested these for AMR-informativeness in the Act I reassessment (§6, refuted, 17.86% vs 18.81%) — a
good test — but never for dark-family informativeness.

**Who is misled if this is wrong.** (i) Structural-biology and wet-lab groups who take the 93- or
3,343-family list as a target set: F1 says the prevalence criterion that ranked those lists was
partly counting deposits, so effort would be spent on well-sequenced clones. (ii) Database and tool
authors (PlasAnn, MOB-suite, Pfam) told that 57.4% of a compartment is invisible to them — an
overstated figure sets an unreachable benchmark and may misdirect curation effort toward novelty
when F4 says a chunk of it is sub-threshold backbone. (iii) Anyone citing "71.1% of small cryptic
plasmids are untypeable" as motivation for a new typing scheme, when F6 says a majority of that gap
is gene count and ~9 points is the payload definition. (iv) Reviewers of the eventual paper, if the
93-family list is presented without the dereplication test.

---

## THE "SO WHAT?" VERDICT

Suppose C1 is entirely true as stated. What follows?

The honest answer is: less than the framing implies, because **the dark plasmidome is not a property
of small cryptic plasmids.** From the project's own table:

```
$ cat data/plasmidscope_primary/smallcryptic_dark_matter.tsv
grp                 n_plasmids  total_cds  total_orf  pooled_pct_dark
small-cryptic            43450     234112     210523             89.9
small-cargo              20775     130740      67386             51.5
large-conjugative        25222    3753993    1614486             43.0
large-other              46735    5000639    3667413             73.3
```

Large non-conjugative plasmids are **73.3% dark** and contribute **3,667,413 dark ORFs** — 9.7× the
378,552 this project studied. Summing all four groups gives 5,559,808 dark ORFs in the analysis set,
of which the studied set is **6.8%**. The compartment was chosen for tractability, which is a fine
reason, but the resulting claim is being written as though small cryptic plasmids are where the dark
matter lives. They are not; they are where it is *concentrated per plasmid* (89.9% vs 73.3%), which
is a much weaker statement, and the deck's own "beyond <20 kb" slide concedes the test of whether
"dark = small" is a real property has not been run.

So the strongest version of the "so what" is: *"there exist thousands of uncharacterised protein
families on plasmids" restates microbial dark matter, which has been known for twenty years and has
its own Nature-2024 catalogue (FESNov) that this project cites.* Scoping it to small cryptic plasmids
narrows the population but does not make the finding new, and 6.8% coverage of the available dark
proteome makes the boundary look arbitrary.

**What would make it matter, in ascending order of cost:**

1. Show the residue is *enriched* in this compartment relative to a matched control — e.g. that
   small-plasmid dark families are more lineage-crossing, more conserved, or more compositionally
   distinctive than the 3.67 M dark ORFs on large plasmids. Right now there is no comparison at all.
   This is a day's work and it is the missing keystone.
2. Show it is *bounded* — that the family accumulation curve saturates, so the space is finite and
   worth cataloguing. The project shows the opposite ("half of it needs 3,009 families... there is no
   small core"), which as written argues the space is *not* tractable.
3. Show *selection*: dN/dS on the families with enough nucleotide diversity, or the observation the
   90% notebook already gestures at — ≥90% amino-acid identity across 42 independent lineages is a
   selection argument, and it is the best one in the project. Make it quantitative.
4. Show *expression*: any metatranscriptome hit on any of the five.
5. Show *function*: the cassette knockout.

Items 1 and 3 are cheap and would move this from "we measured how bad annotation is" to "we found
something". Until then the deck's own closing question — *"Is 57.4% a discovery or a measurement of
how bad plasmid annotation is?"* — has to be answered "predominantly the latter", and the project's
answer ("we think it is both") is currently unearned in the first half.

---

## CONFIRMATION-BIAS AND FRAMING AUDIT

**The asymmetry is real and it is large.**

*How the negative result was treated.* The anti-defense hypothesis got: a Pfam screen; a
purpose-built database (dbAPIS, 290 HMMs) reported at **three** stringency tiers so the conclusion
could not depend on the cutoff; a **deliberately aggressive adjudication regex** that the author
notes will discard true positives, making the answer a floor; a per-hit audit against independent
Pfam calls that caught `acriia21` → RepA_N; and then a **third fully independent tool**
(DefenseFinder/AntiDefenseFinder) with a different evidence rule (genomic context), which reproduced
the same two false positives. Result: 460 → 158 → 91.

*How the positive result was treated.* The recurrence claim got: one control (≥2 MOB clusters), no
dereplication, no null model, no sampling-depth analysis, no alternative lineage proxy. The threshold
sweep — which is genuinely excellent work — was run only after the fact and only on the clustering
knob, not on any of these.

That is three independent screens plus false-positive adjudication for the hypothesis they wanted to
kill, and one weak criterion for the hypothesis they wanted to keep. It is the textbook signature.

**Specific framing findings, deck vs reports.**

| # | Location | Finding |
|---|---|---|
| A | Deck Part 5, "the cryptic set" | "Replicon call from mob_typer on only 35.0%; 65% are untypeable" — pooled. Stratified it is IMG-PR-only **76.5%** vs isolate-derived **34.6%** (my recompute). The sibling report carries a provenance stratification for the 47,031 set; the deck carries none for the 71,414 set, and the speaker note escalates it to "the standard typing toolkit is already blind to two thirds of this set". |
| B | Deck Part 7 verdict + speaker note | "**46.9%** of what Pfam rescued... That is **18.5%** of the entire dark proteome." `reports/pfam_dark_validation.md` §4.3 says **47.3%** and **18.6%**. Small, but the deck is the artefact that circulates. |
| C | Deck Part 8 ledger vs Part 7 verdict | The ledger slide totals to "still unexplained **59.3%**"; the verdict slide says **57.4%**. Two different partitions (39.4% family-size-weighted Pfam+ vs the four-cell cross-tab) presented in one deck without reconciling. Both are correct; a reader cannot tell why they differ. |
| D | Deck errata slide, Defect 2 | Deck: "The index was computed over all features rather than over CDS only." Report §6.2: the index was numbered *within the emitted subset*, advancing only on dark CDS in `--mode dark`. These are different bugs. The errata slide is the most credibility-earning slide in the deck and it misdescribes the defect. |
| E | Deck errata slide, Defect 1 | Deck: "Fixed by de-duplicating on plasmid ID before extraction." Report: first-occurrence-wins dedup **on the FASTA header** in `scripts/cluster_dark_orfs.sh`. Again a different fix. |
| F | Deck Part 11 | The recommended candidate is "**dORF0356**". `widespread_dark_orfs.ipynb` calls the same sequence "**dORF0350**". The sequence matches; the identifier does not. |
| G | `reports/small_cryptic_methodology.md` §4 | "for small-cryptic the rate *falls* with increasing distance (17.4% in the >0.1 bin)". The project's own stored table peaks at **42.3%** in the ≤0.01 bin (F8). Endpoint selection. |
| H | `reports/small_cryptic_methodology.md` §3 | Pooled dispersal result (significant, cryptic narrower) leads; isolate-only null follows and is then converted into "near-parity" (F7). The K = 10 ratio is quoted; the K = 40 ratio, twice as large, is not. |
| I | `PROJECT_OVERVIEW.md` §3 | "there are **no duplicate plasmids** across source databases" — true only for 100% nucleotide identity; 29.9% of carriers share an exact proteome with another (F9). |

**Where the framing is better than the field standard, and should be said.** The deck leads Act II
with a refuted manuscript thesis and a rejected hypothesis. It devotes an entire slide to a bug that
inflated its own novelty estimate by ten points, with the line "failing loudly would have been better
than the answer we liked". It presents the 93 → 5 collapse on its own slide and calls 93 an upper
bound in three separate places. Its "Not established" column is longer and more damaging than most
papers' limitations sections, and it correctly places "nothing here is expression evidence" last, as
the thing that should bother them most. `reports/clustering_threshold_sensitivity.md` §4 has an
explicit "What does not survive" section. That is unusually honest work, and findings A–I are errors
of emphasis inside a fundamentally candid document — not spin.

**Two things the deck does *not* do that I expected to find, and did not:** it does not present the
71.1% typing blind spot, the dispersal paradox (C3), or the host-range inversion (C5) at all. Those
claims live only in `reports/small_cryptic_methodology.md`, which carries their caveats in the same
section. So the cherry-picking charge on C3 and C5 **does not land against the deck** — there was
nothing to cherry-pick from. F7 and F8 are criticisms of the report, not of the presentation.

---

## THE COLLAPSING SHORTLIST (attack line 6)

The sequence is 93 → 61 (criterion enforced) → 39 (UniRef50) → 19 (NMPFamsDB) → 5 (UniRef90), and
from 50% upward the most prevalent family is one Pfam can name. Add my F1 result and the 30%-identity
93 becomes ~30 under clonal dereplication.

**The case that this is a permissive-threshold artefact.** A finding that loses 95% of its members
when the family definition is tightened along a standard axis, and whose flagship object (88 MOB
lineages, 2,468 plasmids) becomes 5 proteins in 1 lineage at 70%, was a property of the clustering,
not of a protein. The project says exactly this, in those words. The 30% + `-c 0.8 --cov-mode 0`
combination is *loose on identity and strict on coverage* — an unusual pairing that the project
itself identifies as the cause of the non-nesting problem. Under a hostile reading, 93 was never a
result; it was a shortlist generated by a greedy set-cover at a threshold that maximises apparent
family size, and the honest number is 5.

**The counter-case, which I think is stronger than it first looks.** `zoom_two_dark_orfs.ipynb` §6
demonstrates that the clusterings are **not nested**: ORF-A is a 20-protein family at 30% and a
411-protein family at 90%, because symmetric 80% coverage forbids a 208 aa hub and a 155 aa protein
from ever sharing a family, and at 30% the hub absorbs everything it covers. So 93 → 5 is not monotone
attrition of a single object — it is two different partitions of the same sequence space, and the
90% run *discovers* families the 30% run hid. The project found this itself and describes it as "the
methodological finding I would most like feedback on". That is the correct reading and it is a real
contribution.

**My verdict:** the honest number is neither 93 nor 5 — it is that the family-count deliverable is
not well defined under this method, and the project should stop reporting one. What *is* well defined
is the per-candidate evidence, and the five UniRef90 survivors carry it (F14). Report the five,
report the cassette, report the sensitivity table, and drop "3,343 families / 115,811 ORFs" from the
abstract position it currently occupies.

---

## ATTACK LINES NOT TESTED

Stated plainly so the panel can see the gaps in my coverage.

- **C6 (the surviving Arc-1 One Health claim).** **NOT TESTED.** I read
  `manuscript/reassessment_round1.md` and `review_round1.md` but re-fit no model and recomputed no
  odds ratio. My only reading-based observation, which the project already makes itself: host
  coverage is 99.2% for clinical against 20.1% for wildlife, so the ~50% attenuation is estimated on
  a subset over-representing the highest-effect compartments, and the project's own hedge ("the
  precise fraction uncertain") is the right one. I add nothing to it.
- **Whether the circular contigs are plasmids at all**, by a purpose-built classifier (geNomad,
  CheckV, PlasClass). **NOT TESTED** — no such tool or database is installed. My Pfam-hallmark screen
  (F12) bounds the phage signal but cannot see anything Pfam has no model for.
- **DIAMOND/BLAST against UniRef90 or NR.** **NOT TESTED** — no such database on disk. This is the
  most consequential of my gaps (F4).
- **The DefenseFinder and dbAPIS adjudications**, and the FlashFold/candidate selection. **NOT
  TESTED** — I read the notebook narratives and did not re-run any search.
- **Figures.** **NOT TESTED** — I read no figure and verified no panel.
- **Anything upstream of the master table** (PlasAnn gene calls, RGI, MOB-suite internals, BioSample
  scraping) — taken on trust, as the verification pack also flags.
- **`--cov-mode` / coverage sensitivity.** Only identity was swept by the project; I did not sweep
  coverage either, though F14's non-nesting result suggests it is the more consequential knob.

---

## WHAT WOULD CHANGE MY MIND

My strongest counter-argument is F1 + F2: *the recurrence result is a readout of clonal
over-sampling at a particular sampling depth, and the control that would show this was never run.*

Four analyses would defeat it, in order of decisiveness:

1. **Rebuild the entire dark-family pipeline on a dereplicated carrier set** — one plasmid per MOB
   cluster, or better, ANI-dereplicated at 99% with skani/fastANI since `mob_cluster` is not fit for
   this purpose (F5). Re-cluster, re-roll-up Pfam, re-derive the shortlist. If the Pfam-negative
   substantial-family count lands materially above my ~40 and the residue fraction stays near 30%,
   F1 is answered and C1/C4 stand as written. If it lands near 40, the headline numbers must be
   restated. **This is the one experiment I would demand before publication.**
2. **A sampling-depth null.** Simulate World A explicitly: junk sequence evolving on the observed
   plasmid phylogeny/cluster structure, sampled at the observed depth, clustered identically. If it
   produces recurrence well below 84.7%, the World A/World B test becomes a real test and F2
   dissolves. If it produces ~80%, the framing must be abandoned.
3. **DIAMOND against UniRef90 (and ideally NR) for the 3,343 representatives**, reported alongside
   the Pfam result. If the residue is still unnamed at 30% identity / 50% coverage against UniRef90,
   F4's second half is answered and the novelty claim becomes much harder to attack. Hours of
   compute. Its absence is currently the easiest thing for a hostile reviewer to point at.
4. **A matched comparison against the large-plasmid dark proteome** (the 3.67 M ORFs already
   extracted, `data/fam_3300022589/big_dark_orfs.faa`). If small-plasmid dark families are
   demonstrably more lineage-crossing or more conserved than large-plasmid ones, the "so what"
   collapses and the compartment choice becomes a finding rather than a convenience.

Two smaller things would move me on individual claims: a **TOST equivalence test** on the isolate-only
dispersal data would convert F7's "near-parity is unsupported" into either a real parity result or a
clean concession; and reporting the **overlap-region** typing comparison (F6) instead of the crude
9.1× would make C2 defensible as stated.

Finally, what would *not* change my mind: another HMM database. The ledger slide already shows that
every database screen run so far moves the explained fraction by 1.3 points. The binding constraints
are the carrier population and the family definition, not the model set.

---

*Ends. Findings F11, F12 and F13 are attacks that failed under test and should be weighted as
evidence in the project's favour.*
