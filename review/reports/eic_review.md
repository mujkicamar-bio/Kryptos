# Editor-in-Chief — internal rigour audit

**Reviewer seat:** EIC, microbial genomics (resources + methods). One of five independent reviewers;
written without sight of, or reference to, the other four.
**Scope:** the whole research programme, both arcs, plus the group-meeting deck.
**Standard applied:** *is this true, reproducible, and does it support what it claims* — not venue fit.
**Date:** 2026-09-02
**Built on:** `review/VERIFICATION_PACK.md` (46/52 headline quantities reproduce exactly). I have not
re-litigated arithmetic. Where I ran my own computations they are new quantities the project has
never computed, and each is shown with the command context that produced it.

---

## 1 · Summary of the programme (150 words)

A single 208,248-plasmid resource — complete, non-redundant plasmids from PlasmidScope, folded into a
67-column master table across seven annotation layers — has been used to chase three successive
theses. The founding question was multireplicon architecture, inherited from de Quinto et al. Attention
moved to plasmid-borne AMR across One Health compartments; that produced a full NAR draft whose
central claim ("mobility rivals provenance") the project's own simulated peer review and reassessment
then refuted, showing it to be a plasmid-size artefact. Attention moved again to the "dark
plasmidome": 71,414 payload-free small plasmids carrying 378,552 proteins PlasAnn can only call `ORF`.
That arc established that the dark proteome is recurrent rather than idiosyncratic, partitioned it
against Pfam-A into ~28.5% annotation failure and ~57.4% genuine residue, rejected an anti-defense
hypothesis from three independent directions, and produced a shortlist that collapsed from 93 families
to 5 under its own sensitivity analysis.

---

## 2 · Overall assessment (300–400 words)

This is not a sequence of competent analyses in search of a thesis, and it is not thesis drift. It is
a research programme with one durable object — the master table — and a genuinely healthy
self-correction record on top of it. The claim history reads as escalating honesty, not escalating
ambition: each pivot was forced by a *specific refutation the project itself produced and wrote down*.
The AMR arc did not fade out, it was killed by `reassess_review_findings.py` and the death certificate
is in the repository. The anti-defense hypothesis was killed by three independent screens, two of
which reproduced each other's false positives. The 93-family shortlist was killed by a compliance
audit the project ran on itself after reading MMseqs2's own help text. Programmes that drift do not
generate this kind of paper trail. This one converged: the question got narrower and better-posed at
every step.

The problem is not the science. The problem is that **the correction pipeline stops at the point of
discovery.** Defects are found, quantified, and written up beautifully — in a *new* document. The old
document is left standing, uncorrected and un-annotated, and it keeps being cited. Three separate
instances of this are documented below (F1, F3, F8), and the most serious is that
`reports/amr_onehealth_methodology.md` §5 still asserts "**Mobility rivals provenance**" as an
established finding twenty-seven minutes' worth of file-timestamps after the reassessment refuted it —
and the deck cites that report as the source for the AMR results it presents. A reader arriving at any
single document in this repository cannot tell whether they are holding a live number or a dead one.
For a programme whose defining virtue is self-correction, this is the one defect that undermines the
virtue itself.

The second-order problem is a systematic miscalibration in *which direction* the project doubts
itself. It over-doubts the resource layer (presenting the Arc-2 catalogue as if the "93 → 5 collapse"
threatened it, when it does not) and under-doubts the resource layer's actual defect (`pf_n_inc`,
which is documented as something it is not). I recomputed both. The catalogue is far more robust than
the programme claims; the multireplicon column is far more broken than the programme admits.

There is a publishable paper here. It is not the one the project has been drafting.

---

## 3 · The publishable unit

**What is publishable right now is the Arc-2 catalogue paper, not the Arc-1 AMR paper.** Arc 1's
manuscript has a refuted title, a refuted headline claim, a section its own reviewers told it to cut,
and a resource-layer column defect (F2) sitting under its founding question. Arc 2 has a positive,
external-database-validated, threshold-robust result and a concrete deliverable. Below is the title
and abstract I would defend.

### Proposed title

**A recurrent, lineage-spanning protein space in the small plasmidome is invisible to both
plasmid-specific annotation and Pfam**

### Proposed abstract (every number below reproduces; provenance noted)

> Small plasmids are the modal element in complete-plasmid databases and the least characterised. We
> assembled 208,248 complete, non-redundant bacterial plasmids from PlasmidScope and annotated them
> uniformly (PlasAnn, MOB-suite, CARD/RGI). Within the non-artefact analysis set, 71,414 plasmids are
> under 20 kb and carry no annotated resistance, virulence or metal/biocide gene; these contribute
> 378,552 proteins that PlasAnn can only label "open reading frame" — 82.8% of their translated
> proteome. Clustering this dark proteome against itself (MMseqs2, 30% identity, 80% mutual coverage)
> yields 92,752 families, of which only 15.3% are singletons: 84.7% of dark ORFs sit in families with
> at least one other member, and 63.9% sit in families spanning two or more independent MOB lineages,
> so the recurrence is not clonal redundancy. Searching every family representative against Pfam-A
> 38.2 at Pfam's curated gathering thresholds — with the 78,397 PlasAnn-named proteins from the same
> plasmids as a positive control, which hit at 76.7% against 26.3% for dark representatives —
> partitions the dark proteome exactly: 28.5% is a PlasAnn annotation failure that Pfam corrects, of
> which roughly half carries a domain in a category PlasAnn already annotates, and 57.4% is named by
> neither tool. After Pfam, 3,343 families of ≥10 members retain no Pfam domain at all, holding
> 115,811 dark ORFs (30.6% of the dark proteome), 72.9% of them spanning ≥2 lineages. This catalogue
> is robust to the clustering choice: re-clustered with the criterion strictly enforced
> (`--cluster-reassign`) at the same identity it comprises 3,490 families and 112,413 ORFs, sharing
> 94.8% of its ORF content with the original (Jaccard 0.928), and it retains 2,878 families at
> UniRef90 stringency. Three independent screens (Pfam-A, dbAPIS, DefenseFinder/AntiDefenseFinder)
> reject anti-defense as an explanation for this residue, which accounts for at most 0.042% of it;
> defense systems, which the payload-free definition never excluded, account for 1.32%. The residue is
> therefore neither annotation noise nor a known functional class, and profile-HMM search has reached
> its limit on it.

**Numbers that survive and carry the paper:** 208,248 · 71,414 · 378,552 · 92,752 · 15.3% / 84.7% ·
63.9% · 76.7% vs 26.3% · 57.4/28.5/10.9/3.2 · 3,343 / 115,811 / 30.6% / 72.9% · 70.3% recurrence at
UniRef90 · 38.6–39.1% Pfam-named across all thresholds · 0.042% and 1.32%.

**Numbers that must NOT appear:** 24,583 multireplicon (F2) · 95,442 replicon coverage (F9) ·
"mobility rivals provenance" and everything downstream of it (F1) · the 93-family shortlist presented
as a result rather than as a worked example of threshold-dependence (F5) · "88 MOB lineages" for
`IMGPR_plasmid_3300022589_000018|2` (a property of the threshold, as the project itself established).

**The one new sentence I would add, which the project has not written:** the *catalogue* is
threshold-robust while the *shortlist* is not, and conflating the two has made Arc 2 look weaker than
its evidence. See F4.

---

## 4 · Findings

### [CRITICAL] F1 — Refuted claims still stand as established findings in the document the deck cites

**Evidence anchor:** `reports/amr_onehealth_methodology.md:243–252`, §5 "What this establishes":

> "2. **Mobility rivals provenance.** Conjugative plasmids carry AMR at 60.1% overall, and the
> conjugative fraction of the natural environment out-carries the non-mobilizable fraction of the
> clinic."
> "5. The **multireplicon–AMR association is universal but strongest outside the clinic**, extending
> de Quinto et al. into exactly the compartments their cultured/clinical dataset could not reach."

**What is wrong.** `manuscript/reassessment_round1.md:44` states the opposite of item 2 — "At every
matched size the clinical non-mobilizable plasmid is the more likely carrier, by roughly 1.5–1.8-fold
… The sentence the paper is named after is an artefact of size composition" — and §5 of that same file
recommends *cutting* item 5 outright. The reassessment is timestamped 2026-07-27 15:04; the
methodology report 2026-07-27 14:37. It was never touched again. `grep -n "reassess\|review_round1\|
refut" reports/amr_onehealth_methodology.md` returns **nothing**: the file carries no forward pointer,
no banner, no strikethrough. Its §3.3 heading still reads "The two axes are **independent**", which the
reassessment showed to be false (61% of the raw conjugative effect is plasmid size). This is the file
the group-meeting deck cites as its source for the AMR gradient ("analysis · reports/
amr_onehealth_methodology.md §3.1"). The deck's own peer-review slide is honest; the document
underneath it is not.

**Fix.** Add a dated header block to `reports/amr_onehealth_methodology.md`: *"Superseded in part on
2026-07-27 — see `manuscript/reassessment_round1.md`."* Then edit §3.3's heading to "separable but
correlated", strike §5 items 2 and 5, and replace them with the reassessment's three surviving
results. Do not delete the originals; mark them. The project's whole credibility model is that
corrections are visible.

---

### [CRITICAL] F2 — `pf_n_inc` is documented as something it is not, corrupting the founding question

**Evidence anchor:** `reports/master_table_data_dictionary.md:129`:

> "| 52 | `pf_n_inc` | float | 100% | Count of distinct Inc **families** (`0` = untyped/none).
> **`pf_n_inc ≥ 2` = multireplicon** — the paper's definition. |"

Against `scripts/harvest_typing.py:78`, which computes `"pf_n_inc": len(incs)` where `incs` is the
sorted list of **allele-level** names, while the family-collapsed `fams` is written to a separate
column and never counted.

**What is wrong.** I verified this directly against the master table rather than the code alone. For
**every one of the 208,248 rows**, `pf_n_inc` equals the allele count exactly and never equals the
family count. 8,354 rows — **20.1% of all PlasmidFinder-typed plasmids** — have an inflated value. The
headline multireplicon figure is 24,583 by allele and **22,422 by family (−8.8%)**. Worked examples:

```
GenBank_MN182750.1   pf_n_inc = 24   families = IncFIA;IncFIB;IncFII;IncHI1B;repB_KLEB_VIR_AP006726;repE   (6)
GenBank_CP035124.1   pf_n_inc = 25   families = IncFIA;IncFIB;IncFIC;IncFII;IncQ1;repB;repE                (7)
```

This is worse than an inflated statistic. It is a **resource-integrity defect**: the authoritative
column reference actively misdescribes the column, and instructs a reuser to threshold on it as a
family-level multireplicon flag. Anyone who downloads this table and follows its own documentation
will compute a wrong number and have no way to discover it. The defect sits directly under the
question `PROJECT_OVERVIEW.md:14–16` names as "the central biological question driving the design".

It also contradicts a locked decision. `EXECUTION_PLAN.md` "Decisions locked in (2026-06-30)" item 3
states mob_typer "is the primary replicon … typing tool — **not** PlasmidFinder", and the deck states
the standing rule as "every replicon / Inc number in this deck comes from mob_typer's
`mob_rep_types`". Yet `PROJECT_OVERVIEW.md:92` and `reports/typing_methodology.md:44,94` both publish
the PlasmidFinder allele-based 24,583 without qualification. The deck's own multireplicon figure
(18,685 / 24.7% of typed, mob-based) is a different number entirely, and the two are nowhere
reconciled.

**Fix.** Three edits, in this order. (a) Correct `master_table_data_dictionary.md:129` to "Count of
distinct Inc **allele names**; use `pf_inc_families` for family-level analysis" and remove "the
paper's definition". (b) Add a `pf_n_fam` column in `harvest_typing.py` (`len(fams)` — the value is
already computed and discarded) and rebuild the master; this is a one-line change. (c) Replace 24,583
with 22,422 in `PROJECT_OVERVIEW.md:92` and `typing_methodology.md:44,94`, or delete the PlasmidFinder
multireplicon figure from both in favour of the mob-based 18,685, consistent with the locked rule.

---

### [CRITICAL] F3 — The dark-clustering corrections were never propagated; three documents still call the sweep un-run

*This is the specific check the audit asked for: were `reports/dark_orf_clustering.md` numbers
regenerated after the `--cluster-reassign` finding? **No — and the decision not to regenerate them is
defensible, but the failure to annotate them is not.***

**Evidence anchor:** `reports/dark_orf_clustering.md:216–218`, Limitations item 2:

> "2. **One clustering threshold.** 30% identity / 80% coverage only. Family counts are threshold-
> dependent; the qualitative split (recurrent vs. singleton) should be checked at 40% and 50%
> **before publication**."

and `reports/pfam_dark_validation.md:210`:

> "5. **One clustering threshold** upstream (30% identity / 80% coverage); untested at 40% / 50%."

and `notebooks_my/dark_plasmidome.ipynb` cell 21:

> "- **One clustering threshold** (30% id / 80% cov); untested at 40% and 50%."

**What is wrong.** The sweep *was* run, on 2026-09-01, at 30/50/70/90% with the criterion enforced,
and reported in `reports/clustering_threshold_sensitivity.md`. All three statements above are now
false. Worse, none of the three documents mentions the *other*, more serious finding — that the
published 30% clustering was enforcing its stated criterion on only **89.12%** of member–representative
pairs. `dark_orf_clustering.md` §2 (Method) still presents "30% identity with 80% bidirectional
coverage (`--cov-mode 0`)" as what was enforced, which the project's own audit has since shown it was
not. Verified by cross-reference sweep:

```
grep -rln "clustering_threshold_sensitivity" reports/ plans/ *.md manuscript/ notebooks_my/
  → dark_families.ipynb, dark_families_90pct.ipynb, plans/conservative_reclustering.md, dark_plasmidome_deck.html
grep -rn "cluster.reassign" reports/dark_orf_clustering.md reports/pfam_dark_validation.md
  → (no matches)
```

The two reports that *define* the dark proteome are the only Arc-2 documents that do not know the
sensitivity analysis exists. `dark_plasmidome.ipynb` was modified on 2026-09-02 — the day *after* the
sweep — and still was not updated.

The decision itself is right and well-argued: `clustering_threshold_sensitivity.md:18–20` deliberately
preserves the original outputs so the earlier reports "remain reproducible". Preserving the artefacts
is correct. Leaving a stale "should be checked before publication" in a report that has since been
checked, and a Method section describing a criterion that was not enforced, is not.

**Fix.** Add a three-line dated note at the top of `dark_orf_clustering.md` and
`pfam_dark_validation.md`: *"2026-09-01 — the threshold sweep called for in Limitations has been run,
and an audit found this clustering enforced its stated criterion on 89.12% of sampled pairs. See
`reports/clustering_threshold_sensitivity.md`. Numbers below are unchanged and remain reproducible
from the preserved `dark30_*`/`mix30_*` outputs."* Rewrite the stale limitation in all three
documents to point at the sweep instead of requesting it. Amend `dark_orf_clustering.md` §2 to state
the enforced-vs-documented criterion gap.

---

### [MAJOR] F4 — The flagship deliverable was never re-derived under the criterion fix, and it is far more robust than the programme claims

**Evidence anchor:** `reports/pfam_dark_validation.md:229` — "`dark_novel_families.tsv` | **the
deliverable** — 3,343 substantial Pfam-negative families" — built entirely on the 30% clustering that
F3 shows was under-enforcing its criterion. `plans/conservative_reclustering.md` §5 "Out of scope"
does not list re-deriving it, and `clustering_threshold_sensitivity.md` §4 tracks only the *core*
shortlist (93 → 61 → 39 → 19 → 5), never the catalogue.

**What is wrong.** The single most important Arc-2 output has no reported sensitivity, while the
much less important shortlist has a whole report devoted to its collapse. The result is that the
programme's own narrative — deck Part 10, "Two conclusions are threshold-proof. The shortlist is
not." — leaves the catalogue's status ambiguous, and the deck's summary slide lists "The shortlist
moves 93 → 5 across the sweep" under **Not established**, immediately below the 3,343 catalogue under
**Supported**, with nothing said about whether the collapse touches the catalogue.

I computed it. Applying the identical criterion (≥10 members, no Pfam-A hit at `--cut_ga`) to the
re-clustered `data/dark_orf_run/recluster/dark{30,50,90}_fampfam.tsv`:

| clustering | families | dark ORFs | % of 378,552 |
|---|---:|---:|---:|
| 30%, as published | 3,343 | 115,811 | 30.6% |
| **30% + `--cluster-reassign`** | **3,490** | **112,413** | **29.7%** |
| 50% (UniRef50) + RA | 3,352 | 103,593 | 27.4% |
| 90% (UniRef90) + RA | 2,878 | 77,031 | 20.3% |

And the identity of the catalogue is preserved, not just its size — comparing ORF membership between
the published catalogue and the criterion-enforced rebuild at the same identity:

```
ORF-level overlap 109,829 · 94.8% of the original catalogue's ORFs retained
                            97.7% of the corrected catalogue already present
                            Jaccard 0.928 · 3,053 of 3,343 representative IDs shared
```

**The catalogue is stable to within ~5% under the very fix that destroyed the shortlist.** This is a
materially stronger result than anything the project currently claims for Arc 2, and it exists only
because the shortlist and the catalogue were never separated in the sensitivity analysis. The
programme is under-claiming here — which is the mirror image of F2, where it over-claims.

**Fix.** Add these two tables to `reports/clustering_threshold_sensitivity.md` §4 as a new subsection
"Sensitivity of the catalogue (as distinct from the shortlist)". Regenerate
`dark_novel_families.tsv` at 30%+`--cluster-reassign` as the canonical deliverable (3,490 rows),
retaining the original alongside it. Then rewrite the deck's Part 13 ledger so "catalogue: robust" and
"shortlist: threshold-dependent" are two separate bullets rather than adjacent ambiguities.

---

### [MAJOR] F5 — The shortlist has a blind spot the project documented and then did not act on

**Evidence anchor:** `notebooks_my/zoom_two_dark_orfs.ipynb` §6, and deck Part 12 "the sting":

> "The 30% view understated ORF-A twentyfold — and it would have been missed entirely by ranking
> families at that threshold, which is exactly how the 93-family shortlist was built."

**What is wrong.** This is an excellent and generalisable methodological finding — greedy set-cover
under a symmetric 80% coverage constraint makes long proteins into hubs that strand short ones — and
the project states it plainly and then draws no consequence from it. Two consequences follow that are
never drawn. First, the "93 → 61 → 39 → 19 → 5" sequence is presented as a *collapse of one set*, but
if the clusterings are not nested (as §6 proves), the sweep can also **surface** core-eligible families
that the 30% ranking never saw; the project never checked whether it did. Second, the 3,343 catalogue
inherits the same construction and the same blind spot.

I checked the first of these, because the programme's headline "93 → 5" depends on it. Mapping the
members of each 90% core family back through `dark30_cluster.tsv`:

```
90%-core COMPASS_KP718939.1|8       → 30% family GenBank_CP096868.1|3        (138/140)  in_93_core=True
90%-core IMGPR_..._2551306362|6     → 30% family IMGPR_..._2551306362|6      (126/126)  in_93_core=True
90%-core IMGPR_..._2563366806|1     → 30% family IMGPR_..._2660238460|5      (228/230)  in_93_core=True
90%-core IMGPR_..._2700989144|10    → 30% family IMGPR_..._2700989144|10     (105/105)  in_93_core=True
90%-core IMGPR_..._2903362858|1     → 30% family DDBJ_AP027707.1|1           (187/187)  in_93_core=True
```

**All five 90% survivors are genuinely five of the ninety-three.** The "collapse" framing is therefore
correct at the endpoints — but the project asserted it without ever testing it, on the same page where
it proved that clusterings at different thresholds are not nested. That is an unearned claim that
happens to be true.

**Fix.** State the nesting check (above) explicitly in `clustering_threshold_sensitivity.md` §4 —
it takes ten lines and converts an assertion into a result. Then add one sentence to
`dark_families.ipynb` §5 acknowledging that the shortlist's construction inherits the hub effect from
§6 of `zoom_two_dark_orfs.ipynb`, so 93 is not only an upper bound on stringency but a possibly
incomplete draw on short proteins.

---

### [MAJOR] F6 — The front-door document does not know half the programme exists

**Evidence anchor:** `PROJECT_OVERVIEW.md:3` — "*Last updated 2026-07-27. **This is the front-door
document**: what the project is, what has been built so far, where everything lives*". Verified:

```
grep -in "dark\|cryptic\|pfam\|notebooks_my" PROJECT_OVERVIEW.md   →   (no matches)
```

**What is wrong.** The entire dark-plasmidome arc (2026-08-26 → 2026-09-02: four reports, one plan,
nine notebooks, the 56-slide deck, `data/dark_orf_run/`, `data/pfam_run/`) is invisible from the
document that declares itself the entry point. §7's repository map lists "reports/ 6 methodology
write-ups" (there are now ten plus a deck) and "notebooks/ 5 investigation notebooks" (there is also
`notebooks_my/` with nine). §2's status table still ends at Phase 3, and its manuscript row —
"⚠️ draft, needs reframing" — is the *only* place in the repository where the reader is warned that
`manuscript/amr_onehealth_NAR.md` is refuted. A reader who opens the manuscript directly gets no such
warning (F8).

**Fix.** Update `PROJECT_OVERVIEW.md` to cover both arcs, with a status table row per Arc-2 stage and
a repository-map entry for `notebooks_my/`, `data/dark_orf_run/` and `data/pfam_run/`. Add an
explicit "Claim status" section listing which headline claims are live, which are superseded, and by
what — this is the structural fix for F1, F3 and F8 together.

---

### [MAJOR] F7 — No version control, no archive plan, and a resource that cannot yet be handed to anyone

**Evidence anchor:** `ls -d .git` → **no git repository**. `grep -rin "licen\|zenodo\|figshare\|
deposit"` across `reports/*.md`, `manuscript/*.md` and the root docs returns no archiving, licensing
or DOI plan for the master table or the catalogue.

**What is wrong.** The programme's stated virtue is that "every number is produced by a named script"
(`reports/dark_orf_clustering.md:4`, and the standing rule in `EXECUTION_PLAN.md` decision 6). That
guarantee is only as good as the ability to say *which version* of the script. Without version
control there is no way to reconstruct that `harvest_typing.py` as it stands today is the one that
produced the master table folded in July, or to see the two dark-ORF bug fixes as diffs rather than as
prose. The corrections in `pfam_dark_validation.md` §6 are described in exemplary detail and are
nonetheless unverifiable, because the pre-fix code no longer exists anywhere.

Separately, for the resource contribution assessed under §5 below: the master table and the catalogue
are currently a path on one HPC filesystem. There is no licence, no checksum manifest, no archived
copy, and — for `dark_novel_families.tsv`, the declared "deliverable" — no column dictionary and no
documented route from its `rep` IDs to the actual sequences (they are joinable against
`data/dark_orf_run/dark30_rep_seq.fasta`, but this is stated nowhere).

**Fix.** (a) `git init` now and commit the current state as a baseline, even retroactively — the
history is lost but the future is not. (b) Write a six-line column dictionary for
`dark_novel_families.tsv` into `pfam_dark_validation.md` §8 and state the FASTA join. (c) Before any
submission, deposit master table + catalogue + representative FASTA on Zenodo with a licence and a
checksum manifest, and cite the DOI from `PROJECT_OVERVIEW.md`.

---

### [MAJOR] F8 — The refuted manuscript carries no warning in its own file

**Evidence anchor:** `manuscript/amr_onehealth_NAR.md:1` — the title is still "**Mobility rivals
provenance**: plasmid-borne antimicrobial resistance across One Health compartments in 104 169
complete plasmids", and the abstract still closes:

> "Conjugative plasmids from natural environments carry resistance more often than non-mobilizable
> plasmids from the clinic, indicating that transfer capability deserves parity with sampling origin
> in genomic AMR surveillance."

The file was last modified 2026-08-13, i.e. **seventeen days after** `reassessment_round1.md` refuted
exactly this sentence, and it contains no reference to the review, the reassessment, or its own
status.

**What is wrong.** This is a lower-severity instance of F1 (the manuscript is at least flagged in
`PROJECT_OVERVIEW.md`'s status table, whereas the methodology report is flagged nowhere), but it is
the document most likely to be circulated to a collaborator, and it is the one that would do the most
damage if it escaped. It also still contains §"Multireplicon architecture and resistance outside the
clinic", which both reviewers and the reassessment said to cut, and Table 4 *P* values the
reassessment showed to be overstated 2–9×.

**Fix.** Add a header block immediately under the title: *"SUPERSEDED 2026-07-27. The central claim of
this draft is refuted by `manuscript/reassessment_round1.md`. Retained for provenance; do not
circulate."* Do not silently retitle — the refutation is part of the record.

---

### [MINOR] F9 — "95,442 (46%) any-standard-threshold replicon coverage" does not reproduce and is published twice

**Evidence anchor:** `PROJECT_OVERVIEW.md:95` and `reports/typing_methodology.md:93`, both stating
"**Any standard-threshold replicon (PlasmidFinder ∪ MOB-suite ∪ PLSDB): 95,442 (46%)**". The
verification pack §7.3 recomputes the natural reconstruction as 97,668 (46.9%).

**What is wrong.** The exact set expression behind 95,442 is recorded nowhere, so the number cannot be
audited or regenerated — which is a direct violation of the project's own standing rule that every
report cite the script producing it.

**Fix.** Record the exact expression in `typing_methodology.md` §3 or replace both figures with the
reproducible 97,668 (46.9%).

---

### [MINOR] F10 — Two "payload-free small plasmidome" definitions share one vocabulary, and the deck misattributes its own threshold

**Evidence anchor:** `reports/small_cryptic_methodology.md:36` — "**small** — `size_bp < 10,000`" with
a conjugation filter → 47,031. Against `reports/dark_orf_clustering.md:51` — "Small = `size_bp <
20,000`", no conjugation filter → 71,414. And `reports/dark_plasmidome_deck.html` speaker note,
Part 5:

> "If asked why 20 kb: it is the size threshold in the locked small-cryptic methodology, and the size
> distribution is genuinely bimodal around it."

**What is wrong.** The locked small-cryptic methodology specifies **10 kb**, not 20 kb, and explicitly
says the 10 kb cut "is *not* the antimode: it is a deliberately conservative line inside the small
mode's shoulder". The deck's fallback answer to the most obvious question about its scoping is
factually wrong about its own source document. The 20 kb choice is separately defensible — the antimode
is at 18,757 bp — but that is not the justification given. Credit where due: the deck otherwise handles
this cleanly, using the 71,414 object consistently and never quoting the 47,031 scale figures beside
it, which is more discipline than the reports show.

**Fix.** Correct the speaker note to justify 20 kb on the measured antimode (18,757 bp) rather than by
citing a report that specifies 10 kb. Give the two sets distinct names in both reports —
e.g. *small-cryptic set* (<10 kb, conjugation-filtered, n=47,031) and *dark-ORF working set* (<20 kb,
unfiltered, n=71,414) — and use them consistently.

---

### [MINOR] F11 — The proposed wet-lab target is named by two different IDs

**Evidence anchor:** `notebooks_my/widespread_dark_orfs.ipynb` cell 16 — "**The one actually worth
working on** — `dORF0350`, **75 aa, 41 plasmids, 11 MOB lineages, 12 habitats**". Deck Part 11 — "The
one worth working on: `dORF0356` — 75 aa, 41 plasmids, 11 MOB lineages, 12 habitats". Identical
statistics, identical amino-acid sequence, different identifier.

**What is wrong.** Small, but this is a candidate the deck asks the room to fold and synthesise; the
ID is its handle. `grep -o "dORF[0-9]*"` returns `dORF0350` from the notebook and `dORF0356` from the
deck, with `dORF0001` and `dORF0105` matching in both.

**Fix.** Reconcile against the notebook's output and correct the deck.

---

### [MINOR] F12 — Denominator error on lineage spread, and two remaining set-size ambiguities

**Evidence anchor:** `reports/dark_orf_clustering.md:179` — "The top family spans **224 of the 4,174
lineages (5.4%)**".

**What is wrong.** 4,174 is the number of distinct `mob_cluster` values across all 71,414 plasmids in
the working set; but only **63,993** of those plasmids actually carry a dark ORF, and they span
**4,004** clusters. A dark family can only occur on a carrier, so the correct denominator is 4,004 and
the figure is 5.6%. Verified:

```
distinct mob_cluster in the 71,414 subset : 4,174
dark-ORF carrier plasmids                 : 63,993
distinct mob_cluster among carriers       : 4,004   ← the deck uses this; the report uses the other
```

Both numbers are individually correct and both documents say "in this subset", meaning different
subsets. Two smaller instances of the same family of problem remain from the verification pack and
should be closed at the same time: the 143,503 / 143,590 analysis-set pair (§7.1 — 143,590 survives
only in `METADATA_ENRICHMENT_PLAN.md:95`, so this is nearly closed), and the 562 / 694 PlasAnn-AMR
definitional leak (§7.4), quoted in two reports from two different sources.

**Fix.** Correct the denominator at `dark_orf_clustering.md:179`; define "lineages in this subset" once,
against carriers, and use it everywhere; recompute the AMR leak from the master column (694) and
propagate.

---

### [MINOR] F13 — The foundational Arc-2 notebook has no narrative at all

**Evidence anchor:** `notebooks_my/cryptic_plasmids.ipynb` — 24 cells, **0 markdown cells**.

**What is wrong.** Every other notebook in `notebooks_my/` carries 6–9 markdown cells of genuinely
high-quality reasoning. This one — which establishes the 71.5% unannotatable figure, the 43.5%
dark-end-to-end figure, and the payload-free definition that all of Arc 2 rests on, and which is cited
as "upstream" by four other notebooks and two reports — has none. Its results are reconstructable only
from `reports/small_cryptic_methodology.md`, which describes a *different* set (F10). The known
Acinetobacter clonal-collapse lesson, referenced by name in three other notebooks as "the trap that
caught the Acinetobacter signal in cell 10", is documented nowhere in prose.

**Fix.** Add markdown cells matching the standard of the sibling notebooks, at minimum: the
payload-free definition and its known leak, the 71.5% / 43.5% / 54.2% distinction (the deck's Part 5
speaker note already contains the right explanation and can be lifted), and the cell-10 Acinetobacter
result the rest of the programme cites as its founding control.

---

## 5 · Scores (0–100)

| Dimension | Score | Justification |
|---|---:|---|
| **Internal validity** | **72** | Controls are real and repeatedly decisive — the positive control on the Pfam assay, the MOB-lineage clonal control, the three-way anti-defense adjudication — but the resource layer carries an undetected column defect (F2) under the founding question, and the shortlist's construction has a blind spot the project identified and did not act on (F5). |
| **External validity** | **58** | The collection is 62% IMG/PR metagenome-only with 46% host coverage and a lineage skew so severe that one MOB cluster holds 19,427 of 71,414 plasmids; the programme states these limits candidly everywhere but cannot escape them, and the Arc-2 claims are bounded by "invisible to Pfam-A 38.2", not by novelty. |
| **Reproducibility** | **68** | Every report names its scripts, tool versions are pinned, seeds are fixed, an independent recompute reproduced 46/52 numbers, and the preserved `dark30_*` outputs were deliberately protected — but there is **no version control** (F7), one published figure has no recorded expression (F9), and the manuscript's own "single deterministic script" guarantee cannot be checked against the code that produced it. |
| **Claim–evidence calibration** | **55** | Bimodal rather than uniformly poor: the deck and the newer Arc-2 documents are exemplary, while the older Arc-1 documents assert refuted claims as established (F1, F8) and the catalogue is materially under-claimed (F4). The programme doubts the wrong things. |
| **Contribution** | **74** | The 67-column master table and the Pfam-negative family catalogue are genuine, reusable objects, and the Pfam partition of a plasmid dark proteome against a proper positive control is, as far as this audit can tell, not otherwise available; the multireplicon and AMR arcs contribute methodological cautionary tales rather than findings. |
| **Documentation quality** | **63** | Individually, several documents are among the best internal methods writing I have reviewed — `pfam_dark_validation.md` §6 and the deck especially. Collectively they are not a corpus: the front door is blind to half the programme (F6), three documents request an analysis that has been run (F3), and no document states which claims are still live. |

---

## 6 · Genuine strengths

These are real and I would not want them lost in a revision.

1. **The self-correction is substantive, not performative.** `reports/pfam_dark_validation.md` §6.2
   documents a CDS-indexing bug that silently dropped 10.7% of dark ORFs into the "no Pfam match"
   bucket — i.e. **the bug inflated the project's own novelty claim**, and it is written up in full,
   with the corrupted cross-tab (59.2/26.8/8.7/5.3) printed beside the correct one
   (57.4/28.5/10.9/3.2). The verification pack confirms the fix holds with zero join loss. Publishing
   the number that made your own result look better is the strongest single credibility signal in this
   repository.

2. **The positive control on the Pfam assay is the right experiment.** Searching the 78,397
   PlasAnn-*named* proteins from the same plasmids alongside the dark representatives converts an
   uninterpretable "26.3% hit rate" into evidence, because the control hits at 76.7%. Most dark-matter
   surveys do not do this. It is also what makes the "28.5% is annotation failure" finding possible at
   all, and that finding is against the project's own interest.

3. **Three independent lines killed the anti-defense hypothesis, and two of them reproduced each
   other's false positives.** dbAPIS and AntiDefenseFinder independently called `AcrIIA21` on RepA_N
   replication initiators. Recognising that the concordance of *artefacts* is itself evidence — and
   reporting the adjudicated 158 as "a floor, not an estimate" — is careful work.

4. **The clustering-compliance audit was self-initiated.** Nobody made the project read MMseqs2's help
   text and discover that `--cluster-reassign` was missing; it then measured the damage
   (89.12% → 99.26% compliance), swept four thresholds against literature landmarks, and reported that
   its own shortlist collapses 93 → 5. That is the opposite of thesis protection.

5. **The deck is the most honest document in the repository** and should be the template for the
   corpus. It labels three negative results and two methodological corrections as deliberate content
   ("Those are the parts I most want this room to attack"), separates unit-of-analysis claims
   explicitly ("Every 'X% of the dark plasmidome' claim needs its unit stated"), and its closing
   ledger's right-hand column ("Nothing here is expression evidence. We have not shown a single one of
   these proteins is made") is a better Limitations section than either manuscript has.

6. **Definitional discipline where it was hardest.** The habitat taxonomy was locked before analysis;
   the standing rule to build One Health compartments from `hab_sub` + `is_clinical` rather than
   `hab_top` is enforced by importing the map from the Phase-3 script so notebooks cannot drift; the
   simulated-community exclusion is filtered on `hab_top` rather than `geo_source` for a stated reason.

---

## 7 · Disposition

**State of the work.** Two arcs at very different maturities, in one repository with no mechanism for
telling them apart. **Arc 1 is closed and should be declared closed** — its central claim is refuted
by its own analysis, its resource layer has an unfixed column defect under its founding question, and
its remaining value is a cautionary methods note about size confounding in habitat-stratified plasmid
comparisons. **Arc 2 is a genuine paper roughly one month of work from submittable** — the analysis is
substantially done, the key result is external-database-validated, and my recompute in F4 shows the
central deliverable is more robust than the project has claimed. What blocks it is not analysis. It is
that the repository currently cannot tell a reader which of its own numbers are alive.

**Nothing should be submitted anywhere until items 1 and 2 below are done.** They are both editing,
not science, and together they are perhaps two days.

### The three things that must happen next, ranked

**1. Reconcile the claim record across the corpus — before any other work.** (F1, F3, F6, F8)
Add dated supersession headers to `reports/amr_onehealth_methodology.md` and
`manuscript/amr_onehealth_NAR.md`; add forward pointers to `clustering_threshold_sensitivity.md` from
`dark_orf_clustering.md`, `pfam_dark_validation.md` and `dark_plasmidome.ipynb`, and rewrite the three
stale "untested at 40% and 50%" limitations; rebuild `PROJECT_OVERVIEW.md` to cover both arcs with an
explicit **Claim status** table listing every headline claim as live / superseded / withdrawn, and by
what. This is the single highest-value action in the programme, because every other strength here is
downstream of the reader being able to trust a document at face value. `git init` in the same pass.

**2. Fix the resource layer, and promote the catalogue to its corrected form.** (F2, F4, F7, F9)
Correct `master_table_data_dictionary.md:129`; add `pf_n_fam` to `harvest_typing.py` (the value is
already computed and thrown away) and rebuild; replace 24,583 with 22,422 or drop the PlasmidFinder
multireplicon figure in favour of the mob-based 18,685 per the locked rule; record or correct the
95,442 expression; regenerate `dark_novel_families.tsv` at 30%+`--cluster-reassign` as the canonical
3,490-family deliverable and fold the stability table from F4 into
`clustering_threshold_sensitivity.md`. Then deposit master table, catalogue and representative FASTA
with a licence and checksums.

**3. Write the Arc-2 catalogue paper to the abstract in §3 — and run the structure search first.**
Every Arc-2 document reaches the same wall: `pfam_dark_validation.md` Limitation 1, "the honest next
step is `foldseek` against the AlphaFold DB"; the deck's Part 13 "Structure first — it is cheap, and it
is the test we keep deferring". They are right, and the paper is meaningfully weaker without it,
because "invisible to Pfam-A 38.2 at its gathering thresholds" is a statement about a database whereas
"no structural match in AFDB/PDB" is a statement about the proteins. At a median 130 aa over 3,490
representatives this is days of compute. Do it, then write. Do **not** lead the paper with the 93- or
5-family shortlist: lead with the catalogue and the partition, and present the shortlist collapse
where it belongs — as the methods result about threshold-dependence that it actually is.
