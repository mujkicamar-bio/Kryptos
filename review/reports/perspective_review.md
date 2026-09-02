# Reviewer 3 — novel-protein-family / structural-bioinformatics perspective

*Scope: the dark-plasmidome protein-family work. Read `review/VERIFICATION_PACK.md` first; the
arithmetic is confirmed and I did not re-check totals. Everything below attacks definition, design
and interpretation, and every quantitative claim I make was computed by me from the files on disk.
Commands and outputs are in §Q and §S. I modified nothing; all my scratch output lives outside the
project tree.*

---

## Summary

The dark-plasmidome pipeline is unusually honest for this kind of work — it ships a positive
control, a threshold sweep, a clustering-compliance audit and a public defect log, all of which most
novel-family papers omit. But it establishes novelty against **one database at one precision-tuned
threshold**. In the FESNov / NMPFamsDB world, "no Pfam-A hit at `--cut_ga`" is the first of five or
six mandatory filters; this project published the first filter's output as "the genuine residue". I
tested that. Re-searching the same local Pfam-A at `E<1e-3` recovers 853 of the 3,343 "novel"
families (25.5%) — mostly relaxases, Rep initiators, RHH/HTH regulators and antitoxins — and
re-clustering under the *Nature* 2024 criterion the project itself cites merges 246 more onto
Pfam-named families. **28.4% of the catalogue, and 55% of the 93-family core, falls with no new
database.** Two further problems are decisive: the flagship "224 lineages" table names five
representatives that do not exist in the current output, and a sampling null — which I built,
because nobody had — inverts the recurrence argument. The ORF-call worry I largely cleared: three
independent tests exonerate the calls, with a bounded ~15–20% residue.

---

## Findings

### [CRITICAL] F1 — Novelty is established against Pfam alone; the field requires a homology cascade, and I can show what it would remove

**Evidence anchor:** `reports/pfam_dark_validation.md` §3.2 (`hmmsearch --cut_ga` vs Pfam-A 38.2, the
only external search); §4.2 cross-tab "57.4% — the genuine residue"; §7 limitation 1.
`data/pfam_run/dark_novel_families.tsv` (3,343 families) is the deliverable derived from that single
test.

**What is wrong.** The two catalogues this project is implicitly competing with do not accept a
single-database negative:

| | FESNov (Rodríguez del Río et al., *Nature* 2024) | NMPFamsDB (*NAR* 2024) | this project |
|---|---|---|---|
| sequence DB screen | UniRef / UniProtKB **and** eggNOG orthologous groups | UniRef, Pfam, PDB | Pfam-A only |
| domain/profile screen | Pfam **plus** profile-level (HMM-vs-HMM) search | profile–profile (HHsearch-class) on the SSN | Pfam sequence-vs-profile only |
| structure screen | AlphaFold/ESMFold + structural comparison on candidates | structural comparison on candidates | **none** |
| ORF/quality filter | explicit gene-call and length quality control | ≥100-member families only | **none** |
| output of "novel" | family + HMM + MSA + accession | family + HMM + MSA + accession | family id + counts |

Concretely missing, each of which is standard and cheap relative to what has already been spent:
UniRef50/UniRef90 or UniProtKB via MMseqs2; eggNOG/COG/KEGG orthology assignment (eggNOG-mapper);
NCBI nr; profile–profile remote homology (HHblits/HHsearch against Pfam + UniRef30, or MMseqs2
profile search); and structure-based search (Foldseek vs AFDB / PDB / ESM Atlas). Pfam-A covers
~30k families and, crucially, **`--cut_ga` is tuned for precision, not recall** — it is the wrong
operating point for a novelty claim, where you want maximum recall on the negative side. The
project treats a precision-tuned threshold as if it bounded detectable homology.

**What I computed.** Two independent demonstrations that the cascade would bite, both run without
any new database:

1. **Clustering alone already de-novels 7.4%.** Re-clustering the 92,752 dark family
   representatives at 30% identity under FESNov's own coverage rule (`-c 0.5 --cov-mode 1`) puts
   **246 of the 3,343 "novel" families (7.4%) into a group that also contains a Pfam-**named** dark
   family** (F4). Those are, by a published standard, homologs of things Pfam names. No new data
   was needed to find this — only the criterion the project cites in its own literature table.
2. **Relaxing only the threshold, not the database, de-novels a quarter of the catalogue.**
   I re-searched all 3,343 deliverable representatives against the *same local* Pfam-A 38.2, same
   HMMER, changing only `--cut_ga` to `-E 1e-3`:

   ```
   hmmsearch --noali -E 1e-3 --domE 1e-3 --cpu 24 --domtblout relaxed3343.domtbl \
             data/refs/pfam/Pfam-A.hmm data/flashfold_run/targets_top3343.faa
   ```

   ```
   queries 3343 -> 853 with a hit at full-sequence E < 1e-3   (25.5%)
                   444 at E < 1e-4  (13.3%)
                   226 at E < 1e-5  ( 6.8%)
   recovered representatives median 138 aa vs 102 aa for the non-recovered
   most-recovered Pfam families: HTH_17 (29), NikA-like (22), HTH_23 (18), RHH_1 (18),
     DUF3847 (17), Rep_1 (15), HTH_36 (13), HTH_11 (13), PhdYeFM_antitox (12),
     RepB-RCR_reg (10), CopG_antitoxin (9), Relaxase (9)
   ```

   These are *the plasmid backbone again* — relaxases, Rep initiators, ribbon-helix-helix and HTH
   regulators, antitoxins. The project's §4.3 finding that the Pfam-rescued half is mostly
   machinery PlasAnn already has categories for **extends into the supposedly novel residue**.

**The combined effect, with no new database at all:**

```
of the 3,343 'novel' families:
  relaxed Pfam (same DB, E<1e-3)           :  853  (25.5%)
  merges with a Pfam-named family (FESNov) :  246  ( 7.4%)
  UNION                                    :  948  (28.4%)   -> 2,395 survivors
within core_dark_families.tsv (93)         :   51  (54.8%)   -> 42 survivors
```

**28.4% of the novelty catalogue, and 55% of the 93-family core, falls to two operations that use
only files already on this disk.** That is before UniRef, eggNOG, nr, profile–profile or structure —
i.e. before the actual field-standard cascade begins.

**Concrete fix.** Do not publish 57.4% / 3,343 as novelty. Run, in this order, and report the
survivor count after each step as a funnel (this is exactly how FESNov presents its catalogue):
Pfam `--cut_ga` → Pfam `E<1e-3` → MMseqs2 vs UniRef50 (`-s 7.5`) → eggNOG-mapper → MMseqs2 profile
or HHblits vs UniRef30 → Foldseek. Until then the defensible phrasing is the one the report already
uses once and then abandons: *"invisible to Pfam-A 38.2 at its gathering thresholds"* — never
"genuinely uncharacterized", never "novel families".

---

### [CRITICAL] F2 — The flagship lineage-spread table cites five family representatives that do not exist in the output

**Evidence anchor:** `reports/dark_orf_clustering.md` §R2, table "Most lineage-widespread families",
rows `GenBank_CP048556.1|10` (1,089 proteins / 1,079 plasmids / **224 lineages** / 26 habitats),
`IMGPR_plasmid_3300038695_000007|9`, `IMGPR_plasmid_3300022496_000003|3`,
`IMGPR_plasmid_3300047678_000089|1`, `IMGPR_plasmid_3300007123_000008|1`. The same 224 figure carries
the report's Interpretation section ("conservation across 224 lineages is not the signature of
neutral junk"), and it is the sentence this project's argument rests on.

**What is wrong.** None of those five ids is a cluster representative in the shipped clustering.
Each is an ordinary *member* of a different family. Checked with
`awk -F'\t' -v r="<rep>" '$1==r{c++}END{print c+0}' data/dark_orf_run/dark30_cluster.tsv` (0 for all five)
and the reciprocal `$2==r` lookup:

```
GenBank_CP048556.1|10             -> member of IMGPR_plasmid_3300024992_000013|1
IMGPR_plasmid_3300038695_000007|9 -> member of IMGPR_plasmid_3300014621_000006|4
IMGPR_plasmid_3300022496_000003|3 -> member of IMGPR_plasmid_2568526604_000001|8
IMGPR_plasmid_3300047678_000089|1 -> member of IMGPR_plasmid_3300037153_000622|11
IMGPR_plasmid_3300007123_000008|1 -> member of IMGPR_plasmid_3300007271_000003|3
```

`darkfam_top200.tsv` — the file the report says produced that table — leads with
`IMGPR_plasmid_3300047172_000056|15` (1,123 / 1,115 / **234** / 31), a different family. The
published table is pre-fix output that was never regenerated after the 2026-08-28 dedup and
CDS-reindex corrections. The corrections notice at the bottom of the report claims "every number
above is post-fix"; §R2 is not. The verification pack checked family *counts* (§4) but not these
rows, so this survived.

**Concrete fix.** Regenerate §R2 from `darkfam_top200.tsv`, and grep the whole repo, the deck and
the memory notes for `224` and for those five accessions. Then re-read F3 before deciding what the
corrected number is worth.

---

### [CRITICAL] F3 — The recurrence-across-lineages argument inverts under the correct null, which nobody built

**Evidence anchor:** `reports/dark_orf_clustering.md` §R2 and Interpretation ("neutral junk does not
stay conserved across 224 lineages"); `notebooks_my/dark_families_90pct.ipynb` ("Sequence
conservation that strong across that many lineages is the signature of purifying selection");
limitation 3 in both reports ("`mob_cluster` is a proxy, not a phylogeny") acknowledges the concern
but no test follows.

**What is wrong.** "N lineages" has no meaning without the sampling expectation. The right first
null is: *if the same number of carrier plasmids had been drawn at random from the corpus, how many
mob_clusters would they have touched?* That is analytic (occupancy under sampling without
replacement) and takes seconds. I computed it (§Q4). Every result runs the wrong way:

| family set | median plasmids | median observed lineages | median **expected** | obs/exp |
|---|---:|---:|---:|---:|
| all families ≥10 members | — | 4 | 13.5 | **0.24** |
| `dark_novel_families.tsv` (3,343) | — | 3 | 12.8 | **0.20** |
| `core_dark_families.tsv` (93) | — | 22 | 114.5 | **0.18** |
| top family `IMGPR_..._3300022589_000018\|2` | 2,468 | 88 | 892.4 | **0.10** |
| the five 90%-identity survivors | 101–227 | 14–42 | 67–139 | **0.12–0.47** |

Not one family in the 93-family core exceeds its random expectation; **0.4% of the 3,343 do.** A
random draw of 1,115 plasmids from this corpus touches ~515 mob_clusters, so the corrected flagship
figure (234 lineages) is a **2.2-fold deficit**, not a surfeit. Real plasmid gene families are
lineage-*restricted* relative to random sampling, which is exactly what vertical inheritance plus a
concentrated corpus predicts (the largest mob_cluster alone holds 29.3% of the 63,993 carriers;
Simpson index 0.087).

**The one version of the argument that survives is the calibrated one, and it points the other
way.** Using Pfam-named dark families as a "known real gene" reference:

- Pfam-**positive** families (≥10 members): median obs/exp = **0.312**
- Pfam-**negative** families (the deliverable): median obs/exp = **0.197** (Mann-Whitney p = 2.2e-72)
- individual known backbone families for scale: `Rep_1` 0.228, `Rep3_N` 0.289, `Mob_Pre` 0.213,
  `MobA_MobL` 0.308, `HTH_3` 0.397

So the candidate-novel families are *significantly less* lineage-spread than the known plasmid
backbone. Whatever else is true of them, they are not distinguished from junk by travelling further.

**Concrete fix.** (i) Delete every "not the signature of neutral junk" inference. (ii) Report
obs/exp, not raw lineage counts, everywhere lineage spread appears — including in the shortlist
criterion, which currently ranks on a sampling-dominated statistic. (iii) For a real selection
claim, use the right instrument: dN/dS on within-family codon alignments. The nucleotide sequence is
already on disk (`data/plasann_run/gbk/`, `ORIGIN` block — I used it in §Q), so this is a day's
work, and for the five 90%-identity survivors it is the *only* argument that would carry. "≥90%
amino-acid identity across 14–42 lineages" is compatible with purifying selection and equally
compatible with recent transfer; dN/dS separates them and lineage counting cannot.

---

### [MAJOR] F4 — The coverage axis is the dominant clustering axis and it was deliberately not swept; "stricter" is not "conservative" for a novelty catalogue

**Evidence anchor:** `reports/clustering_threshold_sensitivity.md` §2 ("The coverage setting was
already stricter than any of them … So the sweep is on identity, with coverage held at 80%") and
§6 ("One coverage setting … because nothing in the literature reviewed uses a stricter value");
`plans/conservative_reclustering.md` D3. Contradicted inside the project by
`notebooks_my/zoom_two_dark_orfs.ipynb` §6.

**What is wrong.** The project's own zoom notebook already proves the point and the sweep report
does not absorb it: with `-c 0.8 --cov-mode 0`, a 208 aa and a 155 aa protein *can never* share a
family at any identity (0.75 < 0.80), so the clusterings across the sweep are **not nested** and
ORF-A is a 20-member family at 30% and a 400-member family at 90%. That is a coverage artifact, and
it means the sweep varied the axis that mattered less.

More importantly, the framing is backwards. Strict bidirectional coverage is conservative when you
are *asserting* homology (UniRef's use case: don't merge fragments). It is anti-conservative when
you are *counting novel families*, because fragmentation splits one family into several and every
split adds a row to the novelty catalogue. I measured the size of this effect:

- re-clustering the 92,752 representatives at the project's own criterion merges only 73 (0.08%) —
  a good self-consistency check;
- re-clustering them at FESNov's criterion (`-c 0.5 --cov-mode 1`, 30% id) yields **72,335
  clusters — a 22.0% reduction**;
- the deliverable itself: **3,343 → 2,989 families (−10.6%)**; the 93-family core → 88.

So "3,343 novel families" is ~11% inflated relative to the *Nature* 2024 standard the project cites
as looser than its own, and the whole-catalogue family count is ~22% inflated.

**Second, on whether a 30% MMseqs2 `easy-cluster` family is a "protein family" at all.** In this
field, no — not on its own. `easy-cluster` is greedy set-cover around a single representative, so
membership is a star topology anchored on one sequence, and the *pairwise* relationship between two
members is never tested. At 30% identity that is below the twilight zone, where sequence identity is
not a reliable homology signal, and the project's own compliance audit shows why this matters:
`reports/clustering_threshold_sensitivity.md` §3 measured **89.12% compliance** on the published
clustering — one member in nine did not meet the stated criterion against its own representative,
almost always on coverage. NMPFamsDB builds a sequence-similarity network and clusters with HipMCL;
FESNov validates with profile-level search. Both are asserting something about the *family*, not
about distances to a hub.

**Concrete fix.** (i) Report the deliverable under two criteria — the project's and FESNov's — and
lead with the lower count. (ii) Move from identity clustering to profile clustering for the
deliverable only (3,343 families is tiny): build a per-family MSA, then MMseqs2 profile search or
HHsearch all-vs-all, and merge families whose profiles match. This does double duty — it is also
step 5 of the F1 cascade and it is what makes the catalogue citable (F5). (iii) Retire the
non-`--cluster-reassign` clustering. `reports/dark_orf_clustering.md`,
`reports/pfam_dark_validation.md` and three notebooks still stand on the 89.1%-compliant run; the
sweep report's own numbers say `--cluster-reassign` alone moves 12.00% of members and takes the
core set 93 → 61. A manuscript cannot cite a clustering whose stated criterion it is known to
violate for 11% of members.

---

### [MAJOR] F5 — The deliverable is not a resource another lab can consume

**Evidence anchor:** `data/pfam_run/dark_novel_families.tsv` — six columns,
`rep / proteins / plasmids / lineages / habitats / pfam_hit`, keyed on
`IMGPR_plasmid_3300022589_000018|2`-style ids. `data/dark_orf_run/darkfam_stats.tsv` — seven
columns, same key. `data/dark_orf_run/core_dark_families.tsv` (93),
`recluster/core_dark_families_90.tsv` (5).

**What is wrong.** Against what NMPFamsDB and FESNov ship, this is a table of counts, not a
catalogue. Missing, in descending order of importance:

1. **Sequences.** `dark_novel_families.tsv` contains no sequence at all. The representatives exist
   in `data/flashfold_run/targets_top3343.faa` under *different* identifiers (`fam00001`…) which
   appear nowhere else, so the two files cannot be joined without the intermediate
   `targets_top3343.tsv`. Nobody can use the deliverable as shipped.
2. **HMM profiles and MSAs.** These are the unit of exchange in this field. Without a profile,
   another lab cannot ask "is my protein in your family?" — which is the entire point of publishing
   a family catalogue.
3. **A stable accession scheme.** `IMGPR_plasmid_3300022589_000018|2` is a coordinate into *this*
   run: it changes when the clustering changes (F2 is exactly this failure mode). Families need
   opaque, versioned accessions (`PDF00001`, plus a catalogue version) with representative ids kept
   as a separate column.
4. **Per-family confidence.** One boolean (`pfam_hit`) is not confidence. The field expects, per
   family: member count, aligned-fraction / MSA depth and coverage, internal identity distribution,
   the propagation-concordance figure (`pfam_dark_validation.md` §4.5 measures 96.1% overall and
   23/200 families internally mixed — that should be *per family*, not a global caveat), and, once
   run, mean pLDDT and best Foldseek hit.
5. **Taxonomic and ecological metadata.** `lineages` and `habitats` are integers. Users need the
   host taxon distribution (LCA), the habitat vector, and — given F3 — obs/exp rather than raw
   counts. The columns exist upstream in `plasmid_metadata_master.tsv`; they are simply not joined.
6. **Length/quality columns.** `rep_len_aa` is in `darkfam_stats.tsv` but dropped from
   `dark_novel_families.tsv`, which is the file a reader is pointed at. Given §Q, length is the
   single most decision-relevant column and it is the one missing.
7. **A submission target.** None is named anywhere. Realistic options: deposit profiles + MSAs in
   Zenodo with a DOI and offer them to InterPro as candidate unintegrated signatures; or submit
   families to Pfam as new entries (the project's own top rescued families are Pfam entries, so the
   route is familiar); or ask NMPFamsDB to ingest.

**Concrete fix.** Ship one directory: `families.tsv` (accession, version, representative id, length,
n_members, n_plasmids, obs/exp lineage, habitat vector, host LCA, propagation concordance, all
external-search results as separate columns), `representatives.faa` keyed on accession, `msa/`,
`hmm/` (one pressed concatenation), and a `README` stating the exact criterion. That is a day's work
and it is the difference between a table and a resource.

---

### [MAJOR] F6 — Representative-only search is used to label 378,552 proteins, and partial-domain rescues are counted as full

**Evidence anchor:** `reports/pfam_dark_validation.md` §3.1 (`DARKREP` = one representative per
family), §4.1 ("Weighted by family size, 149,160 of 378,552 dark ORFs (39.4%)"), §4.5 (96.1%
propagation, 23/200 families internally mixed), §7 limitation 2 ("Median Pfam alignment covers ~52%
of the dark protein").

**What is wrong.** Two compounding approximations. (i) Family-level Pfam status is imputed from one
sequence and then multiplied by family size, so ~4% of 378,552 ≈ 15,000 ORFs are mislabelled and the
error is not random — it concentrates in exactly the large, heterogeneous, 30%-identity families the
novelty claim depends on. The `VALID` check is good design, but a check is not a fix when the
searched set is only 92,752 sequences — three-quarters the size of the `NAMED` set that *was*
searched exhaustively. (ii) A hit covering 52% of the query is scored as "Pfam names it outright"
in the §4.2 cross-tab. In this field a half-covered protein is a **partial** assignment: a known
domain plus an uncharacterized region, which is a different and often more interesting object than
either "known" or "novel".

**Concrete fix.** Search all 378,552 (the earlier run did 174,827 in ~6 minutes on 88×4 cores — this
is ~2× that, i.e. minutes, not a reason to approximate). Then replace the 2×2 cross-tab with a
three-way split: fully-covered hit / partial hit (report the uncovered residue length) / no hit. The
"partial" class is a publishable finding in its own right — novel domains attached to known plasmid
backbone proteins.

---

### [MINOR] F7 — The structure test is not deferred, it is blocked, and the tooling is absent

**Evidence anchor:** `data/flashfold_run/` contains `targets_top250.faa`, `targets_top3343.faa`,
their TSVs and `install.log`. Nothing folded. `install.log` is a conda failure log:
`InvalidArchiveError … [Errno 122] Disk quota exceeded` repeated across `ncurses`, `libgfortran`,
`liblapack`, `pdbfixer`, `openblas`, `python_abi`. I also checked the machine: **no `foldseek`, no
`hhblits`/`hhsearch`, no `colabfold_batch` anywhere** on the system or in any of the six conda envs;
only `mmseqs` and `diamond` (both in `panaroo`).

**What is wrong.** Five separate documents name Foldseek-vs-AlphaFold-DB as "the honest next step"
(`pfam_dark_validation.md` §7.1, `dark_plasmidome.ipynb`, `dark_families.ipynb`, `antidefense.ipynb`,
`zoom_two_dark_orfs.ipynb`). The reader is left believing it is a scheduling choice. It is a home-
directory quota failure that nobody recorded as a blocker, and Foldseek — the tool actually named —
was never installed at all.

**Concrete fix.** State the blocker in the limitations, install into `/scratch` (6.1 TB free) rather
than `~/.conda`, and see §S.

---

### [MINOR] F8 — Two orphaned framing claims worth one sentence each

**Evidence anchor / fix:**
- `reports/clustering_threshold_sensitivity.md` §2 lists FESNov as "MMseqs2 30% / 50% `--cov-mode 1`
  / min size ≥3" and concludes the project is stricter. That table compares only the *clustering*
  parameters and omits FESNov's novelty cascade entirely, which is what makes FESNov a novel-family
  catalogue rather than a clustering. As written, the table reads as evidence that this project's
  methods exceed the *Nature* 2024 standard. They do not; they are a subset of it. Add the missing
  columns (as in F1's table).
- `notebooks_my/dark_families_90pct.ipynb` presents the five 90%-identity survivors as the strongest
  per-candidate evidence, on 70/49/121/105/78 aa representatives. Given §Q, the shortlist that
  carries a "these are doing something" claim is the one most exposed to the short-ORF and the
  sampling problems simultaneously. It needs both the dN/dS test (F3) and the per-candidate
  coding-potential numbers in §Q8 attached to it before it is presented that way. In fairness:
  two of the five appear in the 3,343 set and **neither was recovered by the relaxed Pfam search**,
  so they survive Tier 0 intact — they are the best candidates in the project, and they deserve the
  two tests (dN/dS, independent gene calling) that would actually settle them.

---

## Q — Quantitative section: the ORF-call problem (point 2)

**Verdict up front, because it cuts against the reviewer's prior: the length skew is real and large,
but the ORF calls are mostly sound. Three independent tests exonerate them, and the residual risk is
bounded at roughly 15–20% of the Pfam-negative dark set, concentrated below 100 aa.**

Interpreter throughout: `/gorilla/home/amujkic/.conda/envs/genesis/bin/python`. Scripts written to
scratch, not to the project tree.

### Q1 — The raw length distributions

```python
# lengths parsed from data/dark_orf_run/dark_orfs.faa and all_cds.faa (|N suffix = PlasAnn-named)
```

```
all dark ORFs              n= 378552  min= 29  q25= 86  med=138.0  q75=249  <50aa= 3.94%  <80= 20.52%  <100aa= 34.59%
NAMED CDS (control)        n=  78397  min= 29  q25=100  med=185.0  q75=312  <50aa= 1.05%  <80= 12.31%  <100aa= 24.99%
```

Both sets floor at **29 aa**, i.e. a ~90 nt caller minimum. The dark histogram *rises* from that
floor to a mode at 80–89 aa rather than decaying from it — random-sequence ORFs would decay
monotonically from the floor. So the dark set as a whole is not a spurious-ORF pile. The report's
"only 3.9% under 50 aa" is true; it is also the wrong statistic, because it is computed on the whole
dark set rather than on the Pfam-negative residue that carries the novelty claim.

### Q2 — The decisive split: Pfam-negative dark families are half the length of Pfam-positive ones

```python
st = pd.read_csv("data/dark_orf_run/darkfam_stats.tsv", sep="\t")      # rep_len_aa
pf = pd.read_csv("data/pfam_run/darkfam_pfam.tsv",  sep="\t")          # pfam_hit
d  = st.merge(pf[["rep","pfam_hit"]], on="rep", validate="1:1")        # 92,752 rows
```

```
DARKREP Pfam-POSITIVE (rep len)   n= 24357  med=225.0  q25=129  q75=367  <50aa= 0.77%  <100aa= 14.94%
DARKREP Pfam-NEGATIVE (rep len)   n= 68395  med=102.0  q25= 69  q75=166  <50aa= 8.37%  <100aa= 48.75%
   Mann-Whitney p = 0 ; median ratio neg/pos = 0.453

size-weighted (per dark ORF):
  Pfam-POS ORFs    n=149160  med=212.0   <50aa=0.25%   <100aa=18.45%
  Pfam-NEG ORFs    n=229392  med=108.0   <50aa=6.28%   <100aa=45.26%

the deliverable:
  dark_novel_families.tsv (3,343)  med=109.0  <50aa= 5.56%  <100aa= 44.48%
  fam>=10 Pfam-POS (2,280)         med=226.0  <50aa= 0.26%  <100aa= 15.96%
```

**The skew the reviewer predicted is present and it is a factor of 2.2 in median length.** Nearly
half the deliverable's representatives are under 100 aa. This is diagnostic and it is not reported
anywhere in the project: `dark_orf_clustering.md` §R1 reassures the reader with the whole-set
median (138 aa) and never splits by Pfam status.

### Q3 — But how much of that is HMM detection power? A control, and its limits

Same test on the `NAMED` set, where every protein is known to be real:

```
NAMED Pfam-POSITIVE   n= 60124  med=214.0
NAMED Pfam-NEGATIVE   n= 18273  med=143.0     median ratio neg/pos = 0.668
```

Pfam hit rate as a function of length:

```
len bin      DARKREP n  DARK hit%    NAMED n  NAMED hit%
29-49             5913       3.2%        825      21.0%
50-59             5540       6.6%       1223      57.9%
60-79            13882       9.9%       7603      73.8%
80-99            11644      14.7%       9944      66.8%
100-149          18591      21.4%      13956      75.2%
300-499           9478      60.9%      17766      92.1%
```

Detection power is real (0.668 vs 0.453 in the ratio), but it does not explain the gap: Pfam finds
58–74% of *real* 50–79 aa plasmid proteins and only 7–10% of dark ones in the same band. **Caveat I
must state: the `NAMED` control is biased.** PlasAnn named those proteins by homology, so short
NAMED proteins are pre-selected for database detectability. The control therefore over-states
achievable power. It bounds the effect rather than measuring it — which is why I went to
composition and coding potential instead.

### Q4 (note) — the lineage-spread null is in F3; the commands are:

```python
# occupancy under sampling without replacement, mob_cluster sizes from plasmid_metadata_master.tsv
# restricted to the 63,993 dark-ORF carrier plasmids, 4,004 distinct clusters
def E_lin(n):  # expected distinct clusters when n carriers are drawn from N
    a = gammaln(N-vc+1) - gammaln(N-vc-n+1); b = gammaln(N+1) - gammaln(N-n+1)
    return float((1-np.exp(a-b)).sum())
```

### Q5 — Amino-acid composition against a random-ORF null

Expected composition of a translated random ORF at GC = 43.6% (the small-cryptic median), from the
standard code; total-variation distance from it:

```
                                        TVdist_to_randomORF
NAMED CDS (real proteins, control)                   0.1670
all dark ORFs                                        0.1547
dark reps, Pfam-POSITIVE                             0.1447
dark reps, Pfam-NEGATIVE                             0.1381
dark reps Pfam-NEG, <100 aa                          0.1388
dark reps Pfam-NEG, <50 aa                           0.1498
dark_novel_families.tsv 3,343 reps                    0.1544

Cys:  null 0.0326 | NAMED 0.0086 | darkNEG 0.0128 | darkNEG<100 0.0145
Ser:  null 0.0978 | NAMED 0.0581 | darkNEG 0.0663
Glu:  null 0.0326 | NAMED 0.0787 | darkNEG 0.0679
```

The dark Pfam-negative set carries every selection signature the real proteins carry — 2.5-fold Cys
depletion, Ser depletion, 2-fold Glu enrichment — so it is not translated random DNA in bulk. It
sits marginally *closer* to the null than NAMED does, consistent with a minority admixture. I tried
to size that admixture with a two-component fit (`obs = (1-f)·NAMED + f·null`, length-matched) and
**the control failed**: Pfam-*positive* dark ORFs at 50–79 aa returned f = 33.7% where it should
have returned ~0. Composition cannot resolve this; I am reporting the failure rather than the
number. That is itself the argument for §S.

### Q6 — Genomic context: overlap and strand, all 456,949 CDS

Parsed `Start/End/Strand/Category` from all 1,072 `data/plasann_run/annot/*shard_*.tsv.gz`,
sorted within plasmid, computed neighbour overlaps.

```
dark (Open reading frame)  n=378552  ovl>30bp= 4.15%  ovl>60bp=0.90%  opposite-strand to prev=24.4%  antisense-ovl>30bp=1.40%
named CDS                  n= 78397  ovl>30bp= 7.05%  ovl>60bp=0.97%  opposite-strand to prev=26.7%  antisense-ovl>30bp=1.15%
  dark 0-50 aa             n= 14905  ovl>30bp= 4.92%                  opposite-strand=25.1%          antisense-ovl=1.95%
 named 0-50 aa             n=   825  ovl>30bp= 3.88%                  opposite-strand=28.0%          antisense-ovl=0.24%
```

**Dark ORFs overlap their neighbours *less* than named CDS do (4.15% vs 7.05%) and are slightly
*more* often co-oriented with them (24.4% vs 26.7% opposite-strand).** Antisense "shadow ORF"
overlap — the classic spurious-call signature — is 1.40% vs 1.15%, a negligible excess. This test
finds no artifact.

### Q7 — Coding potential vs an uncalled-ORF null, the decisive test

I parsed nucleotide sequence from the `ORIGIN` blocks of 30 randomly chosen GBK tarballs (1,634
plasmids from the 71,414 set, seed 20260902), built a codon log-likelihood-ratio model **from named
CDS only** (null = mononucleotide composition, **stop codons excluded from both model and scoring**
so the ORF-definition circularity is removed), and scored four classes. The negative control is the
right one: **66,013 ORFs ≥29 aa in all six frames of the same plasmids that PlasAnn did *not* call**
(≤30 bp same-strand overlap with any annotated CDS).

```
class                        n    meanS     medS      q10     %S>0    SD+ %
named                     2312    0.173    0.176    0.060    94.8%    18.1%
dark_pfam+                3559    0.143    0.147    0.025    93.1%    21.1%
dark_pfam-                5633    0.122    0.125   -0.044    83.0%    21.9%
uncalled                 66013   -0.136   -0.139   -0.313    16.5%     1.3%
```

The separation from the spurious null is enormous: **83.0% of Pfam-negative dark ORFs score above
zero versus 16.5% of uncalled ORFs from the same replicons.** Shine-Dalgarno-like motif in the 20 nt
upstream: 21.9% for dark Pfam-negative versus 1.3% for uncalled (random 16-mers from the same
upstream pool: 22.3% under a looser motif). *Caveat: both tests are partly circular — Prodigal-class
callers score coding statistics and RBS motifs when choosing starts, so called ORFs are enriched for
both by construction. The comparison bounds the worst case; it does not prove independence.*

The bound, using the length-matched 90th percentile of the uncalled distribution as a threshold:

```
  29-50 aa:  uncalled P90 S=0.047 -> 32.6% of dark_pfam- below it (named control: 41.7%, n=24)
  50-80 aa:  uncalled P90 S=0.022 -> 24.5% of dark_pfam- below it (named control:  8.4%)
  80-100 aa: uncalled P90 S=0.015 -> 21.9% of dark_pfam- below it (named control:  7.5%)
  100+ aa:   uncalled P90 S=0.086 -> 37.4% of dark_pfam- below it (named control: 15.2%)
```

**Excess over the known-real control: 14–22 percentage points.** That is my quantitative answer to
point 2: on the order of **one in six of the Pfam-negative dark ORFs has coding-potential
indistinguishable from the uncalled-ORF background**, versus one in twelve for proteins known to be
real. Not the wholesale artifact a naive reading of the length skew suggests; not negligible either,
and it lands squarely on the sub-100 aa half of the deliverable.

### Q8 — The actual deliverable representatives, individually

Targeted extraction of nucleotide context for the 812 plasmids carrying the 93-family core, the five
90% survivors and the 250 FlashFold targets; same scoring, null rebuilt from the uncalled ORFs of
those same plasmids (61,753 of them; a looser model, trained on 1,122 named CDS, so treat this as
indicative):

```
uncalled-ORF null on the same 812 plasmids: mean S=-0.103  P90=0.115  P99=0.292
  core_dark_families.tsv (93)     medlen=130  mean S=0.116  %S>null_P90=49.5%  SD+=24.7%
  core_dark_families_90 (5)       medlen= 78  mean S=0.150  %S>null_P90=60.0%  SD+=20.0%
  flashfold targets_top250        medlen=124  mean S=0.119  %S>null_P90=51.6%  SD+=22.8%

the five 90%-identity survivors:
  IMGPR_plasmid_2903362858_000001|1    78 aa   S= 0.283   SD+ no
  IMGPR_plasmid_2563366806_000011|1    49 aa   S= 0.210   SD+ no
  COMPASS_KP718939.1|8                 70 aa   S= 0.157   SD+ no
  IMGPR_plasmid_2551306362_000001|6   105 aa   S= 0.054   SD+ no
  IMGPR_plasmid_2700989144_000001|10  121 aa   S= 0.046   SD+ yes
```

Two of the five shortlisted candidates score barely above their own plasmids' uncalled-ORF mean.
Half the 93-family core sits below the null's 90th percentile. **The shortlist that carries the
strongest per-candidate claim in the project is the one that most needs an independent gene call.**

### Q9 — What is still not tested, and how to close it

Not done here, and each is a day or less:

1. **Independent gene calling.** Prodigal-2.6.3 (`-p meta` for the IMG-PR contigs, `-c` for closed
   circular) and/or Pyrodigal on all 71,414 replicons, then report for the deliverable: % of
   representatives recovered with identical start and stop, % with a shifted start, % not called at
   all. This is limitation #4 in both reports and it is the cheapest thing on the list. A family
   whose representative Prodigal does not call should be flagged in the catalogue, not silently
   shipped.
2. **Per-plasmid codon usage against the host.** I used a pooled model. A per-replicon model
   (named CDS of the same plasmid, falling back to the host genome's codon table via
   `mob_host_range`) is the sharper version and would also flag recently-acquired genes.
3. **Provenance stratification.** Dark ORFs are 78% IMG-PR-derived (296,329 of 378,552) and
   `small_cryptic_methodology.md` §7 shows provenance is a genuine confounder. I found only a mild
   length difference by source (IMG-PR median 141 aa vs GenBank 127 aa, 33.5% vs 38.6% under
   100 aa), which is reassuring, but the coding-potential test should be re-run stratified by
   provenance before the manuscript claims the metagenomic material behaves like the isolate
   material.

---

## S — Protocol: the structure test (point 5), sized

**Is folding the right next move? No — not first.** Folding 3,343 proteins to discover that some of
them are UniRef50 members would be an expensive way to skip a cheap step. And Foldseek against
AlphaFold-DB is *not* an alternative to folding: AFDB is a structure database, so you need predicted
structures for your queries either way (Foldseek's `--search-type` accepts sequence input against
structure DBs via 3Di prediction, but for 100 aa novel proteins that is weaker than folding first).
The correct order is: exhaust sequence, then fold what survives, then search structures. Each tier
below is gated on the previous one.

**Environment prerequisite (F7).** Install into `/scratch` (6.1 TB free), never `~/.conda` — the
existing `install.log` failure is a home quota error. `foldseek`, `hhsuite` and either
`fair-esm`/ESMFold or `colabfold` are all absent from this system today.

### Tier 0 — finish the job with what is already on disk (2 h, CPU, 48 cores)

- Re-run `hmmsearch` on `targets_top3343.faa` vs the local `data/refs/pfam/Pfam-A.hmm` at
  `-E 1e-3` (no `--cut_ga`) — **I ran this; 853/3,343 = 25.5% recovered, F1.**
- Re-cluster the 3,343 under FESNov's criterion and against the Pfam-positive families — **I ran
  this; −10.6% family count and 246 further families de-novelled, F4.**
- Build per-family MSAs (MAFFT or `mmseqs result2msa`) and run all-vs-all HHsearch on the 3,343
  profiles. This merges remote families *and* produces the MSAs the catalogue needs (F5).

**Proves:** how much of the residue is Pfam-detectable below the gathering threshold, and how much
is internal duplication. **Fails to prove:** anything about databases other than Pfam.
**Measured attrition: 28.4% (948 of 3,343), before HHsearch.** Tier 0 is therefore already done
except for the profile step, and it has removed more than a quarter of the catalogue.

### Tier 1 — sequence exhaustion (1 day, CPU)

- MMseqs2 `easy-search -s 7.5 -e 1e-3` against **UniRef50** (~65 GB indexed) and, for the
  survivors, **UniRef90**; then **NCBI nr** via DIAMOND `--ultra-sensitive` (already installed).
- **eggNOG-mapper** (diamond mode) against eggNOG 6 for orthologous-group and COG/KEGG assignment.
- MMseqs2 profile search (Tier 0's profiles) against UniRef30 — remote homology at profile level.

**Proves:** whether these families are unknown to the protein universe or only to plasmid tooling.
**Fails to prove:** anything about fold, and nothing about the sub-50 aa material, where all
sequence methods are underpowered. **Expected attrition: this is where the largest cut happens.**
Given that Tier 0 already removed 28.4% using Pfam alone, and that the 2,395 survivors have a median
representative length of ~102 aa on 78%-metagenomic carriers, my estimate is that **a further
30–50% of the survivors will acquire a UniRef50 or eggNOG assignment**. The plausible landing zone
for a defensible novel-family count is therefore **~1,200–1,700, not 3,343** — and that is the set
worth folding.

### Tier 2 — fold the survivors, with controls (3–8 GPU-hours)

- **ESMFold**, not AlphaFold2. Single-sequence, no MSA — which is the correct choice here precisely
  because these proteins have no deep MSA to give AF2 an advantage, and it removes the MSA-depth
  confound from the pLDDT comparison. The full 3,343 set is 485,962 residues, median 109 aa; on one
  A100 that is a few GPU-hours, so folding everything is affordable even before attrition.
- **Fold three sets in the same batch — this is the part the project has not planned and the part
  that makes the result interpretable:**
  1. the Tier-1 survivors (the test set);
  2. **1,000 length-matched Pfam-*positive* dark family representatives** (positive control: known
     real plasmid proteins, so their pLDDT distribution calibrates "foldable");
  3. **1,000 length-matched uncalled ORFs** from the same plasmids (negative control — I generated
     61,753 of these in §Q8; the extraction is already written).
- Report the pLDDT distribution per set, not a mean.

**Proves:** whether the candidates are foldable at all, *relative to real and to spurious proteins
of the same length*. **Fails to prove:** function. And pLDDT on a 70 aa protein is weak evidence
either way — which is exactly why control set (3) is non-negotiable. Without it, a mean pLDDT of 65
is uninterpretable; with it, "candidates 68, real controls 71, spurious controls 42" is a result and
"candidates 45, spurious controls 44" kills the catalogue.

### Tier 3 — structural search (4–8 h, CPU + ~250 GB disk)

- `foldseek easy-search` of the predicted structures against **AFDB50** (or AFDB/Proteome), **PDB**,
  and **ESM Atlas** (the last is the one that matters most for metagenomic proteins).
- Report per query: best E-value, best TM-score (`--alignment-type 1`), and the taxon of the hit.
  Call structural homology at **E < 1e-3 and TM ≥ 0.5**; report the 0.4–0.5 band separately.
- Restrict interpretation to queries with mean pLDDT ≥ 70. Below that, a Foldseek non-hit is not
  evidence of a novel fold.

**Proves:** for high-pLDDT queries, whether the fold is known. A hit with TM ≥ 0.5 to a
characterized protein is a function hypothesis and removes the family from the novel catalogue. A
high-pLDDT, no-hit family is a genuine novel-fold candidate — the strongest claim available without
experiments, and the one that would make this publishable in the FESNov/NMPFamsDB conversation.
**Fails to prove:** that a no-hit family is functional. And a low-pLDDT no-hit family proves nothing
at all — that is the honest ceiling of the whole approach, and it should be stated in the abstract
rather than discovered by a referee.

### The one-line summary of the protocol

**Tier 1 is the decisive test, not Tier 2.** Sequence exhaustion is ~1 day of CPU and will change the
headline number more than folding will. Folding is what turns the survivors into a *contribution*;
sequence exhaustion is what makes the contribution honest. Doing Tier 2 before Tier 1 would be
publishing a structure survey of proteins that are in UniProt.

---

## Scores

| axis | score | justification |
|---|---:|---|
| **Novelty-claim validity** | **35** / 100 | The single external test is well designed — `--cut_ga` avoids threshold-shopping, and the `NAMED` positive control is genuinely good practice that most papers in this field skip. But novelty rests on one database at a precision-tuned threshold, with no UniRef, no orthology, no profile–profile and no structure, and the cost is now measured rather than speculated: **28.4% of the deliverable and 54.8% of the 93-family core fall to a relaxed threshold on the same database plus the project's own cited clustering criterion** (F1, F4) — before the field-standard cascade even starts. The recurrence argument that was supposed to carry the functional interpretation inverts under a sampling null (F3). The claim is not wrong, but as published it is an upper bound presented as a result. Points retained for the self-critical framing ("measurably both", the 28.5% artifact half) — the project found and published its own previous biggest deflation, which is why I expect it to absorb this one. |
| **Catalogue usability** | **25** / 100 | `dark_novel_families.tsv` has no sequences, no profiles, no MSAs, no stable accessions, no per-family confidence, no host taxonomy, and drops the length column that matters most; its representatives are keyed on ids that F2 shows are not stable across reruns, and the only FASTA of them uses a disjoint identifier space in a different directory. Another lab cannot ask "is my protein in your family?", which is the function of a family catalogue. Points for the fact that the underlying data — lineage, habitat, host, Pfam accession, length — all exists upstream and only needs joining, so the gap is a day's work, not a redesign. |
| **Methods currency (vs field standard)** | **40** / 100 | MMseqs2 + HMMER + Pfam is a competent 2015-era stack. The field standard as of the *Nature* 2024 / *NAR* 2024 papers is a homology cascade terminating in structure, plus profile-based family definition. Absent: Foldseek, HHsuite, any structure predictor, any protein language model, eggNOG, UniRef, dN/dS. Points earned back — and this is not a small amount — for the threshold sweep, the `--cluster-reassign` compliance audit (measuring your own clustering's non-compliance and publishing 89.12% is rare and admirable), the dbAPIS cross-hit adjudication that caught `AcrIIA21` firing on `RepA_N`, and the DefenseFinder run that used genomic context as an independent evidence rule. Those are the instincts of a careful lab using out-of-date tools. |
| **Cross-disciplinary awareness** | **55** / 100 | High awareness, low execution. The literature table in `clustering_threshold_sensitivity.md` §2 correctly identifies FESNov, NMPFamsDB, UniRef and PlasX as the comparison set; Foldseek is named as the next step in five documents; `zoom_two_dark_orfs.ipynb` independently discovered the coverage-non-nesting problem and the ColE1 gene-neighbourhood signal, which is exactly the right cross-field instinct. But the borrowing stops at citation: the FESNov comparison omits FESNov's actual novelty pipeline (F8), the neighbourhood insight is never generalised beyond two ORFs (§X2), and no protein language model appears anywhere. |

---

## Genuine strengths

These are real and I would say so in a referee report.

1. **The `NAMED` positive control.** Searching the 78,397 PlasAnn-named proteins in the same pass is
   the single best design decision in this project. Without it, 26.3% would be uninterpretable. Most
   published novel-family screens have no such control, and the 76.7% figure is what makes the dark
   set's depletion a measurement rather than an assertion.
2. **Publishing the defect log instead of quietly fixing it.** `pfam_dark_validation.md` §6.2
   documents a CDS-indexing bug that silently mislabelled 10.7% of dark ORFs, explains that it was
   caught only because two independently-computed totals disagreed (39.4% vs 35.5%), gives the
   corrupted cross-tab alongside the corrected one, and adds an assert. That is better scientific
   hygiene than most published pipelines.
3. **Measuring your own clustering's non-compliance.** Reading MMseqs2's help, noticing
   `--cluster-reassign` was missing, then aligning 1,066 members back to their own representatives
   and publishing **89.12%** — with the honest note that homology is not in doubt (worst E-value
   1.2e-43) but the stated criterion was not the enforced one — is a level of self-audit I rarely
   see. It is also why F4's fix is available: the work is done, it just has not been promoted.
4. **The anti-defense negative, and how it was reached.** Rejecting an attractive hypothesis at
   0.042%, catching that the raw dbAPIS signal was 3× inflated, and specifically identifying that
   the largest "anti-CRISPR" signal was `AcrIIA21` firing on replication initiators, then
   corroborating with an independent tool that produced *the same artifacts* — that is a properly
   conducted negative result. The `psiA`/`psiB` observation (classic anti-defense genes are in the
   *named* proteome, not the dark one) is a genuinely informative side finding.
5. **The ColE1 cassette analysis.** `zoom_two_dark_orfs.ipynb` is the best piece of biology in the
   set: two "independent" families shown to be one cassette (Jaccard 0.80), positioned at a fixed
   ±2-gene offset around a relaxase on 87% of carriers, on a Col156 replicon, with the
   per-carrier-vs-per-pairing denominator trap caught and reported. It also produced the coverage
   non-nesting insight that F4 is built on.
6. **Threshold honesty.** Reporting that the core target set collapses 93 → 61 → 39 → 19 → **5**,
   and that the most prevalent dark family is a 30%-threshold object whose "88 lineages" is a
   property of the threshold rather than of a protein, is the opposite of how this result would
   usually be presented.

---

## X — The three highest-value cross-disciplinary moves, ranked

### X1 — Finish the homology cascade before anything else (1 day CPU; changes the headline)

Not glamorous, and it is the highest-value move by a wide margin. UniRef50 + eggNOG + profile–profile
will, on my estimate, reassign 40–60% of the 3,343. Every downstream investment — folding, structural
search, experiments, the manuscript's framing — is mis-sized until that number is known. This is
also the move that converts the project from "we found dark matter" to "we bounded dark matter",
which is the claim the field will accept. Tools: MMseqs2 (already installed), DIAMOND (already
installed), eggNOG-mapper, HHsuite.

### X2 — Generalise the gene-neighbourhood signal from two ORFs to the whole catalogue (2 days; free data)

The ColE1 result *is* a genomic-context function inference, and the project treats it as a one-off
anecdote. It is a method, and everything it needs is already on disk. I ran it (validating my CDS
index reconstruction at 100% precision / 0% false positives against `dark30_cluster.tsv`) on the 12
most prevalent Pfam-negative "novel" families:

```
family                                  members   modal flanking annotated category   share   (% flanks annotated)
IMGPR_plasmid_3300042343_000004|14         747     Plasmid Maintenance/Replication      73%    21%
GenBank_CP027661.1|11                      637     Plasmid Maintenance/Replication      88%    52%
IMGPR_plasmid_2914076620_000004|7          539     Plasmid Maintenance/Replication      99%    39%
IMGPR_plasmid_3300037293_000337|6          505     Toxin-Antitoxin System               93%    36%
IMGPR_plasmid_2630968543_000001|3          493     Conjugation                          59%    79%
IMGPR_plasmid_3300022589_000018|2         2487     Toxin-Antitoxin System               77%     1%
```

`GenBank_CP027661.1|11` — a deliverable family on 621 plasmids and 76 lineages — sits beside a
replication/maintenance gene on 88% of flanks with half its flanks annotated. That is a testable
function hypothesis for a "novel" protein, obtained from data already computed, and it is exactly the
evidence class (guilt by association / conserved gene neighbourhood) that the antidefense notebook
names as "the one test that could support the hypothesis where sequence search structurally cannot"
and then does not run. Build it properly: per family, the flanking-gene category distribution and its
entropy, an intergenic-distance histogram (operon evidence), and a family–family co-occurrence network
over all 71,414 replicons. Ship those columns in the catalogue. Borrow from the operon/synteny and
gene-neighbourhood literature (the STRING neighbourhood score, `mmseqs`-based context conservation,
GeneFuncPred-style approaches) — none of it needs a GPU.

### X3 — ESMFold + Foldseek with matched positive *and* negative controls (1 week; §S Tiers 2–3)

Third, not first, and only on Tier-1 survivors. The design point that makes it worth doing is the
control structure in §S Tier 2: fold length-matched Pfam-positive dark families and length-matched
uncalled ORFs alongside the candidates. Without those, pLDDT on 100 aa proteins is not evidence and
a referee will say so. With them, this becomes the decisive experiment and the paper's centrepiece.

**Runner-up, worth naming:** protein language models. ESM2/ProtT5 embeddings plus k-NN label transfer
from Swiss-Prot, or DeepFRI / ProteInfer for direct GO prediction, would give per-family function
hypotheses and — more usefully here — an embedding space in which the 3,343 families can be
positioned relative to the known plasmid proteome. Embedding-space clustering is also the modern
answer to F4's twilight-zone problem: it groups remote homologs that 30% identity cannot. It is
cheap (hours on one GPU for 3,343 sequences) and it is the most conspicuous absence from a 2026
dark-matter project.

---

## Two things I could not check

- Everything upstream of the master table (PlasAnn gene calls' internals, MOB-suite, RGI, the
  BioSample/GOLD environment scraping) is tool output taken on trust, as the verification pack notes.
  My §Q7–Q9 tests the gene calls' *output* against nulls, which is not the same as validating the
  caller.
- No figure was inspected.
