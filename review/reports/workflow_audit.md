# Workflow and research-practice audit

*Produced by the orchestrating reviewer (not one of the five panel seats), 2026-09-02. Every claim
below was checked by a command run against the repository; the commands are given.*

This audit covers the dimension the five reviewers were not assigned: **how the work is done**, as
distinct from whether the science or the statistics are right. It is deliberately separate because
the project's own standards are unusually high — every report names its producing scripts, defects
are documented rather than quietly corrected, a locked exclusion decision is honoured across
phases — and the failures below are failures *against that standard*, not against a generic one.

---

## W1 [CRITICAL] The project is not under version control

```
$ git -C . status
fatal: not a git repository (or any parent up to mount point /gorilla)
$ ls -d .git   ->  (absent)
```

62 scripts, 15 notebooks, 13 reports, and a 25 GB data tree, with **no repository, no history, no
diff, no blame, no tags**.

Why this is critical here specifically, rather than a generic best-practice complaint:

- The project has already had **two silent pipeline defects that changed published numbers**
  (`reports/pfam_dark_validation.md` §6.1, §6.2). One of them — the CDS-indexing bug — "affected no
  number *in this report*" but "silently corrupted the first Pfam cross-tab", changing
  59.2/26.8/8.7/5.3 to 57.4/28.5/10.9/3.2. With no VCS there is **no way to establish which
  documents, figures or notebook narratives were written against the corrupted numbers** and which
  were written after the fix. The deck asserts every figure was regenerated
  (`reports/dark_plasmidome_deck.html:846`); that assertion is currently unfalsifiable.
- `plans/conservative_reclustering.md` §3 declares a **"Non-negotiable constraint": nothing
  overwrites `data/dark_orf_run/dark30_*` or `mix30_*`**, because those files back three notebooks
  and two reports. That constraint is enforced by nothing but memory and care. A single mistyped
  output prefix silently invalidates two published reports with no way to detect or undo it.
- Reports carry dates like "*numbers revised the same day after two pipeline defects were fixed*"
  — a statement that only means something if the prior state is recoverable. It is not.

**Fix (30 minutes, no risk):** `git init`; a `.gitignore` excluding `data/`, `envs/`, `tools/`
(25 GB + 5.4 GB + 89 MB); commit `scripts/ notebooks/ notebooks_my/ reports/*.md plans/ manuscript/
*.md environment.yml`. Then tag the current state `dark-plasmidome-2026-09-02` so the deck's claims
are anchored to a recoverable tree. Add a checksum manifest for the frozen `dark30_*`/`mix30_*`
files so the non-negotiable constraint becomes enforceable rather than aspirational.

---

## W2 [CRITICAL] The entire Arc-2 narrative lives in notebooks saved with zero outputs

```
$ for nb in notebooks_my/*.ipynb; do python3 -c "...count outputs..."; done
notebooks_my/analysis.ipynb              code_cells_with_output=0
notebooks_my/antidefense.ipynb           code_cells_with_output=0
notebooks_my/cryptic_plasmids.ipynb      code_cells_with_output=0
notebooks_my/dark_families_90pct.ipynb   code_cells_with_output=0
notebooks_my/dark_families.ipynb         code_cells_with_output=0
notebooks_my/dark_plasmidome.ipynb       code_cells_with_output=0
notebooks_my/defense_systems.ipynb       code_cells_with_output=0
notebooks_my/widespread_dark_orfs.ipynb  code_cells_with_output=0
notebooks_my/zoom_two_dark_orfs.ipynb    code_cells_with_output=0
```

Versus the Arc-1 notebooks, which retain theirs:

```
notebooks/onehealth_amr_investigation.ipynb   code_cells_with_output=55
notebooks/small_cryptic_investigation.ipynb   code_cells_with_output=43
notebooks/plasmid_data_investigation.ipynb    code_cells_with_output=42
```

Every headline number in the dark-plasmidome story — "158 dark ORFs, 0.042%", "3,647 defense
systems on 3,362 plasmids", "dORF0001, 83 aa, on 820 plasmids", "336 plasmids carry both",
"93 families", "five survivors" — exists **only as hand-typed prose in a markdown cell**, with the
code beside it but no preserved output to check the prose against.

This is the exact mechanism by which a documented number can drift from the code that produced it,
and it is not hypothetical: the "562-plasmid AMR definitional leak" quoted in two reports
recomputes as **694** from the master column (`review/VERIFICATION_PACK.md` §7.4), and the leak
figure originates in one of these output-free notebooks.

The regression is the notable part. Arc 1 used a disciplined generator pattern —
`scripts/build_*_notebook.py` emits the notebook, the notebook is executed with
`jupyter nbconvert --execute --inplace`, outputs are committed, and the report cites the script.
Arc 2 abandoned that pattern for hand-authored notebooks **exactly when the science became more
contested and more corrected**.

**Fix:** either (a) execute the `notebooks_my/` notebooks and save with outputs — this is the
cheap fix and it makes every prose number checkable; or (b) return to the Arc-1 generator pattern
for anything that will be cited in a report or manuscript. (a) is enough for now. Note the
standing constraint that these notebooks are edited cell-by-cell in VS Code, so the execute step
must be done deliberately rather than as a side effect.

---

## W3 [MAJOR] Arc 2's compute environment is outside the project and is not captured anywhere

Every dark-plasmidome result was produced by interpreters in `$HOME`:

```
$ grep -rho "conda/envs/[a-z_]*" scripts/ reports/*.md notebooks_my/*.ipynb | sort | uniq -c
     10 conda/envs/panaroo
      4 conda/envs/genesis
      3 conda/envs/plasann_env
```

But the project's captured environment specs are only:

```
$ ls envs/live_env_specs/
rgi_env.yml  segmantx_env.yml  stats_env.yml
```

`genesis` and `panaroo` — which ran MMseqs2 18.8cc5c, HMMER 3.4, and all the pandas/numpy analysis
behind Arc 2 — are **not exported, not listed in `PROJECT_OVERVIEW.md` §6, and not described in
`tools/README.md`**. `PROJECT_OVERVIEW.md` mentions only a `genesis_nb` *kernel* (line 200), never
the environment. They live in `$HOME`, which `tools/README.md` itself warns has "a small disk quota
that a single full bioinformatics env blows through" — i.e. they are in the location the project
documented as unsafe, and they are the location the newest science depends on.

This directly violates the project's own provenance standard, which is otherwise well kept: all 13
`reports/*.md` cite their producing scripts (checked — none missing).

**Fix:** `micromamba env export` (or `conda env export --no-builds`) both envs into
`envs/live_env_specs/`, add them to `PROJECT_OVERVIEW.md` §6 with the tool versions already recorded
in `reports/dark_orf_clustering.md` (MMseqs2 18.8cc5c, HMMER 3.4), and state which env produced
which report.

---

## W4 [MAJOR] `PROJECT_OVERVIEW.md` — declared "the front-door document" — has no knowledge of Arc 2

```
$ grep -n -i "dark\|cryptic\|pfam" PROJECT_OVERVIEW.md
(no matches)
```

It is dated 2026-07-27 and stops there. It still lists as "⏳ next" a Phase-3 programme (ARG
co-occurrence, Inc-resolved carriage, geography) that was abandoned, still describes the NAR
manuscript as the live output, and contains no mention of the 47,031/71,414 payload-free
compartment, the 378,552 dark ORFs, the Pfam adjudication, the defense-system finding, the
threshold sweep, or the six weeks of work that produced them. Five of the thirteen files in
`reports/` are invisible to it.

Anyone opening this repository — a collaborator, a reviewer, the author in six months — is told the
project's current state is something it stopped being in July.

**Fix:** one editing pass. Add an Arc-2 section mirroring §4's structure (stage → script → output →
report), update §2's status table, update §8 to list all 13 reports and both notebook directories,
and correct §6 per W3.

---

## W5 [MAJOR] Two incompatible definitions of "the small payload-free plasmidome" circulate freely

Documented in full at `review/VERIFICATION_PACK.md` §7.5. In short: `< 10 kb` **with** a
conjugation-machinery filter → 47,031 (`reports/small_cryptic_methodology.md`), versus `< 20 kb`
**without** it → 71,414 (`reports/dark_orf_clustering.md` and all nine `notebooks_my/`). Verified
that the first is a strict subset of the second (overlap = 47,031 exactly).

The workflow problem, as opposed to the scientific one: the two sets have no distinct names. Both
are called "the small cryptic/payload-free plasmidome", the ID file for the larger one is named
`cryptic_small_ids.txt`, and the deck moves between the two objects' numbers without signalling the
change of denominator. This is a naming failure that will produce a wrong sentence in a manuscript
with near-certainty.

**Fix:** name them (e.g. `SC10` and `PF20`), state the definition once in a shared definitions
block, and grep every report, notebook and slide for sentences that cross between them.

---

## W6 [MINOR] No tests anywhere, in a codebase whose defects were both silent

```
$ find . -maxdepth 3 -iname "*test*" -not -path "./data/*" -not -path "./envs/*" -not -path "./tools/*"
(nothing)
```

Both documented defects were caught by *ad hoc* consistency checks — "two independently computed
totals disagreed" — not by anything systematic. The project already knows what its invariants are
and states them in prose: cluster members must sum to the input count; dark IDs must be a strict
subset of all-CDS IDs; the two PlasAnn shard series must both be globbed and then deduplicated;
family stats must sum to 378,552. Each of those is a two-line assertion.

**Fix:** a single `scripts/check_invariants.py` asserting the five or six known invariants, run
after any pipeline stage. This is a ~50-line file and it would have caught both historical defects.

---

## W7 [MINOR] Portability limits, worth recording rather than fixing

- 12 sbatch scripts hardcode `-A uppmax2025-2-42` and `-p pelle`. Fine now; a reproduction attempt
  after the allocation expires will fail confusingly. One `#SBATCH -A ${SLURM_ACCOUNT:-...}` pattern
  or a note in each report's Reproduce block resolves it.
- 22 of 62 scripts hardcode `/gorilla/proj/...` or `/gorilla/home/amujkic/...`. Acceptable for a
  single-site project; it should simply be stated as a known limitation rather than discovered by
  whoever tries to run it elsewhere.
- `prodigal` was installed in the (now decommissioned) `plasmid_phase0` env. Independent gene
  calling — named as an unaddressed limitation in three separate reports — was available in-project
  and was never run. Worth noting that the blocker was never tooling.

---

## What is genuinely good here, and should not be lost

These are real, and unusual:

1. **Every one of the 13 `reports/*.md` names the scripts that produced it.** Checked
   programmatically; none missing. This is the single practice that made the verification pack
   possible at all, and it is why 46 of 52 headline numbers could be independently recomputed.
2. **Defects are published, not patched.** `reports/pfam_dark_validation.md` §6 documents both bugs
   including the one that "affected no number in this report", with the corrupted and corrected
   cross-tabs shown side by side. Most projects delete that section.
3. **A refuted thesis was written up as a first-class artifact.** `manuscript/reassessment_round1.md`
   works through six reviewer objections, upholds four, qualifies one, and **refutes one in the
   paper's own favour** — then recommends reframing rather than patching. That is the correct
   response and it is documented in a form that survives.
4. **A locked decision, honoured across phases.** The Simulated-artifact exclusion was made before
   the analysis and is applied identically in every downstream stage (verified: 143,503 everywhere
   it matters).
5. **The sensitivity analysis was run even though it destroyed the headline.** The threshold sweep
   took the core family set from 93 to 5 and was published anyway, with "what does not survive" as
   a named section.

W1–W3 matter *because* of this. A project that documents its own defects this carefully but cannot
reconstruct which version of the data a given document was written against has built a strong
provenance culture on top of an unversioned foundation.

---

## Ranked actions

| # | Action | Cost | What it buys |
|---|---|---|---|
| 1 | `git init` + `.gitignore` + tag current state (**W1**) | 30 min | Makes every provenance claim in every report verifiable rather than asserted |
| 2 | Execute + save `notebooks_my/*.ipynb` with outputs (**W2**) | 1–2 h | Makes the whole Arc-2 narrative checkable against its own code |
| 3 | Export `genesis` + `panaroo` specs into `envs/live_env_specs/` (**W3**) | 15 min | Arc 2 becomes reproducible at all |
| 4 | Rewrite `PROJECT_OVERVIEW.md` to cover Arc 2 (**W4**) | 1 h | The front door stops lying about what the project is |
| 5 | Name the two payload-free sets and audit crossings (**W5**) | 2 h | Prevents a wrong sentence in the manuscript |
| 6 | `scripts/check_invariants.py` (**W6**) | 1 h | Would have caught both historical defects |
