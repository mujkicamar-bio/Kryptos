# Plan — conservative re-clustering of the dark plasmid proteome

Status: **awaiting approval**. Written 2026-09-01.

## 1. Why

The dark proteome was clustered once, by `scripts/cluster_dark_orfs.sh`:

```
mmseqs easy-cluster --min-seq-id 0.3 -c 0.8 --cov-mode 0
```

Two problems, both now measured rather than assumed:

1. **Single threshold.** Already listed as limitation #2 in `reports/dark_orf_clustering.md`.
   Every headline number in `dark_plasmidome.ipynb`, `dark_families.ipynb` and
   `widespread_dark_orfs.ipynb` is conditional on 30%.
2. **Criterion violations.** MMseqs2's own help states cascaded clustering "can cluster
   sequence that do not fulfill the clustering criteria", corrected by `--cluster-reassign`,
   which we did not pass. Audit of family `IMGPR_plasmid_3300022589_000018|2`: **13.2%** of its
   991 distinct members fail the stated criterion against the representative (123 of 131 fail on
   coverage, 8 on identity). Homology is not in doubt (worst E-value 1.2e-43) but the stated
   criterion is not what was actually enforced.

## 2. Decisions needed before execution

- **D1 — threshold set.** Recommend a sweep, not one value (precedent: PlasX/MobMess ran 12):
  `0.3` (rebuild, comparable), `0.5` (UniRef50), `0.7` (NMPFamsDB), `0.9` (UniRef90),
  all at `-c 0.8 --cov-mode 0 --cluster-reassign`.
- **D2 — input scope.** Small-plasmid dark proteome only (378,552 proteins), or the union with
  the big-plasmid dark ORFs extracted on 2026-09-01 (5,135,708) = 5,514,260? Recommend small
  first, union as a separate follow-on.
- **D3 — coverage.** Keep `-c 0.8 --cov-mode 0` throughout. It is already stricter than the
  Nature FESNov paper (`-c 0.5 --cov-mode 1`) and matches UniRef's 80% overlap rule. Raising it
  further has no precedent in the literature reviewed.

## 3. Non-negotiable constraint

**Nothing overwrites `data/dark_orf_run/dark30_*` or `mix30_*`.** Those files back three
notebooks and two reports. All output goes to `data/dark_orf_run/recluster/dark{ID}_*`.

## 4. Tasks

### T1 — parameterised clustering script
Write `scripts/recluster_dark_orfs.sh <min_seq_id>`, mirroring `cluster_dark_orfs.sh` but with
`--cluster-reassign` and a threshold argument, writing to `data/dark_orf_run/recluster/`.
**Verify:** runs at 0.9 on a 10,000-sequence subset in under 2 minutes and emits a cluster TSV.

### T2 — run the sweep on the small dark proteome
SLURM job, 24 cores, over the 378,552 proteins at each threshold in D1.
**Verify:** four cluster TSVs; for each, `sum(members) == 378,552`, no member appears twice, and
family count increases monotonically with identity.

### T3 — criterion-compliance audit
For each threshold, sample 200 families, `mmseqs easy-search` every member against its own
representative, report the fraction meeting identity AND coverage directly.
**Verify:** reassigned runs reach >=99% compliance, against the **86.8%** measured on the
current 30% clustering. This is the test that proves the rebuild fixed the defect.

### T4 — rebuild family stats and Pfam roll-up per threshold
Reuse `scripts/summarize_dark_orf_clusters.py` and `scripts/summarize_pfam.py`. The per-sequence
Pfam results in `data/pfam_run/pfam_per_seq.tsv` are threshold-independent, so **no new HMM
search is needed** — only the family->Pfam roll-up changes.
**Verify:** per-threshold `darkfam_stats.tsv` and `darkfam_pfam.tsv`; dark ORF totals match T2.

### T5 — threshold-sensitivity report
`reports/clustering_threshold_sensitivity.md`: how every headline number moves across
thresholds — % of dark ORFs in recurrent families, singleton share, core-family count, the
identity of the most prevalent family, and whether the "88 lineages" result survives.
**Verify:** report names the exact scripts that produced it (project standing rule).

### T6 — notebook reconciliation
Add one threshold-sensitivity cell to `dark_families.ipynb` rather than rewriting it, so the
30% results stay reproducible and the sensitivity is visible beside them.
**Verify:** notebook executes end to end; existing cells unchanged.

### T7 — optional follow-on
Cluster the 5,514,260-protein union (small + big dark ORFs) at the chosen threshold. Separate
job; ~15x the input, so size from T2's measured runtime rather than guessing.

## 5. Out of scope

Re-running Pfam (per-sequence results are threshold-independent), re-running DefenseFinder or
dbAPIS, and anything touching the FlashFold target sets.
