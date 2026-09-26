# Full-run estimate: wall time, CPU-hours, memory and storage

Scope: one submission of `workflow/run_pipeline.sbatch` (one Slurm job, 96 cores, 600 GB,
`MEM_MB=550000`, 10-day limit) over the non-simulated analysis set: 143,503 plasmids,
8,700,021,996 bp, 9,317,050 ORFs, 3,497,616 unique proteins (`hmmer_z` documents the last
two). Code state: contract revision 6 of the labels-build work (narrow_at 0.7, AntiFam skip,
`dark_cooccurrence`, label databases, CONJScan, structure search on CPU inside the job).

Produced by leaf 1.8 on 2026-09-26. The per-stage numbers were computed by a scratch script
(`estimate.py`, reproduced in the Basis section by its rules and constants) from
`results_test/benchmarks/*.tsv` of test job 6990867 plus the measurements listed in Basis.
Labels used below: **MEASURED** (at full scale or on a large sample), **SCALED** (a test-run
measurement multiplied by a driver), **UNBENCHMARKED** (no measurement of this configuration
exists), **UNCERTAIN** (a stated assumption dominates the figure).

## Totals

| Quantity | Estimate | Range / note |
|---|---|---|
| Wall time on one 96-core node | **~72 h (3.0 days)** | ~90 h if the cascade selects as many proteins as the old broad-family proxy; ~105 h if in addition ClusteredNR costs per residue what full nr did (see Uncertainties) |
| CPU-hours actually used | **~3,050 CPU-h** | T5 1,455 + ISEScan 987 + structure search 467 = 95% of it |
| Core-hours billed (96 x wall) | **~6,900 core-h** | 96 x 72.1 h; ~8,700 at 90 h, ~10,100 at 105 h |
| 10-day limit (240 h) | **safe** | the worst case above (~105 h) uses 44% of the limit |
| Memory against `MEM_MB` 550 GB | **sufficient for every measured stage**; largest measured peak is T5 at 185 GB | UNCERTAIN for the Python table stages (see Per stage, memory) |
| Final output | **~81 GB** | plus ~1.5 GB Mash distances |
| Peak disk during the run | **~115 GB** in `results/` | output plus the largest concurrent temporaries |
| Fits on the project filesystem | **yes** | 26 TiB free of the 30,000 GiB quota |

How the wall time is composed. A rule that declares `threads: workflow.cores` takes the whole
node, so those rules (orf_call, tier_search T1-T5, family_network, defence_systems,
conjugation_systems, integrons, is_elements) run one at a time and nothing runs beside them:
their walls add up to **37.5 h** (T5 19.3 h, is_elements 17.4 h, the rest 0.8 h). The other
rules run between them, several at once, and their critical chain is
pre-cascade **2.1 h** (dereplicate, PlasmidScope import, AntiFam screen, the three family
clusterings, the search clustering) plus post-cascade **32.5 h** (cascade_resolve, dark_set,
protein_families, **structure_search 30.7 h** on 16 threads, context, co-occurrence, report).
plasmid_lineage (6.9 h, 16 threads), label_databases (1.2 h), orthology (1.6 h) and the
evolution rules fit inside those windows. Total 37.5 + 2.1 + 32.5 = **72.1 h**.

With the GPU split (`STRUCTURE_ON_GPU=1`, then `structure_gpu.sbatch`) the CPU job loses the
30.7 h structure search: ~41 h of CPU-node wall and ~4,000 core-h billed, plus one GPU job.
The GPU time itself is **UNBENCHMARKED** (the only GPU test, job 6876033, was cancelled before
it ran); `structure_gpu.sbatch` requests one L40S for up to 12 h.

## Per stage

Columns: declared threads; full-run wall (h) at those threads; CPU-h; peak memory (GB); the
scaling driver; the basis. `node` = rule takes the whole node (`workflow.cores`). Test values
are from job 6990867 (100 plasmids, 16 cores). Drivers (full / test): bp x1,607;
ORFs x1,644 (9,317,050 / 5,666); unique proteins x660 (3,497,616 / 5,297); selected proteins
x895 (254,155 / 284); tier-1 queries x323 ((152,324 + 1,000) / 475). Where no larger
measurement exists a multi-threaded stage is assumed to reach 50% parallel efficiency
(UNCERTAIN); single-threaded stages scale their test wall linearly (an upper bound, because
the test wall is dominated by start-up).

| Stage | Threads | Wall h | CPU-h | Peak GB | Driver | Basis |
|---|---|---|---|---|---|---|
| analysis_set | 1 | 0.8 | 0.5 | 6.2 | bp | SCALED: test 1.8 s / 1.1 CPU-s |
| orf_call | node | 0.1 | 4.8 | 32 (declared) | bp | SCALED: test 10.1 s / 10.8 CPU-s |
| orf_index, dereplicate, sweep_cohort, quality_gate, dark_set, extract_cds, cascade_resolve | 1 | 0.03-0.23 each | < 0.1 | < 0.1 | ORFs or unique proteins | SCALED: each <= 1.3 s in the test |
| plasmidscope_import | 1 | 0.08 | 0.1 | ~1 | fixed: reads the whole 5.2 GB PlasmidScope table | MEASURED: test 193.5 s; proxy scan of the same table 233 s |
| artefact_screen | 16 | 0.1 | 0.8 | 17.7 | unique proteins | SCALED: test 4.9 s / 4.5 CPU-s |
| protein_clustering | 32 | 0.84 | 13.4 (upper) | not measured (declared 128) | unique proteins | MEASURED at 3,312,766 proteins, 16 threads: intermediate 492 s, broad 1,880 s; close assumed <= intermediate; x1.056 |
| cascade_selection | 32 | 0.03 | 0.5 | ~5 | selected proteins | MEASURED: search clustering of 254,155 proteins 7 s (16 threads) |
| plasmid_lineage | 16 | 6.9 | 49 | 0.6 | plasmids squared (all-vs-all Mash) | MEASURED on 19,999 plasmids: dist 481 s / 3,423 CPU-s, sketch 6.7 s; x51.5 (n^2) |
| tier_search T1 (Pfam GA) | node | 0.07 | 3.4 | 0.1 | 153,324 queries | SCALED: 0.079 CPU-s/query |
| tier_search T2 (Pfam 1e-5) | node | 0.06 | 3.0 | 0.1 | 127,178 queries (394/475 of T1) | SCALED: 0.084 CPU-s/query |
| tier_search T3 (pharokka) | node | 0.25 | 11.9 | 1.9 | 126,533 queries (392/475) | SCALED: 0.337 CPU-s/query |
| tier_search T4 (Swiss-Prot) | node | 0.04 | 2.0 | 1.8 | 126,533 queries | SCALED: 0.058 CPU-s/query |
| **tier_search T5 (ClusteredNR)** | node | **19.3** | **1,455** | 185 | 119,028 queries, 18.6 M residues | **UNBENCHMARKED**: full-nr fit x0.636 (residues); CPU at the measured 75.5 busy cores |
| orthology (eggNOG-mapper) | 16 | 1.6 | 12.7 | 7.4 | 59,959 cascade-named queries | MEASURED: job 6990877, 2,000 queries 4,658 s, 20,000 queries 4,965 s |
| label_databases | 32 | 1.2 | 16.1 | 3.7 | unique proteins | MEASURED: job 6994476, 50,000 proteins 124.5 s / 938.6 CPU-s |
| protein_labels | 1 | 0.7 | < 0.1 | < 0.1 | unique proteins | SCALED: test 3.7 s |
| protein_families | 1 | 0.2 | 0.1 | 66 (linear upper bound) | unique proteins | SCALED: test 1.1 s, 152 MB |
| family_network | node | 0.06 | 0.3 | 1.6 | network nodes | MEASURED: all-vs-all MMseqs2 over 116,546 node representatives 200 s / 973 CPU-s (16 threads) |
| rarity | 1 | 0.9 | 0.8 | 10 | selected proteins | SCALED: test 3.8 s |
| synteny | 1 | 0.5 | 0.2 | 0.3 | ORFs | SCALED: test 1.0 s |
| recurrence | 1 | 0.9 | 0.7 | **251 (linear upper bound; declared 16)** | ORFs | SCALED: test 2.1 s, 206 MB |
| family_evolution | 16 | 0.05 | 0.4 | 0.5-431 (see memory note) | selected proteins | SCALED: test 3.1 s, 543 MB over 16 workers |
| consensus_recheck | 8 | 1.7 | 6.9 | < 0.1 | selected proteins | SCALED: test 16.7 s / 27.9 CPU-s |
| defence_search | 16 | 0.14 | 1.1 | 79 (linear upper bound) | unique proteins | SCALED: test 53.6 s / 6.2 CPU-s |
| defence_gembase | 1 | 0.3 | < 0.1 | < 0.1 | ORFs | SCALED: test 0.7 s |
| defence_systems | node | 0.04 | 1.9 | **220 (linear upper bound; declared 48)** | ORFs | SCALED: test 145.5 s / 4.2 CPU-s, 187 MB |
| conjugation_systems | node | 0.04 | 2.0 | ~7 | ORFs | SCALED: leaf 1.3, 0.77 ms CPU per ORF; test 92.7 s / 1.8 CPU-s |
| integrons | node | 0.15 | 7.1 | 1.5 | bp | SCALED: test 41.6 s / 15.9 CPU-s |
| **is_elements (ISEScan)** | node | **17.4** | **987** | 5.2 | bp | MEASURED: job 6994476, 1,000 random plasmids (62.1 Mbp) 1,342 s / 25,364 CPU-s at 32 threads (utilisation 0.59) |
| **structure_search (CPU)** | 16 | **30.7** | **467** | ~12 | 177,193 dark representatives | SCALED: two-point fit, 0.623 s and 9.48 CPU-s per query |
| context_features | 1 | 0.9 | 0.7 | < 0.1 | ORFs | SCALED: test 2.0 s |
| dark_cooccurrence | 1 | 0.02 | < 0.1 | ~1 (ceiling 16) | dark-family pairs per plasmid | SCALED: test 0.2 s; leaf 1.9: ~0.73 M pair occurrences, < 1 GB |
| annotation_report | 1 | 0.3 | 0.2 | **167 (linear upper bound; declared 32)** | ORFs | SCALED: test 0.8 s, 154 MB |
| feature_files, annotate_plasmids | 1 | 0.8, 1.3 | 0.5, 0.1 | 21, < 0.1 | ORFs | SCALED |
| preflight, prepare_control, prepare_decoys, clonal_registry, check_hmmer_z | 1 | < 0.02 | < 0.1 | < 0.1 | fixed | test values; independent of the set size |

**Memory against 550 GB.** Every measured peak is below its declared `mem_mb` and far below
550 GB: T5 185 GB (full nr, `-b 16 -c 1`, 96 threads; ClusteredNR with the same block size is
expected to be similar, declared 220 GB), protein clustering is declared 128 GB, and every
other measured stage is below 8 GB. **UNCERTAIN:** for the Python table stages the only
memory measurement is the 100-plasmid test, and the column above scales its excess over a
50 MB interpreter baseline linearly with the driver. That linear bound exceeds the declared
`mem_mb` for recurrence (251 GB vs 16), defence_systems (220 vs 48), annotation_report (167
vs 32) and protein_families (66 vs 64). family_evolution's 543 MB is 16 worker processes,
so it does not scale with the data. Snakemake schedules by declared `mem_mb`, so if these
bounds were real, two of these rules running together could exceed the 600 GB allocation.
The bound is pessimistic (much of a small test's RSS is fixed), but it is not refuted by any
measurement. A 1,000-plasmid run would settle it; until then the risk is recorded here.

## Basis

Everything below was measured; the source files named here (`results_bench/`,
`results_test_small/`, `results_test_ps/`, `results/logs/bench_nr.*`, `results/bench/nr/`,
`data/bench/`) are to be deleted, so every number used is copied here.

**Test run, job 6990867** (this leaf; `config/test/config.yaml`, 100 plasmids, 16 cores,
64 GB, all 45 rules incl. label_databases, conjugation_systems, structure_search on CPU,
dark_cooccurrence and the report): wall 1:41:17, 15:37:13 CPU (56,233 CPU-s), batch MaxRSS
7.9 GB, COMPLETED; `check_testrun.py` prints TESTRUN_OK with labels from amrfinder, bacmet,
card, dbapis, mobileog, oritdb and tadb. Counts: 5,666 ORFs; 5,297 unique proteins; 284
selected (275 representatives + 9 members), 2,729 Tier 0, 2,280 not selected, 4 AntiFam;
tier queries T1 475 (275 + 100 controls + 100 decoys), T2 394, T3 392, T4 392; T4 left 346
unresolved = 242 real + 4 controls + 100 decoys, of which 28 are named by Pfam or Swiss-Prot
and skip T5, so **214 of 275 representatives (77.8%) would reach T5** (binomial 95% interval
about 73-83%); 209 dark proteins, 198 structure queries; orthology queried 67; dark families
199. Per-rule benchmarks (s, CPU-s, max RSS MB) are those cited in the table.

**Cascade selection re-measured at the intermediate level** (proxy on PlasmidScope's proteins
of the analysis set, 2026-09-25): 8,753,749 ORF rows, 3,312,766 unique proteins, 1,239,766 not
annotated by Tier 0 (`plasmidscope.ps_class`), 148,284 of them on a small plasmid, 270,644
unique proteins on a small plasmid. All proteins clustered as the pipeline does (MMseqs2
easy-cluster, 50% identity, 80% coverage, cov-mode 1, cluster-mode 2, cluster-reassign):
1,007,019 families, 78,088 selected, **254,155 selected proteins**, search clustering (90/80,
cov-mode 0) **152,324 representatives** with 23,777,414 residues (mean 156 aa). The same
procedure at the broad level (30/50) gives 580,964 families, 334,383 selected, 203,193
representatives, 34,780,713 residues. The earlier figure (356,959 selected / 216,546
searched; broad families clustered from the unannotated proteins only) is therefore an upper
bound. Selection happens before the cascade, so **narrow_at 0.7 does not change these
counts**; it changes only how many representatives reach the deeper tiers, which the test
run measures at 0.7 (above). The AntiFam skip changes nothing measurable: the test's 4
AntiFam proteins were all outside the selection already (leaf 1.9); the proxy was not
AntiFam-screened (UNCERTAIN, expected negligible).

**nr benchmark, job 6928733** (96 threads, 600 GB, `--very-sensitive -b 16 -c 1`,
`--max-target-seqs 5`, full nr 2025-03: 707 M sequences, 273 G residues). test5 cold: 152
queries, 32,327 residues, 4,158.83 s, max RSS 135.8 GB, 151 with a hit; test5 warm 2,700.21 s;
test100: 3,808 queries (the smoke set's T4-unresolved proteins, 877,792 residues, mean 230.5
aa), 8,942.70 s, max RSS 185.3 GB, 3,689 with a hit. Job TotalCPU 13-19:16:46 (1,192,606
CPU-s) over 15,802 s of runs = 75.5 busy cores. Fit: 5.658 ms per query residue after a
fixed 3,976 s. An earlier attempt with DIAMOND's default `-b 2 -c 4` (job 6823981, 2,000
queries, 16 threads) did not finish in 12 h.
T5 now: 0.778 x 152,324 + 500 decoys + ~20 controls = 119,028 queries x 156 aa = 18.6 M
residues -> full nr 109,100 s (30.3 h); **ClusteredNR (173,589,472,897 residues, 545,597,271
sequences, `data/refs/clustered_nr/metadata.json`) is 0.636 of full nr's residues -> 19.3 h.
UNBENCHMARKED** (user decision not to benchmark ClusteredNR): the scaling assumes DIAMOND's
cost is proportional to database residues; the benchmark used 5 target sequences and
production uses 25 (unmeasured cost); linear extrapolation from 0.88 M to 18.6 M residues is
conservative for DIAMOND, whose per-letter cost falls as the query block grows.

**Earlier tier benchmarks on the same 100 plasmids** (CPU-s per query, used for T1-T4):
results_test_ps (16 cores): T1 2,780 queries 118.71 s / 220.95 CPU-s; T2 2,696, 116.02 s /
226.12; T3 2,694, 235.26 s / 908.28; T4 2,694, 15.09 s / 155.96. results_bench (96 cores, all
proteins searched): T1 5,504, 394.99 s / 6,834.88; T2 5,061, 391.74 s / 7,053.14; T3 5,055,
145.34 s / 577.55; T4 5,055, 6.26 s / 299.16; T5 (full nr) 3,910, 10,485.64 s / 755,905 CPU-s,
max RSS 180 GB. results_bench's T1/T2 CPU per query is 15x higher (hmmsearch at 96 threads);
using it would add ~100 CPU-h and < 1 h of wall.

**Structure search** (ProstT5 + Foldseek on CPU, 16 threads, scope representatives): test
job 6990867 198 queries 157.7 s / 1,770.9 CPU-s, 10.9 GB; results_test_ps 1,855 of 1,947
dark proteins 1,190.55 s / 17,486.41 CPU-s, 12.0 GB; results_test_small 166.56 s / 1,686.13
CPU-s. Fit 0.6233 s + 34.3 s fixed, 9.48 CPU-s per query. Full: 198 x 895 = 177,193 queries
(selected-protein driver; UNCERTAIN, the dark fraction of the full selection is unknown
before the cascade).

**eggNOG-mapper, job 6990877** (16 cores, `-m diamond`, the rule's command, random samples of
the proxy's searched representatives): 2,000 queries 4,657.89 s / 45,619 CPU-s, 7.24 GB, 850
annotated; 20,000 queries 4,965.08 s / 41,913 CPU-s, 7.40 GB, 8,327 annotated. Slope 0.0171
s/query after ~4,624 s fixed. Earlier 100-plasmid runs: 72 queries 3,543 s / 44,513 CPU-s;
587 queries 3,564 s / 45,135; 3,493 queries 2,383 s / 32,563. Full queries: 67 x 895 = 59,959.

**Scale benchmarks, job 6994476** (32 cores, the pipeline's own rules via Snakemake):
label_databases on 50,000 random proxy proteins 124.49 s / 938.59 CPU-s, 3.69 GB, 8,336 label
rows on 7,433 proteins; is_elements on 1,000 random analysis-set plasmids (62,110,556 bp)
1,342.26 s / 25,364.18 CPU-s, 5.22 GB, mean load 18.5 of 32. Linear in bp:
4.08e-4 CPU-s per bp (the test run gives 4.00e-4: 2,163.6 CPU-s for 5,412,798 bp).

**Other full-scale measurements** (login node, 16 threads): family clustering of 3,312,766
proteins, intermediate 492 s, broad 1,880 s, broad temporary+output peak 2.55 GB; Mash on a
seeded 19,999-plasmid sample (1,230,140,713 bp): sketch 6.66 s / 39.3 CPU-s, dist 481.05 s /
3,422.6 CPU-s, 0.57 GB, 456,013 rows (29.8 MB); network search over 116,546 intermediate
families with a small-plasmid protein 200.38 s / 972.9 CPU-s, 1.56 GB, 2,476,757 hits
(218 MB).

**Earlier 100-plasmid runs for reference** (wall s): results_test_small (job 6970802, 1:23:02,
14:05:12 CPU): orthology 3,543.01, is_elements 441.03, structure_search 166.56,
defence_systems 99.19, tier T3 180.33; results_test_ps (jobs 6957648/6963002, ~1:40 each):
orthology 3,564.24, structure_search 1,190.55, is_elements 432.93, T3 235.26. Previous test
job 6985208 (before revisions 5-6): 1:33:10, 14:54:38 CPU.

## Storage

**Filesystem.** The run writes `results/` on `/gorilla/proj` (= `/proj/h-mel-phylo`, shared
network filesystem, 15 PB, 4.7 PB free). Project quota (`uquota`, 2026-09-26):
3,790.6 GiB used of 30,000 GiB, 489,107 files of 2,000,000. DIAMOND spills to
`resources.tmpdir` = node-local `/scratch`.

**(a) Final output per stage**, from `du` of results_test (job 6990867) with each file scaled
by its stage's driver; files that already list the whole analysis set (clonal_registry.tsv,
analysis_set.txt, small_plasmids.txt) and the label-database copies are fixed:

| Stage directory | Test | Driver | Full run |
|---|---|---|---|
| 01_analysis_set | 27.9 MB | bp x1,607 (22 MB fixed) | 8.7 GB |
| 02_orf_calling | 3.7 MB | ORFs x1,644 | 6.1 GB |
| 03_dereplication | 3.1 MB | unique proteins x660 | 2.1 GB |
| 04_orf_qc | 0.2 MB | unique proteins | 0.2 GB |
| 05_annotation_cascade | 1.2 MB | unique proteins | 0.8 GB |
| 06_annotation_tables | 9.5 MB | ORFs | 15.6 GB |
| 07_orthology | 0.4 MB | unique proteins | 0.3 GB |
| 08_protein_labels | 662.1 MB | unique proteins (653 MB database copies fixed) | 6.9 GB |
| 09_quality_gate | 0.3 MB | unique proteins | 0.2 GB |
| 10_clustering | 6.5 MB | unique proteins; + Mash distances measured, 29.8 MB x51.5 = 1.5 GB | 5.8 GB |
| 11_distribution_and_evolution | 0.3 MB | selected proteins x895 | 0.3 GB |
| 12_context_and_structure | 18.4 MB | ORFs (DefenseFinder/MacSyFinder/CONJScan run directories are kept) | 30.3 GB |
| 13_synteny, 14_rarity | 0.2 MB | ORFs, selected proteins | 0.2 GB |
| 15_report | 1.7 MB | ORFs | 2.8 GB |
| logs, benchmarks | 0.5 MB | ORFs | 0.8 GB |
| **Total** | **736 MB** | | **~81 GB** |

File count does not scale with the data: 9,712 files in the test, 9,526 of them one per HMM
model in the kept MacSyFinder `hmmer_results` directories.

**(b) Peak temporaries.** The pipeline creates working directories with
`plasmidann.scratch.scratch_dir` and deletes them on success (`scratch.release`): the MMseqs2
tmp of each family clustering, the search clustering, the network search (`network_tmp`, 218
MB of hits at full scale), the Foldseek tmp (`foldseek_tmp`), family_evolution's tmp and each
cascade tier's directory. They are kept only when a stage fails. NOT deleted: the
label_databases work directory (653 MB of DIAMOND databases plus hits), the DefenseFinder
`phase1`/`phase2` and CONJScan `conjscan/run` directories and `10_clustering/mash` (all counted
in (a)). Measured peaks: a 15-second disk sampler on job 6990867 saw 05_annotation_cascade at
7.8 MB (final 1.2), 10_clustering 8.4 MB (final 6.5) and 12_context_and_structure 50.1 MB
(final 18.4); at full scale the broad clustering of 3.31 M proteins peaked at 2.55 GB
(temporaries + outputs). Scaled: cascade tier temporaries ~2 GB (x323), Foldseek/defence
temporaries ~28 GB (x895, UNCERTAIN), MMseqs2 ~3 GB. DIAMOND's T5 spill on node-local
`/scratch` is **not measured** (the nr benchmark wrote to node-local `$TMPDIR` without
recording its size); with 119 k queries and 25 targets each it holds a few GB of alignments
(UNCERTAIN), and it is not in the project quota.

**(c) Reference data and environments already on disk** (`du`, 2026-09-26): data/refs
720 GB in total - clustered_nr 658 GB (DIAMOND database 208 GB + the downloaded BLAST volumes
450 GB, which DIAMOND no longer needs once `clustered_nr.dmnd` exists), eggnog 48 GB,
foldseek 6.5 GB, pfam 4.5 GB, pharokka 1.8 GB, labels 0.87 GB, macsyfinder 0.34 GB,
amrfinder 0.24 GB, antifam 0.05 GB, conjscan 0.02 GB, control 1.7 MB. envs 7.8 GB
(plasmidann 3.0, pharokka 1.7, amrfinder 1.0, conjscan 0.66).

**(d) Total and peak.** Final output ~81 GB; peak during the run ~115 GB in `results/`
(output plus the largest set of concurrent temporaries, ~33 GB). Against the quota:
3,791 + 115 = ~3,906 GiB of 30,000, so **the full run fits** with a wide margin, and removing
`data/refs/clustered_nr/blastdb` would return 450 GB.

## Uncertainties and unbenchmarked items

- **T5 on ClusteredNR is UNBENCHMARKED.** 19.3 h assumes cost proportional to database
  residues (0.636 of full nr). If ClusteredNR costs as much per query as full nr: 30.3 h.
  The 5 -> 25 target-sequence change is not measured.
- **Selection size.** 152,324 searched representatives is a proxy on PlasmidScope's proteins,
  not on the pipeline's own gene calls. Using the broad-level proxy instead (203,193
  representatives, 334,383 selected, 171 aa mean) raises T5 to ~28 h and structure search to
  ~40 h: wall ~90 h. The old broad-from-unannotated figure (216,546 representatives) gives T5
  ~30 h (~47 h at full-nr cost).
- **Share reaching T5** (77.8%) and the dark share (198 structure queries per 284 selected)
  come from 275 representatives of the test set (binomial 95% interval ~73-83%).
- **50% parallel efficiency** is assumed for multi-threaded stages without a larger
  measurement (T1-T4, orf_call, integrons, defence, conjugation): together < 1 h of wall.
- **Memory of the Python table stages** at full scale (see Per stage). This is the one item
  that could make the run fail rather than run longer.
- **GPU structure search** has no measurement.
