# Tool provenance and install notes (Phase 0.1)

All environments are built with **micromamba** (not the system `conda`, which hangs/times out
on this node — see "Why micromamba" below), rooted at:

```
envs/.micromamba_root/
```

(project-local, not `$HOME`, because `$HOME` has a small disk quota that a single full
bioinformatics env blows through — see "Why project-local root" below.)

To use any environment in a shell:

```bash
export MAMBA_ROOT_PREFIX=/gorilla/proj/h-mel-phylo/h-mel-phylo/private/amar/plasmid_analyis/envs/.micromamba_root
ENV=$MAMBA_ROOT_PREFIX/envs/<env_name>
PYTHONNOUSERSITE=1 PATH="$ENV/bin:$PATH" $ENV/bin/<tool> ...
```

**Both `PYTHONNOUSERSITE=1` and prepending `$ENV/bin` to `PATH` are required every time**, not
optional convenience — see "Two recurring gotchas" below.

## Environments built

| env name | purpose | key packages | exported spec |
|---|---|---|---|
| `plasmid_phase0` | general analysis (seqkit, csvtk, prodigal, hmmer, blast, mmseqs2, Python stack) | python=3.10, seqkit, csvtk, prodigal, hmmer, blast, mmseqs2, pandas, scipy, statsmodels, networkx, seaborn, biopython | `envs/plasmid_phase0.environment.yml` |
| `mob_suite_env` | **MOB-suite** (`mob_typer`) — isolated because it needs an old pandas API | mob_suite 3.1.9, blast, mash 2.3, gsl=2.6 (pinned) | `envs/mob_suite_env.environment.yml` |
| `plasann_env` | **PlasAnn** annotation pipeline | plasann 1.1.6 (bioconda), prodigal, blast, infernal | `envs/plasann_env.environment.yml` |
| `segmantx_env` | **SegMantX** (Hanke & Dagan 2025) duplication detection — run from `tools/SegMantX/` (cloned, not packaged) | pandas<2.2, python<3.13 (pinned), scipy, biopython, plotly, streamlit, blast | `envs/segmantx_env.environment.yml` |

mob_suite, PlasAnn, and SegMantX were each put in their **own** environment rather than one
shared env, after a shared-env solve produced an incompatible pandas for `mob_suite` (see below).
Isolating them avoids any repeat of that problem and lets each tool's own dependency solve pick
versions that actually work with that tool's code.

## Why micromamba, not `conda create`

The system `conda` (`/sw/apps/conda/latest/rackham_stage/bin/conda`, UPPMAX/Rackham-style install)
hung indefinitely on every `conda create`/`conda search` call here, even with
`--override-channels -c conda-forge -c bioconda` and a 60s `timeout` wrapped around the actual
binary. CPU time consumed was tiny (~10s) relative to wall time (>60s), so it isn't a slow
solve — something in this conda install's channel/index handling stalls. The system `.condarc`
also injects a large set of `file:///sw/apps/conda/latest/rackham/local_repo/...` local mirror
channels that get merged in even with `--override-channels` unless you also pass `--no-rc`-style
isolation, compounding the problem.

Fix: downloaded the micromamba 2.8.1 static binary directly
(`https://micro.mamba.pm/api/micromamba/linux-64/latest`) to `~/bin/micromamba` and used it for
every env build with `--no-rc --override-channels -c conda-forge -c bioconda -c defaults`. Solves
that hung indefinitely under `conda` completed in 10-15s under micromamba.

## Why project-local micromamba root, not `$HOME`

First env-creation attempt rooted at `$HOME/.micromamba` failed mid-transaction with
`Disk quota exceeded` while linking `csvtk` — `$HOME` (`/gorilla/home/amujkic`) has a per-user
quota that a ~350-package bioinformatics env (852MB) exceeds. `/proj/.../plasmid_analyis` has no
such quota (7.2 PB free on `/gorilla/proj` at audit time). Root moved to
`plasmid_analyis/envs/.micromamba_root/` and all four envs were built there without issue.

## Two recurring gotchas (apply to every tool call below)

1. **`~/.local` pandas 3.0.4 shadows every conda env's own pandas.** A pre-existing user-level
   `pip install --user pandas` at `~/.local/lib/python3.1x/site-packages/pandas` sits ahead of (or
   gets merged ahead of) each conda env's site-packages on `sys.path`, because Python's user-site
   mechanism is enabled by default. This silently replaced each env's deliberately-pinned pandas
   with pandas 3.0.4, breaking `mob_typer` (`EmptyDataError` / `line_terminator` removed from the
   pandas 2.0+/3.0 API) and SegMantX (stricter pandas 3.0 dtype assignment rules). **Always run
   tools with `PYTHONNOUSERSITE=1`** to disable user-site lookup, or the env's own pinned pandas is
   not actually the one running.
2. **Conda/micromamba envs' own `bin/` must be prepended to `PATH`, not just invoked by absolute
   path.** `mob_typer`, `PlasAnn`, and SegMantX's `test_modules` all shell out to sibling binaries
   (`blastn`, `prodigal`, `mash`, or `SegMantX` itself) by bare name via `subprocess`/`shutil.which`,
   which only finds them if the env's `bin/` is on `$PATH` — invoking the entry-point script by
   absolute path alone is not sufficient.

## MOB-suite (`mob_typer`)

- Source: bioconda `mob_suite` package (upstream: Robertson & Nash 2018, github.com/phac-nml/mob-suite).
- Installed version: **3.1.9** (bioconda's current build at install time; the `mob_suite` package
  pulled into the shared `plasmid_phase0` solve had resolved to `1.4.9.1`, a much older build —
  isolating it in its own env let the solver pick the current 3.1.9 instead).
- Reference databases (NCBI plasmid finder DB, relaxase/MOB-family proteins, MPF proteins, oriT
  sequences, ~451MB) auto-download from Zenodo (`zenodo.org/records/10304948/files/data.tar.gz`)
  on first run via `mob_init` — this happened automatically as part of the bioconda post-link step
  and again as a runtime check; took ~70s once over this node's network.
- **Fix applied:** the shared-solve pandas (whatever the unconstrained solver picked) broke
  `mob_typer`'s `blast_df.to_csv(..., line_terminator=...)` call (`line_terminator` removed in
  pandas 2.0+). Re-solving `mob_suite_env` in isolation still produced pandas 3.0.4 via the
  `~/.local` shadow issue above — once `PYTHONNOUSERSITE=1` was applied, the env's own pinned
  pandas (1.5.3, which *does* still have `line_terminator`) took effect correctly.
- **Fix applied:** `mash` (used for nearest-neighbor/cluster typing) failed with
  `error while loading shared libraries: libgsl.so.25: cannot open shared object file` — the
  solver had picked `gsl 2.7.1` (providing `libgsl.so.27`) alongside a `mash` build linked against
  `libgsl.so.25`. Pinned `gsl=2.6` (provides `.so.25`), which also pulled in a matching `mash`
  rebuild (2.3, build `he348c14_1`). Without this fix `mob_typer` ran to completion but **silently
  wrote zero output rows** (no exception surfaced to the exit code) — this is worth remembering
  for any future "mob_typer completed but the TSV only has a header" symptom.
- **Verified on 5-plasmid test set** (`data/processed/test5/test5.fna`, real IMG/PR records,
  3.0-5.7 kb each): produced one row per contig with `rep_type(s)`, `relaxase_type(s)`, and
  `predicted_mobility` populated for 4/5 (1/5 had no replicon hit, `-`, which is expected for a
  reference-based tool on environmental plasmids — this untypeable fraction is exactly the
  caveat EXECUTION_PLAN.md Phase 1 calls out as a finding to report explicitly).
- Usage: `mob_typer --multi --infile <fasta> --out_file <tsv>` (use `--multi` to type every contig
  in a multi-FASTA as a separate sample; without it, `mob_typer` summarizes the whole input as one
  sample).

## PlasAnn

- Source: bioconda `plasann` package (upstream: github.com/ajlopatkin/PlasAnn, Lopatkin Lab).
- Installed version: **1.1.6**.
- Dependencies (BLAST+, Prodigal, Infernal) are **not** bundled by the bioconda recipe as
  hard pins that resolve automatically into the same env in all cases — `--check-deps` reported
  all three missing until `prodigal`, `blast`, and `infernal` were explicitly installed into
  `plasann_env`. Confirmed present and on `PATH` via `PlasAnn --check-deps` → "All external
  dependencies are installed!".
- **Capability surface** (from `--check-deps`/`--help` and a real run): gene prediction (Prodigal)
  + functional annotation (BLAST against a bundled reference set, optionally UniProt with
  `--uniprot-blast`), mobile-element detection (oriC/oriV, oriT, transposons, replicons), ncRNA
  detection (Infernal/Rfam), intergenic-region gene discovery, GenBank + annotation CSV + circular
  plasmid-map PNG output per plasmid.
- **Gotcha:** PlasAnn's single-file mode (`-i <file.fasta>`) requires the FASTA to contain exactly
  one record (`Bio.SeqIO.read` raises `ValueError: More than one record found in handle` on a
  multi-FASTA). For batches, use **folder mode**: split into one-file-per-plasmid first (e.g.
  `seqkit split -i <multi.fasta> -O <dir>`) and pass `-i <dir> -t fasta`.
- **Verified on 5-plasmid test set** (split into individual files, folder mode): 5/5 processed
  successfully in 49.3s, 5-11 annotated features per plasmid, GenBank + CSV + map PNG written per
  plasmid.

## SegMantX

- Source: cloned from `github.com/DMH-biodatasci/SegMantX` (upstream: Hanke & Dagan 2025, *MBE*
  42(10):msaf242) into `tools/SegMantX/` — **not** a packaged/installable tool, it is run in place
  as `tools/SegMantX/SegMantX <module> ...` with `tools/SegMantX/` on `$PATH` (its own
  `test_modules` and `chain_alignments` steps shell out to the sibling `SegMantX` entry script by
  bare name).
- Dependencies installed manually into `segmantx_env` per `SegMantX.yml` in the repo (pandas,
  scipy, biopython, plotly, streamlit, blast) — did **not** use `conda env create -f SegMantX.yml`
  directly because that file's `prefix: SegMantX` line and `anaconda::`-channel pins fight the
  no-system-conda / micromamba setup here; replicated the same package list manually instead.
- **Fix applied:** as installed, the unconstrained solve picked **pandas 3.0.3**, which broke
  `fetch_nucleotide_chains` (`TypeError: Invalid value '[ True False]' for dtype 'str'` — pandas
  3.0's stricter dtype assignment rejects assigning a bool array into what it now infers as a
  string-dtype column). Pinned `pandas<2.2` and `python<3.13` (re-solved to pandas 2.1.4), after
  which SegMantX's own `test_modules` self-test passed end-to-end (alignment generation, chaining,
  fasta extraction, HTML dotplot visualization, all 4 bundled real-genome test cases: "All modules
  are fine using the test data.").
- **Verified independently on our own 5-plasmid IMG/PR test set** (pipe-delimited
  `IMGPR_plasmid_...|taxon|scaffold` headers, not the tool's own NCBI-accession test set):
  `generate_alignments` (blastn self-alignment, circular topology) → `chain_self_alignments`
  (gap=5000bp) ran cleanly and produced a real chained duplication call (16 chained alignment
  hits, 96.4% identity, spanning ~27.5kb) confirming the tool handles IMG/PR's header format and
  ID scheme without modification.
- Core modules: `generate_alignments` (blastn-based local alignments, `-Q`/`-SA` for self-vs-self
  on circular sequences) → `chain_self_alignments` (chains local hits into duplication calls,
  `-G` max gap bp, `-SG` scaled gap) → `fetch_nucleotide_chains` (extract chained regions as
  FASTA) / `visualize_chains` (interactive HTML dotplot).

## Quick reference: invoking each tool

```bash
export MAMBA_ROOT_PREFIX=/gorilla/proj/h-mel-phylo/h-mel-phylo/private/amar/plasmid_analyis/envs/.micromamba_root

# MOB-suite
ENV=$MAMBA_ROOT_PREFIX/envs/mob_suite_env
PYTHONNOUSERSITE=1 PATH="$ENV/bin:$PATH" $ENV/bin/mob_typer --multi --infile <in.fna> --out_file <out.tsv>

# PlasAnn (folder mode, one file per plasmid)
ENV=$MAMBA_ROOT_PREFIX/envs/plasann_env
PYTHONNOUSERSITE=1 PATH="$ENV/bin:$PATH" $ENV/bin/PlasAnn -i <dir_of_single_fastas> -o <outdir> -t fasta

# SegMantX (run from tools/SegMantX/, repo dir must also be on PATH)
ENV=$MAMBA_ROOT_PREFIX/envs/segmantx_env
SEGDIR=/gorilla/proj/h-mel-phylo/h-mel-phylo/private/amar/plasmid_analyis/tools/SegMantX
PYTHONNOUSERSITE=1 PATH="$ENV/bin:$SEGDIR:$PATH" $ENV/bin/python $SEGDIR/SegMantX generate_alignments -q <in.fna> -b <out.blast.x7> -a <out.alignment_coordinates.tsv> -Q -SA
```
