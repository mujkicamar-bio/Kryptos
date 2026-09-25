#!/bin/bash -l
#SBATCH -A uppmax2025-2-42
#SBATCH -p pelle
#SBATCH -n 1
#SBATCH -c 96
#SBATCH --mem=600G
#SBATCH -t 10-00:00:00
#SBATCH -J dark_plasmidome
#SBATCH --dependency=singleton
#SBATCH -o results/logs/pipeline_%j.out
#
# The whole pipeline: 143,503 plasmids -> a complete, evidence-rich annotation of every
# plasmid ORF and its dark subset. It produces no ranking and no candidate list; spec
# section 76 puts experimental prioritisation in a separate downstream workflow.
#
# ONE SUBMISSION ANNOTATES EVERYTHING. `sbatch workflow/run_pipeline.sbatch` runs every
# stage in this one job: the cascade, the plasmid label databases and AMRFinderPlus
# (label_databases), DefenseFinder, CONJScan (conjugation_systems), IntegronFinder,
# ISEScan, synteny, the context terms, the structure search and the report. Calibrating
# the context terms (tools/calibrate_context.py) is a post-run tool and is not a rule.
#
# ENVIRONMENT
#
# envs/plasmidann is a materialised copy of workflow/envs/plasmidann.yaml and holds all
# eleven executables the workflow invokes. Putting it on PATH is what makes the conda:
# directives meaningful: an earlier version of this script activated nothing, so all 26
# directives were inert and the run silently depended on whatever the login shell had.
#
# To let Snakemake build the environment itself instead, drop the PATH line and add
# --use-conda. Do not do both: two environments means the pre-flight rule checks one and
# the searches run in the other, which is how a missing tool reaches hour 45 unnoticed.
#
# RESUMING
#
# Every stage is one job. Re-submitting this script resumes: --rerun-incomplete discards
# the partial output of the job that was killed and Snakemake rebuilds only what is
# missing, so a failure costs one stage rather than the run. DIAMOND cannot checkpoint, so
# a tier killed by the time limit restarts from zero: -t is the partition maximum (10 days)
# because the serial estimate is ~120-230 h.
#
# --dependency=singleton, with the same job name in structure_gpu.sbatch: both invocations
# below run --nolock, so a second submission would otherwise write the same outputs at the
# same time. SLURM now starts it only when the first has ended.
#
# STRUCTURE SEARCH (S8d) runs HERE, on CPU, by default (~15-35 h on 16 CPU threads for the
# expected dark representatives), so that one submission produces the whole annotation.
# ProstT5 is far faster on a GPU, so a GPU split remains available as an option:
#
#   STRUCTURE_ON_GPU=1 sbatch workflow/run_pipeline.sbatch   # stops before structure_search
#   sbatch workflow/structure_gpu.sbatch                     # S8d on the GPU partition
#   sbatch workflow/run_pipeline.sbatch                      # the rest, resuming
#
# With STRUCTURE_ON_GPU=1 the job stops before structure_search (--omit-from) until
# structure_hits.tsv exists. The split trades queue waits for CPU hours; the default
# trades CPU hours for a single run, which is what was asked for.
#
# TOOLS IN THEIR OWN ENVIRONMENTS. AMRFinderPlus (envs/amrfinder) and the MacSyFinder 2.1.6
# that runs CONJScan (envs/conjscan) are not on PATH: config names them by path
# (amrfinder.executable, conjugation.exe), as it names pharokka, and pre-flight checks them
# there - including that the MacSyFinder is new enough for CONJScan's model grammar.
set -euo pipefail

# MEMORY. Every rule declares mem_mb, and in local execution Snakemake ignores all of it
# unless a ceiling is supplied here: it schedules on the core budget alone. Without this the
# 600 GB allocation is respected by luck rather than by accounting - orf_call, integrons and
# the nr tier can all become schedulable at once. 550 GB of the 600 leaves
# headroom for the OS and for tool temporary files.
MEM_MB=550000
CORES=${SLURM_CPUS_PER_TASK:-96}
OMIT=()
if [ -n "${STRUCTURE_ON_GPU:-}" ] && [ ! -e results/12_context_and_structure/structure_hits.tsv ]; then
    OMIT=(--omit-from structure_search)
fi

# The repository root is derived from this script's own location, so the pipeline runs
# from wherever it is checked out rather than from one person's home directory.
ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
cd "$ROOT"

# WHICH PYTHON RUNS THE SCRIPTS
#
# Snakemake executes a `script:` directive with the interpreter that is running Snakemake,
# NOT with whatever is first on PATH. Invoking a Snakemake from one environment while
# putting a different environment on PATH therefore splits the run in two: the shell tools
# (diamond, hmmsearch, pharokka) resolve through PATH, and every Python script runs
# somewhere else. That split is silent until a script imports something only one of the
# two environments has - it surfaced as ModuleNotFoundError: scipy in context_features,
# after 34 of 35 rules had already succeeded.
#
# So the environment that holds the tools also runs Snakemake, and PATH agrees with it.
SNAKEMAKE="$ROOT/envs/plasmidann/bin/snakemake"
export PATH="$ROOT/envs/plasmidann/bin:$PATH"
# Tool caches off the home directory, whose file quota is 94% used.
export XDG_CACHE_HOME="${XDG_CACHE_HOME:-$ROOT/.cache}"

if [ ! -x "$SNAKEMAKE" ]; then
    echo "no snakemake at $SNAKEMAKE - see README.md, section Installation" >&2
    exit 1
fi

# Pre-flight runs first in its own invocation so a missing tool or database fails in
# seconds, before the whole node is spent on gene calling.
#
# The target is placed AFTER --nolock deliberately. Both --resources and --configfile
# accept one or more values, so a target written after either is consumed as one of them:
# after --resources it fails with "dictionary update sequence element #1 has length 1",
# and after --configfile it is opened as a second config file. --nolock is a boolean, so
# a positional argument following it is parsed as a target.
"$SNAKEMAKE" -s workflow/Snakefile \
  --resources mem_mb=$MEM_MB -j 1 --nolock \
  results/05_annotation_cascade/preflight.tsv

"$SNAKEMAKE" -s workflow/Snakefile \
  -j "$CORES" \
  --resources mem_mb=$MEM_MB \
  --rerun-incomplete \
  --nolock \
  --keep-going \
  --printshellcmds \
  --rerun-triggers mtime \
  "${OMIT[@]}"
