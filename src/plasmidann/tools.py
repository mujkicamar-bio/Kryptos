"""Every external executable the pipeline invokes, and which stage needs it.

WHY THIS FILE EXISTS

v1's most expensive failure was a job that ran for 45 hours and then died because DIAMOND
was not on PATH. The identical failure happened twice. The pre-flight rule was written to
end that class of failure - and then it checked only the tools named in the cascade tier
list, four of the eleven the workflow actually runs. A missing mafft still killed S7b after
S3 had run for days; a missing prodigal made every one of the 600 IntegronFinder shards
fail; a missing foldseek emptied the novel_fold stratum.

The registry is therefore declared once, here, and pre-flight walks all of it. A test in
tests/test_tools_registry.py scans workflow/scripts/ for subprocess invocations and fails
if a script runs anything this list does not name, so the list cannot silently fall behind
the code again.

`optional_when` names a config predicate: foldseek is only required when structural
evidence is required, because a run without it is a declared, honest degradation rather
than a broken run.
"""

REQUIRED_TOOLS = [
    {"name": "hmmsearch", "stage": "S2b, S3",
     "why": "AntiFam artefact screen and every hmmer cascade tier"},
    {"name": "diamond", "stage": "S3",
     "why": "the Swiss-Prot and nr cascade tiers"},
    {"name": "hhblits_omp", "stage": "Stage 4 tier 3 (PHROGs)",
     "why": "the PHROGs tier; PHROGs ships HH-suite HHM profiles and documents HH-suite "
            "as the way to search them, so nothing else can read that database"},
    {"name": "ffindex_from_fasta", "stage": "Stage 4 tier 3 (PHROGs)",
     "why": "hhblits_omp reads its queries as an ffindex, not as a FASTA"},
    {"name": "tantan", "stage": "S2b",
     "why": "low-complexity masking; without it every protein reports 0 masked"},
    {"name": "mash", "stage": "Stage 6",
     "why": "plasmid lineage clustering - sequence independence, separate from MOB class"},
    {"name": "mmseqs", "stage": "S6b",
     "why": "deep-homology clustering into families; nothing downstream has families"},
    {"name": "mafft", "stage": "S7b",
     "why": "codon alignments for dN/dS; without it every family reports ALIGNMENT_FAILED"},
    {"name": "defense-finder", "stage": "S8a phase 1",
     "why": "defence component hits; without it the defence_island stratum is empty"},
    {"name": "macsyfinder", "stage": "S8a phase 2",
     "why": "calls systems from gene adjacency; components alone are not systems"},
    {"name": "integron_finder", "stage": "S8b",
     "why": "cassette arrays; without it the integron_cassette stratum is empty"},
    {"name": "prodigal", "stage": "S8b (invoked by integron_finder)",
     "why": "IntegronFinder calls genes with prodigal and exits non-zero without it"},
    {"name": "cmsearch", "stage": "S8b (invoked by integron_finder)",
     "why": "IntegronFinder locates attC sites with an Infernal covariance model"},
    {"name": "RNAcode", "stage": "S7b",
     "why": "coding-potential signal independent of the gene caller, on both strands"},
    {"name": "emapper.py", "stage": "S4b",
     "why": "COG and KEGG terms for the annotated fraction; S8's pathway axis needs them",
     "optional_when": "orthology.required is false"},
    {"name": "foldseek", "stage": "S8d",
     "why": "structural homology; without it the folds evidence line can never fire",
     "optional_when": "structure.required is false"},
]


def tool_names():
    """The set of executable names, for membership checks."""
    return {t["name"] for t in REQUIRED_TOOLS}


# Which config predicate governs each conditional tool.
_CONDITIONAL = {"foldseek": "structure", "emapper.py": "orthology"}


def required_tools(structure_required=True, orthology_required=True):
    """The registry entries that must be present for THIS run's configuration.

    Two tools are conditional, and both follow the same contract: when the stage is
    declared required a missing tool or database fails pre-flight in seconds, and when it
    is not, the stage is skipped and its evidence is recorded as ABSENT rather than as
    searched and not found. Everything else is unconditional, because a pipeline that
    quietly skips one of the other stages produces a table built on missing evidence rather
    than on absent evidence, and those are not the same claim.
    """
    enabled = {"structure": structure_required, "orthology": orthology_required}
    return [t for t in REQUIRED_TOOLS
            if "optional_when" not in t or enabled[_CONDITIONAL[t["name"]]]]
