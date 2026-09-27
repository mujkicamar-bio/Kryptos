"""Every external executable the pipeline finds on PATH, and which stage needs it.

Pre-flight checks every entry before any search. An entry with `only_if_required` names a
config section (`structure`, `orthology`): the tool is needed only when that section's
`required` is true, and otherwise the stage is skipped and recorded as NOT_RUN. Tools that
live in their own environment - pharokka, AMRFinderPlus and the MacSyFinder that runs
CONJScan - are named by path in config and checked by path in pre-flight, so they are not
listed here.
"""
import pathlib
import re

REQUIRED_TOOLS = [
    {"name": "hmmsearch", "stage": "S2b, S3",
     "why": "AntiFam artefact screen and every hmmer cascade tier"},
    {"name": "diamond", "stage": "S3, S4d",
     "why": "the Swiss-Prot and nr cascade tiers and the plasmid label databases"},
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
    {"name": "isescan.py", "stage": "S8e",
     "why": "IS elements; without it no ORF can be placed inside or beside an IS element"},
    {"name": "RNAcode", "stage": "S7b",
     "why": "coding-potential signal independent of the gene caller, on both strands"},
    {"name": "emapper.py", "stage": "S4b",
     "why": "COG and KEGG terms for the annotated fraction; S8's pathway axis needs them",
     "only_if_required": "orthology"},
    {"name": "foldseek", "stage": "S8d",
     "why": "structural homology; without it the folds evidence line can never fire",
     "only_if_required": "structure"},
]


def tool_names():
    """The set of executable names, for membership checks."""
    return {t["name"] for t in REQUIRED_TOOLS}


def required_tools(structure_required=True, orthology_required=True):
    """The registry entries that must be present for this run's configuration.

    Every other stage is unconditional: skipping it quietly would give a table built on
    missing evidence rather than on evidence recorded as absent.
    """
    enabled = {"structure": structure_required, "orthology": orthology_required}
    return [t for t in REQUIRED_TOOLS
            if "only_if_required" not in t or enabled[t["only_if_required"]]]


# ------------------------------------------------------------------------------------
# MacSyFinder model grammar. CONJScan 2.0.2 and later write their definitions in grammar
# 2.1 (the `vers` attribute of <model>), which MacSyFinder reads from 2.1.6 with MacSyLib
# 1.0.4 (CONJScan README, changelog 2.0.2). MacSyFinder 2.1.4 - the version DefenseFinder
# pins - stops with "has not the right version. version supported is '2.0'".
# ------------------------------------------------------------------------------------
MIN_MACSYFINDER = {"2.0": (2, 0), "2.1": (2, 1, 6)}


def macsyfinder_version(text):
    """(major, minor, patch) from `macsyfinder --version` output, None when absent.

    2.1.6 prints 'MacSyFinder 2.1.6', 2.1.4 prints 'Macsyfinder 2.1.4'.
    """
    match = re.search(r"macsyfinder\s+(\d+(?:\.\d+)+)", text, re.IGNORECASE)
    return tuple(int(x) for x in match.group(1).split(".")) if match else None


def model_grammars(models_dir):
    """Every grammar version declared by a model definition under `models_dir`."""
    grammars = set()
    for xml in pathlib.Path(models_dir).rglob("definitions/**/*.xml"):
        match = re.search(r"<model\b[^>]*\bvers=\"([^\"]+)\"", xml.read_text())
        if match:
            grammars.add(match.group(1))
    return grammars


def grammar_problem(grammars, version):
    """Why a MacSyFinder `version` cannot read definitions of `grammars`, or ''."""
    unknown = sorted(g for g in grammars if g not in MIN_MACSYFINDER)
    if unknown:
        return f"model grammar {', '.join(unknown)} is not known to plasmidann.tools"
    need = max(MIN_MACSYFINDER[g] for g in grammars)
    if version < need:
        return (f"the models use grammar {max(grammars)}, which needs MacSyFinder >= "
                f"{'.'.join(map(str, need))}, and the executable is MacSyFinder "
                f"{'.'.join(map(str, version))}")
    return ""
