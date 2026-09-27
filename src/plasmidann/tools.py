"""Every external executable the pipeline finds on PATH, and the rules that run it.

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
    {"name": "hmmsearch", "stage": "artefact_screen, tier_search, consensus_recheck",
     "why": "AntiFam artefact screen and every hmmer cascade tier"},
    {"name": "diamond", "stage": "tier_search, label_databases",
     "why": "the Swiss-Prot and nr cascade tiers and the plasmid label databases"},
    {"name": "tantan", "stage": "artefact_screen",
     "why": "low-complexity masking; without it every protein reports 0 masked"},
    {"name": "mash", "stage": "plasmid_lineage",
     "why": "plasmid lineage clustering - sequence independence, separate from MOB class"},
    {"name": "mmseqs", "stage": "protein_clustering, cascade_selection, family_network",
     "why": "deep-homology clustering into families; nothing downstream has families"},
    {"name": "mafft", "stage": "family_evolution",
     "why": "codon alignments for dN/dS; without it every family reports ALIGNMENT_FAILED"},
    {"name": "defense-finder", "stage": "defence_search",
     "why": "defence component hits; without it the defence_island stratum is empty"},
    {"name": "macsyfinder", "stage": "defence_systems",
     "why": "calls systems from gene adjacency; components alone are not systems"},
    {"name": "integron_finder", "stage": "integrons",
     "why": "cassette arrays; without it the integron_cassette stratum is empty"},
    {"name": "prodigal", "stage": "integrons (invoked by integron_finder)",
     "why": "IntegronFinder calls genes with prodigal and exits non-zero without it"},
    {"name": "cmsearch", "stage": "integrons (invoked by integron_finder)",
     "why": "IntegronFinder locates attC sites with an Infernal covariance model"},
    {"name": "isescan.py", "stage": "is_elements",
     "why": "IS elements; without it no ORF can be placed inside or beside an IS element"},
    {"name": "yn00", "stage": "family_evolution",
     "why": "PAML yn00 pairwise dN/dS; without it every family reports YN00_FAILED"},
    {"name": "RNAcode", "stage": "family_evolution",
     "why": "coding-potential signal independent of the gene caller, on both strands"},
    {"name": "emapper.py", "stage": "orthology",
     "why": "COG and KEGG terms for the annotated fraction; dark ORF neighbourhoods are described with them",
     "only_if_required": "orthology"},
    {"name": "foldseek", "stage": "structure_search",
     "why": "structural homology; without it the folds evidence line can never fire",
     "only_if_required": "structure"},
]


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
