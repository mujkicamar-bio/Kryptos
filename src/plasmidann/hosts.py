"""The host of a plasmid, from whichever of three sources records it.

Sources in order of preference:

  PLSDB species (master table `plsdb_species`)
  PlasmidScope per-record host (provenance/complete_provenance.tsv, `host`), the host
      each source database reports, IMG/PR's among them
  GenBank/RefSeq source organism (provenance/working_set.tsv, `organism`)

A name is kept only if it names an organism: "human gut metagenome", "uncultured
bacterium", "Plasmid pDS56", "Tomato aspermy virus" and "Enterobacteriaceae bacterium" do
not. The species is the binomial; a name resolved only to the genus ("Acidovorax sp.
DW039", "Wolbachia endosymbiont of Drosophila") counts toward the genus but not the
species. A phytoplasma recorded by its common name ("Onion yellows phytoplasma") is of the
genus Candidatus Phytoplasma.
"""
import re

_NOT_AN_ORGANISM = re.compile(
    r"^(NA|none|unknown|unidentified|uncultured|unclassified|environmental|synthetic|"
    r"cloning|vector|expression|shuttle|mixed|bacterium|bacteria|prokaryote|"
    r"archaeon)\b|metagenome|\bplasmid\b|virus\b", re.I)
_NOT_A_SPECIES = {"sp.", "sp", "spp.", "cf.", "genomosp.", "endosymbiont", "symbiont"}
# A genus is one capitalised Latin word; "[Clostridium]" marks a misplaced genus.
_GENUS = re.compile(r"^(\[[A-Z][a-z]+\]|[A-Z][a-z]+)$")
# Suffixes of the family, order and phylum names (International Code of Nomenclature of
# Prokaryotes, Rule 8: Parker et al. 2019, IJSEM 69:S1; Oren & Garrity 2021, IJSEM
# 71:005056 for the phylum).
_HIGHER_RANK = ("aceae", "ales", "ota")
# The noun that follows a higher taxon in "Mollicutes bacterium", "Nostocales cyanobacterium".
_UNNAMED_MEMBER = {"bacterium", "cyanobacterium", "archaeon"}


def normalise(name):
    """(species, genus) for one recorded host name; ("", "") when it names no organism.

    A name that lists several hosts separated by commas gives the first that names one.
    """
    for part in (name or "").split(","):
        species, genus = _normalise_one(part)
        if genus:
            return species, genus
    return "", ""


def _normalise_one(name):
    name = name.replace("_", " ").strip()
    if not re.search(r"[A-Za-z]", name) or _NOT_AN_ORGANISM.search(name):
        return "", ""
    if "phytoplasma" in name.lower() and not name.startswith("Candidatus"):
        return "", "Candidatus Phytoplasma"
    words = name.split()
    prefix = ""
    if words[0] == "Candidatus" and len(words) > 1:
        prefix, words = "Candidatus ", words[1:]
    # "Enterobacteriaceae bacterium" names a family, class or order, not a genus; a
    # free-text first word ("endosymbiont", "'Brassica") names no genus; and in "Plautia
    # stali symbiont" the binomial is the symbiont's animal host.
    if (not _GENUS.match(words[0]) or words[0].endswith(_HIGHER_RANK)
            or words[1:2] and words[1] in _UNNAMED_MEMBER
            or len(words) > 2 and words[-1] in ("symbiont", "endosymbiont")):
        return "", ""
    genus = prefix + words[0]
    if len(words) < 2 or words[1] in _NOT_A_SPECIES:
        return "", genus
    return f"{genus} {words[1]}", genus


def resolve(candidates):
    """The first source that names an organism: (species, genus, source).

    `candidates` is [(source, name), ...] in order of preference.
    """
    for source, name in candidates:
        species, genus = normalise(name)
        if genus:
            return species, genus, source
    return "", "", ""
