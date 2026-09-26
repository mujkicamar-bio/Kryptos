"""The host of a plasmid, from whichever of three sources records it.

Sources in order of preference, with the plasmids of the 143,503-plasmid analysis set for
which each names an organism under the rules below:

  PLSDB species (master table `plsdb_species`)                           45,642 (31.8%)
  PlasmidScope per-record host (provenance/complete_provenance.tsv, `host`), the host
      each source database reports, IMG/PR's among them                  26,214 (18.3%)
  GenBank/RefSeq source organism (provenance/working_set.tsv, `organism`)  70,421 (49.1%)

Together they name a host for 87,567 plasmids (61.0%): 98.4% of the isolate plasmids but
only 23.9% of the metagenomic ones. Where PLSDB and the organism field both name one they
agree on the genus for 45,230 of 45,640 (99.1%).

A name is kept only if it names an organism: "human gut metagenome", "uncultured
bacterium", "Plasmid pDS56" and "Enterobacteriaceae bacterium" do not. The species is the
binomial; a name resolved only to the genus ("Acidovorax sp. DW039") counts toward the
genus but not the species.
"""
import re

_NOT_AN_ORGANISM = re.compile(
    r"^(NA|none|unknown|unidentified|uncultured|unclassified|environmental|synthetic|"
    r"plasmid|cloning|vector|expression|shuttle|mixed|bacterium|bacteria|prokaryote|"
    r"archaeon)\b|metagenome", re.I)
_NOT_A_SPECIES = {"sp.", "sp", "spp.", "cf.", "genomosp."}
# A genus is one capitalised Latin word; "[Clostridium]" marks a misplaced genus.
_GENUS = re.compile(r"^\[?[A-Z][a-z]+\]?$")
# Suffixes of the family, order and phylum names (International Code of Nomenclature of
# Prokaryotes, Rule 8: Parker et al. 2019, IJSEM 69:S1; Oren & Garrity 2021, IJSEM
# 71:005056 for the phylum).
_HIGHER_RANK = ("aceae", "ales", "ota")
# The noun that follows a higher taxon in "Mollicutes bacterium", "Nostocales cyanobacterium".
_UNNAMED_MEMBER = {"bacterium", "cyanobacterium", "archaeon"}


def normalise(name):
    """(species, genus) for one recorded host name; ("", "") when it names no organism."""
    name = (name or "").replace("_", " ").strip()
    if not re.search(r"[A-Za-z]", name) or _NOT_AN_ORGANISM.search(name):
        return "", ""
    words = name.split()
    prefix = ""
    if words[0] == "Candidatus" and len(words) > 1:
        prefix, words = "Candidatus ", words[1:]
    # "Enterobacteriaceae bacterium" names a family, class or order, not a genus, and a
    # free-text first word ("endosymbiont", "'Brassica") names no genus.
    if (not _GENUS.match(words[0]) or words[0].endswith(_HIGHER_RANK)
            or words[1:2] and words[1] in _UNNAMED_MEMBER):
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
