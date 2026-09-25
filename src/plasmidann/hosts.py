"""The host of a plasmid, from whichever of three sources records it.

PLSDB species names a host for 45,733 of the 143,503 analysis-set plasmids (31.9%). Two more
sources add to it (measured 2026-09-24, after the filter below):

  PlasmidScope per-record host (provenance/complete_provenance.tsv, `host`), which carries
      the host each source database reports - IMG/PR's among them        26,267 (18.3%)
  GenBank/RefSeq source organism (provenance/working_set.tsv, `organism`)  70,571 (49.2%)

Taken in that order of preference, 87,718 plasmids (61.1%) have a host: 98.6% of the
isolate plasmids but only 24.0% of the metagenomic ones, most of which no source assigns
to a host. Where PLSDB and the organism field both exist they agree on the genus for
45,313 of 45,733 (99.1%).

A name is kept only if it names an organism: "human gut metagenome", "uncultured
bacterium" and "Plasmid pDS56" do not, and would each count as a host. The species is the
binomial; a name resolved only to the genus ("Acidovorax sp. DW039") counts toward the
genus but not the species.
"""
import re

_NOT_AN_ORGANISM = re.compile(
    r"^(NA|none|unknown|unidentified|uncultured|unclassified|environmental|synthetic|"
    r"plasmid|cloning|vector|expression|shuttle|mixed|bacterium|bacteria|prokaryote|"
    r"archaeon)\b|metagenome", re.I)
_NOT_A_SPECIES = {"sp.", "sp", "spp.", "cf.", "bacterium", "genomosp."}


def normalise(name):
    """(species, genus) for one recorded host name; ("", "") when it names no organism."""
    name = (name or "").replace("_", " ").strip()
    if not re.search(r"[A-Za-z]", name) or _NOT_AN_ORGANISM.search(name):
        return "", ""
    words = name.split()
    prefix = ""
    if words[0] == "Candidatus" and len(words) > 1:
        prefix, words = "Candidatus ", words[1:]
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
