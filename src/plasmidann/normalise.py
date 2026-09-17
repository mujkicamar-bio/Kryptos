"""Collapsing free-text product names to a canonical form.

WHY ONLY THE FREE-TEXT KINDS NEED THIS

Pfam family names, COG identifiers, GO terms, EC numbers and KO identifiers are controlled
vocabularies: a finite set of exact strings that can be enumerated and reviewed. Product
names from nr and Swiss-Prot are not. The same protein is deposited as 'Plasmid replication
initiator protein RepA', 'replication initiation protein' and 'putative replication
initiator protein' by three different submitters.

Left as written, one concept becomes dozens of categories, each with too few members to
test, and the frequency ranking that makes the review tractable does not work. Normalised,
the surface forms collapse and the head of the distribution is a bounded list.

WHAT IS DELIBERATELY NOT DONE HERE

No stemming, no synonym dictionary, no fuzzy matching. Those would merge products that a
reviewer should see separately, and the merge would be invisible in the output. This
function only removes tokens that carry no functional information: the submitter's hedging,
the generic suffixes, and inconsistent whitespace and case.
"""
import re

# Confidence hedges. They describe the submitter's certainty, not the protein, and keeping
# them splits every product into a confident and a hedged form.
_PREFIXES = ("putative ", "probable ", "predicted ", "possible ", "conserved ")

# Generic suffixes applied inconsistently across submissions. 'MobA/MobL family protein'
# and 'MobA/MobL protein' are one product.
#
# Longest first, so 'domain-containing protein' is removed whole rather than leaving a
# stranded 'domain-containing' behind the shorter ' protein' match.
_SUFFIXES = (" domain-containing protein", " domain containing protein",
             " superfamily protein", " family protein", " containing protein",
             " like protein", " protein")


def product_name(text):
    """Canonical form of a free-text product name.

    Lower-cased, hedges and generic suffixes removed, whitespace collapsed. Stripping never
    empties the string: a product that is nothing but a suffix is returned whole, because
    an empty result would merge it with every parse failure in the collection.
    """
    value = re.sub(r"\s+", " ", (text or "").strip().lower())
    if not value:
        return ""

    changed = True
    while changed:
        changed = False
        for prefix in _PREFIXES:
            if value.startswith(prefix) and len(value) > len(prefix):
                value, changed = value[len(prefix):].strip(), True
        for suffix in _SUFFIXES:
            if value.endswith(suffix) and len(value) > len(suffix):
                value, changed = value[:-len(suffix)].strip(), True
    return value
