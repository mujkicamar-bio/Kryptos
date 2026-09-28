"""Context terms (rule context_features, family_context_terms.tsv): what a family's
neighbours are, counted per independent lineage.

WHAT A TERM IS

A term names one thing an ORF's neighbourhood holds, as '<type>:<name>':

    amr:<label>          CARD protein homolog model (card_amr_family), or an AMRFinderPlus
                         element of type AMR (amrfinder_gene)
    metal:<label>        BacMet 2.0 experimentally confirmed (bacmet_compound), or an
                         AMRFinderPlus STRESS element of subtype METAL or BIOCIDE
    ta:<label>           TADB 3.0 validated (tadb_ta)
    conj_role:<label>    oriTDB 2.0 relaxase / auxiliary / T4CP (oritdb_role)
    mge:<label>          mobileOG-db major category, every evidence class (mobileog_category)
    antidefence:<label>  dbAPIS (dbapis_family) and Anti-CRISPRdb (acrdb_family), every
                         evidence class; the class is in the label's sub_label
    defence:<system>     DefenseFinder system
    conj:<system>        CONJScan system

Labels are used verbatim, as the label step wrote them. KEGG is not a context source, and
pharokka's CARD rows are not either: amr: comes from the CARD search of the label step.

TWO NEIGHBOUR RULES

A gene label counts only from a neighbour within +-3 genes that is in the same directon as
the ORF - same strand, intergenic gaps of at most context.max_operon_gap - which is the
FESNov rule (Rodriguez del Rio et al. 2024, Nature 626:377). A system is a multi-gene call
whose components sit on both strands, so a system term counts when the ORF is itself a
component or a component lies within +-3 genes on either strand.

A neighbour that belongs to the focal family is a tandem paralogue: its label says what
the family is, not what surrounds it, so it is excluded. The ORF's own system membership
is never excluded - it is where the ORF sits, not a neighbour.

THE UNIT IS THE LINEAGE

Clonal copies of one plasmid are one observation. A family's conservation of a term is the
fraction of the plasmid lineages (plasmid_lineage.tsv) it occurs in where at least one
occurrence has the term.
Below two lineages there is nothing to conserve across: the row has status
TOO_FEW_LINEAGES and no conservation.
"""
import re

from darkorf.status import SUCCESS, TOO_FEW_LINEAGES

# The kind of a label (plasmidann.labeldb) and the term type it gives. amrfinder_gene is
# absent on purpose: its type depends on the element type, see label_term.
TERM_PREFIX = {
    "card_amr_family": "amr",
    "tadb_ta": "ta",
    "bacmet_compound": "metal",
    "oritdb_role": "conj_role",
    "mobileog_category": "mge",
    "dbapis_family": "antidefence",
    "acrdb_family": "antidefence",
}
AMRFINDER_KIND = "amrfinder_gene"

# Conservation across one lineage is a single observation, not a conservation.
MIN_LINEAGES = 2

COLUMNS = ["family_id", "term_type", "term", "n_lineages",
           "n_lineages_with_term", "conservation", "status",
           "window_covers_plasmid_fraction"]


def label_term(kind, label, sub_label=""):
    """The context term one protein label gives, or None when the kind gives none.

    AMRFinderPlus reports an element type (AMR, STRESS, VIRULENCE) and subtype (for STRESS:
    METAL, BIOCIDE, ACID, HEAT ...), which the label step writes into sub_label. AMR
    elements are amr:, STRESS METAL and BIOCIDE are metal:, and the rest have no term. A
    row without a recognisable element type raises: guessing would file a mercury
    reductase under amr:.
    """
    if kind in TERM_PREFIX:
        return f"{TERM_PREFIX[kind]}:{label}"
    if kind != AMRFINDER_KIND:
        return None
    tokens = set(re.split(r"[^A-Z]+", (sub_label or "").upper()))
    if "AMR" in tokens:
        return f"amr:{label}"
    if "STRESS" in tokens:
        return f"metal:{label}" if tokens & {"METAL", "BIOCIDE"} else None
    if "VIRULENCE" in tokens:
        return None
    raise ValueError(
        f"AMRFinderPlus label {label!r} has sub_label {sub_label!r}, which names no "
        "element type (AMR, STRESS or VIRULENCE); the term type cannot be decided.")


def system_term(prefix, system):
    """'defence:Clover' from a MacSyFinder model path or a bare model name."""
    return f"{prefix}:{system.rsplit('/', 1)[-1]}"


def orf_term_sources(orf_id, near, directon, label_terms, system_terms):
    """Every (term, source ORF) in one ORF's context; the source is None for its own system.

    `near` holds the orf_ids within the +-window neighbourhood (plasmidann.context.flanks),
    `directon` the ORF's transcriptional unit (plasmidann.context.directons), and
    `label_terms` and `system_terms` map an orf_id to the terms it carries. The source is
    kept so that the family step can drop tandem paralogues, which only it can recognise.
    """
    near = set(near)
    out = {(t, None) for t in system_terms.get(orf_id, ())}
    for n in near:
        out |= {(t, n) for t in system_terms.get(n, ())}
    for n in near & set(directon):
        out |= {(t, n) for t in label_terms.get(n, ())}
    return sorted(out, key=lambda s: (s[0], s[1] or ""))


def window_covers_plasmid(near, n_genes):
    """Whether the neighbourhood `near` holds every other gene of a plasmid of `n_genes`.

    On such a plasmid everything is everything's neighbour, so a context term there says
    little about the ORF; the family rows report how often that is the case.
    """
    return len(set(near)) == n_genes - 1


def family_term_rows(family_id, occurrences, family_of_orf, lineage_of):
    """One row per term in a family's context, in COLUMNS order.

    `occurrences` holds (orf_id, sources, covers) for every ORF of the family's members,
    with sources from orf_term_sources and covers from window_covers_plasmid. orf_id is
    '<plasmid_id>|<ordinal>'. A plasmid missing from `lineage_of` raises KeyError: rule
    plasmid_lineage assigns every analysed plasmid a lineage, so a gap is an input mismatch.

    conservation is counted over lineages, and is empty under TOO_FEW_LINEAGES;
    window_covers_plasmid_fraction is the fraction of the family's occurrences (ORFs)
    whose window covers their plasmid.
    """
    lineages = set()
    with_term = {}
    n_covers = 0
    for orf_id, sources, covers in occurrences:
        plasmid_id = orf_id.rsplit("|", 1)[0]
        if plasmid_id not in lineage_of:
            raise KeyError(f"plasmid {plasmid_id} of family {family_id} has no "
                           "plasmid lineage")
        lineage = lineage_of[plasmid_id]
        lineages.add(lineage)
        n_covers += bool(covers)
        for term, source in sources:
            if source is not None and family_of_orf.get(source) == family_id:
                continue
            with_term.setdefault(term, set()).add(lineage)
    n = len(lineages)
    covers_fraction = round(n_covers / len(occurrences), 6) if occurrences else 0.0
    measured = n >= MIN_LINEAGES
    return [{
        "family_id": family_id,
        "term_type": term.split(":", 1)[0], "term": term,
        "n_lineages": n, "n_lineages_with_term": len(with_term[term]),
        "conservation": round(len(with_term[term]) / n, 6) if measured else "",
        "status": SUCCESS if measured else TOO_FEW_LINEAGES,
        "window_covers_plasmid_fraction": covers_fraction,
    } for term in sorted(with_term)]
