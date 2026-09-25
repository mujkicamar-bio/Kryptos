"""Layer B evidence description: what the record SAYS about a family.

This module describes evidence. It does not select on it. Spec section 76 draws that
boundary: the core pipeline produces no top_1000, no candidate_score, no novelty_score and
no experimental_rank, and section 2.3 forbids collapsing evidence into a single number
because evidence must remain named and interpretable.

The selection machinery that used to live here - eligibility, a Pareto front, per-stratum
quotas, a synthesis order - is deleted with the rest of Layer C. Section 77 lists what a
downstream prioritisation workflow may use, and every item is a column this pipeline
produces. The decision is made on the tables, not inside them.

WHAT REMAINS, AND WHY IT IS DESCRIPTION RATHER THAN SELECTION

  reality_lines     COUNTED and NAMED. Independent lines of evidence that a family is a
                    real protein rather than a gene-calling artefact. The names matter as
                    much as the count: "3 of 4, missing purifying_selection" tells a reader
                    what to go and check, which a score never can.

                    Absence of evidence never counts against a family. A dN/dS status of
                    NO_DIVERGENCE simply fails to fire a line; it is not read as evidence
                    of neutrality, because identical sequences are what a highly conserved
                    family looks like.

  darkness_state    NAMED, two states. A dark protein whose fold is recognisable and one
                    whose fold is not are different situations, and averaging them into a
                    number would discard the distinction that matters most.

Both statements are made ABOUT a family and recorded beside it. Neither ranks anything.
"""


# The four independent lines of evidence that a dark ORF is a real protein. Each is a
# named boolean with a stated rationale, so a reader can disagree with any one of them
# specifically rather than with an opaque number.
#
# EVERY THRESHOLD ARRIVES AS AN ARGUMENT. They were literals inside these lambdas - 0.5, 2
# and 3 - and two of them duplicated values that already lived in config/targets.yaml and
# were read correctly by S7b. A change to either config value would have silently diverged
# from what S9 eligibility actually applied: no error, no schema violation, and no entry in
# the provenance ledger to notice by. That is the same defect the cascade fixed once
# already, where three of six thresholds were module constants outside config and schema.
REALITY_TESTS = [
    {
        "name": "purifying_selection",
        # The strongest available evidence, and the only one that catches a shadow ORF -
        # which is conserved and multi-species, but whose locus is under selection for the
        # gene on the OTHER strand.
        "test": lambda f, thr: (f.get("dnds_status") == "MEASURED"
                                and _num(f.get("dnds_median")) is not None
                                and _num(f.get("dnds_median")) < thr["dnds_purifying_max"]),
    },
    {
        "name": "multi_lineage",
        # Present on two or more Stage 6 lineages (plasmid_lineage, Mash), so not a
        # single-lineage accident: one clone sequenced forty times is one lineage. Not MOB-suite
        # clusters, which assign the nearest reference however distant and so put unrelated
        # novel plasmids into one cluster.
        "test": lambda f, thr: (_int(f.get("independent_plasmid_cluster_count"))
                                >= thr["min_lineages"]),
    },
    {
        "name": "is_family",
        # FESNov's minimum, and enough to align. Deliberately the SAME number as
        # evolution.min_members_for_dnds - see check_reality_config.
        "test": lambda f, thr: _int(f.get("n_members")) >= thr["min_members"],
    },
    {
        "name": "folds",
        # A recognisable structure means there is something for the sequence to be.
        "test": lambda f, thr: bool(f.get("structural_match")),
    },
]


def _num(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def _int(v):
    try:
        return int(float(v))
    except (TypeError, ValueError):
        return 0


# Which lines are NOT independent of which others.
#
# `purifying_selection` requires a MEASURED dN/dS, and dN/dS is only computed for families
# with at least min_members_for_dnds (3) members - which is exactly the `is_family` test.
# One strictly implies the other, so a candidate firing only those two was reported as
# carrying "2 of 4 independent lines" and cleared an eligibility floor of 2 on a single
# line of evidence. The whole library could have been admitted that way, which is the
# opposite of what the floor is for.
#
# The nesting is declared here rather than fixed by deleting a test, because both numbers
# are worth reporting: what fired, and how much of it was independent.
IMPLIED_BY = {
    "is_family": ("purifying_selection",),
}


def reality_lines(family, thresholds):
    """Independent lines of evidence that this is a real protein.

    `thresholds` carries dnds_purifying_max, min_lineages and min_members, all from
    config/targets.yaml. Nothing here reads a module constant.

    Returns (n_independent, fired_independent, fired_implied).

    The names matter as much as the count: "3 of 4, missing purifying_selection" tells a
    reader what to go and check, which a score never can. `fired_implied` holds the lines
    that DID fire but add nothing, because another line that fired already entails them -
    see IMPLIED_BY. They are reported rather than dropped, so the table still shows
    everything that was true about the family.

    Absence of evidence never counts AGAINST a family here. A dN/dS status of
    NO_DIVERGENCE simply fails to fire this line; it is not scored as evidence of
    neutrality, because identical sequences are what a highly conserved family looks like.
    """
    fired = [t["name"] for t in REALITY_TESTS if t["test"](family, thresholds)]
    fired_set = set(fired)
    implied = [name for name in fired
               if any(other in fired_set for other in IMPLIED_BY.get(name, ()))]
    independent = [name for name in fired if name not in implied]
    return len(independent), independent, implied


def darkness_state(family):
    """Two named states, because the dark set is already dark by construction.

    DARK_NO_STRUCTURE   nothing named it and nothing recognises its fold - maximal novelty,
                        maximal risk, and the open-discovery arm of the portfolio.
    DARK_FOLD_KNOWN     no sequence homolog, but the fold is recognisable - a remote
                        homolog, so the function is guessable and the assay is choosable.

    These are different experiments, which is exactly why they must not be averaged into
    one number.
    """
    return "DARK_FOLD_KNOWN" if family.get("structural_match") else "DARK_NO_STRUCTURE"




def reality_thresholds(evolution, min_lineages=2):
    """Assemble the reality-test thresholds.

    min_members comes from evolution.min_members_for_dnds rather than being restated, so
    the nesting IMPLIED_BY relies on cannot drift apart.

    min_lineages is how many Stage 6 lineages make recurrence independent rather than
    clonal (spec section 2.6): two, the smallest count that is more than one.
    """
    return {"dnds_purifying_max": evolution["dnds_purifying_max"],
            "min_members": evolution["min_members_for_dnds"],
            "min_lineages": min_lineages}
