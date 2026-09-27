"""Evidence description per dark family: the reality lines and the darkness state.

This module describes evidence and does not select on it: no score, rank or shortlist.

  reality_lines     Independent lines of evidence that a family is a real protein rather
                    than a gene-calling artefact, counted AND named, so that a reader sees
                    which line is missing. Absence of evidence never counts against a
                    family: a dN/dS status of NO_DIVERGENCE fails to fire a line and is
                    not read as neutral evolution.

  darkness_state    Whether the fold of a dark family is recognisable, as a named state.
"""

# Present on this many Stage 6 lineages makes recurrence independent rather than clonal:
# two, the smallest count that is more than one.
MIN_LINEAGES = 2

# The four lines of evidence that a dark ORF is a real protein, each a named boolean over a
# family row and the thresholds from reality_thresholds.
REALITY_TESTS = [
    {
        "name": "purifying_selection",
        # Median dN/dS below evolution.dnds_purifying_max.
        "test": lambda f, thr: (f.get("dnds_status") == "MEASURED"
                                and _num(f.get("dnds_median")) is not None
                                and _num(f.get("dnds_median")) < thr["dnds_purifying_max"]),
    },
    {
        "name": "multi_lineage",
        # Present on two or more Stage 6 lineages (plasmid_lineage, Mash): one clone
        # sequenced forty times is one lineage. Not MOB-suite clusters, which assign the
        # nearest reference however distant and so put unrelated novel plasmids together.
        "test": lambda f, thr: (_int(f.get("independent_plasmid_cluster_count"))
                                >= thr["min_lineages"]),
    },
    {
        "name": "is_family",
        # FESNov's minimum family size, and enough to align: the same number as
        # evolution.min_members_for_dnds (see reality_thresholds).
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


# Lines implied by another line. purifying_selection needs a MEASURED dN/dS, which is
# computed only for families of at least min_members_for_dnds members - the is_family test -
# so is_family adds nothing when purifying_selection fired. Both are reported: what fired,
# and how much of it was independent.
IMPLIED_BY = {
    "is_family": ("purifying_selection",),
}


def reality_lines(family, thresholds):
    """Independent lines of evidence that this is a real protein.

    `thresholds` comes from reality_thresholds. Returns (n_independent, fired_independent,
    fired_implied): `fired_implied` holds the lines that fired but are entailed by another
    line that fired (IMPLIED_BY), reported rather than dropped.
    """
    fired = [t["name"] for t in REALITY_TESTS if t["test"](family, thresholds)]
    fired_set = set(fired)
    implied = [name for name in fired
               if any(other in fired_set for other in IMPLIED_BY.get(name, ()))]
    independent = [name for name in fired if name not in implied]
    return len(independent), independent, implied


def darkness_state(family):
    """DARK_FOLD_KNOWN when the family has a structural match, else DARK_NO_STRUCTURE.

    DARK_NO_STRUCTURE   nothing names it and nothing recognises its fold.
    DARK_FOLD_KNOWN     no sequence homolog, but a recognisable fold: a remote homolog, so a
                        function can be proposed.
    """
    return "DARK_FOLD_KNOWN" if family.get("structural_match") else "DARK_NO_STRUCTURE"


def reality_thresholds(evolution):
    """The reality-test thresholds from the `evolution` config section.

    min_members is evolution.min_members_for_dnds, so the nesting IMPLIED_BY relies on
    cannot drift apart; min_lineages is MIN_LINEAGES.
    """
    return {"dnds_purifying_max": evolution["dnds_purifying_max"],
            "min_members": evolution["min_members_for_dnds"],
            "min_lineages": MIN_LINEAGES}
