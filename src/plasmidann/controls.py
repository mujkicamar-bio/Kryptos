"""The positive control for the annotation cascade (success criterion SC2).

Following ECLIPSE (Bioinformatics 42:8), which recovered 99.2-100% of 246 virulence, 42
AMR and 75 essential genes as annotated. Anything KNOWN that emerges dark is a recall
failure in the cascade, and the run should stop rather than produce a target list built on
a broken annotation step. Its absence from v1 was the single point of unanimous reviewer
criticism.

The control set is independent by construction: reviewed Swiss-Prot proteins of known
function, spiked into the query set at S2c and carried through every tier exactly as a real
protein is. An earlier version drew the control from proteins the cascade itself had
labelled, which cannot detect the failure that matters, because a protein the cascade
MISSED never enters a self-drawn control set.

This lived in plasmidann.backbone alongside a curated list of plasmid-backbone family
names. It was never part of that list and does not depend on it, so it moved here when the
list was deleted.
"""


def control_recall(rows):
    """Fraction of a known-function control set that the cascade called FUNCTIONAL.

    The gate for success criterion SC2. Raises on an empty control set rather than
    returning 1.0 or 0.0: an empty control silently passing is exactly how a gate stops
    being a gate, and it is the failure mode that would reproduce v1's missing control.
    """
    if not rows:
        raise ValueError(
            "the positive control set is empty - a gate that cannot fail is not a gate. "
            "Check that the control ids resolved against the annotation table.")
    n = sum(1 for r in rows if r.get("functional_class") == "FUNCTIONAL")
    return round(n / len(rows), 4)
