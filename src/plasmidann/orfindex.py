"""ORF identifiers: {plasmid_id}|{n}, assigned once over every ORF of a plasmid."""


def assign_orf_ids(orfs):
    """Number every ORF per plasmid in coordinate order as {plasmid_id}|{n}, n from 1.

    Ids are assigned once, over the complete ORF set for a plasmid. Already-indexed input
    is refused, so a subset is never renumbered from 1.
    """
    if any("orf_id" in o for o in orfs):
        raise ValueError("orfs are already indexed; assign_orf_ids must not renumber them")
    ordered = sorted(orfs, key=lambda o: (o["plasmid_id"], o["start"], o["end"]))
    n = {}
    for o in ordered:
        pid = o["plasmid_id"]
        n[pid] = n.get(pid, 0) + 1
        o["orf_id"] = f"{pid}|{n[pid]}"
    return ordered
