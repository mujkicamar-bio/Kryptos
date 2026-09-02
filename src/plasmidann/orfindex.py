def assign_orf_ids(orfs):
    """Number every ORF per plasmid in coordinate order as {plasmid_id}|{n}, n from 1.

    Ids are assigned once, over the complete ORF set for a plasmid. Re-indexing an
    already-indexed collection is refused: numbering a subset from 1 is how the same
    protein came to hold two different ids in 2026-08.
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
