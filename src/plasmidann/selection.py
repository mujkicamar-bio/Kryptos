"""Which proteins the cascade annotates (rule cascade_selection).

A protein is annotated when it is not skipped - already annotated by PlasmidScope (Tier 0),
or flagged by AntiFam as a probable non-protein - and its family (clustering.primary)
holds a small-plasmid protein that is not skipped either. Large-plasmid members of those
families are annotated with the same cascade, so "known" and "unknown" mean the same on
both sides of every such family.
"""


def select(clusters, on_small, skipped):
    """The proteins to annotate: {seq_id}.

    clusters  {representative: [member, ...]} at the selecting resolution; every protein
              is a member of exactly one cluster
    on_small  seq_ids that occur on at least one small plasmid
    skipped   seq_ids the cascade never searches: those PlasmidScope annotates (Tier 0
              ANNOTATED) and those AntiFam flags
    """
    selected = set()
    for members in clusters.values():
        open_ = [m for m in members if m not in skipped]
        if any(m in on_small for m in open_):
            selected.update(open_)
    return selected
