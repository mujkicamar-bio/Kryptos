"""Which proteins the cascade annotates (S2s, cascade_selection).

The study is about small plasmids. A protein is annotated when it is not already annotated
by PlasmidScope (Tier 0) and its family holds a small-plasmid protein that is not
annotated by Tier 0 either - that is, a family in which a small plasmid carries a protein
still to be explained. Large-plasmid members of those families are annotated with the same
cascade, so "known" and "unknown" mean the same on both sides of every such family.

The family (the primary clustering, intermediate: 50% identity, 80% coverage) decides WHICH
proteins; the search clustering (90% identity; UniRef90) decides which of them are actually
SEARCHED. No annotation is transferred at the family level: 50% identity says two proteins
are related, not that they do the same thing.
"""


def select(clusters, on_small, ps_annotated):
    """The proteins to annotate: {seq_id}.

    clusters      {representative: [member, ...]} at the selecting resolution; every
                  protein is a member of exactly one cluster
    on_small      seq_ids that occur on at least one small plasmid
    ps_annotated  seq_ids PlasmidScope annotates (Tier 0 ANNOTATED)
    """
    selected = set()
    for members in clusters.values():
        open_ = [m for m in members if m not in ps_annotated]
        if any(m in on_small for m in open_):
            selected.update(open_)
    return selected
