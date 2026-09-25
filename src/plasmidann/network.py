"""The family network: 50%-identity clusters linked by sequence similarity.

Built as Durairaj et al. built theirs over UniRef50 (Nature 2023, 622:646, Methods), so
the two sets of numbers can be compared directly:

  nodes   cluster representatives (ours: the `intermediate` clustering, 50% identity,
          80% coverage - the UniRef50 analogue)
  edges   MMseqs2 all-against-all, "an alignment that covers at least 50% of one of the
          proteins with an E value < 1 x 10^-4", at most four outbound edges per node
  dark    "functional brightness" - full-length annotation coverage - at most 5%
  groups  asynchronous label propagation (networkx asyn_lpa_communities)

A family network is a map for reading, not a classification: nothing here changes which
proteins are dark or which family a protein belongs to.
"""
import math

import networkx as nx


def edges_from_hits(hits, max_out, min_cov, max_evalue):
    """Undirected weighted edges from MMseqs2 hit rows.

    `hits` are dicts with query, target, evalue, qcov, tcov (coverages as fractions). A hit
    counts when it covers at least `min_cov` of EITHER protein - a small protein wholly
    inside a larger one is a relationship - and its E-value is below `max_evalue`. Each
    query contributes its `max_out` best targets by E-value. An edge found in both
    directions is kept once, with the better E-value.
    """
    by_query = {}
    for h in hits:
        if h["query"] == h["target"]:
            continue
        ev = float(h["evalue"])
        if ev >= max_evalue or max(float(h["qcov"]), float(h["tcov"])) < min_cov:
            continue
        by_query.setdefault(h["query"], []).append((ev, h["target"]))

    edges = {}
    for q, targets in by_query.items():
        for ev, t in sorted(targets)[:max_out]:
            key = tuple(sorted((q, t)))
            if key not in edges or ev < edges[key]:
                edges[key] = ev
    return edges


def member_brightness(row):
    """Annotation coverage of one protein, 0..1, from its protein_annotation row.

    A FUNCTIONAL protein whose span was never measured - a PlasmidScope Tier-0 transfer or
    a family-level pharokka assignment - is fully known, so it counts as 1.0; reading its
    empty or zero explained_fraction as 'dark' would say the opposite of the evidence.
    """
    if row.get("functional_class") == "FUNCTIONAL" and \
            row.get("annot_completeness") == "NOT_MEASURED":
        return 1.0
    try:
        return float(row.get("explained_fraction") or 0.0)
    except ValueError:
        return 0.0


def node_brightness(member_rows):
    """The brightness a cluster REACHES: the best-annotated member's coverage."""
    return max((member_brightness(r) for r in member_rows), default=0.0)


def communities(edges, nodes, seed):
    """{node: community index}. Asynchronous label propagation, seeded so it repeats.

    Label propagation is stochastic; the seed makes the partition reproducible, which is
    what lets a community id be quoted. Isolated nodes are communities of one.
    """
    g = nx.Graph()
    g.add_nodes_from(nodes)
    for (a, b), ev in edges.items():
        # Durairaj weight edges by E-value; a smaller E-value is a stronger link, so the
        # weight is its negative log, floored for E = 0.
        g.add_edge(a, b, weight=-math.log10(ev) if ev > 0 else 300.0)
    out = {}
    for i, comm in enumerate(sorted(nx.community.asyn_lpa_communities(
            g, weight="weight", seed=seed), key=lambda c: (-len(c), min(c)))):
        for n in comm:
            out[n] = i
    return out
