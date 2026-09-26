"""Plasmid lineages: the unit of independence, separate from MOB class.

Mash (Ondov et al. 2016, Genome Biol 17:132) distances between whole plasmids. Pairs with
distance <= max_distance and p-value <= max_pvalue are linked, and the connected components
(single linkage) are the lineages. Single linkage can chain distant plasmids through
near-identical intermediates, which enlarges lineages and so undercounts independence. The
thresholds are configuration values (config/targets.yaml) recorded on every output row.
The lineage id is its smallest member id, so it does not depend on input order.
"""
import collections


def parse_mash_dist(lines, max_distance, max_pvalue=1e-10):
    """(a, b) edges from `mash dist` output lines, self-comparisons dropped.

    Mash writes: reference<TAB>query<TAB>distance<TAB>p-value<TAB>shared-hashes. The
    p-value filter rejects a small distance computed from very few shared hashes.
    """
    for line in lines:
        a, b, distance, pvalue = line.split("\t")[:4]
        if a != b and float(distance) <= max_distance and float(pvalue) <= max_pvalue:
            yield a, b


def connected_components(names, edges):
    """{root: [members]} for `names` given (a, b) edges; a name with no edge is alone."""
    parent = {name: name for name in names}

    def find(x):
        root = x
        while parent[root] != root:
            root = parent[root]
        while parent[x] != root:
            parent[x], x = root, parent[x]
        return root

    for a, b in edges:
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb

    groups = collections.defaultdict(list)
    for name in parent:
        groups[find(name)].append(name)
    return groups


def lineage_clusters(names, edges):
    """{plasmid: lineage id}, the id being the smallest member of the component."""
    assignment = {}
    for members in connected_components(names, edges).values():
        cluster_id = min(members)
        for member in members:
            assignment[member] = cluster_id
    return assignment
