"""Stage 6: plasmid lineage clusters, an independence measure separate from MOB class.

WHY THIS EXISTS AT ALL

Spec section 33 opens with the point: "MOB classification and plasmid sequence similarity
are separate concepts." MOB typing describes the relaxase a plasmid carries. It says
nothing about whether two plasmids are the same molecule sequenced twice.

Section 33.2 gives the worked example, and the three numbers are the whole argument:

    plasmid_occurrences          2,143
    MOB_clusters                     8
    independent_plasmid_clusters    47

A protein on 2,143 plasmid records is not 2,143 observations. It may be one clinical
plasmid deposited 2,143 times. MOB breadth of 8 does not fix that, because a single
widespread lineage can carry one relaxase type. The 47 is the number that answers "how many
independent times has evolution done this", and section 2.6 makes distinguishing these a
design principle: "recurrence is not independence".

THE METHOD

Mash, which section 33.1 names ("a configured sequence-based approach such as Mash/sketch
or an appropriate ANI-like clustering method"). Each plasmid is sketched to a small set of
hashed k-mers, all pairs within a distance threshold are linked, and the connected
components are the lineage clusters.

SINGLE LINKAGE, AND WHAT IT COSTS

Connected components mean single linkage: A and C land in one cluster if both are close to
B, even when A and C are far apart. That is the wrong choice for taxonomy and the right one
here, because the question is "could these have descended from one recent ancestor", and a
chain of near-identical intermediates is evidence that they could. The cost is that single
linkage can chain distant things together through a dense middle, which INFLATES cluster
size and therefore UNDER-counts independence. That direction is the conservative one: it
makes recurrence look less independent than it is, never more.

THE THRESHOLD IS A PARAMETER, NOT A TRUTH

Mash distance 0.05 is roughly 95% ANI, the conventional species-level boundary and the
usual bar for "same plasmid lineage". It is declared in config and recorded on every row,
because the cluster count is a direct function of it and any statement about independence
inherits that choice.
"""
import collections


def parse_mash_dist(text, max_distance, max_pvalue=1e-10):
    """Edges from `mash dist` output, as (a, b) pairs below the distance threshold.

    Mash writes: reference<TAB>query<TAB>distance<TAB>p-value<TAB>shared-hashes

    Self-comparisons are dropped; every pair otherwise appears twice and the caller's
    union-find is idempotent, so no deduplication is needed here.

    The p-value filter is separate from the distance filter and both are needed. A small
    distance computed from two shared hashes out of a thousand is not evidence of anything,
    and Mash reports exactly that case with a large p-value. Filtering on distance alone
    would link plasmids that share almost nothing.
    """
    for line in text.splitlines():
        if not line.strip():
            continue
        fields = line.split("\t")
        if len(fields) < 4:
            continue
        a, b, distance, pvalue = fields[0], fields[1], fields[2], fields[3]
        if a == b:
            continue
        try:
            if float(distance) <= max_distance and float(pvalue) <= max_pvalue:
                yield a, b
        except ValueError:
            continue


def connected_components(names, edges):
    """Group `names` into connected components given an iterable of (a, b) edges.

    Union-find with path compression. Every name gets a component even with no edges, so a
    plasmid resembling nothing is a lineage cluster of one rather than absent from the
    table - which is the honest record, and is also the common case for a novel plasmid.
    """
    parent = {name: name for name in names}

    def find(x):
        root = x
        while parent[root] != root:
            root = parent[root]
        while parent[x] != root:
            parent[x], x = root, parent[x]
        return root

    for a, b in edges:
        # An edge naming a plasmid that is not in `names` means the sketch set and the
        # analysis set disagree; ignoring it silently would drop a real link, so the
        # plasmid is admitted and the disagreement shows up as a cluster.
        for node in (a, b):
            if node not in parent:
                parent[node] = node
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb

    groups = collections.defaultdict(list)
    for name in parent:
        groups[find(name)].append(name)
    return groups


def lineage_clusters(names, edges):
    """Assign each plasmid a lineage cluster id.

    The id is the lexicographically smallest member of its component, not an ordinal.
    Cluster numbering that depended on iteration order would change between runs on
    identical input, and every count derived from it would be uncomparable - the same
    defect the family ids had (spec section 5.4).
    """
    groups = connected_components(names, edges)
    assignment = {}
    for members in groups.values():
        cluster_id = min(members)
        for member in members:
            assignment[member] = cluster_id
    return assignment
