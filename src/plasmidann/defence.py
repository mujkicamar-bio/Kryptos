"""S8a: splitting DefenseFinder's two phases across the representation each one needs.

THE PROBLEM THIS SOLVES

DefenseFinder runs in two phases with incompatible requirements.

  Phase 1  HMM search, 1,887 profiles. Per protein. Order-independent.
           Answers: "does this protein look like a defence component?"

  Phase 2  MacSyFinder system calling, 711 model definitions. Each carries a quorum rule
           AND a co-localisation constraint. A real one reads:

               <model inter_gene_max_space="3" min_mandatory_genes_required="2" ...>

           Answers: "do these components sit together in a way that works?"

v2 ran both on unique_proteins.faa - dereplicated and ordered by SHA-256 hash - with the
default --db-type ordered_replicon, which asserts the file is in genomic order. Phase 1 was
unaffected. Phase 2 was reading a shuffled deck: proteins that were neighbours on a plasmid
were scattered, and proteins from unrelated plasmids became "adjacent".

Demonstrated in review: the same two proteins adjacent produced 1 system; separated by five
decoys, 0 systems - with identical HMM hits.

THE SPLIT

  1. SEARCH      1,887 profiles against the 3.5M dereplicated proteins. Cheap, no
                 redundancy, and order is irrelevant to what Phase 1 asks.
  2. PROPAGATE   component labels back to every ORF sharing that sequence. A protein
                 identical on forty plasmids is searched once and labelled forty times.
  3. PRUNE       keep only plasmids carrying at least one component. A plasmid with none
                 cannot produce a system, so it never needs the expensive ordered
                 representation. Measured effect: roughly a 4-5x reduction in what
                 reaches Phase 2.
  4. CALL        MacSyFinder in gembase mode over ordered per-plasmid gene lists, with
                 --replicon-topology circular, since 94% of these plasmids are closed.
"""


def propagate_components(component_hits, orf_to_seq):
    """Spread per-unique-protein component labels onto every ORF that shares the sequence.

    `component_hits` maps seq_id -> component name; `orf_to_seq` maps orf_id -> seq_id.
    Returns orf_id -> component name for the ORFs that inherit one.

    This is what makes dereplication safe for Phase 1: identity is exactly the property
    the HMM search depends on, so one search result is valid for every copy.
    """
    return {orf: component_hits[seq]
            for orf, seq in orf_to_seq.items() if seq in component_hits}


def candidate_plasmids(orf_labels):
    """Plasmids carrying at least one defence component, and therefore worth ordering.

    A plasmid with no component cannot satisfy any model's quorum, so writing it out in
    genomic order for Phase 2 would be pure waste. orf_id is `<plasmid_id>|<ordinal>`, so
    the plasmid is the part before the last pipe.
    """
    return {orf.rsplit("|", 1)[0] for orf in orf_labels}


def order_orfs(orfs):
    """Sort one plasmid's ORFs into genomic order.

    An origin-spanning gene - reconstructed at S1 across the cut point of a circular
    plasmid - runs start..length then 1..end, so its start is GREATER than its end. Sorting
    naively on start would place it last when on the circle it is adjacent to the first
    gene, and MacSyFinder would measure the wrong number of intervening genes across the
    junction. 94% of these plasmids are circular, so this is the common case.

    Such genes are placed first: on a circular replicon written from coordinate 1, the gene
    straddling the origin is the one preceding position 1.
    """
    def key(o):
        spans = str(o.get("spans_origin", "0")) == "1"
        return (0 if spans else 1, int(o["start"]))
    return sorted(orfs, key=key)


def gembase_id(plasmid_id, position):
    """Build a MacSyFinder gembase identifier: <replicon>_<zero-padded position>.

    Gembase mode lets one file hold many replicons and still treat each separately, which
    turns tens of thousands of per-plasmid runs into one. MacSyFinder recovers the replicon
    by splitting on the LAST underscore.

    Plasmid ids here contain underscores of their own - COMPASS_AB007909.1 - so they are
    rewritten with hyphens. Without that, every plasmid would parse as a different,
    malformed replicon and co-localisation would break in a new way.

    The position is zero-padded so that lexical order matches genomic order; unpadded,
    position 10 would sort before position 2.
    """
    return f"{plasmid_id.replace('_', '-')}_{position:05d}"
