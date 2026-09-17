"""S8a: splitting DefenseFinder's two phases across the representation each one needs.

DefenseFinder runs in two phases. Phase 1 is an HMM search over 1,887 profiles: per
protein, order-independent. Phase 2 is MacSyFinder system calling over 711 model
definitions, each carrying a co-localisation constraint - a real one reads

    <model inter_gene_max_space="3" min_mandatory_genes_required="2" ...>

so it is entirely about gene adjacency.

v2 ran the whole thing on unique_proteins.faa: dereplicated and ordered by SHA-256 hash,
with the default --db-type ordered_replicon. Phase 1 was fine. Phase 2 was told that a
cryptographic hash ordering was genomic order. Demonstrated in review: the same two
proteins adjacent gave 1 system, separated by 5 decoys gave 0, with identical HMM hits.

The fix searches the dereplicated set once, propagates component labels to every ORF that
shares the sequence, and calls systems on ordered per-plasmid gene lists.
"""
from plasmidann.defence import (propagate_components, candidate_plasmids,
                                order_orfs, gembase_id)


ORF_TO_SEQ = {
    "p1|1": "aaa", "p1|2": "bbb", "p1|3": "ccc",
    "p2|1": "aaa", "p2|2": "zzz",
    "p3|1": "yyy",
}


def test_a_component_hit_propagates_to_every_orf_sharing_the_sequence():
    """This is the whole point of dereplicating first: a protein identical on forty
    plasmids is searched once and labelled forty times."""
    hits = {"aaa": "RM_type_II"}

    labels = propagate_components(hits, ORF_TO_SEQ)

    assert labels["p1|1"] == "RM_type_II"
    assert labels["p2|1"] == "RM_type_II"
    assert "p1|2" not in labels


def test_propagation_does_not_invent_labels_for_unhit_sequences():
    assert propagate_components({}, ORF_TO_SEQ) == {}


# --- pruning: only plasmids that could possibly carry a system ------------------------

def test_only_plasmids_carrying_a_component_become_candidates():
    """System calling needs ordered proteins, which is the expensive representation. A
    plasmid with no component hit cannot produce a system, so it never needs to be
    written out in genomic order."""
    labels = {"p1|1": "RM_type_II", "p2|1": "RM_type_II"}

    assert candidate_plasmids(labels) == {"p1", "p2"}


def test_a_plasmid_with_no_components_is_pruned():
    assert "p3" not in candidate_plasmids({"p1|1": "RM_type_II"})


# --- ordering: the thing the hash ordering destroyed ---------------------------------

def test_orfs_are_returned_in_genomic_order():
    """MacSyFinder counts intervening genes. If the order is wrong, inter_gene_max_space
    is measuring nothing."""
    orfs = [
        {"orf_id": "p1|3", "start": 900, "end": 1200, "spans_origin": "0"},
        {"orf_id": "p1|1", "start": 10, "end": 300, "spans_origin": "0"},
        {"orf_id": "p1|2", "start": 400, "end": 800, "spans_origin": "0"},
    ]

    assert [o["orf_id"] for o in order_orfs(orfs)] == ["p1|1", "p1|2", "p1|3"]


def test_an_origin_spanning_gene_sorts_to_the_start_of_the_replicon():
    """A gene reconstructed across the cut runs start..length then 1..end, so start > end.
    Sorting naively on start would place it last, when on the circle it is adjacent to the
    first gene. 94% of these plasmids are circular, so this is the common case, not an
    edge case."""
    orfs = [
        {"orf_id": "p1|2", "start": 400, "end": 800, "spans_origin": "0"},
        {"orf_id": "p1|1", "start": 4800, "end": 200, "spans_origin": "1"},
        {"orf_id": "p1|3", "start": 900, "end": 1200, "spans_origin": "0"},
    ]

    assert [o["orf_id"] for o in order_orfs(orfs)] == ["p1|1", "p1|2", "p1|3"]


# --- gembase naming ------------------------------------------------------------------

def test_gembase_ids_encode_replicon_and_position():
    """MacSyFinder's gembase mode splits an id on the last underscore to recover the
    replicon, so one run can hold many plasmids and still treat each separately."""
    assert gembase_id("COMPASS_AB007909.1", 7) == "COMPASS-AB007909.1_00007"


def test_gembase_ids_are_sortable_in_genomic_order():
    """Zero padding matters: position 2 must sort before position 10."""
    ids = [gembase_id("p1", i) for i in (10, 2, 1)]
    assert sorted(ids) == [gembase_id("p1", 1), gembase_id("p1", 2), gembase_id("p1", 10)]


def test_underscores_in_a_plasmid_id_do_not_break_the_replicon_split():
    """Plasmid ids here contain underscores (COMPASS_AB007909.1). MacSyFinder splits on
    the LAST underscore, so an unescaped id would make every plasmid its own malformed
    replicon name."""
    gid = gembase_id("COMPASS_AB007909.1", 1)
    replicon, position = gid.rsplit("_", 1)
    assert replicon == "COMPASS-AB007909.1"
    assert position == "00001"
