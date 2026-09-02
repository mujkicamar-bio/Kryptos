from plasmidann.cascade import explained_fraction


def test_overlapping_domains_are_merged_not_summed():
    """Real case: MFS_1 hits a 392 aa protein at 8-337 and 263-387.
    Best-hit coverage says 0.84; the union says 0.97."""
    assert explained_fraction(392, [(8, 337), (263, 387)]) == 0.9694


def test_disjoint_domains_add_up():
    assert explained_fraction(100, [(1, 10), (51, 60)]) == 0.2


def test_a_protein_with_no_hits_is_entirely_unexplained():
    assert explained_fraction(300, []) == 0.0
