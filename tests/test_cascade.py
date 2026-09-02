from plasmidann.cascade import classify


def test_hypothetical_protein_hit_is_not_an_annotation():
    """A homolog that is itself unnamed leaves the protein unknown - the prize class."""
    hits = [{"tier": "T6", "label": "hypothetical protein", "coverage": 0.95}]

    r = classify(hits)

    assert r["functional_class"] == "UNCHARACTERIZED_HOMOLOG"
    assert r["annot_tier"] == "T6"


def test_partial_coverage_is_domain_only_not_functional():
    """Median Pfam alignment covered ~52% of these proteins; a third is not an annotation."""
    hits = [{"tier": "T1", "label": "HTH_3", "coverage": 0.33}]

    assert classify(hits)["functional_class"] == "DOMAIN_ONLY"


def test_no_hits_anywhere_is_none():
    assert classify([])["functional_class"] == "NONE"


def test_shallowest_functional_hit_wins_regardless_of_input_order():
    """Cascade order is authority order (design 6.3 rule 6)."""
    hits = [
        {"tier": "T6", "label": "conjugal transfer protein TraD", "coverage": 0.9},
        {"tier": "T1", "label": "Rep_1", "coverage": 0.85},
    ]

    r = classify(hits)

    assert r["annot_tier"] == "T1"
    assert r["annot_label"] == "Rep_1"


def test_a_deep_real_name_beats_a_shallow_uninformative_one():
    """A 'hypothetical' hit at T4 must not block a real name found at T6."""
    hits = [
        {"tier": "T4", "label": "hypothetical protein", "coverage": 0.99},
        {"tier": "T6", "label": "relaxase MobA", "coverage": 0.8},
    ]

    r = classify(hits)

    assert r["functional_class"] == "FUNCTIONAL"
    assert r["annot_tier"] == "T6"


def test_the_winning_hits_coverages_and_evalue_are_preserved():
    """Coverage decides the class AND must survive into the table, for stratification."""
    hits = [{"tier": "T4", "label": "relaxase MobA", "coverage": 0.91,
             "target_coverage": 0.42, "evalue": "1e-40"}]

    r = classify(hits)

    assert r["annot_qcov"] == 0.91
    assert r["annot_tcov"] == 0.42
    assert r["annot_evalue"] == "1e-40"


def test_unresolved_proteins_carry_no_coverage():
    r = classify([])
    assert r["annot_qcov"] is None and r["annot_tcov"] is None
