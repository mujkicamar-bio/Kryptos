from plasmidann.cascade import classify

TIERS = ["T1", "T2", "T3", "T4", "T5", "T6", "T7", "T8"]


def test_hypothetical_protein_hit_is_not_an_annotation():
    """A homolog that is itself unnamed leaves the protein unknown - the prize class."""
    hits = [{"tier": "T6", "label": "hypothetical protein", "coverage": 0.95}]

    r = classify(hits, explained=0.0, min_coverage=0.5, tier_order=TIERS)

    assert r["functional_class"] == "UNCHARACTERIZED_HOMOLOG"
    assert r["annot_tier"] == "T6"


def test_partial_coverage_is_domain_only_not_functional():
    """Median Pfam alignment covered ~52% of these proteins; a third is not an annotation."""
    hits = [{"tier": "T1", "label": "HTH_3", "coverage": 0.33}]

    r = classify(hits, explained=0.33, min_coverage=0.5, tier_order=TIERS)

    assert r["functional_class"] == "DOMAIN_ONLY"


def test_no_hits_anywhere_is_none():
    r = classify([], explained=0.0, min_coverage=0.5, tier_order=TIERS)
    assert r["functional_class"] == "NONE"


def test_shallowest_functional_hit_wins_regardless_of_input_order():
    """Cascade order is authority order (design 6.3 rule 6)."""
    hits = [
        {"tier": "T6", "label": "conjugal transfer protein TraD", "coverage": 0.9},
        {"tier": "T1", "label": "Rep_1", "coverage": 0.85},
    ]

    r = classify(hits, explained=0.9, min_coverage=0.5, tier_order=TIERS)

    assert r["annot_tier"] == "T1"
    assert r["annot_label"] == "Rep_1"


def test_a_deep_real_name_beats_a_shallow_uninformative_one():
    """A 'hypothetical' hit at T4 must not block a real name found at T6."""
    hits = [
        {"tier": "T4", "label": "hypothetical protein", "coverage": 0.99},
        {"tier": "T6", "label": "relaxase MobA", "coverage": 0.8},
    ]

    r = classify(hits, explained=0.8, min_coverage=0.5, tier_order=TIERS)

    assert r["functional_class"] == "FUNCTIONAL"
    assert r["annot_tier"] == "T6"


def test_the_winning_hits_coverages_and_evalue_are_preserved():
    """Coverage decides the class AND must survive into the table, for stratification."""
    hits = [{"tier": "T4", "label": "relaxase MobA", "coverage": 0.91,
             "target_coverage": 0.42, "evalue": "1e-40"}]

    r = classify(hits, explained=0.91, min_coverage=0.5, tier_order=TIERS)

    assert r["annot_qcov"] == 0.91
    assert r["annot_tcov"] == 0.42
    assert r["annot_evalue"] == "1e-40"


def test_unresolved_proteins_carry_no_coverage():
    r = classify([], explained=0.0, min_coverage=0.5, tier_order=TIERS)
    assert r["annot_qcov"] is None and r["annot_tcov"] is None


# --- the class is a property of the protein, not of one alignment -------------------

def test_a_protein_explained_by_two_partial_domains_is_functional():
    """RepA_N and Bac_RepA_C cover 92% of a replication initiator between them. Neither
    single domain clears min_coverage, so per-hit classification calls this DOMAIN_ONLY
    and sends a fully annotated backbone protein to the screening pool. Measured on
    17% of all Pfam-hit proteins."""
    hits = [
        {"tier": "T1", "label": "RepA_N", "coverage": 0.40, "evalue": "1e-30"},
        {"tier": "T1", "label": "Bac_RepA_C", "coverage": 0.52, "evalue": "1e-25"},
    ]

    r = classify(hits, explained=0.92, min_coverage=0.5, tier_order=TIERS)

    assert r["functional_class"] == "FUNCTIONAL"


def test_one_small_domain_on_a_long_protein_stays_domain_only():
    """The genuine case the DOMAIN_ONLY class exists for: most of the protein is
    unexplained, so it remains a screening target."""
    hits = [{"tier": "T1", "label": "HTH_3", "coverage": 0.15, "evalue": "1e-12"}]

    r = classify(hits, explained=0.15, min_coverage=0.5, tier_order=TIERS)

    assert r["functional_class"] == "DOMAIN_ONLY"


def test_the_label_comes_from_the_strongest_evalue_not_the_widest_alignment():
    """Measured: 15.4% of labels change. ABC_membrane aligns further, Peptidase_C39 is
    seventeen orders of magnitude more significant and is the right name."""
    hits = [
        {"tier": "T1", "label": "ABC_membrane", "coverage": 0.62, "evalue": "1e-23"},
        {"tier": "T1", "label": "Peptidase_C39", "coverage": 0.28, "evalue": "6.5e-40"},
    ]

    r = classify(hits, explained=0.75, min_coverage=0.5, tier_order=TIERS)

    assert r["annot_label"] == "Peptidase_C39"
    assert r["annot_evalue"] == "6.5e-40"


def test_how_many_domains_built_the_explanation_is_reported():
    """0.9 explained by one domain and by six fragments are different claims."""
    hits = [
        {"tier": "T1", "label": "RepA_N", "coverage": 0.40, "evalue": "1e-30"},
        {"tier": "T1", "label": "Bac_RepA_C", "coverage": 0.52, "evalue": "1e-25"},
        {"tier": "T3", "label": "hypothetical protein", "coverage": 0.9, "evalue": "1e-8"},
    ]

    r = classify(hits, explained=0.92, min_coverage=0.5, tier_order=TIERS)

    assert r["n_informative_hits"] == 2


def test_homology_depth_follows_the_configured_tiers_not_a_fixed_list():
    """A four-tier cascade must not report depths from an eight-tier constant."""
    hits = [{"tier": "T3", "label": "relaxase MobA", "coverage": 0.9, "evalue": "1e-40"}]

    r = classify(hits, explained=0.9, min_coverage=0.5,
                 tier_order=["T1", "T2", "T3", "T4"])

    assert r["homology_depth"] == 3


# --- homology_depth means how deep we had to dig, not where the best E-value was -------

def test_homology_depth_is_the_shallowest_tier_that_named_the_protein():
    """`homology_depth` asks "how far did we have to dig before anything named it?", and
    that is a property of the CASCADE, not of one hit. A protein Pfam named at T1 is a
    shallow, well-characterised protein even if a deeper tier produced a better E-value
    for the same thing - and depth is used downstream to describe how obscure a protein
    is."""
    hits = [
        {"tier": "T1", "label": "RepA_N", "coverage": 0.4, "evalue": "1e-20",
         "start": 1, "end": 40},
        {"tier": "T4", "label": "replication initiator protein RepA", "coverage": 0.9,
         "evalue": "1e-80", "start": 1, "end": 90},
    ]
    out = classify(hits, explained=0.9, min_coverage=0.5,
                   tier_order=["T1", "T2", "T3", "T4"])

    assert out["homology_depth"] == 1, (
        "T1 named this protein, so the cascade dug one tier - not four")


def test_homology_depth_reports_the_uninformative_depth_when_nothing_named_it():
    """A protein only ever called 'hypothetical' still has a depth: how far we searched
    before anything at all recognised it. That is the number the dark set is described by."""
    hits = [
        {"tier": "T3", "label": "hypothetical protein", "coverage": 0.8,
         "evalue": "1e-30", "start": 1, "end": 80},
        {"tier": "T4", "label": "hypothetical protein", "coverage": 0.9,
         "evalue": "1e-60", "start": 1, "end": 90},
    ]
    out = classify(hits, explained=0.0, min_coverage=0.5,
                   tier_order=["T1", "T2", "T3", "T4"])
    assert out["functional_class"] == "UNCHARACTERIZED_HOMOLOG"
    assert out["homology_depth"] == 3


def test_a_family_assignment_without_coordinates_is_functional_not_domain_only():
    """The pharokka tier reports a family, an annotation and an E-value, and no span:
    pharokka deletes its alignment tables on exit. The families are whole-protein
    clusters, so the hit is a family-level assignment, not a domain. Classing it
    DOMAIN_ONLY because explained is 0 would say "a fragment of this protein matched",
    which is the opposite of what the tool reported."""
    hits = [{"tier": "T3", "label": "ParA-like partition protein", "coverage": "",
             "evalue": "1e-42", "start": "", "end": ""}]

    r = classify(hits, explained=0.0, min_coverage=0.5, tier_order=TIERS)

    assert r["functional_class"] == "FUNCTIONAL"
    assert r["annot_tier"] == "T3"
    assert r["span_measured"] == 0, (
        "the row must say its completeness was not measured, so that a reader does not "
        "take explained_fraction=0 for a fragment")


def test_a_measured_span_still_governs_when_one_is_present():
    """A Pfam domain covering a third of the protein AND a family assignment: FUNCTIONAL,
    because the family says the whole protein is known. But the span WAS measured, and
    the row says so."""
    hits = [{"tier": "T1", "label": "HTH_3", "coverage": 0.33, "evalue": "1e-10",
             "start": 1, "end": 40},
            {"tier": "T3", "label": "ParA-like partition protein", "coverage": "",
             "evalue": "1e-42", "start": "", "end": ""}]

    r = classify(hits, explained=0.33, min_coverage=0.5, tier_order=TIERS)

    assert r["functional_class"] == "FUNCTIONAL"
    assert r["span_measured"] == 1


# --- the label comes from the most authoritative tier (spec section 20) ----------------

def test_a_curated_tier_keeps_the_label_when_a_deeper_tier_aligns_better():
    """Automated transfer must not outrank curated evidence. Ranked on E-value alone, an
    nr title took the label from an informative Swiss-Prot hit on 39% of the proteins that
    had one. The strongest label is still reported beside it."""
    hits = [
        {"tier": "T4", "label": "Replication initiator protein RepA", "coverage": 0.9,
         "evalue": "1e-40", "start": 1, "end": 90},
        {"tier": "T5", "label": "WP_1.1 replication protein [Escherichia coli]",
         "coverage": 0.95, "evalue": "1e-90", "start": 1, "end": 95},
    ]
    out = classify(hits, explained=0.95, min_coverage=0.5,
                   tier_order=["T1", "T2", "T3", "T4", "T5"])

    assert out["annot_label"] == "Replication initiator protein RepA"
    assert out["annot_tier"] == "T4"
    assert out["best_evalue_label"] == "WP_1.1 replication protein [Escherichia coli]"
    assert out["best_evalue_tier"] == "T5"


def test_a_protein_named_only_by_its_domain_is_domain_only():
    """"X domain-containing protein" is PGAP's name for a domain-level assignment: it
    explains its span but does not make the protein FUNCTIONAL, however much it covers.
    Any full name beside it does."""
    domain = {"tier": "T5", "label": "WP_2.1 GNAT domain-containing protein [Bacillus]",
              "coverage": 0.95, "evalue": "1e-60", "start": 1, "end": 95}
    out = classify([domain], explained=0.95, min_coverage=0.5,
                   tier_order=["T1", "T2", "T3", "T4", "T5"])
    assert out["functional_class"] == "DOMAIN_ONLY"
    assert out["named_by_domain_only"] == 1

    full = {"tier": "T5", "label": "WP_3.1 N-acetyltransferase [Bacillus]",
            "coverage": 0.9, "evalue": "1e-50", "start": 1, "end": 90}
    out = classify([domain, full], explained=0.95, min_coverage=0.5,
                   tier_order=["T1", "T2", "T3", "T4", "T5"])
    assert out["functional_class"] == "FUNCTIONAL"
    assert out["named_by_domain_only"] == 0
