from plasmidann.cascade import informative_spans, explained_fraction


def test_a_hypothetical_hit_explains_nothing_however_well_it_aligns():
    """A 95%-coverage hit to 'hypothetical protein' must not satisfy the narrowing
    threshold - nothing has been explained, so the search must continue."""
    hits = [{"label": "MULTISPECIES: hypothetical protein", "start": 5, "end": 195}]

    assert informative_spans(hits) == []
    assert explained_fraction(200, informative_spans(hits)) == 0.0


def test_only_real_function_names_contribute_to_explained_fraction():
    hits = [
        {"label": "hypothetical protein", "start": 1, "end": 100},
        {"label": "relaxase MobA", "start": 101, "end": 200},
    ]

    assert informative_spans(hits) == [(101, 200)]
    assert explained_fraction(200, informative_spans(hits)) == 0.5


def test_uninformative_hits_are_still_reported_not_discarded():
    """Someone else has seen this protein. That is worth recording."""
    from plasmidann.cascade import uninformative_hits

    hits = [
        {"label": "hypothetical protein", "tier": "T5", "start": 1, "end": 50},
        {"label": "relaxase MobA", "tier": "T1", "start": 60, "end": 90},
    ]

    kept = uninformative_hits(hits)

    assert len(kept) == 1
    assert kept[0]["label"] == "hypothetical protein"
    assert kept[0]["tier"] == "T5"
