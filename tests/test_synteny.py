"""Stage 9: synteny and context conservation.

Stage 8 asks what a dark ORF sits next to, once. This asks whether the arrangement RECURS:

    A - B - DARK - C - D
    A - B - DARK - C - D
    A - B - DARK - C
    A - B - DARK - C - D

Not four observations of "near B" but one conserved gene order seen four times, which is a
far stronger claim: order survives because the arrangement matters.

The counting unit is the plasmid LINEAGE (Stage 6), not the occurrence: forty copies of one
redeposited plasmid are one observation, not forty.
"""
import pytest
from conftest import FakeSnakemake, read_tsv, run_script, write_tsv

from plasmidann import synteny


def one_lineage_each(occurrences):
    """Give every occurrence its own lineage: the case where lineages and copies coincide."""
    return [dict(o, lineage=f"L{i}") for i, o in enumerate(occurrences)]


def conserved_example():
    """The four occurrences of the module docstring, each on its own lineage."""
    return one_lineage_each([
        {"left": ["B", "A"], "right": ["C", "D"], "operon": True},
        {"left": ["B", "A"], "right": ["C", "D"], "operon": True},
        {"left": ["B", "A"], "right": ["C"], "operon": True},
        {"left": ["B", "A"], "right": ["C", "D"], "operon": True},
    ])


def test_the_example_is_fully_conserved():
    result = synteny.conservation(conserved_example())

    assert result["lineage_left_conservation"] == 1.0
    assert result["lineage_right_conservation"] == 1.0
    assert result["lineage_synteny_conservation"] == 1.0
    assert result["modal_synteny"] == "B|C"
    assert result["status"] == "SUCCESS"


def test_a_conserved_left_and_a_variable_right_are_reported_separately():
    """A real and common arrangement - the left gene may be the promoter-sharing partner -
    that a single averaged context score would hide."""
    occurrences = one_lineage_each([
        {"left": ["B"], "right": ["C"]},
        {"left": ["B"], "right": ["X"]},
        {"left": ["B"], "right": ["Y"]},
        {"left": ["B"], "right": ["Z"]},
    ])

    result = synteny.conservation(occurrences)

    assert result["lineage_left_conservation"] == 1.0
    assert result["lineage_right_conservation"] == 0.25
    assert result["lineage_synteny_conservation"] == 0.25, (
        "ordered synteny cannot exceed its weakest side")


def test_ordered_synteny_is_stricter_than_the_unordered_neighbourhood():
    """The same two genes in swapped positions are one neighbourhood but not one synteny."""
    occurrences = one_lineage_each([
        {"left": ["B"], "right": ["C"]},
        {"left": ["C"], "right": ["B"]},
    ])

    result = synteny.conservation(occurrences)

    assert result["lineage_neighborhood_conservation"] == 1.0
    assert result["lineage_synteny_conservation"] == 0.5


def test_a_missing_neighbour_is_not_counted_as_a_failure_to_conserve():
    """An ORF at the end of a contig has no left neighbour. Counting that against the
    family would penalise it for where the assembler cut, not for its biology."""
    occurrences = one_lineage_each([
        {"left": [], "right": ["C"]},
        {"left": ["B"], "right": ["C"]},
        {"left": ["B"], "right": ["C"]},
    ])

    result = synteny.conservation(occurrences)

    assert result["lineage_left_conservation"] == 1.0, (
        "the truncated occurrence was counted as a conservation failure")
    assert result["lineage_right_conservation"] == 1.0


def test_a_single_lineage_is_not_measured():
    """Conservation over one lineage is 1.0 by construction and means nothing. Emitting it
    would put every singleton at the top of a ranking of conserved context."""
    result = synteny.conservation(one_lineage_each([{"left": ["B"], "right": ["C"]}]))

    assert result["status"] == "TOO_FEW_LINEAGES"
    assert result["lineage_synteny_conservation"] == ""


def test_a_family_with_no_context_at_all_is_reported_as_such():
    """Distinct from 'measured and found variable', and from 'too few lineages'."""
    result = synteny.conservation(one_lineage_each([{"left": [], "right": []},
                                                    {"left": [], "right": []}]))

    assert result["context_recurrence"] == 0
    assert result["n_lineages"] == 0
    assert result["status"] == "NO_CONTEXT"
    assert result["lineage_synteny_conservation"] == ""


def test_the_occurrence_and_lineage_counts_travel_with_the_measurement():
    result = synteny.conservation(conserved_example())

    assert result["n_occurrences"] == 4
    assert result["context_recurrence"] == 4
    assert result["n_lineages"] == 4
    assert result["n_lineages_discordant"] == 0


def test_operon_conservation_is_its_own_measurement():
    occurrences = one_lineage_each([
        {"left": ["B"], "right": ["C"], "operon": True},
        {"left": ["B"], "right": ["C"], "operon": False},
    ])

    result = synteny.conservation(occurrences)

    assert result["lineage_operon_like_conservation"] == 0.5
    assert result["lineage_synteny_conservation"] == 1.0, (
        "operon membership and gene order are different measurements")


# ---------------------------------------------------------------------------------------
# The lineage rule (fractional vote)
# ---------------------------------------------------------------------------------------

def test_forty_clonal_copies_are_one_vote():
    """The redeposition case: forty copies of one plasmid are one biological event, and
    one lineage is conserved with itself by construction."""
    occurrences = [{"left": ["B"], "right": ["C"], "operon": True, "lineage": "L1"}
                   for _ in range(40)]

    result = synteny.conservation(occurrences)

    assert result["n_occurrences"] == 40
    assert result["n_lineages"] == 1
    assert result["status"] == "TOO_FEW_LINEAGES"
    assert result["lineage_synteny_conservation"] == ""
    assert result["modal_synteny"] == ""


def test_one_copy_per_lineage_gives_the_share_among_occurrences():
    """With one copy per lineage, a conservation is the share of the most common non-empty
    value among the occurrences."""
    occurrences = one_lineage_each([
        {"left": ["B"], "right": ["C"], "operon": True},
        {"left": ["B"], "right": ["C"], "operon": False},
        {"left": ["B"], "right": ["X", "C"], "operon": True},
        {"left": ["Y"], "right": [], "operon": False},
        {"left": [], "right": ["C", "D"], "operon": True},
        {"left": ["B", "A"], "right": ["C"], "operon": True},
        {"left": [], "right": [], "operon": True},
    ])

    result = synteny.conservation(occurrences)

    # Six occurrences have context. Left: B in 4 of the 5 non-empty; right: C in 4 of 5;
    # neighbourhood {B,C} in 2 of 6; synteny B|C in 3 of 6; operon-like in 4 of 6.
    assert result["lineage_left_conservation"] == 0.8
    assert result["lineage_right_conservation"] == 0.8
    assert result["lineage_neighborhood_conservation"] == 0.3333
    assert result["lineage_synteny_conservation"] == 0.5
    assert result["lineage_operon_like_conservation"] == 0.6667
    assert result["n_lineages"] == 6


def test_a_side_with_a_single_voting_lineage_has_no_value():
    """Two lineages have context, but only L2 has a left neighbour: one lineage conserved
    with itself is not a measurement, so left conservation is empty while right, voted by
    both, is measured."""
    occurrences = [{"left": [], "right": ["R1"], "lineage": "L1"},
                   {"left": ["X"], "right": ["R1"], "lineage": "L2"}]

    result = synteny.conservation(occurrences)

    assert result["status"] == "SUCCESS"
    assert (result["lineage_left_conservation"], result["modal_left"]) == ("", "")
    assert (result["lineage_right_conservation"], result["modal_right"]) == (1.0, "R1")


def test_a_discordant_lineage_splits_its_vote():
    """L1 = {A|B, C|D}, L2 = {A|B}: L1 gives 1/2 to each of its arrangements, so A|B scores
    1.5 of 2 lineages. The within-lineage rearrangement is counted, not averaged away."""
    occurrences = [
        {"left": ["A"], "right": ["B"], "lineage": "L1"},
        {"left": ["C"], "right": ["D"], "lineage": "L1"},
        {"left": ["A"], "right": ["B"], "lineage": "L2"},
    ]

    result = synteny.conservation(occurrences)

    assert result["lineage_synteny_conservation"] == 0.75
    assert result["modal_synteny"] == "A|B"
    assert result["n_lineages"] == 2
    assert result["n_lineages_discordant"] == 1


def test_many_copies_in_one_lineage_do_not_outvote_two_lineages():
    """L1 = ten copies of X|Y, L2 = {P|Q}, L3 = {P|Q}. Over occurrences X|Y would win with
    10/12 = 0.8333; over lineages P|Q wins with 2/3."""
    occurrences = ([{"left": ["X"], "right": ["Y"], "lineage": "L1"}] * 10
                   + [{"left": ["P"], "right": ["Q"], "lineage": "L2"},
                      {"left": ["P"], "right": ["Q"], "lineage": "L3"}])

    result = synteny.conservation(occurrences)

    assert result["lineage_synteny_conservation"] == 0.6667
    assert result["modal_synteny"] == "P|Q"
    assert result["n_lineages_discordant"] == 0


def test_a_lineage_with_only_empty_values_does_not_vote_on_that_side():
    """L2 has no left neighbour anywhere (its copies sit at a record end): it does not vote
    for left, but its right neighbour still votes."""
    occurrences = [
        {"left": ["B"], "right": ["C"], "lineage": "L1"},
        {"left": [], "right": ["C"], "lineage": "L2"},
        {"left": [], "right": ["X"], "lineage": "L2"},
        {"left": ["B"], "right": ["X"], "lineage": "L3"},
    ]

    result = synteny.conservation(occurrences)

    assert result["lineage_left_conservation"] == 1.0, "L2 voted with an empty value"
    # right: L1 -> C (1), L2 -> C (1/2) + X (1/2), L3 -> X (1): a tie at 1.5 of 3.
    assert result["lineage_right_conservation"] == 0.5
    assert result["modal_right"] == "C"
    assert result["n_lineages"] == 3


def test_a_tie_is_broken_by_the_smallest_value_whatever_the_input_order():
    occurrences = [
        {"left": ["Z"], "right": ["C"], "lineage": "L1"},
        {"left": ["A"], "right": ["C"], "lineage": "L2"},
        {"left": ["M"], "right": ["C"], "lineage": "L3"},
        {"left": ["Z"], "right": ["C"], "lineage": "L4"},
        {"left": ["A"], "right": ["C"], "lineage": "L5"},
    ]

    forward = synteny.conservation(occurrences)
    backward = synteny.conservation(list(reversed(occurrences)))

    assert forward["modal_left"] == backward["modal_left"] == "A"
    assert forward["modal_synteny"] == backward["modal_synteny"] == "A|C"
    assert forward["lineage_left_conservation"] == 0.4


def test_fractional_ties_are_exact():
    """Ten copies worth 1/10 each make exactly one vote. Summed in floating point they make
    0.9999999999999999, which lost the tie to a single-copy lineage and broke the
    lexicographic rule."""
    occurrences = ([{"left": ["A"], "right": ["C"], "lineage": "L1"}] * 10
                   + [{"left": ["B"], "right": ["C"], "lineage": "L2"}])

    result = synteny.conservation(occurrences)

    assert result["modal_left"] == "A"
    assert result["lineage_left_conservation"] == 0.5


def test_operon_conservation_is_the_mean_of_within_lineage_fractions():
    """L1 = 1 operon-like of 2 copies, L2 = 1 of 1: (0.5 + 1) / 2, not 2 of 3 copies."""
    occurrences = [
        {"left": ["B"], "right": ["C"], "operon": True, "lineage": "L1"},
        {"left": ["B"], "right": ["C"], "operon": False, "lineage": "L1"},
        {"left": ["B"], "right": ["C"], "operon": True, "lineage": "L2"},
    ]

    result = synteny.conservation(occurrences)

    assert result["lineage_operon_like_conservation"] == 0.75


def test_the_minimum_number_of_lineages_is_a_parameter():
    occurrences = one_lineage_each([{"left": ["B"], "right": ["C"]}] * 2)

    assert synteny.conservation(occurrences)["status"] == "SUCCESS"
    assert synteny.conservation(occurrences, min_lineages=3)["status"] == "TOO_FEW_LINEAGES"


# ---------------------------------------------------------------------------------------
# The script: two levels, lineages, strand, origin
# ---------------------------------------------------------------------------------------

def _run(fixture_dir, genes, pmap, clusters, lineages, small, dark, dark_families,
         topology=None, window=3, levels=("close", "intermediate"), primary="intermediate"):
    """Write a Stage 9 fixture and run synteny.py on it; returns rows by family_id.

    genes: (orf_id, start, end, strand); the plasmid is the orf_id before '|'.
    clusters: {level: {member: representative}}.
    """
    plasmids = sorted({g[0].rsplit("|", 1)[0] for g in genes})
    annotation = fixture_dir / "plasmid_annotation.tsv"
    write_tsv(annotation, ["orf_id", "plasmid_id", "start", "end", "strand"],
              [[o, o.rsplit("|", 1)[0], s, e, st] for o, s, e, st in genes])
    mapping = fixture_dir / "protein_map.tsv"
    mapping.write_text("".join(f"{sid}\t{','.join(o)}\n" for sid, o in pmap.items()))
    cluster_files = []
    for level in levels:
        path = fixture_dir / f"families_{level}_cluster.tsv"
        path.write_text("".join(f"{rep}\t{mem}\n" for mem, rep in clusters[level].items()))
        cluster_files.append(str(path))
    lineage = fixture_dir / "plasmid_lineage.tsv"
    write_tsv(lineage, ["plasmid_id", "plasmid_lineage_cluster"], sorted(lineages.items()))
    small_ids = fixture_dir / "small_plasmids.txt"
    small_ids.write_text("".join(f"{p}\n" for p in small))
    dark_ids = fixture_dir / "dark_ids.txt"
    dark_ids.write_text("".join(f"{s}\n" for s in dark))
    fams = fixture_dir / "dark_families.tsv"
    write_tsv(fams, ["family_id", "members"], [[f, ""] for f in dark_families])
    registry = fixture_dir / "clonal_registry.tsv"
    write_tsv(registry, ["plasmid_id", "topology"],
              [[p, (topology or {}).get(p, "linear")] for p in plasmids])
    lengths = fixture_dir / "plasmid_lengths.tsv"
    write_tsv(lengths, ["plasmid_id", "length_bp"], [[p, 1200] for p in plasmids])
    out = fixture_dir / "synteny.tsv"

    run_script("synteny.py", FakeSnakemake(
        input={"annotation": str(annotation), "families": str(fams), "map": str(mapping),
               "clusters": cluster_files, "lineage": str(lineage),
               "dark_ids": str(dark_ids), "small_ids": str(small_ids),
               "registry": str(registry), "lengths": str(lengths)},
        output={"tsv": str(out)},
        params={"context": {"neighbourhood_window": window, "max_operon_gap": 100},
                "primary": primary,
                "synteny": {"levels": list(levels), "min_lineages": 2}}))
    return out


def _three_plasmids(fixture_dir, dark_families=("intermediate:d",), lineages=None):
    """x - d - y on P1, P2 (one lineage) and P3 (another). The left neighbours x1, x2, x3
    share an intermediate family but no close cluster; the right neighbours y1, y2, y3
    share one close cluster."""
    genes, pmap = [], {}
    for plasmid, left, right in [("P1", "x1", "y1"), ("P2", "x2", "y2"), ("P3", "x3", "y3")]:
        for k, (sid, start) in enumerate([(left, 1), ("d", 400), (right, 800)], 1):
            orf = f"{plasmid}|{k}"
            genes.append((orf, start, start + 300, 1))
            pmap.setdefault(sid, []).append(orf)
    clusters = {
        "close": {"x1": "x1", "x2": "x2", "x3": "x3", "y1": "y1", "y2": "y1", "y3": "y1",
                  "d": "d"},
        "intermediate": {"x1": "x1", "x2": "x1", "x3": "x1", "y1": "y1", "y2": "y1",
                         "y3": "y1", "d": "d"},
    }
    return _run(fixture_dir, genes, pmap, clusters,
                lineages or {"P1": "LA", "P2": "LA", "P3": "LB"},
                small=["P1", "P2", "P3"], dark=["d"], dark_families=dark_families)


def test_the_script_measures_every_dark_cluster_at_both_levels(fixture_dir):
    rows = {r["family_id"]: r for r in read_tsv(_three_plasmids(fixture_dir))}

    assert set(rows) == {"close:d", "intermediate:d"}, (
        "only clusters holding a dark small-plasmid member are measured")
    close, inter = rows["close:d"], rows["intermediate:d"]
    assert (close["level"], inter["level"]) == ("close", "intermediate")
    assert close["intermediate_family_ids"] == "intermediate:d"
    assert close["synteny_min_lineages"] == inter["synteny_min_lineages"] == "2"
    for row in (close, inter):
        assert row["n_occurrences"] == "3"
        assert row["context_recurrence"] == "3"
        assert row["n_lineages"] == "2"
        assert row["status"] == "SUCCESS"

    # Same family on the left at the intermediate level; three different close clusters.
    assert inter["lineage_left_conservation"] == "1.0"
    assert inter["modal_left"] == "intermediate:x1"
    assert inter["lineage_synteny_conservation"] == "1.0"
    assert inter["n_lineages_discordant"] == "0"
    # Close: LA splits 1/2 + 1/2 over x1 and x2, LB gives 1 to x3.
    assert close["modal_left"] == "close:x3"
    assert close["lineage_left_conservation"] == "0.5"
    assert close["lineage_right_conservation"] == "1.0"
    assert close["n_lineages_discordant"] == "1"
    assert float(close["lineage_synteny_conservation"]) < float(
        inter["lineage_synteny_conservation"])
    # Every plasmid is small here, so the small_ block repeats the measurement.
    assert close["small_lineage_left_conservation"] == "0.5"
    assert inter["small_status"] == "SUCCESS"


def test_the_script_counts_clonal_copies_as_one_lineage(fixture_dir):
    rows = {r["family_id"]: r for r in read_tsv(
        _three_plasmids(fixture_dir, lineages={"P1": "LA", "P2": "LA", "P3": "LA"}))}

    assert rows["intermediate:d"]["n_occurrences"] == "3"
    assert rows["intermediate:d"]["n_lineages"] == "1"
    assert rows["intermediate:d"]["status"] == "TOO_FEW_LINEAGES"
    assert rows["intermediate:d"]["lineage_synteny_conservation"] == ""


def test_the_script_fails_when_the_dark_set_disagrees_with_dark_families(fixture_dir):
    """Two definitions of 'dark family' must not drift apart silently."""
    with pytest.raises(SystemExit, match="dark_families"):
        _three_plasmids(fixture_dir, dark_families=("intermediate:d", "intermediate:x1"))


def test_the_script_fails_when_a_plasmid_has_no_lineage(fixture_dir):
    with pytest.raises(SystemExit, match="lineage"):
        _three_plasmids(fixture_dir, lineages={"P1": "LA", "P2": "LA"})


def test_the_small_variant_uses_small_plasmid_occurrences_only(fixture_dir):
    """S1 and S2 (small) agree on the right; L1 (large) does not."""
    genes, pmap = [], {}
    for plasmid, left, right in [("S1", "x1", "y1"), ("S2", "x2", "y2"), ("L1", "x3", "z1")]:
        for k, (sid, start) in enumerate([(left, 1), ("d", 400), (right, 800)], 1):
            orf = f"{plasmid}|{k}"
            genes.append((orf, start, start + 300, 1))
            pmap.setdefault(sid, []).append(orf)
    same = {"x1": "x1", "x2": "x1", "x3": "x1", "y1": "y1", "y2": "y1", "z1": "z1", "d": "d"}
    out = _run(fixture_dir, genes, pmap, {"close": same, "intermediate": same},
               {"S1": "A", "S2": "B", "L1": "C"}, small=["S1", "S2"], dark=["d"],
               dark_families=["intermediate:d"])

    row = {r["family_id"]: r for r in read_tsv(out)}["intermediate:d"]
    assert row["n_occurrences"] == "3" and row["small_n_occurrences"] == "2"
    assert row["lineage_left_conservation"] == "1.0"
    assert float(row["lineage_right_conservation"]) == round(2 / 3, 4)
    assert row["small_lineage_right_conservation"] == "1.0"
    assert row["small_n_lineages"] == "2"


def test_left_and_right_are_read_on_the_genes_own_strand_and_wrap_the_origin(fixture_dir):
    """x-d-y on the plus strand and y-d-x on the minus strand are one arrangement written
    in two orientations: upstream of d is x in both. On a circular record the window also
    wraps, so a gene at the record start still has a left neighbour."""
    genes = [("P1|1", 1, 300, 1), ("P1|2", 400, 700, 1), ("P1|3", 800, 1100, 1),
             ("P2|1", 1, 300, -1), ("P2|2", 400, 700, -1), ("P2|3", 800, 1100, -1)]
    pmap = {"x": ["P1|1", "P2|3"], "d": ["P1|2", "P2|2"], "y": ["P1|3", "P2|1"]}
    same = {"x": "x", "d": "d", "y": "y"}
    out = _run(fixture_dir, genes, pmap, {"close": same, "intermediate": same},
               {"P1": "A", "P2": "B"}, small=["P1", "P2"], dark=["d", "x"],
               dark_families=["intermediate:d", "intermediate:x"],
               topology={"P1": "circular", "P2": "circular"}, window=1)

    rows = {r["family_id"]: r for r in read_tsv(out)}
    assert rows["intermediate:d"]["modal_left"] == "intermediate:x"
    assert rows["intermediate:d"]["lineage_synteny_conservation"] == "1.0"
    # x is first on P1 and last on P2; on a circle its upstream neighbour is y in both.
    assert rows["intermediate:x"]["modal_left"] == "intermediate:y"
    assert rows["intermediate:x"]["lineage_left_conservation"] == "1.0"
    assert rows["close:x"]["modal_left"] == "close:y"
