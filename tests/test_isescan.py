"""ISEScan output parsing, against rows taken from a real run on the test plasmids."""
import pathlib

from plasmidann.isescan import COLUMNS, parse_isescan

SAMPLE = pathlib.Path(__file__).parent / "data" / "isescan_sample.tsv"


def test_every_element_is_parsed_with_its_family_and_boundaries():
    rows = parse_isescan(SAMPLE)
    assert [r["family"] for r in rows] == ["IS5", "IS21", "IS3", "IS21"]
    assert (rows[0]["start"], rows[0]["end"]) == (33252, 35386)
    assert all(set(r) == set(COLUMNS) for r in rows)


def test_complete_and_partial_are_kept_apart():
    """Partial elements are reported, as ISEScan's defaults do: a partial IS on a plasmid
    is still an IS-derived region, and dropping it would hide exactly the degenerate
    transposase fragments that otherwise look like dark ORFs."""
    assert [r["complete"] for r in parse_isescan(SAMPLE)] == [1, 1, 0, 0]


def test_absent_values_stay_absent():
    """ISEScan writes '-:-' for an element with no TIR pair and leaves strand empty for
    some elements. Both must read as absence, never as a TIR called '-:-'."""
    rows = parse_isescan(SAMPLE)
    assert rows[2]["tir"] == ""
    assert rows[3]["strand"] == ""
    assert rows[1]["strand"] == "-"


def test_ids_are_numbered_per_plasmid():
    ids = [r["is_id"] for r in parse_isescan(SAMPLE)]
    assert ids[0] == "COMPASS_NZ_CP010683.1|IS1"
    assert len(set(ids)) == len(ids)
