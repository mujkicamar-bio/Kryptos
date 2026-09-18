"""Stage 13: explicit backgrounds (spec section 52).

Section 52.2 is an absolute: "No normalized prevalence field is valid without its
normalization definition." Section 53 requires context enrichment to be compared against an
appropriate background rather than a flat corpus rate, because on a 5 kb plasmid carrying
six genes a plus or minus three neighbourhood is the whole molecule - and small cryptic
plasmids are a stratum of interest here, so a flat background errs exactly where it does
most damage.
"""
from plasmidann import background


def test_a_normalised_value_never_travels_without_its_definition():
    """Section 52.2, asserted. A normalised number on its own is not a result: its meaning
    depends entirely on a choice the reader cannot see."""
    result = background.normalise(
        observed_count=8, observed_total=10,
        background_count=100, background_total=1000,
        method=background.METHOD_STRATIFIED,
        background_definition="stratified on plasmid_length = 4kb; 40 plasmids")

    for required in ("normalization_method", "normalization_version",
                     "background_definition"):
        assert result[required], f"{required} is empty beside a normalised value"
    assert result["normalized_prevalence"] == 8.0


def test_the_raw_values_are_always_retained():
    """Section 52.1. Normalisation must never be the only thing recorded, or the reader
    cannot recompute it under a different background."""
    result = background.normalise(3, 12, 50, 500, "m", "d")

    assert result["raw_count"] == 3
    assert result["raw_prevalence"] == 0.25
    assert result["background_prevalence"] == 0.1


def test_a_zero_background_gives_no_normalised_value_rather_than_infinity():
    """An infinite ratio is not a measurement. A reader ranking on it would put every
    never-before-seen association above every well-supported one."""
    result = background.normalise(5, 10, 0, 1000, "m", "d")

    assert result["normalized_prevalence"] == ""
    assert result["raw_prevalence"] == 0.5, "the raw value must survive"


def test_length_bands_are_logarithmic():
    """Plasmid size spans four orders of magnitude. A linear band wide enough to be useful
    at 200 kb would put every cryptic plasmid in one bucket."""
    small = background.length_band(3_000)
    also_small = background.length_band(3_500)
    large = background.length_band(200_000)

    assert small == also_small, "two cryptic plasmids landed in different bands"
    assert small != large


def test_a_missing_covariate_does_not_merge_a_plasmid_into_another_stratum():
    """Dropping an absent covariate would pool plasmids with missing metadata into the
    strata of plasmids that have it, which makes metadata completeness look like biology."""
    known = background.stratum_key({"size_bp": 5000, "hab_top": "Host-associated"},
                                   ["plasmid_length", "hab_top"])
    unknown = background.stratum_key({"size_bp": 5000},
                                     ["plasmid_length", "hab_top"])

    assert known != unknown
    assert "unknown" in unknown


def test_plasmids_group_into_strata_by_the_configured_covariates():
    plasmids = {
        "p1": {"size_bp": 3000, "hab_top": "Host-associated"},
        "p2": {"size_bp": 3200, "hab_top": "Host-associated"},
        "p3": {"size_bp": 200000, "hab_top": "Host-associated"},
    }

    strata = background.build_background(plasmids, ["plasmid_length", "hab_top"])

    sizes = sorted(len(v) for v in strata.values())
    assert sizes == [1, 2], f"expected a pair and a singleton, got {sizes}"


def test_the_background_definition_says_whether_stratification_was_used():
    """A reader must be able to tell whether a comparison rests on forty plasmids or four,
    and whether the stratified background was used at all."""
    stratified = background.describe_background(
        ["plasmid_length"], ("4kb",), n_plasmids=40, fell_back=False)
    pooled = background.describe_background(
        ["plasmid_length"], ("4kb",), n_plasmids=1000, fell_back=True)

    assert "stratified" in stratified and "40" in stratified
    assert "pooled" in pooled and "below min_stratum_size" in pooled


def test_the_normalization_version_is_recorded():
    """Two runs whose stratification differs are not comparable, and nothing else in the
    row would say so."""
    result = background.normalise(1, 2, 1, 2, "m", "d")

    assert result["normalization_version"] == background.NORMALIZATION_VERSION
