"""The configuration file is the single source of every threshold (spec §6, §15.3).

These tests exist because three of six thresholds in an earlier version of this project
lived as Python constants, which meant they could not be validated, swept, or recorded in
the output. A threshold that is not in config is a threshold nobody can audit.
"""
import pathlib
import yaml
import jsonschema

ROOT = pathlib.Path(__file__).resolve().parents[1]


def test_config_validates_against_its_schema():
    config = yaml.safe_load((ROOT / "config" / "config.yaml").read_text())
    schema = yaml.safe_load((ROOT / "config" / "schemas" / "config.schema.yaml").read_text())
    jsonschema.validate(config, schema)


def test_the_two_length_parameters_are_distinct_keys():
    """Spec §11.1: the gene-caller floor and the discovery cutoff are different concepts
    that happen to share a value. Conflating them is how a discovery flag silently becomes
    a deletion."""
    config = yaml.safe_load((ROOT / "config" / "config.yaml").read_text())
    assert config["orf"]["min_call_length_aa"] == 20
    assert config["discovery"]["min_dark_candidate_length_aa"] == 20


def test_cascade_thresholds_match_the_spec():
    """Spec §15.3. narrow_at controls compute; min_explained controls science."""
    config = yaml.safe_load((ROOT / "config" / "config.yaml").read_text())
    assert config["cascade"]["narrow_at"] == 0.9
    assert config["cascade"]["min_explained"] == 0.5
    assert config["cascade"]["narrow_at"] >= config["cascade"]["min_explained"]


def test_the_background_stratification_is_configured():
    """Spec section 53: context enrichment is compared against an appropriate background,
    and a stratified one is what this collection needs.

    A flat corpus background divides by one number for the whole collection. On a 5 kb
    cryptic plasmid a plus or minus three neighbourhood is the entire molecule, so a flat
    background under-corrects for small plasmids and over-corrects for large ones - and
    small cryptic plasmids are a stratum of interest here, which puts the error exactly
    where it does most damage.

    Which covariates to stratify on is a run-time choice, not a constant: each added
    covariate makes every stratum smaller, and below min_stratum_size the background is
    too thin to estimate anything from.
    """
    targets = yaml.safe_load((ROOT / "config" / "targets.yaml").read_text())
    schema = yaml.safe_load(
        (ROOT / "config" / "schemas" / "targets.schema.yaml").read_text())
    jsonschema.validate(targets, schema)

    background = targets["background"]
    assert background["covariates"], (
        "an empty covariate list is a pooled background wearing a stratified label")
    assert background["min_stratum_size"] >= 2, (
        "a stratum of one plasmid is that plasmid; it cannot be its own background")
