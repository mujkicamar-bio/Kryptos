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
