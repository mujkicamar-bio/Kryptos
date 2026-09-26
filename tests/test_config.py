"""The configuration file is the single source of every threshold (spec §6, §15.3).

These tests exist because three of six thresholds in an earlier version of this project
lived as Python constants, which meant they could not be validated, swept, or recorded in
the output. A threshold that is not in config is a threshold nobody can audit.
"""
import pathlib

import jsonschema
import yaml

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
    """Spec §15.3. narrow_at controls compute; min_explained controls science.

    Read from the cascade files the workflow loads (Snakefile: cascade_config), in
    production and in the smoke and benchmark configurations alike. narrow_at 0.7 is a
    user decision (2026-09-25)."""
    for path in (ROOT / "config" / "cascade.yaml", ROOT / "config" / "test" / "cascade.yaml",
                 ROOT / "config" / "bench" / "cascade.yaml"):
        cascade = yaml.safe_load(path.read_text())
        assert cascade["narrow_at"] == 0.7, path
        assert cascade["min_explained"] == 0.5, path
        assert cascade["narrow_at"] >= cascade["min_explained"], path


def test_the_cascade_thresholds_live_only_in_the_cascade_files():
    """A second copy of narrow_at in config.yaml, which the workflow never read, once
    stated 0.9 while the run used 0.7. Only the cascade files may carry these keys."""
    for path in (ROOT / "config" / "config.yaml", ROOT / "config" / "test" / "config.yaml",
                 ROOT / "config" / "bench" / "config.yaml"):
        assert "cascade" not in yaml.safe_load(path.read_text()), path


def test_targets_match_their_schema():
    targets = yaml.safe_load((ROOT / "config" / "targets.yaml").read_text())
    schema = yaml.safe_load(
        (ROOT / "config" / "schemas" / "targets.schema.yaml").read_text())
    jsonschema.validate(targets, schema)


def test_the_test_configuration_validates_against_the_same_schema():
    """config/test is what the smoke run uses; a key the production config gained and the
    test config lacks fails the test run at load time."""
    config = yaml.safe_load((ROOT / "config" / "test" / "config.yaml").read_text())
    schema = yaml.safe_load((ROOT / "config" / "schemas" / "config.schema.yaml").read_text())
    jsonschema.validate(config, schema)


def test_the_gpu_targets_differ_from_targets_only_in_structure_gpu():
    """targets.gpu.yaml files are copies of targets.yaml with structure.gpu enabled. A key
    added to one and not the others makes the GPU run use different thresholds."""
    targets = yaml.safe_load((ROOT / "config" / "targets.yaml").read_text())
    for path in (ROOT / "config" / "targets.gpu.yaml",
                 ROOT / "config" / "test" / "targets.gpu.yaml"):
        gpu = yaml.safe_load(path.read_text())
        assert gpu["structure"].pop("gpu") is True
        expected = {**targets, "structure": {k: v for k, v in targets["structure"].items()
                                             if k != "gpu"}}
        assert gpu == expected, f"{path.name} drifted from targets.yaml"


def test_the_label_databases_and_conjscan_are_configured():
    config = yaml.safe_load((ROOT / "config" / "config.yaml").read_text())
    targets = yaml.safe_load((ROOT / "config" / "targets.yaml").read_text())
    assert config["labels"]["dir"] == "data/refs/labels"
    assert config["amrfinder"]["executable"] == "envs/amrfinder/bin/amrfinder"
    assert config["references"]["conjscan_models"] == "data/refs/conjscan"
    # The KO list maps Tier 0 KOs to symbols for the disagreement table only.
    assert config["references"]["kegg_ko_list"].startswith("data/refs/labels/kegg_ko/")
    assert targets["conjugation"]["exe"] == "envs/conjscan/bin/macsyfinder"
    assert targets["conjugation"]["version"] == "2.1.0"
    # The primary resolution must be a synteny level (the family table reads its rows),
    # and close is the level the ORF table reads.
    assert targets["clustering"]["primary"] in targets["synteny"]["levels"]
    assert "close" in targets["synteny"]["levels"]
    assert targets["synteny"]["min_lineages"] == 2
