"""The configuration files: schema validity, and relations between values that the schemas
cannot express."""
import pathlib
import re

import jsonschema
import pytest
import yaml
from snakemake.utils import update_config

ROOT = pathlib.Path(__file__).resolve().parents[1]
CONFIG = ROOT / "config"


def load(path):
    return yaml.safe_load(path.read_text())


def test_the_three_configurations_validate_against_the_schema():
    """The Snakefile always loads config/config.yaml and merges a --configfile over it, so
    the test and benchmark configurations are validated as merged."""
    schema = load(CONFIG / "schemas" / "config.schema.yaml")
    base = load(CONFIG / "config.yaml")
    jsonschema.validate(base, schema)
    for path in (CONFIG / "test" / "config.yaml", CONFIG / "bench" / "config.yaml"):
        merged = load(CONFIG / "config.yaml")
        update_config(merged, load(path))
        jsonschema.validate(merged, schema)


def test_a_misspelt_config_key_is_refused():
    schema = load(CONFIG / "schemas" / "config.schema.yaml")
    config = load(CONFIG / "config.yaml")
    config["input"]["max_plasmid_size"] = config["input"].pop("max_plasmid_size_bp")
    try:
        jsonschema.validate(config, schema)
    except jsonschema.ValidationError:
        return
    raise AssertionError("a misspelt input key passed the schema")


def test_the_smoke_and_bench_cascades_use_the_production_thresholds():
    production = load(CONFIG / "cascade.yaml")
    assert production["narrow_at"] >= production["min_coverage"]
    for path in (CONFIG / "test" / "cascade.yaml", CONFIG / "bench" / "cascade.yaml"):
        cascade = load(path)
        for key in ("narrow_at", "min_coverage", "full_at", "partial_at",
                    "max_target_seqs", "search_clustering",
                    "artefact_screen"):
            assert cascade[key] == production[key], f"{path}: {key}"
        # Tiers T1-T4 are searched identically; only their notes may differ.
        def tiers(c):
            return [{k: v for k, v in t.items() if k != "note"} for t in c["tiers"][:4]]
        assert tiers(cascade) == tiers(production), f"{path}: tiers T1-T4"


def test_the_targets_file_validates_against_the_schema():
    jsonschema.validate(load(CONFIG / "targets.yaml"),
                        load(CONFIG / "schemas" / "targets.schema.yaml"))


def test_a_misspelt_targets_or_cascade_key_is_refused():
    """An optional key read with a default would otherwise be ignored when misspelt."""
    for name, section, key in (("targets", "structure", "scop"),
                               ("targets", "rarity", "rare_max_lineage"),
                               ("cascade", "search_clustering", "min_seqid")):
        config = load(CONFIG / f"{name}.yaml")
        config[section][key] = 1
        with pytest.raises(jsonschema.ValidationError):
            jsonschema.validate(config, load(CONFIG / "schemas" / f"{name}.schema.yaml"))


def test_the_report_levels_are_synteny_levels():
    """The family table reads the primary-level synteny rows, the ORF table the close
    rows."""
    targets = load(CONFIG / "targets.yaml")
    assert targets["clustering"]["primary"] in targets["synteny"]["levels"]
    assert "close" in targets["synteny"]["levels"]


def test_the_installer_installs_the_configured_conjscan_release():
    (installed,) = re.findall(r'^CONJSCAN_VERSION = "([^"]+)"',
                              (ROOT / "tools" / "install_tool_envs.py").read_text(), re.M)
    assert installed == load(CONFIG / "targets.yaml")["conjugation"]["version"]
