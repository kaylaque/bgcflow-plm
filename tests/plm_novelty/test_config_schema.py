"""Tests for workflow/schemas/config.schema.yaml — plm-novelty boolean toggle."""
from pathlib import Path

import jsonschema
import pytest
import yaml

SCHEMA_PATH = Path(__file__).parents[2] / "workflow" / "schemas" / "config.schema.yaml"


@pytest.fixture
def schema():
    with open(SCHEMA_PATH) as fh:
        return yaml.safe_load(fh)


def _minimal_config(extra_rules=None):
    """Return a minimal config dict that satisfies required fields."""
    cfg = {
        "projects": [{"name": "test_proj", "samples": "samples.csv", "rules": "rules.yaml"}],
        "resources_path": {"antismash_db": "resources/antismash_db"},
    }
    if extra_rules:
        cfg["rules"] = extra_rules
    return cfg


def test_plm_novelty_true_validates(schema):
    """plm-novelty: true should pass schema validation."""
    cfg = _minimal_config({"plm-novelty": True})
    jsonschema.validate(instance=cfg, schema=schema)  # must not raise


def test_plm_novelty_false_validates(schema):
    """plm-novelty: false should also pass."""
    cfg = _minimal_config({"plm-novelty": False})
    jsonschema.validate(instance=cfg, schema=schema)


def test_plm_novelty_string_fails(schema):
    """plm-novelty: 'yes' (string) must fail validation."""
    cfg = _minimal_config({"plm-novelty": "yes"})
    with pytest.raises(jsonschema.ValidationError):
        jsonschema.validate(instance=cfg, schema=schema)
