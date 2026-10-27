"""nextflow_schema.json must stay in lock-step with the pipeline.

Seqera Platform builds its launch form from nextflow_schema.json, so drift is silent and
user-facing: a parameter missing from the schema can't be set from the form, one that is
in the schema but not the pipeline is ignored, and a wrong default misleads. These tests
tie the schema to nextflow.config, main.nf's validation and the web launcher.
"""

import ast
import csv
import json
import re

import pytest
from jsonschema import Draft202012Validator

from conftest import REPO

SCHEMA = json.loads((REPO / "nextflow_schema.json").read_text())
INPUT_SCHEMA = json.loads((REPO / "assets" / "schema_input.json").read_text())


def _config_params():
    """name -> default from the top-level `params { }` block of nextflow.config."""
    text = (REPO / "nextflow.config").read_text()
    block = re.search(r"^params \{\n(.*?)^\}", text, re.S | re.M).group(1)
    params = {}
    for line in block.splitlines():
        m = re.match(r"\s*(\w+)\s*=\s*(.+?)\s*(//.*)?$", line)
        if not m:
            continue
        raw = m.group(2)
        if raw[0] in "'\"":
            value = raw[1:-1]
        elif raw in ("true", "false"):
            value = raw == "true"
        elif raw == "null":
            value = None
        else:
            value = float(raw) if "." in raw else int(raw)
        params[m.group(1)] = value
    return params


def _schema_props():
    props = {}
    for group in SCHEMA["$defs"].values():
        props.update(group["properties"])
    return props


CONFIG = _config_params()
PROPS = _schema_props()


def test_schemas_are_valid_json_schema():
    Draft202012Validator.check_schema(SCHEMA)
    Draft202012Validator.check_schema(INPUT_SCHEMA)


def test_every_group_is_referenced():
    refs = {r["$ref"].split("/")[-1] for r in SCHEMA["allOf"]}
    assert refs == set(SCHEMA["$defs"])


def test_schema_and_config_declare_the_same_params():
    assert set(PROPS) == set(CONFIG), (
        f"only in schema: {set(PROPS) - set(CONFIG)}; only in nextflow.config: {set(CONFIG) - set(PROPS)}")


@pytest.mark.parametrize("name", sorted(CONFIG))
def test_schema_default_matches_config(name):
    if CONFIG[name] is None:
        assert "default" not in PROPS[name], f"{name}: config default is null, schema has one"
    else:
        assert PROPS[name].get("default") == CONFIG[name]


@pytest.mark.parametrize("name", sorted(PROPS))
def test_every_param_has_type_and_description(name):
    assert PROPS[name].get("type") and PROPS[name].get("description")


def test_variant_caller_choices_match_main_nf():
    main_nf = (REPO / "main.nf").read_text()
    allowed = ast.literal_eval(re.search(r"params\.variant_caller in (\[[^\]]+\])", main_nf).group(1))
    assert sorted(PROPS["variant_caller"]["enum"]) == sorted(allowed)


def test_default_image_is_this_release():
    """The default container is the image published for this release (pushed when the vX.Y.Z
    git tag is pushed). Bumping manifest.version without it would run the previous image."""
    version = re.search(r"^\s*version\s*=\s*'([^']+)'", (REPO / "nextflow.config").read_text(),
                        re.M).group(1)
    assert CONFIG["strainqc_container"] == f"ghcr.io/digibio-ut/strainqc:{version}"


def test_launcher_only_forwards_known_params():
    """The web launcher passes UI options as --<key>; an unknown key would be silently ignored."""
    tree = ast.parse((REPO / "launcher" / "app" / "runner.py").read_text())
    forwarded = next(ast.literal_eval(n.value) for n in ast.walk(tree)
                     if isinstance(n, ast.Assign) and any(getattr(t, "id", "") == "_FORWARD_OPTS"
                                                          for t in n.targets))
    assert forwarded <= set(PROPS), f"launcher forwards unknown params: {forwarded - set(PROPS)}"


# ── samplesheet schema ─────────────────────────────────────────────────────────

def _rows(name):
    with open(REPO / "assets" / "test" / name, newline="") as fh:
        return [{k: v for k, v in r.items() if v} for r in csv.DictReader(fh)]   # empty cell = absent


def test_required_columns_match_main_nf_check():
    assert sorted(INPUT_SCHEMA["items"]["required"]) == ["fastq_1", "sample"]
    assert "containsKey('sample')" in (REPO / "main.nf").read_text()


@pytest.mark.parametrize("sheet", ["samplesheet.csv", "samplesheet_two.csv", "samplesheet_three.csv"])
def test_bundled_samplesheets_validate(sheet):
    rows = _rows(sheet)
    Draft202012Validator(INPUT_SCHEMA).validate(rows)
    ids = [r["sample"] for r in rows]
    assert len(ids) == len(set(ids))                         # uniqueEntries: ["sample"]


def test_duplicate_ids_violate_unique_entries():
    ids = [r["sample"] for r in _rows("samplesheet_dupes.csv")]
    assert INPUT_SCHEMA["uniqueEntries"] == ["sample"] and len(ids) != len(set(ids))


def test_sample_ids_with_spaces_are_rejected():
    bad = [{"sample": "my sample", "fastq_1": "a.fastq.gz"}]
    assert list(Draft202012Validator(INPUT_SCHEMA).iter_errors(bad))
