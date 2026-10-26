"""Shared setup for the bin/ unit + contract tests.

The bin/ scripts are plain executables that Nextflow puts on PATH; here they are imported
as modules (bin/ is added to sys.path) and, for the contract tests, also run as a CLI
exactly as the EMIT_RESULTS / CNVKIT_CALL processes run them.
"""

import json
import runpy
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
BIN = REPO / "bin"
FIXTURES = Path(__file__).resolve().parent / "fixtures"

sys.path.insert(0, str(BIN))


@pytest.fixture
def fx():
    """Path to tests/python/fixtures."""
    return FIXTURES


@pytest.fixture
def run_bin(tmp_path, monkeypatch):
    """Run a bin/ script through its CLI entry point (`__main__`, real argv), the way the
    pipeline invokes it, with tmp_path as the working directory. In-process rather than a
    subprocess so coverage.py sees the code it runs."""
    def _run(script, *args):
        monkeypatch.chdir(tmp_path)
        monkeypatch.setattr(sys, "argv", [script, *map(str, args)])
        runpy.run_path(str(BIN / script), run_name="__main__")
    return _run


@pytest.fixture
def load_json():
    def _load(path):
        return json.loads(Path(path).read_text())
    return _load
