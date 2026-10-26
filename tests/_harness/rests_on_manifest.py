"""pytest plugin: write every collected item's node id and its ``rests_on`` targets to a JSON file.

The ``rests_on`` resolution gate (``tests/gates/test_test_layers.py``) runs
``pytest --collect-only -p tests._harness.rests_on_manifest tests`` in a
subprocess and reads the file. The marker arguments are read from pytest's
resolved marks (decorators, class marks, module ``pytestmark``), so a target
spelled as an f-string or a module constant arrives evaluated. The plugin
writes nothing unless ``ORPHEUS_RESTS_ON_MANIFEST`` names the output path.

It lives under ``tests/_harness/`` and not under ``tests/gates/``, because the
test-layer gate scans every file under ``tests/gates/``.
"""
from __future__ import annotations

import json
import os
import pathlib

import pytest

#: The environment variable naming the output path.
MANIFEST_ENV = "ORPHEUS_RESTS_ON_MANIFEST"


def pytest_collection_finish(session: pytest.Session) -> None:
    out = os.environ.get(MANIFEST_ENV)
    if not out:
        return
    rows = [
        {
            "nodeid": item.nodeid,
            "rests_on": [str(arg) for mark in item.iter_markers("rests_on") for arg in mark.args],
        }
        for item in session.items
    ]
    pathlib.Path(out).write_text(json.dumps(rows), encoding="utf-8")
