"""The ``.claude/`` harness view is GENERATED from ``docs/development/`` by
``tools/docs/generate_harness.py``. This gate is the reader of its ``--check``:
a check nobody runs is not an instrument. It fails on drift between a generated
file and its source, a core over its token budget, a link to a heading MyST
would not mint, a manifest entry whose target lies outside ``.claude/`` (the
repo-root ``CLAUDE.md`` excepted: its on-boarding block is generated), or a
stamped file no manifest entry produces. Run the generator to repair drift; fix
the docs page for everything else.
"""
from __future__ import annotations

import pathlib
import subprocess
import sys

import pytest

pytestmark = pytest.mark.foundation

REPO_ROOT = pathlib.Path(__file__).resolve().parents[1]


def test_generated_harness_is_current() -> None:
    run = subprocess.run(
        [sys.executable, "-m", "tools.docs.generate_harness", "--check"],
        cwd=REPO_ROOT, capture_output=True, text=True, check=False,
    )
    report = run.stdout + run.stderr
    assert run.returncode == 0, f"tools.docs.generate_harness --check failed:\n{report}"
