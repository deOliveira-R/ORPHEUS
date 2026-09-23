"""The harness view (``.claude/``, ``CLAUDE.md``) is GENERATED from
``docs/development/`` by ``tools/harness``. This gate is the reader of its
``--check``: a check nobody runs is not an instrument. It fails on drift between
a generated file and its source, a core over its token budget, a link to a
heading MyST would not mint, a source page whose ``harness:`` block is missing
or malformed, or a stamped file no source page produces. Run the generator to
repair drift; fix the docs page for everything else.
"""
from __future__ import annotations

import pathlib
import subprocess
import sys

import pytest

pytestmark = pytest.mark.foundation

REPO_ROOT = pathlib.Path(__file__).resolve().parents[2]


def test_generated_harness_is_current() -> None:
    run = subprocess.run(
        [sys.executable, "-m", "tools.harness", "--check"],
        cwd=REPO_ROOT, capture_output=True, text=True, check=False,
    )
    report = run.stdout + run.stderr
    assert run.returncode == 0, f"tools.harness --check failed:\n{report}"


def test_every_rule_file_is_generated_or_installed() -> None:
    """Every ``.claude/rules/*.md`` is written from ``docs/development/rules/`` or
    installed by a tool that recorded it in its manifest (``nexus setup``, the
    routing rule ``nexus-tools`` since 2026-09-21). A hand-maintained rule would
    load always-on while no source page owned it and no installer could update
    it; the last four came under generation on 2026-09-20 (K3b), so an unstamped,
    unrecorded file here is a new hand-written one. An installed file's drift
    instrument is ``nexus setup --check``, not this gate.
    """
    from tools.harness.render import STAMP_PREFIX
    from tools.harness.targets.claude_code import ClaudeCode

    installed = ClaudeCode(repo_root=REPO_ROOT).installed()
    unowned = [
        path.name
        for path in sorted((REPO_ROOT / ".claude" / "rules").glob("*.md"))
        if STAMP_PREFIX not in path.read_text(encoding="utf-8")[:1200] and path not in installed
    ]
    assert not unowned, f"hand-maintained rules (no GENERATED stamp, no installer manifest entry): {unowned}"
