r"""G00a: the stored bytes. Every array of ``snapshots/MANIFEST.txt`` recomputed on the live tree equals its stored bytes.

This file reads bit identity against bytes stored on one host (the stamp in the first line of the manifest). On any other host the
gate FAILS and says so, because a skip would be a silent green. It is not listed in the CI step: CI runs on ``ubuntu-latest`` and
the snapshots were taken on macOS arm64 with Accelerate. Its sibling ``test_snapshot_controls.py`` compares a capture with a
capture in one process and runs anywhere.

Every assertion is here and none in ``_snapshot.py``: a helper's ``assert`` is removed by ``python -O``.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from tests.gates.cross_method.w1 import _snapshot as snap

MANIFEST = snap.SNAPSHOT_DIR / "MANIFEST.txt"


def _require(cond: object, msg: str) -> None:
    if not cond:
        pytest.fail(msg)


def check_every_snapshot_equal(manifest: Path = MANIFEST) -> None:
    """Fail unless the live capture equals the stored arrays and the manifest lines, on the stamped host."""
    stamp, lines = snap.read_manifest(manifest)
    live_stamp = snap.host_stamp()
    _require(
        stamp == live_stamp,
        f"bit identity cannot be read on this host: the manifest was taken on {stamp} and this host is {live_stamp}. "
        f"Re-take the snapshots on this host at the last green commit (a commit of its own) before comparing.",
    )
    live = snap.capture()
    stored = snap.load_snapshots(Path(manifest).parent)
    differing = snap.differing_keys(live, stored)
    _require(not differing, f"{len(differing)} of {len(live)} arrays differ from the stored bytes: {differing[:10]}")
    _require(
        snap.manifest_lines(live) == lines,
        "the manifest lines (shape, dtype, sha1 per array) differ from the live capture although the stored arrays agree",
    )


def test_every_snapshot_equal():
    check_every_snapshot_equal()
