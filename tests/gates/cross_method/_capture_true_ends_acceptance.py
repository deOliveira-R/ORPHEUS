r"""Freeze the #707 Step A acceptance answers into ``true_ends_acceptance.npz`` (run ONCE, on the pre-carve tree).

    .venv/bin/python -O -m tests.gates.cross_method._capture_true_ends_acceptance

Refuses to overwrite an existing fixture (the capture is the anchor that predates the carve; a move is a declared
event in ``true_ends_acceptance_rebaselines.json``, never a re-capture). Records the commit, ``git diff --stat`` of
``orpheus/`` against ``main`` (empty, or derivations only, on the capture tree), numpy and the platform.
"""

from __future__ import annotations

import platform
import subprocess
import sys
import time
from typing import Any, cast

import numpy as np

from tests.gates.cross_method import _true_ends_acceptance_cases as cases


def _git(*args: str) -> str:
    return subprocess.check_output(["git", *args], text=True).strip()


def main() -> None:
    if cases.FIXTURE.exists() and "--force" not in sys.argv:
        raise SystemExit(f"{cases.FIXTURE} exists: the acceptance anchor is never re-captured (declare an event)")
    import orpheus

    payload: dict[str, np.ndarray] = {}
    for case in cases.cases():
        t0 = time.perf_counter()
        for field, value in case.solve().items():
            payload[f"{case.name}/{field}"] = np.asarray(value)
        print(f"{case.name}: {time.perf_counter() - t0:.2f} s", flush=True)
    payload["meta/commit"] = np.array(_git("rev-parse", "HEAD"))
    payload["meta/orpheus_diff_vs_main"] = np.array(_git("diff", "--stat", "main", "--", "orpheus"))
    payload["meta/orpheus_file"] = np.array(orpheus.__file__)
    payload["meta/numpy"] = np.array(np.__version__)
    payload["meta/platform"] = np.array(platform.platform())
    np.savez(cases.FIXTURE, **cast("dict[str, Any]", payload))
    print(f"wrote {cases.FIXTURE}: {len(payload)} keys")


if __name__ == "__main__":
    main()
