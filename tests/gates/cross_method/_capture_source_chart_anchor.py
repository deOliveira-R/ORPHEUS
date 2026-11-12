r"""Freeze the density-chart anchor of the #707 Step A chart law (run ONCE, on the pre-carve tree).

    .venv/bin/python -O -m tests.gates.cross_method._capture_source_chart_anchor

Refuses to overwrite: the anchor predates the carve, which is what makes ``new == m ⊙ anchor`` a law and not a
re-capture (``retirement-audit`` D.14: a snapshot re-captured by the carve it gates pins nothing).
"""

from __future__ import annotations

import platform
import subprocess
import sys
import warnings
from typing import Any, cast

import numpy as np

from tests.gates.cross_method import _source_chart_cases as cases


def main() -> None:
    if cases.ANCHOR.exists() and "--force" not in sys.argv:
        raise SystemExit(f"{cases.ANCHOR} exists: the anchor is never re-captured")
    import orpheus

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        payload = cases.everything()
    payload["meta/commit"] = np.array(subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip())
    payload["meta/orpheus_diff_vs_main"] = np.array(
        subprocess.check_output(["git", "diff", "--stat", "main", "--", "orpheus"], text=True).strip())
    payload["meta/orpheus_file"] = np.array(orpheus.__file__)
    payload["meta/numpy"] = np.array(np.__version__)
    payload["meta/platform"] = np.array(platform.platform())
    np.savez(cases.ANCHOR, **cast("dict[str, Any]", payload))
    print(f"wrote {cases.ANCHOR}: {len(payload)} keys")


if __name__ == "__main__":
    main()
