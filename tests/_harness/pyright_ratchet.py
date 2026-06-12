"""Pyright error ratchet — the #226 burn-down instrument.

The package carries a large pre-existing type-error surface
(ORPHEUS#226: 734 errors at filing). Reaching zero is a per-module
burn-down; until then the enforceable invariant is monotonicity:
**no change may increase any module's pyright error count.** This
module is the single source for both sides of that contract:

- ``collect_module_counts`` — run pyright, bucket errors by top-level
  ``orpheus`` subpackage (the burn-down axis from the issue).
- ``python -m tests._harness.pyright_ratchet --update`` — regenerate
  the checked-in baseline after an intentional improvement.
- ``tests/test_pyright_ratchet.py`` — compare live counts against the
  baseline (both directions: increase = regression, decrease =
  baseline must be tightened so the ratchet stays taut).

The baseline records the pyright version it was measured with:
version drift can move counts without any code change, so the test
surfaces the recorded version in its failure message as a triage
hint (it does NOT skip on mismatch — a silent skip would un-ratchet
every machine with a newer pyright).

Faithful-analysis recipe (validated in #226): per-checkout venv +
run pyright FROM the checkout root, so worktree sessions never see
main-rooted cross-tree artifacts.
"""

from __future__ import annotations

import json
import shutil
import subprocess
import sys
from collections import Counter
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
BASELINE_PATH = Path(__file__).with_name("pyright_baseline.json")
_ROOT_BUCKET = "(root)"
UPDATE_CMD = "python -m tests._harness.pyright_ratchet --update"


def find_pyright() -> str | None:
    """The host pyright (npm/homebrew per README setup), if any."""
    return shutil.which("pyright")


def collect_module_counts(pyright: str) -> tuple[dict[str, int], str]:
    """Run pyright over ``orpheus/`` and bucket errors by top-level
    subpackage. Returns ``(counts, pyright_version)``.

    Buckets only diagnostics of severity ``error`` — warnings are not
    ratcheted (they fluctuate with pyright releases far more than
    errors do).
    """
    proc = subprocess.run(
        [pyright, "--outputjson", "orpheus/"],
        cwd=REPO_ROOT,
        capture_output=True,
        text=True,
        timeout=600,
    )
    # pyright exits 1 when it finds errors — that IS the expected
    # state during the burn-down; only refuse unparseable output.
    report = json.loads(proc.stdout)

    counts: Counter[str] = Counter()
    for diag in report["generalDiagnostics"]:
        if diag.get("severity") != "error":
            continue
        rel = Path(diag["file"]).resolve().relative_to(REPO_ROOT / "orpheus")
        bucket = rel.parts[0] if len(rel.parts) > 1 else _ROOT_BUCKET
        counts[bucket] += 1
    return dict(sorted(counts.items())), report["version"]


def read_baseline() -> tuple[dict[str, int], str]:
    """Returns ``(per-module counts, pyright version at measurement)``."""
    data = json.loads(BASELINE_PATH.read_text())
    return data["modules"], data["pyright_version"]


def write_baseline(counts: dict[str, int], version: str) -> None:
    BASELINE_PATH.write_text(json.dumps(
        {
            "_comment": (
                "Pyright error ratchet baseline (ORPHEUS#226). "
                "Counts must only go DOWN. Regenerate after an "
                f"intentional improvement with: {UPDATE_CMD}"
            ),
            "pyright_version": version,
            "total": sum(counts.values()),
            "modules": counts,
        },
        indent=2,
    ) + "\n")


def main(argv: list[str] | None = None) -> int:
    args = argv if argv is not None else sys.argv[1:]
    pyright = find_pyright()
    if pyright is None:
        print("pyright not found on PATH (npm install -g pyright)", file=sys.stderr)
        return 1
    counts, version = collect_module_counts(pyright)
    if "--update" in args:
        write_baseline(counts, version)
        print(f"baseline updated: {sum(counts.values())} errors "
              f"across {len(counts)} modules (pyright {version})")
        return 0
    print(json.dumps({"pyright_version": version, "modules": counts}, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
