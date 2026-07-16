"""Foundation tests for the capability-matrix meta-generator.

Pins three software invariants on
:mod:`tools.verification.generate_capability_matrices`:

1. ``--check`` mode exits 0 when matrices are in sync with the
   ``cases.py:capability_rows()`` registries (regression gate
   for hand-edited matrices vs registry drift).
2. Generation is deterministic — running the meta-generator twice
   produces byte-identical output.
3. At least the two existing packages (``peierls_nystrom`` and
   ``fn_method``) are auto-discovered. Catches accidental breakage
   of the discovery loop.

These tests are tagged ``@pytest.mark.foundation`` because they
verify a software invariant of the documentation infrastructure
(the meta-generator), not a numerical solver claim — so they do not
carry a ``vv_level`` or ``verifies(...)`` mark.
"""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pytest

from tools.verification.generate_capability_matrices import (
    MATRIX_GLOB,
    MATRIX_OUTPUT_DIR,
)

REPO_ROOT = Path(__file__).resolve().parents[2]


def _run_generator(*extra_args: str) -> subprocess.CompletedProcess:
    return subprocess.run(
        [
            sys.executable,
            "-m",
            "tools.verification.generate_capability_matrices",
            *extra_args,
        ],
        cwd=REPO_ROOT,
        check=False,
        capture_output=True,
        text=True,
    )


@pytest.mark.foundation
def test_check_mode_exits_zero_when_in_sync():
    """``--check`` must report exit 0 when matrices match cases.py.

    First force a fresh write so the on-disk matrices match the
    registry, then run ``--check``. Any drift between the two would
    surface as a non-zero exit — that is exactly what this test
    forbids in the green state.
    """
    write_result = _run_generator()
    assert write_result.returncode == 0, (
        f"meta-generator write failed:\n"
        f"stdout={write_result.stdout!r}\nstderr={write_result.stderr!r}"
    )
    check_result = _run_generator("--check")
    assert check_result.returncode == 0, (
        f"--check failed after a clean write — meta-generator is "
        f"non-deterministic or there is on-disk corruption:\n"
        f"stdout={check_result.stdout!r}\nstderr={check_result.stderr!r}"
    )


@pytest.mark.foundation
def test_generation_is_deterministic():
    """Running the meta-generator twice must produce identical output.

    Reads each ``_<pkg>_capability_matrix.inc.rst`` after the first
    write, runs the generator again, and asserts byte equality.
    Catches non-deterministic dict iteration / set ordering bugs.
    """
    first_write = _run_generator()
    assert first_write.returncode == 0, first_write.stderr

    first_contents = {
        path.name: path.read_text(encoding="utf-8")
        for path in sorted(MATRIX_OUTPUT_DIR.glob(MATRIX_GLOB))
    }
    assert first_contents, "no capability matrices found after first write"

    second_write = _run_generator()
    assert second_write.returncode == 0, second_write.stderr

    second_contents = {
        path.name: path.read_text(encoding="utf-8")
        for path in sorted(MATRIX_OUTPUT_DIR.glob(MATRIX_GLOB))
    }

    assert first_contents.keys() == second_contents.keys(), (
        f"matrix file set changed between runs:\n"
        f"  first:  {sorted(first_contents)}\n"
        f"  second: {sorted(second_contents)}"
    )
    for name in first_contents:
        assert first_contents[name] == second_contents[name], (
            f"matrix {name} differs between two consecutive runs — "
            f"meta-generator is non-deterministic"
        )


@pytest.mark.foundation
def test_at_least_peierls_nystrom_and_fn_method_discovered():
    """The two existing packages must be auto-discovered.

    Asserts that running the meta-generator emits matrices for at
    least ``peierls_nystrom`` and ``fn_method``. This is the
    minimum-discovery floor — any future regression of the
    discovery loop (broken ``pkgutil.iter_modules`` traversal,
    accidental rename of ``cases.capability_rows``, etc.) trips
    this test.
    """
    write_result = _run_generator()
    assert write_result.returncode == 0, write_result.stderr

    found = {
        path.name for path in MATRIX_OUTPUT_DIR.glob(MATRIX_GLOB)
    }
    expected = {
        "_peierls_nystrom_capability_matrix.inc.rst",
        "_fn_method_capability_matrix.inc.rst",
    }
    missing = expected - found
    assert not missing, (
        f"meta-generator failed to discover expected packages: "
        f"missing={sorted(missing)}; found={sorted(found)}"
    )
