r"""G01i: the declared signature ``(n - kernel.dim, kernel.dim, 0)`` holds on the formal path and is checked on the raw path (N01, spec 2.1 S2, 3.3, 7.2).

Leg 1, the raw path: ``DenseMetric(matrix, kernel=Coord(n, ()))`` on an indefinite matrix (``diag([1, -1])``) or one with a kept exact
zero (``diag([1, 0])``) is refused with "kept eigenvalue is not positive" ONLY under ``if __debug__:``, so the refusal is read in a
subprocess WITHOUT ``-O`` (the canonical ``python -O`` strips the block, so it cannot be asserted in the gate's own process). Run
with ``-O`` the same code prints ``accepted`` for both, and the faces of the second are NaN (a kept exact zero divides by zero: the
raw path's own business; the check is a misuse check, not a guarantee, ruling of 2026-10-30).

Leg 2, the formal path, in the gate's own process under ``-O``: for each of the 17 frame-built Grams, with ``kernel.dim`` from
``frame_gram_kernel`` LIVE, the sign of ``u_i . vt_i`` of ``np.linalg.svd(g_sym, hermitian=True)`` over the kept directions (all but
``kernel.dim``) has no entry ``<= 0``. The sign test is computed HERE, not by the production body. The guarantee assumes
non-negative quadrature weights, which N03 enforces (#706) and N01 does not.
"""

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

from tests.gates.numerics.w1 import _stored_grams as sg

pytestmark = pytest.mark.foundation

_ROOT = Path(__file__).resolve().parents[4]

CODE = r"""
import numpy as np
from orpheus.numerics.declaration import Coord, Spanned
from orpheus.numerics.metric import DenseMetric
for m, k in (
    (np.diag([1.0, -1.0]), Coord(2, ())),
    (np.diag([1.0, 0.0]), Coord(2, ())),
    (np.diag([1.0, 0.0]), Spanned(2, ("k",), 1)),
):
    try:
        d = DenseMetric(m, kernel=k)
    except ValueError as e:
        print("refused: " + str(e))
    else:
        print("accepted nan=" + str(bool(np.isnan(np.asarray(d.inverse_matrix)).any())))
"""


def _require(cond: object, msg: str) -> None:
    if not cond:
        pytest.fail(msg)


def _run(*flags: str) -> list[str]:
    env = {k: v for k, v in os.environ.items() if k != "PYTHONOPTIMIZE"}
    env["PYTHONPATH"] = os.pathsep.join([str(_ROOT), *([env["PYTHONPATH"]] if env.get("PYTHONPATH") else [])])
    r = subprocess.run([sys.executable, *flags, "-c", CODE], env=env, cwd=_ROOT, capture_output=True, text=True, timeout=120)
    _require(r.returncode == 0, f"the subprocess ({flags}) exited {r.returncode}: {r.stderr[-400:]}")
    return r.stdout.strip().splitlines()


def refusal_reads(lines: list[str]) -> None:
    """Raise unless the first two lines are the refusal "kept eigenvalue is not positive" and the third (a legal deficient form: the
    kernel is declared) is ``accepted nan=False`` (reads the printed messages, whole)."""
    if len(lines) != 3 or not all(l.startswith("refused: ") and "kept eigenvalue is not positive" in l for l in lines[:2]):
        raise AssertionError(f"the raw path did not refuse both matrices with 'kept eigenvalue is not positive': {lines}")
    if lines[2] != "accepted nan=False":
        raise AssertionError(f"the legal deficient form diag([1, 0]) with its kernel declared was not accepted: {lines[2]}")


def test_raw_path_refusal_under_debug() -> None:
    refusal_reads(_run())  # no -O: the check stands


def test_without_debug_the_raw_path_accepts_and_a_kept_zero_gives_nan_faces() -> None:
    lines = _run("-O")
    _require(lines == ["accepted nan=False", "accepted nan=True", "accepted nan=False"], f"under -O: {lines}")


def test_the_refusal_reading_reds_on_the_O_run() -> None:
    """The mutation of 7.2 as a control: the refusal assertion applied to the ``-O`` run's output must red."""
    with pytest.raises(AssertionError):
        refusal_reads(_run("-O"))


def kept_directions_have_positive_sign(g_sym: np.ndarray, dim: int) -> bool:
    """True iff ``u_i . vt_i > 0`` on every kept direction of the hermitian SVD (the sign of the i-th eigenvalue)."""
    u, _, vt = np.linalg.svd(g_sym, hermitian=True)
    kept = g_sym.shape[0] - dim
    return bool(np.all(np.einsum("ji,ij->i", u, vt)[:kept] > 0))


@pytest.mark.parametrize("family", [f for f in sg.families() if sg.stored(f).frame_built])
@pytest.mark.rests_on("tests/gates/numerics/w1/test_dense_kernel_declared.py::test_declared_kernel_is_exact_for_the_frames")
def test_formal_path_has_no_non_positive_kept_direction(family: str) -> None:
    e = sg.stored(family)
    dim = sg.declared_for(e).dim
    _require(kept_directions_have_positive_sign(e.g_sym, dim), f"[{family}] a kept direction (all but the declared {dim}) has a non-positive eigenvalue")


def test_the_formal_path_population() -> None:
    frames = [f for f in sg.families() if sg.stored(f).frame_built]
    _require(len(frames) == 17, f"{len(frames)} frame-built Grams, expected 17")


@pytest.mark.rests_on("tests/gates/numerics/w1/test_dense_kernel_declared.py::test_declared_kernel_is_exact_for_the_frames")
def test_the_raw_path_check_is_not_a_rank_detector() -> None:
    """The documented limit: a kernel declared ONE DIMENSION SHORT on a deficient Gram is accepted by the sign test on at least one of
    the 14 deficient harness Grams (6 of 14 on the recording host, platform round-off), and the inverse face then reaches ``> 1e12``.
    The raw-path check is a misuse check, not a guarantee, so an accepted under-declaration is the limit, not a defect to be read as
    caught (ruling of 2026-10-30)."""
    accepted: list[tuple[str, float]] = []
    for family in sg.families():
        e = sg.stored(family)
        dim = sg.declared_for(e).dim
        if dim < 1 or not kept_directions_have_positive_sign(e.g_sym, dim - 1):
            continue
        u, s, vt = np.linalg.svd(e.g_sym, hermitian=True)
        kept = e.g_sym.shape[0] - (dim - 1)
        with np.errstate(all="ignore"):
            inv = vt.T @ (np.concatenate([1.0 / s[:kept], np.zeros(len(s) - kept)])[:, None] * u.T)
        accepted.append((family, float(np.max(np.abs(inv)))))
    _require(len(accepted) >= 1, "no under-declared deficient Gram was accepted by the sign test: the documented limit is not exercised")
    _require(all(m > 1e12 or not np.isfinite(m) for _, m in accepted), f"an accepted under-declaration with a small inverse: {accepted}")


def test_the_controls_read_non_positive() -> None:
    _require(not kept_directions_have_positive_sign(np.diag([1.0, -1.0]), 0), "diag([1, -1]) declared full rank must read non-positive")
    _require(not kept_directions_have_positive_sign(np.diag([1.0, 0.0]), 0), "diag([1, 0]) declared full rank must read non-positive")
