r"""#707 Step A acceptance: the STATE side does not move when the data side V* is stored integrated (R6).

Kind: RECORD, deliberately. Each row compares today's solve with the frozen pre-carve answer of the same solver
(``true_ends_acceptance.npz``, captured once by ``_capture_true_ends_acceptance.py`` at the commit in ``meta/commit``),
bit for bit. It certifies that the carve changed no k and no flux; it is NOT evidence that either is correct, which is
the business of the reference rows these problems already have (the regression snapshots' closed-form corroboration,
the carry rows, the MMS rows). R6's claim "k and the fluxes do not move" is exactly a bit-identity claim on an
unchanged solver, which is why ``np.array_equal`` and not a band.

A move is admitted only as a DECLARED event in ``true_ends_acceptance_rebaselines.json`` (read through
``expected()``), of one of two kinds, each gated here:

* ``gauge`` — the adjoint-eigen importance under the restated normalisation (census C8: the coordinate sum of
  ``F†ψ*`` becomes the pairing). The event states the factor ``c``; the row requires ``new ≈ c·old`` to
  ``10 × flux_tol`` (per-outer renormalisation changes the rounding, so the trajectory is not exactly proportional)
  and the adjoint ``k`` unmoved.
* a DECLARED BAND (``_true_ends_acceptance_cases.PRINCIPLED``): the rows Step A is predicted to move off
  bit-identity (multi-D DAG, LD, diffusion), each with its rtol and reason under the three criteria; bit-identity is
  still required everywhere else, and a banded row that stays bit-equal passes as bit-equal.
* ``principled`` — a non-bit move on any OTHER state row, with the three criteria of ``vv-principles`` ("Bit-identity vs
  principled-equivalence") written out and the measured ULP distance. Predicted candidates (spec row ACC-3): the
  multi-D DAG kernel (census C3 gains a ÷V on an integrated source), the LD DAG (the ``/mass_diag`` retirement), and
  diffusion k (``solver.py:431`` restated as a pairing).

On any red the message prints ``max|a−b|/max|b|`` before the nulp count (a nulp figure alone is uninterpretable).
"""

from __future__ import annotations

import subprocess
from functools import cache

import numpy as np
import pytest

from tests.gates.cross_method import _true_ends_acceptance_cases as acc

pytestmark = [pytest.mark.foundation, pytest.mark.filterwarnings("ignore::Warning")]

_CASES = {c.name: c for c in acc.cases()}


@cache
def _solved(name: str) -> dict[str, np.ndarray]:
    return _CASES[name].solve()


def _distance(a: np.ndarray, b: np.ndarray) -> str:
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    rel = float(np.max(np.abs(a - b)) / max(float(np.max(np.abs(b))), np.finfo(float).tiny))
    nulp = int(np.max(np.abs(a.view(np.int64) - b.view(np.int64)))) if a.shape == b.shape else -1
    return f"max|a-b|/max|b| = {rel:.3e}, max nulp = {nulp}, differing {int(np.sum(a != b))} of {a.size}"


def _event_for(key: str) -> dict | None:
    for event in reversed(acc.rebaseline_events()):
        if key in event["entries"]:
            return event
    return None


def test_the_fixture_is_the_pre_carve_anchor() -> None:
    """The capture holds every case's keys, no source value, and was taken at an ancestor of HEAD."""
    head = acc.head()
    keys = [k for k in head.files if not k.startswith("meta/")]
    assert len(_CASES) == 19 and len(keys) == 43, (len(_CASES), len(keys))
    for name in _CASES:
        assert any(k.startswith(f"{name}/") for k in keys), f"{name} has no frozen key"
    fields = {k.rsplit("/", 1)[1] for k in keys}
    assert fields <= {"k", "scalar_flux", "angular_flux"}, f"a non-state field is frozen: {fields}"
    # The capture commit must be an ancestor of HEAD: the anchor predates the carve. A squash or a rebase of the branch
    # REDS this leg on purpose (the recorded commit then names nothing in history): re-point ``meta/commit`` through a
    # declared event naming the rewrite, as the carry fixture's banner does, never by re-capturing.
    commit = str(head["meta/commit"])
    ancestor = subprocess.run(["git", "merge-base", "--is-ancestor", commit, "HEAD"], capture_output=True)
    assert ancestor.returncode == 0, f"the anchor {commit} is not an ancestor of HEAD"


@pytest.mark.parametrize("name", sorted(_CASES))
def test_state_is_unmoved(name: str) -> None:
    got = _solved(name)
    frozen = [k for k in acc.head().files if k.startswith(f"{name}/")]
    assert frozen, name
    for key in frozen:
        field = key.rsplit("/", 1)[1]
        event = _event_for(key)
        new = got[field]
        if event is not None and event["kind"] == "gauge":
            continue  # gated by test_declared_gauge_events
        want = acc.expected(key)
        assert new.shape == want.shape, f"{key}: shape {new.shape} != {want.shape}"
        declared = acc.band(name, field)
        if declared is None or np.array_equal(new, want):
            assert np.array_equal(new, want), f"{key} moved: {_distance(new, want)}"
            continue
        rtol, reason = declared
        rel = float(np.max(np.abs(new - want)) / max(float(np.max(np.abs(want))), np.finfo(float).tiny))
        assert rel <= rtol, f"{key} moved beyond its declared band {rtol:.1e} ({reason}): {_distance(new, want)}"


def test_declared_gauge_events() -> None:
    """Every ``gauge`` event: adjoint-eigen family only, k unmoved, ``new = c·old`` to 10 × flux_tol."""
    n = 0
    for event in acc.rebaseline_events():
        if event["kind"] != "gauge":
            continue
        for key, entry in event["entries"].items():
            name = key.rsplit("/", 1)[0]
            assert _CASES[name].family == "adjoint_eigen", f"{key}: a gauge event outside the adjoint-eigen family"
            assert not key.endswith("/k"), f"{key}: a gauge cannot move k"
            old = acc.head()[key]
            new = _solved(name)[key.rsplit("/", 1)[1]]
            c = float(entry["factor"])
            tol = 10.0 * acc.EIGEN["flux_tol"]
            rel = float(np.max(np.abs(new - c * old)) / np.max(np.abs(new)))
            assert rel <= tol, f"{key}: new != c*old (c = {c}): {rel:.3e} > {tol:.0e}"
            n += 1
    for event in acc.rebaseline_events():
        assert event["kind"] in {"gauge", "principled"}, event["kind"]
        if event["kind"] == "principled":
            for criterion in ("principled", "independent_reference", "fp_explanation", "max_ulp"):
                assert criterion in event, f"a principled event lacks {criterion!r}"
            for key, entry in event["entries"].items():
                # the carry pattern: the recorded old value IS the capture, and the key did move (else no event)
                np.testing.assert_array_equal(np.asarray(entry["old"], float), acc.head()[key],
                                              err_msg=f"{key}: the event's old value is not the capture")
                assert not np.array_equal(np.asarray(entry["old"], float), np.asarray(entry["new"], float)), key
    assert n >= 0


@pytest.mark.parametrize("body", ["slab_het_dd", "2d_het_dd"])
def test_adjoint_k_is_forward_k_in_the_capture(body: str) -> None:
    """REFERENCE between the two daggered solves, free from the capture: ``|k − k†| ≤ 10 × keff_tol`` (no solve)."""
    k, k_adj = float(acc.expected(f"sn/{body}/eigen/k")), float(acc.expected(f"sn/{body}/adjoint_eigen/k"))
    assert abs(k - k_adj) <= 10.0 * acc.EIGEN["keff_tol"], (body, k, k_adj)
