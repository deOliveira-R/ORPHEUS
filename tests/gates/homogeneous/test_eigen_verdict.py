r"""Homogeneous's eigen entry under the Fredholm verdict (U0): q is the fission emission chi.

Spec ``scratch/membership_covector/u0/spec_u0.md`` rung 4 (H1-H4). On a singular loss the eigen posing is solvable iff
no fission emission reaches an absorption-free closed class; the solve inverts the deflated loss ``A + Pi`` (the user's
ruling, 2026-10-31), so an exact zero pivot of ``A`` is never factorised. Before U0 (`[M]` 2026-10-31) both singular
fissile fixtures raised ``LinAlgError ... exactly singular (zero LU pivot)``, by pivot luck.
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.homogeneous.solver import HomogeneousProblem, solve_homogeneous_infinite
from orpheus.numerics.evidence import Measured
from orpheus.numerics.fredholm import Gauge, NoSolution, Solve

from tests.gates.transport._singular_loss_fixtures import (
    FIS_FEEDS, FIS_ISOLATED, FIS_ISOLATED_KINF, FIS_REGULAR, PURE_1G, dense_of,
)

_FRAGMENT = "no steady solution"
_EPS = np.finfo(float).eps

pytestmark = [pytest.mark.l1, pytest.mark.rests_on(
    "tests/gates/transport/test_absorbing_chain_predictor.py",
    "tests/gates/numerics/test_fredholm_verdict.py",
)]


def _dense_loss(mix):
    p = HomogeneousProblem(mix)
    n = int(np.prod(p.space.shape))
    return p, dense_of(lambda e: np.ravel(p.loss.apply(e.reshape(p.space.shape))), n)


@pytest.mark.catches("ERR-123")
def test_H1_emission_into_an_absorption_free_class_is_refused_by_the_verdict_not_the_pivot():
    """DISCRIMINATION: before U0 the same input raised LinAlgError('... zero LU pivot'). The left co-vector is (0.5, 1),
    so chi = (1, 0) pairs at 0.5 (the absorbing-chain probability of reaching the class)."""
    with pytest.raises(NoSolution, match=_FRAGMENT) as info:
        solve_homogeneous_infinite(FIS_FEEDS)
    assert "zero LU pivot" not in str(info.value)


def test_H2_an_isolated_absorption_free_class_is_gauged():
    """Gauge: k_inf = 1.25 (closed form nu Sf1 / (St1 - Ss11)); the loss is diag(0.1, 0), an EXACT zero pivot, so the
    member comes from (A + Pi)^-1 and its class component is zero (the minimum-norm member) to round-off."""
    r = solve_homogeneous_infinite(FIS_ISOLATED)
    assert isinstance(r.solvability, Gauge)
    assert r.k_inf == pytest.approx(FIS_ISOLATED_KINF, rel=1e-13)
    phi = np.ravel(np.asarray(r.flux, dtype=float))
    assert abs(phi[1]) <= 64 * _EPS * np.abs(phi).max()


def test_H3_a_regular_medium_is_solved_and_carries_its_rcond():
    """k_inf unchanged (`[M]` 0.625 before and after); the result carries the Solve verdict with the LU's rcond, at least
    1/cond_1 and at most 10/cond_1 (`[M]` 0.067)."""
    r = solve_homogeneous_infinite(FIS_REGULAR)
    assert r.k_inf == pytest.approx(0.625, rel=1e-13)
    verdict = r.solvability
    assert isinstance(verdict, Solve) and isinstance(verdict.conditioning, Measured)
    _, A = _dense_loss(FIS_REGULAR)
    true = 1.0 / np.linalg.cond(A, 1)
    assert true * (1 - 64 * _EPS) <= verdict.conditioning.value <= 10 * true


def test_H4_a_non_fissile_singular_medium_is_still_refused():
    """Not U0's question (ruled 2026-10-31): an eigen posing with F = 0 has no eigenpair, whatever A is; the refusal stays
    (the ScaleGauge's, at any message). The LU of A + Pi before it is legitimate (the deflated loss is nonsingular)."""
    with pytest.raises(ValueError):
        solve_homogeneous_infinite(PURE_1G)
