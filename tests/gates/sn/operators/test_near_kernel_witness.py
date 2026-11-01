r"""U3b: the near-kernel witness against SN's predicted kernel and the dense SVD of the production loss.

Spec ``scratch/membership_covector/u3b/gate_spec.md`` §2.3-2.4 (rows S-W1..S-W3, B-W1..B-W3). The witness
(:mod:`orpheus.numerics.near_kernel`) never decides in production; these rows check that what SN PREDICTS from declared
structure (:attr:`~orpheus.sn.problem.SNProblem.predicted_kernel`) is what the operator HAS. Two oracles, independent of
each other: the dense SVD of the loss assembled through the production matvec (``_dense``), and the witness's own routes.

The carrier-valued splitting reaches the witness through ``FlattenedOperator`` (the witness acts on flat coordinates).
"""
from __future__ import annotations

from functools import reduce
from operator import add

import numpy as np
import pytest

from orpheus.numerics.flat_operator import FlattenedOperator
from orpheus.numerics.near_kernel import NearKernel
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.splitting import Splitting, resolve_schedule
from tests.gates.sn.operators import test_loss_kernel_gauge as T
from tests.gates.numerics._near_kernel_oracles import dense_kernels, derived_containment_tol, residual
from tests.gates.sn.operators.test_predicted_loss_kernel import FEEDER, _GL8, _R, _columns, _dense, _problem

pytestmark = [pytest.mark.l1, pytest.mark.rests_on(
    "tests/gates/numerics/test_near_kernel.py",
    "tests/gates/sn/operators/test_transpose_sweep_multid.py",
)]

_LS2 = Quadrature.level_symmetric(sn_order=2)


def _parts(problem, schedule="jacobi"):
    r"""The production splitting :math:`A = M - N` of the problem's loss, flattened onto one chart."""
    system, template = T._loss_system(problem)
    splitting = Splitting.from_schedule(system, resolve_schedule(problem, schedule))
    return (FlattenedOperator(splitting.implicit, template),
            FlattenedOperator(reduce(add, splitting.explicit), template))


def _dense_witness(problem):
    implicit, lagged = _parts(problem)
    return NearKernel.dense(implicit - lagged)


def _block_witness(problem, schedule="jacobi"):
    return NearKernel.of_splitting(*_parts(problem, schedule))


def _metric(problem):
    return np.asarray(problem.full_field_space.metric_diagonal(), dtype=float).ravel()


# ── S-W1: the union counts and spans, dense route ───────────────────────────────────────────────────────────

_UNION = [("LS2 ng=1 c=0.9", dict(ng=1, quad=_LS2, c=0.9), 2), ("LS2 ng=1 c=1", dict(ng=1, quad=_LS2, c=1.0), 3),
          ("LS4 ng=2 c=0.9", dict(ng=2, c=0.9), 12), ("LS4 ng=2 c=1", dict(ng=2, c=1.0), 14)]


@pytest.mark.parametrize("label,kw,dim", _UNION, ids=[r[0] for r in _UNION])
def test_SW1_the_witness_counts_the_union_and_the_predicted_spans_contain_it(label, kw, dim):
    """Measured dim == dense dim == predicted dim; the predicted right span equals the measured one; the predicted left
    co-vectors' Riesz representatives (``kernel.H.right``) contain the measured left vectors; gap >= 1e10. `[M]` the
    test-architect's p2: gap >= 2.29e12, left containment <= 7.1e-14."""
    problem = T._mesh((3, 4), [(T._R, T._R)] * 2, **kw)
    A = _dense(problem)
    sv, R, LH, k = dense_kernels(A, _metric(problem))
    w = _dense_witness(problem)
    kernel = problem.predicted_kernel
    n = A.shape[0]
    assert w.dimension == k == dim == kernel.dimension, (w.dimension, k, kernel.dimension)
    predicted_right = _columns(kernel.right, n)
    assert residual(predicted_right, w.right) <= 1e-12 and residual(w.right, predicted_right) <= 1e-12
    assert residual(_columns(kernel.H.right, n), w.left) <= 1e-12
    assert w.singular_gap is not None and w.singular_gap >= 1e10


# ── S-W2: the witness is blind to the prediction ────────────────────────────────────────────────────────────

def test_SW2_the_witness_is_predictor_blind(monkeypatch):
    """A 12-only predictor (the physical half made trivial) on LS4 ng=2 c=1: both witness routes still measure 14, and the
    comparison S-W1 asserts would red on this input."""
    import sys
    from orpheus.numerics.fredholm import PredictedKernel
    import orpheus.sn.operators.loss_kernel_gauge as M

    def trivial(problem):
        return PredictedKernel.trivial(problem.full_field_space, "mutated: no physical kernel")

    rebound = 0
    for module in list(sys.modules.values()):
        if getattr(module, "predicted_physical_kernel", None) is M.predicted_physical_kernel:
            monkeypatch.setattr(module, "predicted_physical_kernel", trivial)
            rebound += 1
    assert rebound >= 1
    problem = T._mesh((3, 4), [(T._R, T._R)] * 2, ng=2, c=1.0)
    assert problem.predicted_kernel.dimension == 12          # activation: the predictor IS mutated
    for w in (_dense_witness(problem), _block_witness(problem)):
        assert w.dimension == 14
        assert w.dimension != problem.predicted_kernel.dimension


# ── S-W3: the tangential slots are quotiented ───────────────────────────────────────────────────────────────

@pytest.mark.parametrize("c,expected", [(0.9, 0), (1.0, 1)])
def test_SW3_the_tangential_slots_are_quotiented(c, expected):
    """product(4,4): the raw nullity is at least 112 (the tangential trace slots, ker G), modulo ker G it is 0 (c = 0.9)
    or 1 (c = 1). `[M]` the test-architect's p1/p1b. At c = 1 the predictor refuses today (the scope door's tangential
    entry), so the witness measures alone; the row is its quotient's only gate on a production loss."""
    problem = T._mesh((3, 4), [(T._R, T._R)] * 2, c=c, ng=1, quad=Quadrature.product(n_mu=4, n_phi=4))
    A = _dense(problem)
    g = _metric(problem)
    sv, R, LH, k = dense_kernels(A, g)
    assert k == expected and int(np.sum(sv < 1e-10 * sv[0])) >= 112
    for w in (_dense_witness(problem), _block_witness(problem)):
        assert w.dimension == expected
        if expected:
            on_metric = g > 0
            assert residual(R[on_metric], w.right[on_metric]) <= 1e-10


# ── B-W: the block route against the dense route ────────────────────────────────────────────────────────────

def _thick1d():
    from orpheus.mesh import AxisMesh
    from tests.gates._exterior import conditions_for
    axes = (AxisMesh(edges=np.linspace(0, 32.0, 65), role_low=T._R, role_high=T._R),)
    return T._as_problem(axes, Quadrature.gauss_legendre(n_ordinates=16), {0: T._mixture(1, 1.0)}, scheme=None,
                         boundary_conditions=conditions_for(axes, T._CONDITIONS))


def _thick2d():
    from orpheus.mesh import AxisMesh
    from tests.gates._exterior import conditions_for
    axes = tuple(AxisMesh(edges=T._graded(4.0, 8), role_low=T._R, role_high=T._R) for _ in range(2))
    return T._as_problem(axes, Quadrature.level_symmetric(sn_order=4), {0: T._mixture(1, 1.0)}, scheme=None,
                         boundary_conditions=conditions_for(axes, T._CONDITIONS))


from orpheus.numerics.operator import MissingAdjoint

_GS = pytest.mark.xfail(strict=True, raises=MissingAdjoint,
                        reason="#681: the Gauss-Seidel splitting's transposed solve, so its left route")

#: (id, builder, schedule, k, the averaged map's gap `[M]` the test-architect's p7/p9 or None, marks)
_BLOCK = [
    ("ls4c09-jacobi", lambda: T._mesh((3, 4), [(T._R, T._R)] * 2, c=0.9, ng=2), "jacobi", 12, 0.02543, ()),
    ("ls4c09-gs", lambda: T._mesh((3, 4), [(T._R, T._R)] * 2, c=0.9, ng=2), "gauss_seidel", 12, 0.03199, (_GS,)),
    ("ls4c1-jacobi", lambda: T._mesh((3, 4), [(T._R, T._R)] * 2, c=1.0, ng=2), "jacobi", 14, 0.07729, ()),
    ("ls4c1-gs", lambda: T._mesh((3, 4), [(T._R, T._R)] * 2, c=1.0, ng=2), "gauss_seidel", 14, 0.08726, (_GS,)),
    ("feeder1d-jacobi", lambda: _problem((6,), [(_R, _R)], _GL8, {0: FEEDER}), "jacobi", 1, 0.2019, ()),
    ("thick1d-jacobi", _thick1d, "jacobi", 1, 0.002395, (pytest.mark.slow,)),
    ("thick2d-3.2mfp-jacobi", _thick2d, "jacobi", 7, None, (pytest.mark.slow,)),
]


@pytest.mark.parametrize("label,build,schedule,k,gap",
                         [pytest.param(*r[:5], marks=r[5], id=r[0]) for r in _BLOCK])
def test_BW_the_block_route_agrees_with_the_dense_route(label, build, schedule, k, gap):
    """Right AND left spans (the left through the dual split's map; d >= 2 rests on the multi-D transposed sweep, #616)
    within the Davis-Kahan tolerance derived from the route's own stopping residual; each right vector is null to that
    residual; the averaged map's gap equals the dense one to rtol 1e-3."""
    problem = build()
    A = _dense(problem)
    sv, R, LH, kd = dense_kernels(A, _metric(problem))
    w = _block_witness(problem, schedule)
    assert w.dimension == kd == k
    tol = derived_containment_tol(w.stop_residual, sv, k)
    assert residual(R, w.right) <= tol and residual(LH, w.left) <= tol
    for z in w.right.T:
        assert np.linalg.norm(A @ z) <= (w.stop_residual * (1 + 1e-6) + 64 * np.finfo(float).eps) * sv[0]
    if gap is not None:
        assert w.splitting_gap is not None
        np.testing.assert_allclose(w.splitting_gap, gap, rtol=1e-3)


# ── the flattened sweep the block route reads ───────────────────────────────────────────────────────────────

def test_the_flattened_sweep_inverts_the_flattened_loss_and_transposes():
    r"""``FlattenedOperator(L + C).inverse()`` is the flattened sweep (the chart is an isometry, so
    :math:`(\flat A\flat^{-1})^{-1} = \flat A^{-1}\flat^{-1}`): it inverts the flattened :math:`L + C` to round-off, and
    its dual arrow is its Euclidean transpose (the multi-D transposed sweep, #616). The block route of the witness
    reads exactly these two maps; this row pins them directly, on the d = 2 mirror box."""
    problem = T._mesh((3, 4), [(T._R, T._R)] * 2, c=0.9, ng=2)
    implicit, _ = _parts(problem)
    assert implicit.is_invertible
    sweep = implicit.inverse()
    rng = np.random.default_rng(616)
    n = implicit.n_flat
    x, y = rng.standard_normal(n), rng.standard_normal(n)
    back = np.asarray(implicit.apply(sweep.apply(x)))
    assert np.linalg.norm(back - x) <= 1e-12 * np.linalg.norm(x), np.linalg.norm(back - x) / np.linalg.norm(x)
    lhs = float(np.dot(np.asarray(sweep.apply(x)), y))
    rhs = float(np.dot(x, np.asarray(sweep.dual().apply(y))))
    assert abs(lhs - rhs) <= 1e-12 * (abs(lhs) + abs(rhs)), (lhs, rhs)
