r"""#616: the transposed sweep :math:`(L+C)^{-\mathsf T}` on the multi-D loss representations.

Spec ``scratch/membership_covector/u3b/gate_spec.md`` §1, rows T1-T4, T6, T7 (T5 merged into T2-T4, 2026-10-31). The 1-D reverse scan (``CumprodScan``, #280) is
the foundation; each multi-D representation (``ScanMarch``, ``MovingFrontierWindow``, ``FullFieldWavefront``) is FORCED
on every fixture it supports, with an activation leg asserting the forced class is the one the operator holds.

Every reference is the dense :math:`M = L + C` assembled column by column through the operator's own ``apply`` (the
forward matvec, which #616 does not touch), then ``np.linalg.solve(M.T, b)``: independent of the transposed walk.

Tolerances: the FORWARD one-sided identity reads 0.1-0.2 eps, normalised by ``‖M‖₁‖x‖∞``, on these fixtures with
``cond₁(M) <= 441`` (`[M]` 2026-10-31, ``scratch/membership_covector/u3b/ta/p8_616_baseline.py``); a reverse
substitution has the same depth, so the transposed identity is held to 64 eps, and the agreement with the dense
transposed solve to ``10 cond₁(M) eps``.

Before #616 landed, every multi-D row was red (``has_transpose_sweep`` False, ``(L+C).H.inverse()`` raising
``NotInvertible``) and the 1-D rows were the positive control of the same assertions. ``[M]`` 2026-10-31, on the landed
carve: 55 of 55 green; dropping the written-row zeroing of the trace contract reddens 39 of 55, and dropping the DD
kernel's coupling term reddens 33 of 55.
"""
from __future__ import annotations

from contextlib import contextmanager
from functools import lru_cache
from unittest.mock import patch

import numpy as np
import pytest

from orpheus.numerics.quadrature import Quadrature
from tests.gates.sn.operators import test_g_adjoint_reciprocity as G
from tests.gates.sn.operators import test_loss_kernel_gauge as T
from tests.gates.transport._singular_loss_fixtures import dense_of

pytestmark = [
    pytest.mark.l0,
    pytest.mark.rests_on(
        "tests/gates/sn/solve/test_scan_march_end_to_end.py",
        "tests/gates/sn/operators/test_g_adjoint_reciprocity.py",
    ),
]

_EPS = np.finfo(float).eps


def _ld():
    from orpheus.transport.spatial.linear_discontinuous import LinearDiscontinuous
    return LinearDiscontinuous()


def _first(built):
    return built[0] if isinstance(built, tuple) else built


#: fixture id -> builder of an SNProblem.
_FIXTURES = {
    "slab_2g": lambda: _first(G._make_slab(ng=2)),
    "cart2d_2g": lambda: _first(G._make_cart2d(ng=2)),
    "d2_vacuum": lambda: T._mesh((3, 4), [(T._V, T._V)] * 2, c=0.5, ng=2),
    "d2_graded_mirror": lambda: T._mesh((3, 4), [(T._R, T._R)] * 2, c=0.9, ng=2, stretch=1.5),
    "ld_2d_2g": lambda: _first(G._make_ld_2d(ng=2)),
    "d3_mirror": lambda: T._mesh((2, 3, 2), [(T._R, T._R)] * 3, c=0.9, ng=1),
    "cart2d_product44_mirror": lambda: T._mesh((3, 4), [(T._R, T._R)] * 2, c=0.9, ng=1,
                                               quad=Quadrature.product(n_mu=4, n_phi=4)),
    #: Lebedev order 5 on the 2-D box: 2 pure-z ordinates (no in-plane streaming: the cell solve is sigma psi = q) and
    #: 4 tangential ones (zero trace metric on one axis's faces); 56 of 364 dofs in ker G (`[M]` qa p1).
    "lebedev5_mirror": lambda: T._mesh((3, 4), [(T._R, T._R)] * 2, c=0.9, ng=1, quad=Quadrature.lebedev(order=5)),
    "ld_product44_mirror": lambda: T._mesh((3, 4), [(T._R, T._R)] * 2, c=0.9, ng=1,
                                           quad=Quadrature.product(n_mu=4, n_phi=4), scheme=_ld()),
}

#: (fixture, representation) cells: every multi-D representation on every fixture it supports (`[M]` p8 supports()).
_CELLS = [
    ("slab_2g", "CumprodScan"),
    ("cart2d_2g", "ScanMarch"), ("cart2d_2g", "MovingFrontierWindow"),
    ("cart2d_2g", "FullFieldWavefront"),
    ("d2_vacuum", "ScanMarch"), ("d2_vacuum", "MovingFrontierWindow"),
    ("d2_vacuum", "FullFieldWavefront"),
    ("d2_graded_mirror", "ScanMarch"), ("d2_graded_mirror", "MovingFrontierWindow"),
    ("d2_graded_mirror", "FullFieldWavefront"),
    ("ld_2d_2g", "MovingFrontierWindow"), ("ld_2d_2g", "FullFieldWavefront"),
    ("d3_mirror", "FullFieldWavefront"),
    ("cart2d_product44_mirror", "ScanMarch"), ("cart2d_product44_mirror", "MovingFrontierWindow"),
    ("cart2d_product44_mirror", "FullFieldWavefront"),
    ("lebedev5_mirror", "ScanMarch"), ("lebedev5_mirror", "MovingFrontierWindow"),
    ("lebedev5_mirror", "FullFieldWavefront"),
]

#: The singular-metric fixtures, for the full-space reciprocity row T6 (ker G slots included).
_SINGULAR = [("cart2d_product44_mirror", "MovingFrontierWindow"), ("lebedev5_mirror", "ScanMarch"),
             ("lebedev5_mirror", "MovingFrontierWindow"), ("ld_product44_mirror", "FullFieldWavefront")]


def _params(cells):
    return [pytest.param(f, r, id=f"{f}-{r}") for f, r in cells]


@contextmanager
def _forced(name: str):
    """Make ``default_for`` return representation ``name`` (the ``window_forced`` idiom of the end-to-end gates)."""
    from orpheus.sn import loss_representation as lr

    cls = getattr(lr, name)

    def forced(problem, spatial_closure, angular_closure):
        if not cls.supports(problem, spatial_closure).ok:
            raise AssertionError(f"{name} does not support this fixture: the cell is mis-tabled")
        return cls(problem, spatial_closure, angular_closure)

    with patch.object(lr, "default_for", forced):
        yield cls


@lru_cache(maxsize=None)
def _case(fixture: str, representation: str):
    """``(LC, template, M dense, metric diagonal)`` with ``representation`` forced and its activation checked."""
    with _forced(representation) as cls:
        problem = _FIXTURES[fixture]()
        system, template = T._loss_system(problem)
        lc = system.factors.streaming_collision
        rep = lc.loss_representation
    if type(rep) is not cls:
        raise AssertionError(f"activation: the operator holds {type(rep).__name__}, not the forced {representation}")
    n = template.to_flat().size
    M = dense_of(lambda e: lc.apply(type(template).from_flat(e, template)).to_flat(), n)
    g = np.asarray(problem.full_field_space.metric_diagonal(), float).ravel()
    return lc, template, M, g


def _quotient(M, g):
    r"""The G > 0 coordinates and :math:`M_{pp}`. On a singular trace metric the ker-G slots carry identically zero rows
    and columns of :math:`M` (`[M]` rank deficiency = dim ker G on product(4,4) and Lebedev-5, qa p1), so :math:`M_{pp}`
    is :math:`M` on the Hilbert space the question is posed on; on a full metric it is :math:`M` itself."""
    pos = np.flatnonzero(g > 0.0)
    Mpp = M[np.ix_(pos, pos)]
    zero = np.flatnonzero(g == 0.0)
    if zero.size and (np.abs(M[zero]).max() > 0.0 or np.abs(M[:, zero]).max() > 0.0):
        raise AssertionError("premise: M is not zero on the ker-G slots; the quotient form does not apply")
    return pos, Mpp


def _flat(op_apply, template, x):
    return op_apply(type(template).from_flat(x, template)).to_flat()


def _source(n: int) -> np.ndarray:
    return np.random.default_rng(7).standard_normal(n)


# ── T1: the flags ──────────────────────────────────────────────────────────────────────────────────────────

@pytest.mark.parametrize("fixture,representation", _params(_CELLS))
def test_T1_the_representation_declares_its_transposed_sweep(fixture, representation):
    """THEOREM of the type: ``has_transpose_sweep`` and ``is_transpose_invertible``. First red: today's False."""
    lc, *_ = _case(fixture, representation)
    assert lc.loss_representation.has_transpose_sweep is True
    assert lc.is_transpose_invertible is True


# ── T2/T3: the one-sided transposed identity, and the dense transposed solve ────────────────────────────────

@pytest.mark.parametrize("fixture,representation", _params(_CELLS))
def test_T2_the_transposed_sweep_inverts_the_transpose(fixture, representation):
    r"""One-sided: :math:`\lVert M^{\mathsf T}x - b\rVert_\infty / (\lVert M\rVert_1\lVert x\rVert_\infty) \le 64\varepsilon`
    for :math:`x = (L+C)^{-\mathsf T}b` on the G > 0 block (all of it on a full metric), b random over bulk AND trace rows (a reciprocity row is blind to a symmetric
    regression, ``vv-principles`` #37). Arms: ``forward_schedule``, ``drop_outflow_rows``, ``skip_trace_rows``."""
    lc, template, M, g = _case(fixture, representation)
    pos, Mpp = _quotient(M, g)
    b = np.zeros(M.shape[0])
    b[pos] = _source(pos.size)
    x = _flat(lc.solve_transpose, template, b)[pos]
    residual = np.abs(Mpp.T @ x - b[pos]).max() / (np.abs(Mpp).sum(axis=0).max() * np.abs(x).max())
    assert residual <= 64 * _EPS, residual / _EPS


@pytest.mark.parametrize("fixture,representation", _params(_CELLS))
def test_T3_the_transposed_sweep_equals_the_dense_transposed_solve(fixture, representation):
    """REFERENCE: ``np.linalg.solve(M.T, b)`` within ``10 cond₁(M) eps`` relative, through ``(L+C).inverse()``'s
    transpose (the operator door, a second spelling beside T2's verb). On Lebedev-5 it is the only row that reaches the
    pure-z branch (no in-plane streaming): arm ``pure_z_times`` (multiply by sigma instead of divide) reds it
    (`[M]` qa p2: 1.7e-2 on the G > 0 block)."""
    lc, template, M, g = _case(fixture, representation)
    pos, Mpp = _quotient(M, g)
    b = np.zeros(M.shape[0])
    b[pos] = _source(pos.size)
    x = _flat(lc.inverse().apply_transpose, template, b)[pos]
    exact = np.linalg.solve(Mpp.T, b[pos])
    error = np.abs(x - exact).max() / np.abs(exact).max()
    assert error <= 10 * np.linalg.cond(Mpp, 1) * _EPS, error


# ── T4: the metric adjoint ─────────────────────────────────────────────────────────────────────────────────

@pytest.mark.parametrize("fixture,representation", _params(_CELLS))
def test_T4_the_adjoint_inverse_is_the_metric_transpose_and_reciprocal(fixture, representation):
    r""":math:`((L+C)^{-1})^{\mathsf H} = G^{-1}(L+C)^{-\mathsf T}G` on the G > 0 rows (64 eps relative), and
    :math:`\langle (L+C)^{-1}p, q\rangle_G = \langle p, (L+C)^{\mathsf H\,-1}q\rangle_G`. The grid holds no uniform
    square fixture, where a Euclidean ``.H`` would coincide (arm ``euclidean_H``)."""
    lc, template, M, g = _case(fixture, representation)
    n = M.shape[0]
    pos, Mpp = _quotient(M, g)
    kappa = np.linalg.cond(Mpp, 1)
    rng = np.random.default_rng(11)
    p, q = rng.standard_normal(n), rng.standard_normal(n)
    p[g == 0.0] = 0.0
    q[g == 0.0] = 0.0
    adjoint_inverse = lc.inverse().H   # the swap law's door: (M^-1)^H = (M^H)^-1
    y = _flat(adjoint_inverse.apply, template, q)
    expected = np.linalg.solve(Mpp.T, (g * q)[pos]) / g[pos]
    np.testing.assert_allclose(y[pos], expected, rtol=0, atol=64 * _EPS * kappa * np.abs(expected).max())
    forward = _flat(lc.inverse().apply, template, p)
    lhs, rhs = float(np.dot(g * forward, q)), float(np.dot(g * p, y))
    assert abs(lhs - rhs) <= 64 * _EPS * kappa * (abs(lhs) + abs(rhs)), (lhs, rhs)


# ── T6: the trace contract on the ker-G slots ──────────────────────────────────────────────────────────────

@pytest.mark.parametrize("fixture,representation", _params(_SINGULAR))
def test_T6_the_transposed_sweep_is_the_transpose_of_the_sweep_on_the_full_space(fixture, representation):
    r"""Euclidean reciprocity over the WHOLE coordinate space, ker-G slots included:
    :math:`|\langle S b, c\rangle - \langle b, S_{\mathsf T} c\rangle| \le 64\varepsilon\,\lVert Sb\rVert\lVert c\rVert`
    for :math:`S` the forward sweep and :math:`S_{\mathsf T}` the transposed one, 12 random pairs. The forward sweep
    WRITES into the tangential slots (`[M]` p1b: max 1.75), so the transposed one must read them as their transpose; no
    metric enters, so the rows T2-T4 (posed on the G > 0 block) are blind to this clause by construction. `[M]` qa p2:
    3.8e-17 to 5.1e-17 honest; arms ``tangential_wrong_face`` 2.1e-2 to 2.9e-2 and ``tangential_not_zeroed`` 1.7e-2 to
    2.1e-2 (the G > 0 block reads 1e-17 under both: the clause lives only on ker G)."""
    lc, template, M, g = _case(fixture, representation)
    assert np.any(g == 0.0), "the fixture lost its ker-G slots: the row would test nothing"
    rng = np.random.default_rng(5)
    worst = 0.0
    for _ in range(12):
        b, c = rng.standard_normal(g.size), rng.standard_normal(g.size)
        Sb = _flat(lc.solve, template, b)
        STc = _flat(lc.solve_transpose, template, c)
        worst = max(worst, abs(float(Sb @ c) - float(b @ STc)) / (np.linalg.norm(Sb) * np.linalg.norm(c)))
    assert worst <= 64 * _EPS, worst / _EPS


# ── T7: the two wavefront schedules agree bitwise ──────────────────────────────────────────────────────────

@pytest.mark.parametrize("fixture", ["cart2d_2g", "ld_2d_2g", "lebedev5_mirror"])
@pytest.mark.parametrize("seed", range(5))
def test_T7_the_window_and_the_full_field_wavefront_transpose_agree_bitwise(fixture, seed):
    """``array_equal`` over 5 seeds (`[M]` qa p2: identical at seed 9 on all three). The two classes share the octant
    walk's transposed arithmetic and differ in storage (a rolling window against the full field), so this row pins the
    window's bookkeeping, not the arithmetic: a disagreement means the window dropped or reordered a cell. It is not
    an independent reference for the transpose (T2-T4 are)."""
    window, template_w, *_ = _case(fixture, "MovingFrontierWindow")
    full, template_f, *_ = _case(fixture, "FullFieldWavefront")
    c = np.random.default_rng(seed).standard_normal(template_w.to_flat().size)
    np.testing.assert_array_equal(_flat(window.solve_transpose, template_w, c), _flat(full.solve_transpose, template_f, c))
