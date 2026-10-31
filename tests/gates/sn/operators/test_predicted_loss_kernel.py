r"""SN's predicted loss kernel (U0): :math:`\ker A` and :math:`\ker A^{\mathsf T}` from structure, against a dense SVD.

:attr:`~orpheus.sn.problem.SNProblem.predicted_kernel` is the sum of two predictors on the full field (bulk then
trace, flat):

* :func:`~orpheus.sn.operators.loss_kernel_gauge.predicted_discretization_kernel`, the diamond-difference face modes:
  right generators :math:`\Phi` on the trace, left :math:`W\Phi` with :math:`W = -w_n(\Omega_n\cdot\hat n_f)A_f`
  (the numerics-investigator's derivation, ``scratch/membership_covector/u0/left_basis/MEMO.md``);
* :func:`~orpheus.sn.operators.loss_kernel_gauge.predicted_physical_kernel`, one mode per absorption-free closed
  energy class: right the class's energy null vector :math:`e` on every bulk and trace slot, left :math:`u\,W` with
  :math:`u` the energy left vector (1 on the class, the absorbing-chain probability on its feeders), bulk weight
  :math:`w_n V_{\rm cell}` (``scratch/membership_covector/u0/physical/MEMO.md``).

Every row's reference is the dense SVD of the production loss assembled column by column (``_apply_loss``), which
shares nothing with the predictors; the left weights of S16b are rebuilt here from the quadrature, the trace layout
and the cell edges, never from the module's ``_balance_weights`` (X4). Spec
``scratch/membership_covector/u0/spec_u0.md`` rows S8, S16a, S16b and the physical-predictor rows.

Teeth (`[M]` 2026-10-31, ``scratch/membership_covector/u0/battery/``, a ``-p`` plugin per arm under ``python -O``):
``left_is_right`` (the balance weights replaced by ones), ``trace_weights_zeroed`` (their trace part zeroed) and
``e_for_u`` (the energy right vector handed as the left one); the arm table is in the battery's log.
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.data.macro_xs.synthetic import make_mixture
from orpheus.geometry.faces import face_normal
from orpheus.mesh import AxisMesh
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.operators.loss_kernel_gauge import (
    predicted_discretization_kernel,
    predicted_kernel_dimension,
    predicted_physical_kernel,
)
from orpheus.sn.solver import _as_problem
from orpheus.transport.boundary_conditions import WallCondition
from tests.gates.sn.operators.test_loss_kernel_gauge import _R, _apply_loss, _loss_system, _mesh

pytestmark = [
    pytest.mark.l0,
    pytest.mark.rests_on(
        "tests/gates/sn/operators/test_loss_kernel_gauge.py::test_EVERY_basis_vector_is_annihilated_by_the_production_matvec",
        "tests/gates/sn/operators/test_loss_kernel_gauge.py::test_the_dimension_matches_a_DENSE_SVD_of_the_assembled_operator",
        "tests/gates/transport/test_absorbing_chain_predictor.py",
    ),
]

_SPAN_TOL = 1e-10


# ── the oracle ───────────────────────────────────────────────────────────────────────────────────────────

def _dense(problem) -> np.ndarray:
    system, template = _loss_system(problem)
    n = template.to_flat().size
    A = np.empty((n, n))
    e = np.zeros(n)
    for j in range(n):
        e[j] = 1.0
        A[:, j] = _apply_loss(system, template, e)
        e[j] = 0.0
    return A


def _dense_kernels(A):
    U, s, Vt = np.linalg.svd(A)
    k = int(np.sum(s < 1e-10 * s[0]))
    n = A.shape[0]
    return U[:, n - k:], Vt[n - k:].T


def _columns(summands, n) -> np.ndarray:
    """Summands embedded as full-field columns, by hand (never the production ``_embedded_span``)."""
    cols = []
    for s in summands:
        g = np.asarray(s.generators, dtype=float)
        for j in range(g.shape[1]):
            v = np.zeros(n)
            v[np.asarray(s.support)] = g[:, j]
            cols.append(v)
    return np.column_stack(cols) if cols else np.zeros((n, 0))


def _residual(span, vecs) -> float:
    """``|vecs - P_span vecs| / |vecs|``, P the Euclidean projector onto span's numerical range."""
    if vecs.shape[1] == 0:
        return 0.0
    if span.shape[1] == 0:
        return 1.0
    Q, sv, _ = np.linalg.svd(span, full_matrices=False)
    Q = Q[:, sv > 1e-10 * sv[0]]
    return float(np.linalg.norm(vecs - Q @ (Q.T @ vecs)) / np.linalg.norm(vecs))


def _assert_spans_agree(label, predicted, measured):
    a, b = _residual(predicted, measured), _residual(measured, predicted)
    assert a <= _SPAN_TOL and b <= _SPAN_TOL, f"{label}: measured in predicted {a:.2e}, predicted in measured {b:.2e}"


# ── fixtures ─────────────────────────────────────────────────────────────────────────────────────────────

def _mix(sig_t, sig_s, sig_2=None):
    sig_t = np.asarray(sig_t, float)
    sig_s = np.asarray(sig_s, float)
    ng = sig_t.size
    s2 = np.zeros((ng, ng)) if sig_2 is None else np.asarray(sig_2, float)
    return make_mixture(sig_t=sig_t, sig_c=sig_t - sig_s.sum(axis=1) - s2.sum(axis=1), sig_f=np.zeros(ng),
                        nu=np.zeros(ng), chi=np.zeros(ng), sig_s=sig_s, sig_2=s2)


#: Every value is DYADIC, so the class's declared removal is exactly zero in floats (the predicate reads it with
#: ``math.fsum``; the investigator's 1.2 / 0.5 / 0.7 nets to round-off, which the data-as-stored rule reads as REGULAR,
#: spec section 7, #665). Class {1, 2} absorption-free; group 0 an absorbing FEEDER (absorption 0.125, u_0 = 0.625/0.75).
FEEDER = _mix([1.0, 1.25, 1.5], [[0.25, 0.5, 0.125], [0.0, 0.5, 0.75], [0.0, 0.5, 1.0]])
#: (n,2n) inside the class (1 -> 2, 0.125, counted twice by the loss), s12 lowered by 0.25 to keep the class balanced.
N2N = _mix([1.0, 1.25, 1.5], [[0.25, 0.5, 0.125], [0.0, 0.5, 0.5], [0.0, 0.5, 1.0]],
           sig_2=[[0, 0, 0], [0, 0, 0.125], [0, 0, 0]])
#: A density-scaled copy of FEEDER (x 2, exact; the same e and u): a two-region class the closed forms still cover.
FEEDER_X2 = _mix(2.0 * np.array([1.0, 1.25, 1.5]), 2.0 * np.array([[0.25, 0.5, 0.125], [0.0, 0.5, 0.75], [0.0, 0.5, 1.0]]))
#: A different class block (e proportional to (1, 0.75) against FEEDER's (2/3, 1)): the closed form's scope boundary.
OTHER_CLASS = _mix([1.5, 1.0, 2.0], [[0.0, 0.75, 0.5], [0.0, 0.25, 0.75], [0.0, 1.0, 1.0]])


def test_premise_the_classes_are_exactly_balanced_in_floats():
    import math
    for m in (FEEDER, N2N, FEEDER_X2, OTHER_CLASS):
        S, S2 = m.SigS[0].toarray(), m.Sig2[0].toarray()
        for g in (1, 2):
            assert math.fsum([m.SigT[g], *(-S[g]), *(-2.0 * S2[g])]) == 0.0


def _graded(extent, n, stretch):
    w = stretch ** np.arange(n)
    return np.concatenate([[0.0], np.cumsum(w / w.sum() * extent)])


def _problem(cells, roles, quad, mats, *, conds=None, mat_map=None, stretch=1.0):
    axes = tuple(AxisMesh(edges=_graded(e, n, stretch), role_low=lo, role_high=hi)
                 for e, n, (lo, hi) in zip((1.0, 2.0, 3.0), cells, roles))
    return _as_problem(axes, quad, mats, mat_map=mat_map, boundary_conditions=conds or {})


def _two_region(cells):
    mm = np.zeros(cells, int)
    mm[cells[0] // 2:] = 1
    return mm


_GL8 = Quadrature.gauss_legendre(n_ordinates=8)
_LS2 = Quadrature.level_symmetric(sn_order=2)
_WHITE = {0: WallCondition.white(1.0)}

#: (label, builder, expected physical dimension). d = 1 has no DD freedom, so its kernel is the physical one alone.
_PHYSICAL = [
    ("d1 (6,) GL8 graded1.5 mirror, feeder", lambda: _problem((6,), [(_R, _R)], _GL8, {0: FEEDER}, stretch=1.5), 1),
    ("d1 (6,) GL8 white(1.0), feeder", lambda: _problem((6,), [(0, 0)], _GL8, {0: FEEDER}, conds=_WHITE), 1),
    ("d1 (6,) GL8 two regions, density-scaled", lambda: _problem((6,), [(_R, _R)], _GL8, {0: FEEDER, 1: FEEDER_X2},
                                                                  mat_map=_two_region((6,))), 1),
    ("d2 (3,4) LS2 mirror, n2n in the class", lambda: _problem((3, 4), [(_R, _R)] * 2, _LS2, {0: N2N}), 1),
]


# ── S16a / S16b: the diamond-difference left basis ───────────────────────────────────────────────────────

_DD = [("LS2 ng=1", dict(ng=1, quad=_LS2)), ("LS4 ng=2", dict(ng=2))]


@pytest.mark.parametrize("label,kw", _DD, ids=[r[0] for r in _DD])
def test_S16a_the_DD_left_basis_spans_the_dense_left_null_space(label, kw):
    """REFERENCE: dense SVD at c = 0.9 (the DD kernel alone). `[M]` the left null vectors are trace-only and lie outside
    span(Phi), span(G Phi), span(G^-1 Phi) (residual 1.00 each), so no re-dressing of the right basis passes; the
    predicted left summands' supports are trace dofs only (exact zeros on the bulk)."""
    problem = _mesh((3, 4), [(_R, _R)] * 2, c=0.9, **kw)
    A = _dense(problem)
    L, R = _dense_kernels(A)
    kernel = predicted_discretization_kernel(problem)
    n, nb = A.shape[0], int(np.prod(problem.full_field_space.interior_space.shape))
    left, right = _columns(kernel.left, n), _columns(kernel.right, n)
    assert L.shape[1] == R.shape[1] == predicted_kernel_dimension(problem) > 0
    _assert_spans_agree(f"{label} left", left, L)
    _assert_spans_agree(f"{label} right", right, R)
    assert min(int(np.min(s.support)) for s in kernel.left) >= nb
    np.testing.assert_array_equal(left[:nb], 0.0)


def _independent_trace_weights(problem) -> np.ndarray:
    r"""``-w_n (Omega_n . n_f) A_f`` per flat trace slot, from the quadrature, the trace layout and the cell edges."""
    quad = problem.quad
    mu = np.stack([np.asarray(quad.mu_x, float), np.asarray(quad.mu_y, float), np.asarray(quad.mu_z, float)], axis=1)
    w = np.asarray(quad.weights, float)
    ndim = int(problem.ndim)
    out = np.zeros(int(np.prod(problem.angular_trace.shape)))
    for face, slot in problem.angular_trace.layout.faces.items():
        axis, sign = face_normal(face)
        area = np.ones(())
        for c in (c for c in range(ndim) if c != axis):
            area = np.multiply.outer(area, np.diff(np.asarray(problem.axes[c].edges, float)))
        per_ordinate = -w * sign * mu[:, axis]
        values = per_ordinate[:, None, *([None] * area.ndim)] * np.asarray(area)[None, None]
        out[slot.offset:slot.offset + int(np.prod(slot.shape))] = np.broadcast_to(values, slot.shape).ravel()
    return out


@pytest.mark.parametrize("label,kw", _DD, ids=[r[0] for r in _DD])
def test_S16b_the_DD_left_generators_are_J_times_the_right_ones(label, kw):
    r"""The MEMO's law, per summand: left generators = :math:`J\Phi`, :math:`J = -w_n(\Omega_n\cdot\hat n_f)A_f`, with J
    rebuilt here. Arms ``left_is_right`` and ``trace_weights_zeroed`` red it."""
    problem = _mesh((3, 4), [(_R, _R)] * 2, c=0.9, **kw)
    kernel = predicted_discretization_kernel(problem)
    nb = int(np.prod(problem.full_field_space.interior_space.shape))
    J = _independent_trace_weights(problem)
    assert len(kernel.left) == len(kernel.right) > 0
    for right, left in zip(kernel.right, kernel.left):
        np.testing.assert_array_equal(left.support, right.support)
        expected = J[np.asarray(right.support) - nb][:, None] * np.asarray(right.generators)
        np.testing.assert_allclose(left.generators, expected, rtol=0,
                                   atol=64 * np.finfo(float).eps * np.abs(expected).max())


# ── the physical predictor ───────────────────────────────────────────────────────────────────────────────

@pytest.mark.parametrize("label,build,dim", _PHYSICAL, ids=[r[0] for r in _PHYSICAL])
def test_the_physical_predictor_spans_the_dense_kernels(label, build, dim):
    """REFERENCE: dense SVD; the DD half added where present (d = 2). Covers a feeder group (u_0 = 0.625/0.75, the
    absorbing-chain probability), a white(1.0) wall (conserving, not a mirror), two regions of density-scaled copies and
    an (n,2n) transfer inside the class (counted twice). Arm ``e_for_u`` reds every left leg; ``trace_weights_zeroed``
    and ``left_is_right`` red the left legs."""
    problem = build()
    A = _dense(problem)
    L, R = _dense_kernels(A)
    n = A.shape[0]
    physical = predicted_physical_kernel(problem)
    dd = predicted_discretization_kernel(problem)
    assert physical.dimension == dim
    right = np.column_stack([_columns(dd.right, n), _columns(physical.right, n)])
    left = np.column_stack([_columns(dd.left, n), _columns(physical.left, n)])
    assert R.shape[1] == L.shape[1] == dd.dimension + dim
    _assert_spans_agree(f"{label} right", right, R)
    _assert_spans_agree(f"{label} left", left, L)


def test_the_physical_right_mode_is_annihilated_by_the_loss_to_round_off():
    """Object leg beside the spans (a count is blind to a rotation): ``|A r| <= 64 eps ||A||_1 ||r||_inf``."""
    problem = _PHYSICAL[0][1]()
    A = _dense(problem)
    (r,) = _columns(predicted_physical_kernel(problem).right, A.shape[0]).T
    assert np.abs(A @ r).max() <= 64 * np.finfo(float).eps * np.abs(A).sum(0).max() * np.abs(r).max()


_SCOPE_FRAGMENT = "outside the closed-form kernel's derived scope"


def _scope_witnesses():
    from orpheus.transport.spatial.linear_discontinuous import LinearDiscontinuous
    return [
        ("different energy null vectors", lambda: _problem((6,), [(_R, _R)], _GL8, {0: FEEDER, 1: OTHER_CLASS},
                                                           mat_map=_two_region((6,))), "different energy null vectors"),
        ("tangential ordinate on a closed face", lambda: _mesh((3, 4), [(_R, _R)] * 2, c=1.0, ng=1,
                                                               quad=Quadrature.product(n_mu=4, n_phi=4)), "tangential"),
        ("moment closure", lambda: _mesh((3, 4), [(_R, _R)] * 2, c=1.0, ng=1, scheme=LinearDiscontinuous()),
         "moment closure"),
    ]


@pytest.mark.parametrize("index", range(3), ids=["heterogeneous_class", "tangential_ordinate", "moment_closure"])
def test_the_closed_form_kernel_refuses_outside_its_derived_scope_at_one_door(index, monkeypatch):
    """SCOPE-BOUNDARY (ruled 2026-10-31): ONE door, ``_require_closed_form_scope``, four entries; three witnessed here
    (the non-Cartesian entry has no cheap fixture). The entry is named in the message. A spy on the door counts 1 call,
    so the refusal is the door's and no second door exists for these inputs. `[M]` physical/MEMO.md: the constant
    right form fails at 2.7e-1 on the heterogeneous class."""
    import orpheus.sn.operators.loss_kernel_gauge as M
    label, build, entry = _scope_witnesses()[index]
    problem = build()
    calls = []
    honest = M._require_closed_form_scope
    monkeypatch.setattr(M, "_require_closed_form_scope", lambda *a, **k: (calls.append(1), honest(*a, **k))[1])
    with pytest.raises(ValueError, match=_SCOPE_FRAGMENT) as info:
        predicted_physical_kernel(problem)
    assert entry in str(info.value), (label, str(info.value))
    assert calls == [1]


# ── S8: the union count ──────────────────────────────────────────────────────────────────────────────────

_UNION = [("LS2 ng=1 c=0.9", dict(ng=1, quad=_LS2, c=0.9), 2, 0), ("LS2 ng=1 c=1", dict(ng=1, quad=_LS2, c=1.0), 2, 1),
          ("LS4 ng=2 c=0.9", dict(ng=2, c=0.9), 12, 0), ("LS4 ng=2 c=1", dict(ng=2, c=1.0), 12, 2)]


@pytest.mark.parametrize("label,kw,dd,classes", _UNION, ids=[r[0] for r in _UNION])
def test_S8_the_union_count_and_spans(label, kw, dd, classes):
    """THEOREM, dense: dim ker A = dim ker A^T = DD law + |absorption-free classes| (`[M]` 2, 3, 12, 14), and the hub's
    predicted kernel spans both dense kernels. At c = 1 the diagonal-only scattering makes each group its own class."""
    problem = _mesh((3, 4), [(_R, _R)] * 2, **kw)
    A = _dense(problem)
    L, R = _dense_kernels(A)
    kernel = problem.predicted_kernel
    n = A.shape[0]
    assert predicted_kernel_dimension(problem) == dd
    assert R.shape[1] == L.shape[1] == kernel.dimension == dd + classes
    _assert_spans_agree(f"{label} right", _columns(kernel.right, n), R)
    _assert_spans_agree(f"{label} left", _columns(kernel.left, n), L)
