r"""The Fredholm verdict on a posing A psi = q, decided from a predicted kernel: a numerics type, no neutrons.

Space-membership campaign, unit U0 (the user's ruling R3, 2026-10-30: "Fredholm as the decision rule, kernel measurement
as witness"; spec ``scratch/membership_covector/u0/spec_u0.md`` rows V1-V7, by the test-architect). The fixture breaks
both stabilisers of the claim: the left kernel differs from the right one (a coupling into the kernel rows), and the
metric G is non-uniform (ratio 16), so a verdict that pairs through the metric, or a Euclidean projector, answers
differently. The oracles (a dense SVD for the left kernel, the bordered system for the minimum-G-norm member) are
independent of the projector under test.
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.numerics.axis import Axis, BasisKind, Declared
from orpheus.numerics.assembled_operator import SparseAssembledOperator
from orpheus.numerics.evidence import NotApplicable
from orpheus.numerics.fredholm import (
    CONSISTENCY_RTOL, Gauge, NoSolution, PredictedKernel, Solve, fredholm_alternative,
)
from orpheus.numerics.kernel_projection_gauge import KernelSummand
from orpheus.numerics.matrix_inverse_operator import MatrixInverseOperator
from orpheus.numerics.space import FunctionSpace

pytestmark = [pytest.mark.l0]

_N = 9
_G = np.array([0.5, 1.0, 2.0, 0.25, 3.0, 1.5, 0.75, 4.0, 1.25])
SPACE = FunctionSpace.of_axes(Axis("dof", (_N,), weights=_G, origin=Declared(BasisKind.NODAL)))
_K = np.array([2, 7])
_R = np.setdiff1d(np.arange(_N), _K)
_rng = np.random.default_rng(20261031)
_A = np.zeros((_N, _N))
_A[np.ix_(_R, _R)] = _rng.standard_normal((_R.size, _R.size)) + 4.0 * np.eye(_R.size)
_A[np.ix_(_K, _R)] = _rng.standard_normal((_K.size, _R.size))   # kernel ROWS coupled: left kernel != right kernel
_RIGHT = np.eye(_N)[:, _K]                                         # A e_K = 0 exactly (columns K are zero)
_LEFT = np.linalg.svd(_A)[0][:, -_K.size:]                         # dense left null vectors (independent oracle)
_ALL = np.arange(_N)
_SOLVED = NotApplicable("the fixture factorises nothing")

_LABEL = "left_covector_{}"


def _kernel(left_columns=None):
    """The predicted kernel: the right summand e_K, one left summand per dense left null vector (labelled)."""
    columns = range(_K.size) if left_columns is None else left_columns
    return PredictedKernel(
        SPACE, (KernelSummand(_K, np.eye(_K.size), "right"),),
        tuple(KernelSummand(_ALL, _LEFT[:, [j]], _LABEL.format(j)) for j in columns), "the fixture",
    )


def _verdict(q, kernel=None):
    return fredholm_alternative(_kernel() if kernel is None else kernel, q, conditioning=_SOLVED)


def _balanced(seed):
    q = _A @ np.random.default_rng(seed).standard_normal(_N)          # in range(A): every left co-vector pairs to 0
    return q


def test_premise_left_and_right_kernels_differ_and_G_is_not_uniform():
    assert np.linalg.norm(_LEFT.T @ _RIGHT) < 0.999 * _K.size
    assert _G.max() / _G.min() == 16.0
    np.testing.assert_array_equal(_A @ _RIGHT, 0.0)


@pytest.mark.verifies("fredholm-posing-alternative")
def test_V1_a_trivial_kernel_is_Solve():
    v = _verdict(np.ones(_N), PredictedKernel.trivial(SPACE, "no kernel"))
    assert isinstance(v, Solve) and v.because == "no kernel" and v.conditioning is _SOLVED


@pytest.mark.verifies("fredholm-posing-alternative")
def test_V2_q_orthogonal_to_every_left_covector_is_Gauge_and_the_pairing_is_canonical():
    """Arm ``verdict_pairs_with_metric`` (pairs G w with q) reds: the balanced q is canonically orthogonal, not G-orth."""
    q = _balanced(1)
    assert max(abs(float(_LEFT[:, j] @ q)) for j in range(_K.size)) < 1e-12 * np.linalg.norm(q)
    assert max(abs(float((_G * _LEFT[:, j]) @ q)) for j in range(_K.size)) > 1e-3 * np.linalg.norm(q)  # the arm's reading
    v = _verdict(q)
    assert isinstance(v, Gauge) and v.projector.dimension == _K.size and v.pairing.value <= CONSISTENCY_RTOL


@pytest.mark.verifies("fredholm-posing-alternative")
def test_V3_one_nonzero_pairing_is_Refuse_and_names_the_covector():
    q = _balanced(2) + 0.1 * np.linalg.norm(_balanced(2)) * _LEFT[:, 1]
    with pytest.raises(NoSolution, match="no steady solution") as info:
        _verdict(q)
    assert repr(_LABEL.format(1)) in str(info.value) and repr(_LABEL.format(0)) not in str(info.value)


@pytest.mark.verifies("fredholm-posing-alternative")
def test_V4_the_gauged_member_is_the_minimum_G_norm_member():
    """REFERENCE: one dense bordered solve [[A, R], [(G R)^T, 0]] (never the SUT's projector). Arm ``euclidean_projector``."""
    q = _balanced(3)
    v = _verdict(q)
    x_any = np.linalg.lstsq(_A, q, rcond=None)[0] + _RIGHT @ np.array([3.0, -2.0])
    B = np.block([[_A, _RIGHT], [(_G[:, None] * _RIGHT).T, np.zeros((2, 2))]])
    x_ref = np.linalg.solve(B, np.concatenate([q, [0.0, 0.0]]))[:_N]
    np.testing.assert_allclose(v.projector.gauge(x_any), x_ref, rtol=0, atol=64 * np.finfo(float).eps * np.linalg.cond(B) * np.abs(x_ref).max())


def _cosine(w, q, support):
    """The cut's statistic, by hand: the coordinate cosine on the annihilator's declared support; 0 when q vanishes there."""
    den = float(np.linalg.norm(w[support]) * np.linalg.norm(q[support]))
    return 0.0 if den == 0.0 else abs(float(w[support] @ q[support])) / den


@pytest.mark.parametrize("t,admitted", [(0.0, True), (1e-14, True), (1e-6, False), (1e-2, False)],
                         ids=["balanced", "round_off", "1e-6", "1e-2"])
@pytest.mark.verifies("fredholm-posing-alternative")
def test_V5_the_cut_reads_the_cosine_on_the_declared_support(t, admitted):
    """Arm ``cut_1e-2`` reds the 1e-6 row; arm ``cosine_with_metric`` (G in the pairing) reds the balanced row.
    [M] U0 window (probes/p_cosine.py): largest Gauge 2.4e-15, smallest Refuse 0.092."""
    q = _balanced(4) + t * np.linalg.norm(_balanced(4)) * _LEFT[:, 0]
    full = np.arange(_N)
    c = max(_cosine(_LEFT[:, j], q, full) for j in range(_K.size))
    assert (c <= 1e-10) is admitted
    if admitted:
        assert isinstance(_verdict(q), Gauge)
    else:
        with pytest.raises(NoSolution):
            _verdict(q)
    assert CONSISTENCY_RTOL == 1e-10


@pytest.mark.verifies("fredholm-posing-alternative")
def test_V5b_a_source_off_the_support_does_not_dilute_the_cosine():
    """A source component outside an annihilator's declared support leaves c unchanged (the denominator reads q_S only).
    Arm ``cosine_full_norm`` (||q||_2 over every dof) reds it: the large off-support component would hide the pairing."""
    supp = np.array([0, 1, 3, 4, 5, 6, 8])
    w = np.zeros(_N); w[supp] = np.random.default_rng(9).standard_normal(supp.size)
    q = np.zeros(_N); q[supp] = w[supp]                          # c = 1 on the support
    q_off = q.copy(); q_off[[2, 7]] = 1e12                       # off the support
    assert _cosine(w, q_off, supp) == pytest.approx(1.0)
    kernel = PredictedKernel(SPACE, (KernelSummand(np.array([2]), np.ones((1, 1)), "r"),),
                             (KernelSummand(supp, w[supp, None], "w"),), "the dilution fixture")
    assert kernel.consistency(q_off)[0] == pytest.approx(1.0)


@pytest.mark.verifies("fredholm-posing-alternative")
def test_V6_unequal_left_and_right_counts_are_refused_at_construction():
    """A square A has Fredholm index 0: as many left co-vectors as right kernel vectors."""
    with pytest.raises(ValueError, match="index"):
        _kernel(left_columns=[0])


def _projector(R, G):
    """The U3 projector onto span(R), G-orthogonal, built by hand: R (R^T G R)^-1 R^T G."""
    return R @ np.linalg.solve(R.T @ (G[:, None] * R), (G[:, None] * R).T)


@pytest.mark.verifies("fredholm-deflated-loss")
def test_V7a_the_deflated_inverse_returns_the_minimum_G_norm_member():
    """RULED (the user, 2026-10-31): a factorising method inverts A + Pi. For a consistent q, (A + Pi)^-1 q is the
    minimum-G-norm member. REFERENCE: the bordered system (test oracle only). [M] orchestrator (deflated_probe.py,
    3 diffusion fixtures): cond(A + Pi) 35-42 against 1e17, agreement 1e-15. Arm deflate_euclidean (R R^T) reds."""
    q = _balanced(5)
    A = SparseAssembledOperator(_A, domain=SPACE, codomain=SPACE)
    x = np.ravel(MatrixInverseOperator(_kernel().deflate(A)).apply(q))
    B = np.block([[_A, _RIGHT], [(_G[:, None] * _RIGHT).T, np.zeros((2, 2))]])
    x_ref = np.linalg.solve(B, np.concatenate([q, [0.0, 0.0]]))[:_N]
    np.testing.assert_allclose(x, x_ref, rtol=0, atol=64 * np.finfo(float).eps * np.linalg.cond(_A + _projector(_RIGHT, _G)) * np.abs(x_ref).max())


@pytest.mark.verifies("fredholm-deflated-loss")
def test_V7b_a_non_semisimple_kernel_is_refused_at_construction():
    """A + Pi is nonsingular iff L^T R is nonsingular. Jordan fixture: A = [[0, 1], [0, 0]] (+) I, ker A = e_0, ker A^T =
    e_1, L^T R = 0: A + Pi is singular. The SN DD kernel is of this kind [M] (L^T R DD columns 0), which is why SN keeps
    iterate-then-project. Fragment pinned at the carve."""
    A = np.eye(4); A[0, 0] = 0.0; A[1, 1] = 0.0; A[0, 1] = 1.0
    sp = FunctionSpace.of_axes(Axis("dof", (4,), weights=np.ones(4), origin=Declared(BasisKind.NODAL)))
    jordan = PredictedKernel(sp, (KernelSummand(np.array([0]), np.ones((1, 1)), "right"),),
                             (KernelSummand(np.array([1]), np.ones((1, 1)), "left"),), "a Jordan chain")
    with pytest.raises(ValueError, match="semisimple"):
        jordan.deflate(SparseAssembledOperator(A, domain=sp, codomain=sp))
