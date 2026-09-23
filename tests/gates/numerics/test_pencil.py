r"""``OperatorPencil`` — its DEFINING laws (consumers campaign step 2, C3b).

Gate specification: ``scratch/_consumers/planning/test_architect_step2_delta.md``
§B2 (rows a–h, j; row i dropped — ``[M]`` the shipped A mixtures carry no
(n,2n) channel, so "a polluted rhs raises the rank" has no shipped witness).

⛔ A limit: ``[M]`` ``at(σ).as_matrix()`` RAISES on the SN composite
(``OperatorSum.as_matrix`` probes with bare arrays; ``CoupledOperator.apply``
refuses), so the matrix-form laws (d)/(f)/(h) are 0-D-only; the SN arm gets the
ends law (a)/(b), the AC-a gate and the subcritical-source witness.
⛔ Declared blindness for (f): the 0-D pose is a rank-1 (point) spatial axis, so
every metric weight is a scalar and a scalar Gram commutes with everything — (f)
cannot see a metric error; the metric-loaded partner is the SN composite's
reciprocity suite (``tests/gates/sn/operators/test_g_adjoint_reciprocity.py``).
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.derivations.common.xs_library import get_mixture
from orpheus.homogeneous.solver import HomogeneousProblem, solve_homogeneous_infinite
from orpheus.numerics.operator import LinearOperator
from orpheus.numerics.pencil import OperatorPencil

pytestmark = pytest.mark.foundation


def _require(cond: object, msg: str) -> None:
    if not cond:
        raise AssertionError(msg)


def _zero_d(fam: str = "A", ng: str = "2g") -> tuple[HomogeneousProblem, OperatorPencil]:
    problem = HomogeneousProblem(get_mixture(fam, ng))
    return problem, OperatorPencil(problem.loss, problem.production)


def _matrices(pencil: OperatorPencil) -> tuple[np.ndarray, np.ndarray]:
    return np.asarray(pencil.lhs.as_matrix(), dtype=float), np.asarray(pencil.rhs.as_matrix(), dtype=float)


class TestLawTheEndsLaw:
    @pytest.mark.parametrize("pair", ["the sharp pair (equal shapes, different spaces)", "the coarse pair (different shapes)"])
    def test_law_a_mismatched_pair_is_REFUSED(self, pair: str) -> None:
        """``[M]`` the sharp pair: ``(rec.loss, hub.fission)`` — ``CoupledSpace(160,)``
        vs ``FullFieldSpace(160,)``, shapes EQUAL — a SHAPE-keyed ends check would
        admit it; the coarse pair is the plan's named witness."""
        from orpheus.sn.coupled_system import build_within_group_system
        from tests.gates.sn.operators.test_step2_posed_fission_anchors import _slab_hub

        hub = _slab_hub(); rec = build_within_group_system(hub, hub.mat_xs)
        rhs = hub.fission if pair.startswith("the sharp") else hub.fission.isotropic_energy
        with pytest.raises(ValueError, match="ends law"):
            OperatorPencil(rec.loss, rhs)

    def test_law_the_posed_pair_is_ADMITTED(self) -> None:
        from orpheus.sn.coupled_system import build_within_group_system
        from tests.gates.sn.operators.test_step2_posed_fission_anchors import _slab_hub

        hub = _slab_hub(); rec = build_within_group_system(hub, hub.mat_xs)
        pencil = OperatorPencil(rec.loss, rec.production)
        _require(pencil.lhs is rec.loss and pencil.rhs is rec.production, "the pencil holds the posed pair by identity")
        _require(hub.pencil.lhs is hub.system.loss, "the hub's pencil is over the hub's record")


class TestLawTheAffineFamily:
    def test_law_at_zero_is_the_lhs(self) -> None:
        """``[M]`` a mutation reds by RAISING (``ScaledOperator(0.0, ·)`` is refused) — paired with the affine row."""
        _, pencil = _zero_d()
        _require(pencil.at(0.0) is pencil.lhs, "at(0) must be the lhs object itself")

    @pytest.mark.parametrize("sigma", [0.5, 1.0, 1.875, 2.0])
    def test_law_at_sigma_is_the_affine_combination(self, sigma: float) -> None:
        """``[M]`` ``array_equal`` True, ``max|Δ| = 0.0`` at all four σ (matrix form)."""
        _, pencil = _zero_d()
        A, F = _matrices(pencil)
        _require(np.array_equal(np.asarray(pencil.at(sigma).as_matrix(), dtype=float), A - sigma * F),
                 "at(σ) must be A − σ·M exactly in matrix form")

    def test_law_the_increment_is_minus_tau_M(self) -> None:
        """The APPLIED form is banded, not bit-identical (``[M]`` 144/200 draws differ,
        draw-stable statistic ≤ 9.470e-16 over 300 draws ≈ 4.3 ε): gate at 8·ε·scale."""
        _, pencil = _zero_d()
        A, F = _matrices(pencil)
        rng = np.random.default_rng(0)
        eps = np.finfo(float).eps
        for _ in range(300):
            x = rng.random(A.shape[1]); sigma, tau = rng.random(2)
            d = np.asarray(pencil.at(sigma + tau).as_matrix() @ x) - np.asarray(pencil.at(sigma).as_matrix() @ x)
            ref = -tau * (F @ x)
            scale = max(np.max(np.abs(A @ x)), np.max(np.abs(sigma * (F @ x))), 1e-300)
            _require(np.max(np.abs(d - ref)) / scale <= 8 * eps, "the increment at(σ+τ) − at(σ) must be −τ·M to 8 ε")


class TestLawTheAdjoint:
    def test_law_H_is_the_pair_of_adjoints(self) -> None:
        _, pencil = _zero_d()
        A, F = _matrices(pencil)
        _require(np.array_equal(np.asarray(pencil.at(1.0).H.as_matrix(), dtype=float), (A - F).T),
                 "the adjoint pencil's member at σ=1 is (A − F)ᵀ (0-D: a scalar Gram commutes — see the module docstring)")
        _require(np.array_equal(np.asarray(pencil.H.rhs.as_matrix(), dtype=float), F.T), ".H daggers the rhs too")


class TestLawTheRank:
    @pytest.mark.parametrize("fam,rank", [("A", 1), ("B", 0), ("C", 0), ("D", 0)])
    @pytest.mark.parametrize("ng", ["1g", "2g", "4g"])
    def test_law_rhs_rank_counts_the_finite_spectrum(self, fam: str, ng: str, rank: int) -> None:
        """Weierstrass–Kronecker: exactly rank(M) finite eigenvalues — ``[M]`` A: 1, B/C/D: 0 (a two-sided gate)."""
        _, pencil = _zero_d(fam, ng)
        _require(pencil.rhs_rank == rank, f"{fam}/{ng}: rhs_rank {pencil.rhs_rank} ≠ {rank}")

    @pytest.mark.parametrize("ng", ["1g", "2g", "4g"])
    def test_law_k_inf_is_the_trace_when_rhs_is_rank_one(self, ng: str) -> None:
        """REFERENCE (closed form): ``k_inf == trace(A⁻¹F)`` — ``[M]`` exact at 1g/2g, 2.2e-16 at 4g ⟹ rtol 1e-13, not array_equal."""
        problem, pencil = _zero_d("A", ng)
        A, F = _matrices(pencil)
        k_trace = float(np.trace(np.linalg.solve(A, F)))
        k_inf = float(solve_homogeneous_infinite(problem.mixture).k_inf)
        _require(abs(k_trace - k_inf) <= 1e-13 * abs(k_inf), f"k_inf {k_inf!r} vs trace(A⁻¹F) {k_trace!r}")


class _Dense(LinearOperator):
    """An OPAQUE operator — a matrix and nothing else (no factors, no ends)."""

    def __init__(self, m: np.ndarray) -> None:
        self._m = np.asarray(m, dtype=float)

    @property
    def domain(self):  # noqa: ANN201
        return None

    @property
    def codomain(self):  # noqa: ANN201
        return None

    def apply(self, x):  # noqa: ANN001, ANN201
        return self._m @ np.asarray(x, dtype=float)

    def as_matrix(self, *a, **k):  # noqa: ANN001, ANN002, ANN003, ANN201
        return self._m


class TestTheCPOpaquePair:
    def test_law_the_constructor_takes_two_opaque_operators(self) -> None:
        """The pencil reaches for no ``.factors``/``.a``/``.b`` — CP's ``(I − K, F)`` is admitted as-is."""
        rng = np.random.default_rng(1)
        A, M = _Dense(np.eye(3) + 0.1 * rng.random((3, 3))), _Dense(rng.random((3, 3)))
        pencil = OperatorPencil(A, M)
        x = rng.random(3)
        _require(np.allclose(pencil.at(0.7).apply(x), A.apply(x) - 0.7 * M.apply(x)), "at(σ) composes through the algebra alone")
