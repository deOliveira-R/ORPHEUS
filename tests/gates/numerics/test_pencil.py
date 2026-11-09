r"""``OperatorPencil`` — its DEFINING laws (consumers campaign step 2, C3b; re-based layer DAG Phase 2).

Gate specification: ``scratch/_consumers/planning/test_architect_step2_delta.md``
§B2 (rows a–h, j). Re-based (``scratch/layer_dag/refiling/specs_2r/numerics_laws.md``)
onto :mod:`tests.gates.numerics._energy_pencil`, a pencil built from numerics
objects alone.

The ends law's SHARP pair is reproduced in numerics: two spaces of EQUAL shape
and different content (a synthetic and a partitioned energy axis), so a
shape-keyed ends check still admits it and this gate still refuses it. The SN
instance (``CoupledSpace(160,)`` vs ``FullFieldSpace(160,)``) and the hub
identity ``hub.pencil.lhs is hub.system.loss`` are SN claims, asserted at
``tests/gates/sn/architecture/test_posing.py``.

The adjoint row is no longer metric-blind: the 0-D counting pose made every
Gram a scalar (the old module docstring's declared blindness for (f)); the
weighted-axis row poses the pencil on a non-scalar diagonal Gram, where
:math:`\mathcal{A}(\sigma)^\dagger = G^{-1}\mathcal{A}(\sigma)^{\mathsf T} G \ne \mathcal{A}(\sigma)^{\mathsf T}`.
"""
from __future__ import annotations


import re

import numpy as np
import pytest

from tests.gates.numerics._energy_pencil import FISSILE, POINT, EnergyPencil, GroupData, bound, energy_space, non_fissile, two_spectra, weighted_space

from orpheus.numerics.axis import EnergyAxis
from orpheus.numerics.operator import LinearOperator
from orpheus.numerics.pencil import OperatorPencil
from orpheus.numerics.space import FunctionSpace

pytestmark = pytest.mark.foundation


def _require(cond: object, msg: str) -> None:
    if not cond:
        raise AssertionError(msg)


def _fissile(ng: int = 2) -> EnergyPencil:
    return EnergyPencil.of(FISSILE[ng])


def _partitioned_space(ng: int) -> FunctionSpace:
    """Equal SHAPE to ``energy_space(ng)``, different CONTENT (a real group partition)."""
    edges = np.geomspace(2.0e7, 1.0e-5, ng + 1)
    return FunctionSpace.of_axes(EnergyAxis.of_edges(edges), POINT)


class TestLawTheEndsLaw:
    @pytest.mark.parametrize("pair", ["the sharp pair (equal shapes, different spaces)", "the coarse pair (different shapes)"])
    def test_law_a_mismatched_pair_is_REFUSED(self, pair: str) -> None:
        """The sharp pair: the same matrix bound to ``energy(2)⊗point`` and to a
        partitioned ``energy(2)⊗point`` — shapes EQUAL, so a SHAPE-keyed ends
        check would admit it; the coarse pair differs in shape."""
        p = _fissile(2)
        other = _partitioned_space(2) if pair.startswith("the sharp") else energy_space(4)
        _require(other != p.space, "the fixture's two spaces are distinct")
        if pair.startswith("the sharp"):
            _require(other.shape == p.space.shape, "the sharp pair has EQUAL shapes")
        rhs = bound(np.eye(other.shape[0]), other)
        with pytest.raises(ValueError, match="ends law"):
            OperatorPencil(p.loss, rhs)

    def test_law_the_posed_pair_is_ADMITTED(self) -> None:
        p = _fissile(2)
        pencil = OperatorPencil(p.loss, p.production)
        _require(pencil.lhs is p.loss and pencil.rhs is p.production, "the pencil holds the posed pair by identity")


class TestLawTheAffineFamily:
    def test_law_at_zero_is_the_lhs(self) -> None:
        """``[M]`` the short-circuit ``at(0.0) is lhs`` — paired with the affine row."""
        pencil = _fissile().pencil
        _require(pencil.at(0.0) is pencil.lhs, "at(0) must be the lhs object itself")

    @pytest.mark.parametrize("sigma", [0.5, 1.0, 1.875, 2.0])
    def test_law_at_sigma_is_the_affine_combination(self, sigma: float) -> None:
        p = _fissile()
        A, F = p.matrices
        _require(np.array_equal(np.asarray(p.pencil.at(sigma).as_matrix(), dtype=float), A - sigma * F),
                 "at(σ) must be A − σ·M exactly in matrix form")

    def test_law_the_increment_is_minus_tau_M(self) -> None:
        """The APPLIED form is banded, not bit-identical: gate at 8·ε·scale over 300 draws."""
        p = _fissile()
        A, F = p.matrices
        rng = np.random.default_rng(0)
        eps = np.finfo(float).eps
        for _ in range(300):
            x = rng.random(p.space.shape); sigma, tau = rng.random(2)
            d = np.asarray(p.pencil.at(sigma + tau).apply(x)).ravel() - np.asarray(p.pencil.at(sigma).apply(x)).ravel()
            x = x.ravel(); ref = -tau * (F @ x)
            scale = max(np.max(np.abs(A @ x)), np.max(np.abs(sigma * (F @ x))), 1e-300)
            _require(np.max(np.abs(d - ref)) / scale <= 8 * eps, "the increment at(σ+τ) − at(σ) must be −τ·M to 8 ε")


class TestLawTheAdjoint:
    def test_law_H_is_the_pair_of_adjoints(self) -> None:
        p = _fissile()
        A, F = p.matrices
        _require(np.array_equal(np.asarray(p.pencil.at(1.0).H.as_matrix(), dtype=float), (A - F).T),
                 "the adjoint pencil's member at σ=1 is (A − F)ᵀ on the counting pose")
        _require(np.array_equal(np.asarray(p.pencil.H.rhs.as_matrix(), dtype=float), F.T), ".H daggers the rhs too")

    def test_law_H_reads_the_metric_on_a_weighted_space(self) -> None:
        """NEW rung (closes the 0-D declared blindness): on a diagonal Gram G ≠ cI,
        ``at(σ).H = G⁻¹ (A − σF)ᵀ G``, which is NOT the Euclidean transpose."""
        w = np.array([1.0, 3.0])
        p = EnergyPencil(FISSILE[2], weighted_space(w))
        A, F = p.matrices
        G = np.diag(w)
        metric_adjoint = np.linalg.solve(G, (A - F).T @ G)
        H = np.asarray(p.pencil.at(1.0).H.as_matrix(), dtype=float)
        _require(not np.allclose(metric_adjoint, (A - F).T), "non-vacuity: the metric adjoint differs from the transpose here")
        _require(np.allclose(H, metric_adjoint, rtol=1e-14, atol=0.0), f"at(1).H must be G⁻¹(A − F)ᵀG; got {H}")


_RANK_CASES = [(ng, fam, rank) for ng in (1, 2, 4) for fam, rank in (("fissile", 1), ("non-fissile", 0), ("two spectra", 2)) if not (fam == "two spectra" and ng == 1)]
_FAMILY = {"fissile": lambda ng: FISSILE[ng], "non-fissile": non_fissile, "two spectra": two_spectra}


class TestLawTheRank:
    @pytest.mark.parametrize("ng,family,rank", _RANK_CASES, ids=[f"{f}-{n}g" for n, f, _ in _RANK_CASES])
    def test_law_rhs_rank_counts_the_finite_spectrum(self, ng: int, family: str, rank: int) -> None:
        """Weierstrass–Kronecker: exactly rank(M) finite eigenvalues — a three-sided gate (0, 1, 2; rank 2 needs ng ≥ 2)."""
        pencil = EnergyPencil.of(_FAMILY[family](ng)).pencil
        _require(pencil.rhs_rank == rank, f"{family}/{ng}g: rhs_rank {pencil.rhs_rank} ≠ {rank}")


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


# ── the fixture refuses what it cannot pose (each refusal a ValueError, each input failing ONE guard) ──────────


_BASE = FISSILE[2]


@pytest.mark.parametrize(
    "build, fragment",
    [
        (lambda: GroupData(_BASE.sigma_t, np.zeros((2, 3)), _BASE.spectra, _BASE.productions), "must be (ng, ng)"),
        (lambda: GroupData(_BASE.sigma_t, np.array([[0.18, 0.0], [-0.07, 0.70]]), _BASE.spectra, _BASE.productions), "is negative"),
        (lambda: GroupData(_BASE.sigma_t, np.array([[0.18, 0.01], [0.07, 0.70]]), _BASE.spectra, _BASE.productions), "not lower-triangular"),
        (lambda: GroupData(np.array([0.18, 0.80]), _BASE.transfer, _BASE.spectra, _BASE.productions), "diagonal sigma_t - S_gg is not positive"),
        (lambda: GroupData(np.array([0.24, 0.80]), _BASE.transfer, _BASE.spectra, _BASE.productions), "absorption sigma_t - out-scatter is negative"),
        (lambda: GroupData(_BASE.sigma_t, _BASE.transfer, (np.array([0.9, 0.0]),), _BASE.productions), "is not a probability"),
        (lambda: EnergyPencil(_BASE, energy_space(4)), "is not the 2-group space"),
    ],
    ids=["transfer-shape", "negative-transfer", "up-scatter", "singular-diagonal", "negative-absorption", "chi-not-a-probability", "wrong-space"],
)
def test_negative_the_fixture_refuses_an_illegal_medium_or_a_foreign_space(build, fragment: str) -> None:
    with pytest.raises(ValueError, match=re.escape(fragment)):
        build()
