r"""``EigenPosing`` / ``SourcePosing`` — their DEFINING laws (consumers campaign step 2, C3b; re-based layer DAG Phase 2).

Gate specification: ``scratch/_consumers/planning/test_architect_step2_delta.md``
§B3 and the C3 delta §D.3; re-based onto :mod:`tests.gates.numerics._energy_pencil`
(``scratch/layer_dag/refiling/specs_2r/numerics_laws.md``).  ⚠ Row (d) is an ARITY
gate, not an adjoint-correctness gate (``k† = k`` is a theorem the tree asserts
elsewhere, and ``eig(Mᵀ) = eig(M)`` puts the factor-order family inside the shared
spectrum's stabiliser).

The leg "the fixture's k IS the hub's k_inf" left with the homogeneous fixture:
the homogeneous k∞ is pinned against its exact rational value by
``tests/gates/homogeneous/test_kinf_exact_reference.py``.
"""
from __future__ import annotations


import numpy as np
import pytest

from tests.gates.numerics._energy_pencil import FISSILE, EnergyPencil, bound

from orpheus.numerics.pencil import OperatorPencil
from orpheus.numerics.posing import ALPHA_MAP, K_MAP, EigenPosing, SourcePosing

pytestmark = pytest.mark.foundation


def _require(cond: object, msg: str) -> None:
    if not cond:
        raise AssertionError(msg)


def _exact_pair(ng: int = 2):
    """The exact eigenpair of the numerics pencil: k and the dominant eigenvector of A⁻¹F."""
    p = EnergyPencil.of(FISSILE[ng])
    k, phi = p.exact_pair
    return p, p.pencil, k, phi


@pytest.mark.rests_on("tests/gates/numerics/test_pencil.py::TestLawTheAffineFamily::test_law_at_zero_is_the_lhs", "tests/gates/numerics/test_pencil.py::TestLawTheAffineFamily::test_law_at_sigma_is_the_affine_combination")
class TestLawTheEigenPosing:
    def test_law_the_rayleigh_at_the_solution_is_the_eigenvalue(self) -> None:
        """At the exact eigenpair the Rayleigh quotient with ANY weight is the eigenvalue."""
        _, pencil, k, phi = _exact_pair()
        posing = EigenPosing(pencil, K_MAP)
        _require(abs(posing.rayleigh(phi, w=1.0) - k) <= 1e-12 * k, "rayleigh(φ*, w=1) must be k")
        _require(abs(posing.rayleigh(phi, w=phi) - k) <= 1e-12 * k, "…and with the stationary weight too")

    def test_law_the_eigen_residual_vanishes_at_the_solution(self) -> None:
        _, pencil, k, phi = _exact_pair()
        posing = EigenPosing(pencil, K_MAP)
        r = np.asarray(posing.residual(phi, k)); scale = np.max(np.abs(np.asarray(pencil.lhs.apply(phi))))
        _require(np.max(np.abs(r)) <= 1e-12 * scale, "A φ* − F φ*/k* must vanish")

    def test_law_the_balance_functional_is_ONE_functional(self) -> None:
        """β(ψ, λ) = ⟨1, 𝒜(μ(λ))ψ − q⟩: the eigen kind with q = 0, the source kind
        with the member at μ and q — the SAME number, bit-identically."""
        _, pencil, k, phi = _exact_pair()
        rng = np.random.default_rng(0); psi = rng.random(phi.shape); lam = 0.7 * k
        eigen = EigenPosing(pencil, K_MAP).balance(psi, lam)
        source = SourcePosing(pencil.at(K_MAP.inverse(lam)), np.zeros_like(psi)).balance(psi)
        _require(eigen == source, f"one functional, two kinds: {eigen!r} vs {source!r}")


class TestLawTheAdjointArity:
    @pytest.mark.rests_on("tests/gates/numerics/test_pencil.py::TestLawTheAdjoint::test_law_H_is_the_pair_of_adjoints")
    def test_law_the_eigen_H_is_NULLARY(self) -> None:
        _, pencil, k, _ = _exact_pair()
        adj = EigenPosing(pencil, K_MAP).H()
        A_H = np.asarray(adj.pencil.lhs.as_matrix(), dtype=float); F_H = np.asarray(adj.pencil.rhs.as_matrix(), dtype=float)
        k_adj = float(np.max(np.linalg.eigvals(np.linalg.solve(A_H, F_H)).real))
        _require(abs(k_adj - k) < 1e-12, "k† = k (an ARITY gate — the adjoint's correctness lives in the reciprocity suite)")

    def test_law_the_source_H_is_UNARY(self) -> None:
        _, pencil, _, phi = _exact_pair()
        posing = SourcePosing(pencil.lhs, np.ones_like(phi))
        with pytest.raises(TypeError):
            posing.H()  # type: ignore[call-arg]  # Python's own TypeError: the PRODUCTION signature demands a detector
        adj = posing.H(np.ones_like(phi))
        _require(np.array_equal(np.asarray(adj.source), np.ones_like(phi)), "H(detector) poses the detector as the adjoint source")


@pytest.mark.rests_on("tests/gates/numerics/test_pencil.py::TestLawTheAffineFamily::test_law_at_zero_is_the_lhs", "tests/gates/numerics/test_pencil.py::TestLawTheAffineFamily::test_law_at_sigma_is_the_affine_combination")
class TestLawTheSourcePosing:
    def test_law_the_balance_and_the_residual_close_at_the_solution(self) -> None:
        p, pencil, _, phi = _exact_pair()
        A, _ = p.matrices
        q = np.ones(p.space.shape); psi = np.linalg.solve(A, q.ravel()).reshape(p.space.shape)
        posing = SourcePosing(pencil.lhs, q)
        _require(np.max(np.abs(np.asarray(posing.residual(psi)))) <= 1e-12 * np.max(np.abs(q)), "Aψ* − q vanishes")
        _require(np.all(np.asarray(posing.residual(2.0 * psi)) > 0), "the residual is SIGNED with the loss convention: too much flux reads Aψ − q > 0")
        _require(abs(posing.balance(psi)) <= 1e-12 * np.sum(np.abs(q)), "⟨1, Aψ* − q⟩ vanishes")
        _require(posing.balance(2.0 * psi) > 0, "the balance is SIGNED: too much flux reads positive (Aψ − q > 0)")


class TestLawEachSpectralMapReadsThePencilsEigenvalue:
    r"""A spectral map is the bijection between the PENCIL's eigenvalue μ,
    :math:`A\psi = \mu M\psi`, and the physical one; what pins a map is the
    physical equation it names, solved here from the definition and never from
    the map (ERR-089).
    """

    @staticmethod
    def _alpha_pencil():
        r"""The prompt-α pencil :math:`(A - F)\psi = -\alpha\,T\psi`, :math:`T = 1/v` (speeds 3 and 1, any unit)."""
        p = EnergyPencil.of(FISSILE[2])
        prompt = p.pencil.at(1.0)
        inverse_speed = 1.0 / np.array([3.0, 1.0])
        T = bound(np.diag(inverse_speed), p.space)  # 1/v, bound to the pencil's space on both ends
        A = np.asarray(prompt.as_matrix(), dtype=float)
        mus, vecs = np.linalg.eig(np.linalg.solve(np.diag(inverse_speed), A))
        i = int(np.argmin(np.real(mus)))  # the fundamental mode: the largest α
        psi = np.real(vecs[:, i]); psi = psi / psi.sum()
        alpha = -(A @ psi) / (inverse_speed * psi)  # α from the physical equation alone
        _require(np.ptp(alpha) <= 1e-9 * abs(alpha[0]), "ψ is an eigenvector: one α in every group")
        return OperatorPencil(prompt, T), psi.reshape(p.space.shape), float(alpha[0])

    @pytest.mark.catches("ERR-089")
    def test_law_the_alpha_map_solves_the_alpha_equation(self) -> None:
        pencil, psi, alpha = self._alpha_pencil()
        posing = EigenPosing(pencil, ALPHA_MAP)
        _require(alpha > 0.0, "the fixture is supercritical (k = 1.39), so the fundamental α is positive")
        r = np.asarray(posing.residual(psi, alpha)); scale = np.max(np.abs(np.asarray(pencil.lhs.apply(psi))))
        _require(np.max(np.abs(r)) <= 1e-12 * scale, "(A − μ(α) T) ψ must vanish at the physical α")
        _require(abs(posing.rayleigh(psi) - alpha) <= 1e-12 * abs(alpha), "the Rayleigh quotient must be the physical α")

    def test_law_the_k_map_solves_the_k_equation(self) -> None:
        p, pencil, k, phi = _exact_pair()
        A, F = p.matrices
        k_def = float((F @ phi).sum() / (A @ phi).sum())  # A φ = F φ / k
        _require(abs(K_MAP(K_MAP.inverse(k_def)) - k_def) <= 1e-14 * k_def, "the round trip")
        _require(abs(EigenPosing(pencil, K_MAP).rayleigh(phi) - k_def) <= 1e-12 * k_def, "the Rayleigh quotient is k")
