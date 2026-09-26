r"""``EigenPosing`` / ``SourcePosing`` — their DEFINING laws (consumers campaign step 2, C3b).

Gate specification: ``scratch/_consumers/planning/test_architect_step2_delta.md``
§B3 and the C3 delta §D.3.  ⚠ Row (d) is an ARITY gate, not an adjoint-correctness
gate (``k† = k`` is a theorem the tree asserts elsewhere, and ``eig(Mᵀ) = eig(M)``
puts the factor-order family inside the shared spectrum's stabiliser).
⚠ Row (a) at the SN tier is a REFERENCE-class agreement: ``[M]`` (C3 delta §D.3)
``rayleigh(ψ, w=1)`` and ``compute_keff`` differ by the CONVERGENCE RESIDUAL
(rel 7.2e-10 at the default tolerances, 6.6e-14 with more outers) — they are one
functional at the solution, not one body on an iterate.
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.derivations.common.xs_library import get_mixture
from orpheus.homogeneous.solver import HomogeneousProblem, solve_homogeneous_infinite
from orpheus.numerics.pencil import OperatorPencil
from orpheus.numerics.operator import DiagonalOperator
from orpheus.numerics.posing import ALPHA_MAP, K_MAP, EigenPosing, SourcePosing

pytestmark = pytest.mark.foundation


def _require(cond: object, msg: str) -> None:
    if not cond:
        raise AssertionError(msg)


def _exact_pair(fam: str = "A", ng: str = "2g"):
    """The 0-D exact eigenpair: k_inf and the dominant eigenvector of A⁻¹F."""
    problem = HomogeneousProblem(get_mixture(fam, ng))
    pencil = OperatorPencil(problem.loss, problem.production)
    A = np.asarray(pencil.lhs.as_matrix(), dtype=float); F = np.asarray(pencil.rhs.as_matrix(), dtype=float)
    w, V = np.linalg.eig(np.linalg.solve(A, F))
    i = int(np.argmax(np.real(w))); phi = np.real(V[:, i]); phi = phi / phi.sum()
    # the 0-D operators are bound on the (ng, 1) energy × point-spatial carrier
    return problem, pencil, float(np.real(w[i])), phi.reshape(-1, 1)


class TestLawTheEigenPosing:
    def test_law_the_rayleigh_at_the_solution_is_the_eigenvalue(self) -> None:
        """At the exact eigenpair the Rayleigh quotient with ANY weight is the eigenvalue."""
        problem, pencil, k, phi = _exact_pair()
        posing = EigenPosing(pencil, K_MAP)
        _require(abs(posing.rayleigh(phi, w=1.0) - k) <= 1e-12 * k, "rayleigh(φ*, w=1) must be k_inf")
        _require(abs(posing.rayleigh(phi, w=phi) - k) <= 1e-12 * k, "…and with the stationary weight too")
        _require(abs(k - float(solve_homogeneous_infinite(problem.mixture).k_inf)) <= 1e-12 * k, "the fixture's k IS the hub's k_inf")

    def test_law_the_eigen_residual_vanishes_at_the_solution(self) -> None:
        _, pencil, k, phi = _exact_pair()
        posing = EigenPosing(pencil, K_MAP)
        r = np.asarray(posing.residual(phi, k)); scale = np.max(np.abs(np.asarray(pencil.lhs.apply(phi))))
        _require(np.max(np.abs(r)) <= 1e-12 * scale, "A φ* − F φ*/k* must vanish")

    def test_law_the_balance_functional_is_ONE_functional(self) -> None:
        """β(ψ, λ) = ⟨1, 𝒜(μ(λ))ψ − q⟩: the eigen kind evaluates it with q = 0, the source
        kind with the member at μ and q — the SAME number, bit-identically."""
        _, pencil, k, phi = _exact_pair()
        rng = np.random.default_rng(0); psi = rng.random(phi.shape); lam = 0.7 * k
        eigen = EigenPosing(pencil, K_MAP).balance(psi, lam)
        source = SourcePosing(pencil.at(K_MAP.inverse(lam)), np.zeros_like(psi)).balance(psi)
        _require(eigen == source, f"one functional, two kinds: {eigen!r} vs {source!r}")


class TestLawTheAdjointArity:
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


class TestLawTheSourcePosing:
    def test_law_the_balance_and_the_residual_close_at_the_solution(self) -> None:
        _, pencil, _, phi = _exact_pair()
        A = np.asarray(pencil.lhs.as_matrix(), dtype=float)
        q = np.ones((A.shape[0], 1)); psi = np.linalg.solve(A, q)
        posing = SourcePosing(pencil.lhs, q)
        _require(np.max(np.abs(np.asarray(posing.residual(psi)))) <= 1e-12 * np.max(np.abs(q)), "Aψ* − q vanishes")
        _require(np.all(np.asarray(posing.residual(2.0 * psi)) > 0), "the residual is SIGNED with the loss convention: too much flux reads Aψ − q > 0 (step 3 U1 — it read q − Aψ before)")
        _require(abs(posing.balance(psi)) <= 1e-12 * np.sum(np.abs(q)), "⟨1, Aψ* − q⟩ vanishes")
        _require(posing.balance(2.0 * psi) > 0, "the balance is SIGNED: too much flux reads positive (Aψ − q > 0)")


class TestLawEachSpectralMapReadsThePencilsEigenvalue:
    r"""A spectral map is the bijection between the PENCIL's eigenvalue μ,
    :math:`A\psi = \mu M\psi` (``at(μ)`` singular), and the physical one.

    Its round trip and its one-division form are internal laws a wrong map
    also satisfies; what pins a map is the physical equation it names, solved
    here from the definition and never from the map. ⚠ The corpus also writes μ
    for the eigenvalue of :math:`A^{-1}M` (the power iteration's), where
    :math:`k = \mu` and :math:`\alpha = -1/\mu`; in the pencil's convention
    :math:`k = 1/\mu` and :math:`\alpha = -\mu`. ``ALPHA_MAP`` once mixed the
    two (ERR-089).
    """

    @staticmethod
    def _alpha_pencil():
        r"""The 0-D prompt-α pencil: :math:`(L+C-S-F)\psi = -\alpha\,T\psi`, :math:`T = 1/v`.

        The speeds are not physical (3 and 1 in any unit): the law holds for any
        positive :math:`T`, and physical speeds four orders apart would make
        :math:`T^{-1}A` ill-scaled enough that ``np.linalg.eig`` itself limits
        the fixture to about 1e-12 (``[M]`` residual 2.5e-12 relative at
        2e9 and 2.2e5 cm/s), which is the fixture's error, not the map's.
        """
        problem = HomogeneousProblem(get_mixture("A", "2g"))
        prompt = OperatorPencil(problem.loss, problem.production).at(1.0)
        inverse_speed = 1.0 / np.array([3.0, 1.0])  # fast, thermal (any unit)
        T = DiagonalOperator(inverse_speed.reshape(-1, 1), broadcast_axes=())
        A = np.asarray(prompt.as_matrix(), dtype=float)
        mus, vecs = np.linalg.eig(np.linalg.solve(np.diag(inverse_speed), A))
        i = int(np.argmin(np.real(mus)))  # the fundamental mode: the largest α
        psi = np.real(vecs[:, i]); psi = psi / psi.sum()
        # α from the physical equation alone: −(Aψ)_g / (Tψ)_g, the same in every group
        alpha = -(A @ psi) / (inverse_speed * psi)
        _require(np.ptp(alpha) <= 1e-9 * abs(alpha[0]), "ψ is an eigenvector: one α in every group")
        return OperatorPencil(prompt, T), psi.reshape(-1, 1), float(alpha[0])

    @pytest.mark.catches("ERR-089")
    def test_law_the_alpha_map_solves_the_alpha_equation(self) -> None:
        pencil, psi, alpha = self._alpha_pencil()
        posing = EigenPosing(pencil, ALPHA_MAP)
        _require(alpha > 0.0, "mixture A 2g is supercritical (k_inf = 1.875), so the fundamental α is positive")
        r = np.asarray(posing.residual(psi, alpha)); scale = np.max(np.abs(np.asarray(pencil.lhs.apply(psi))))
        _require(np.max(np.abs(r)) <= 1e-12 * scale, "(A − μ(α) T) ψ must vanish at the physical α")
        _require(abs(posing.rayleigh(psi) - alpha) <= 1e-12 * abs(alpha), "the Rayleigh quotient must be the physical α")

    def test_law_the_k_map_solves_the_k_equation(self) -> None:
        _, pencil, k, phi = _exact_pair()
        A = np.asarray(pencil.lhs.as_matrix(), dtype=float); F = np.asarray(pencil.rhs.as_matrix(), dtype=float)
        k_def = float((F @ phi.ravel()).sum() / (A @ phi.ravel()).sum())  # A φ = F φ / k
        _require(abs(K_MAP(K_MAP.inverse(k_def)) - k_def) <= 1e-14 * k_def, "the round trip")
        _require(abs(EigenPosing(pencil, K_MAP).rayleigh(phi) - k_def) <= 1e-12 * k_def, "the Rayleigh quotient is k")
