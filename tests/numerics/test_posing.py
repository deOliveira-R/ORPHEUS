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
from orpheus.numerics.posing import K_MAP, EigenPosing, SourcePosing

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
    i = int(np.argmax(w.real)); phi = np.real(V[:, i]); phi = phi / phi.sum()
    # the 0-D operators are bound on the (ng, 1) energy × point-spatial carrier
    return problem, pencil, float(w[i].real), phi.reshape(-1, 1)


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
        _require(np.max(np.abs(np.asarray(posing.residual(psi)))) <= 1e-12 * np.max(np.abs(q)), "q − Aψ* vanishes")
        _require(abs(posing.balance(psi)) <= 1e-12 * np.sum(np.abs(q)), "⟨1, Aψ* − q⟩ vanishes")
        _require(posing.balance(2.0 * psi) > 0, "the balance is SIGNED: too much flux reads positive (Aψ − q > 0)")
