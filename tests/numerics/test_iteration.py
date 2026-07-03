"""Foundation + L1 tests for ``orpheus.numerics.iteration``.

Wave E Round 1 (Issue #163) ships :class:`SourceIteration` and
:class:`KEigenvalue` as stand-alone iteration primitives that consume
the Wave A :class:`LinearOperator` Protocol triple :math:`(L, S, F)`.

Tests in this file:

* **Foundation (synthetic):** L0 dense-matrix fixtures where the
  ground truth is :func:`numpy.linalg.solve` /
  :func:`numpy.linalg.eig`.  Pin the algorithmic correctness of the
  primitives in isolation from any transport solver.
* **Foundation (apply-guards):** the constructors raise
  :class:`TypeError` when their argument operators lack the
  required Protocol surface.
* **L1 (SN integration gate):** build an actual SN operator triple
  (:class:`~orpheus.sn.operators.streaming.InvertibleOperator`
  (= ``A = L + C``) /
  :class:`ScatteringOperator` / :class:`FissionOperator`) for a
  2-group homogeneous slab and assert that :class:`KEigenvalue`
  recovers the same :math:`k_{\\rm eff}` as :func:`solve_sn`.  This
  is the gate test that the new primitives compose with the existing
  SN operator algebra.
"""

from __future__ import annotations

import warnings

import numpy as np
import pytest

from orpheus.derivations.common.xs_library import get_mixture
from orpheus.geometry import Mesh1D
from orpheus.numerics.iteration import (
    KEigenvalue,
    KrylovAcceleration,
    SourceIteration,
)
from orpheus.numerics.operator import (
    InverseOperator,
    LinearOperator,
    NotInvertible,
    ZeroOperator,
)
from orpheus.transport.fields.angular_boundary_flux import AngularBoundaryFlux


# ───────────────────────────────────────────────────────────────────────
# Synthetic fixtures (L0 ground truth)
# ───────────────────────────────────────────────────────────────────────


class MatrixOperator(LinearOperator):
    """Test operator backed by a dense numpy matrix.

    Same shape as the fixture in ``test_operator.py`` — kept independent
    here so tests in this file are self-contained.
    """

    def __init__(
        self,
        matrix: np.ndarray,
        *,
        can_solve: bool = False,
        can_transpose: bool = False,
    ) -> None:
        self.matrix = np.asarray(matrix, dtype=float)
        self._can_solve = bool(can_solve)
        self._can_transpose = bool(can_transpose)

    def apply(self, x: np.ndarray) -> np.ndarray:
        return self.matrix @ x

    def solve(self, b: np.ndarray) -> np.ndarray:
        return np.linalg.solve(self.matrix, b)

    def apply_transpose(self, x: np.ndarray) -> np.ndarray:
        return self.matrix.T @ x

    @property
    def is_invertible(self) -> bool:
        return self._can_solve

    @property
    def is_adjointable(self) -> bool:
        return self._can_transpose

    def inverse(self) -> InverseOperator:
        # The #226 step-3 driver contract: the caller builds the inverse
        # OPERATOR and SourceIteration applies it.  The generic
        # InverseOperator delegates apply → this leaf's solve
        # (bit-identical to the pre-step-3 ``L.solve`` step).
        return InverseOperator(self)


@pytest.fixture
def rng() -> np.random.Generator:
    return np.random.default_rng(20260509)


# ───────────────────────────────────────────────────────────────────────
# Foundation: synthetic SourceIteration
# ───────────────────────────────────────────────────────────────────────


@pytest.mark.foundation
def test_source_iteration_recovers_direct_solve(rng):
    """L0 ground truth: SourceIteration on (A − S) matches np.linalg.solve.

    Build a 4×4 SPD-ish ``A`` and a contraction ``S`` (spectral radius
    well below 1).  Solve ``(A − S)·ψ = q`` two ways:

    * Directly: ``np.linalg.solve(A − S, q)``.
    * By SourceIteration with ``F = ZeroOperator()``.

    The fixed-point iteration converges geometrically at rate
    :math:`\\rho(A^{-1}\\,S)`; the two answers must agree to
    1e-10 absolute.
    """
    n = 4
    # Diagonal-dominant matrix → well-conditioned solve.
    A_mat = np.eye(n) * 4.0 + 0.1 * rng.standard_normal((n, n))
    A_mat = 0.5 * (A_mat + A_mat.T) + n * np.eye(n)  # symmetric, dominant
    # Small contraction: spectral radius of A^{-1} S < 0.1.
    S_mat = 0.05 * rng.standard_normal((n, n))

    A = MatrixOperator(A_mat, can_solve=True)
    S = MatrixOperator(S_mat)
    F = ZeroOperator()

    q = rng.standard_normal(n)
    expected = np.linalg.solve(A_mat - S_mat, q)

    si = SourceIteration(A.inverse(), S, F, max_iter=1000, tol=1e-14)
    psi, residuals = si.solve(q)

    np.testing.assert_allclose(psi, expected, atol=1e-10, rtol=1e-10)
    # Residual history must be monotonically (or near-monotonically)
    # decreasing — a contraction map produces geometric decay.
    assert residuals[-1] < 1e-10, (
        f"Residual did not reach 1e-10; final={residuals[-1]:.2e}"
    )


@pytest.mark.foundation
def test_source_iteration_with_fission_term(rng):
    """SourceIteration on full (A − S − F) recovers direct solve.

    Same construction as above but with a small ``F`` term.  Solves
    ``(A − S − F)·ψ = q`` by fixed-point iteration; compares to
    ``np.linalg.solve(A − S − F, q)``.
    """
    n = 4
    A_mat = np.eye(n) * 5.0 + 0.05 * rng.standard_normal((n, n))
    A_mat = 0.5 * (A_mat + A_mat.T) + n * np.eye(n)
    S_mat = 0.05 * rng.standard_normal((n, n))
    F_mat = 0.05 * rng.standard_normal((n, n))

    A = MatrixOperator(A_mat, can_solve=True)
    S = MatrixOperator(S_mat)
    F = MatrixOperator(F_mat)

    q = rng.standard_normal(n)
    expected = np.linalg.solve(A_mat - S_mat - F_mat, q)

    si = SourceIteration(A.inverse(), S, F, max_iter=2000, tol=1e-14)
    psi, _ = si.solve(q)

    np.testing.assert_allclose(psi, expected, atol=1e-10, rtol=1e-10)


@pytest.mark.foundation
def test_source_iteration_with_explicit_solve_realisation():
    r"""The caller controls the inverse step by BUILDING the inverse operator.

    #226 taxonomy step 3 (superseding the R-1 Step B ``solve`` contract):
    the solver layer builds ``A.inverse()`` — whose ``apply`` delegates to
    the leaf's ``solve``, bit-identical — and :class:`SourceIteration`
    only APPLIES it.  This test exercises the pattern with a dense matrix
    wrapped in a ``MatrixOperator``: the caller's chosen inverse action
    (``np.linalg.solve`` under the hood) reaches the driver as the
    inverse-application operator, not as a duck-typed ``.solve`` surface.
    """
    n = 3
    A_mat = np.diag([5.0, 6.0, 7.0])
    S_mat = 0.1 * np.array([[0.0, 1.0, 0.0],
                            [1.0, 0.0, 1.0],
                            [0.0, 1.0, 0.0]])

    # The caller-controlled inverse: the leaf carries the solve
    # realisation; .inverse() lifts it to the operator the driver applies.
    A = MatrixOperator(A_mat, can_solve=True)
    S = MatrixOperator(S_mat)
    F = ZeroOperator()

    q = np.array([1.0, 2.0, 3.0])
    expected = np.linalg.solve(A_mat - S_mat, q)

    si = SourceIteration(A.inverse(), S, F, max_iter=500, tol=1e-14)
    psi, _ = si.solve(q)
    np.testing.assert_allclose(psi, expected, atol=1e-12)


# ───────────────────────────────────────────────────────────────────────
# Foundation: synthetic KrylovAcceleration (R-1 sibling of SourceIteration)
# ───────────────────────────────────────────────────────────────────────


@pytest.mark.foundation
def test_krylov_acceleration_recovers_direct_solve(rng):
    """L0 ground truth: KrylovAcceleration on (A − S) matches np.linalg.solve.

    Same algebraic setup as
    :func:`test_source_iteration_recovers_direct_solve`.  GMRES on
    the composed (A − S) matvec, with ``A.inverse().apply`` as the
    default preconditioner.  Convergence to 1e-10 should take far fewer
    matvecs than source iteration because A − S is well-conditioned.
    """
    n = 4
    A_mat = np.eye(n) * 4.0 + 0.1 * rng.standard_normal((n, n))
    A_mat = 0.5 * (A_mat + A_mat.T) + n * np.eye(n)
    S_mat = 0.05 * rng.standard_normal((n, n))

    A = MatrixOperator(A_mat, can_solve=True)
    S = MatrixOperator(S_mat)
    F = ZeroOperator()

    q = rng.standard_normal(n)
    expected = np.linalg.solve(A_mat - S_mat, q)

    krylov = KrylovAcceleration(A, S, F, max_iter=200, tol=1e-12)
    psi, residuals = krylov.solve(q)

    np.testing.assert_allclose(psi, expected, atol=1e-10, rtol=1e-10)
    assert residuals, "GMRES callback never fired"
    assert residuals[-1] < 1e-8, (
        f"GMRES residual did not reach 1e-8; final={residuals[-1]:.2e}"
    )


@pytest.mark.foundation
def test_krylov_acceleration_with_fission_term(rng):
    """KrylovAcceleration on full (A − S − F) recovers direct solve."""
    n = 4
    A_mat = np.eye(n) * 5.0 + 0.05 * rng.standard_normal((n, n))
    A_mat = 0.5 * (A_mat + A_mat.T) + n * np.eye(n)
    S_mat = 0.05 * rng.standard_normal((n, n))
    F_mat = 0.05 * rng.standard_normal((n, n))

    A = MatrixOperator(A_mat, can_solve=True)
    S = MatrixOperator(S_mat)
    F = MatrixOperator(F_mat)

    q = rng.standard_normal(n)
    expected = np.linalg.solve(A_mat - S_mat - F_mat, q)

    krylov = KrylovAcceleration(A, S, F, max_iter=200, tol=1e-12)
    psi, _ = krylov.solve(q)

    np.testing.assert_allclose(psi, expected, atol=1e-10, rtol=1e-10)


@pytest.mark.foundation
def test_krylov_acceleration_explicit_preconditioner():
    """Caller-supplied ``preconditioner`` shadows the default inverse choice.

    R-1 Step B (2026-05-19) — the parameter name is ``preconditioner``
    (not ``inverter``).  Pass an ``A`` that is NOT invertible and
    supply ``preconditioner`` — construction must succeed and GMRES
    must converge using the supplied preconditioner.
    """
    n = 3
    A_mat = np.diag([5.0, 6.0, 7.0])
    S_mat = 0.1 * np.array([[0.0, 1.0, 0.0],
                            [1.0, 0.0, 1.0],
                            [0.0, 1.0, 0.0]])

    A = MatrixOperator(A_mat, can_solve=False)
    S = MatrixOperator(S_mat)
    F = ZeroOperator()

    inv_A = np.linalg.inv(A_mat)
    preconditioner = lambda q: inv_A @ q

    q = np.array([1.0, 2.0, 3.0])
    expected = np.linalg.solve(A_mat - S_mat, q)

    krylov = KrylovAcceleration(
        A, S, F, preconditioner=preconditioner, max_iter=100, tol=1e-12,
    )
    psi, _ = krylov.solve(q)
    np.testing.assert_allclose(psi, expected, atol=1e-10)


@pytest.mark.foundation
def test_krylov_acceleration_works_without_preconditioner():
    """KrylovAcceleration runs unpreconditioned when A is not invertible.

    No ``preconditioner`` supplied, ``A`` not invertible — GMRES
    still converges, just with more iterations (M = I, the identity
    preconditioner).
    """
    n = 5
    # Well-conditioned diagonal-dominant A so unpreconditioned GMRES
    # still converges quickly.
    A_mat = np.eye(n) * 10.0
    A = MatrixOperator(A_mat, can_solve=False)
    S = ZeroOperator()
    F = ZeroOperator()

    q = np.arange(1.0, n + 1.0)
    expected = q / 10.0

    krylov = KrylovAcceleration(A, S, F, max_iter=50, tol=1e-12)
    assert krylov._preconditioner is None, (
        "Expected no preconditioner when A is not invertible and no "
        "preconditioner is supplied."
    )
    psi, _ = krylov.solve(q)
    np.testing.assert_allclose(psi, expected, atol=1e-10)


@pytest.mark.foundation
def test_krylov_acceleration_high_scattering_beats_source_iteration():
    """At c → 1, GMRES converges in many fewer matvecs than SI.

    The whole point of the KrylovAcceleration sibling is the
    spectral-radius win when the scattering ratio is high.  Pin the
    qualitative win at c ≈ 0.9 (~ρ(A⁻¹S) ≈ 0.9 ⇒ SI needs
    ~log(tol)/log(0.9) ≈ 220 iterations to reach 1e-10 vs GMRES at
    well under that).
    """
    n = 8
    # Diagonal A; nearly-uniform S to make ρ(A⁻¹·S) ≈ 0.9.
    A_mat = np.eye(n) * 1.0
    # All entries = 0.9/n so the row-sum is 0.9 and A⁻¹·S has spectral
    # radius exactly 0.9.
    S_mat = np.full((n, n), 0.9 / n)
    A = MatrixOperator(A_mat, can_solve=True)
    S = MatrixOperator(S_mat)
    F = ZeroOperator()
    q = np.arange(1.0, n + 1.0)

    si = SourceIteration(A.inverse(), S, F, max_iter=500, tol=1e-10)
    _, si_residuals = si.solve(q)

    krylov = KrylovAcceleration(A, S, F, max_iter=500, tol=1e-10)
    _, kr_residuals = krylov.solve(q)

    # GMRES should converge in well under SI's iteration count.  The
    # exact ratio is problem-dependent; pin the qualitative gap at 5×.
    assert len(kr_residuals) < len(si_residuals) / 5, (
        f"KrylovAcceleration ({len(kr_residuals)} iters) was not "
        f"meaningfully faster than SourceIteration "
        f"({len(si_residuals)} iters) at c=0.9 — the algorithmic win "
        f"that motivates the sibling primitive is missing."
    )


@pytest.mark.foundation
def test_krylov_acceleration_requires_apply_on_A():
    class BrokenA:
        pass  # genuinely no apply — the eager guard rejects

    with pytest.raises(TypeError, match=r"requires 'apply' on A"):
        KrylovAcceleration(BrokenA(), ZeroOperator(), ZeroOperator())


@pytest.mark.foundation
def test_krylov_acceleration_requires_apply_on_first_coupling():
    """A coupling gain without apply → TypeError at construction.

    Wave O #208 O.2a: the drivers take the variadic ``(A, *gains)`` shape; the
    per-gain apply check names the offending gain by index (the legacy
    ``S``/``F`` named slots are retired).
    """
    class BrokenS:
        pass

    A = MatrixOperator(np.eye(3), can_solve=True)
    with pytest.raises(
        TypeError,
        match=r"requires 'apply' on every coupling operator; gain 0",
    ):
        KrylovAcceleration(A, BrokenS(), ZeroOperator())


@pytest.mark.foundation
def test_krylov_acceleration_requires_apply_on_later_coupling():
    """A broken gain at a non-zero index is caught and named by its index."""
    class BrokenF:
        pass

    A = MatrixOperator(np.eye(3), can_solve=True)
    with pytest.raises(
        TypeError,
        match=r"requires 'apply' on every coupling operator; gain 1",
    ):
        KrylovAcceleration(A, ZeroOperator(), BrokenF())


# ───────────────────────────────────────────────────────────────────────
# Foundation: synthetic KEigenvalue
# ───────────────────────────────────────────────────────────────────────


@pytest.mark.foundation
def test_keigenvalue_recovers_dominant_eigenvalue(rng):
    """L0 ground truth: KEigenvalue matches numpy.linalg.eig dominant root.

    Build ``A`` (full apply + solve) and ``F`` (apply only).  Solve
    the generalised eigenvalue problem :math:`A\\,\\psi = (1/k)\\,F\\,\\psi`
    two ways:

    * Directly: largest eigenvalue of ``A^{-1}·F`` from
      :func:`numpy.linalg.eig`.
    * By :class:`KEigenvalue` with ``S = ZeroOperator``.

    Power iteration converges to the dominant eigenvalue; the two
    must agree to 1e-9 absolute.
    """
    n = 4
    # Diagonally-dominant L with positive diagonal.
    A_mat = np.diag([3.0, 5.0, 7.0, 11.0]) + 0.05 * rng.standard_normal((n, n))
    A_mat = 0.5 * (A_mat + A_mat.T) + 5.0 * np.eye(n)
    # F: small dominance ratio (k_0 / k_1 ≈ 2 for fast convergence).
    F_mat = np.diag([2.0, 1.0, 0.5, 0.25])

    A = MatrixOperator(A_mat, can_solve=True)
    S = ZeroOperator()
    F = MatrixOperator(F_mat)

    # Reference: largest k = eig(A^{-1} F) (the K = A⁻¹F multiplication
    # operator's dominant root).
    A_inv_F = np.linalg.solve(A_mat, F_mat)
    eigvals = np.linalg.eigvals(A_inv_F)
    expected_keff = float(np.max(np.real(eigvals)))

    initial = np.ones(n)
    ke = KEigenvalue(
        A, S, F,
        max_outer=500, keff_tol=1e-12, flux_tol=1e-12,
        max_inner=500, inner_tol=1e-14,
    )
    keff, keff_history, psi = ke.solve(initial_guess=initial)

    assert abs(keff - expected_keff) < 1e-9, (
        f"KEigenvalue keff={keff!r} expected≈{expected_keff!r}; "
        f"history[-3:]={keff_history[-3:]}"
    )

    # The eigenvector should satisfy (A⁻¹F)·ψ = k·ψ to high precision.
    np.testing.assert_allclose(A_inv_F @ psi, keff * psi, atol=1e-7)


# ───────────────────────────────────────────────────────────────────────
# Foundation: eager apply-guards (carve P4 — TypeError, no registry)
# ───────────────────────────────────────────────────────────────────────


@pytest.mark.foundation
def test_source_iteration_requires_apply_on_A_inv():
    """A_inv without apply → TypeError at construction."""
    class BrokenAInv:
        pass  # genuinely no apply

    A_inv = BrokenAInv()
    S = MatrixOperator(np.eye(2))
    F = ZeroOperator()

    with pytest.raises(TypeError, match="apply"):
        SourceIteration(A_inv, S, F)


@pytest.mark.foundation
def test_source_iteration_requires_apply_on_S():
    """S without apply → TypeError at construction."""
    class BrokenS:
        pass

    A = MatrixOperator(np.eye(2), can_solve=True)
    S = BrokenS()
    F = ZeroOperator()

    with pytest.raises(TypeError, match="apply"):
        SourceIteration(A.inverse(), S, F)


@pytest.mark.foundation
def test_source_iteration_requires_apply_on_F():
    """F without apply → TypeError at construction."""
    class BrokenF:
        pass

    A = MatrixOperator(np.eye(2), can_solve=True)
    S = MatrixOperator(np.eye(2))
    F = BrokenF()

    with pytest.raises(TypeError, match="apply"):
        SourceIteration(A.inverse(), S, F)


@pytest.mark.foundation
def test_invertibility_obligation_lives_at_the_inverse_builder():
    """The R-1 Step B "L must solve" gate MIGRATED to the builder (#226 step 3).

    :class:`SourceIteration` no longer demands ``CAP_SOLVE`` — its step
    operator arrives pre-inverted, and an APPLY-ONLY step operator is
    acceptable BY DESIGN (the windowed product ``P @ A.inverse()`` is
    exactly that shape: no round-trip promise, just the family's canonical
    seeded-apply signature).  The "can this be inverted?" obligation now
    discharges where the inverse is BUILT: ``.inverse()`` on a
    non-invertible leaf raises with the domain message.
    """
    # The obligation fires at the builder …
    with pytest.raises(NotInvertible, match="invertible"):
        MatrixOperator(np.eye(2), can_solve=False).inverse()

    # … and the driver runs an apply-only, seeded-signature step operator
    # END-TO-END (the windowed-product shape): apply-only, the
    # inverse action baked into ``apply``.  Zero gains → one exact step.
    class ApplyOnlyStep:
        def apply(self, rhs, *, initial_guess=None):
            return rhs / 2.0  # the exact inverse of L = 2·I

    si = SourceIteration(ApplyOnlyStep(), max_iter=5, tol=1e-14)
    psi, _ = si.solve(np.array([2.0, 4.0]))
    np.testing.assert_allclose(psi, np.array([1.0, 2.0]))


@pytest.mark.foundation
def test_keigenvalue_requires_invertible_A():
    """KEigenvalue (the posing layer that BUILDS A.inverse()) guards
    invertibility at construction with a domain message — not an
    AttributeError from a missing ``.inverse`` (#226 step 3)."""
    A = MatrixOperator(np.eye(2), can_solve=False)
    S = MatrixOperator(np.eye(2))
    F = MatrixOperator(np.eye(2))

    with pytest.raises(NotInvertible, match="INVERTIBLE"):
        KEigenvalue(A, S, F)


@pytest.mark.foundation
def test_keigenvalue_rejects_non_power_method():
    """eigenvalue_method != 'power' raises NotImplementedError."""
    A = MatrixOperator(np.eye(2), can_solve=True)
    S = ZeroOperator()
    F = MatrixOperator(np.eye(2))

    with pytest.raises(NotImplementedError, match="FEAST"):
        KEigenvalue(A, S, F, eigenvalue_method="feast")


# ───────────────────────────────────────────────────────────────────────
# L1: SN integration gate
# ───────────────────────────────────────────────────────────────────────


@pytest.mark.l1
@pytest.mark.verifies("multigroup")
def test_keigenvalue_matches_solve_sn_2g_slab():
    """L1 gate: KEigenvalue on the SN operator triple matches solve_sn.

    Build a 2-group homogeneous-material 1-D slab, run both:

    * :func:`solve_sn` (the legacy power_iteration-based path), and
    * :class:`KEigenvalue` directly on
      ``(InvertibleOperator (= L + C), ScatteringOperator, FissionOperator)``
      with adapter shims that present scalar-flux shapes consistent
      across operators.

    Assert recovered keff agrees to 1e-9.  Both paths use the SAME
    underlying operators — the only difference is whether the iteration
    primitive is the legacy ``power_iteration(solver)`` (Wave-E pre-
    Round-1) or the new ``KEigenvalue`` (Wave-E Round 1).  Bit-identity
    is not required (the new primitive's keff_estimator may differ
    from the legacy ``compute_keff`` at the last few digits); equality
    to 1e-9 is the operational test that the algebra is consistent.
    """
    # Suppress the eigenvalue.py deprecation warning for this test —
    # we are deliberately exercising the legacy path as a reference.
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", DeprecationWarning)

        from orpheus.transport.operators.fission import FissionOperator
        from orpheus.sn.mesh.augmented_mesh import SNMesh
        from orpheus.numerics.quadrature import Quadrature
        from orpheus.transport.operators.scattering import ScatteringOperator
        from orpheus.sn.solver import SNSolver, solve_sn
        from orpheus.sn.loss_representation import transport_sweep

    # 2-group homogeneous 1-D slab — the same canonical fixture
    # ``test_solver_components.py`` uses for component checks.
    mix = get_mixture("A", "2g")
    materials = {0: mix}
    mesh = Mesh1D(
        edges=np.linspace(0.0, 5.0, 11),
        mat_ids=np.zeros(10, dtype=int),
    )
    quad = Quadrature.gauss_legendre(n_ordinates=8)

    # Reference: solve_sn (legacy power_iteration path).
    ref = solve_sn(
        materials, mesh, quad,
        inner_solver="source_iteration",
        scattering_order=0,
        max_outer=500,
        keff_tol=1e-9,
        flux_tol=1e-8,
        max_inner=500,
        inner_tol=1e-10,
    )
    expected_keff = ref.keff

    # Build the SN operator triple from the same precomputed solver
    # data used for the reference run.  ``A_inv_adapter`` (defined
    # below) wraps :func:`transport_sweep` directly; the
    # :class:`InvertibleOperator` (= ``L + C``) on the SNSolver is
    # unused here.  Solver instance retained to provide
    # ``solver.scattering_op`` / ``solver.fission_op`` / ``solver.mat_xs``.
    sn_mesh = SNMesh(mesh, quad, materials)
    solver = SNSolver(sn_mesh, scattering_order=0)
    # The canonical S, F operators built directly from solver state.
    S = solver.scattering_op
    F = solver.fission_op

    # ── Adapter shims to keep the iteration primitive scalar-flux-only.
    #
    # The SN operator triple has an internal shape-mismatch for
    # historical reasons (Wave D Round 3 docstring §"Vector layout"):
    #   • F.apply takes scalar phi (nx,ny,ng), returns scalar Q (nx,ny,ng)
    #   • S.apply takes angular psi (N,nx,ny,ng), returns angular Q (N,nx,ny,ng)
    #   • L.solve takes Q+psi_bc+Q_aniso, returns (psi, phi) tuple.
    # Round 2 normalises this; for the L1 gate test we wrap each
    # operator into a thin scalar-in/scalar-out facade.
    # Issue #197 PR-TYPED-2 — typed AngularBoundaryFlux replaces psi_bc: dict.
    boundary_flux = AngularBoundaryFlux.zeros_on(sn_mesh)

    class A_inv_adapter(LinearOperator):
        """Adapter: rhs (ng, nx, ny) → phi via the unified sweep.

        Issue #196 PR-INDEX-5: principled layout throughout.  Returns
        scalar flux (drops angular flux for shape consistency with
        F.apply / S.apply scalar facade).

        #226 step 3: carries the invertibility pair (``is_invertible`` +
        ``inverse()``) so :class:`KEigenvalue` — which now BUILDS the
        inverse and hands it to the inner driver — can lift this leaf's
        ``solve`` through the generic :class:`InverseOperator`.
        """

        @property
        def is_invertible(self) -> bool:
            return True

        def inverse(self) -> InverseOperator:
            return InverseOperator(self)

        def apply(self, phi):  # not used by the iteration primitive
            return phi

        def solve(self, rhs):
            # R-1 Step 4 A1: single per-ordinate source carrier.
            # ``rhs`` is bare ndarray (ng, nx, ny) — wrap via the
            # canonical iso → per-ord factory at the adapter boundary.
            from orpheus.transport.source_sinks import AngularSourceSink
            from orpheus.sn.solver import _reflect_outflow_into_inflow
            source = AngularSourceSink.from_isotropic(rhs, sn_mesh)
            # Wave O (#208) O.4a.2 — the bare ``transport_sweep`` no longer
            # re-applies the reflective BC at entry; drive the −B coupling
            # explicitly (reflect the persisted outflow — ``boundary_flux``
            # is the closure-scoped partner-flux carrier — into the inflow
            # slots) before each sweep, exactly as ``_solve_fixed_source_si``
            # does in production.
            _reflect_outflow_into_inflow(boundary_flux, sn_mesh)
            _angular, scalar = transport_sweep(
                source, solver.mat_xs.total_cross_section, sn_mesh,
                boundary_flux,
            )
            return scalar

    class S_scalar_adapter(LinearOperator):
        """Adapter: phi (ng, nx, ny) → P0 scattering source (ng, nx, ny).

        Issue #196 PR-INDEX-5: principled end-to-end.
        """

        def apply(self, phi):
            Q = np.zeros_like(phi)
            S.add_iso_source(Q, phi)
            S.add_n2n_source(Q, phi)
            return Q

    class F_scalar_adapter(LinearOperator):
        """Adapter: phi (ng, nx, ny) → fission source (ng, nx, ny).

        Issue #196 PR-INDEX-5: principled end-to-end.
        """

        def apply(self, phi):
            return F.apply(phi)

    A_adapt = A_inv_adapter()
    S_adapt = S_scalar_adapter()
    F_adapt = F_scalar_adapter()

    # SN-aware keff estimator — matches SNSolver.compute_keff so the
    # outer-iteration math is the same as the legacy reference path.
    def sn_keff_estimator(A_, S_, F_, phi):
        return solver.compute_keff(phi)

    ng = solver.ng
    # Issue #196 PR-INDEX-5: principled initial guess.
    initial = np.ones((ng, *sn_mesh.spatial_shape))

    ke = KEigenvalue(
        A_adapt, S_adapt, F_adapt,
        max_outer=500, keff_tol=1e-9, flux_tol=1e-8,
        max_inner=500, inner_tol=1e-10,
        keff_estimator=sn_keff_estimator,
    )
    keff, keff_history, _phi = ke.solve(initial_guess=initial)

    assert abs(keff - expected_keff) < 1e-9, (
        f"KEigenvalue keff={keff!r} differs from solve_sn "
        f"keff={expected_keff!r} by {abs(keff-expected_keff):.2e}; "
        f"history[-3:]={keff_history[-3:]}"
    )
