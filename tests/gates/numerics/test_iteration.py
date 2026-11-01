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
  (:class:`~orpheus.sn.operators.streaming.StreamingCollisionOperator`
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

from orpheus.numerics.green_operator import GreenOperator
from orpheus.numerics.iteration import (
    KEigenvalue,
    KrylovAcceleration,
    SourceIteration,
)
from orpheus.numerics.pencil import OperatorPencil
from orpheus.numerics.posing import K_MAP, EigenPosing
from orpheus.numerics.operator import ZeroOperator
from orpheus.numerics.operator import (
    IdentityOperator,
    InverseOperator,
    LinearOperator,
    NotInvertible,
)


# ───────────────────────────────────────────────────────────────────────
# Synthetic fixtures (L0 ground truth)
# ───────────────────────────────────────────────────────────────────────


class MatrixOperator(LinearOperator):
    """Test operator backed by a dense numpy matrix.

    Same shape as the fixture in ``test_operator.py`` — kept independent
    here so tests in this file are self-contained.
    """
    # S4-amendment: the base DEMANDS an answer from every subclass; this
    # double is a deliberately-unbound probe, so it DECLARES the unbound
    # state instead of inheriting a silent default (which no longer exists).
    domain = None
    codomain = None

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
    psi, record = si.solve(q)
    residuals = _trajectory(record)

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
    the composed (A − S) matvec, preconditioned by ``A.inverse()`` (stated:
    the preconditioner has no default).  Convergence to 1e-10 should take far fewer
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

    krylov = KrylovAcceleration(
        A, S, F, preconditioner=A.inverse(), max_iter=200, tol=1e-12,
    )
    psi, record = krylov.solve(q)
    residuals = _trajectory(record)

    np.testing.assert_allclose(psi, expected, atol=1e-10, rtol=1e-10)
    assert record.iterated, "GMRES callback never fired"
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

    krylov = KrylovAcceleration(
        A, S, F, preconditioner=A.inverse(), max_iter=200, tol=1e-12,
    )
    psi, _ = krylov.solve(q)

    np.testing.assert_allclose(psi, expected, atol=1e-10, rtol=1e-10)


@pytest.mark.foundation
def test_krylov_acceleration_explicit_preconditioner():
    """The preconditioner is an operator independent of ``A``'s invertibility.

    R-1 Step B (2026-05-19) — the parameter name is ``preconditioner``
    (not ``inverter``).  Pass an ``A`` that is NOT invertible and
    supply a preconditioner OPERATOR — construction must succeed and GMRES
    must converge using it.
    """
    n = 3
    A_mat = np.diag([5.0, 6.0, 7.0])
    S_mat = 0.1 * np.array([[0.0, 1.0, 0.0],
                            [1.0, 0.0, 1.0],
                            [0.0, 1.0, 0.0]])

    A = MatrixOperator(A_mat, can_solve=False)
    S = MatrixOperator(S_mat)
    F = ZeroOperator()

    preconditioner = MatrixOperator(np.linalg.inv(A_mat))

    q = np.array([1.0, 2.0, 3.0])
    expected = np.linalg.solve(A_mat - S_mat, q)

    krylov = KrylovAcceleration(
        A, S, F, preconditioner=preconditioner, max_iter=100, tol=1e-12,
    )
    psi, _ = krylov.solve(q)
    np.testing.assert_allclose(psi, expected, atol=1e-10)


@pytest.mark.foundation
def test_krylov_acceleration_works_without_preconditioner():
    """Plain GMRES is the STATED identity preconditioner, ``M = I``.

    ``A`` not invertible and ``preconditioner=IdentityOperator()``: GMRES
    still converges, just with more iterations.  Until 2026-10-04 the same
    run was reached by OMITTING the preconditioner, which silently chose
    ``A⁻¹`` when ``A`` was invertible and the identity when it was not;
    the choice is now always written at the call site (the user's ruling).
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

    identity = IdentityOperator()
    krylov = KrylovAcceleration(
        A, S, F, preconditioner=identity, max_iter=50, tol=1e-12,
    )
    assert krylov.preconditioner is identity
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
    _, si_record = si.solve(q)

    krylov = KrylovAcceleration(
        A, S, F, preconditioner=A.inverse(), max_iter=500, tol=1e-10,
    )
    _, kr_record = krylov.solve(q)

    # GMRES should converge in well under SI's iteration count.  The
    # exact ratio is problem-dependent; pin the qualitative gap at 5×.
    # ⭐ Counts, not trajectory lengths.  The two drivers use DIFFERENT
    # conventions (SI measures differences, so its pass count exceeds its
    # trajectory by one; GMRES gets one callback per iteration), and this
    # comparison was silently mixing them until each driver began stating
    # its own count (#340 F11).
    assert kr_record.n_iterations < si_record.n_iterations / 5, (
        f"KrylovAcceleration ({kr_record.n_iterations} iters) was not "
        f"meaningfully faster than SourceIteration "
        f"({si_record.n_iterations} iters) at c=0.9 — the algorithmic win "
        f"that motivates the sibling primitive is missing."
    )


@pytest.mark.foundation
def test_krylov_acceleration_requires_apply_on_A():
    class BrokenA:
        pass  # genuinely no apply — the eager guard rejects

    with pytest.raises(TypeError, match=r"requires 'apply' on A"):
        KrylovAcceleration(
            BrokenA(), ZeroOperator(), ZeroOperator(),
            preconditioner=IdentityOperator(),
        )


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
        KrylovAcceleration(
            A, BrokenS(), ZeroOperator(), preconditioner=IdentityOperator(),
        )


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
        KrylovAcceleration(
            A, ZeroOperator(), BrokenF(), preconditioner=IdentityOperator(),
        )


@pytest.mark.foundation
def test_krylov_acceleration_requires_a_stated_preconditioner():
    """The preconditioner has no default: omitting it is refused.

    The user's ruling (2026-10-04): no default is right for every model
    (``A⁻¹`` is SN's sweep but diffusion's whole solve), and a silent
    identity default hands the next model's Krylov path iteration counts
    that grow with the mesh, which is what issue #200 removed from SN.
    Plain GMRES is spelled ``preconditioner=IdentityOperator()``.
    """
    A = MatrixOperator(np.eye(3), can_solve=True)
    with pytest.raises(TypeError, match=r"preconditioner"):
        KrylovAcceleration(A, ZeroOperator())  # type: ignore[call-arg]


@pytest.mark.foundation
def test_krylov_acceleration_requires_apply_on_the_preconditioner():
    """A bare callable (the contract until 2026-10-04) is refused at construction.

    The preconditioner is a ``LinearOperator``: GMRES applies it to typed
    fields, and a callable typed ``ndarray -> ndarray`` misstated that.
    """
    A = MatrixOperator(np.eye(3), can_solve=True)
    with pytest.raises(TypeError, match=r"requires 'apply' on the preconditioner"):
        KrylovAcceleration(
            A, ZeroOperator(), preconditioner=lambda q: q,  # type: ignore[arg-type]
        )


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
        EigenPosing(OperatorPencil(A - S, F), K_MAP), GreenOperator(A, S),
        max_outer=500, keff_tol=1e-12, flux_tol=1e-12,
        max_inner=500, inner_tol=1e-14,
    )
    _o = ke.solve(initial_guess=initial)
    keff, keff_history, psi = _o.keff, _o.keff_history, _o.flux_distribution

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
    """KEigenvalue's invertibility obligation is its Green's (S4, ruling 68):
    the Green the Strategy hands it refuses a non-invertible implicit part at
    construction with a domain message, before any KEigenvalue exists — not an
    AttributeError from a missing ``.inverse`` (#226 step 3). Until S4
    KEigenvalue held the split itself and raised its own ``INVERTIBLE``."""
    A = MatrixOperator(np.eye(2), can_solve=False)
    S = MatrixOperator(np.eye(2))

    with pytest.raises(NotInvertible, match="invertible implicit part M"):
        GreenOperator(A, S)


@pytest.mark.foundation
def test_keigenvalue_rejects_non_power_method():
    """eigenvalue_method != 'power' raises NotImplementedError."""
    A = MatrixOperator(np.eye(2), can_solve=True)
    S = ZeroOperator()
    F = MatrixOperator(np.eye(2))

    with pytest.raises(NotImplementedError, match="FEAST"):
        KEigenvalue(EigenPosing(OperatorPencil(A - S, F), K_MAP), GreenOperator(A, S), eigenvalue_method="feast")


# ───────────────────────────────────────────────────────────────────────
# L1: SN integration gate
# ───────────────────────────────────────────────────────────────────────




def _trajectory(record) -> list[float]:
    """The single criterion's trajectory — what ``solve`` returned as a bare
    list before #340 N2a gave the drivers an
    :class:`~orpheus.numerics.convergence.IterationRecord`.

    Not a compatibility shim: these rows are ABOUT the residual sequence, so
    naming the extraction keeps them reading that way.  Rows that are about
    the iteration COUNT read ``record.n_iterations`` directly instead —
    the two are not interchangeable across drivers.
    """
    return list(record.criteria[0].trajectory)








# ───────────────────────────────────────────────────────────────────────
# Foundation: the GMRES exact-breakdown carve-out (the ERR-053 warn
# boundary) — B.2d d3 NIT-2 anchor
# ───────────────────────────────────────────────────────────────────────
#
# The warn branch in :meth:`KrylovAcceleration.solve` carries a PERMANENT
# exact-breakdown carve-out: a final preconditioned residual of LITERAL
# ``0.0`` means the Krylov space collapsed AT the solution (``M⁻¹(b − Ax)
# = 0`` with a nonsingular preconditioner implies ``Ax = b`` exactly) —
# that is CONVERGENCE, the opposite of the ERR-053 restart-truncation
# stall the warning surfaces, so it must NOT warn even when scipy stamps
# ``info > 0``.  The guard's first production caller (the B.2d
# transitional dead-ray GMRES padding) dissolved at the d2 eviction,
# leaving it caller-less-but-permanent — these anchors keep both branch
# arms exercised.  scipy's ``info`` stamping on breakdown is
# version-dependent (current scipy stamps ``info == 0`` on every small
# reachable case because its convergence test is ``<=``), so the branch
# pair is pinned through a deterministic gmres STUB while the real-scipy
# arm pins the caller-visible contract.


def _singular_consistent() -> tuple[MatrixOperator, np.ndarray]:
    r"""``A = diag(1, 0)`` with ``b ∈ range(A)`` — the minimal
    singular-but-CONSISTENT system.  GMRES from ``x0 = 0`` solves it
    exactly in one iteration (``K₁ = span{b}`` contains the solution), so
    the ``pr_norm`` history ends at LITERAL ``0.0`` — the exact-breakdown
    signature."""
    return MatrixOperator(np.diag([1.0, 0.0])), np.array([3.0, 0.0])


def _krylov_warnings(krylov: KrylovAcceleration, b: np.ndarray):
    """Run ``krylov.solve(b)`` recording only this boundary's RuntimeWarnings.

    Returns the pr_norm TRAJECTORY (not the record): all three callers are
    ERR-053 gates asserting on the residual tail, which is the trajectory's
    job.  ``solve`` began returning an
    :class:`~orpheus.numerics.convergence.IterationRecord` at #340 N2a.
    """
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        x, record = krylov.solve(b)
        history = _trajectory(record)
    ours = [
        w for w in caught
        if issubclass(w.category, RuntimeWarning)
        and "KrylovAcceleration.solve" in str(w.message)
    ]
    return x, history, ours


@pytest.mark.foundation
def test_singular_consistent_exact_breakdown_solves_clean() -> None:
    r"""REAL scipy: the singular-consistent exact-breakdown solve returns
    the EXACT solution with a literal-``0.0`` residual tail and NO
    ERR-053 warning — the caller-visible contract of the carve-out
    (``tol=0.0`` = "solve to exactness" makes the breakdown the only
    possible exit).

    Carries no ``catches("ERR-053")``: this row guards the carve-out
    against OVER-warning, the opposite direction of the defect, and the
    defect re-dropped (the ``restart=min(50, n_dof)`` clamp, the discarded
    ``info``, and both together) leaves it green (``[M]`` 2026-10-05,
    ``scratch/reference_architecture/p3/gates200b/arm53_*.log``). The
    catchers are :func:`test_info_warning_fires_on_genuine_nonconvergence`
    and :func:`test_the_gmres_nonconvergence_warning_is_ESCALATABLE` (the
    discarded ``info``), and the restart gates of
    ``tests/gates/sn/solve/`` (the clamp)."""
    A_op, b = _singular_consistent()
    krylov = KrylovAcceleration(
        A_op, preconditioner=IdentityOperator(), tol=0.0, max_iter=30, restart=2,
    )
    x, history, ours = _krylov_warnings(krylov, b)
    if not history or history[-1] != 0.0:
        pytest.fail(
            f"the singular-consistent probe no longer ends at a LITERAL 0.0 "
            f"pr_norm (tail = {history[-3:] if history else []}) — the "
            f"exact-breakdown fixture premise broke; redesign the anchor."
        )
    np.testing.assert_allclose(
        A_op.matrix @ x, b, atol=5e-16,
        err_msg="exact breakdown did not return the exact solution",
    )
    if ours:
        pytest.fail(
            f"exact breakdown WARNED ({ours[0].message}) — the ERR-053 "
            f"carve-out regressed: a converged solve is being reported as "
            f"a restart stall."
        )


def _stub_gmres_returning(info: int, pr_norm_tail: float, solution=None):
    """A deterministic ``spla.gmres`` stand-in: feeds the callback one
    ``pr_norm`` value, stamps ``info`` and returns ``solution`` (zeros when
    none is given) — pins the guard BRANCH independent of scipy's
    version-dependent breakdown stamping."""

    def _stub(A, b, x0=None, M=None, rtol=None, atol=None, maxiter=None,
              restart=None, callback=None, callback_type=None):
        if callback is not None:
            callback(pr_norm_tail)
        x = (
            np.zeros_like(np.asarray(b)) if solution is None
            else np.asarray(solution, dtype=float)
        )
        return x, info

    return _stub


@pytest.mark.foundation
def test_exact_breakdown_guard_suppresses_the_info_warning(monkeypatch) -> None:
    r"""GUARD arm: ``info > 0`` WITH a literal-``0.0`` tail (the d1-observed
    breakdown stamping) AND an iterate that solves the system must NOT warn
    — the carve-out recognizes the collapsed-at-the-solution Krylov space as
    convergence.  Since 2026-10-05 the stub returns the true solution
    ``(3, 0)``: the carve-out is confirmed against the true residual, so a
    stub returning a non-solution now warns, correctly
    (:func:`test_breakdown_on_an_inconsistent_system_warns_and_is_not_accepted`).

    Carries no ``catches("ERR-053")``, for the reason given on
    :func:`test_singular_consistent_exact_breakdown_solves_clean`: the
    re-dropped defect leaves it green (``[M]`` 2026-10-05)."""
    import orpheus.numerics.iteration as iteration_mod

    monkeypatch.setattr(
        iteration_mod.spla, "gmres",
        _stub_gmres_returning(info=7, pr_norm_tail=0.0, solution=[3.0, 0.0]),
    )
    A_op, b = _singular_consistent()
    krylov = KrylovAcceleration(
        A_op, preconditioner=IdentityOperator(), tol=1e-12, max_iter=30, restart=2,
    )
    _x, history, ours = _krylov_warnings(krylov, b)
    if history != [0.0]:
        pytest.fail(f"stub plumbing broke: history = {history}")
    if ours:
        pytest.fail(
            f"info=7 with a literal-0.0 tail WARNED ({ours[0].message}) — "
            f"the exact-breakdown guard branch regressed."
        )


@pytest.mark.foundation
@pytest.mark.catches("ERR-053")
def test_info_warning_fires_on_genuine_nonconvergence(monkeypatch) -> None:
    r"""TEETH arm: the SAME stub with a NONZERO tail (a genuine ERR-053
    stall signature) MUST warn — proving the guard arm's silence above is
    the carve-out biting, not dead warn machinery."""
    import orpheus.numerics.iteration as iteration_mod

    monkeypatch.setattr(
        iteration_mod.spla, "gmres",
        _stub_gmres_returning(info=7, pr_norm_tail=0.37),
    )
    A_op, b = _singular_consistent()
    krylov = KrylovAcceleration(
        A_op, preconditioner=IdentityOperator(), tol=1e-12, max_iter=30, restart=2,
    )
    _x, _history, ours = _krylov_warnings(krylov, b)
    if not ours:
        pytest.fail(
            "info=7 with a NONZERO residual tail did not warn — the "
            "ERR-053 non-convergence surface went silent (the guard has "
            "no teeth)."
        )
    if "ERR-053" not in str(ours[0].message):
        pytest.fail(f"the warning lost its ERR-053 pointer: {ours[0].message}")


@pytest.mark.foundation
@pytest.mark.catches("ERR-098")
def test_breakdown_on_an_inconsistent_system_warns_and_is_not_accepted() -> None:
    r"""qa's F2 (2026-10-05): scipy sets the preconditioned residual to a
    literal 0.0 at EVERY Arnoldi breakdown, including on an INCONSISTENT
    singular system, where ``A x = b`` has no solution.  ``A = diag(1, 0)``
    with ``b = (1, 1)``: real scipy GMRES breaks down with a 0.0 tail and
    ``info > 0`` at a true relative residual of 0.71.  Until 2026-10-05 the
    carve-out accepted that tail alone, so the solve was SILENT and its
    record read converged.  The carve-out now confirms the breakdown
    against the true residual.  Reddens on the tail-only carve-out: no
    warning, ``accepted`` True."""
    A_op = MatrixOperator(np.diag([1.0, 0.0]))
    b = np.array([1.0, 1.0])
    krylov = KrylovAcceleration(
        A_op, preconditioner=IdentityOperator(), tol=1e-8, max_iter=5, restart=2,
    )
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        x, record = krylov.solve(b)
    if not record.criteria[0].trajectory or record.criteria[0].trajectory[-1] != 0.0:
        pytest.fail(
            "premise broke: scipy no longer ends this solve on a literal-0.0 "
            f"tail ({record.criteria[0].trajectory}); the row is blind"
        )
    true_relative = float(np.linalg.norm(b - A_op.matrix @ x) / np.linalg.norm(b))
    assert true_relative > 0.5, f"the system became solvable: {true_relative:.2e}"
    assert any(issubclass(w.category, RuntimeWarning) for w in caught), (
        "a breakdown on an inconsistent system was silent: the carve-out "
        "trusted the 0.0 tail without the true residual"
    )
    assert record.accepted is False and not record.converged, (
        f"the record of a refused solve reads accepted={record.accepted}, "
        f"converged={record.converged}"
    )


@pytest.mark.foundation
@pytest.mark.catches("ERR-097")
@pytest.mark.parametrize("seed", range(4))
def test_krylov_record_is_not_converged_when_scipy_refuses_a_cleared_criterion(seed: int) -> None:
    r"""The user's ruling (2026-10-05, ERR-097): a record cannot read
    converged when its producer's own acceptance test failed.  scipy's GMRES
    stops its inner loop on ``‖P r‖ ≤ rtol·‖P b‖`` and accepts on the TRUE
    residual.  With a strongly anisotropic preconditioner
    (``P = diag(logspace(-3, 3))``, qa's F1 fixture, n = 60, one cycle of
    the full size) the recorded criterion clears (about 5e-9 against
    ``tol = 1e-8``) while the true relative residual is 2e-4 to 7e-4 and
    scipy returns ``info = 1`` (``[M]`` 2026-10-05,
    ``scratch/reference_architecture/p3/operator_contract/probe_fold.py``).
    The row asserts that premise (it is blind without it), then that the
    record reads ``accepted=False`` and NOT converged.  Reddens when the
    record ignores scipy's acceptance (``accepted`` left ``None``)."""
    n = 60
    rng = np.random.default_rng(seed)
    A = MatrixOperator(np.eye(n) + 0.5 * rng.standard_normal((n, n)) / np.sqrt(n))
    P = MatrixOperator(np.diag(np.logspace(-3, 3, n)))
    b = rng.standard_normal(n)
    krylov = KrylovAcceleration(A, preconditioner=P, tol=1e-8, max_iter=1, restart=n)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        x, record = krylov.solve(b)
    true_relative = float(np.linalg.norm(b - A.matrix @ x) / np.linalg.norm(b))
    if not (record.criteria[0].cleared and true_relative > 10 * krylov.tol):
        pytest.fail(
            f"premise broke: criterion cleared={record.criteria[0].cleared}, "
            f"true residual {true_relative:.1e}; the row no longer exercises "
            f"a criterion that clears on a refused solve"
        )
    assert record.accepted is False, f"scipy refused, the record says accepted={record.accepted}"
    assert not record.converged, (
        f"the record reads converged on a solve scipy refused (true relative "
        f"residual {true_relative:.1e})"
    )


# ───────────────────────────────────────────────────────────────────────
# Foundation: the inner GMRES record reads scipy's stopping quantity under a
# preconditioner (qa's #200 review, item 3; the user's ruling 2026-10-04)
# ───────────────────────────────────────────────────────────────────────
#
# scipy's left-preconditioned GMRES steers its inner loop on ‖P r_k‖ ≤
# rtol·‖P b‖, reports ‖P r_k‖/‖b‖ to a ``pr_norm`` callback, and accepts the
# solve on the true residual ‖b − A x‖ ≤ rtol·‖b‖ (``gmres`` in
# ``scipy/sparse/linalg/_isolve/iterative.py``). The record's one criterion,
# ``pr_residual``, is the callback reading rescaled by ‖b‖/‖P b‖, so it IS
# the loop's own test, ‖P r_k‖/‖P b‖. Until 2026-10-04 it recorded the raw
# reading, which is relative only when ‖P b‖ = ‖b‖: true of the identity and
# false of the SN sweep (‖P b‖/‖b‖ from 0.04 to 21.7 on qa's slabs).
#
# The fixture: a 40-unknown non-symmetric tridiagonal with a varying diagonal
# d, preconditioned by s·diag(1/d). Scaling P by a scalar s leaves every
# GMRES iterate unchanged in exact arithmetic (the Arnoldi basis is
# normalised; the Hessenberg entries, the rotated right-hand side and the
# stopping threshold all scale by s), and for s a power of two the scaling is
# exact in binary floating point, so the iterates are bit-identical. The two
# scales 1/16 and 16 put ‖P b‖/‖b‖ at 0.020 and 5.1 (``[M]`` 2026-10-05,
# ``scratch/reference_architecture/p3/gates200b/probe_l0.py``): far from 1 in
# both directions, which is where the raw reading and the ruled one part.

_PRECONDITIONER_SCALES = (1.0 / 16.0, 16.0)


def _preconditioned_tridiagonal(
    scale: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """``(A, M, b, x0)``: the non-symmetric tridiagonal, ``scale`` times its
    Jacobi inverse, a source and a non-zero initial guess (the SN driver
    always warm-starts, so the reference row must too)."""
    rng = np.random.default_rng(5)
    n = 40
    diagonal = 2.0 + rng.uniform(0.0, 3.0, n)
    A = (
        np.diag(diagonal)
        + np.diag(np.full(n - 1, -1.3), -1)
        + np.diag(np.full(n - 1, -0.4), 1)
    )
    b = rng.normal(size=n)
    x0 = 0.1 * rng.normal(size=n)
    return A, scale * np.diag(1.0 / diagonal), b, x0


def _gmres_spy(monkeypatch) -> list[int]:
    """Wrap scipy's ``gmres`` as ``KrylovAcceleration.solve`` reaches it and
    capture every returned ``info`` (the acceptance flag the record is held
    to; the Mode-11 sentinel that the row ran the solve it reads)."""
    import orpheus.numerics.iteration as iteration_mod

    infos: list[int] = []
    honest = iteration_mod.spla.gmres

    def spy(*args, **kwargs):
        x, info = honest(*args, **kwargs)
        infos.append(int(info))
        return x, info

    monkeypatch.setattr(iteration_mod.spla, "gmres", spy)
    return infos


def _solve_preconditioned(
    scale: float, *, tol: float, max_iter: int, restart: int, warm: bool,
):
    A, M, b, x0 = _preconditioned_tridiagonal(scale)
    krylov = KrylovAcceleration(
        MatrixOperator(A), preconditioner=MatrixOperator(M),
        tol=tol, max_iter=max_iter, restart=restart,
    )
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        x, record = krylov.solve(b, initial_guess=x0 if warm else None)
    return (A, M, b), x, record


@pytest.mark.foundation
@pytest.mark.catches("ERR-097")
def test_krylov_record_is_invariant_under_preconditioner_scaling() -> None:
    r"""THEOREM: the recorded ``pr_residual`` trajectory does not depend on a
    scalar scaling of the preconditioner. ``M`` and ``s·M`` drive GMRES
    through the same iterates (exactly, for ``s`` a power of two), and the
    ruled criterion ``‖M r_k‖/‖M b‖`` is homogeneous of degree zero in ``M``,
    so the trajectories are bit-identical: ``[M]`` 21 iterations, ``array_equal``
    at ``s`` = 1/16, 1 and 16. Reddens on the raw callback reading
    ``‖M r_k‖/‖b‖`` (the trajectory scales with ``s``) and on any normaliser
    that is not homogeneous of degree one in ``M``."""
    _system, _x, reference = _solve_preconditioned(
        1.0, tol=1e-8, max_iter=50, restart=40, warm=True,
    )
    expected = _trajectory(reference)
    if not len(expected) > 5:
        pytest.fail(f"the fixture converged in {len(expected)} iterations; "
                    f"a trajectory that short gates nothing")
    for scale in _PRECONDITIONER_SCALES:
        _system, _x, record = _solve_preconditioned(
            scale, tol=1e-8, max_iter=50, restart=40, warm=True,
        )
        np.testing.assert_array_equal(
            np.asarray(_trajectory(record)), np.asarray(expected),
            err_msg=f"s = {scale}: the inner record's pr_residual moved under "
                    f"a scalar scaling of the preconditioner — it is not the "
                    f"preconditioned RELATIVE residual ‖M r‖/‖M b‖",
        )


@pytest.mark.foundation
@pytest.mark.catches("ERR-097")
@pytest.mark.parametrize("scale", _PRECONDITIONER_SCALES)
def test_krylov_record_last_is_the_preconditioned_relative_residual(
    scale: float,
) -> None:
    r"""REFERENCE: the record's last ``pr_residual`` equals
    ``‖M (b − A x)‖ / ‖M b‖`` recomputed by dense numpy at the RETURNED
    iterate ``x`` (a warm-started solve, so the normaliser ``‖M b‖`` is told
    apart from the initial residual ``‖M r_0‖``). Left-preconditioned GMRES
    minimises exactly that norm, so the two agree to the rounding of the
    recomputation, bounded by ``n ε (‖M‖‖A‖‖x‖ + ‖M‖‖b‖)/‖M b‖``: ``[M]``
    3.7e-17 against a bound of 5.2e-14, the value 7.0e-9. Reddens on the raw
    reading (off by ‖M b‖/‖b‖), on the inverted factor, and on a normaliser
    ‖M r_0‖."""
    (A, M, b), x, record = _solve_preconditioned(
        scale, tol=1e-8, max_iter=50, restart=40, warm=True,
    )
    last = _trajectory(record)[-1]
    reference = np.linalg.norm(M @ (b - A @ x)) / np.linalg.norm(M @ b)
    eps = np.finfo(float).eps
    rounding = b.size * eps * np.linalg.norm(M, 2) * (
        np.linalg.norm(A, 2) * np.linalg.norm(x) + np.linalg.norm(b)
    ) / np.linalg.norm(M @ b)
    if not reference < 1e-6:
        pytest.fail(f"s = {scale}: the solve did not converge "
                    f"(‖M r‖/‖M b‖ = {reference:.2e}); the row's premise broke")
    if abs(last - reference) > rounding:
        pytest.fail(
            f"s = {scale}: the record's last pr_residual {last:.6e} is not the "
            f"preconditioned relative residual at the returned iterate "
            f"{reference:.6e} (|Δ| = {abs(last - reference):.2e} > rounding "
            f"{rounding:.2e})"
        )


#: ``(label, scale, tol, max_iter, restart, scipy accepts)``. The ACCEPTED rows
#: converge in one restart cycle; the EXHAUSTED rows stop after three Arnoldi
#: steps (one cycle of restart 3) with ``‖M r‖/‖M b‖`` = 0.08, ``info = 1``.
#: ``[M]`` Under the raw reading the s = 16 accepted row reads 3.4e-8 > tol (the false
#: "not converged" qa measured on the thin slab) and the s = 1/1024 exhausted
#: row reads 2.6e-5 < tol = 1e-3 (the false "converged"); the other two rows
#: are the controls that read the same either way.
_VERDICT_ROWS = (
    ("accepted_s16", 16.0, 1e-8, 50, 40, True),
    ("accepted_s1_16", 1.0 / 16.0, 1e-8, 50, 40, True),
    ("exhausted_s1_1024", 1.0 / 1024.0, 1e-3, 1, 3, False),
    ("exhausted_s1024", 1024.0, 1e-3, 1, 3, False),
)

#: The two rows on which the raw reading is wrong (ERR-097); the other two are
#: its controls and catch nothing.
_ERR_097_VERDICT_ROWS = frozenset({"accepted_s16", "exhausted_s1_1024"})


@pytest.mark.foundation
@pytest.mark.parametrize(
    "label, scale, tol, max_iter, restart, accepted",
    [
        pytest.param(
            *row, id=row[0],
            marks=pytest.mark.catches("ERR-097") if row[0] in _ERR_097_VERDICT_ROWS else (),
        )
        for row in _VERDICT_ROWS
    ],
)
def test_krylov_record_verdict_agrees_with_scipy_acceptance(
    label: str, scale: float, tol: float, max_iter: int, restart: int,
    accepted: bool, monkeypatch,
) -> None:
    r"""The record reads converged exactly when scipy returned ``info = 0``, on
    four rows: two accepted and two budget-exhausted, at preconditioner scales
    on both sides of 1. Honest scope: the two predicates are different
    quantities (the loop's ‖M r‖/‖M b‖ and scipy's true ‖r‖/‖b‖), and they
    agree here because the fixture's residual is not concentrated where M is
    small; a solve that clears the first and not the second is the claim
    check's (``orpheus.sn.solver._check_convergence_claim``), not this
    record's. Reddens on the raw reading at ``accepted_s16`` (false
    not-converged) and at ``exhausted_s1_1024`` (false converged)."""
    infos = _gmres_spy(monkeypatch)
    _system, _x, record = _solve_preconditioned(
        scale, tol=tol, max_iter=max_iter, restart=restart, warm=False,
    )
    if len(infos) != 1:
        pytest.fail(f"{label}: expected one scipy gmres call, saw {len(infos)}")
    if (infos[0] == 0) is not accepted:
        pytest.fail(f"{label}: the fixture's premise moved — scipy returned "
                    f"info = {infos[0]}")
    if record.converged is not accepted:
        pytest.fail(
            f"{label}: scipy returned info = {infos[0]} but the inner record "
            f"reads converged = {record.converged} (last pr_residual "
            f"{_trajectory(record)[-1]:.3e}, tol {tol:.0e})"
        )


# ═══════════════════════════════════════════════════════════════════════
# Step 5 (R-5.2/R-5.3) — the ρ-honest free-identity stop (C1/C5/r5)
# ═══════════════════════════════════════════════════════════════════════
#
# The SI stop is the equation residual via the free identity
# ``r_n = rhs_{n−1} − rhs_n = Σ g_i (ψ_{n−1} − ψ_n)`` normalized by
# ``‖q_ext‖`` — EXACT (= ``Aψ_n − q``) when the step operator is an exact
# inverse of M. These rows pin the identity ELEMENT-WISE against
# explicitly-computed references on a dense toy (Mode-12-safe: the object,
# not just a norm), the zero-gain honest exit, and the q ≈ 0 guard.


class _RecordingInverse:
    """A dense exact inverse that RECORDS every (rhs, ψ) pair — the
    reference sequence the free-identity rows recompute against."""

    def __init__(self, A: np.ndarray) -> None:
        self._inv = np.linalg.inv(A)
        self.rhs_seen: list[np.ndarray] = []
        self.psi_out: list[np.ndarray] = []

    def apply(self, rhs, /, *, initial_guess=None):
        del initial_guess
        psi = self._inv @ np.asarray(rhs, dtype=float)
        self.rhs_seen.append(np.asarray(rhs, dtype=float).copy())
        self.psi_out.append(psi.copy())
        return psi


@pytest.mark.foundation
def test_stop_is_the_free_identity_residual_elementwise(rng):
    """C1 — the recorded history ≡ ‖rhs_{n−1} − rhs_n‖/‖q‖ ≡
    ‖(S+F)·Δψ‖/‖q‖ ≡ ‖Aψ − q‖/‖q‖ (three independently-computed
    spellings; the last is the HONEST residual — exact-M toy).
    Discrimination is BY VALUE: the fixture's ‖ψ‖ and ‖q‖ scales differ
    O(10²), so a ψ-normalized mutant (the historical stop) sits O(10²)
    off every reference; a gain dropped from the rhs bookkeeping breaks
    the (S+F)·Δψ match O(1)."""
    n = 6
    A_mat = np.diag(rng.uniform(2.0, 3.0, n))
    S_mat = 0.35 * rng.random((n, n))
    F_mat = 0.25 * rng.random((n, n))
    A_inv = _RecordingInverse(A_mat)
    S, F = MatrixOperator(S_mat), MatrixOperator(F_mat)
    q = 200.0 * (rng.random(n) + 0.5)  # ‖q‖ scale ≫ tol, ≠ ‖ψ‖ scale
    si = SourceIteration(A_inv, S, F, max_iter=200, tol=1e-10)
    psi, record = si.solve(q)
    history = _trajectory(record)

    if not history or not history[-1] < 1e-10:
        pytest.fail("fixture did not converge — the rows below assume the "
                    "break path")
    q_norm = np.linalg.norm(q)
    # The BREAK pass computes one final rhs that never reaches the inverse
    # (the stop fires before the apply) — reconstruct it from the returned
    # iterate; the recorded applies plus that tail ARE the comparison pairs.
    # Reconstructed with the loop's own accumulation ORDER
    # ((q + S·ψ) + F·ψ — associativity matters at the 1e-12 bar):
    psi_arr = np.asarray(psi)
    rhs_seq = A_inv.rhs_seen + [(q + S_mat @ psi_arr) + F_mat @ psi_arr]
    psi_seq = [np.zeros(n)] + A_inv.psi_out
    if len(history) != len(rhs_seq) - 1:
        pytest.fail(
            f"history length {len(history)} ≠ #rhs-comparisons "
            f"{len(rhs_seq) - 1} — the stop is not the rhs-difference"
        )
    # The residual is a DIFFERENCE of O(‖q‖)-scale vectors, so near
    # convergence every independent re-spelling is cancellation-limited
    # (absolute noise ~ eps·‖rhs‖). The identity rows therefore run on
    # the iterations where r is MEASURABLE (well above that floor); the
    # stop-VALUE row (the same subtraction the production loop performs)
    # is bit-safe and runs on every iteration.
    cancellation_floor = 1e-12 * q_norm
    measurable = 0
    for k, res in enumerate(history):
        r_rhs = rhs_seq[k] - rhs_seq[k + 1]
        np.testing.assert_allclose(
            res, np.linalg.norm(r_rhs) / q_norm, rtol=1e-12,
            err_msg=f"iter {k}: the stop value is not ‖r‖/‖q_ext‖ — "
                    f"a ψ-normalized mutant would sit O(‖ψ‖/‖q‖) off",
        )
        if np.linalg.norm(r_rhs) < 1e4 * cancellation_floor:
            continue
        measurable += 1
        r_gain = (S_mat + F_mat) @ (psi_seq[k] - psi_seq[k + 1])
        np.testing.assert_allclose(
            r_rhs, r_gain, rtol=1e-9, atol=cancellation_floor,
            err_msg=f"iter {k}: rhs-difference ≠ (S+F)·Δψ — a gain "
                    f"dropped from the rhs bookkeeping",
        )
        # The HONEST residual of the iterate the comparison certifies:
        r_true = A_mat @ psi_seq[k + 1] - (
            q + (S_mat + F_mat) @ psi_seq[k + 1]
        )
        np.testing.assert_allclose(
            np.linalg.norm(r_rhs), np.linalg.norm(r_true),
            rtol=1e-6, atol=cancellation_floor,
            err_msg=f"iter {k}: the free identity is not the true "
                    f"residual on an exact-M toy",
        )
    if measurable < 10:
        pytest.fail(
            f"only {measurable} iterations above the cancellation floor — "
            f"the identity rows lost their teeth (fixture drift)"
        )
    # And the exit claim is honest: the converged iterate's true residual.
    r_final = A_mat @ psi - (q + (S_mat + F_mat) @ psi)
    if not np.linalg.norm(r_final) / q_norm < 1e-10:
        pytest.fail("the exit claim is dishonest on an exact-M toy")


@pytest.mark.foundation
def test_zero_gain_exits_after_one_apply_with_exact_residual(rng):
    """C5 — zero gains ⟹ rhs ≡ q is constant ⟹ the first comparison sees
    r = 0 EXACTLY: one inverse apply, history == [0.0], and the returned
    iterate is the exact solve (a zero-gain path that iterates or exits
    non-zero is a bookkeeping bug)."""
    n = 5
    A_mat = np.diag(rng.uniform(1.0, 2.0, n))
    A_inv = _RecordingInverse(A_mat)
    q = rng.random(n) + 0.5
    si = SourceIteration(A_inv, max_iter=50, tol=1e-12)
    psi, record = si.solve(q)
    history = _trajectory(record)
    if len(A_inv.rhs_seen) != 1:
        pytest.fail(
            f"zero-gain SI ran {len(A_inv.rhs_seen)} inverse applies — "
            f"expected exactly ONE (ψ₁ = A⁻¹q is exact when A = M)"
        )
    if history != [0.0]:
        pytest.fail(f"zero-gain history {history} — expected [0.0] exactly")
    np.testing.assert_allclose(psi, np.linalg.solve(A_mat, q), rtol=1e-13)


@pytest.mark.foundation
def test_zero_source_zero_start_exits_clean():
    """r5 — ``q_ext = 0`` with the zero cold start: the zero solution is
    found at the first comparison (res = 0/1e-30 = 0.0 exactly), no
    division blow-up, no spurious non-convergence."""
    n = 4
    A_inv = _RecordingInverse(np.eye(n) * 2.0)
    S = MatrixOperator(0.3 * np.eye(n))
    si = SourceIteration(A_inv, S, max_iter=50, tol=1e-12)
    psi, record = si.solve(np.zeros(n))
    history = _trajectory(record)
    np.testing.assert_array_equal(np.asarray(psi), 0.0)
    if history != [0.0]:
        pytest.fail(f"zero-source history {history} — expected [0.0]")
    if not np.isfinite(history).all():
        pytest.fail("the q≈0 guard leaked a non-finite residual")


@pytest.mark.foundation
@pytest.mark.catches("ERR-053")
def test_the_gmres_nonconvergence_warning_is_ESCALATABLE(monkeypatch) -> None:
    r"""#340 R3 — the tree's only non-convergence announcement from inside
    ``numerics`` must answer to the PUBLISHED escalation flag.

    ``ESCALATION_FLAG`` is
    ``-W error::orpheus.numerics.convergence.ConvergenceWarning``, and a
    filter names a category: a bare ``RuntimeWarning`` does not match it.
    So until 2026-08-10 a CI run could be configured to make truncation
    fatal, PASS, and still have swallowed this one — a gate not covering
    what its own recipe claims, which is the defect class #340 exists to
    remove.

    ⭐ The teeth are in the SECOND leg, and only there.  Asserting that the
    warning is caught by a ``RuntimeWarning`` filter is satisfied by the
    pre-fix code (``ConvergenceWarning`` IS a ``RuntimeWarning``), so that
    leg cannot discriminate and is not written.  Asserting the category is
    ``ConvergenceWarning`` reds the moment anyone widens it back — which is
    the regression this row exists to catch, and the one a reviewer would
    wave through as "it's only a warning category".

    The narrowing is one-directional and safe by construction: every
    existing consumer filters on the BASE class (``_krylov_warnings`` above
    matches ``issubclass(w.category, RuntimeWarning)``), so it keeps
    matching.  Widening back would break only this row.
    """
    import orpheus.numerics.iteration as iteration_mod
    from orpheus.numerics.convergence import ConvergenceWarning

    monkeypatch.setattr(
        iteration_mod.spla, "gmres",
        _stub_gmres_returning(info=7, pr_norm_tail=0.37),
    )
    A_op, b = _singular_consistent()

    def _fresh_krylov():
        return KrylovAcceleration(
            A_op, preconditioner=IdentityOperator(),
            tol=1e-12, max_iter=30, restart=2,
        )

    # leg 1 — the category IS the escalatable one, not merely a RuntimeWarning.
    with pytest.warns(ConvergenceWarning) as caught:
        _fresh_krylov().solve(b)
    assert any("ERR-053" in str(w.message) for w in caught), (
        "the escalatable warning must still carry its ERR-053 pointer — the "
        "category is what CI filters on, the pointer is what a human follows"
    )

    # leg 2 — the escalation the published flag performs, done in-process:
    # the SAME solve must RAISE once that category is an error.  This is the
    # leg with teeth; it reds on any widening of the category.
    with warnings.catch_warnings():
        warnings.simplefilter("error", ConvergenceWarning)
        with pytest.raises(ConvergenceWarning):
            _fresh_krylov().solve(b)


# ───────────────────────────────────────────────────────────────────────
# ERR-079 / #349 — the budget and the trajectory count the SAME thing
# ───────────────────────────────────────────────────────────────────────
#
# ⭐ These are the gates whose ABSENCE let ERR-079 ship.  `[M]` 2026-08-13:
# every one of the 58 `IterationRecord`s in the suite was hand-built, so no
# gate had ever read `budget` / `n_iterations` / `exhausted_budget` off a
# record produced by a real `KrylovAcceleration.solve` — including the
# fixture LABELLED `inner(gmres)`, which omits `iterations_run` and so does
# not even reproduce the producer's shape.  A synthetic record cannot catch a
# producer's units error, because the test author picks both numbers.


@pytest.mark.foundation
@pytest.mark.catches("ERR-079")
def test_a_CONVERGED_gmres_solve_did_not_exhaust_its_budget():
    """A healthy Krylov solve must not report that it ran out (#349).

    ⛔ The defect, verbatim: ``exhausted_budget`` was
    ``n_iterations >= budget`` over a pair in different units — scipy's
    ``maxiter`` counts restart CYCLES while the ``pr_norm`` callback fires per
    inner ARNOLDI STEP.  Raising the knob did not help, because the two
    numbers were never commensurable.

    ⭐ The fixture is sized so the OLD spelling would say the opposite, and
    that precondition is ASSERTED rather than assumed — a row that silently
    stopped discriminating would otherwise go on reading green forever
    (``vv`` #19: only the reading under the wrong structure carries
    information).  `[M]` 30 distinct eigenvalues need 29 Arnoldi steps
    against a budget of 5 cycles, so ``n_iterations >= limit`` holds by a
    factor of ~6 while the honest ceiling is 150.
    """
    n = 30
    A = MatrixOperator(np.diag(np.linspace(1.0, 10.0, n)), can_solve=True)
    max_iter, restart = 5, 30

    krylov = KrylovAcceleration(
        A, ZeroOperator(), preconditioner=IdentityOperator(),
        max_iter=max_iter, tol=1e-12, restart=restart,
    )
    _, record = krylov.solve(np.ones(n))
    budget = record.budget

    # The producer states the exchange rate, and it is scipy's own.
    assert budget.iterations_per_unit == min(restart, n)
    assert budget.limit == max_iter
    assert budget.in_iterations == max_iter * min(restart, n)

    # ⚠ The anti-dud precondition: without this, the row below could pass on a
    # fixture that never entered the defect's regime at all.
    assert record.n_iterations >= budget.limit, (
        f"fixture no longer discriminates — it ran "
        f"{record.n_iterations} steps against a limit of {budget.limit}, so "
        f"the retired `n_iterations >= budget` spelling would have agreed "
        f"with the correct one and this gate would be vacuous"
    )

    assert record.converged is True
    assert record.exhausted_budget is False, (
        "a solve whose criterion cleared did not stop because it ran out"
    )
    assert record.truncated is False
    assert record.status == "CONVERGED"


@pytest.mark.foundation
@pytest.mark.catches("ERR-079")
def test_a_gmres_solve_that_really_DID_run_out_still_says_so():
    """The negative leg — without it, ``exhausted_budget = False`` passes.

    ``vv`` #11: the row above validates that the property stops crying wolf;
    only this one validates that it can still cry.  ``restart=1`` collapses
    the cycle and the step to the same thing, so the honest ceiling IS the
    knob and a hard system exhausts it.
    """
    n = 30
    rng_local = np.random.default_rng(0)
    hard = np.eye(n) + 0.9 * rng_local.standard_normal((n, n)) / np.sqrt(n)
    A = MatrixOperator(hard, can_solve=True)
    max_iter = 4

    krylov = KrylovAcceleration(
        A, ZeroOperator(), preconditioner=IdentityOperator(),
        max_iter=max_iter, tol=1e-14, restart=1,
    )
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        _, record = krylov.solve(np.ones(n))

    assert record.budget.iterations_per_unit == 1
    assert record.budget.in_iterations == max_iter
    assert record.n_iterations == max_iter
    assert record.converged is False
    assert record.exhausted_budget is True
    assert record.truncated is True
    assert record.status == "TRUNCATED"


@pytest.mark.foundation
@pytest.mark.catches("ERR-079")
def test_the_advice_names_a_setting_in_the_KNOBs_units_not_the_trajectorys():
    """*"set max_inner=N"* must be typeable into ``max_inner`` (#349).

    ``projected_iterations`` fits the observed rate over the TRAJECTORY, so
    its answer is in Arnoldi steps; the knob takes cycles.  Before the fix
    both halves of *"needs about N iterations: set X=N"* printed the same N,
    which on this arm over-states the required setting by ``restart``.

    The two numbers below must DIFFER for the gate to mean anything, so that
    is asserted too — on a unit-consistent producer they coincide by design
    and the row would be a tautology.
    """
    n = 30
    rng_local = np.random.default_rng(1)
    hard = np.eye(n) + 0.9 * rng_local.standard_normal((n, n)) / np.sqrt(n)
    A = MatrixOperator(hard, can_solve=True)

    krylov = KrylovAcceleration(
        A, ZeroOperator(), preconditioner=IdentityOperator(),
        max_iter=2, tol=1e-14, restart=8,
    )
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        _, record = krylov.solve(np.ones(n))

    projected = record.projected_iterations()
    assert projected is not None, "the fixture must be projectable"
    setting = record.budget.covering(projected)

    assert setting != projected, (
        "this arm's knob and trajectory are in different units, so a gate "
        "where the two coincide is not testing the conversion"
    )
    # The recommendation must actually BUY what was projected...
    assert setting * record.budget.iterations_per_unit >= projected
    # ...and not a cycle more than needed.
    assert (setting - 1) * record.budget.iterations_per_unit < projected
