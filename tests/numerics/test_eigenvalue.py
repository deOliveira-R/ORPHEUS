r"""``DirectEigenvalue`` — the direct (non-iterative) dominant eigenpair (L0/L1).

#276 P4-D. ``DirectEigenvalue`` is the DIRECT sibling of
:func:`~orpheus.numerics.eigenvalue.power_iteration`: the EXACT dominant
eigenpair of the generalized eigenproblem

.. math::

    A\,\varphi \;=\; \frac{1}{k}\,F\,\varphi
    \qquad\Longleftrightarrow\qquad
    A^{-1}F\,\varphi \;=\; k\,\varphi ,

returned in one LAPACK call (``eig(solve(A, F))``) for SMALL densifiable
operators. ``power_iteration`` converges to the SAME eigenpair iteratively
for large sweep-only operators.

This file is the **headline pure-math gate**: it verifies ``DirectEigenvalue``
against TRANSPORT-UNRELATED linear algebra with **hand-derived** closed-form
eigenpairs. Structural independence is by construction — every reference
eigenvalue/eigenvector is either (a) hand-verified ``M @ v == k·v``, (b) built
from chosen eigenpairs via ``M = V diag(λ) V⁻¹`` (recovery asserted), or (c) a
closed-form rank-1 identity ``k = vᵀA⁻¹u``. **No reference value is produced by
calling the same ``np.linalg.eig`` the SUT uses.** Once verified this way,
``DirectEigenvalue`` becomes trusted machinery whose structural independence as
a *verification tool* comes from the consumer assembling ``A``/``F`` differently
(e.g. the homogeneous solver), not from the eigensolver.

vv claim layers (1.5 gate)
    * Convergence-order: n/a (direct solver, no discretisation).
    * Flux-shape: the eigenVECTOR recovery (``TestStandardEigenproblem``,
      ``TestGeneralizedEigenproblem``) — a model/shape claim, pinned by the
      intrinsic generalized-eigenpair residual law ``F φ = k A φ``.
    * Eigenvalue: the dominant-VALUE recovery — a closed-form (pillar 1) claim,
      the ONLY pillar that proves eigenvalues (MMS cannot). The references here
      ARE closed-form, so the pillar matches the claim.

vv Mode-8 (``-O`` strips bare ``assert``): structural gates route through
:func:`_require` (a ``pytest.fail`` call, fires under ``-O``); value gates
through ``np.testing.*``. The canonical invocation is ``python -O -m pytest``.

vv Mode-10/teeth: ``TestGatesHaveTeeth`` carries an executable reference of the
SUT body parametrised by the documented mutations (argmax→argmin, drop the sign
flip, drop ``real()``, ``A⁻¹F``→``F⁻¹A``) and proves each mutation reddens the
gate that claims to catch it — confirmed NOW, before the SUT lands.
"""
from __future__ import annotations

import warnings
from typing import Any

import numpy as np
import pytest
from numpy.typing import ArrayLike

import orpheus.numerics.eigenvalue as _eig

pytestmark = pytest.mark.l1


# ───────────────────────────────────────────────────────────────────────────
# Binding adapter (PRE-IMPL tolerant + API-shape agnostic).
#
# The final spelling of the primitive is open (see the spec memo): a class
# ``DirectEigenvalue(A, F).solve() -> (k, phi)`` OR a free function
# ``direct_eigenvalue(A, F) -> (k, phi)``. This adapter binds to whichever
# ships and SKIPS cleanly while neither exists, so the file collects green
# pre-impl (consistent with the ``importorskip`` pre-impl pattern elsewhere).
# When the SUT lands, every gate below exercises the REAL primitive.
# ───────────────────────────────────────────────────────────────────────────


def _solve(A: np.ndarray, F: np.ndarray) -> tuple[float, np.ndarray]:
    """Return ``(k, phi)`` from the new direct dominant-eigenpair primitive."""
    cls = getattr(_eig, "DirectEigenvalue", None)
    if cls is not None:
        return cls(A, F).solve()
    fn = getattr(_eig, "direct_eigenvalue", None)
    if fn is not None:
        return fn(A, F)
    pytest.skip(
        "#276 P4-D PRE-IMPL: neither DirectEigenvalue nor direct_eigenvalue is "
        "on orpheus.numerics.eigenvalue yet."
    )


def _require(condition: object, message: str) -> None:
    """Fail with ``message`` if false. Fires under ``-O`` (NOT a bare assert).

    ``condition`` is typed ``object`` (not ``bool``) so numpy scalar/0-d
    boolean results — ``np.bool_`` from a comparison, a 0-d ``NDArray[bool_]``
    from a matmul-then-compare — bind cleanly; ``not condition`` resolves them
    via ``__bool__`` at runtime.
    """
    if not condition:
        pytest.fail(message)


def _unit(v: ArrayLike) -> np.ndarray:
    """Unit-norm copy (for direction-up-to-scale comparisons)."""
    arr = np.asarray(v, dtype=float).ravel()
    return arr / np.linalg.norm(arr)


def _assert_generalized_eigenpair(A, F, k, phi, *, rtol=1e-10) -> None:
    r"""The intrinsic LAW: ``F φ = k A φ`` (i.e. ``A⁻¹F φ = k φ``).

    The STRONGEST structurally-independent check — it needs no a-priori
    eigenvector and catches ``A⁻¹F``→``F⁻¹A``, transposes, and a wrong
    dominant selection that is nonetheless self-consistent only if it is a
    genuine eigenpair. Equivalent to the brief's ``‖A φ − (1/k) F φ‖`` form,
    written division-free.
    """
    A = np.asarray(A, dtype=float)
    F = np.asarray(F, dtype=float)
    phi = np.asarray(phi, dtype=float).ravel()
    lhs = F @ phi
    rhs = k * (A @ phi)
    scale = np.linalg.norm(rhs) + np.linalg.norm(lhs)
    resid = np.linalg.norm(lhs - rhs)
    _require(
        resid <= rtol * scale + 1e-12,
        f"generalized-eigenpair law violated: ‖Fφ − kAφ‖={resid:.3e} "
        f"(scale={scale:.3e}); (k, φ) is not an eigenpair of A⁻¹F.",
    )


# ═══════════════════════════════════════════════════════════════════════════
# PILLAR 1 — standard eigenproblem (A = I, F = M): dominant eigenpair of M.
# ═══════════════════════════════════════════════════════════════════════════


class TestStandardEigenproblem:
    """``A = I`` ⟹ ``A⁻¹F = F``; the dominant eigenpair of a known matrix."""

    def test_symmetric_2x2_sanity(self):
        """L1 sanity: ``[[2,1],[1,2]]`` → k=3, φ ∝ [1,1] (hand-derived).

        Hand proof: ``M @ [1,1] = [3,3] = 3·[1,1]`` and ``M @ [1,-1] =
        [1,-1]`` so the spectrum is {3, 1}, dominant 3. Symmetric, so this is
        a WEAK eigenvector-teeth case (Mᵀ = M; a transpose bug is invisible) —
        it is here as a floor, the teeth live in ``test_nonsymmetric_*``.
        """
        M = np.array([[2.0, 1.0], [1.0, 2.0]])
        k, phi = _solve(np.eye(2), M)
        np.testing.assert_allclose(k, 3.0, rtol=1e-12, atol=0)
        np.testing.assert_allclose(np.abs(_unit(phi)), _unit([1.0, 1.0]), rtol=1e-10)
        _assert_generalized_eigenpair(np.eye(2), M, k, phi)

    def test_nonsymmetric_eigenvector_recovered(self):
        """L1 eigenvector teeth: non-symmetric ``M = V diag(λ) V⁻¹``.

        λ = [4, 2, 1] (distinct, dominant 4); dominant eigenvector is the
        chosen ``V[:,0] = [2, 1, 1]`` (non-trivial: not axis-aligned, not
        all-equal). Recovery of THIS vector is what a variable-swap /
        transpose in the eigenvector extraction would break — a symmetric or
        axis-aligned case cannot see that bug.
        """
        lam = np.array([4.0, 2.0, 1.0])
        V = np.array([[2.0, 1.0, 0.0], [1.0, -1.0, 1.0], [1.0, 0.0, -1.0]])
        M = V @ np.diag(lam) @ np.linalg.inv(V)  # built from CHOSEN eigenpairs
        k, phi = _solve(np.eye(3), M)
        np.testing.assert_allclose(k, 4.0, rtol=1e-12, atol=0)
        # direction up to sign+scale; the chosen v0 = [2,1,1].
        cos = abs(_unit(phi) @ _unit(V[:, 0]))
        _require(cos > 1 - 1e-10, f"dominant eigenvector misaligned: |cos|={cos:.12f}.")
        _assert_generalized_eigenpair(np.eye(3), M, k, phi)

    def test_argmax_selects_dominant_in_middle_of_lapack_order(self):
        """L1 dominant-selection teeth: diagonal ``F`` → eig order is the

        diagonal order. ``F = diag([2, 5, 1])`` (A = I) returns eigenvalues
        in EXACTLY [2, 5, 1] order, so the dominant (5) sits in the MIDDLE —
        ``argmax`` must pick index 1, distinct from index 0 (→2), index −1
        (→1), and ``argmin`` (→1). A diagonal matrix's eigenvalues ARE its
        diagonal (no eig call needed for the ground truth); numpy returns them
        in order with eigenvectors = I. φ = e₁ here (axis-aligned ON PURPOSE —
        the job is the SELECTION value, the vector teeth are elsewhere).
        """
        F = np.diag([2.0, 5.0, 1.0])
        k, phi = _solve(np.eye(3), F)
        np.testing.assert_allclose(k, 5.0, rtol=1e-12, atol=0)
        np.testing.assert_allclose(_unit(phi), _unit([0.0, 1.0, 0.0]), rtol=1e-10)


# ═══════════════════════════════════════════════════════════════════════════
# PILLAR 1 — genuine generalized problem (A ≠ I, F ≠ I): exercises A⁻¹.
# ═══════════════════════════════════════════════════════════════════════════


class TestGeneralizedEigenproblem:
    r"""``A`` non-diagonal AND ``F`` ≠ I — the form the homogeneous solver uses.

    The homogeneous ``F = χ ⊗ νΣ_f`` is rank-1 ``u vᵀ``; then ``A⁻¹F =
    (A⁻¹u) vᵀ`` is rank-1 with a SINGLE nonzero eigenvalue and closed forms

        k = vᵀ A⁻¹ u   (the unique nonzero eigenvalue),
        φ = A⁻¹ u      (its eigenvector direction),

    both hand-computable. ``A⁻¹``→genuine solve (A non-diagonal), so this is
    the case that pins the ``np.linalg.solve(A, F)`` step (a standard ``A = I``
    case leaves the solve trivial).
    """

    @staticmethod
    def _rank1_case():
        A = np.array([[3.0, 1.0, 0.0], [1.0, 4.0, 1.0], [0.0, 1.0, 5.0]])
        u = np.array([1.0, 2.0, 3.0])  # χ
        v = np.array([0.5, 1.5, 0.7])  # νΣ_f
        F = np.outer(u, v)
        Ainv_u = np.linalg.solve(A, u)  # reference built independently of eig
        k_closed = float(v @ Ainv_u)  # k = vᵀ A⁻¹ u
        return A, F, u, v, Ainv_u, k_closed

    def test_rank1_closed_form_eigenvalue(self):
        """L1 eigenvalue (closed form): k = vᵀ A⁻¹ u for rank-1 F = u vᵀ.

        The reference uses ``np.linalg.solve`` (a DIFFERENT primitive than the
        SUT's ``eig``) on a vector — structurally independent of the
        eigensolver. The other two eigenvalues are exactly 0, so ``argmax``
        must pick the nonzero one (an argmax→argmin mutation would return ~0).
        """
        A, F, _u, _v, _Ainv_u, k_closed = self._rank1_case()
        k, _phi = _solve(A, F)
        np.testing.assert_allclose(k, k_closed, rtol=1e-11, atol=0)

    def test_rank1_eigenvector_is_Ainv_u(self):
        """L1 flux-shape: the dominant eigenvector direction is A⁻¹u."""
        A, F, _u, _v, Ainv_u, _k = self._rank1_case()
        _k_got, phi = _solve(A, F)
        cos = abs(_unit(phi) @ _unit(Ainv_u))
        _require(cos > 1 - 1e-10, f"rank-1 eigenvector ≠ A⁻¹u: |cos|={cos:.12f}.")

    def test_rank1_residual_law(self):
        """L1 intrinsic: ``F φ = k A φ`` (the generalized-eigenpair law)."""
        A, F, _u, _v, _Ainv_u, _k = self._rank1_case()
        k, phi = _solve(A, F)
        _assert_generalized_eigenpair(A, F, k, phi)

    def test_full_rank_generalized_via_construction(self):
        r"""L1 full-rank generalized: ``A⁻¹F = M`` with chosen eigenpairs.

        Build ``M = V diag(λ) V⁻¹`` (known dominant), pick an invertible non-
        diagonal ``A``, set ``F = A M`` → ``A⁻¹F = M`` exactly. This exercises
        the solve with a FULL-rank F (the rank-1 case has n−1 zero eigenvalues;
        here all three are nonzero), so a bug that only shows on full-rank F is
        caught.
        """
        lam = np.array([6.0, 3.0, 2.0])
        V = np.array([[1.0, 1.0, 0.0], [0.0, 1.0, 1.0], [1.0, 0.0, 1.0]])
        M = V @ np.diag(lam) @ np.linalg.inv(V)
        A = np.array([[2.0, 0.5, 0.0], [0.5, 3.0, 0.5], [0.0, 0.5, 4.0]])
        F = A @ M  # A⁻¹F = M, eigenpairs of M are known by construction
        k, phi = _solve(A, F)
        np.testing.assert_allclose(k, 6.0, rtol=1e-10, atol=0)
        cos = abs(_unit(phi) @ _unit(V[:, 0]))
        _require(cos > 1 - 1e-9, f"generalized dominant eigenvector misaligned: |cos|={cos}.")
        _assert_generalized_eigenpair(A, F, k, phi)


# ═══════════════════════════════════════════════════════════════════════════
# Sign normalisation — φ.sum() ≥ 0 (the φ output convention).
# ═══════════════════════════════════════════════════════════════════════════


class TestSignNormalisation:
    r"""The SUT flips φ so ``φ.sum() ≥ 0`` (a physical, non-negative-spectrum
    convention). The teeth case is hand-verified AND lands on the negative-sum
    orientation as the SUT's eig backend returns it.
    """

    # M = Mᵀ, hand-verified eigenpair: M @ [1,1,0] = [4,4,0] = 4·[1,1,0].
    # Spectrum {−2, 4, 2} (trace = 4 = −2+4+2 ✓), dominant 4. As numpy returns
    # it the dominant eigenvector is the NEGATIVE-sum orientation [−1,−1,0]/√2,
    # so the φ.sum() ≥ 0 flip is genuinely exercised here.
    _M = np.array([[1.0, 3.0, -2.0], [3.0, 1.0, -2.0], [-2.0, 2.0, 2.0]])

    def test_dominant_eigenpair_is_hand_value(self):
        """L1: hand-derived eigenpair (4, [1,1,0]) is recovered."""
        k, phi = _solve(np.eye(3), self._M)
        np.testing.assert_allclose(k, 4.0, rtol=1e-12, atol=0)
        cos = abs(_unit(phi) @ _unit([1.0, 1.0, 0.0]))
        _require(cos > 1 - 1e-10, f"eigenvector ≠ [1,1,0]: |cos|={cos:.12f}.")

    def test_returned_sum_is_nonnegative(self):
        """L1 sign contract: φ.sum() ≥ 0 (the output convention)."""
        _k, phi = _solve(np.eye(3), self._M)
        s = float(np.asarray(phi, dtype=float).sum())
        _require(s >= 0.0, f"φ.sum()={s:.6e} < 0 — the sign normalisation did not fire.")

    def test_sign_flip_is_exercised_here(self):
        """Teeth precondition (HONEST, never vacuous): on THIS matrix the SUT's
        eig backend returns a NEGATIVE-sum dominant eigenvector, so the flip is
        load-bearing. If a future LAPACK changes its sign convention this
        precondition reddens (telling us to repick), rather than the φ.sum()≥0
        gate going silently vacuous.
        """
        raw_vals, raw_vecs = np.linalg.eig(self._M)  # the SUT's backend, raw
        dom = int(np.argmax(np.real(raw_vals)))
        raw = np.real(raw_vecs[:, dom])
        _require(
            raw.sum() < 0,
            f"PRECONDITION BROKEN: raw dominant eigenvector sum={raw.sum():.4e} ≥ 0 "
            f"on this LAPACK — the sign-flip gate no longer has teeth here; pick a "
            f"new negative-sum matrix.",
        )


# ═══════════════════════════════════════════════════════════════════════════
# Complex spectrum — real()-extraction of a REAL dominant eigenpair.
# ═══════════════════════════════════════════════════════════════════════════


class TestComplexSpectrum:
    r"""A real matrix can carry a complex-conjugate eigenvalue pair; when the
    DOMINANT (max-real) eigenvalue is real, the SUT's ``real()`` must return a
    real ``float`` k and a real-dtype φ. The presence of complex eigenvalues
    forces ``eig`` to a complex128 dtype, so a dropped ``real()`` leaks a
    complex k / φ — caught by the dtype assertions.
    """

    # block diag(5, rot(1, 3)): eigenvalues {5, 1+3i, 1-3i}; dominant 5 (real),
    # eigenvector e₀.
    _M = np.array([[5.0, 0.0, 0.0], [0.0, 1.0, -3.0], [0.0, 3.0, 1.0]])

    def test_real_dominant_returns_real_float(self):
        """L1: dominant real eigenvalue 5.0 returned as a real ``float``."""
        k, _phi = _solve(np.eye(3), self._M)
        _require(
            not isinstance(k, complex) and float(np.imag(k)) == 0.0,
            f"k must be a real float; got {k!r} (real() dropped?).",
        )
        np.testing.assert_allclose(float(np.real(k)), 5.0, rtol=1e-12, atol=0)

    def test_returned_eigenvector_is_real_dtype(self):
        """L1: φ is a real-dtype array (not complex128) — the real() teeth."""
        _k, phi = _solve(np.eye(3), self._M)
        phi = np.asarray(phi)
        _require(
            np.isrealobj(phi),
            f"φ must be real-dtype; got dtype={phi.dtype} (real() dropped?).",
        )
        np.testing.assert_allclose(_unit(phi), _unit([1.0, 0.0, 0.0]), rtol=1e-10)


# ═══════════════════════════════════════════════════════════════════════════
# Intrinsic properties — the DEFINING laws of the type (project standard:
# every math-bearing type ships a test of its invariants).
# ═══════════════════════════════════════════════════════════════════════════


class TestIntrinsicProperties:
    """Determinism + dominance + the residual law, on a generic case."""

    @staticmethod
    def _generic():
        lam = np.array([5.0, 3.0, 1.5])
        V = np.array([[2.0, 1.0, 0.0], [1.0, 1.0, 1.0], [1.0, 0.0, 2.0]])  # det 3
        M = V @ np.diag(lam) @ np.linalg.inv(V)
        A = np.array([[2.0, 0.3, 0.0], [0.3, 2.5, 0.4], [0.0, 0.4, 3.0]])
        F = A @ M
        return A, F, V, lam

    def test_determinism_bit_identical(self):
        """Two solves of the same (A, F) return bit-identical (k, φ)."""
        A, F, _V, _lam = self._generic()
        k1, p1 = _solve(A, F)
        k2, p2 = _solve(A, F)
        _require(k1 == k2, f"non-deterministic k: {k1!r} vs {k2!r}.")
        _require(
            np.array_equal(np.asarray(p1), np.asarray(p2)),
            "non-deterministic φ (not bit-identical across two solves).",
        )

    def test_returns_the_dominant_max_real_eigenvalue(self):
        """The returned k is ≥ the real part of every eigenvalue of A⁻¹F.

        The bound is computed from the CONSTRUCTED spectrum (max(λ)=5), not
        from an eig call — so it independently pins that the SELECTED eigenvalue
        is the dominant one (argmax→argmin would return 1.5 < 5).
        """
        A, F, _V, lam = self._generic()
        k, _phi = _solve(A, F)
        np.testing.assert_allclose(k, float(lam.max()), rtol=1e-10, atol=0)

    def test_residual_law_holds(self):
        """The generalized-eigenpair law ``F φ = k A φ`` (intrinsic)."""
        A, F, _V, _lam = self._generic()
        k, phi = _solve(A, F)
        _assert_generalized_eigenpair(A, F, k, phi)


# ═══════════════════════════════════════════════════════════════════════════
# Edge / robustness contracts (see the spec memo for the RECOMMENDED contract).
# ═══════════════════════════════════════════════════════════════════════════


class TestEdgeContracts:
    """Singular A, degenerate dominant, complex dominant — pin the contract."""

    def test_singular_A_raises(self):
        """RECOMMENDED contract: a singular ``A`` propagates
        ``np.linalg.LinAlgError`` (a singular loss matrix is a malformed
        problem — fail loud, do not return garbage). ``np.linalg.solve``
        already raises; the SUT must NOT swallow it.
        """
        A = np.array([[1.0, 2.0], [2.0, 4.0]])  # rank-1, singular
        F = np.eye(2)
        with pytest.raises(np.linalg.LinAlgError):
            _solve(A, F)

    def test_degenerate_dominant_returns_a_valid_eigenpair(self):
        """RECOMMENDED contract: a repeated dominant eigenvalue returns the
        FIRST argmax; the eigenvector is implementation-defined WITHIN the
        degenerate eigenspace but MUST satisfy the residual law. We pin only k
        and the law (NOT a specific direction — that is genuinely ambiguous).
        """
        M = np.diag([4.0, 4.0, 1.0])  # dominant 4 with multiplicity 2
        k, phi = _solve(np.eye(3), M)
        np.testing.assert_allclose(k, 4.0, rtol=1e-12, atol=0)
        _assert_generalized_eigenpair(np.eye(3), M, k, phi)

    def test_complex_dominant_contract(self):
        r"""RECOMMENDED contract: when the max-real eigenvalue is itself COMPLEX
        (a conjugate pair dominates), the spec'd ``real()`` would SILENTLY
        return a non-eigenvalue (residual ``‖Fφ − kAφ‖`` ≠ 0). The principled
        choice is to RAISE (the dominant of a Perron–Frobenius resolvent must
        be real; a complex dominant means malformed input). This test asserts
        the raise IF the guard is implemented; if the implementer keeps the
        silent ``real()`` fallback, it asserts the residual law EXPOSES the
        corruption instead (so the corruption is never invisible).
        """
        # eigenvalues {3+2i, 3-2i, 2}: the dominant max-real is the complex pair.
        M = np.array([[3.0, -2.0, 0.0], [2.0, 3.0, 0.0], [0.0, 0.0, 2.0]])
        try:
            k, phi = _solve(np.eye(3), M)
        except (ValueError, np.linalg.LinAlgError):
            return  # guard implemented (the RECOMMENDED fail-loud) — contract met.
        # No guard: the silent real() fallback MUST be detectable as a non-eigenpair.
        A = np.eye(3)
        resid = np.linalg.norm(M @ np.asarray(phi, float).ravel()
                               - k * (A @ np.asarray(phi, float).ravel()))
        _require(
            resid > 1e-6,
            "complex-dominant fallback returned something that passes the "
            "eigenpair law — impossible unless the guard silently fixed it; "
            "review the contract.",
        )


# ═══════════════════════════════════════════════════════════════════════════
# Mode-10 TEETH PROOF — every gate is provably ABLE to red.
#
# The SUT is not implemented yet, so this class carries an EXECUTABLE REFERENCE
# of the documented SUT body, parametrised by each named mutation, and proves
# the gate that CLAIMS to catch the mutation actually reddens. This runs NOW
# (no SUT needed), doubles as a precise executable spec, and is the §0.5
# standing-discipline evidence that no gate is a false green.
#
# The reference is NOT the production primitive (the real gates bind via
# ``_solve``); it exists ONLY to validate the teeth.
# ═══════════════════════════════════════════════════════════════════════════


def _ref_solve(A, F, *, mut: str | None = None) -> tuple[Any, np.ndarray]:
    """Executable reference of the spec'd SUT body, with optional mutation.

    Clean (``mut=None``) == the documented behaviour. Each ``mut`` is one of
    the §3 mutations the spec memo enumerates.
    """
    if mut == "invert_swap":  # A⁻¹F → F⁻¹A
        M = np.linalg.solve(F, A)
    else:
        M = np.linalg.solve(A, F)
    eigvals, eigvecs = np.linalg.eig(M)
    if mut == "argmin":  # argmax → argmin
        dom = int(np.argmin(np.real(eigvals)))
    else:
        dom = int(np.argmax(np.real(eigvals)))
    if mut == "drop_real_k":  # real(eigvals[dom]) → eigvals[dom]
        k: Any = eigvals[dom]
    else:
        k = float(eigvals[dom].real)
    if mut == "drop_real_phi":  # real(eigvecs[:,dom]) → eigvecs[:,dom]
        phi = eigvecs[:, dom]
    else:
        phi = np.real(eigvecs[:, dom])
    if mut != "drop_sign":  # the φ.sum() ≥ 0 flip
        if np.real(np.sum(phi)) < 0:
            phi = -phi
    return k, phi


class TestGatesHaveTeeth:
    """Each documented mutation reddens the gate that claims to catch it."""

    def test_argmin_mutation_reddens_value_gate(self):
        """argmax→argmin: PM-4 value gate (k==5) reddens (returns 1, the min)."""
        F = np.diag([2.0, 5.0, 1.0])
        k_clean, _ = _ref_solve(np.eye(3), F)
        k_mut, _ = _ref_solve(np.eye(3), F, mut="argmin")
        _require(np.isclose(k_clean, 5.0), "clean ref should select 5.")
        _require(
            not np.isclose(k_mut, 5.0),
            f"argmin mutation did NOT move k off 5.0 (got {k_mut}); gate has no teeth.",
        )

    def test_argmin_mutation_reddens_rank1_value_gate(self):
        """argmax→argmin also reddens the rank-1 closed-form gate (returns ~0)."""
        A = np.array([[3.0, 1.0, 0.0], [1.0, 4.0, 1.0], [0.0, 1.0, 5.0]])
        u, v = np.array([1.0, 2.0, 3.0]), np.array([0.5, 1.5, 0.7])
        F = np.outer(u, v)
        k_closed = float(v @ np.linalg.solve(A, u))
        k_mut, _ = _ref_solve(A, F, mut="argmin")
        _require(
            not np.isclose(k_mut, k_closed, rtol=1e-6),
            f"argmin did not move k off the closed form {k_closed:.6f} (got {k_mut}).",
        )

    def test_drop_sign_mutation_reddens_sign_gate(self):
        """Drop the flip: φ.sum() < 0 on the PM-5 matrix → the ≥0 gate reddens."""
        M = TestSignNormalisation._M
        _k, phi = _ref_solve(np.eye(3), M, mut="drop_sign")
        _require(
            float(np.real(np.sum(phi))) < 0,
            "drop-sign mutation left φ.sum() ≥ 0 on the PM-5 matrix; the sign "
            "gate would not catch it (LAPACK returned positive-sum — repick).",
        )

    def test_drop_real_mutations_redden_dtype_gates(self):
        """Drop real(): k goes complex and φ goes complex-dtype on PM-6."""
        M = TestComplexSpectrum._M
        k_mut, _ = _ref_solve(np.eye(3), M, mut="drop_real_k")
        _require(
            isinstance(k_mut, complex) or np.asarray(k_mut).dtype.kind == "c",
            f"drop_real_k left k real ({k_mut!r}); dtype gate has no teeth.",
        )
        _k, phi_mut = _ref_solve(np.eye(3), M, mut="drop_real_phi")
        _require(
            not np.isrealobj(np.asarray(phi_mut)),
            "drop_real_phi left φ real-dtype; dtype gate has no teeth.",
        )

    def test_invert_swap_mutation_reddens_residual_and_value(self):
        """A⁻¹F→F⁻¹A: the dominant-value AND the residual law redden.

        Uses a FULL-RANK F (= A·M, det ≠ 0) so ``F⁻¹A`` exists — a rank-1 F is
        singular, so the swap would merely raise (a different, weaker signal).
        """
        lam = np.array([6.0, 3.0, 2.0])
        V = np.array([[1.0, 1.0, 0.0], [0.0, 1.0, 1.0], [1.0, 0.0, 1.0]])
        M = V @ np.diag(lam) @ np.linalg.inv(V)
        A = np.array([[2.0, 0.5, 0.0], [0.5, 3.0, 0.5], [0.0, 0.5, 4.0]])
        F = A @ M  # A⁻¹F = M, dominant 6; F is full rank
        k_mut, phi_mut = _ref_solve(A, F, mut="invert_swap")
        phi_mut = np.asarray(phi_mut, dtype=float).ravel()
        # value gate reddens (the swapped dominant is ~0.5, not 6)
        _require(
            not np.isclose(float(np.real(k_mut)), 6.0, rtol=1e-6),
            f"invert_swap did not move k off the dominant 6.0 (got {k_mut}).",
        )
        # residual law reddens: F φ ≠ k A φ for the swapped pair
        resid = np.linalg.norm(F @ phi_mut - float(np.real(k_mut)) * (A @ phi_mut))
        scale = np.linalg.norm(float(np.real(k_mut)) * (A @ phi_mut)) + np.linalg.norm(F @ phi_mut)
        _require(
            resid > 1e-8 * scale + 1e-10,
            f"invert_swap pair satisfied the eigenpair law (resid={resid:.3e}); "
            f"the residual-law gate has no teeth against it.",
        )

    def test_clean_reference_passes_every_gate(self):
        """Control: the UN-mutated reference satisfies all the gates it should
        (so the teeth tests above are rejecting the MUTATION, not a broken
        reference).
        """
        # PM-4 value
        k, _ = _ref_solve(np.eye(3), np.diag([2.0, 5.0, 1.0]))
        _require(np.isclose(k, 5.0), "clean ref failed PM-4 value.")
        # PM-5 sign
        _k, phi = _ref_solve(np.eye(3), TestSignNormalisation._M)
        _require(float(np.real(np.sum(phi))) >= 0, "clean ref failed sign gate.")
        # PM-6 real dtype
        kk, pp = _ref_solve(np.eye(3), TestComplexSpectrum._M)
        _require(not isinstance(kk, complex) and np.isrealobj(np.asarray(pp)),
                 "clean ref failed real() gate.")
        # rank-1 residual law
        A = np.array([[3.0, 1.0, 0.0], [1.0, 4.0, 1.0], [0.0, 1.0, 5.0]])
        F = np.outer([1.0, 2.0, 3.0], [0.5, 1.5, 0.7])
        kc, phic = _ref_solve(A, F)
        phic = np.asarray(phic, float).ravel()
        resid = np.linalg.norm(F @ phic - float(np.real(kc)) * (A @ phic))
        _require(resid < 1e-8, f"clean ref failed residual law (resid={resid:.3e}).")


# ═══════════════════════════════════════════════════════════════════════════
# RQI — bordered-matrix Rayleigh-Quotient Iteration (#276 P4-E, #277)
#
# The third eigenvalue engine, verified against the now-trusted
# ``direct_eigenvalue`` oracle: RQI must converge to its EXACT eigenpair. The
# references stay structurally independent — the oracle is verified above
# against hand-derived closed forms, and the nearest-eigenvalue case uses
# CONSTRUCTED eigenpairs (``M = V diag(λ) V⁻¹``), not an eig call. The
# distinguishing teeth are the NEAREST-eigenvalue contract (not the dominant)
# and the QUADRATIC convergence (vs power iteration's linear rate).
# ═══════════════════════════════════════════════════════════════════════════


def _rqi(A, F, **kw) -> tuple[float, np.ndarray]:
    """Bind the RQI primitive (skips cleanly pre-impl)."""
    fn = getattr(_eig, "rayleigh_quotient_iteration", None)
    if fn is None:
        pytest.skip(
            "#276 P4-E PRE-IMPL: rayleigh_quotient_iteration not on "
            "orpheus.numerics.eigenvalue yet."
        )
    return fn(A, F, **kw)


class TestRayleighQuotientIteration:
    r"""Bordered-RQI polishes an estimate to the nearest eigenpair, superlinearly."""

    @staticmethod
    def _generalized():
        # A non-diagonal, F full-rank; A⁻¹F = M has the known spectrum {6, 3, 2}.
        lam = np.array([6.0, 3.0, 2.0])
        V = np.array([[1.0, 1.0, 0.0], [0.0, 1.0, 1.0], [1.0, 0.0, 1.0]])
        M = V @ np.diag(lam) @ np.linalg.inv(V)
        A = np.array([[2.0, 0.5, 0.0], [0.5, 3.0, 0.5], [0.0, 0.5, 4.0]])
        F = A @ M
        return A, F, V, lam

    def test_polishes_to_dominant_matches_direct_eigenvalue(self):
        """Warm-started near the dominant, RQI recovers ``direct_eigenvalue``'s
        EXACT (k, φ) and satisfies the generalized-eigenpair law."""
        A, F, _V, _lam = self._generalized()
        k_dir, phi_dir = _solve(A, F)  # the verified oracle (dominant)
        v0 = phi_dir + 0.1 * np.array([1.0, -0.5, 0.3])  # warm start, off-axis
        k, phi = _rqi(A, F, v0=v0)
        np.testing.assert_allclose(k, k_dir, rtol=1e-12, atol=0)
        cos = abs(_unit(phi) @ _unit(phi_dir))
        _require(cos > 1 - 1e-12, f"RQI eigenvector ≠ dominant: |cos|={cos:.14f}.")
        _assert_generalized_eigenpair(A, F, k, phi)

    def test_converges_to_nearest_eigenvalue_not_dominant(self):
        """The DEFINING RQI contract: warm-started near a SUBDOMINANT eigenvector,
        RQI converges to THAT eigenvalue (nearest the Rayleigh quotient), NOT the
        dominant. The reference eigenpair is CONSTRUCTED (``V diag(λ) V⁻¹``),
        independent of the eigensolver."""
        A, F, V, _lam = self._generalized()  # spectrum {6, 3, 2}; subdominant 3 = V[:,1]
        v0 = V[:, 1] + 0.05 * np.array([0.2, 0.1, -0.3])
        k, phi = _rqi(A, F, v0=v0)
        np.testing.assert_allclose(k, 3.0, rtol=1e-10, atol=0)  # the nearest, not 6
        _assert_generalized_eigenpair(A, F, k, phi)

    def test_residual_law_standard_problem(self):
        """A = I: RQI on a plain non-symmetric matrix satisfies ``M φ = k φ``."""
        lam = np.array([4.0, 2.0, 1.0])
        V = np.array([[2.0, 1.0, 0.0], [1.0, -1.0, 1.0], [1.0, 0.0, -1.0]])
        M = V @ np.diag(lam) @ np.linalg.inv(V)
        _k_dir, phi_dir = _solve(np.eye(3), M)
        k, phi = _rqi(np.eye(3), M, v0=phi_dir + 0.1 * np.array([0.3, -0.2, 0.5]))
        np.testing.assert_allclose(k, 4.0, rtol=1e-12, atol=0)
        _assert_generalized_eigenpair(np.eye(3), M, k, phi)

    def test_quadratic_convergence(self):
        """Superlinear teeth: from a moderate perturbation the error roughly SQUARES
        each step (locally quadratic) and reaches machine precision in a handful of
        steps — what distinguishes RQI from power iteration's linear rate.

        Re-runs from the SAME ``v0`` with ``max_iter`` = 1, 2, 3 (deterministic)
        and checks ``e_{n+1} ≤ C·e_n²`` with a generous ``C``, plus that 4 steps
        reach < 1e-10 while 2 steps do NOT (genuinely iterating, not one-shot).
        """
        A, F, _V, _lam = self._generalized()
        k_exact, phi_dir = _solve(A, F)
        v0 = phi_dir + 0.2 * np.array([1.0, -1.0, 0.5])

        with warnings.catch_warnings():  # low max_iter ⟹ expected non-converged warning
            warnings.simplefilter("ignore", RuntimeWarning)
            e = [abs(_rqi(A, F, v0=v0, max_iter=m)[0] - k_exact) for m in (1, 2, 3)]
            k4 = _rqi(A, F, v0=v0, max_iter=4)[0]

        C = 50.0  # generous quadratic constant
        _require(
            e[1] <= C * e[0] ** 2 + 1e-15,
            f"not quadratic at step 2: e2={e[1]:.3e} > C·e1²={C * e[0] ** 2:.3e}.",
        )
        _require(
            e[2] <= C * e[1] ** 2 + 1e-15,
            f"not quadratic at step 3: e3={e[2]:.3e} > C·e2²={C * e[1] ** 2:.3e}.",
        )
        _require(
            abs(k4 - k_exact) < 1e-10,
            f"4 RQI steps did not reach 1e-10 (|Δk|={abs(k4 - k_exact):.3e}).",
        )
        _require(
            e[1] > 1e-10,
            f"2 steps already at machine precision (e2={e[1]:.3e}) — perturbation too "
            f"small to demonstrate iteration; pick a larger v0 offset.",
        )

    def test_sign_normalised_output(self):
        """φ.sum() ≥ 0 output convention (matches direct_eigenvalue). Non-vacuous:
        warm-starting near ``−φ_dom`` makes RQI converge to the negative-sum
        orientation, so the sign flip is load-bearing here."""
        A, F, _V, _lam = self._generalized()
        _k_dir, phi_dir = _solve(A, F)
        v0 = -phi_dir + 0.1 * np.array([0.4, -0.3, 0.2])  # negative-sum hemisphere
        _k, phi = _rqi(A, F, v0=v0)
        s = float(np.asarray(phi).sum())
        _require(s >= 0.0, f"φ.sum()={s:.3e} < 0 — the sign normalisation did not fire.")

    def test_shape_validation_raises(self):
        """Non-square / mismatched shapes raise ValueError (boundary guard)."""
        with pytest.raises(ValueError):
            _rqi(np.ones((2, 3)), np.eye(2))
        with pytest.raises(ValueError):
            _rqi(np.eye(2), np.eye(3))
