r"""GreenOperator L0 gates — taxonomy §12 step 4 (#226; verification spec
PART III §20–§23).

**What is pinned.** The preconditioned-splitting inverse of a general
operator sum: ``(D − N).inverse()`` is a
:class:`~orpheus.numerics.green_operator.GreenOperator` whose ``apply``
returns :math:`(D-N)^{-1}q` with TRUE relative residual below ``tol`` —
verified against the structurally-independent dense LU
(:func:`numpy.linalg.solve` of the materialized sum), never against the
iteration's own bookkeeping.  The name-earning invariants (taxonomy §13):
**G-Neumann** (the partial sums of :math:`\sum_k (D^{-1}N)^k D^{-1} q` ARE
the Green application — a generic inverse has no splitting to satisfy) and
**G-reciprocity** (:math:`\langle\phi_2, G\phi_1\rangle = \langle
G^{\mathsf T}\phi_2, \phi_1\rangle` with :math:`G^{\mathsf T}` built
MANUALLY over the transposed operands — no ``.H`` minting, #280).  G-kernel
folds into the anchor via the :math:`\delta_j` input basis (spec §21
ruling).

**The fixture** (config-blindness §0.6): ``D`` = non-uniform diagonal, ``N``
= scaled NON-INVOLUTION cyclic permutation — NON-symmetric (a symmetric
``A`` makes G-reciprocity vacuous), non-commuting, ``N ≠ 0`` (a zero gain
nulls every gain mutation), with an EXACTLY computable spectral radius:
:math:`(D^{-1}\alpha P)^4 = \alpha^4/\prod d_i \cdot I` for the 4-cycle, so
:math:`\rho = \alpha/(\prod d_i)^{1/4}` and the Neumann 4-step decay ratio
is pinned TIGHTLY, not to a 5% band.

**Tolerance discipline (spec §20.0):** Green is ITERATIVE — its gates
assert at ``SAFETY × tol`` (driver-tol rule, lessons L7), NOT nulp; the one
row that breaks the §12.0 count-the-ops rule.  ``OperatorSum.solve`` was
DELETED at carve P4 (a generic sum's inverse action IS the GreenOperator;
solving is spelled ``sum.inverse().apply(b)``), so no solve-vs-inverse
tautology remains to (not) assert; the mixin back-half ``green.solve``
(the un-invert face) is what the back-half gate pins.

**ConvergenceFailure teeth** (spec §18.A + §22.2): the promise reads the
TRUE equation residual, not the driver's ρ-blind increment (Signature 9) —
pinned by the near-critical ``ρ=0.99`` case where the increment converges
while the equation is ~1e2 off; the refinement loop's positive control
(generous budget → the promise IS met) pins that the check-only false-raise
for ``ρ > 1/2`` is designed out.  The divergent split raises LOUDLY
(with a convergent control — a bare ``pytest.raises`` alone is Mode-10
blind).

**Teeth (spec §25, mutation-verified 2026-07-02, in-process patches):**
M-GRN-SIGN (gain sign `+t`) REDs the anchor + G-Neumann; M-GRN-SWAP
(precond↔gain) REDs the anchor; M-GRN-TOL (promise check deleted) REDs the
divergence gate; M-GRN-INCREMENT (increment instead of true residual) REDs
the near-critical gate; M-GRN-SEED (drop ``initial_guess`` to the driver
start) REDs ONLY the §23 spy here — the landed §14 spy is BLIND to it;
M-MIXIN-CODOM / M-MIXIN-INV RED the involution/round-trip on all three
siblings.
"""
from __future__ import annotations

import warnings

import numpy as np
import pytest

from orpheus.numerics.green_operator import ConvergenceFailure, GreenOperator
from orpheus.numerics.iteration import SourceIteration
from orpheus.numerics.operator import (
    DiagonalOperator,
    NotInvertible,
    PermutationOperator,
    ScaledOperator,
    ZeroOperator,
)

pytestmark = pytest.mark.foundation

_RNG = np.random.default_rng(226)

#: §20.0 driver-tol rule: iterative inverse ⇒ assertion rtol = SAFETY × tol
#: (each gate reads the operator's own ``green.tol``).
_SAFETY = 10.0

_D_DIAG = np.array([4.0, 5.0, 6.0, 7.0])
_P4 = np.roll(np.arange(4), 1)          # 4-cycle (non-involution)
_ALPHA = 2.0                            # gain scale: ρ = 2/(4·5·6·7)^¼ ≈ 0.372
_RHO = _ALPHA / float(np.prod(_D_DIAG)) ** 0.25
_Q = _RNG.uniform(1.0, 2.0, 4)


def _perm_matrix(perm: np.ndarray) -> np.ndarray:
    return np.eye(perm.size)[perm]


def _split():
    """``A = D − αP`` (leading invertible diagonal, non-symmetric gain)."""
    D = DiagonalOperator(_D_DIAG)
    N = ScaledOperator(_ALPHA, PermutationOperator(_P4))
    A = D - N
    A_dense = np.diag(_D_DIAG) - _ALPHA * _perm_matrix(_P4)
    return A, A_dense


def _green_of(A) -> GreenOperator:
    """``A.inverse()`` narrowed — the §24.4 factory-dispatch pin at every use."""
    green = A.inverse()
    assert isinstance(green, GreenOperator)  # plain sum → Green (not Sweep)
    return green


def test_g_i1_roundtrip_anchor_and_kernel_columns():
    """G-I1 both ways at driver tol + the dense-LU anchor, INCLUDING the
    δ_j basis (G-kernel fold: ``G.apply(δ_j)`` = column j of the inverse =
    the unit-point-source response)."""
    A, A_dense = _split()
    green = _green_of(A)
    x = _RNG.standard_normal(4)
    np.testing.assert_allclose(
        green.apply(A.apply(x)), x, rtol=_SAFETY * green.tol,
        err_msg="G-I1: A⁻¹(Ax) ≠ x",
    )
    np.testing.assert_allclose(
        A.apply(green.apply(_Q)), _Q, rtol=_SAFETY * green.tol,
        err_msg="G-I1: A(A⁻¹q) ≠ q",
    )
    inputs = [_Q] + [np.eye(4)[j] for j in range(4)]  # random + δ_j basis
    for j, q in enumerate(inputs):
        np.testing.assert_allclose(
            green.apply(q), np.linalg.solve(A_dense, q),
            rtol=_SAFETY * green.tol, atol=_SAFETY * green.tol,
            err_msg=f"Green ≠ dense LU of the materialized sum (input {j})",
        )


def test_g_i2_involution_object_identity():
    """``(G)⁻¹`` is the SUM ITSELF, by object identity (mixin back-half)."""
    A, _ = _split()
    assert A.inverse().inverse() is A


def test_back_half_solve_is_the_forward_matvec():
    """The mixin's ``solve`` on the INVERSE object is the FORWARD action —
    anchored against the dense sum (independent), NOT via ``sum.apply``
    (the mixin back-half is DEFINED as ``inner.apply`` — asserting that
    delegation against itself would be the definition-tautology; note
    ``OperatorSum.solve`` itself is retired, carve P4).  One sibling suffices:
    the line is single-sourced on the mixin, so a delegation bug here is
    a bug on all three (M-MIXIN-CODOM/solve teeth)."""
    A, A_dense = _split()
    green = _green_of(A)
    x = _RNG.standard_normal(4)
    np.testing.assert_allclose(
        green.solve(x), A_dense @ x, rtol=1e-12,
        err_msg="inverse.solve ≠ the forward matvec (mixin back-half)",
    )


def test_g_neumann_partial_sums():
    """THE name-earner (§21.1, the §17 falsifier-5 shape run as a gate):
    the Neumann partial sums :math:`\\sum_k (D^{-1}N)^k D^{-1} q` converge
    to ``green.apply(q)`` AND to the independent dense LU, with the EXACT
    4-step geometric decay :math:`\\rho = \\alpha/(\\prod d_i)^{1/4}` (the
    permutation's spectrum makes per-step norms rotate, so the pin is the
    exact 4-cycle contraction, not a fuzzy successive-ratio band)."""
    A, A_dense = _split()
    green = _green_of(A)
    D_inv = DiagonalOperator(_D_DIAG).inverse()
    N = ScaledOperator(_ALPHA, PermutationOperator(_P4))

    term = D_inv.apply(_Q)          # k = 0
    acc = term.copy()
    norms = [float(np.linalg.norm(term))]
    for _k in range(1, 40):
        term = D_inv.apply(N.apply(term))   # (D⁻¹N)^k D⁻¹ q
        acc = acc + term
        norms.append(float(np.linalg.norm(term)))

    np.testing.assert_allclose(
        acc, green.apply(_Q), rtol=_SAFETY * green.tol,
        err_msg="Neumann partial sums ≠ Green application (G-Neumann)",
    )
    np.testing.assert_allclose(
        acc, np.linalg.solve(A_dense, _Q), rtol=1e-10,
        err_msg="Neumann partial sums ≠ dense LU (independent oracle)",
    )
    # (D⁻¹αP)⁴ = α⁴/∏d · I exactly ⇒ the 4-step norm ratio IS ρ⁴.
    four_step = (norms[20] / norms[16]) ** 0.25
    np.testing.assert_allclose(
        four_step, _RHO, rtol=1e-10,
        err_msg="Neumann decay ratio ≠ the splitting's spectral radius",
    )


def test_g_reciprocity_euclidean_transposed_operands():
    """G-reciprocity (§21.3): ``⟨φ₂, Gφ₁⟩ = ⟨Gᵀφ₂, φ₁⟩`` with ``Gᵀ`` built
    MANUALLY as the Green over the TRANSPOSED operands (Dᵀ = D diagonal;
    (αP)ᵀ = αP⁻¹, the argsort permutation) — not ``.H`` (#280).  Pins the
    split derivation on a DIFFERENT operand config without a second dense
    oracle; ``A`` is non-symmetric so the identity is not vacuous."""
    A, _ = _split()
    A_T = DiagonalOperator(_D_DIAG) - ScaledOperator(
        _ALPHA, PermutationOperator(np.argsort(_P4)),
    )
    gA, gAT = _green_of(A), _green_of(A_T)
    for seed in (1, 2, 3):
        rng = np.random.default_rng(seed)
        phi1, phi2 = rng.standard_normal(4), rng.standard_normal(4)
        lhs = float(np.dot(phi2, gA.apply(phi1)))
        rhs = float(np.dot(gAT.apply(phi2), phi1))
        np.testing.assert_allclose(
            lhs, rhs, rtol=_SAFETY * gA.tol,
            err_msg="Green reciprocity ⟨φ₂,Gφ₁⟩ ≠ ⟨Gᵀφ₂,φ₁⟩",
        )


def test_three_term_left_spine_flatten():
    """``(D − N₁) − N₂`` flattens to precond D + gains [N₁, N₂] — anchored
    against the 3-term dense LU (the M-GRN-FLATTEN L0 leg)."""
    D = DiagonalOperator(_D_DIAG)
    N1 = ScaledOperator(1.0, PermutationOperator(_P4))
    N2 = DiagonalOperator(np.full(4, 0.5))
    A = (D - N1) - N2
    A_dense = (
        np.diag(_D_DIAG) - _perm_matrix(_P4) - np.diag(np.full(4, 0.5))
    )
    green = _green_of(A)
    np.testing.assert_allclose(
        green.apply(_Q), np.linalg.solve(A_dense, _Q),
        rtol=_SAFETY * green.tol,
        err_msg="3-term left-spine flatten produced the wrong splitting",
    )


def test_divergent_split_raises_loudly_convergent_control_does_not():
    """§22.2: ``ρ(D⁻¹N) > 1`` (the C+L trap's L0 analog) raises
    :class:`ConvergenceFailure`; the convergent control MUST NOT raise (a
    bare raises-gate with no control is Mode-10 blind)."""
    D = DiagonalOperator(_D_DIAG)
    div = D + ScaledOperator(8.0, PermutationOperator(_P4))  # ρ ≈ 1.49
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        with pytest.raises(ConvergenceFailure, match="residual"):
            GreenOperator(div, max_iter=200).apply(_Q)
    conv, _ = _split()
    _green_of(conv).apply(_Q)  # control: no raise


def test_near_critical_increment_lies_true_residual_governs():
    """§18.A (Signature 9, ρ-blind stopping): at ``ρ = 0.99`` the driver's
    INCREMENT converges below tol while the TRUE equation residual is
    ~``ρ/(1−ρ) ≈ 1e2`` larger.  A tight budget must RAISE on the true
    residual (M-GRN-INCREMENT's gate — an increment-checking Green returns
    the wrong iterate silently); a generous budget lets the refinement
    loop MEET the true promise (the positive control that the promise is
    driven, not merely checked — a check-only design would falsely raise
    for every ρ > 1/2)."""
    n = 6
    D = DiagonalOperator(np.full(n, 1.0))
    N = DiagonalOperator(np.full(n, 0.99))
    q = _RNG.uniform(1.0, 2.0, n)
    A_dense = np.diag(np.full(n, 1.0)) - np.diag(np.full(n, 0.99))

    # Budget calibration (M-GRN-INCREMENT teeth): at ρ=0.99/tol=1e-4 the
    # increment-stop lands at ~460 steps (0.99ⁿ < 1e-2·tol/1e-4) and the
    # true-residual refinement closes at ~920, so a 600 budget sits
    # BETWEEN them — the honest Green raises (promise unmet at ~2e-3),
    # while an increment-checking mutant returns the wrong iterate
    # silently at ~460 (this gate then reds on the missing raise).
    with pytest.raises(ConvergenceFailure, match="TRUE relative residual"):
        GreenOperator(D - N, max_iter=600, tol=1e-4).apply(q)

    psi = GreenOperator(D - N, max_iter=5000, tol=1e-4).apply(q)
    true_res = np.linalg.norm(A_dense @ psi - q) / np.linalg.norm(q)
    assert true_res < 1e-4, (
        f"refinement returned with the promise unmet: {true_res:.3e}"
    )
    np.testing.assert_allclose(
        psi, np.linalg.solve(A_dense, q), rtol=1e-2,
        err_msg="near-critical refinement diverged from the LU anchor",
    )


def test_green_threads_initial_guess_to_driver_start(monkeypatch):
    """§23 Mode-11 pin: ``Green.apply(q, initial_guess=x₀)`` threads x₀ as
    the INTERNAL driver's start.  The landed §14 spy (a level below —
    per-iterate threading inside the driver) is BLIND to this drop, and a
    warm start is value-invisible (the converged answer is
    seed-independent), so this spy is the only catcher (M-GRN-SEED)."""
    A, _ = _split()
    green = _green_of(A)
    x0 = _RNG.standard_normal(4)
    seen = []
    real_solve = SourceIteration.solve

    def spy(self, q_ext, initial_guess=None):
        seen.append(initial_guess)
        return real_solve(self, q_ext, initial_guess=initial_guess)

    monkeypatch.setattr(SourceIteration, "solve", spy)
    green.apply(_Q, initial_guess=x0)
    if not seen or seen[0] is None:
        pytest.fail("Green dropped initial_guess to the driver start")
    np.testing.assert_array_equal(
        seen[0], x0,
        err_msg="Green threaded something other than the caller's start",
    )


def test_leading_non_invertible_refuses_at_construction():
    """§22.3: the canonical-ordering refusal — a sum whose LEADING term is
    not invertible cannot designate a preconditioner; the message names the
    canonical spelling (the §18.B contract pin)."""
    with pytest.raises(NotInvertible, match="canonical ordering"):
        (ZeroOperator() + DiagonalOperator(_D_DIAG)).inverse()


def test_tol_zero_is_unsatisfiable():
    """``tol=0`` demands an exact iterative inverse — always raises (the
    documented edge; a silent best-effort return would be a wrong answer)."""
    A, _ = _split()
    with pytest.raises(ConvergenceFailure):
        GreenOperator(A, max_iter=50, tol=0.0).apply(_Q)
