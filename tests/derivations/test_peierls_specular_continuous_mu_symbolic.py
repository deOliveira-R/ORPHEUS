r"""Paired symbolic-vs-textbook contract test for the Phase-5
continuous-:math:`\mu` multi-bounce specular kernel.

Math-origin pattern: the SymPy derivation in
:mod:`orpheus.derivations.peierls_specular.continuous_mu` is the
**source of truth** for the Sanchez 1986 [SanchezTTSP1986]_ Eq. (A6)
↔ ORPHEUS M1 sketch equivalence and the diagonal-singularity finding
that blocks production wiring of ``boundary="specular_continuous_mu"``.

The reference implementation
:func:`compute_K_bc_specular_continuous_mu_sphere` (in
:mod:`peierls_geometry`) is preserved for the textbook reference even
though the production multi-bounce specular path uses
``closure="specular_multibounce"``. This test pins the four SymPy
verifications that gate any future production wiring:

V1. Multi-bounce factor :math:`T(\mu) = 1/(1 - e^{-2a\mu})` has
    :math:`\mu \cdot T(\mu) \to 1/(2a)` as :math:`\mu \to 0^+`
    (finite, by L'Hôpital).
V2. The M1 sketch differs from Sanchez Eq. (A6) by a factor of
    :math:`\mu`; the corrected M1 form (no :math:`\mu` in the
    numerator of the multi-bounce factor) is **algebraically
    equivalent** to Sanchez. This is the µ-weight convention question
    that R1 closure required.
V3. The :math:`2\alpha` BC-prefactor reduces to 0 at
    :math:`\alpha \to 0`, recovering the vacuum kernel.
V4. The diagonal-singularity finding (documentary, non-gating): the
    Sanchez integrand has a :math:`1/\mu^2` singularity at the surface
    diagonal :math:`\rho = a` (non-integrable) and a :math:`1/\mu`
    integrable singularity at interior diagonals.

The verified equation is :eq:`peierls-phase5-sphere-Kbc` (Phase-5
sphere :math:`K_{\rm bc}^{\rm cmu}` formula, Sanchez Eq. A6).

The canonical derivation source is
``orpheus/derivations/peierls_specular/continuous_mu.py:1`` (full
SymPy module docstring with the four verification narratives).
"""
from __future__ import annotations

import pytest
import sympy as sp

from orpheus.derivations.peierls_specular import (
    derive_diagonal_singularity,
    derive_m1_equivalence,
    derive_multi_bounce_factor,
    derive_vacuum_reduction,
)


# ═══════════════════════════════════════════════════════════════════════
# V1 — Multi-bounce factor: µ · T(µ) finite at µ → 0
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.foundation
def test_v1_mu_T_limit_is_one_over_2a():
    r"""V1: :math:`\lim_{\mu \to 0^+} \mu \cdot T(\mu) = 1/(2a)`.

    The Sanchez :math:`T(\mu) = 1/(1 - e^{-2a\mu})` has a simple pole
    at :math:`\mu = 0`; the product :math:`\mu \cdot T(\mu)` is the
    form that appears in the M1 sketch and must be finite to admit
    smooth quadrature.

    Reference: workbench (now lifted) lines 129-158
    (:func:`derive_multi_bounce_factor`).
    """
    result = derive_multi_bounce_factor()
    a = sp.Symbol("a", positive=True)
    expected = sp.Rational(1, 2) / a
    diff = sp.simplify(result["limit"] - expected)
    assert diff == 0, (
        f"V1 fails: lim_{{µ→0+}} µ·T(µ) = {result['limit']}, "
        f"expected 1/(2a); residual {diff}"
    )
    assert result["pass"] is True


# ═══════════════════════════════════════════════════════════════════════
# V2 — Sanchez Eq. (A6) ↔ M1-corrected algebraic equivalence
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.foundation
def test_v2_sanchez_eq_corrected_M1():
    r"""V2: the load-bearing identity of the lift —
    :math:`I_{\rm Sanchez}(\mu) = I_{M1,{\rm corrected}}(\mu)` after
    dropping the spurious :math:`\mu` from the numerator of the
    multi-bounce factor in the M1 sketch.

    The naive M1 sketch carries an extra :math:`\mu` factor in the
    multi-bounce numerator; this is wrong by a factor of
    :math:`1/\mu`. The corrected form (no extra :math:`\mu`)
    is bit-equal to Sanchez Eq. (A6) under the explicit
    :math:`G_{\rm in}, F_{\rm out}` identification.

    Reference: workbench (now lifted) lines 161-283
    (:func:`derive_m1_equivalence`).
    """
    result = derive_m1_equivalence()
    # The corrected diff must be exactly zero
    assert result["diff_corrected"] == 0, (
        f"V2 corrected form not bit-equal to Sanchez: "
        f"diff = {result['diff_corrected']}"
    )
    # The naive ratio is documentary — should be 1/µ
    mu = sp.Symbol("mu", positive=True)
    expected_ratio = 1 / mu
    diff_ratio = sp.simplify(result["ratio_naive"] - expected_ratio)
    assert diff_ratio == 0, (
        f"V2 naive M1 / Sanchez ratio drifted: got "
        f"{result['ratio_naive']}, expected 1/µ"
    )
    assert result["pass"] is True


# ═══════════════════════════════════════════════════════════════════════
# V3 — Vacuum-BC reduction: α → 0 prefactor → 0
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.foundation
def test_v3_alpha_zero_kills_BC_kernel():
    r"""V3: the :math:`2\alpha` BC-prefactor reduces to 0 at
    :math:`\alpha \to 0`, switching off the BC kernel and recovering
    the vacuum form.

    Reference: workbench (now lifted) lines 370-415
    (:func:`derive_vacuum_reduction`).
    """
    result = derive_vacuum_reduction()
    assert result["limit_alpha_zero"] == 0, (
        f"V3 fails: lim_{{α→0}} 2α = {result['limit_alpha_zero']}, "
        f"expected 0"
    )
    assert result["pass"] is True


# ═══════════════════════════════════════════════════════════════════════
# V4 — Diagonal singularity finding (documentary, non-gating)
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.foundation
def test_v4_diagonal_singularity_documented():
    r"""V4 (documentary): the Sanchez integrand on the diagonal
    :math:`\rho' \to \rho` has a leading :math:`1/\mu` form for
    interior :math:`\rho < a` (integrable, log singularity) and a
    :math:`1/\mu^2` form at the surface :math:`\rho = a`
    (non-integrable). This test does **not** gate production — it
    asserts only that the leading-order analysis returns a non-zero
    expression, witnessing the production-blocker that prevents
    wiring ``boundary="specular_continuous_mu"`` directly.

    Future investigators wiring the continuous-µ closure must
    address the singularity (adaptive µ-quadrature, singularity
    subtraction, or a change-of-variables to absorb it). See
    :mod:`orpheus.derivations.peierls_specular.continuous_mu` V4
    docstring for the three resolution options.

    Reference: workbench (now lifted) lines 286-367
    (:func:`derive_diagonal_singularity`).
    """
    result = derive_diagonal_singularity()
    leading = result["leading_at_mu_zero"]
    # The leading expression should not be zero — the singularity is real.
    assert leading != 0, (
        "V4 leading-order expansion is unexpectedly zero — the "
        "diagonal singularity should be present in the integrand."
    )
    # Always passing — this is documentary
    assert result["pass"] is True
