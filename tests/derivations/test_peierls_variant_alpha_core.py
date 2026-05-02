"""Foundation tests for the Phase-2 unified Variant α core primitives.

This file pins ONE foundation-tagged test per primitive in
:mod:`orpheus.derivations.continuous.peierls_greens_function.variant_alpha_core`,
verifying the shared functions reproduce the inlined formulas that
were previously duplicated in the sphere and cylinder solvers.

Per ``algebra-of-record`` discipline these are software invariants
(not L0/L1 verification claims about the underlying physics), hence
the ``@pytest.mark.foundation`` tag and absence of any ``verifies``
annotation.
"""
from __future__ import annotations

import math

import numpy as np
import pytest

from orpheus.derivations.continuous.peierls_greens_function.variant_alpha_core import (
    apply_variant_alpha_closure,
    compute_resolvent_T,
)


@pytest.mark.foundation
def test_resolvent_matches_inlined_formula_scattered_points():
    """``compute_resolvent_T(α, τ)`` reproduces the closed-form
    ``1/(1 - α·exp(-τ))`` to machine precision across a representative
    sweep of :math:`(\\alpha, \\tau)` points spanning the regimes the
    sphere and cylinder solvers exercise.
    """
    # Mix of tau and alpha covering: vacuum-ish (small alpha), partial
    # reflection, near-closed (alpha near 1), and short/long periods.
    alphas = [0.0, 0.1, 0.5, 0.9, 0.99, 1.0]
    taus = [1e-3, 0.1, 0.5, 1.0, 2.5, 10.0]
    for alpha in alphas:
        for tau in taus:
            expected = 1.0 / (1.0 - alpha * math.exp(-tau))
            actual = compute_resolvent_T(alpha, tau)
            assert actual == pytest.approx(expected, rel=0.0, abs=1e-15), (
                f"Resolvent mismatch at alpha={alpha}, tau={tau}: "
                f"expected {expected}, got {actual}"
            )


@pytest.mark.foundation
def test_resolvent_alpha_zero_returns_one():
    """Vacuum BC should give :math:`T = 1` exactly for any
    :math:`\\tau_{\\rm period}`.
    """
    for tau in [0.0, 1e-3, 1.0, 100.0]:
        assert compute_resolvent_T(0.0, tau) == 1.0


@pytest.mark.foundation
def test_closure_alpha_zero_returns_F():
    """At :math:`\\alpha = 0` the closure must collapse to the
    first-leg integral exactly: :math:`\\psi_{\\rm new} = F`.
    """
    F = np.array([1.0, 2.0, 3.5, 0.0])
    B = np.array([0.5, 1.0, 2.0, 5.0])
    tau_first_leg = np.array([0.1, 0.5, 1.0, 2.0])
    tau_period = np.array([0.2, 1.0, 2.0, 4.0])
    out = apply_variant_alpha_closure(F, B, tau_first_leg, tau_period, alpha=0.0)
    np.testing.assert_array_equal(out, F)


@pytest.mark.foundation
def test_closure_matches_inlined_sphere_formula():
    """``apply_variant_alpha_closure`` reproduces the inlined sphere
    formula

    .. code::

        denom = 1.0 - alpha * exp(-sigma_t * L_p)
        psi_surf = alpha * B / denom
        psi_new = F + exp(-sigma_t * L_back) * psi_surf

    at scattered points, scalar-valued.
    """
    rng = np.random.default_rng(42)
    for _ in range(20):
        F = float(rng.uniform(0.1, 5.0))
        B = float(rng.uniform(0.1, 5.0))
        sigma_t = float(rng.uniform(0.1, 2.0))
        L_back = float(rng.uniform(0.1, 5.0))
        L_p = float(rng.uniform(0.1, 5.0))
        alpha = float(rng.uniform(0.0, 1.0))

        # Inlined formula (sphere greens_function.py, pre-refactor).
        denom = 1.0 - alpha * math.exp(-sigma_t * L_p)
        psi_surf = alpha * B / denom
        psi_new_inlined = F + math.exp(-sigma_t * L_back) * psi_surf

        # Unified closure.
        psi_new_unified = apply_variant_alpha_closure(
            F=F, B=B,
            tau_first_leg=sigma_t * L_back,
            tau_period=sigma_t * L_p,
            alpha=alpha,
        )

        assert psi_new_unified == pytest.approx(psi_new_inlined, rel=0.0, abs=1e-15)


@pytest.mark.foundation
def test_closure_matches_inlined_cylinder_formula():
    """``apply_variant_alpha_closure`` reproduces the inlined cylinder
    formula

    .. code::

        exp_period = exp(-sigma_t * L_period_3D)
        denom = 1.0 - alpha * exp_period
        psi_surf = alpha * B / denom
        psi_new = F + exp(-sigma_t * L_first_3D) * psi_surf

    at scattered points spanning the (mu_axial, b) phase-space the
    cylinder solver explores.
    """
    rng = np.random.default_rng(7)
    for _ in range(20):
        F = float(rng.uniform(0.0, 5.0))
        B = float(rng.uniform(0.0, 5.0))
        sigma_t = float(rng.uniform(0.1, 2.0))
        # 3D path lengths can be large for grazing-axial rays.
        L_first_3D = float(rng.uniform(0.1, 50.0))
        L_period_3D = float(rng.uniform(0.1, 50.0))
        alpha = float(rng.uniform(0.0, 1.0))

        exp_period = math.exp(-sigma_t * L_period_3D)
        denom = 1.0 - alpha * exp_period
        psi_surf = alpha * B / denom
        psi_new_inlined = F + math.exp(-sigma_t * L_first_3D) * psi_surf

        psi_new_unified = apply_variant_alpha_closure(
            F=F, B=B,
            tau_first_leg=sigma_t * L_first_3D,
            tau_period=sigma_t * L_period_3D,
            alpha=alpha,
        )

        assert psi_new_unified == pytest.approx(psi_new_inlined, rel=0.0, abs=1e-15)


@pytest.mark.foundation
def test_closure_vectorised_matches_scalar_loop():
    """Vectorised application produces identical results to a
    scalar-by-scalar loop.
    """
    rng = np.random.default_rng(11)
    n = 15
    F = rng.uniform(0.1, 5.0, n)
    B = rng.uniform(0.1, 5.0, n)
    tau_first_leg = rng.uniform(0.1, 5.0, n)
    tau_period = rng.uniform(0.1, 5.0, n)
    alpha = 0.7

    vec = apply_variant_alpha_closure(F, B, tau_first_leg, tau_period, alpha)
    scalar = np.array([
        apply_variant_alpha_closure(
            F[i], B[i], tau_first_leg[i], tau_period[i], alpha,
        )
        for i in range(n)
    ])
    np.testing.assert_array_equal(vec, scalar)
