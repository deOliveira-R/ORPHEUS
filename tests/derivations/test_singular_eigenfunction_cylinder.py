r"""Foundation + L1 tests for the Westfall-Metcalf 1973 singular-eigenfunction
expansion of the bare-critical infinite cylinder.

Test breakdown
--------------

* **Branch-1 SymPy gates** — one ``@pytest.mark.foundation`` test per
  ``derive_*()`` in
  :mod:`orpheus.derivations.continuous.singular_eigenfunction.origins.cylinder_derivations`.
  Eight tests pinning the algebraic identities that justify the
  WM-72 reduction:

  - V_se-cyl.1: dispersion function
    :math:`\Lambda(\nu) = 1 - c\nu\,\mathrm{atanh}(1/\nu)` is identical
    to the slab/sphere form (medium property, geometry-independent).
  - V_se-cyl.2: discrete pseudo-eigenfunction
    :math:`\eta_0(\mu) = c\nu_0^2/(\nu_0^2 - \mu^2)` satisfies WM-72
    Eq. 15. **Note**: catches a typo in printed Eq. 17 (single :math:`\nu_0`
    in numerator → should be :math:`\nu_0^2`).
  - V_se-cyl.3: Bessel-Wronskian identity
    :math:`K_1(z)\,I_0(z) + I_1(z)\,K_0(z) = 1/z` (load-bearing for
    the Eq. 9 integrodifferential reduction).
  - V_se-cyl.4: bare-cylinder reduction (corrected)
    :math:`c_1 = c_2 \implies D(\nu) = 0`, :math:`d_0 = 0`,
    :math:`A(\nu) = B(\nu)` (NOT zero — the iterative Fredholm
    coupling).
  - V_se-cyl.5: bare-cylinder criticality structure (corrected) —
    documents the **Fredholm-iterated** Eqs. 30 ↔ 31 + Eq. 32
    structure with the corrected q-formula
    (:math:`q = (R/\nu)\,K_0(R/\mu)\,I_1(R/\nu) + (R/\mu)\,K_1(R/\mu)\,
    I_0(R/\nu)`; second-term denominator is :math:`\mu`, NOT :math:`R`,
    forced by the Wronskian identity :math:`q(\mu, \mu) = 1`).
  - V_se-cyl.6: discrete eigenfunction normalisation matches WM-72
    Eq. 21d.
  - V_se-cyl.7: bare-cylinder neutron density profile is
    :math:`\rho(r) = J_0(r/u_0)`.
  - V_se-cyl.8: Mitsis-Zweifel singular-subtraction identity for
    Eq. 31 — the load-bearing algebra behind the production solver's
    M_Aφ diagonal handling.

* **Branch-2 production** — bare-cylinder solver tests:

  - ``test_solver_returns_valid_result`` — solver returns finite,
    positive :math:`r_c` and reasonable :math:`u_0`.
  - ``test_dispersion_root_consistency`` — returned :math:`u_0`
    satisfies the dispersion relation
    :math:`1 - c\,u_0\,\mathrm{atan}(1/u_0) = 0`.
  - ``test_scalar_flux_at_origin_is_unity`` — V_se-cyl.7 flux
    reconstruction gives :math:`\rho(0) = J_0(0) = 1`.
  - ``test_scalar_flux_monotone_decreasing`` — the bare-cylinder
    profile is monotone-decreasing on :math:`(0, R_c)` (since
    :math:`R_c < j_{0,1} u_0`).
  - ``test_pseudo_angular_flux_positive`` — discrete-mode
    pseudo-flux :math:`\Phi_1^{(0)}(r, \mu) > 0` everywhere in
    :math:`(0, R_c) \times (0, 1)`.

* **L1 reference-value gate (Sood ``Ua-1-0-CY``)** — solver reproduces
  the published Sood :math:`r_c = 1.72500292` mfp at :math:`c = 1.30`
  to **≤ 1e-5 relative** at :math:`n_{\rm grid} = 24`. The hardened
  WM-72 Fredholm method with Mitsis-Zweifel singular subtraction +
  Lagrangian-derivative diagonal handling reaches ~3e-7 in practice,
  comfortably exceeding the 1e-5 brief target. The Variant α cylinder
  cross-check (already shipped at 8.5e-6 in
  :mod:`tests.derivations.test_trajectory_resolvent_cylinder_xverif_sood2003`)
  is now joined by this WM-72 path as a **second, structurally-
  independent** anchor at the same precision.

* **L1 WM-72 Table II benchmark** — six-configuration parametrized
  test asserting agreement to ≤ 1e-5 relative against WM-72's
  published Table II values for :math:`c \in \{1.05, 1.10, 1.20, 1.30,
  1.40, 2.00\}` at :math:`n_{\rm grid} = 24`. Empirically each test
  case reaches ≤ 4e-7.

References
----------

* Westfall & Metcalf 1973, *Nucl. Sci. Eng.* **52**, 1-11.
* Sood, Forster & Parsons 1999, LA-13511 (Table 13, case ``Ua-1-0-CY``).
"""
from __future__ import annotations

import math

import numpy as np
import pytest
import scipy.special

from orpheus.derivations.continuous.singular_eigenfunction import (
    solve_singular_eigenfunction_cylinder_bare_critical,
)
from orpheus.derivations.continuous.singular_eigenfunction.origins.cylinder_derivations import (
    derive_bare_cylinder_criticality_condition,
    derive_bare_cylinder_reduction,
    derive_bessel_wronskian_identity,
    derive_discrete_eigenfunction_normalization,
    derive_discrete_pseudo_eigenfunction,
    derive_dispersion_function,
    derive_flux_reconstruction_bare_cylinder,
    derive_singular_subtraction_eq31,
)
from orpheus.derivations.continuous.sood_registry import LA13511_CASES


# Suppress numpy / scipy integration warnings from the kernel quadrature
# (the weakly-singular kernel triggers "roundoff" warnings that are
# accuracy-bounded, not bugs).
pytestmark = [
    pytest.mark.filterwarnings(
        "ignore:.*roundoff error.*:scipy.integrate.IntegrationWarning"
    ),
    pytest.mark.filterwarnings("ignore::RuntimeWarning"),
]


# ═══════════════════════════════════════════════════════════════════
# Branch-1 SymPy gates (one test per derive_*() function)
# ═══════════════════════════════════════════════════════════════════


@pytest.mark.foundation
def test_v_se_cyl_1_dispersion_function():
    """V_se-cyl.1: cylindrical dispersion function ≡ slab/sphere form."""
    r = derive_dispersion_function()
    assert r["pass"], f"V_se-cyl.1 FAIL: {r}"


@pytest.mark.foundation
def test_v_se_cyl_2_discrete_pseudo_eigenfunction():
    r"""V_se-cyl.2: η_0(μ) = c·ν_0²/(ν_0² - μ²) satisfies WM-72 Eq. 15.

    This test catches the typo in printed WM-72 Eq. 17: the published
    form ``c·ν_0/(ν_0² - μ²)`` does NOT satisfy Eq. 15; the correct
    form is ``c·ν_0²/(ν_0² - μ²)``.
    """
    r = derive_discrete_pseudo_eigenfunction()
    assert r["pass"], f"V_se-cyl.2 FAIL: {r}"


@pytest.mark.foundation
def test_v_se_cyl_3_bessel_wronskian():
    """V_se-cyl.3: Wronskian K_1·I_0 + I_1·K_0 = 1/z."""
    r = derive_bessel_wronskian_identity()
    assert r["pass"], f"V_se-cyl.3 FAIL: {r}"


@pytest.mark.foundation
def test_v_se_cyl_4_bare_cylinder_reduction():
    r"""V_se-cyl.4 (corrected): c_1 = c_2 ⇒ D(ν) = 0, d_0 = 0,
    A(ν) = B(ν) (NOT zero).

    Catches the original Phase B1 documentation error which claimed
    A=B=0 by omitting Eq. 33's middle (c_1/c_2)·N_2 term.
    """
    r = derive_bare_cylinder_reduction()
    assert r["pass"], f"V_se-cyl.4 FAIL: {r}"


@pytest.mark.foundation
def test_v_se_cyl_5_bare_cylinder_criticality():
    r"""V_se-cyl.5 (corrected): bare-cylinder Fredholm criticality
    structure with the corrected q-formula

    .. math::

       q(\nu, \mu) = (R/\nu)\,K_0(R/\mu)\,I_1(R/\nu)
                  + (R/\mu)\,K_1(R/\mu)\,I_0(R/\nu)

    The original Phase B1 SymPy used :math:`R\,K_1\,I_0` for the
    second term; the Wronskian identity :math:`q(\mu, \mu) = 1`
    forces the corrected :math:`(R/\mu)` denominator. This test pins
    that correction symbolically.
    """
    r = derive_bare_cylinder_criticality_condition()
    assert r["pass"], f"V_se-cyl.5 FAIL: {r}"


@pytest.mark.foundation
def test_v_se_cyl_6_discrete_normalization():
    """V_se-cyl.6: N_0 integral matches WM-72 Eq. 21d."""
    r = derive_discrete_eigenfunction_normalization()
    assert r["pass"], f"V_se-cyl.6 FAIL: {r}"


@pytest.mark.foundation
def test_v_se_cyl_7_flux_reconstruction():
    """V_se-cyl.7: bare-cylinder ρ(r) = J_0(r/u_0)."""
    r = derive_flux_reconstruction_bare_cylinder()
    assert r["pass"], f"V_se-cyl.7 FAIL: {r}"


@pytest.mark.foundation
def test_v_se_cyl_8_singular_subtraction():
    r"""V_se-cyl.8: Mitsis-Zweifel singular subtraction for Eq. 31.

    Verifies the structural identity

    .. math::

       \int_0^1 \mu^2 \eta_{2\nu}(\mu) \Phi'(\mu) d\mu
       = \int_0^1 \frac{c\,\nu^2 [\mu^2\Phi'(\mu) - \nu^2\Phi'(\nu)]}
                       {\nu^2 - \mu^2} d\mu + \nu^2 \Phi'(\nu)

    via the dispersion-identity collapse :math:`(1-\lambda)+\lambda = 1`
    of the PV residue and δ contributions. This is the load-bearing
    algebraic backbone of the production solver's
    :func:`...one_group._build_M_A_phi` diagonal handling.
    """
    r = derive_singular_subtraction_eq31()
    assert r["pass"], f"V_se-cyl.8 FAIL: {r}"


# ═══════════════════════════════════════════════════════════════════
# Branch-2 production gates
# ═══════════════════════════════════════════════════════════════════


@pytest.mark.foundation
def test_solver_returns_valid_result():
    """Solver returns finite, positive r_c and a sensible u_0."""
    res = solve_singular_eigenfunction_cylinder_bare_critical(c=1.30, n_grid=32)
    assert res.converged, "Solver did not converge"
    assert math.isfinite(res.r_c_mfp), "r_c is not finite"
    assert res.r_c_mfp > 0, f"r_c = {res.r_c_mfp} should be positive"
    assert res.u_0 > 0, f"u_0 = {res.u_0} should be positive"
    # u_0 for c=1.30 should be ~0.946.
    assert 0.8 < res.u_0 < 1.0, f"u_0 = {res.u_0} outside expected range"


@pytest.mark.foundation
def test_dispersion_root_consistency():
    """Returned u_0 satisfies 1 - c·u_0·atan(1/u_0) = 0."""
    res = solve_singular_eigenfunction_cylinder_bare_critical(c=1.30, n_grid=32)
    Lam = 1.0 - 1.30 * res.u_0 * math.atan(1.0 / res.u_0)
    assert abs(Lam) < 1e-12, (
        f"Dispersion residual {Lam:.3e} too large; u_0 = {res.u_0}"
    )


@pytest.mark.foundation
def test_scalar_flux_at_origin_is_unity():
    r"""ρ(0) = J_0(0) = 1 by V_se-cyl.7."""
    res = solve_singular_eigenfunction_cylinder_bare_critical(c=1.30, n_grid=32)
    rho_0 = res.compute_scalar_flux(0.0)
    np.testing.assert_allclose(rho_0, 1.0, atol=1e-14)


@pytest.mark.foundation
def test_scalar_flux_monotone_decreasing():
    r"""ρ(r) = J_0(r/u_0) is monotone-decreasing on (0, R_c).

    Since R_c < j_{0,1}·u_0 (the bare-critical radius is always
    smaller than the first zero of J_0 by physics — neutrons leak
    before the diffusion-theory critical limit), the dominant
    eigenfunction is monotone-decreasing across the cylinder.
    """
    res = solve_singular_eigenfunction_cylinder_bare_critical(c=1.30, n_grid=32)
    r_grid = np.linspace(0.0, res.r_c_mfp, 50)
    rho = res.compute_scalar_flux(r_grid)
    assert np.all(np.diff(rho) < 0), (
        "ρ(r) should be monotone-decreasing on (0, R_c)"
    )
    # Surface flux must still be positive (no zero before R_c).
    assert rho[-1] > 0, f"ρ(R_c) = {rho[-1]} should be positive"


@pytest.mark.foundation
def test_pseudo_angular_flux_positive():
    r"""Φ_1^{(0)}(r, μ) > 0 everywhere on (0, R_c) × (0, 1)."""
    res = solve_singular_eigenfunction_cylinder_bare_critical(c=1.30, n_grid=32)
    r_grid = np.linspace(0.01, res.r_c_mfp * 0.99, 10)
    mu_grid = np.linspace(0.05, 0.95, 10)
    psi = res.compute_pseudo_angular_flux(r_grid, mu_grid)
    assert psi.shape == (len(r_grid), len(mu_grid))
    assert np.all(psi > 0), "Pseudo-angular flux should be strictly positive"


@pytest.mark.foundation
def test_solver_rejects_subcritical():
    """Solver raises ValueError for c <= 1."""
    with pytest.raises(ValueError, match="c > 1"):
        solve_singular_eigenfunction_cylinder_bare_critical(c=0.9)
    with pytest.raises(ValueError, match="c > 1"):
        solve_singular_eigenfunction_cylinder_bare_critical(c=1.0)


# ═══════════════════════════════════════════════════════════════════
# L1 reference-value gate — Sood Ua-1-0-CY
# ═══════════════════════════════════════════════════════════════════


@pytest.mark.l1
def test_solver_matches_sood_ua_1_0_cy_to_1e5():
    r"""Sood ``Ua-1-0-CY`` (c=1.30): WM-72 cylinder reproduces the
    published :math:`r_c = 1.72500292` mfp to **≤ 1e-5 relative** at
    :math:`n_{\rm grid} = 24`.

    This is the brief's hardening target. The hardened WM-72 Fredholm
    method with Mitsis-Zweifel singular subtraction reaches ~3e-7
    relative in practice; the 1e-5 assertion has 30× slack to absorb
    platform variation. Together with the Variant α cylinder
    cross-check (already at 8.5e-6 in
    :mod:`tests.derivations.test_trajectory_resolvent_cylinder_xverif_sood2003`),
    Sood ``Ua-1-0-CY`` now has TWO structurally-independent anchors at
    the same precision:

    * Variant α via bouncing-characteristic integration with analytical
      bounce-period summation;
    * WM-72 via singular-eigenfunction expansion with
      Mitsis-Zweifel-subtracted Fredholm coupling.

    Both close the Sood reference to ≤ 1e-5 — a strong V&V duo.
    """
    case = LA13511_CASES["Ua-1-0-CY"]
    truth_mfp = case.truth.critical_dimension_mfp

    res = solve_singular_eigenfunction_cylinder_bare_critical(
        c=1.30, n_grid=24, sigma_t=case.materials[0].SigT[0],
    )
    err_rel = abs(res.r_c_mfp - truth_mfp) / truth_mfp
    assert err_rel < 1.0e-5, (
        f"Sood Ua-1-0-CY: WM-72 R_c = {res.r_c_mfp:.9f} mfp, "
        f"truth = {truth_mfp}, rel err = {err_rel:.3e} > 1e-5."
    )
    # Also check: r_c_cm derived correctly.
    assert res.r_c_cm is not None
    expected_r_c_cm = res.r_c_mfp / case.materials[0].SigT[0]
    np.testing.assert_allclose(res.r_c_cm, expected_r_c_cm, rtol=1e-12)
    # Eq. 32 residual at convergence should be small:
    assert res.criticality_residual < 1e-8, (
        f"Eq. 32 residual at converged R = {res.criticality_residual:.3e}, "
        f"expected ≤ 1e-8."
    )


# ── WM-72 Table II benchmark configurations ──────────────────────────
# Six bare-cylinder critical-radius values from the original Westfall-
# Metcalf 1973 paper (Nucl. Sci. Eng. 52, 1-11), Table II — published
# to 7 significant figures. Sood ``Ua-1-0-CY`` (c=1.30, R=1.72500292)
# refines the c=1.30 entry to 8 figures; we use the Sood value where
# available.

# (c, R_truth_mfp, source)
_WM72_TABLE_II = [
    (1.05, 5.411288, "WM-72 Table II"),
    (1.10, 3.577391, "WM-72 Table II"),
    (1.20, 2.287209, "WM-72 Table II"),
    (1.30, 1.72500292, "Sood Ua-1-0-CY (8-digit refinement of WM-72 c=1.30)"),
    (1.40, 1.396979, "WM-72 Table II"),
    (2.00, 0.668613, "WM-72 Table II"),
]


@pytest.mark.l1
@pytest.mark.parametrize("c, R_truth, source", _WM72_TABLE_II,
                         ids=lambda x: f"c={x}" if isinstance(x, float) else str(x))
def test_wm72_table_ii_six_configurations(c, R_truth, source):
    r"""WM-72 Table II — six bare-cylinder configurations at ≤ 1e-5
    relative agreement at :math:`n_{\rm grid} = 24`.

    Empirically each test case reaches ≤ 4e-7 (matching WM-72's
    7-significant-figure quoted precision with their original 24-GL
    quadrature). The 1e-5 assertion has 25× slack for platform
    variation.
    """
    res = solve_singular_eigenfunction_cylinder_bare_critical(
        c=c, n_grid=24,
    )
    err_rel = abs(res.r_c_mfp - R_truth) / R_truth
    assert err_rel < 1.0e-5, (
        f"{source}: c = {c}, WM-72 R_c = {res.r_c_mfp:.9f} mfp, "
        f"truth = {R_truth}, rel err = {err_rel:.3e} > 1e-5."
    )


@pytest.mark.l1
def test_convergence_with_n_grid():
    r"""WM-72 hardened method achieves spectral-rate convergence at
    Sood ``Ua-1-0-CY``.

    Pins the convergence pattern: errors at increasing :math:`n_{\rm grid}`
    monotone-decrease; at :math:`n = 12` the error is ≤ 1e-5; at
    :math:`n = 24` it reaches ~3e-7. The Mitsis-Zweifel singular
    subtraction + Lagrangian-derivative diagonal is the load-bearing
    machinery: both convergence rate AND absolute floor are 4-6
    orders of magnitude better than the original Phase B1 prototype's
    :math:`O(1/n)` algebraic floor.
    """
    truth = 1.72500292
    errs = {}
    for n in [12, 16, 24, 32]:
        res = solve_singular_eigenfunction_cylinder_bare_critical(c=1.30, n_grid=n)
        errs[n] = abs(res.r_c_mfp - truth)

    # Strict monotone-decrease across every refinement.
    ns = sorted(errs.keys())
    for i in range(len(ns) - 1):
        n_lo, n_hi = ns[i], ns[i + 1]
        assert errs[n_hi] < errs[n_lo], (
            f"Convergence failure: n={n_lo} err={errs[n_lo]:.3e}, "
            f"n={n_hi} err={errs[n_hi]:.3e}; expected monotone-decrease."
        )
    # n=12 must already be below 1e-5 absolute (the brief's target).
    assert errs[12] < 1.0e-5, (
        f"n=12 err={errs[12]:.3e} > 1e-5; the hardened method should "
        f"reach the 1e-5 target at n=12."
    )
    # n=24 must reach ~1e-6 absolute (matches WM-72 Table II precision).
    assert errs[24] < 1.0e-6, (
        f"n=24 err={errs[24]:.3e} > 1e-6; expected 7-digit precision "
        f"as quoted in WM-72 Table II at the same quadrature order."
    )
