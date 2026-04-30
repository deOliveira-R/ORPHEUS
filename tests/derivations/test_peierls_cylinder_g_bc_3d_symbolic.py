r"""Paired symbolic-vs-production contract test for the 3-D
:math:`G_{\rm bc}^{\rm cyl}` cylinder white-BC kernel.

Math-origin pattern: the SymPy derivation in
:mod:`orpheus.derivations.peierls_cylinder_g_bc_3d` is the **source of
truth** for the corrected 3-D :math:`G_{\rm bc}^{\rm cyl}` form
(Issue #112 Phase C). This test pins:

1. The :math:`z \to \mathrm{Ki}_2` reduction
   (:eq:`peierls-cyl-Gbc-3d-derivation` derivation).
2. The numerical match between the surface-centric reference
   evaluator :func:`G_bc_cyl_correct` and the production
   observer-centric :func:`compute_G_bc_cylinder_3d`.
3. The thin-cell limits :math:`r \to 0,\,\Sigma_t \to 0`:
   :math:`G_{\rm correct}(0) = 4\,\mathrm{Ki}_2(0) = 4` and
   :math:`G_{\rm legacy}(0) = 2\,\mathrm{Ki}_1(0) = \pi`.
4. A frozen quantification of the legacy/correct discrepancy at
   :math:`(r=0,\,\Sigma_t R = 1)` to lock the historical
   3-bug fix permanently.

The verified equation is :eq:`peierls-cyl-Gbc-3d-final` (cylinder
mode-0 :math:`G_{\rm bc}^{\rm cyl}` final form, Issue #112 Phase C).

The canonical derivation source is
``orpheus/derivations/peierls_cylinder_g_bc_3d.py:1`` (full SymPy
docstring with the :math:`z` → :math:`\mathrm{Ki}_2` reduction) and
the shipped production primitive lives at
``orpheus/derivations/peierls_geometry.py:1725``
(:func:`compute_G_bc_cylinder_3d`).
"""
from __future__ import annotations

import numpy as np
import pytest
import sympy as sp

from orpheus.derivations._kernels import ki_n_mp
from orpheus.derivations.peierls_cylinder_g_bc_3d import (
    G_bc_cyl_correct,
    G_bc_cyl_legacy_three_bug,
    Ki_n,
)
from orpheus.derivations.peierls_geometry import (
    CYLINDER_1D,
    compute_G_bc_cylinder_3d,
)


# ═══════════════════════════════════════════════════════════════════════
# Algebraic reduction: ∫ exp(-σd_3D)/d_3D³ dz · cosθ · 2πR ⇒ Ki_2 form
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.l1
@pytest.mark.verifies("peierls-cyl-Gbc-3d-final")
def test_z_integral_to_Ki2():
    r"""Symbolic identity (the load-bearing reduction):

    .. math::

       \int_{-\infty}^{+\infty} \frac{e^{-\Sigma_t\,d_{\rm 3D}}}
                                       {d_{\rm 3D}^3}\,\mathrm dz
       \;=\; \frac{2}{d^2}\,\mathrm{Ki}_2(\Sigma_t\,d).

    Verified by SymPy substitution :math:`z = d \tan\alpha` and the
    Bickley convention :math:`\mathrm{Ki}_2(x) = \int_0^{\pi/2}
    \cos\alpha\,e^{-x/\cos\alpha}\,\mathrm d\alpha`. Numerically
    checks that the LHS evaluated by ``scipy.integrate.quad`` matches
    :math:`(2/d^2)\,\mathrm{Ki}_2(\Sigma_t d)` from
    :func:`ki_n_mp`.
    """
    from scipy.integrate import quad

    sig_t = 1.0
    for d in (0.1, 0.5, 1.0, 3.0):
        def integrand(z, d=d):
            d_3D = (d * d + z * z) ** 0.5
            return float(np.exp(-sig_t * d_3D) / d_3D ** 3)

        # quad on (-inf, inf): use symmetry, 2 * integral on (0, inf)
        val_z, _ = quad(integrand, 0.0, np.inf, epsabs=1e-13,
                        epsrel=1e-12)
        val_z *= 2.0  # symmetry
        val_ki = (2.0 / d / d) * float(ki_n_mp(2, sig_t * d, 30))
        assert val_z == pytest.approx(val_ki, rel=1e-9), (
            f"z-integral identity fails at d={d}: "
            f"quad={val_z:.16e}, Ki2-form={val_ki:.16e}"
        )


# ═══════════════════════════════════════════════════════════════════════
# Production parity — surface-centric ↔ observer-centric ↔ same number
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.l1
@pytest.mark.verifies("peierls-cyl-Gbc-3d-final")
def test_correct_form_matches_production():
    r"""The surface-centric reference :func:`G_bc_cyl_correct`
    (:eq:`peierls-cyl-Gbc-3d-derivation`) and the production
    observer-centric :func:`compute_G_bc_cylinder_3d`
    (:eq:`peierls-cyl-Gbc-3d-final`) are algebraically equivalent —
    related by a change of variables between surface azimuth
    :math:`\phi` and observer azimuth :math:`\psi`. Their numerical
    values must agree to within quadrature precision at any homogeneous
    cylinder fixture. This test pins that equivalence — any future
    deviation indicates a regression in either kernel.
    """
    radii = np.array([5.0])
    sig_t = np.array([1.0])
    r_nodes = np.array([0.5, 2.5, 4.5])

    g_prod = compute_G_bc_cylinder_3d(
        CYLINDER_1D, r_nodes, radii, sig_t,
        n_surf_quad=64, dps=25,
    )
    g_ref = np.array([
        G_bc_cyl_correct(float(r), float(radii[0]), float(sig_t[0]),
                          n_quad=256, dps=25)
        for r in r_nodes
    ])
    np.testing.assert_allclose(
        g_prod, g_ref, rtol=1e-10, atol=1e-12,
        err_msg=(
            "compute_G_bc_cylinder_3d (observer-centric) drifted from "
            "G_bc_cyl_correct (surface-centric reference)."
        ),
    )


# ═══════════════════════════════════════════════════════════════════════
# Thin-cell asymptotics — locks the 3-bug discrepancy
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.l1
@pytest.mark.verifies("peierls-cyl-Gbc-3d-final")
def test_thin_cell_limit_correct_4():
    r"""At :math:`r \to 0,\,\Sigma_t \to 0,\,R = 1`, the correct
    form gives :math:`G_{\rm correct}(0) \to 4\,\mathrm{Ki}_2(0) = 4`.

    Sanity check that :math:`\mathrm{Ki}_2(0) = 1`
    (closed form by the Wallis integral).
    """
    assert Ki_n(2, 0.0) == pytest.approx(1.0, abs=1e-12)
    g0 = G_bc_cyl_correct(1e-9, 1.0, 1e-9, n_quad=64, dps=25)
    assert g0 == pytest.approx(4.0, rel=1e-3), (
        f"G_correct(r→0, σR→0) should approach 4, got {g0}"
    )


@pytest.mark.l1
@pytest.mark.verifies("peierls-cyl-Gbc-3d-final")
def test_legacy_form_thin_cell_pi():
    r"""At :math:`r \to 0,\,\Sigma_t \to 0,\,R = 1`, the legacy
    3-bug form gives :math:`G_{\rm legacy}(0) \to 2\,\mathrm{Ki}_1(0)
    = \pi`.

    Sanity check that :math:`\mathrm{Ki}_1(0) = \pi/2`
    (Wallis integral with :math:`n - 1 = 0`).
    """
    assert Ki_n(1, 0.0) == pytest.approx(np.pi / 2, rel=1e-12)
    g0 = G_bc_cyl_legacy_three_bug(1e-9, 1.0, 1e-9, n_quad=64, dps=25)
    assert g0 == pytest.approx(np.pi, rel=1e-3), (
        f"G_legacy(r→0, σR→0) should approach π, got {g0}"
    )


@pytest.mark.l1
@pytest.mark.verifies("peierls-cyl-Gbc-3d-final")
def test_three_bug_diff_quantified():
    r"""Frozen documentation of the historical 3-bug fix at
    :math:`(r \approx 0,\,\Sigma_t R = 1)`:
    :math:`G_{\rm legacy}/G_{\rm correct} \approx 0.60` — i.e. the
    legacy form **under-estimated** :math:`G_{\rm bc}^{\rm cyl}` by
    ~40 % at this configuration. The empirical ratio is the witness
    that the three-bug fix landed correctly; any future regression
    that re-introduces any of the three bugs (Bickley order, geometry
    factor, leading coefficient) will trip this gate.

    Reference: see ``peierls_unified.rst:peierls-cyl-Gbc-3d-final``
    section for the cylinder 1G/1R row-sum probe table that drove
    the Phase C fix.
    """
    cor = G_bc_cyl_correct(1e-3, 1.0, 1.0, n_quad=128, dps=25)
    leg = G_bc_cyl_legacy_three_bug(1e-3, 1.0, 1.0, n_quad=128, dps=25)
    ratio = leg / cor
    # Expected ~0.60; allow ±5 % envelope around the historical value.
    assert 0.55 < ratio < 0.65, (
        f"legacy/correct ratio drifted: got {ratio:.4f}, expected ~0.60. "
        f"This may indicate one of the three bugs has been re-introduced."
    )
