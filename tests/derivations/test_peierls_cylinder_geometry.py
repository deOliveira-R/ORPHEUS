"""L0 geometry tests for cylindrical Peierls quadrature and ray walker.

These tests land with the Phase-4.2 C2 scaffolding commit and were
retained after the 2026-04-21 ``peierls_cylinder`` shim cleanup (when
the ``optical_depths_pm`` :math:`\\tau^{\\pm}` chord walker and the
``composite_gl_y`` alias were deleted — the strictly-more-general
directed ray walker :meth:`CurvilinearGeometry.optical_depth_along_ray`
subsumes the chord form). Coverage:

1. :func:`composite_gl_r` — composite Gauss–Legendre quadrature on
   :math:`[0, R]` with breakpoints at each annular radius. Tested
   against closed-form integrals on single-region and multi-annulus
   configurations.

2. :meth:`CurvilinearGeometry.optical_depth_along_ray` (cylinder-1d) —
   integrates :math:`\\Sigma_t(\\cdot)` along a directed ray starting
   at radial coordinate :math:`r_{\\rm obs}` in direction
   :math:`\\mu = \\cos\\omega`. Tested against closed-form optical
   paths for a homogeneous cylinder, two-annulus diameter transits,
   chords that skirt the inner annulus, and a reciprocity check that
   traversing a path forward and backward yields the same optical
   depth.

These are the sanity checks that must pass before the Nyström kernel
assembly in :func:`~orpheus.derivations.peierls_geometry.build_volume_kernel`
is trusted.
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.derivations.continuous.peierls.geometry import CYLINDER_1D, composite_gl_r


# ═══════════════════════════════════════════════════════════════════════
# Composite GL radial quadrature
# ═══════════════════════════════════════════════════════════════════════

@pytest.mark.l0
@pytest.mark.verifies("chord-length")
class TestCylinderCompositeRadialGL:
    """Sanity checks for the radial y-quadrature with breakpoints at r_k."""

    def test_integrates_constant_to_R(self):
        """∫₀ᴿ 1 dy = R (machine precision)."""
        radii = np.array([0.4, 1.0])
        y_pts, y_wts, _ = composite_gl_r(radii, n_panels_per_region=4, p_order=6)
        assert np.isclose(y_wts.sum(), radii[-1], rtol=1e-14)

    def test_integrates_linear_to_half_Rsq(self):
        """∫₀ᴿ y dy = R²/2 (machine precision, GL order ≥ 2 exact)."""
        radii = np.array([0.5, 1.0, 2.0])
        y_pts, y_wts, _ = composite_gl_r(radii, n_panels_per_region=3, p_order=4)
        R = radii[-1]
        integral = np.dot(y_wts, y_pts)
        assert np.isclose(integral, 0.5 * R ** 2, rtol=1e-14)

    def test_panel_bounds_partition_y_axis(self):
        """Concatenation of panel bounds covers [0, R] without gaps."""
        radii = np.array([0.4, 0.9, 1.5])
        _, _, panel_bounds = composite_gl_r(
            radii, n_panels_per_region=2, p_order=4,
        )
        assert panel_bounds[0][0] == pytest.approx(0.0, abs=1e-14)
        assert panel_bounds[-1][1] == pytest.approx(radii[-1], rel=1e-14)
        for p0, p1 in zip(panel_bounds[:-1], panel_bounds[1:]):
            assert p0[1] == pytest.approx(p1[0], rel=1e-14)

    def test_breakpoints_hit_each_annular_radius(self):
        """Each annular radius r_k appears as a panel endpoint."""
        radii = np.array([0.3, 0.7, 1.2])
        _, _, panel_bounds = composite_gl_r(
            radii, n_panels_per_region=2, p_order=4,
        )
        endpoints = sorted({pa for pa, _, _, _ in panel_bounds}
                           | {pb for _, pb, _, _ in panel_bounds})
        for rk in radii:
            assert any(abs(ep - rk) < 1e-13 for ep in endpoints), (
                f"Radius {rk} is not a panel endpoint; endpoints = {endpoints}"
            )


# ═══════════════════════════════════════════════════════════════════════
# Directed-ray optical-path walker
#
# The cylinder ray starting at radial coordinate ``r_obs`` in direction
# ``cos_omega`` crosses the annular shell of radius ``r_k`` at ray
# distances ``s`` solving ``(r_obs + s cos ω)² + (s sin ω)² = r_k²``,
# i.e. ``s² + 2 r_obs cos ω · s + (r_obs² − r_k²) = 0``. The walker
# piecewise integrates ``Σ_t(r_mid(s)) · Δs`` across those crossings.
# ═══════════════════════════════════════════════════════════════════════

@pytest.mark.l0
@pytest.mark.verifies("chord-length")
class TestCylinderOpticalDepthAlongRay:
    """Closed-form verification of :meth:`CurvilinearGeometry.optical_depth_along_ray`
    for the cylinder-1d geometry."""

    def test_homogeneous_equals_sigma_times_rho(self):
        r"""Homogeneous 1-region cylinder: :math:`\tau = \Sigma_t\,\rho`
        regardless of observer position or ray direction."""
        R = 1.0
        radii = np.array([R])
        sig_t = np.array([2.3])

        # Interior observer, arbitrary direction, ray within cell.
        for r_obs, cos_omega, rho in [
            (0.0, 1.0, 0.7),    # radial outward from centre
            (0.3, 0.5, 0.4),    # generic interior
            (0.8, -0.9, 0.1),   # grazing inward
        ]:
            tau = CYLINDER_1D.optical_depth_along_ray(
                r_obs, cos_omega, rho, radii, sig_t,
            )
            assert np.isclose(tau, sig_t[0] * rho, rtol=1e-14, atol=1e-14), (
                f"Homogeneous ray τ ≠ σ·ρ at (r_obs={r_obs}, "
                f"cos_ω={cos_omega}, ρ={rho}): got {tau}, "
                f"expected {sig_t[0] * rho}"
            )

    def test_scales_linearly_with_sigma_t(self):
        """Doubling every region's :math:`\\Sigma_t` doubles :math:`\\tau`."""
        radii = np.array([0.4, 1.0])
        sig_t = np.array([0.9, 2.1])
        r_obs, cos_omega, rho = 0.7, -0.6, 1.3

        tau1 = CYLINDER_1D.optical_depth_along_ray(
            r_obs, cos_omega, rho, radii, sig_t,
        )
        tau2 = CYLINDER_1D.optical_depth_along_ray(
            r_obs, cos_omega, rho, radii, 2 * sig_t,
        )
        assert np.isclose(tau2, 2 * tau1, rtol=1e-14, atol=1e-14)

    def test_two_annulus_diameter_transit(self):
        r"""Observer at :math:`r_{\rm obs} = R`, ray in direction
        :math:`\cos\omega = -1` (straight through the centre), distance
        :math:`\rho = 2R`. The ray crosses inner annulus 1, passes
        through :math:`r = 0`, and exits through inner annulus 1 again
        before leaving the outer annulus:

        .. math::

           \tau(R, -1, 2R) = 2\,\Sigma_{t,1}\,r_1
                           + 2\,\Sigma_{t,2}\,(R - r_1).
        """
        r1, R = 0.4, 1.0
        radii = np.array([r1, R])
        sig_t = np.array([0.9, 2.1])

        tau = CYLINDER_1D.optical_depth_along_ray(
            r_obs=R, cos_omega=-1.0, rho=2 * R,
            radii=radii, sig_t=sig_t,
        )
        expected = 2 * sig_t[0] * r1 + 2 * sig_t[1] * (R - r1)
        assert np.isclose(tau, expected, rtol=1e-14, atol=1e-14), (
            f"Diameter transit τ: got {tau}, expected {expected}"
        )

    def test_chord_skirts_inner_annulus(self):
        r"""Chord at impact parameter :math:`y` with :math:`r_1 < y < R`
        never enters the inner annulus. Boundary observer at
        :math:`r_{\rm obs} = R` heading across the full chord
        (direction angle :math:`\omega` such that
        :math:`\sin\omega = y/R`): the chord length is
        :math:`\rho = 2\sqrt{R^{2} - y^{2}}` and lies entirely in the
        outer annulus.

        .. math::

           \tau = \Sigma_{t,2}\,\rho = 2\,\Sigma_{t,2}\sqrt{R^{2} - y^{2}}.
        """
        r1, R = 0.4, 1.0
        radii = np.array([r1, R])
        sig_t = np.array([0.9, 2.1])
        y = 0.7
        # Observer at (R, 0) heading toward (-√(R²−y²), y). The ray's
        # cos_omega (angle to outward radial, which is +x̂) is the dot
        # of the unit direction with +x̂: direction = (-√(R²-y²)/R, y/R)
        # after normalising — cos_omega = -√(R²-y²)/R.
        chord_half = np.sqrt(R ** 2 - y ** 2)
        cos_omega = -chord_half / R  # direction cosine w.r.t. +x̂
        rho = 2 * chord_half

        tau = CYLINDER_1D.optical_depth_along_ray(
            r_obs=R, cos_omega=cos_omega, rho=rho,
            radii=radii, sig_t=sig_t,
        )
        expected = 2 * sig_t[1] * chord_half
        assert np.isclose(tau, expected, rtol=1e-13, atol=1e-14), (
            f"Outer-annulus chord τ: got {tau}, expected {expected}"
        )

    def test_reverse_traversal_equals_forward(self):
        r"""Physical reciprocity: :math:`\tau` accumulated on the chord
        from :math:`r_A = r_{\rm obs}` in direction :math:`\mu` for
        distance :math:`\rho` equals :math:`\tau` on the reverse trip
        from the end point back to :math:`r_A`.

        Parametrisation: the ray at distance :math:`s` along the
        forward trip has radial coordinate
        :math:`r(s) = \sqrt{r_A^{2} + 2 r_A s \mu + s^{2}}`. Its end
        point sits at :math:`r_B = r(\rho)`, and the (signed) direction
        cosine of the original ray measured from :math:`r_B`'s outward
        radial is :math:`\mu' = (r_A \mu + \rho)/r_B`. Reversing the
        ray flips the sign: :math:`\mu'_{\rm back} = -\mu'`.
        """
        radii = np.array([0.4, 1.0])
        sig_t = np.array([0.9, 2.1])
        r_A, mu, rho = 0.3, 0.7, 1.2

        r_B = np.sqrt(r_A ** 2 + 2 * r_A * rho * mu + rho ** 2)
        mu_prime_back = -(r_A * mu + rho) / r_B

        tau_fwd = CYLINDER_1D.optical_depth_along_ray(
            r_A, mu, rho, radii, sig_t,
        )
        tau_bwd = CYLINDER_1D.optical_depth_along_ray(
            r_B, mu_prime_back, rho, radii, sig_t,
        )
        assert np.isclose(tau_fwd, tau_bwd, rtol=1e-13, atol=1e-14), (
            f"Reciprocity broken: forward τ = {tau_fwd}, reverse τ = {tau_bwd}"
        )

    def test_multi_annulus_ray_entering_inner(self):
        r"""Ray from boundary observer at :math:`r_{\rm obs} = R`,
        aimed at the centre (:math:`\cos\omega = -1`), stopped at
        distance :math:`\rho = R - r_1 + \delta` so it just enters the
        inner annulus by :math:`\delta`. Confirms the walker correctly
        picks up the inner :math:`\Sigma_t` for the :math:`\delta`
        overlap.

        .. math::

           \tau = \Sigma_{t,2}\,(R - r_1) + \Sigma_{t,1}\,\delta.
        """
        r1, R = 0.4, 1.0
        radii = np.array([r1, R])
        sig_t = np.array([0.9, 2.1])
        delta = 0.15
        rho = (R - r1) + delta

        tau = CYLINDER_1D.optical_depth_along_ray(
            r_obs=R, cos_omega=-1.0, rho=rho,
            radii=radii, sig_t=sig_t,
        )
        expected = sig_t[1] * (R - r1) + sig_t[0] * delta
        assert np.isclose(tau, expected, rtol=1e-14, atol=1e-14), (
            f"Partial inner-annulus τ: got {tau}, expected {expected}"
        )
