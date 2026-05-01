r"""Contract tests for :class:`CurvilinearGeometry` (Phase F.1 additions).

Created 2026-04-21 (Phase F.1). This file tests the hollow-core scaffolding
added to :class:`orpheus.derivations.peierls_geometry.CurvilinearGeometry`:

1. ``inner_radius`` field + ``__post_init__`` validation.
2. ``n_surfaces`` property (rank-1 vs rank-2 geometries).
3. ``rho_inner_intersections(r_obs, cos_omega)`` — forward-ρ quadratic
   roots at the inner shell :math:`r = r_0`.
4. ``optical_depth_along_ray`` cavity handling — segments interior to
   :math:`r_0` contribute zero :math:`\tau`.

These are **foundation** tests (software-invariant contracts), not
equation-level verification; no ``verifies()`` label is attached. The
``(None, None)`` sentinel for solid geometry is the regime-A regression
pin: any change that breaks bit-exact recovery of the solid code path
first shows up here.
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.derivations.continuous.peierls.geometry import (
    CYLINDER_1D,
    SLAB_POLAR_1D,
    SPHERE_1D,
    CurvilinearGeometry,
)


# ═══════════════════════════════════════════════════════════════════════
# 1. inner_radius field + validation
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.foundation
class TestInnerRadiusField:
    r"""Contract for the new :attr:`inner_radius` field (Phase F.1)."""

    def test_solid_default_is_zero(self):
        """Existing singletons carry ``inner_radius == 0.0`` by default."""
        assert SLAB_POLAR_1D.inner_radius == 0.0
        assert CYLINDER_1D.inner_radius == 0.0
        assert SPHERE_1D.inner_radius == 0.0

    def test_hollow_cyl_accepts_positive_inner_radius(self):
        g = CurvilinearGeometry(kind="cylinder-1d", inner_radius=0.3)
        assert g.inner_radius == 0.3

    def test_hollow_sph_accepts_positive_inner_radius(self):
        g = CurvilinearGeometry(kind="sphere-1d", inner_radius=0.5)
        assert g.inner_radius == 0.5

    def test_negative_inner_radius_rejected(self):
        with pytest.raises(ValueError, match="inner_radius must be >= 0"):
            CurvilinearGeometry(kind="cylinder-1d", inner_radius=-0.1)
        with pytest.raises(ValueError, match="inner_radius must be >= 0"):
            CurvilinearGeometry(kind="sphere-1d", inner_radius=-0.5)

    def test_slab_with_inner_radius_rejected(self):
        """Slab carries its two faces positionally (x=0, x=L); it must
        not be parameterised with an inner_radius."""
        with pytest.raises(ValueError, match="slab-polar does not carry"):
            CurvilinearGeometry(kind="slab-polar", inner_radius=0.5)

    def test_unsupported_kind_still_rejected(self):
        with pytest.raises(ValueError, match="Unsupported geometry kind"):
            CurvilinearGeometry(kind="annulus-2d")


# ═══════════════════════════════════════════════════════════════════════
# 2. n_surfaces property
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.foundation
class TestNSurfacesProperty:
    r"""The :attr:`n_surfaces` property drives rank-1 vs rank-2 layout
    in the Phase F per-face :class:`BoundaryClosureOperator`."""

    def test_slab_always_two(self):
        """Slab has two faces independent of any radius."""
        assert SLAB_POLAR_1D.n_surfaces == 2

    def test_solid_cyl_one_surface(self):
        assert CYLINDER_1D.n_surfaces == 1

    def test_solid_sph_one_surface(self):
        assert SPHERE_1D.n_surfaces == 1

    def test_hollow_cyl_two_surfaces(self):
        g = CurvilinearGeometry(kind="cylinder-1d", inner_radius=0.2)
        assert g.n_surfaces == 2

    def test_hollow_sph_two_surfaces(self):
        g = CurvilinearGeometry(kind="sphere-1d", inner_radius=0.2)
        assert g.n_surfaces == 2


# ═══════════════════════════════════════════════════════════════════════
# 3. rho_inner_intersections
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.foundation
class TestRhoInnerIntersectionsSolidAndSlab:
    r"""Regime-A sentinel: solid geometries and slabs return
    :math:`(\text{None}, \text{None})` unconditionally."""

    @pytest.mark.parametrize(
        "geometry",
        [SLAB_POLAR_1D, CYLINDER_1D, SPHERE_1D],
    )
    @pytest.mark.parametrize(
        "r_obs, cos_omega",
        [(0.5, 0.5), (1.0, -0.5), (0.1, 0.99)],
    )
    def test_solid_returns_none_pair(self, geometry, r_obs, cos_omega):
        """No hollow shell ⇒ no inner intersections, for any ray."""
        assert geometry.rho_inner_intersections(r_obs, cos_omega) == (None, None)

    def test_slab_with_explicit_zero_inner_radius_still_none(self):
        g = CurvilinearGeometry(kind="slab-polar")
        assert g.rho_inner_intersections(0.5, 0.7) == (None, None)


@pytest.mark.foundation
class TestRhoInnerIntersectionsHollow:
    r"""Closed-form geometric checks for the quadratic roots."""

    def test_radial_inward_ray_single_intersection(self):
        r"""Ray from observer at :math:`r_{\rm obs} = 1` along the radial
        inward direction (:math:`\cos\Omega = -1`) hits the inner shell
        at :math:`\rho^- = r_{\rm obs} - r_0` and re-exits at
        :math:`\rho^+ = r_{\rm obs} + r_0` (through the centre)."""
        g = CurvilinearGeometry(kind="sphere-1d", inner_radius=0.3)
        rho_minus, rho_plus = g.rho_inner_intersections(1.0, -1.0)
        assert rho_minus == pytest.approx(1.0 - 0.3, rel=1e-14)
        assert rho_plus == pytest.approx(1.0 + 0.3, rel=1e-14)

    def test_tangent_ray_has_double_root(self):
        r"""Ray tangent to the inner shell from observer at
        :math:`r_{\rm obs}`: :math:`\sin\Omega_{\rm tan} = r_0/r_{\rm obs}`,
        so :math:`\cos\Omega = -\sqrt{1 - (r_0/r_{\rm obs})^2}` (inward).
        Discriminant is zero ⇒ both roots collapse to
        :math:`\rho = -r_{\rm obs}\cos\Omega`."""
        r_obs = 1.0
        r0 = 0.4
        g = CurvilinearGeometry(kind="cylinder-1d", inner_radius=r0)
        # Inward-tangent direction
        cos_omega = -float(np.sqrt(1.0 - (r0 / r_obs) ** 2))
        rho_minus, rho_plus = g.rho_inner_intersections(r_obs, cos_omega)
        expected = -r_obs * cos_omega  # positive
        assert rho_minus == pytest.approx(expected, rel=1e-12)
        assert rho_plus == pytest.approx(expected, rel=1e-12)

    def test_ray_misses_inner_shell(self):
        """Outward-pointing ray from outside the inner shell: no
        intersection with the cavity."""
        g = CurvilinearGeometry(kind="sphere-1d", inner_radius=0.2)
        # Observer at r=1, outward radial direction
        assert g.rho_inner_intersections(1.0, 1.0) == (None, None)

    def test_ray_parallel_to_axis_misses_when_perpendicular_distance_exceeds_r0(self):
        r"""Cylinder ray with :math:`\cos\Omega = 0` means the ray is
        tangent to the observer's radial circle — equivalently,
        perpendicular distance to the axis equals :math:`r_{\rm obs}`.
        If :math:`r_{\rm obs} > r_0`, the ray stays outside the cavity."""
        g = CurvilinearGeometry(kind="cylinder-1d", inner_radius=0.3)
        assert g.rho_inner_intersections(1.0, 0.0) == (None, None)

    def test_observer_inside_cavity_gives_single_forward_root(self):
        r"""Observer strictly inside the cavity (:math:`r_{\rm obs} < r_0`):
        one of the quadratic roots is behind the observer
        (:math:`\rho < 0`), only the forward root is returned as
        :math:`\rho^+`, :math:`\rho^-` is ``None``.

        (This is a mathematical edge case — the physical solver never
        places observer nodes inside the cavity — but the sentinel must
        still behave correctly for defensive callers.)"""
        g = CurvilinearGeometry(kind="sphere-1d", inner_radius=0.3)
        rho_minus, rho_plus = g.rho_inner_intersections(0.1, 1.0)
        assert rho_minus is None
        assert rho_plus == pytest.approx(0.3 - 0.1, rel=1e-14)

    def test_roots_ordered_rho_minus_le_rho_plus(self):
        """For any forward-hitting ray, the returned pair is ordered."""
        g = CurvilinearGeometry(kind="cylinder-1d", inner_radius=0.25)
        rng = np.random.default_rng(0xF1)
        for _ in range(20):
            r_obs = 0.5 + rng.random()  # (0.5, 1.5)
            cos_omega = 2.0 * rng.random() - 1.0
            rho_minus, rho_plus = g.rho_inner_intersections(r_obs, cos_omega)
            if rho_minus is not None and rho_plus is not None:
                assert rho_minus <= rho_plus


# ═══════════════════════════════════════════════════════════════════════
# 4. optical_depth_along_ray — cavity handling
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.foundation
class TestOpticalDepthCavity:
    r"""Rays that traverse the cavity (:math:`r < r_0`) accumulate τ
    only on the annular-material segments; the cavity segment
    contributes zero.

    Dual-route discipline (plan §9.1): these tests reconstruct τ from
    ``rho_inner_intersections`` + ``rho_max`` by hand algebra and compare
    against ``optical_depth_along_ray``. A Σ_t-sign flip inside the
    cavity segment would be invisible to sum-consistency tests (which
    fire per-surface identically); only hand-algebra reconstruction
    catches it.
    """

    def test_solid_recovery_matches_hollow_at_zero_inner_radius(self):
        """Regime A: ``inner_radius == 0`` reproduces the solid τ byte-exactly."""
        g_solid = CYLINDER_1D
        g_hollow_zero = CurvilinearGeometry(
            kind="cylinder-1d", inner_radius=0.0,
        )
        radii = np.array([1.0])
        sig_t = np.array([1.5])
        rng = np.random.default_rng(0xF2)
        for _ in range(10):
            r_obs = 0.1 + 0.8 * rng.random()
            cos_omega = 2.0 * rng.random() - 1.0
            rho = g_solid.rho_max(r_obs, cos_omega, R=1.0)
            tau_solid = g_solid.optical_depth_along_ray(
                r_obs, cos_omega, rho, radii, sig_t,
            )
            tau_hollow0 = g_hollow_zero.optical_depth_along_ray(
                r_obs, cos_omega, rho, radii, sig_t,
            )
            assert tau_hollow0 == tau_solid

    def test_sphere_radial_through_cavity(self):
        r"""Observer at :math:`r_{\rm obs} = R`, pointing radially
        inward (:math:`\cos\Omega = -1`). Ray hits inner shell at
        :math:`\rho^- = R - r_0`, traverses the cavity to
        :math:`\rho^+ = R + r_0`, and continues out through the
        diametrically opposite annulus until reaching the outer
        boundary at :math:`\rho = 2R`.

        τ over the full chord = Σ_t · (2R − 2r_0) — cavity contributes
        zero, both annular legs contribute Σ_t · (R − r_0)."""
        R, r0 = 1.0, 0.3
        sig_t = 2.0
        g = CurvilinearGeometry(kind="sphere-1d", inner_radius=r0)
        tau = g.optical_depth_along_ray(
            r_obs=R,
            cos_omega=-1.0,
            rho=2.0 * R,
            radii=np.array([R]),
            sig_t=np.array([sig_t]),
        )
        assert tau == pytest.approx(sig_t * (2.0 * R - 2.0 * r0), rel=1e-14)

    def test_cylinder_ray_misses_cavity_matches_solid(self):
        r"""A ray whose perpendicular distance to the axis exceeds
        :math:`r_0` never enters the cavity — τ must match the
        solid-cylinder computation at the same (r_obs, cos_ω)."""
        R, r0 = 1.0, 0.3
        sig_t = np.array([1.5])
        radii = np.array([R])
        g_hollow = CurvilinearGeometry(kind="cylinder-1d", inner_radius=r0)
        # cos_ω = 0 ⇒ perpendicular distance to axis = r_obs = 0.5 > r_0
        r_obs, cos_omega = 0.5, 0.0
        rho = g_hollow.rho_max(r_obs, cos_omega, R)
        tau_hollow = g_hollow.optical_depth_along_ray(
            r_obs, cos_omega, rho, radii, sig_t,
        )
        tau_solid = CYLINDER_1D.optical_depth_along_ray(
            r_obs, cos_omega, rho, radii, sig_t,
        )
        assert tau_hollow == pytest.approx(tau_solid, rel=1e-14)

    def test_cavity_segment_independent_of_sig_t_value(self):
        r"""Doubling Σ_t in the annular material doubles τ on the chord
        — confirming the cavity segment itself contributes zero (else the
        ratio would not be exactly 2). Hand-algebra reconstruction:
        :math:`\tau = \Sigma_t \cdot L_{\rm annulus}` where
        :math:`L_{\rm annulus}` is independent of Σ_t."""
        R, r0 = 1.0, 0.4
        radii = np.array([R])
        g = CurvilinearGeometry(kind="sphere-1d", inner_radius=r0)

        r_obs, cos_omega = 0.8, -0.9  # ray crosses the cavity
        rho = g.rho_max(r_obs, cos_omega, R)
        tau_1 = g.optical_depth_along_ray(
            r_obs, cos_omega, rho, radii, np.array([1.0]),
        )
        tau_2 = g.optical_depth_along_ray(
            r_obs, cos_omega, rho, radii, np.array([2.0]),
        )
        # If cavity contributed nonzero τ, the ratio would deviate from 2
        assert tau_2 == pytest.approx(2.0 * tau_1, rel=1e-14)

    def test_cavity_subtraction_hand_algebra(self):
        r"""Compare hollow-sphere τ to a hand computation built from the
        SOLID τ on the same chord minus the cavity segment's
        Σ_t·(ρ⁺ − ρ⁻) contribution.

        Dual-route catch for cavity-sign bugs: if ``optical_depth_along_ray``
        used ``+ sig_t·(cavity_len)`` instead of ``0``, this would fail.
        """
        R, r0 = 1.0, 0.3
        sig_t_val = 2.5
        radii = np.array([R])
        g_solid = SPHERE_1D
        g_hollow = CurvilinearGeometry(kind="sphere-1d", inner_radius=r0)

        r_obs, cos_omega = 0.9, -0.95  # ray plunges through cavity
        rho = g_solid.rho_max(r_obs, cos_omega, R)
        tau_solid = g_solid.optical_depth_along_ray(
            r_obs, cos_omega, rho, radii, np.array([sig_t_val]),
        )

        rho_in_minus, rho_in_plus = g_hollow.rho_inner_intersections(
            r_obs, cos_omega,
        )
        cavity_length = rho_in_plus - rho_in_minus
        tau_expected = tau_solid - sig_t_val * cavity_length

        tau_hollow = g_hollow.optical_depth_along_ray(
            r_obs, cos_omega, rho, radii, np.array([sig_t_val]),
        )
        assert tau_hollow == pytest.approx(tau_expected, rel=1e-14)
