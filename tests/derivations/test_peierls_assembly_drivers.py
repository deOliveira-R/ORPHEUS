r"""Foundation tests for the per-observer angular-sweep assembler
drivers (`per_observer_angular_assembly` and
`per_surface_centred_angular_assembly`) in
:mod:`orpheus.derivations.peierls_geometry`.

These drivers collapse the ~18 `_per_obs(r_i)` skeletons that
`compute_P_esc_*` / `compute_G_bc_*` / per-face / per-mode primitives
share. The tests pin the contract:

1. The driver builds the right recipe quadrature for each observer.
2. `integrand_at_node(r_i, cos_om, ang_factor, q)` (or
   `(r_i, cos_phi, q)` for the surface-centred sibling) receives the
   correctly-shaped arrays, with `cos_om` and `ang_factor` derived from
   the geometry's `ray_direction_cosine` / `angular_weight` accessors.
3. The driver returns ``prefactor * q.integrate_array(integrand_arr)``
   per observer node — caller-controlled assembly via the callable.
4. ``angular_range`` overrides ``geometry.angular_range`` when supplied.

The drivers are foundation-tier (software-invariant): no associated
theory equation, no V&V level on the L0–L3 ladder. Bit-equivalence of
the migrated `compute_*` consumers is verified by their existing
regression tests in ``test_peierls_*``.
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.derivations.continuous.peierls.geometry import (
    CYLINDER_1D,
    SPHERE_1D,
    per_observer_angular_assembly,
    per_surface_centred_angular_assembly,
)


# ═══════════════════════════════════════════════════════════════════════
# per_observer_angular_assembly — observer-centred ω-sweep driver
# ═══════════════════════════════════════════════════════════════════════

@pytest.mark.foundation
def test_per_observer_constant_integrand_recovers_angular_range():
    r"""Driver applied to integrand ≡ 1 returns
    ``prefactor · ang_factor.sum * (omega_high − omega_low)``.

    For sphere (``angular_weight = sin θ``, ``angular_range = [0, π]``):
    :math:`\int_0^\pi \sin\theta\,d\theta = 2` per observer, so the
    output is ``prefactor · 2`` for every node.
    """
    radii = np.array([5.0])
    sig_t = np.array([0.5])  # noqa: F841 — for completeness
    r_nodes = np.array([0.5, 2.0, 4.5])

    def _integrand(r_i, cos_om, ang_factor, q):
        return ang_factor  # ∫ ang_factor · 1 dω

    result = per_observer_angular_assembly(
        SPHERE_1D, r_nodes, radii, n_per_panel=16,
        integrand_at_node=_integrand,
        prefactor=1.0,
    )
    assert result.shape == (3,)
    np.testing.assert_allclose(result, 2.0, atol=1e-12)


@pytest.mark.foundation
def test_per_observer_prefactor_scales_linearly():
    """Driver multiplies the integrated result by ``prefactor``."""
    radii = np.array([5.0])
    r_nodes = np.array([1.0, 3.0])

    def _integrand(r_i, cos_om, ang_factor, q):
        return ang_factor

    a = per_observer_angular_assembly(
        SPHERE_1D, r_nodes, radii, n_per_panel=16,
        integrand_at_node=_integrand, prefactor=1.0,
    )
    b = per_observer_angular_assembly(
        SPHERE_1D, r_nodes, radii, n_per_panel=16,
        integrand_at_node=_integrand, prefactor=3.5,
    )
    np.testing.assert_allclose(b, 3.5 * a, atol=1e-13)


@pytest.mark.foundation
def test_per_observer_callback_receives_correctly_shaped_args():
    """The integrand callback receives (r_i, cos_om, ang_factor, q)
    with cos_om, ang_factor as ndarrays of length len(q), and r_i as a
    scalar matching the observer node it's evaluating."""
    radii = np.array([2.0, 5.0])  # multi-region → kink subdivision
    r_nodes = np.array([3.0])  # single observer in outer shell
    seen = {}

    def _integrand(r_i, cos_om, ang_factor, q):
        seen["r_i"] = r_i
        seen["cos_om_shape"] = cos_om.shape
        seen["ang_factor_shape"] = ang_factor.shape
        seen["q_len"] = len(q)
        seen["q_panel_count"] = q.n_panels
        return ang_factor

    per_observer_angular_assembly(
        SPHERE_1D, r_nodes, radii, n_per_panel=8,
        integrand_at_node=_integrand,
    )
    assert seen["r_i"] == 3.0
    assert seen["cos_om_shape"] == seen["ang_factor_shape"]
    assert seen["cos_om_shape"][0] == seen["q_len"]
    # r_obs=3.0 with interior shell r=2.0: two tangent angles
    # (forward and backward), so 3 panels.
    assert seen["q_panel_count"] == 3


@pytest.mark.foundation
def test_per_observer_angular_range_override():
    """Supplying ``angular_range`` overrides ``geometry.angular_range``."""
    radii = np.array([5.0])
    r_nodes = np.array([1.0])
    captured = {}

    def _integrand(r_i, cos_om, ang_factor, q):
        captured["interval"] = q.interval
        return ang_factor

    per_observer_angular_assembly(
        SPHERE_1D, r_nodes, radii, n_per_panel=8,
        integrand_at_node=_integrand,
        angular_range=(0.5, 2.5),
    )
    assert captured["interval"] == (0.5, 2.5)


@pytest.mark.foundation
def test_per_observer_matches_hand_coded_skeleton_bit_equally():
    r"""Driver output equals the hand-coded `_per_obs` skeleton for a
    representative integrand. This is the bit-equivalence guard that
    every migrated `compute_*` primitive relies on transitively."""
    from orpheus.derivations.common.quadrature_recipes import (
        observer_angular_quadrature,
    )

    radii = np.array([1.5, 3.0, 5.0])
    sig_t = np.array([0.5, 0.3, 0.1])
    r_nodes = np.array([0.5, 2.0, 4.0])
    R = float(radii[-1])

    # Hand-coded skeleton (the pre-migration pattern).
    def _hand_coded(r_i):
        q = observer_angular_quadrature(
            r_obs=r_i, omega_low=0.0, omega_high=np.pi,
            radii=radii, n_per_panel=12, dps=25,
        )
        cos_om = SPHERE_1D.ray_direction_cosine(q.pts)
        sin_om = SPHERE_1D.angular_weight(q.pts)
        # Synthetic integrand: cos_om² · K_esc(τ) for the longest chord.
        K = np.fromiter(
            (
                float(np.exp(-sig_t[-1] * float(SPHERE_1D.rho_max(r_i, c, R))))
                for c in cos_om
            ),
            dtype=float, count=len(q),
        )
        return 2.0 * q.integrate_array(sin_om * cos_om ** 2 * K)

    expected = np.array([_hand_coded(float(r_i)) for r_i in r_nodes])

    # Driver call with the same integrand.
    def _integrand(r_i, cos_om, ang_factor, q):
        K = np.fromiter(
            (
                float(np.exp(-sig_t[-1] * float(SPHERE_1D.rho_max(r_i, c, R))))
                for c in cos_om
            ),
            dtype=float, count=len(q),
        )
        return ang_factor * cos_om ** 2 * K

    actual = per_observer_angular_assembly(
        SPHERE_1D, r_nodes, radii, n_per_panel=12,
        integrand_at_node=_integrand,
        prefactor=2.0, dps=25,
    )

    np.testing.assert_array_equal(actual, expected)


# ═══════════════════════════════════════════════════════════════════════
# per_surface_centred_angular_assembly — surface-centred φ-sweep driver
# ═══════════════════════════════════════════════════════════════════════

@pytest.mark.foundation
def test_per_surface_centred_constant_integrand_recovers_phi_range():
    r"""Driver with integrand ≡ 1 returns ``prefactor · (phi_high − phi_low)``
    per observer."""
    radii = np.array([5.0])
    r_nodes = np.array([1.0, 3.0])

    def _integrand(r_i, cos_phi, q):
        return np.ones_like(cos_phi)

    result = per_surface_centred_angular_assembly(
        CYLINDER_1D, r_nodes, radii, r_surface=5.0, n_per_panel=16,
        integrand_at_node=_integrand, prefactor=1.0,
    )
    np.testing.assert_allclose(result, np.pi, atol=1e-12)


@pytest.mark.foundation
def test_per_surface_centred_callback_receives_correctly_shaped_args():
    """The integrand callback for the surface-centred driver receives
    (r_i, cos_phi, q) — note the absence of ``ang_factor`` (the
    surface-centred form has no per-geometry angular weight)."""
    radii = np.array([2.0, 5.0])
    r_nodes = np.array([3.0])
    seen = {}

    def _integrand(r_i, cos_phi, q):
        seen["r_i"] = r_i
        seen["cos_phi_shape"] = cos_phi.shape
        seen["q_len"] = len(q)
        seen["q_panel_count"] = q.n_panels
        return np.ones_like(cos_phi)

    per_surface_centred_angular_assembly(
        CYLINDER_1D, r_nodes, radii, r_surface=5.0, n_per_panel=8,
        integrand_at_node=_integrand,
    )
    assert seen["r_i"] == 3.0
    assert seen["cos_phi_shape"][0] == seen["q_len"]
    # r_obs=3.0, r_surface=5.0, interior shell r=2.0: chord-quadratic
    # tangent angles inserted → 3 panels on (0, π).
    assert seen["q_panel_count"] == 3


@pytest.mark.foundation
def test_per_surface_centred_phi_range_override():
    """Supplying ``phi_low/phi_high`` overrides the default ``[0, π]``."""
    radii = np.array([5.0])
    r_nodes = np.array([1.0])
    captured = {}

    def _integrand(r_i, cos_phi, q):
        captured["interval"] = q.interval
        return np.ones_like(cos_phi)

    per_surface_centred_angular_assembly(
        CYLINDER_1D, r_nodes, radii, r_surface=5.0, n_per_panel=8,
        integrand_at_node=_integrand,
        phi_low=0.2, phi_high=2.8,
    )
    assert captured["interval"] == (0.2, 2.8)


@pytest.mark.foundation
def test_per_surface_centred_matches_hand_coded_skeleton_bit_equally():
    r"""Driver output equals the hand-coded surface-centred `_per_obs`
    skeleton for a representative integrand. Bit-equivalence guard for
    the four cylinder G_bc legacy branches that the L3 migration
    routes through this driver."""
    from orpheus.derivations.common.quadrature_recipes import (
        surface_centred_angular_quadrature,
    )

    radii = np.array([2.0, 5.0])
    r_nodes = np.array([0.5, 3.0, 4.5])
    R = float(radii[-1])
    sig_t_val = 0.5

    def _hand_coded(r_i):
        q = surface_centred_angular_quadrature(
            r_obs=r_i, r_surface=R,
            radii=radii, n_per_panel=12, dps=25,
        )
        cos_phi = np.cos(q.pts)
        # Synthetic integrand: 1/d * exp(-sig_t * d)
        kernel = np.empty(len(q))
        for j, c in enumerate(cos_phi):
            d = float(np.sqrt(max(r_i ** 2 + R ** 2 - 2.0 * r_i * R * c, 0.0)))
            kernel[j] = float(np.exp(-sig_t_val * d) / d) if d > 0.0 else 0.0
        return (2.0 / np.pi) * R * q.integrate_array(kernel)

    expected = np.array([_hand_coded(float(r_i)) for r_i in r_nodes])

    def _integrand(r_i, cos_phi, q):
        kernel = np.empty(len(q))
        for j, c in enumerate(cos_phi):
            d = float(np.sqrt(max(r_i ** 2 + R ** 2 - 2.0 * r_i * R * c, 0.0)))
            kernel[j] = float(np.exp(-sig_t_val * d) / d) if d > 0.0 else 0.0
        return kernel

    actual = per_surface_centred_angular_assembly(
        CYLINDER_1D, r_nodes, radii, r_surface=R, n_per_panel=12,
        integrand_at_node=_integrand,
        prefactor=2.0 / np.pi * R, dps=25,
    )

    np.testing.assert_array_equal(actual, expected)
