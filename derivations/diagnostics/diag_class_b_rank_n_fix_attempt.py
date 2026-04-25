"""Test the Issue #132 fix: drop (ρ_max/R)² Jacobian from compute_P_esc_mode.

Per the SymPy derivation in
``derivations/peierls_class_b_rank_n_normalization.py`` the canonical
form is:

    P_esc^{(n)}(r_i) = (1/2) · ∫_0^π exp(-τ) · P̃_n(µ_s) · sin θ dθ

i.e. NO (ρ_max/R)² Jacobian — the area Jacobian dA = ρ²/|µ_s|·dΩ exactly
cancels with the µ_s in the surface basis AND the ρ² in the collision
kernel g = exp(-τ)/(4π ρ²).

This script patches ``compute_P_esc_mode`` in-place to drop the Jacobian,
then runs the failing sphere 1G/1R + 1G/2R rank-N table to verify:

1. Sphere 1R rank-2: should stay near -1.10 % (the legacy mode-0 path
   used compute_P_esc which has no Jacobian, so existing rank-2 already
   benefits from a Jacobian-free mode-0; the rank-2 incremental was
   the *only* path with the spurious Jacobian, and removing it should
   improve, not degrade).
2. Sphere 2R rank-2: should drop from +57 % toward the rank-N convergence
   limit (a few %).
3. Sphere 2R rank-N higher modes: should converge to k_inf within
   a few % (no plateau at +66 %).

If this fixes the catastrophe AND preserves rank-1 + 1R rank-2 behaviour,
ship as the Issue #132 fix.
"""

from __future__ import annotations

import time

import numpy as np

from orpheus.derivations import cp_cylinder, cp_sphere
from orpheus.derivations._kernels import _shifted_legendre_eval
from orpheus.derivations._xs_library import LAYOUTS, get_xs
from orpheus.derivations.peierls_geometry import (
    CYLINDER_1D,
    SPHERE_1D,
    gl_float,
    solve_peierls_mg,
)
from orpheus.derivations import peierls_geometry as _pg


_QUAD_BASE = dict(
    n_panels_per_region=2, p_order=3,
    n_angular=24, n_rho=24, n_surf_quad=24, dps=15,
)


def compute_P_esc_mode_fixed(
    geometry,
    r_nodes,
    radii,
    sig_t,
    n_mode,
    n_angular=32,
    dps=25,
):
    """Patched compute_P_esc_mode: drop (ρ_max/R)² Jacobian.

    Identical to the production
    ``peierls_geometry.compute_P_esc_mode`` (peierls_geometry.py:2642-
    2747) EXCEPT the ``jacobian`` factor is removed from the integrand.
    Per the SymPy derivation in
    ``derivations/peierls_class_b_rank_n_normalization.py`` the canonical
    observer-centred form has no surface-to-observer Jacobian — the
    µ_s weight in the surface basis and the ρ² in the collision kernel
    cancel the area Jacobian, leaving only the Lambertian dΩ measure
    × P̃_n(µ_s).
    """
    if geometry.kind == "slab-polar" and n_mode > 0:
        raise NotImplementedError("slab-polar rank-N not in scope")

    r_nodes = np.asarray(r_nodes, dtype=float)
    radii = np.asarray(radii, dtype=float)
    sig_t = np.asarray(sig_t, dtype=float)
    R = float(radii[-1])

    omega_low, omega_high = geometry.angular_range
    omega_pts, omega_wts = gl_float(n_angular, omega_low, omega_high, dps)
    cos_omegas = geometry.ray_direction_cosine(omega_pts)
    angular_factor = geometry.angular_weight(omega_pts)
    pref = geometry.prefactor

    N = len(r_nodes)
    P = np.zeros(N)
    for i in range(N):
        r_i = r_nodes[i]
        total = 0.0
        for k in range(n_angular):
            cos_om = cos_omegas[k]
            rho_max_val = geometry.rho_max(r_i, cos_om, R)
            if rho_max_val <= 0.0:
                continue
            tau = geometry.optical_depth_along_ray(
                r_i, cos_om, rho_max_val, radii, sig_t,
            )
            K_esc = geometry.escape_kernel_mp(tau, dps)
            mu_exit = (rho_max_val + r_i * cos_om) / R
            p_tilde = float(
                _shifted_legendre_eval(n_mode, np.array([mu_exit]))[0]
            )
            # FIX2 (per literature memory): drop (ρ/R)², insert µ_s
            # weight. The literature-researcher's hollow-sphere diagnosis
            # (.claude/agent-memory/literature-researcher/sanchez_mccormick_rank_n_per_face.md
            # §6c) said the bug is a missing (Ω·n) = cos θ = µ_s
            # µ-weight in the partial-current moment.
            total += (
                omega_wts[k] * angular_factor[k]
                * mu_exit * p_tilde * K_esc
            )
        P[i] = pref * total
    return P


def _build_xs(ng_key, n_regions):
    layout = LAYOUTS[n_regions]
    xs_list = [get_xs(region, ng_key) for region in layout]
    return (
        np.vstack([xs["sig_t"] for xs in xs_list]),
        np.stack([xs["sig_s"] for xs in xs_list], axis=0),
        np.vstack([xs["nu"] * xs["sig_f"] for xs in xs_list]),
        np.vstack([xs["chi"] for xs in xs_list]),
    )


def _solve(geometry, ng_key, n_regions, n_bc_modes):
    cp_module = cp_cylinder if geometry is CYLINDER_1D else cp_sphere
    sig_t, sig_s, nu_sig_f, chi = _build_xs(ng_key, n_regions)
    radii = np.array(cp_module._RADII[n_regions])
    sol = solve_peierls_mg(
        geometry, radii=radii, sig_t=sig_t, sig_s=sig_s,
        nu_sig_f=nu_sig_f, chi=chi,
        boundary="white_rank1_mark", n_bc_modes=n_bc_modes,
        **_QUAD_BASE,
    )
    return sol.k_eff


def _kinf(geometry, ng_key, n_regions):
    cp_module = cp_cylinder if geometry is CYLINDER_1D else cp_sphere
    return cp_module._build_case(ng_key, n_regions).k_inf


def main():
    print("=" * 76)
    print("Issue #132 fix attempt: drop (ρ_max/R)² Jacobian in compute_P_esc_mode")
    print("=" * 76)

    # Capture original
    original_compute_P_esc_mode = _pg.compute_P_esc_mode

    cases = [
        ("sph 1G/1R", SPHERE_1D, "1g", 1, [1, 2, 3, 5, 8]),
        ("sph 1G/2R", SPHERE_1D, "1g", 2, [1, 2, 3, 5, 8]),
        ("sph 2G/2R", SPHERE_1D, "2g", 2, [1, 2, 3]),
        ("cyl 1G/1R", CYLINDER_1D, "1g", 1, [1, 2, 3]),
        ("cyl 1G/2R", CYLINDER_1D, "1g", 2, [1, 2, 3]),
    ]

    results = {}

    for label, geom, ng, nr, ranks in cases:
        kinf = _kinf(geom, ng, nr)
        results[label] = (kinf, [], [])
        print(f"\n--- {label} (k_inf = {kinf:.6f}) ---")
        print(f"  {'rank':>4} {'BEFORE':>14} {'err':>9} {'AFTER':>14} {'err':>9}")
        for n in ranks:
            t0 = time.time()
            # BEFORE (current production)
            _pg.compute_P_esc_mode = original_compute_P_esc_mode
            k_before = _solve(geom, ng, nr, n)
            err_before = (k_before - kinf) / kinf * 100
            # AFTER (patched, no Jacobian)
            _pg.compute_P_esc_mode = compute_P_esc_mode_fixed
            k_after = _solve(geom, ng, nr, n)
            err_after = (k_after - kinf) / kinf * 100
            results[label][1].append(k_before)
            results[label][2].append(k_after)
            print(
                f"  {n:>4} {k_before:>14.10f} {err_before:>+8.2f}% "
                f"{k_after:>14.10f} {err_after:>+8.2f}%   ({time.time()-t0:.1f}s)"
            )

    # Restore original (script-cleanliness)
    _pg.compute_P_esc_mode = original_compute_P_esc_mode

    print()
    print("=" * 76)
    print("SUMMARY")
    print("=" * 76)
    print()
    print("  Sanity gate: sph 1G/1R rank-1 must be UNCHANGED (mode-0 path)")
    print("  Critical gate: sph 1G/2R rank-2 should drop from +57 % to <5 %")
    print("  Convergence gate: sph 1G/2R rank-N should approach k_inf")


if __name__ == "__main__":
    main()
