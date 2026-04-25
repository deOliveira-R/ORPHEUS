"""FIX3 attempt: drop (ρ/R)² Jacobian AND replace reflection_marshak
with identity (drop the Gelbard (2n+1) factors).

Per the SymPy derivation conclusion (HYPOTHESIS path) at
``derivations/peierls_class_b_rank_n_normalization.py``:

  > Replace compute_P_esc_mode integrand with the canonical Lambertian
  > × P̃_n (drop Jacobian) AND
  > Replace reflection_marshak with diag(1) (drop the (2n+1) Gelbard
  > factors) for the rank-N closure on solid Class B — because the new
  > µ-weighted partial-current moments don't need them.

This script monkey-patches BOTH changes simultaneously and re-runs
the failing sphere 1G/2R rank-2 case. Both changes are no-op at n=0
(rank-1 must remain bit-exact).
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


_QUAD = dict(n_panels_per_region=2, p_order=3,
             n_angular=24, n_rho=24, n_surf_quad=24, dps=15)


def compute_P_esc_mode_no_jac(
    geometry, r_nodes, radii, sig_t, n_mode,
    n_angular=32, dps=25,
):
    """Patched: drop (ρ_max/R)² Jacobian. Lambertian × P̃_n form."""
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
            total += omega_wts[k] * angular_factor[k] * p_tilde * K_esc
        P[i] = pref * total
    return P


def reflection_identity(n_modes: int) -> np.ndarray:
    """Replace Gelbard diag(1,3,5,...) with pure identity diag(1,1,...)."""
    return np.eye(n_modes)


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
        **_QUAD,
    )
    return sol.k_eff


def _kinf(geometry, ng_key, n_regions):
    cp_module = cp_cylinder if geometry is CYLINDER_1D else cp_sphere
    return cp_module._build_case(ng_key, n_regions).k_inf


def main():
    print("=" * 76)
    print("Issue #132 FIX3: drop (ρ/R)² Jacobian AND drop Gelbard (2n+1) reflection")
    print("=" * 76)

    orig_P_esc = _pg.compute_P_esc_mode
    orig_reflection_marshak = _pg.reflection_marshak

    cases = [
        ("sph 1G/1R", SPHERE_1D, "1g", 1, [1, 2, 3, 5, 8]),
        ("sph 1G/2R", SPHERE_1D, "1g", 2, [1, 2, 3, 5, 8]),
        ("sph 2G/2R", SPHERE_1D, "2g", 2, [1, 2, 3]),
        ("cyl 1G/1R", CYLINDER_1D, "1g", 1, [1, 2, 3]),
        ("cyl 1G/2R", CYLINDER_1D, "1g", 2, [1, 2, 3]),
    ]

    for label, geom, ng, nr, ranks in cases:
        kinf = _kinf(geom, ng, nr)
        print(f"\n--- {label} (k_inf = {kinf:.6f}) ---")
        print(f"  {'rank':>4} {'BEFORE':>14} {'err':>9} {'FIX3':>14} {'err':>9}")
        for n in ranks:
            t0 = time.time()
            _pg.compute_P_esc_mode = orig_P_esc
            _pg.reflection_marshak = orig_reflection_marshak
            k_before = _solve(geom, ng, nr, n)
            err_before = (k_before - kinf) / kinf * 100
            _pg.compute_P_esc_mode = compute_P_esc_mode_no_jac
            _pg.reflection_marshak = reflection_identity
            k_after = _solve(geom, ng, nr, n)
            err_after = (k_after - kinf) / kinf * 100
            print(
                f"  {n:>4} {k_before:>14.10f} {err_before:>+8.2f}% "
                f"{k_after:>14.10f} {err_after:>+8.2f}%   ({time.time()-t0:.1f}s)"
            )

    _pg.compute_P_esc_mode = orig_P_esc
    _pg.reflection_marshak = orig_reflection_marshak


if __name__ == "__main__":
    main()
