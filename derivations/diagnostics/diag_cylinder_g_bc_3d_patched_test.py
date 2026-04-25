"""Monkey-patch test: does the corrected 3-D G_bc^cyl fix the row-sum bias?

Per derivations/peierls_cylinder_g_bc_3d_derivation.py:

  Correct: G_bc^cyl(r) = (4/π) · ∫_0^π Ki_2(Σ_t · d_2D(r, ψ)) dψ
           with d_2D(r, ψ) = -r cos ψ + sqrt(R² - r² sin² ψ)
                             (backward chord from r in observer direction ψ)

  Current: G_bc^cyl(r) = (2R/π) · ∫_0^π Ki_1(Σ_t · d_surf(r, φ)) / d_surf dφ
           with d_surf(r, φ) = sqrt(r² - 2rR cos φ + R²)
                              (in-plane distance from r to surface point at φ)

Both forms ARE identical at r=0 (limit) but differ at r > 0 by Lambertian
projection factor and Bickley function order. At σ_t·R = 1, current
under-estimates by ~40 %; row-sum K·1/σ_t shows ~11 % under-estimate
empirically.

This script monkey-patches `compute_G_bc` cylinder branch to use the
corrected form, then runs the row-sum probe + Class B k_eff scan to
verify the fix.
"""

from __future__ import annotations

import numpy as np

from orpheus.derivations import cp_cylinder
from orpheus.derivations._kernels import ki_n_mp
from orpheus.derivations._xs_library import LAYOUTS, get_xs
from orpheus.derivations.peierls_geometry import (
    CYLINDER_1D, gl_float, composite_gl_r, build_volume_kernel,
    build_white_bc_correction_rank_n, compute_P_ss_cylinder,
    solve_peierls_mg,
)
from orpheus.derivations import peierls_geometry as _pg


_QUAD_BASE = dict(n_panels_per_region=2, p_order=3,
                  n_angular=24, n_rho=24, n_surf_quad=24, dps=15)


def _compute_G_bc_corrected_cylinder(geometry, r_nodes, radii, sig_t,
                                       n_surf_quad=32, dps=25):
    """Patched G_bc for cylinder: correct 3-D observer-centric form.

    G_bc^cyl(r) = (4/π) · ∫_0^π Ki_2(Σ_t · d_2D(r, ψ)) dψ

    For multi-region σ_t, the 2-D backward chord is integrated piecewise
    via geometry.optical_depth_along_ray, with the backward-direction
    geometry.
    """
    if geometry.kind != "cylinder-1d":
        # Fallback to original for non-cylinder
        from orpheus.derivations.peierls_geometry import compute_G_bc as _orig
        return _orig(geometry, r_nodes, radii, sig_t,
                     n_surf_quad=n_surf_quad, dps=dps)

    r_nodes = np.asarray(r_nodes, dtype=float)
    radii = np.asarray(radii, dtype=float)
    sig_t = np.asarray(sig_t, dtype=float)
    R = float(radii[-1])

    psi_pts, psi_wts = gl_float(n_surf_quad, 0.0, np.pi, dps)

    N = len(r_nodes)
    G_bc = np.zeros(N)
    for i in range(N):
        r_i = float(r_nodes[i])
        total = 0.0
        for k in range(n_surf_quad):
            psi = psi_pts[k]
            cp = float(np.cos(psi))
            sp = float(np.sin(psi))
            # 2-D backward chord from r_i in observer direction ψ
            d_2d = -r_i * cp + np.sqrt(R * R - r_i * r_i * sp * sp)
            if d_2d <= 0:
                continue
            # Multi-region τ along the chord: use the optical_depth_along_ray
            # helper with the ψ direction
            if len(radii) == 1:
                tau_2d = float(sig_t[0]) * float(d_2d)
            else:
                # cos_omega for the BACKWARD direction is -cp (ray direction)
                tau_2d = float(geometry.optical_depth_along_ray(
                    r_i, -cp, float(d_2d), radii, sig_t,
                ))
            total += psi_wts[k] * float(ki_n_mp(2, tau_2d, dps))
        G_bc[i] = float(4.0 / np.pi * total)

    return G_bc


_HEBERT_ENABLED = False  # Module-level toggle


def patched_build_white_bc_correction_rank_n(
    geometry, r_nodes, r_wts, radii, sig_t,
    n_angular=32, n_surf_quad=32, dps=25, n_bc_modes=1,
):
    """Same as the original but uses the corrected cylinder G_bc.

    If ``_HEBERT_ENABLED`` is True at module level, also applies the
    (1 - P_ss^cyl)⁻¹ geometric series factor.
    """
    from orpheus.derivations.peierls_geometry import (
        compute_P_esc, build_closure_operator, BoundaryClosureOperator,
        reflection_marshak,
    )
    if geometry.kind != "cylinder-1d" or n_bc_modes != 1:
        # Use original for non-cylinder or rank > 1
        from orpheus.derivations.peierls_geometry import (
            build_white_bc_correction_rank_n as _orig,
        )
        return _orig(geometry, r_nodes, r_wts, radii, sig_t,
                     n_angular=n_angular, n_surf_quad=n_surf_quad,
                     dps=dps, n_bc_modes=n_bc_modes)

    # Build the rank-1 closure manually with corrected G_bc for cylinder
    R = float(radii[-1])
    sig_t_n = np.array([
        sig_t[geometry.which_annulus(float(r_nodes[i]), radii)]
        for i in range(len(r_nodes))
    ])
    rv = np.array([geometry.radial_volume_weight(float(rj)) for rj in r_nodes])
    divisor = geometry.rank1_surface_divisor(R)

    P_esc_0 = compute_P_esc(geometry, r_nodes, radii, sig_t,
                             n_angular=n_angular, dps=dps)
    G_bc_0 = _compute_G_bc_corrected_cylinder(
        geometry, r_nodes, radii, sig_t,
        n_surf_quad=n_surf_quad, dps=dps,
    )

    P = np.zeros((1, len(r_nodes)))
    G = np.zeros((len(r_nodes), 1))
    P[0, :] = rv * r_wts * P_esc_0
    G[:, 0] = sig_t_n * G_bc_0 / divisor

    op = BoundaryClosureOperator(P=P, G=G, R=np.eye(1))
    K_bc = op.as_matrix()

    if _HEBERT_ENABLED:
        P_ss = compute_P_ss_cylinder(radii, sig_t, n_quad=64, dps=dps)
        K_bc = K_bc / (1.0 - P_ss)

    return K_bc


def row_sum_probe_with_corrected_g_bc(geometry, sig_t_g, radii, with_hebert):
    """Compute K · 1 / σ_t for the corrected G_bc^cyl form."""
    r_nodes, r_wts, panels = composite_gl_r(
        radii, _QUAD_BASE["n_panels_per_region"], _QUAD_BASE["p_order"],
        dps=_QUAD_BASE["dps"],
    )
    K_vol = build_volume_kernel(
        geometry, r_nodes, panels, radii, sig_t_g,
        n_angular=_QUAD_BASE["n_angular"], n_rho=_QUAD_BASE["n_rho"],
        dps=_QUAD_BASE["dps"],
    )
    K_bc = patched_build_white_bc_correction_rank_n(
        geometry, r_nodes, r_wts, radii, sig_t_g,
        n_angular=_QUAD_BASE["n_angular"], n_surf_quad=_QUAD_BASE["n_surf_quad"],
        dps=_QUAD_BASE["dps"], n_bc_modes=1,
    )
    if with_hebert:
        P_ss = compute_P_ss_cylinder(radii, sig_t_g, n_quad=64)
        K_bc = K_bc / (1.0 - P_ss)
    K_total = K_vol + K_bc
    one = np.ones(len(r_nodes))
    rhs = K_total @ one
    sig_t_n = np.array([
        sig_t_g[geometry.which_annulus(float(r_nodes[i]), radii)]
        for i in range(len(r_nodes))
    ])
    rel_dev = (rhs - sig_t_n) / sig_t_n
    return rhs, sig_t_n, rel_dev


def main():
    print("=" * 78)
    print("CORRECTED 3-D G_bc^cyl — row-sum probe + k_eff scan")
    print("=" * 78)

    print("\n--- Row-sum K · 1 / σ_t with patched G_bc^cyl (Hébert closure) ---")
    radii_1r = np.array([1.0])
    sig_t_1g = np.array([1.0])
    rhs, sig_t_n, rel_dev = row_sum_probe_with_corrected_g_bc(
        CYLINDER_1D, sig_t_1g, radii_1r, with_hebert=True,
    )
    print(f"  cyl 1G/1R Hébert: K·1/σ_t mean = {(rhs/sig_t_n).mean():.6f}  "
          f"min = {(rhs/sig_t_n).min():.6f}  max = {(rhs/sig_t_n).max():.6f}")
    print(f"  rel_dev: min = {rel_dev.min():+.4f}  max = {rel_dev.max():+.4f}")

    radii_2r = np.array([0.5, 1.0])
    sig_t_2r = np.array([1.0, 1.0])
    rhs, sig_t_n, rel_dev = row_sum_probe_with_corrected_g_bc(
        CYLINDER_1D, sig_t_2r, radii_2r, with_hebert=True,
    )
    print(f"  cyl 1G/2R-homog Hébert: K·1/σ_t mean = {(rhs/sig_t_n).mean():.6f}")

    # Also without Hébert (just bare Mark + corrected G_bc)
    print("\n--- Row-sum without Hébert (bare Mark + corrected G_bc) ---")
    rhs, sig_t_n, rel_dev = row_sum_probe_with_corrected_g_bc(
        CYLINDER_1D, sig_t_1g, radii_1r, with_hebert=False,
    )
    print(f"  cyl 1G/1R bare Mark: K·1/σ_t mean = {(rhs/sig_t_n).mean():.6f}")

    # k_eff scan with the patched closure (toggle Hébert factor)
    global _HEBERT_ENABLED
    print("\n--- k_eff scan: cylinder Class B with corrected G_bc ---")
    print(f"  {'Config':<12} {'cp k_inf':>10} {'k_bare-Mark':>14} {'err':>9} "
          f"{'k_+Hebert':>14} {'err':>9}")
    print("  " + "-" * 80)
    orig_build = _pg.build_white_bc_correction_rank_n
    _pg.build_white_bc_correction_rank_n = patched_build_white_bc_correction_rank_n

    try:
        for ng_key, n_regions in [("1g", 1), ("1g", 2), ("2g", 1), ("2g", 2)]:
            layout = LAYOUTS[n_regions]
            xs_list = [get_xs(r, ng_key) for r in layout]
            sig_t = np.vstack([xs["sig_t"] for xs in xs_list])
            sig_s = np.stack([xs["sig_s"] for xs in xs_list], axis=0)
            nu_sig_f = np.vstack([xs["nu"] * xs["sig_f"] for xs in xs_list])
            chi = np.vstack([xs["chi"] for xs in xs_list])
            radii = np.array(cp_cylinder._RADII[n_regions])
            kinf = cp_cylinder._build_case(ng_key, n_regions).k_inf

            _HEBERT_ENABLED = False
            sol_bare = solve_peierls_mg(
                CYLINDER_1D, radii=radii, sig_t=sig_t, sig_s=sig_s,
                nu_sig_f=nu_sig_f, chi=chi,
                boundary="white_rank1_mark", n_bc_modes=1,
                **_QUAD_BASE,
            )
            err_bare = (sol_bare.k_eff - kinf) / kinf * 100

            _HEBERT_ENABLED = True
            sol_heb = solve_peierls_mg(
                CYLINDER_1D, radii=radii, sig_t=sig_t, sig_s=sig_s,
                nu_sig_f=nu_sig_f, chi=chi,
                boundary="white_rank1_mark", n_bc_modes=1,
                **_QUAD_BASE,
            )
            err_heb = (sol_heb.k_eff - kinf) / kinf * 100

            print(f"  {ng_key} {n_regions}r          {kinf:>10.6f} "
                  f"{sol_bare.k_eff:>14.10f} {err_bare:>+8.3f}% "
                  f"{sol_heb.k_eff:>14.10f} {err_heb:>+8.3f}%")
    finally:
        _pg.build_white_bc_correction_rank_n = orig_build
        _HEBERT_ENABLED = False


if __name__ == "__main__":
    main()
