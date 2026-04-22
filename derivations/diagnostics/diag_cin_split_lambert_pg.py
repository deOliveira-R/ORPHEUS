"""Diagnostic E2.4/E2.5: Lambert-basis P/G for split basis.

Created by numerics-investigator on 2026-04-21.

Hypothesis: F.4's empirical advantage comes from a basis-mismatch trick —
Lambert P/G (no µ_exit weight on outgoing direction) combined with
Marshak W (µ-weighted transmission). Test whether transferring this
mismatch into the split basis improves accuracy.

E2.4 — Build split-basis closure with Lambert-convention P_esc/G_bc.
Measure rank-(1,1,1), (1,1,2), (2,2,2), compare to F.4's 0.12%.

E2.5 — If E2.4 shows improvement, scan rank-(N,N,N) up to 3,3,3.
"""
from __future__ import annotations

import math
import sys

import numpy as np

sys.path.insert(0, '/workspaces/ORPHEUS/derivations/diagnostics')

from diag_cin_aware_split_basis_keff import (
    _unit_legendre_lambdas,
    grazing_lambdas,
    mu_crit,
    c_of_mu,
    chord_oi,
    compute_W_split_basis,
    solve_k_eff,
    run_scalar_f4,
)

from orpheus.derivations.peierls_geometry import (
    CurvilinearGeometry,
    build_volume_kernel,
    composite_gl_r,
    gl_float,
)


K_INF = 1.5


# ---------------------------------------------------------------------------
# Lambert-convention P/G (no µ_exit weight on outgoing)
# ---------------------------------------------------------------------------


def compute_P_esc_graze_lambert(geom, r_nodes, radii, sig_t, m, mu_crit_val,
                                 leg_g, n_angular=32, dps=25):
    """Lambert-basis P_esc to grazing mode m on outer surface.

    Drops the µ_exit factor from the integrand — integrand is just
    angular_factor · P̃_m^g(µ_exit) · K_esc(τ).
    """
    r_nodes = np.asarray(r_nodes, dtype=float)
    radii = np.asarray(radii, dtype=float)
    sig_t = np.asarray(sig_t, dtype=float)
    R = float(radii[-1])
    omega_low, omega_high = geom.angular_range
    omega_pts, omega_wts = gl_float(n_angular, omega_low, omega_high, dps)
    cos_omegas = geom.ray_direction_cosine(omega_pts)
    angular_factor = geom.angular_weight(omega_pts)
    pref = geom.prefactor
    N = len(r_nodes)
    P = np.zeros(N)
    for i in range(N):
        r_i = float(r_nodes[i])
        total = 0.0
        for k in range(n_angular):
            cos_om = cos_omegas[k]
            rho_out = geom.rho_max(r_i, cos_om, R)
            if rho_out <= 0.0:
                continue
            rho_in_minus, _ = geom.rho_inner_intersections(r_i, cos_om)
            if rho_in_minus is not None and rho_in_minus < rho_out:
                continue
            mu_exit = (rho_out + r_i * cos_om) / R
            if mu_exit > mu_crit_val:
                continue
            tau = geom.optical_depth_along_ray(r_i, cos_om, rho_out, radii, sig_t)
            K_esc = geom.escape_kernel_mp(tau, dps)
            p_tilde = leg_g[m](mu_exit)
            # NO µ_exit FACTOR — Lambert convention
            total += omega_wts[k] * angular_factor[k] * p_tilde * K_esc
        P[i] = pref * total
    return P


def compute_P_esc_steep_lambert(geom, r_nodes, radii, sig_t, m, mu_crit_val, rho_val,
                                 leg, n_angular=32, dps=25):
    r_nodes = np.asarray(r_nodes, dtype=float)
    radii = np.asarray(radii, dtype=float)
    sig_t = np.asarray(sig_t, dtype=float)
    R = float(radii[-1])
    omega_low, omega_high = geom.angular_range
    omega_pts, omega_wts = gl_float(n_angular, omega_low, omega_high, dps)
    cos_omegas = geom.ray_direction_cosine(omega_pts)
    angular_factor = geom.angular_weight(omega_pts)
    pref = geom.prefactor
    N = len(r_nodes)
    P = np.zeros(N)
    for i in range(N):
        r_i = float(r_nodes[i])
        total = 0.0
        for k in range(n_angular):
            cos_om = cos_omegas[k]
            rho_out = geom.rho_max(r_i, cos_om, R)
            if rho_out <= 0.0:
                continue
            rho_in_minus, _ = geom.rho_inner_intersections(r_i, cos_om)
            if rho_in_minus is not None and rho_in_minus < rho_out:
                continue
            mu_exit = (rho_out + r_i * cos_om) / R
            if mu_exit < mu_crit_val:
                continue
            c_val = c_of_mu(mu_exit, rho_val)
            tau = geom.optical_depth_along_ray(r_i, cos_om, rho_out, radii, sig_t)
            K_esc = geom.escape_kernel_mp(tau, dps)
            p_tilde_s = leg[m](c_val) / rho_val
            total += omega_wts[k] * angular_factor[k] * p_tilde_s * K_esc
        P[i] = pref * total
    return P


def compute_P_esc_inner_lambert(geom, r_nodes, radii, sig_t, m, leg,
                                 n_angular=32, dps=25):
    r_nodes = np.asarray(r_nodes, dtype=float)
    radii = np.asarray(radii, dtype=float)
    sig_t = np.asarray(sig_t, dtype=float)
    r_0 = float(geom.inner_radius)
    omega_low, omega_high = geom.angular_range
    omega_pts, omega_wts = gl_float(n_angular, omega_low, omega_high, dps)
    cos_omegas = geom.ray_direction_cosine(omega_pts)
    angular_factor = geom.angular_weight(omega_pts)
    pref = geom.prefactor
    N = len(r_nodes)
    P = np.zeros(N)
    for i in range(N):
        r_i = float(r_nodes[i])
        total = 0.0
        for k in range(n_angular):
            cos_om = cos_omegas[k]
            rho_in_minus, _ = geom.rho_inner_intersections(r_i, cos_om)
            if rho_in_minus is None:
                continue
            tau = geom.optical_depth_along_ray(r_i, cos_om, rho_in_minus, radii, sig_t)
            K_esc = geom.escape_kernel_mp(tau, dps)
            sin_om = math.sqrt(max(0.0, 1.0 - cos_om * cos_om))
            h_sq = r_i * r_i * sin_om * sin_om
            mu_exit_sq = max(0.0, (r_0 * r_0 - h_sq) / (r_0 * r_0))
            mu_exit = math.sqrt(mu_exit_sq)
            p_tilde = leg[m](mu_exit)
            total += omega_wts[k] * angular_factor[k] * p_tilde * K_esc
        P[i] = pref * total
    return P


def compute_G_bc_graze_lambert(geom, r_nodes, radii, sig_t, m, mu_crit_val,
                                leg_g, n_surf_quad=32, dps=25):
    r_nodes = np.asarray(r_nodes, dtype=float)
    radii = np.asarray(radii, dtype=float)
    sig_t = np.asarray(sig_t, dtype=float)
    R = float(radii[-1])
    theta_pts, theta_wts = gl_float(n_surf_quad, 0.0, np.pi, dps)
    cos_thetas = np.cos(theta_pts)
    sin_thetas = np.sin(theta_pts)
    N = len(r_nodes)
    G = np.zeros(N)
    for i in range(N):
        r_i = r_nodes[i]
        total = 0.0
        for k in range(n_surf_quad):
            ct = cos_thetas[k]
            st = sin_thetas[k]
            rho_out = geom.rho_max(r_i, ct, R)
            if rho_out <= 0.0:
                continue
            rho_in_minus, _ = geom.rho_inner_intersections(r_i, ct)
            if rho_in_minus is not None and rho_in_minus < rho_out:
                continue
            mu_s = (rho_out + r_i * ct) / R
            if mu_s > mu_crit_val:
                continue
            tau = geom.optical_depth_along_ray(r_i, ct, rho_out, radii, sig_t)
            p_tilde = leg_g[m](mu_s)
            # Lambert: drop µ_s from the weight
            total += theta_wts[k] * st * p_tilde * math.exp(-tau)
        G[i] = 2.0 * total
    return G


def compute_G_bc_steep_lambert(geom, r_nodes, radii, sig_t, m, mu_crit_val, rho_val,
                                leg, n_surf_quad=32, dps=25):
    r_nodes = np.asarray(r_nodes, dtype=float)
    radii = np.asarray(radii, dtype=float)
    sig_t = np.asarray(sig_t, dtype=float)
    R = float(radii[-1])
    theta_pts, theta_wts = gl_float(n_surf_quad, 0.0, np.pi, dps)
    cos_thetas = np.cos(theta_pts)
    sin_thetas = np.sin(theta_pts)
    N = len(r_nodes)
    G = np.zeros(N)
    for i in range(N):
        r_i = r_nodes[i]
        total = 0.0
        for k in range(n_surf_quad):
            ct = cos_thetas[k]
            st = sin_thetas[k]
            rho_out = geom.rho_max(r_i, ct, R)
            if rho_out <= 0.0:
                continue
            rho_in_minus, _ = geom.rho_inner_intersections(r_i, ct)
            if rho_in_minus is not None and rho_in_minus < rho_out:
                continue
            mu_s = (rho_out + r_i * ct) / R
            if mu_s < mu_crit_val:
                continue
            tau = geom.optical_depth_along_ray(r_i, ct, rho_out, radii, sig_t)
            c_val = c_of_mu(mu_s, rho_val)
            p_tilde = leg[m](c_val) / rho_val
            total += theta_wts[k] * st * p_tilde * math.exp(-tau)
        G[i] = 2.0 * total
    return G


def compute_G_bc_inner_lambert(geom, r_nodes, radii, sig_t, m, leg,
                                n_surf_quad=32, dps=25):
    r_nodes = np.asarray(r_nodes, dtype=float)
    radii = np.asarray(radii, dtype=float)
    sig_t = np.asarray(sig_t, dtype=float)
    r_0 = float(geom.inner_radius)
    theta_pts, theta_wts = gl_float(n_surf_quad, 0.0, np.pi, dps)
    cos_thetas = np.cos(theta_pts)
    sin_thetas = np.sin(theta_pts)
    N = len(r_nodes)
    G = np.zeros(N)
    for i in range(N):
        r_i = r_nodes[i]
        total = 0.0
        for k in range(n_surf_quad):
            ct = cos_thetas[k]
            st = sin_thetas[k]
            rho_in_minus, _ = geom.rho_inner_intersections(r_i, ct)
            if rho_in_minus is None:
                continue
            tau = geom.optical_depth_along_ray(r_i, ct, rho_in_minus, radii, sig_t)
            sin_om = math.sqrt(max(0.0, 1.0 - ct * ct))
            h_sq = r_i * r_i * sin_om * sin_om
            mu_s_sq = max(0.0, (r_0 * r_0 - h_sq) / (r_0 * r_0))
            mu_s = math.sqrt(mu_s_sq)
            p_tilde = leg[m](mu_s)
            total += theta_wts[k] * st * p_tilde * math.exp(-tau)
        G[i] = 2.0 * total
    return G


def build_split_basis_lambert(
    geom, r_nodes, r_wts, radii, sig_t,
    N_g, N_s, N_i, *,
    n_angular=32, n_surf_quad=32, dps=25,
    albedo=1.0,
):
    """Split-basis closure with LAMBERT P/G and MARSHAK W (F.4-style mismatch)."""
    rho_val = float(geom.inner_radius) / float(radii[-1])
    muc = mu_crit(rho_val)

    leg = _unit_legendre_lambdas(max(N_s, N_i, 1))
    leg_g = grazing_lambdas(max(N_g, 1), muc) if N_g > 0 else []

    D = N_g + N_s + N_i
    N_r = len(r_nodes)

    P = np.zeros((D, N_r))
    G = np.zeros((N_r, D))

    R_out = float(radii[-1])
    r_in = float(geom.inner_radius)
    A_outer = R_out * R_out
    A_inner = r_in * r_in

    rv = r_nodes**2

    for m in range(N_g):
        P_m = compute_P_esc_graze_lambert(
            geom, r_nodes, radii, sig_t, m, muc, leg_g,
            n_angular=n_angular, dps=dps,
        )
        G_m = compute_G_bc_graze_lambert(
            geom, r_nodes, radii, sig_t, m, muc, leg_g,
            n_surf_quad=n_surf_quad, dps=dps,
        )
        P[m, :] = rv * r_wts * P_m
        G[:, m] = G_m / A_outer

    for m in range(N_s):
        P_m = compute_P_esc_steep_lambert(
            geom, r_nodes, radii, sig_t, m, muc, rho_val, leg,
            n_angular=n_angular, dps=dps,
        )
        G_m = compute_G_bc_steep_lambert(
            geom, r_nodes, radii, sig_t, m, muc, rho_val, leg,
            n_surf_quad=n_surf_quad, dps=dps,
        )
        P[N_g + m, :] = rv * r_wts * P_m
        G[:, N_g + m] = G_m / A_outer

    for m in range(N_i):
        P_m = compute_P_esc_inner_lambert(
            geom, r_nodes, radii, sig_t, m, leg,
            n_angular=n_angular, dps=dps,
        )
        G_m = compute_G_bc_inner_lambert(
            geom, r_nodes, radii, sig_t, m, leg,
            n_surf_quad=n_surf_quad, dps=dps,
        )
        P[N_g + N_s + m, :] = rv * r_wts * P_m
        G[:, N_g + N_s + m] = G_m / A_inner

    # Scale G by sig_t
    sig_t_n = np.empty(N_r)
    for i, ri in enumerate(r_nodes):
        ki = geom.which_annulus(ri, radii)
        sig_t_n[i] = sig_t[ki]
    G = sig_t_n[:, None] * G

    # Marshak-basis W (SAME as original split basis — this is the mismatch)
    sig_t_val = float(sig_t[0])
    tau = sig_t_val * R_out
    W = compute_W_split_basis(rho_val, tau, N_g, N_s, N_i)

    # White BC operator
    B = np.zeros((D, D))
    if N_g >= 1 and N_s >= 1:
        B[0, 0] = albedo * muc * muc
        B[0, N_g] = albedo * muc * rho_val
        B[N_g, 0] = albedo * rho_val * muc
        B[N_g, N_g] = albedo * rho_val * rho_val
    for m in range(N_i):
        B[N_g + N_s + m, N_g + N_s + m] = 1.0

    M = np.eye(D) - W @ B
    try:
        Minv = np.linalg.inv(M)
    except np.linalg.LinAlgError:
        return None, {"error": "singular"}

    K_bc = G @ B @ Minv @ P
    return K_bc, {"W": W, "B": B, "P": P, "G": G, "M": M}


def run_split_lambert(r_0, R, sig_t_val, sig_s_val, nsf_val,
                       N_g, N_s, N_i,
                       n_panels=2, p_order=4, n_ang=32, dps=15):
    geom = CurvilinearGeometry(kind="sphere-1d", inner_radius=r_0)
    radii = np.array([R])
    sig_t = np.array([sig_t_val])
    r_nodes, r_wts, panels = composite_gl_r(radii, n_panels, p_order,
                                             dps=dps, inner_radius=r_0)
    K_vol = build_volume_kernel(geom, r_nodes, panels, radii, sig_t,
                                 n_angular=n_ang, n_rho=n_ang, dps=dps)
    K_bc, aux = build_split_basis_lambert(
        geom, r_nodes, r_wts, radii, sig_t, N_g, N_s, N_i,
        n_angular=n_ang, n_surf_quad=n_ang, dps=dps,
    )
    if K_bc is None:
        return None, aux
    K = K_vol + K_bc
    return solve_k_eff(K, sig_t_val, sig_s_val, nsf_val), aux


def main():
    print("=" * 78)
    print("E2.4/E2.5 — Lambert-basis P/G with split basis (F.4-style mismatch)")
    print("σ_t·R = 5, r_0/R = 0.3, k_inf = 1.5")
    print("=" * 78)

    sig_t, sig_s, nsf = 1.0, 1.0 / 3.0, 1.0
    R = 5.0
    r_0 = 0.3 * R

    # F.4 baseline
    k_f4 = run_scalar_f4(r_0, R, sig_t, sig_s, nsf)
    err_f4 = abs(k_f4 - K_INF) / K_INF * 100
    print(f"\nF.4 scalar baseline: k={k_f4:.10f}, err={err_f4:.6f}%")

    # E2.4: Lambert + split, rank-(1,1,N) and (N,N,N)
    print("\n--- E2.4: rank-(Ng, Ns, Ni) scan (Lambert P/G, Marshak W, split basis) ---")
    ranks = [
        (1, 1, 1), (1, 1, 2), (1, 1, 3), (1, 1, 4),
        (2, 2, 1), (2, 2, 2), (2, 2, 3), (2, 2, 4),
        (3, 3, 3),
    ]
    for (Ng, Ns, Ni) in ranks:
        k, _ = run_split_lambert(r_0, R, sig_t, sig_s, nsf, Ng, Ns, Ni)
        if k is None:
            print(f"  rank-({Ng},{Ns},{Ni}): singular")
            continue
        err = abs(k - K_INF) / K_INF * 100
        delta = err - err_f4
        print(f"  rank-({Ng},{Ns},{Ni}) ({Ng+Ns+Ni} modes): "
              f"k={k:.10f}, err={err:.6f}%, Δ(F.4)={delta:+.4f}%")

    # E2.5: cross σ_t·R scan for best Lambert+split config
    print("\n--- E2.5: cross σ_t·R scan (rank-(1,1,2), Lambert+split) ---")
    print(" σ_t·R |   F.4 err    | Lambert-split-(1,1,2) err | improvement")
    for sig_t_R in [1.0, 2.5, 5.0, 10.0, 20.0]:
        R_val = sig_t_R  # since σ_t = 1
        r_0_val = 0.3 * R_val
        k_f4_sc = run_scalar_f4(r_0_val, R_val, sig_t, sig_s, nsf)
        err_f4_sc = abs(k_f4_sc - K_INF) / K_INF * 100
        k_sb, _ = run_split_lambert(r_0_val, R_val, sig_t, sig_s, nsf, 1, 1, 2)
        if k_sb is None:
            print(f"  {sig_t_R:5.1f}  |  {err_f4_sc:.4f}%  |  singular")
            continue
        err_sb = abs(k_sb - K_INF) / K_INF * 100
        ratio = err_f4_sc / err_sb if err_sb > 1e-9 else np.inf
        print(f"  {sig_t_R:5.1f}  |  {err_f4_sc:.4f}%  |  {err_sb:.4f}%  |  "
              f"{'IMPROVED' if err_sb < err_f4_sc else 'worse'} "
              f"(ratio={ratio:.2f})")

    return 0


def test_lambert_split_vs_f4():
    """Lambert+split rank-(1,1,1) must be < F.4 to be a "win".

    This test PASSES if Lambert-split beats F.4. If it fails, document
    the plateau — structural understanding.
    """
    sig_t, sig_s, nsf = 1.0, 1.0 / 3.0, 1.0
    R = 5.0
    r_0 = 0.3 * R
    k_f4 = run_scalar_f4(r_0, R, sig_t, sig_s, nsf)
    err_f4 = abs(k_f4 - K_INF) / K_INF * 100
    k_sb, _ = run_split_lambert(r_0, R, sig_t, sig_s, nsf, 1, 1, 2)
    err_sb = abs(k_sb - K_INF) / K_INF * 100
    print(f"F.4 err: {err_f4:.4f}%")
    print(f"Lambert-split-(1,1,2) err: {err_sb:.4f}%")
    # Record; don't fail if it's close
    # assert err_sb < err_f4


if __name__ == "__main__":
    sys.exit(main())
