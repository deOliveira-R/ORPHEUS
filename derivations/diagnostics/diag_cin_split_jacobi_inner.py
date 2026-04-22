"""Diagnostic E2.6: non-Legendre inner basis (Jacobi) for split basis.

Created by numerics-investigator on 2026-04-21.

E2.3 showed the residual lives mostly in INNER mode-1 when we use
rank-(1,1,1), and rank-(1,1,2) captures both significant modes — BUT
the plateau is at 1% even with many inner modes.

The plateau CAN'T come from inner mode resolution (mode-2+ energy is <0.1%).
So the residual must come from the outer representation. But the c_in-aware
split basis on outer is ALREADY a specialized basis that maps the steep
cone to inner-c coordinates. The plateau must be in:
  (A) the grazing sub-basis (shape on [0, µ_crit]) being wrong, OR
  (B) the µ-weight convention itself (formally-consistent but suboptimal).

Goal: test whether a different inner basis (Jacobi P^{(1,0)}) changes
anything. If the plateau is truly "we have the right inner modes just
weighted oddly", Jacobi should give the same number. If the plateau is
"we need different information at inner", Jacobi might help.
"""
from __future__ import annotations

import math
import sys

import numpy as np
import sympy as sp
from scipy import integrate

sys.path.insert(0, '/workspaces/ORPHEUS/derivations/diagnostics')

from diag_cin_aware_split_basis_keff import (
    grazing_lambdas,
    mu_crit,
    c_of_mu,
    chord_oi,
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


def jacobi_c_weighted_symbolic(n_max, alpha, beta, a, b, var, weight_fn):
    """Gram-Schmidt of {1, x, x²,...} with custom weight.

    Returns orthonormal list in given variable.
    """
    basis = []
    for k in range(n_max + 1):
        fk = var**k
        for pj in basis:
            cjk = sp.integrate(fk * pj * weight_fn(var), (var, a, b))
            fk = fk - cjk * pj
        norm_sq = sp.integrate(fk * fk * weight_fn(var), (var, a, b))
        fk = sp.simplify(fk / sp.sqrt(norm_sq))
        basis.append(fk)
    return basis


def jacobi_inner_lambdas(n_max, alpha=0, beta=1):
    """Jacobi-weighted orthonormal basis on c ∈ [0, 1].

    Weight function = c^(α+1) · (1-c)^β, matching P^{(α,β)} on [-1,1] ∩ Gauss
    mapping to [0,1] but with µ-weight convention (extra c factor).
    If alpha=0, beta=0: same as half-range Legendre in c-weight.
    If alpha=1: weight = c² — favors large c (steep angles).
    If beta=1: weight = c*(1-c) — emphasizes mid-range.
    """
    c = sp.Symbol("c", positive=True, real=True)
    if alpha == 0 and beta == 0:
        weight = lambda x: x  # c-weight (matches Legendre half-range)
    elif alpha == 1 and beta == 0:
        weight = lambda x: x**2  # c²-weight
    elif alpha == 0 and beta == 1:
        weight = lambda x: x * (1 - x)
    else:
        weight = lambda x: x ** (alpha + 1) * (1 - x) ** beta
    exprs = jacobi_c_weighted_symbolic(n_max, alpha, beta, 0, 1, c, weight)
    return [sp.lambdify(c, e, modules=["numpy"]) for e in exprs]


def compute_W_jacobi_basis(rho, tau, N_g, N_s, N_i, *,
                            leg_inner_funcs, inner_weight_fn,
                            epsabs=1e-13, epsrel=1e-11):
    """Build W with Jacobi inner basis.

    The outer-steep-to-inner transmission integral changes:
      W_si[m, n] = (1/ρ) ∫_0^1 P̃_m^I(c) P̃_n^s(c_I(µ)) exp(-τ χ) · weight_fn(c) dc
    where the weight is determined by the inner inner-product convention.

    For plain Legendre (half-range), weight = c (same as before).
    For Jacobi (α=1, β=0), weight = c².

    However, the outer steep sub-basis still uses µ-weight (from its
    normalization). The coupling is an overlap between different
    weight conventions — this MIGHT be formally inconsistent, but
    empirical test is cheap.
    """
    from diag_cin_aware_split_basis_keff import _unit_legendre_lambdas

    muc = mu_crit(rho)
    D = N_g + N_s + N_i
    W = np.zeros((D, D))

    leg_s = _unit_legendre_lambdas(max(N_s, 1))  # outer steep: plain Legendre
    leg_g = grazing_lambdas(max(N_g, 1), muc) if N_g > 0 else []

    g0, g1 = 0, N_g
    s0, s1 = N_g, N_g + N_s
    i0, i1 = N_g + N_s, D

    # W_gg unchanged (grazing self)
    for m in range(N_g):
        for n in range(N_g):
            def integrand(mu, m=m, n=n):
                return leg_g[m](mu) * leg_g[n](mu) * math.exp(-2.0 * tau * mu) * mu
            val, _ = integrate.quad(integrand, 0.0, muc, epsabs=epsabs, epsrel=epsrel)
            W[g0 + m, g0 + n] = val

    # W_si: inner-receiving-mode m (Jacobi), steep-outgoing-mode n (Legendre)
    # Integrand: P̃_m^I(c) · [P̃_n(c) / ρ] · exp(-τ χ(c)) · <weight>
    # If we use the inner c-weighted inner product (same as Legendre plain),
    # the inner weight is c. Stick with that.
    for m in range(N_i):
        for n in range(N_s):
            def integrand(c, m=m, n=n):
                return leg_inner_funcs[m](c) * leg_s[n](c) * \
                       math.exp(-tau * chord_oi(c, rho)) * c
            val, _ = integrate.quad(integrand, 0.0, 1.0, epsabs=epsabs, epsrel=epsrel)
            W[i0 + m, s0 + n] = val / rho

    for m in range(N_s):
        for n in range(N_i):
            def integrand(c, m=m, n=n):
                return leg_s[m](c) * leg_inner_funcs[n](c) * \
                       math.exp(-tau * chord_oi(c, rho)) * c
            val, _ = integrate.quad(integrand, 0.0, 1.0, epsabs=epsabs, epsrel=epsrel)
            W[s0 + m, i0 + n] = val / rho

    return W


def compute_P_esc_inner_jacobi(geom, r_nodes, radii, sig_t, m, jacobi_lam,
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
            p_tilde = jacobi_lam[m](mu_exit)
            total += omega_wts[k] * angular_factor[k] * mu_exit * p_tilde * K_esc
        P[i] = pref * total
    return P


def compute_G_bc_inner_jacobi(geom, r_nodes, radii, sig_t, m, jacobi_lam,
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
            p_tilde = jacobi_lam[m](mu_s)
            total += theta_wts[k] * st * mu_s * p_tilde * math.exp(-tau)
        G[i] = 2.0 * total
    return G


def run_split_jacobi(r_0, R, sig_t_val, sig_s_val, nsf_val, N_i,
                      alpha=1, beta=0,
                      n_panels=2, p_order=4, n_ang=32, dps=15):
    """Split basis with N_g=N_s=1 and Jacobi inner basis."""
    from diag_cin_aware_split_basis_keff import (
        _unit_legendre_lambdas,
        compute_P_esc_graze_mode,
        compute_P_esc_steep_mode,
        compute_G_bc_graze_mode,
        compute_G_bc_steep_mode,
    )

    geom = CurvilinearGeometry(kind="sphere-1d", inner_radius=r_0)
    radii = np.array([R])
    sig_t = np.array([sig_t_val])
    r_nodes, r_wts, panels = composite_gl_r(radii, n_panels, p_order,
                                             dps=dps, inner_radius=r_0)
    K_vol = build_volume_kernel(geom, r_nodes, panels, radii, sig_t,
                                 n_angular=n_ang, n_rho=n_ang, dps=dps)

    N_g = N_s = 1
    rho_val = r_0 / R
    muc = mu_crit(rho_val)

    leg_s = _unit_legendre_lambdas(max(N_s, 1))
    leg_g = grazing_lambdas(max(N_g, 1), muc)
    jacobi_lam = jacobi_inner_lambdas(N_i - 1, alpha=alpha, beta=beta)

    D = N_g + N_s + N_i
    N_r = len(r_nodes)
    P = np.zeros((D, N_r))
    G = np.zeros((N_r, D))
    R_out = float(radii[-1])
    A_outer = R_out * R_out
    A_inner = r_0 ** 2
    rv = r_nodes ** 2

    # Grazing
    P_m = compute_P_esc_graze_mode(geom, r_nodes, radii, sig_t, 0, muc, leg_g,
                                     n_angular=n_ang, dps=dps)
    G_m = compute_G_bc_graze_mode(geom, r_nodes, radii, sig_t, 0, muc, leg_g,
                                    n_surf_quad=n_ang, dps=dps)
    P[0, :] = rv * r_wts * P_m
    G[:, 0] = G_m / A_outer

    # Steep
    P_m = compute_P_esc_steep_mode(geom, r_nodes, radii, sig_t, 0, muc, rho_val, leg_s,
                                     n_angular=n_ang, dps=dps)
    G_m = compute_G_bc_steep_mode(geom, r_nodes, radii, sig_t, 0, muc, rho_val, leg_s,
                                    n_surf_quad=n_ang, dps=dps)
    P[1, :] = rv * r_wts * P_m
    G[:, 1] = G_m / A_outer

    # Inner (Jacobi)
    for m in range(N_i):
        P_m = compute_P_esc_inner_jacobi(geom, r_nodes, radii, sig_t, m, jacobi_lam,
                                           n_angular=n_ang, dps=dps)
        G_m = compute_G_bc_inner_jacobi(geom, r_nodes, radii, sig_t, m, jacobi_lam,
                                          n_surf_quad=n_ang, dps=dps)
        P[2 + m, :] = rv * r_wts * P_m
        G[:, 2 + m] = G_m / A_inner

    # sig_t scaling
    sig_t_n = np.empty(N_r)
    for i, ri in enumerate(r_nodes):
        ki = geom.which_annulus(ri, radii)
        sig_t_n[i] = sig_t[ki]
    G = sig_t_n[:, None] * G

    # W (uses Jacobi in inner block)
    tau = sig_t_val * R_out
    W = compute_W_jacobi_basis(rho_val, tau, N_g, N_s, N_i,
                                leg_inner_funcs=jacobi_lam, inner_weight_fn=None)

    # B operator
    B = np.zeros((D, D))
    B[0, 0] = muc * muc
    B[0, N_g] = muc * rho_val
    B[N_g, 0] = rho_val * muc
    B[N_g, N_g] = rho_val * rho_val
    for m in range(N_i):
        B[N_g + N_s + m, N_g + N_s + m] = 1.0

    M = np.eye(D) - W @ B
    K_bc = G @ B @ np.linalg.inv(M) @ P
    K = K_vol + K_bc
    return solve_k_eff(K, sig_t_val, sig_s_val, nsf_val)


def main():
    print("=" * 78)
    print("E2.6 — Jacobi inner basis for split-basis closure")
    print("σ_t·R = 5, r_0/R = 0.3, k_inf = 1.5")
    print("=" * 78)

    sig_t, sig_s, nsf = 1.0, 1.0 / 3.0, 1.0
    R = 5.0
    r_0 = 0.3 * R

    k_f4 = run_scalar_f4(r_0, R, sig_t, sig_s, nsf)
    err_f4 = abs(k_f4 - K_INF) / K_INF * 100
    print(f"\nF.4 err: {err_f4:.6f}%")

    # Baseline Legendre-inner (α=0, β=0) = same as plain
    for (alpha, beta, label) in [(0, 0, "Legendre-inner"),
                                   (1, 0, "Jacobi c²-weighted"),
                                   (0, 1, "Jacobi (1-c)-weighted")]:
        print(f"\n{label} (α={alpha}, β={beta}):")
        for N_i in [1, 2, 3, 4]:
            try:
                k = run_split_jacobi(r_0, R, sig_t, sig_s, nsf, N_i,
                                     alpha=alpha, beta=beta)
                err = abs(k - K_INF) / K_INF * 100
                print(f"  rank-(1,1,{N_i}): k={k:.10f}, err={err:.6f}%")
            except Exception as e:
                print(f"  rank-(1,1,{N_i}): FAILED - {type(e).__name__}: {e}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
