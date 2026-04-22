"""Diagnostic E3.2/E3.3: Physics-asymptote adaptive inner basis.

Created by numerics-investigator on 2026-04-22.

Hypothesis: the optimal inner-surface mode-0 function tracks the
asymptotic shape of ψ^+_inner(c):
    f_0(c; β) = exp(-β · s(c; ρ) / 2)
where s(c; ρ) = sqrt(1 - ρ²(1-c²)) - ρc is the normalized chord from
outer to inner, and β is a calibration parameter (physically β=τ).

Higher modes: Gram-Schmidt({f_0 · P̃_n(c)}) with c-weight on [0,1].

If β=τ gives near-exact k_eff at rank-(1,1,1), we have a universal
principled closure.
"""
from __future__ import annotations

import math
import sys

import numpy as np
import sympy as sp
from scipy import integrate

sys.path.insert(0, '/workspaces/ORPHEUS/derivations/diagnostics')

from diag_cin_aware_split_basis_keff import (
    _unit_legendre_lambdas,
    grazing_lambdas,
    mu_crit,
    c_of_mu,
    chord_oi,
    solve_k_eff,
    run_scalar_f4,
    compute_P_esc_graze_mode,
    compute_P_esc_steep_mode,
    compute_G_bc_graze_mode,
    compute_G_bc_steep_mode,
)
from orpheus.derivations.peierls_geometry import (
    CurvilinearGeometry,
    build_volume_kernel,
    composite_gl_r,
    gl_float,
)


K_INF = 1.5


def build_asymptote_basis_numeric(N_i, beta, rho, *, n_gram=200):
    """Build Gram-Schmidt orthonormal basis for inner flux with:

        mode-0: g_0(c) = exp(-beta * s(c; rho) / 2) / ||.||
        mode-n: Gram-Schmidt of (g_0(c) * c^n, ... ) with c-weight

    Works NUMERICALLY with a fine Gauss-Legendre grid for the inner
    product. Returns a list of callables [f_0, f_1, ..., f_{N_i-1}].
    """
    # Build a fine Gauss-Legendre grid on [0, 1] for c-weighted inner products
    x_gl, w_gl = np.polynomial.legendre.leggauss(n_gram)
    c_grid = 0.5 * (x_gl + 1.0)  # [0, 1]
    w_grid = 0.5 * w_gl * c_grid  # µ-weight (c-weight) on [0, 1]

    def s_of_c(c):
        # s(c; rho) = sqrt(1 - rho^2 (1-c^2)) - rho*c
        return np.sqrt(np.maximum(0.0, 1.0 - rho * rho * (1.0 - c * c))) - rho * c

    # seed function (mode-0 shape)
    g0 = np.exp(-beta * s_of_c(c_grid) / 2.0)

    # Build raw candidates: phi_n = g_0(c) * c^n
    phis = np.stack([g0 * c_grid ** n for n in range(N_i)], axis=0)

    # Gram-Schmidt with c-weighted inner product
    basis_vals = []  # list of numerical arrays on c_grid
    for n in range(N_i):
        f = phis[n].copy()
        for bj in basis_vals:
            cj = np.sum(w_grid * f * bj)
            f = f - cj * bj
        norm_sq = np.sum(w_grid * f * f)
        if norm_sq <= 0:
            raise RuntimeError(f"Degenerate mode at n={n}")
        f = f / math.sqrt(norm_sq)
        basis_vals.append(f)

    # Return callables via interpolation
    basis_funcs = []
    for f_vals in basis_vals:
        # Capture by default arg
        def make_f(f_vals=f_vals.copy()):
            from scipy.interpolate import interp1d
            interp = interp1d(c_grid, f_vals, kind='cubic',
                              bounds_error=False, fill_value="extrapolate")
            return lambda c: float(interp(c)) if np.isscalar(c) else interp(c)
        basis_funcs.append(make_f())
    return basis_funcs


def compute_P_esc_inner_asym(geom, r_nodes, radii, sig_t, m, basis,
                              n_angular=32, dps=25):
    """Projection of volume-escape onto inner surface mode m (asymptote basis)."""
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
            p_tilde = basis[m](mu_exit)
            total += omega_wts[k] * angular_factor[k] * mu_exit * p_tilde * K_esc
        P[i] = pref * total
    return P


def compute_G_bc_inner_asym(geom, r_nodes, radii, sig_t, m, basis,
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
            p_tilde = basis[m](mu_s)
            total += theta_wts[k] * st * mu_s * p_tilde * math.exp(-tau)
        G[i] = 2.0 * total
    return G


def compute_W_asym_basis(rho, tau, N_g, N_s, N_i, *, basis,
                          epsabs=1e-13, epsrel=1e-11):
    """Build W matrix with asymptote-basis inner."""
    muc = mu_crit(rho)
    D = N_g + N_s + N_i
    W = np.zeros((D, D))

    leg_s = _unit_legendre_lambdas(max(N_s, 1))
    leg_g = grazing_lambdas(max(N_g, 1), muc) if N_g > 0 else []

    g0 = 0
    s0 = N_g
    i0 = N_g + N_s

    # W_gg
    for m in range(N_g):
        for n in range(N_g):
            def integrand(mu, m=m, n=n):
                return leg_g[m](mu) * leg_g[n](mu) * math.exp(-2.0 * tau * mu) * mu
            val, _ = integrate.quad(integrand, 0.0, muc, epsabs=epsabs, epsrel=epsrel)
            W[g0 + m, g0 + n] = val

    # W_si: outer-steep → inner-arrive
    for m in range(N_i):
        for n in range(N_s):
            def integrand(c, m=m, n=n):
                return basis[m](c) * leg_s[n](c) * \
                       math.exp(-tau * chord_oi(c, rho)) * c
            val, _ = integrate.quad(integrand, 0.0, 1.0,
                                     epsabs=epsabs, epsrel=epsrel)
            W[i0 + m, s0 + n] = val / rho

    # W_is: inner-emit → outer-steep-arrive
    for m in range(N_s):
        for n in range(N_i):
            def integrand(c, m=m, n=n):
                return leg_s[m](c) * basis[n](c) * \
                       math.exp(-tau * chord_oi(c, rho)) * c
            val, _ = integrate.quad(integrand, 0.0, 1.0,
                                     epsabs=epsabs, epsrel=epsrel)
            W[s0 + m, i0 + n] = val / rho

    return W


def run_asym(r_0, R, sig_t_val, sig_s_val, nsf_val, N_i, beta,
              n_panels=2, p_order=4, n_ang=32, dps=15):
    """Run split-basis rank-(1,1,N_i) with asymptote inner basis."""
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
    basis = build_asymptote_basis_numeric(N_i, beta, rho_val)

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

    # Inner (asymptote basis)
    for m in range(N_i):
        P_m = compute_P_esc_inner_asym(geom, r_nodes, radii, sig_t, m, basis,
                                         n_angular=n_ang, dps=dps)
        G_m = compute_G_bc_inner_asym(geom, r_nodes, radii, sig_t, m, basis,
                                        n_surf_quad=n_ang, dps=dps)
        P[2 + m, :] = rv * r_wts * P_m
        G[:, 2 + m] = G_m / A_inner

    # sig_t scaling on G
    sig_t_n = np.empty(N_r)
    for i, ri in enumerate(r_nodes):
        ki = geom.which_annulus(ri, radii)
        sig_t_n[i] = sig_t[ki]
    G = sig_t_n[:, None] * G

    # W
    tau = sig_t_val * R_out
    W = compute_W_asym_basis(rho_val, tau, N_g, N_s, N_i, basis=basis)

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
    print("E3.2 — Physics-asymptote basis: mode-0 = exp(-β·s(c;ρ)/2)")
    print("k_inf = 1.5")
    print("=" * 78)

    sig_t = 1.0
    sig_s = 1.0 / 3.0
    nsf = 1.0

    # Quick probe at the E2.6 sweet spot first
    sig_t_Rs = [1.0, 2.5, 5.0, 10.0, 20.0, 50.0]
    rhos = [0.1, 0.3, 0.5, 0.7]

    # β values scanned: 0, τ/2, τ, 2τ, and also 3τ and 5τ (over-damped)
    print(f"\n{'σ_t·R':>8} {'ρ':>6} {'F.4':>10} {'β=0 (Leg)':>12} "
          f"{'β=τ/2':>12} {'β=τ':>12} {'β=2τ':>12} {'β=3τ':>12} {'β=5τ':>12}")
    print("-" * 110)
    best_by_point = {}
    for sig_t_R in sig_t_Rs:
        for rho in rhos:
            R = sig_t_R / sig_t
            r_0 = rho * R
            try:
                k_f4 = run_scalar_f4(r_0, R, sig_t, sig_s, nsf)
                err_f4 = abs(k_f4 - K_INF) / K_INF * 100
            except Exception as e:
                err_f4 = float("nan")

            tau = sig_t_R
            betas = [0.0, tau / 2.0, tau, 2.0 * tau, 3.0 * tau, 5.0 * tau]
            errs = []
            for beta in betas:
                try:
                    k = run_asym(r_0, R, sig_t, sig_s, nsf, 1, beta)
                    err = abs(k - K_INF) / K_INF * 100
                except Exception as e:
                    err = float("nan")
                errs.append(err)
            finite = [(b, e) for b, e in zip(betas, errs)
                      if not math.isnan(e) and not math.isinf(e)]
            if finite:
                b_opt, e_opt = min(finite, key=lambda x: x[1])
                best_by_point[(sig_t_R, rho)] = (b_opt, e_opt, err_f4)
            line = f"{sig_t_R:>8.1f} {rho:>6.2f} {err_f4:>10.4f}"
            for e in errs:
                if math.isnan(e):
                    line += f"{'NaN':>13}"
                else:
                    line += f"{e:>12.4f}"
            print(line)

    print("\n" + "=" * 78)
    print("Best β per point (β_opt within {0, τ/2, τ, 2τ, 3τ, 5τ}):")
    print("=" * 78)
    print(f"{'σ_t·R':>8} {'ρ':>6} {'β_opt/τ':>10} {'err_asym':>12} {'err_F4':>12} {'ratio':>8}")
    for (sig_t_R, rho), (b_opt, e_opt, e_f4) in best_by_point.items():
        tau = sig_t_R
        beta_norm = b_opt / tau if tau > 0 else 0
        ratio = e_opt / e_f4 if e_f4 > 0 else float("inf")
        print(f"{sig_t_R:>8.1f} {rho:>6.2f} {beta_norm:>10.2f} "
              f"{e_opt:>12.4f} {e_f4:>12.4f} {ratio:>8.2f}")

    # STOP RULE: if β=τ < 0.01% across σ_t·R ∈ {1, 2.5, 5, 10, 20} at ρ=0.3
    rho_target = 0.3
    beta_tau_errs = []
    for sig_t_R in [1.0, 2.5, 5.0, 10.0, 20.0]:
        if (sig_t_R, rho_target) in best_by_point:
            # find β=τ run
            R = sig_t_R / sig_t
            r_0 = rho_target * R
            try:
                k = run_asym(r_0, R, sig_t, sig_s, nsf, 1, sig_t_R)
                err = abs(k - K_INF) / K_INF * 100
            except Exception:
                err = float("inf")
            beta_tau_errs.append((sig_t_R, err))
    print("\nStop rule check (β=τ at ρ=0.3):")
    for s, e in beta_tau_errs:
        print(f"  σ_t·R={s}: err={e:.4f}%")
    all_under = all(e < 0.01 for s, e in beta_tau_errs)
    if all_under:
        print("  *** STOP RULE TRIGGERED: β=τ gives <0.01% universally at ρ=0.3! ***")
    else:
        print("  Stop rule NOT triggered. Proceed with scan / direction E3.4 or E3.5.")

    return 0


if __name__ == "__main__":
    sys.exit(main())
