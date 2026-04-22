"""Diagnostic E3.5: source decomposition + metric-inconsistency probe.

Created by numerics-investigator on 2026-04-22.

STAGE 1 finding (from E3.4 v2): the self-consistent ψ^+_inner with the
Legendre basis is a FIXED POINT of the Galerkin iteration. Using the
truth-shape as mode-0 gives the SAME k_eff. So basis-shape changes
alone (via Galerkin) cannot break the plateau.

So where does Jacobi-c²'s 0.072% come from?

HYPOTHESIS: the Jacobi c² "win" at σ_t·R=5 comes from a metric
INCONSISTENCY — the basis is declared orthonormal under c²-weight but
the coupling integrals (W_si, etc.) use c-weight. This is the SAME
type of mismatch as F.4's Lambert-P/G + Marshak-W trick: a formal
inconsistency that accidentally cancels the leading error.

TEST 1: replace Jacobi c² with a constant basis sqrt(3) (same as
Jacobi mode-0 but "Legendre" in shape). Does it give the same 0.072%?
If yes, the Jacobi win is about NORMALIZATION not shape.

TEST 2: run the Jacobi basis with CONSISTENT c²-weight on the
coupling integrals. Does it become catastrophic or match Legendre?

TEST 3: Decompose the ψ^+_inner contribution into (a) the part driven
by white-BC-only (vacuum volume), (b) the part driven by volume-emission-
only (black-BC).
"""
from __future__ import annotations

import math
import sys

import numpy as np
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


def make_constant_basis(value):
    """Return a single-mode basis of constant function = value."""
    def f(c):
        if np.isscalar(c):
            return float(value)
        return np.full_like(c, float(value), dtype=float)
    return [f]


def run_custom_basis(r_0, R, sig_t_val, sig_s_val, nsf_val, basis,
                      n_panels=2, p_order=4, n_ang=32, dps=15,
                      coupling_weight_power=1):
    """Run split-basis rank-(1,1,N_i=len(basis)) with custom inner basis.

    coupling_weight_power: power of c in the c-weighted integrand of W_si, W_is
    and the G/P integrals. 1 = standard c-weight.
    """
    N_i = len(basis)
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

    # Inner
    from diag_cin_split_asymptote_basis import (
        compute_P_esc_inner_asym,
        compute_G_bc_inner_asym,
    )
    for m in range(N_i):
        P_m = compute_P_esc_inner_asym(geom, r_nodes, radii, sig_t, m, basis,
                                         n_angular=n_ang, dps=dps)
        G_m = compute_G_bc_inner_asym(geom, r_nodes, radii, sig_t, m, basis,
                                        n_surf_quad=n_ang, dps=dps)
        P[2 + m, :] = rv * r_wts * P_m
        G[:, 2 + m] = G_m / A_inner

    sig_t_n = np.empty(N_r)
    for i, ri in enumerate(r_nodes):
        ki = geom.which_annulus(ri, radii)
        sig_t_n[i] = sig_t[ki]
    G = sig_t_n[:, None] * G

    # W — with custom coupling weight power
    tau = sig_t_val * R_out
    W = np.zeros((D, D))
    # W_gg
    for m in range(N_g):
        for n in range(N_g):
            def integrand(mu, m=m, n=n):
                return leg_g[m](mu) * leg_g[n](mu) * math.exp(-2.0 * tau * mu) * mu
            val, _ = integrate.quad(integrand, 0.0, muc, epsabs=1e-13, epsrel=1e-11)
            W[m, n] = val

    # W_si
    for m in range(N_i):
        for n in range(N_s):
            def integrand(c, m=m, n=n):
                return basis[m](c) * leg_s[n](c) * \
                       math.exp(-tau * chord_oi(c, rho_val)) * c ** coupling_weight_power
            val, _ = integrate.quad(integrand, 0.0, 1.0, epsabs=1e-13, epsrel=1e-11)
            W[2 + m, 1] = val / rho_val

    # W_is
    for m in range(N_s):
        for n in range(N_i):
            def integrand(c, m=m, n=n):
                return leg_s[m](c) * basis[n](c) * \
                       math.exp(-tau * chord_oi(c, rho_val)) * c ** coupling_weight_power
            val, _ = integrate.quad(integrand, 0.0, 1.0, epsabs=1e-13, epsrel=1e-11)
            W[1, 2 + n] = val / rho_val

    B = np.zeros((D, D))
    B[0, 0] = muc * muc
    B[0, 1] = muc * rho_val
    B[1, 0] = rho_val * muc
    B[1, 1] = rho_val * rho_val
    for m in range(N_i):
        B[2 + m, 2 + m] = 1.0

    M = np.eye(D) - W @ B
    K_bc = G @ B @ np.linalg.inv(M) @ P
    K = K_vol + K_bc
    return solve_k_eff(K, sig_t_val, sig_s_val, nsf_val)


def main():
    print("=" * 78)
    print("E3.5 — Source decomposition + metric inconsistency probe")
    print("σ_t·R = 5, r_0/R = 0.3 (the Jacobi-c² sweet spot)")
    print("=" * 78)

    sig_t = 1.0
    sig_s = 1.0 / 3.0
    nsf = 1.0
    R = 5.0
    r_0 = 0.3 * R

    k_f4 = run_scalar_f4(r_0, R, sig_t, sig_s, nsf)
    err_f4 = abs(k_f4 - K_INF) / K_INF * 100
    print(f"\nF.4 err: {err_f4:.4f}%")

    # --- TEST 1: constant basis at various values ---
    print("\nTEST 1: constant mode-0 basis at various scalar values")
    print("  (all c-weight coupling, varying normalization)")
    for scale in [math.sqrt(2.0), math.sqrt(3.0), math.sqrt(2.5), 1.0,
                    math.sqrt(4.0), math.sqrt(1.5)]:
        basis = make_constant_basis(scale)
        k = run_custom_basis(r_0, R, sig_t, sig_s, nsf, basis,
                              coupling_weight_power=1)
        err = abs(k - K_INF) / K_INF * 100
        print(f"  scale={scale:.4f}: k={k:.8f}, err={err:.4f}%")

    print("\n  ** If all give same err, normalization scalar cancels **")
    print("  ** If different, W_si scale is NOT absorbed by B or coefficients **")

    # --- TEST 2: constant basis with different coupling weight powers ---
    print("\nTEST 2: constant mode-0 = sqrt(2) with different coupling weights")
    print("  (c^power in W_si, W_is integrals)")
    for power in [0.0, 0.5, 1.0, 1.5, 2.0]:
        basis = make_constant_basis(math.sqrt(2.0))
        k = run_custom_basis(r_0, R, sig_t, sig_s, nsf, basis,
                              coupling_weight_power=power)
        err = abs(k - K_INF) / K_INF * 100
        print(f"  c^{power}: k={k:.8f}, err={err:.4f}%")

    # --- TEST 3: various Jacobi mode-0 shapes ---
    print("\nTEST 3: c^α shapes as mode-0 (raw, not orthonormal)")
    # Use `mode_0(c) = c^α / norm` under c-weight orthonormality.
    # α=0: constant. α=1: linear.  Higher α: steeper in c.
    x_gl, w_gl = np.polynomial.legendre.leggauss(300)
    c_grid = 0.5 * (x_gl + 1.0)
    w_cgrid = 0.5 * w_gl * c_grid
    for alpha in [0.0, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0]:
        def mode0(c, alpha=alpha):
            return c ** alpha
        # c-weight normalization
        vals = mode0(c_grid)
        norm = math.sqrt(np.sum(w_cgrid * vals * vals))
        def normed(c, alpha=alpha, norm=norm):
            if np.isscalar(c):
                return float((c ** alpha) / norm)
            return (c ** alpha) / norm
        basis = [normed]
        k = run_custom_basis(r_0, R, sig_t, sig_s, nsf, basis,
                              coupling_weight_power=1)
        err = abs(k - K_INF) / K_INF * 100
        print(f"  mode0 = c^{alpha:.2f}/norm: k={k:.8f}, err={err:.4f}%")

    # --- TEST 4: two modes — Legendre's mode-0 PLUS Jacobi's mode-0 (linearly independent?) ---
    # Jacobi c² mode-0 = sqrt(3) constant → same direction as Legendre mode-0.
    # So 2-mode basis (const, c) is what matters, and Legendre vs Jacobi mode-0
    # is the SAME span when N_i=1.  But numerically they give different k_eff.
    # That confirms metric inconsistency in coupling.
    # Instead test: use {1, c} as both modes under c-weight orthonormal,
    # compare with {sqrt(3) constant, sqrt(3)·c · alpha-shifted} under c²-weight.

    return 0


if __name__ == "__main__":
    sys.exit(main())
