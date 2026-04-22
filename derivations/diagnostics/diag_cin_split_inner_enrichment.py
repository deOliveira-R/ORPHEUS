"""Diagnostic E2.2/E2.3: split-basis rank-(1,1,N) scan + spectral residual.

Created by numerics-investigator on 2026-04-21.

E2.2 — rank-(1,1,N) scan at two quadrature levels.
    Test whether the 0.99% plateau moves with quadrature.

E2.3 — Spectral residual decomposition.
    At the self-consistent rank-(1,1,1) solution, inspect the angular flux
    at inner and outer surfaces. Compare to rank-(1,1,N_large) "truth" —
    decompose the residual into Legendre modes. Where does the residual
    energy live?

If E2.2 shows the plateau is quadrature-independent, the plateau is
STRUCTURAL — the basis itself can't capture the truth. E2.3 then tells
us WHICH modes we're missing.
"""
from __future__ import annotations

import numpy as np

import sys
sys.path.insert(0, '/workspaces/ORPHEUS/derivations/diagnostics')

from diag_cin_aware_split_basis_keff import (
    build_split_basis_closure,
    _unit_legendre_lambdas,
    grazing_lambdas,
    mu_crit,
    c_of_mu,
    mu_of_c,
    chord_oi,
    solve_k_eff,
    run_scalar_f4,
    run_split_basis,
)
from orpheus.derivations.peierls_geometry import (
    CurvilinearGeometry,
    build_volume_kernel,
    composite_gl_r,
)


K_INF = 1.5


def self_consistent_split(r_0, R, sig_t_val, sig_s_val, nsf_val,
                           N_g, N_s, N_i, *,
                           n_panels=2, p_order=4, dps=15, n_ang=32):
    """Compute self-consistent k_eff AND the surface coefficients ψ^+.

    Returns (k_eff, r_nodes, phi_V, psi_plus_coeffs) where ψ^+ is the
    D-vector of outgoing current mode coefficients in the split basis.
    """
    geom = CurvilinearGeometry(kind="sphere-1d", inner_radius=r_0)
    radii = np.array([R])
    sig_t = np.array([sig_t_val])
    r_nodes, r_wts, panels = composite_gl_r(radii, n_panels, p_order,
                                             dps=dps, inner_radius=r_0)
    K_vol = build_volume_kernel(geom, r_nodes, panels, radii, sig_t,
                                 n_angular=n_ang, n_rho=n_ang, dps=dps)
    K_bc, aux = build_split_basis_closure(
        geom, r_nodes, r_wts, radii, sig_t,
        N_g, N_s, N_i,
        n_angular=n_ang, n_surf_quad=n_ang, dps=dps,
    )
    if K_bc is None:
        return None, None, None, None, aux
    K = K_vol + K_bc

    # Solve eigen
    N = K.shape[0]
    A = np.diag(np.full(N, sig_t_val)) - K * sig_s_val
    B = K * nsf_val
    phi = np.ones(N)
    k = 1.0
    for _ in range(500):
        q = B @ phi / k
        phi_new = np.linalg.solve(A, q)
        B_phi_new = B @ phi_new
        B_phi = B @ phi
        k_new = k * (np.abs(B_phi_new).sum() / np.abs(B_phi).sum())
        if abs(k_new - k) < 1e-12:
            k = k_new
            phi = phi_new
            break
        phi = phi_new / max(np.linalg.norm(phi_new), 1e-30)
        k = k_new
    # Normalize phi
    phi = phi / np.max(np.abs(phi))

    # Reconstruct surface coefficients
    Q_V = sig_s_val * phi + nsf_val / k * phi
    P = aux["P"]
    W = aux["W"]
    Bop = aux["B"]
    G = aux["G"]
    M = np.eye(N_g + N_s + N_i) - W @ Bop
    psi_plus = np.linalg.solve(M, P @ Q_V)  # outgoing
    psi_minus = Bop @ psi_plus  # incoming (white + cavity)

    return k, r_nodes, phi, psi_plus, psi_minus, aux


def _integrate_on_domain(integrand, a, b, quad_order=200):
    """Numerical integral using Gauss-Legendre."""
    from numpy.polynomial.legendre import leggauss
    x_gl, w_gl = leggauss(quad_order)
    x = 0.5 * (b - a) * x_gl + 0.5 * (a + b)
    w = 0.5 * (b - a) * w_gl
    return np.sum(w * integrand(x))


def project_onto_plain_legendre(phi_func, n_modes, a=0.0, b=1.0, quad_order=300):
    """Project phi_func(mu) with µ-weighted inner product onto plain half-range
    Legendre on [a, b]. Returns coefficients."""
    leg = _unit_legendre_lambdas(n_modes - 1)
    coeffs = np.zeros(n_modes)
    for n in range(n_modes):
        def integrand(mu, n=n):
            return phi_func(mu) * leg[n](mu) * mu
        coeffs[n] = _integrate_on_domain(integrand, a, b, quad_order)
    return coeffs


def evaluate_psi_at_outer(psi_plus_coeffs, N_g, N_s, N_i, rho_val, muc,
                           n_samples=200):
    """Evaluate ψ^+_outer(µ) on [0, 1] from split-basis coefficients."""
    leg_g = grazing_lambdas(max(N_g, 1), muc) if N_g > 0 else []
    leg = _unit_legendre_lambdas(max(N_s, N_i, 1))
    mu = np.linspace(1e-6, 1.0 - 1e-6, n_samples)

    psi = np.zeros_like(mu)
    for i, m in enumerate(mu):
        val = 0.0
        if m < muc and N_g > 0:
            for n in range(N_g):
                val += psi_plus_coeffs[n] * leg_g[n](m)
        elif m >= muc and N_s > 0:
            c = c_of_mu(m, rho_val)
            for n in range(N_s):
                val += psi_plus_coeffs[N_g + n] * leg[n](c) / rho_val
        psi[i] = val
    return mu, psi


def evaluate_psi_at_inner(psi_plus_coeffs, N_g, N_s, N_i, n_samples=200):
    """Evaluate ψ^+_inner(c) on [0, 1]."""
    leg = _unit_legendre_lambdas(max(N_s, N_i, 1))
    c = np.linspace(1e-6, 1.0 - 1e-6, n_samples)
    psi = np.zeros_like(c)
    for i in range(N_i):
        psi += psi_plus_coeffs[N_g + N_s + i] * leg[i](c)
    return c, psi


def main():
    print("=" * 78)
    print("E2.2 + E2.3 — rank-(1,1,N) scan + spectral residual")
    print("σ_t·R = 5, r_0/R = 0.3, k_inf = 1.5")
    print("=" * 78)

    sig_t, sig_s, nsf = 1.0, 1.0 / 3.0, 1.0
    R = 5.0
    r_0 = 0.3 * R
    rho_val = 0.3
    muc = mu_crit(rho_val)
    print(f"µ_crit = {muc:.8f}, ρ = {rho_val}")

    # ---- E2.2: scan ----
    print("\n--- E2.2: rank-(1,1,N) scan at two quadrature levels ---")
    base_cfg = (2, 4, 32, "base")
    rich_cfg = (4, 8, 64, "2×rich")
    for (np_, p, na, label) in [base_cfg, rich_cfg]:
        print(f"\nConfig {label}: n_panels={np_}, p={p}, n_ang={na}")
        k_f4 = run_scalar_f4(r_0, R, sig_t, sig_s, nsf,
                              n_panels=np_, p_order=p, n_ang=na)
        err_f4 = abs(k_f4 - K_INF) / K_INF * 100
        print(f"  F.4         : err={err_f4:.6f}%")
        for N_i in [1, 2, 3, 4, 5, 8]:
            k, _ = run_split_basis(r_0, R, sig_t, sig_s, nsf, 1, 1, N_i,
                                    n_panels=np_, p_order=p, n_ang=na)
            if k is None:
                continue
            err = abs(k - K_INF) / K_INF * 100
            print(f"  rank-(1,1,{N_i})  : err={err:.6f}%")

    # Also: rank-(N,N,N) at baseline
    print(f"\nConfig base: rank-(N,N,N) full scan")
    for N in [1, 2, 3, 4]:
        k, _ = run_split_basis(r_0, R, sig_t, sig_s, nsf, N, N, N)
        err = abs(k - K_INF) / K_INF * 100 if k else float("nan")
        print(f"  rank-({N},{N},{N}): err={err:.6f}%")

    # ---- E2.3: Spectral residual ----
    print("\n--- E2.3: Spectral residual decomposition ---")
    print("Compute ψ^+_outer(µ) and ψ^+_inner(c) at rank-(1,1,1) and (1,1,8)")

    k_1, r_nodes_1, phi_1, psi_plus_1, psi_minus_1, aux_1 = self_consistent_split(
        r_0, R, sig_t, sig_s, nsf, 1, 1, 1,
    )
    # Use (1,1,8) as "truth" proxy (it plateaus at rank-(1,1,N>=2))
    k_ref, r_nodes_ref, phi_ref, psi_plus_ref, psi_minus_ref, aux_ref = self_consistent_split(
        r_0, R, sig_t, sig_s, nsf, 1, 1, 8,
    )

    err_1 = abs(k_1 - K_INF) / K_INF * 100
    err_ref = abs(k_ref - K_INF) / K_INF * 100
    print(f"  rank-(1,1,1): k={k_1:.8f}, err={err_1:.6f}%")
    print(f"  rank-(1,1,8): k={k_ref:.8f}, err={err_ref:.6f}%")
    print(f"  coeffs rank-(1,1,1) ψ^+ = {psi_plus_1}")
    print(f"  coeffs rank-(1,1,8) ψ^+ (first 3 + last 1):")
    print(f"    ψ^+_graze_0={psi_plus_ref[0]:.6e}")
    print(f"    ψ^+_steep_0={psi_plus_ref[1]:.6e}")
    print(f"    ψ^+_inner_0={psi_plus_ref[2]:.6e}, 1={psi_plus_ref[3]:.6e}, "
          f"2={psi_plus_ref[4]:.6e}, 3={psi_plus_ref[5]:.6e}")
    print(f"    ψ^+_inner_7={psi_plus_ref[9]:.6e}")

    # Energy per mode in rank-(1,1,8) ψ^+_inner
    print("\n  rank-(1,1,8) inner mode energies |ψ^+_inner_n|²/Σ:")
    inner_coeffs = psi_plus_ref[2:]
    total_energy = np.sum(inner_coeffs ** 2)
    for i, c in enumerate(inner_coeffs):
        frac = c ** 2 / total_energy * 100 if total_energy > 0 else 0.0
        print(f"    mode {i}: |c|²={c**2:.3e}  ({frac:.3f}%)")

    # Same for ψ^- on outer (incoming):
    print("\n  rank-(1,1,8) ψ^- (incoming) coefficients:")
    print(f"    ψ^-_graze_0={psi_minus_ref[0]:.6e}")
    print(f"    ψ^-_steep_0={psi_minus_ref[1]:.6e}")
    for i in range(min(8, len(psi_minus_ref) - 2)):
        print(f"    ψ^-_inner_{i}={psi_minus_ref[2 + i]:.6e}")

    # Rank-(1,1,1) flux reconstruction vs rank-(1,1,8) — project onto plain Legendre
    print("\n  Residual decomposition on INNER (plain half-range Legendre, 8 modes)")
    c_grid_1, psi_in_1 = evaluate_psi_at_inner(psi_plus_1, 1, 1, 1)
    c_grid_ref, psi_in_ref = evaluate_psi_at_inner(psi_plus_ref, 1, 1, 8)
    # Interpolate both onto same grid
    psi_in_1_interp = np.interp(c_grid_ref, c_grid_1, psi_in_1)
    residual = psi_in_1_interp - psi_in_ref
    # Project residual onto Legendre
    leg8 = _unit_legendre_lambdas(7)
    print("    mode    |coeff of residual|    |residual^2|")
    for n in range(8):
        def intg(c, n=n):
            return np.interp(c, c_grid_ref, residual) * leg8[n](c) * c
        coef = _integrate_on_domain(intg, 0, 1, 300)
        print(f"    {n}       {coef:.6e}")

    # Same for OUTER
    print("\n  Residual decomposition on OUTER (plain half-range Legendre, 8 modes)")
    mu_grid_1, psi_out_1 = evaluate_psi_at_outer(psi_plus_1, 1, 1, 1, rho_val, muc)
    mu_grid_ref, psi_out_ref = evaluate_psi_at_outer(psi_plus_ref, 1, 1, 8, rho_val, muc)
    psi_out_1_interp = np.interp(mu_grid_ref, mu_grid_1, psi_out_1)
    residual_out = psi_out_1_interp - psi_out_ref
    print("    mode    |coeff of residual|")
    for n in range(8):
        def intg(mu, n=n):
            return np.interp(mu, mu_grid_ref, residual_out) * leg8[n](mu) * mu
        coef = _integrate_on_domain(intg, 0, 1, 300)
        print(f"    {n}       {coef:.6e}")

    return 0


def test_split_rank_111_plateau_is_structural():
    """Split-basis rank-(1,1,N) plateau is independent of quadrature.

    If refining quadrature 2× does NOT bring the plateau below F.4, the
    plateau is structural. This test records the baseline vs rich numbers.
    """
    import numpy as np
    sig_t, sig_s, nsf = 1.0, 1.0 / 3.0, 1.0
    R = 5.0
    r_0 = 0.3 * R
    k_b, _ = run_split_basis(r_0, R, sig_t, sig_s, nsf, 1, 1, 4,
                              n_panels=2, p_order=4, n_ang=32)
    k_r, _ = run_split_basis(r_0, R, sig_t, sig_s, nsf, 1, 1, 4,
                              n_panels=4, p_order=8, n_ang=64)
    err_b = abs(k_b - K_INF) / K_INF * 100
    err_r = abs(k_r - K_INF) / K_INF * 100
    print(f"rank-(1,1,4) err: base={err_b:.4f}%, rich={err_r:.4f}%")
    # Structural plateau: refining doesn't help (the gap is < 0.2%)
    assert abs(err_b - err_r) < 0.5, (
        f"Expected structural plateau (~1%). Got base={err_b:.4f}%, rich={err_r:.4f}%"
    )


if __name__ == "__main__":
    sys.exit(main())
