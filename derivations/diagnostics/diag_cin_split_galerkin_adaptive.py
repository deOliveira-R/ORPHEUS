"""Diagnostic E3.4: Self-consistent iterative basis (Galerkin-adaptive).

Created by numerics-investigator on 2026-04-22.

IDEA: the perfect mode-0 inner basis function IS the true ψ^+_inner(c).
We don't know it a priori — but we can iterate:
    1. Start with Legendre (or Jacobi) inner basis.
    2. Solve at rank-(1,1,N_i) — get ψ^+_inner(c; iter=0).
    3. Use the extracted ψ^+_inner(c) as the NEW mode-0 basis function.
    4. Gram-Schmidt for higher modes.
    5. Re-solve. Check k_eff convergence.

If this converges in few iterations and gives < 0.1%, we have a
practical universal closure (with one-time fit per (τ, ρ)).
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
from diag_cin_split_asymptote_basis import (
    compute_P_esc_inner_asym,
    compute_G_bc_inner_asym,
    compute_W_asym_basis,
)
from orpheus.derivations.peierls_geometry import (
    CurvilinearGeometry,
    build_volume_kernel,
    composite_gl_r,
    gl_float,
)


K_INF = 1.5


def build_basis_from_mode0_vals(mode0_vals, c_grid, w_grid, N_i):
    """Given mode-0 values on c_grid, Gram-Schmidt to get N_i modes.

    mode_n seed: mode0(c) * c^n.
    """
    # Ensure positivity of mode-0 for clean normalization
    # (if mode-0 flips sign we lose interpretation but still a valid basis)
    f0 = mode0_vals.copy()
    # Normalize mode-0 first
    norm = math.sqrt(np.sum(w_grid * f0 * f0))
    if norm <= 0:
        raise RuntimeError("Degenerate mode-0")
    f0 = f0 / norm

    basis_vals = [f0]
    for n in range(1, N_i):
        # Seed: mode0 * c^n
        phi = mode0_vals * c_grid ** n
        for bj in basis_vals:
            cj = np.sum(w_grid * phi * bj)
            phi = phi - cj * bj
        nrm = math.sqrt(np.sum(w_grid * phi * phi))
        if nrm <= 0:
            raise RuntimeError(f"Degenerate mode {n}")
        basis_vals.append(phi / nrm)
    return basis_vals


def vals_to_callable(vals, c_grid):
    from scipy.interpolate import interp1d
    interp = interp1d(c_grid, vals, kind='cubic',
                      bounds_error=False, fill_value="extrapolate")
    def f(c):
        if np.isscalar(c):
            return float(interp(c))
        return interp(c)
    return f


def solve_and_extract_psi_inner(r_0, R, sig_t_val, sig_s_val, nsf_val,
                                  basis, N_i, c_grid, w_grid,
                                  n_panels=2, p_order=4, n_ang=32, dps=15):
    """Solve rank-(1,1,N_i) with given inner basis and extract
    ψ^+_inner(c) values on c_grid."""
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

    tau = sig_t_val * R_out
    W = compute_W_asym_basis(rho_val, tau, N_g, N_s, N_i, basis=basis)

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

    # Solve eigenproblem AND extract ψ^+ coefficients
    N = K.shape[0]
    A_eig = np.diag(np.full(N, sig_t_val)) - K * sig_s_val
    B_eig = K * nsf_val
    phi = np.ones(N)
    k = 1.0
    for _ in range(500):
        q = B_eig @ phi / k
        phi_new = np.linalg.solve(A_eig, q)
        B_phi_new = B_eig @ phi_new
        B_phi = B_eig @ phi
        k_new = k * (np.abs(B_phi_new).sum() / np.abs(B_phi).sum())
        if abs(k_new - k) < 1e-12:
            k = k_new
            phi = phi_new
            break
        phi = phi_new / max(np.linalg.norm(phi_new), 1e-30)
        k = k_new
    phi = phi / np.max(np.abs(phi))

    # ψ^+ = (I - W B)^{-1} P Q_V
    Q_V = (sig_s_val + nsf_val / k) * phi
    psi_plus = np.linalg.solve(M, P @ Q_V)

    # ψ^+_inner(c) on c_grid = Σ_m a_inner_m * basis_m(c)
    psi_inner_coeffs = psi_plus[N_g + N_s:]  # len = N_i
    psi_inner_grid = np.zeros_like(c_grid)
    for m in range(N_i):
        psi_inner_grid = psi_inner_grid + psi_inner_coeffs[m] * basis[m](c_grid)

    return k, psi_inner_grid, psi_inner_coeffs


def main():
    print("=" * 78)
    print("E3.4 — Galerkin-adaptive iterative basis")
    print("=" * 78)

    sig_t = 1.0
    sig_s = 1.0 / 3.0
    nsf = 1.0

    # Setup c grid for basis manipulation
    n_gram = 200
    x_gl, w_gl = np.polynomial.legendre.leggauss(n_gram)
    c_grid = 0.5 * (x_gl + 1.0)
    w_grid = 0.5 * w_gl * c_grid  # c-weight on [0, 1]

    test_points = [(5.0, 0.3), (10.0, 0.3), (20.0, 0.3), (2.5, 0.3),
                    (1.0, 0.3), (50.0, 0.3)]

    for sig_t_R, rho in test_points:
        R = sig_t_R / sig_t
        r_0 = rho * R
        print(f"\n--- σ_t·R = {sig_t_R}, r_0/R = {rho} ---")
        try:
            k_f4 = run_scalar_f4(r_0, R, sig_t, sig_s, nsf)
            err_f4 = abs(k_f4 - K_INF) / K_INF * 100
            print(f"F.4 err: {err_f4:.4f}%")
        except Exception:
            err_f4 = float("nan")
            print(f"F.4 err: NaN (solver crashed)")

        # Iteration 0: Legendre basis (half-range orthonormal)
        # In our convention, this is Gram-Schmidt of 1, c, c², ... with c-weight
        # Start with half-range Legendre mode-0 = sqrt(2)
        leg = _unit_legendre_lambdas(4)
        # Build a reasonable seed
        mode0_seed = np.sqrt(2.0) * np.ones_like(c_grid)
        try:
            basis_vals = build_basis_from_mode0_vals(mode0_seed, c_grid, w_grid, 1)
            basis = [vals_to_callable(v, c_grid) for v in basis_vals]
            k0, psi_inner_0, coeffs_0 = solve_and_extract_psi_inner(
                r_0, R, sig_t, sig_s, nsf, basis, 1, c_grid, w_grid,
            )
            err0 = abs(k0 - K_INF) / K_INF * 100
            print(f"  iter 0 (Legendre inner): k={k0:.8f}, err={err0:.4f}%")
        except Exception as e:
            print(f"  iter 0 FAILED: {type(e).__name__}: {e}")
            continue

        # Iterate: use ψ^+_inner from previous iteration as mode-0
        psi_prev = psi_inner_0.copy()
        k_prev = k0
        for iteration in range(1, 6):
            try:
                basis_vals = build_basis_from_mode0_vals(psi_prev, c_grid, w_grid, 1)
                basis = [vals_to_callable(v, c_grid) for v in basis_vals]
                k_new, psi_inner_new, _ = solve_and_extract_psi_inner(
                    r_0, R, sig_t, sig_s, nsf, basis, 1, c_grid, w_grid,
                )
                err_new = abs(k_new - K_INF) / K_INF * 100
                delta_k = abs(k_new - k_prev)
                delta_psi = np.max(np.abs(psi_inner_new - psi_prev))
                print(f"  iter {iteration}: k={k_new:.8f}, err={err_new:.4f}%, "
                      f"Δk={delta_k:.2e}, Δψ={delta_psi:.2e}")
                if delta_k < 1e-8 and iteration > 1:
                    print(f"  converged at iter {iteration}")
                    break
                psi_prev = psi_inner_new.copy()
                k_prev = k_new
            except Exception as e:
                print(f"  iter {iteration} FAILED: {type(e).__name__}: {e}")
                break

    return 0


if __name__ == "__main__":
    sys.exit(main())
