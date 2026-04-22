"""Diagnostic E3.4 v2: Galerkin-adaptive with richer basis seed.

Created by numerics-investigator on 2026-04-22.

LESSON from v1: at rank-(1,1,1) with Gram-Schmidt normalization, the
Legendre iteration is a fixed point (mode-0 = constant is preserved
by Gram-Schmidt with a constant seed). To break this, we need to
ITERATE THE BASIS SHAPE with richer information — i.e., seed with N_i>=2
and use the resulting ψ^+_inner shape as the new mode-0 for a
subsequent rank-(1,1,1) solve.

PROTOCOL:
  1. Solve rank-(1,1,N_large) with Legendre inner. Extract ψ^+_inner.
  2. Use that ψ^+_inner as the MODE-0 of a NEW basis, rank-(1,1,1).
  3. Compare k_eff to rank-(1,1,N_large) result.

If the rank-(1,1,1) with truth-shape mode-0 gives the SAME k_eff as
the rank-(1,1,N_large), we've proven the basis matters, not the rank.

Then iterate: use the new ψ^+_inner (from rank-(1,1,1) with truth-shape)
as the new mode-0 for next iteration. Does it converge to a fixed point?
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
from diag_cin_split_galerkin_adaptive import (
    build_basis_from_mode0_vals,
    vals_to_callable,
    solve_and_extract_psi_inner,
)
from orpheus.derivations.peierls_geometry import (
    CurvilinearGeometry,
    build_volume_kernel,
    composite_gl_r,
    gl_float,
)


K_INF = 1.5


def solve_rank_11N_extract_shape(r_0, R, sig_t_val, sig_s_val, nsf_val, N_i,
                                    c_grid, w_grid,
                                    n_panels=2, p_order=4, n_ang=32, dps=15):
    """Solve rank-(1,1,N_i) with standard Legendre half-range inner basis.

    Extract both k_eff and the self-consistent ψ^+_inner(c) values.
    """
    leg = _unit_legendre_lambdas(max(N_i, 2))
    basis = [leg[m] for m in range(N_i)]
    k, psi_inner_grid, coeffs = solve_and_extract_psi_inner(
        r_0, R, sig_t_val, sig_s_val, nsf_val, basis, N_i, c_grid, w_grid,
        n_panels=n_panels, p_order=p_order, n_ang=n_ang, dps=dps,
    )
    return k, psi_inner_grid, coeffs


def main():
    print("=" * 78)
    print("E3.4 v2 — Rank-(1,1,N) shape → rank-(1,1,1) truth-basis loop")
    print("=" * 78)

    sig_t = 1.0
    sig_s = 1.0 / 3.0
    nsf = 1.0

    n_gram = 300
    x_gl, w_gl = np.polynomial.legendre.leggauss(n_gram)
    c_grid = 0.5 * (x_gl + 1.0)
    w_grid = 0.5 * w_gl * c_grid  # c-weight on [0, 1]

    test_points = [(5.0, 0.3), (10.0, 0.3), (20.0, 0.3), (50.0, 0.3),
                    (2.5, 0.3), (1.0, 0.3),
                    (5.0, 0.1), (5.0, 0.5), (5.0, 0.7),
                    (10.0, 0.1), (10.0, 0.5),
                    (20.0, 0.5),]

    N_large = 8  # capture all significant inner modes

    results = []
    for sig_t_R, rho in test_points:
        R = sig_t_R / sig_t
        r_0 = rho * R
        print(f"\n--- σ_t·R={sig_t_R}, ρ={rho} ---")
        try:
            k_f4 = run_scalar_f4(r_0, R, sig_t, sig_s, nsf)
            err_f4 = abs(k_f4 - K_INF) / K_INF * 100
        except Exception:
            err_f4 = float("nan")

        # Step 1: rank-(1,1,N_large) with Legendre
        try:
            k_large, psi_inner_large, _ = solve_rank_11N_extract_shape(
                r_0, R, sig_t, sig_s, nsf, N_large, c_grid, w_grid,
            )
            err_large = abs(k_large - K_INF) / K_INF * 100
        except Exception as e:
            print(f"  rank-(1,1,{N_large}) FAILED: {type(e).__name__}: {e}")
            continue

        # Step 2: rank-(1,1,1) with psi_inner_large as mode-0
        try:
            basis_vals = build_basis_from_mode0_vals(psi_inner_large, c_grid, w_grid, 1)
            basis = [vals_to_callable(v, c_grid) for v in basis_vals]
            k_adapt, psi_adapt, _ = solve_and_extract_psi_inner(
                r_0, R, sig_t, sig_s, nsf, basis, 1, c_grid, w_grid,
            )
            err_adapt = abs(k_adapt - K_INF) / K_INF * 100
        except Exception as e:
            print(f"  rank-(1,1,1) adaptive FAILED: {type(e).__name__}: {e}")
            continue

        # Step 3: iterate the adaptive basis
        psi_prev = psi_adapt.copy()
        k_prev = k_adapt
        iters_to_converge = 0
        err_iter_final = err_adapt
        for iteration in range(1, 8):
            basis_vals = build_basis_from_mode0_vals(psi_prev, c_grid, w_grid, 1)
            basis = [vals_to_callable(v, c_grid) for v in basis_vals]
            k_new, psi_new, _ = solve_and_extract_psi_inner(
                r_0, R, sig_t, sig_s, nsf, basis, 1, c_grid, w_grid,
            )
            err_new = abs(k_new - K_INF) / K_INF * 100
            if abs(k_new - k_prev) < 1e-10:
                iters_to_converge = iteration
                err_iter_final = err_new
                break
            psi_prev = psi_new.copy()
            k_prev = k_new
            err_iter_final = err_new
            iters_to_converge = iteration

        print(f"  F.4:                      {err_f4:>8.4f}%")
        print(f"  rank-(1,1,{N_large}) Legendre:    {err_large:>8.4f}%")
        print(f"  rank-(1,1,1) truth-mode0: {err_adapt:>8.4f}%")
        print(f"  rank-(1,1,1) converged ({iters_to_converge} iters): {err_iter_final:>8.4f}%")

        results.append((sig_t_R, rho, err_f4, err_large, err_adapt, err_iter_final))

    print("\n" + "=" * 78)
    print("Summary:")
    print("=" * 78)
    print(f"{'σ_t·R':>8} {'ρ':>6} {'F.4':>10} {'rank-(1,1,N)':>14} "
          f"{'adapt':>10} {'conv':>10}")
    for (s, r, ef, el, ea, ec) in results:
        print(f"{s:>8.1f} {r:>6.2f} {ef:>10.4f} {el:>14.4f} {ea:>10.4f} {ec:>10.4f}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
