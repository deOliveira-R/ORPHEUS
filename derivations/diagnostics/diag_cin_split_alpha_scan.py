"""Diagnostic E3.1: brute-force α-scan for Jacobi inner basis.

Created by numerics-investigator on 2026-04-21.
If this test catches a real bug, promote to the matching per-module
folder — tests/cp/ for closure research tests.

Question: does the optimal Jacobi-weight α in the inner-basis weight
function c^(α+1) follow a simple pattern across (σ_t·R, r_0/R)?

We measure α_opt(σ_t·R, r_0/R) as argmin_α |err(rank-(1,1,1), α)|.

If α_opt ≈ f(σ_t·R · (1-ρ)) or similar, we have a principled adaptive
basis. If α_opt is erratic, we need a different approach (see E3.2).
"""
from __future__ import annotations

import math
import sys

import numpy as np
import sympy as sp
from scipy import integrate

sys.path.insert(0, '/workspaces/ORPHEUS/derivations/diagnostics')

from diag_cin_split_jacobi_inner import (
    jacobi_c_weighted_symbolic,
    run_split_jacobi,
)

K_INF = 1.5


def jacobi_inner_lambdas_alpha(n_max, alpha):
    """Jacobi-weighted orthonormal basis on c ∈ [0, 1] with weight c^(α+1).

    (α+1 because the "baseline" c-weight = Legendre has α=0 here in
    our convention.)
    """
    c = sp.Symbol("c", positive=True, real=True)
    weight = lambda x: x ** (float(alpha) + 1.0)
    exprs = jacobi_c_weighted_symbolic(n_max, alpha, 0, 0, 1, c, weight)
    return [sp.lambdify(c, e, modules=["numpy"]) for e in exprs]


def run_with_alpha(r_0, R, sig_t, sig_s, nsf, N_i, alpha,
                   n_panels=2, p_order=4, n_ang=32):
    """Run split-basis rank-(1,1,N_i) with generic Jacobi inner basis (α)."""
    # The run_split_jacobi function accepts (alpha, beta) with discrete cases;
    # we need a monkey-patched version that uses arbitrary α.
    import diag_cin_split_jacobi_inner as mod

    # Save originals
    orig_fn = mod.jacobi_inner_lambdas

    def patched(n_max, alpha=alpha, beta=0):
        return jacobi_inner_lambdas_alpha(n_max, alpha)

    mod.jacobi_inner_lambdas = patched
    try:
        k = run_split_jacobi(
            r_0, R, sig_t, sig_s, nsf, N_i,
            alpha=alpha, beta=0,
            n_panels=n_panels, p_order=p_order, n_ang=n_ang,
        )
    finally:
        mod.jacobi_inner_lambdas = orig_fn
    return k


def scan_alpha_at_point(sig_t_R, rho, alpha_values, sig_t=1.0,
                         n_panels=2, p_order=4, n_ang=32):
    """For a single (σ_t·R, ρ) point, find argmin_α |err|."""
    R = sig_t_R / sig_t
    r_0 = rho * R
    sig_s = sig_t / 3.0
    nsf = sig_t
    # K_inf = nsf / (sig_t - sig_s) = 1 / (1 - 1/3) = 1.5

    results = []
    for alpha in alpha_values:
        try:
            k = run_with_alpha(r_0, R, sig_t, sig_s, nsf, N_i=1, alpha=alpha,
                               n_panels=n_panels, p_order=p_order, n_ang=n_ang)
            err = abs(k - K_INF) / K_INF * 100.0
        except Exception as e:
            err = float("nan")
        results.append((alpha, err))
    return results


def main():
    print("=" * 78)
    print("E3.1 — α-scan for Jacobi inner basis (rank-(1,1,1))")
    print("=" * 78)

    # Fine scan grid — prioritize the key σ_t·R points first
    # Budget: ~6 (σ_t·R) × 3 (ρ) × ~15 (α) = 270 runs at ~1s each = ~4 min
    sig_t_Rs = [1.0, 2.5, 5.0, 10.0, 20.0, 50.0]
    rhos = [0.1, 0.3, 0.5, 0.7]
    alpha_values = [0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.5, 3.0]

    print(f"\nAlpha values: {alpha_values}")
    print(f"σ_t·R values: {sig_t_Rs}")
    print(f"ρ = r_0/R values: {rhos}\n")

    # Store results
    opt_alpha = {}  # (σ_t·R, ρ) -> (α_opt, err_opt)
    all_results = {}

    print("=" * 78)
    print(f"{'σ_t·R':>8} {'ρ':>6}   " + " ".join([f"α={a:<5.2f}" for a in alpha_values]))
    print("=" * 78)

    for sig_t_R in sig_t_Rs:
        for rho in rhos:
            res = scan_alpha_at_point(sig_t_R, rho, alpha_values)
            all_results[(sig_t_R, rho)] = res
            errs = [e for _, e in res]
            # Find best α
            finite = [(a, e) for a, e in res if not math.isnan(e) and not math.isinf(e)]
            if finite:
                a_opt, e_opt = min(finite, key=lambda x: x[1])
                opt_alpha[(sig_t_R, rho)] = (a_opt, e_opt)
            # Print err(α) line
            line = f"{sig_t_R:>8.1f} {rho:>6.2f}   " + " ".join(
                [f"{e:>6.3f}" if not math.isnan(e) else "  NaN " for _, e in res]
            )
            print(line)

    # Summary table
    print("\n" + "=" * 78)
    print("Optimal α per point:")
    print("=" * 78)
    print(f"{'σ_t·R':>8} {'ρ':>6} {'α_opt':>8} {'err_opt':>12}")
    for (sig_t_R, rho), (a_opt, e_opt) in opt_alpha.items():
        print(f"{sig_t_R:>8.1f} {rho:>6.2f} {a_opt:>8.2f} {e_opt:>12.4f}%")

    # Pattern analysis: does α_opt correlate with some combination of (σ_t·R, ρ)?
    print("\n" + "=" * 78)
    print("Pattern analysis:")
    print("=" * 78)
    print(f"{'σ_t·R':>8} {'ρ':>6} {'τ(1-ρ)':>10} {'τρ':>8} {'α_opt':>8}")
    for (sig_t_R, rho), (a_opt, _) in opt_alpha.items():
        print(f"{sig_t_R:>8.1f} {rho:>6.2f} {sig_t_R*(1-rho):>10.2f} "
              f"{sig_t_R*rho:>8.2f} {a_opt:>8.2f}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
