"""Diagnostic E3.6: scale-scan of constant inner mode-0.

Created by numerics-investigator on 2026-04-22.

E3.5 revealed: varying the constant-inner-basis SCALE affects k_eff.
scale=sqrt(3) ≈ 1.732 at σ_t·R=5, ρ=0.3 gives 0.072% (Jacobi-c²'s win)
but Legendre's sqrt(2) ≈ 1.414 gives 1.29%.

TEST: find scale_opt(σ_t·R, ρ) that minimizes err, and see if there's
a principled formula for it. If scale_opt is predictable from τ and ρ,
we have a universal closure.
"""
from __future__ import annotations

import math
import sys

import numpy as np
from scipy import integrate
from scipy.optimize import minimize_scalar

sys.path.insert(0, '/workspaces/ORPHEUS/derivations/diagnostics')

from diag_cin_aware_split_basis_keff import (
    _unit_legendre_lambdas,
    grazing_lambdas,
    mu_crit,
    chord_oi,
    solve_k_eff,
    run_scalar_f4,
)
from diag_cin_split_source_decomposition import (
    make_constant_basis,
    run_custom_basis,
)

K_INF = 1.5


def run_scale(sig_t_R, rho, scale, sig_t=1.0):
    R = sig_t_R / sig_t
    r_0 = rho * R
    sig_s = sig_t / 3.0
    nsf = sig_t
    basis = make_constant_basis(scale)
    k = run_custom_basis(r_0, R, sig_t, sig_s, nsf, basis,
                          coupling_weight_power=1)
    return abs(k - K_INF) / K_INF * 100


def find_optimal_scale(sig_t_R, rho, scale_range=(0.5, 3.0), tol=1e-4):
    """Find scale that minimizes err."""
    def f(scale):
        try:
            return run_scale(sig_t_R, rho, scale)
        except Exception:
            return 1e6

    # Coarse scan first
    scales = np.linspace(scale_range[0], scale_range[1], 20)
    errs = [f(s) for s in scales]
    idx = np.argmin(errs)
    s_best = scales[idx]
    # Bracket around best
    lo = max(scale_range[0], s_best - 0.3)
    hi = min(scale_range[1], s_best + 0.3)
    res = minimize_scalar(f, bracket=(lo, s_best, hi), method='brent',
                           options={'xtol': tol})
    return res.x, res.fun


def main():
    print("=" * 78)
    print("E3.6 — scale scan of constant inner basis")
    print("=" * 78)

    sig_t_Rs = [1.0, 2.5, 5.0, 10.0, 20.0, 50.0, 100.0]
    rhos = [0.1, 0.3, 0.5, 0.7]

    print(f"\n{'σ_t·R':>8} {'ρ':>6} {'F.4':>10} {'scale=√2':>10} {'scale=√3':>10} "
          f"{'scale_opt':>10} {'err_opt':>10}")
    print("-" * 80)

    opts = {}
    for sig_t_R in sig_t_Rs:
        for rho in rhos:
            try:
                R = sig_t_R
                r_0 = rho * R
                k_f4 = run_scalar_f4(r_0, R, 1.0, 1.0 / 3.0, 1.0)
                err_f4 = abs(k_f4 - K_INF) / K_INF * 100
            except Exception:
                err_f4 = float("nan")

            try:
                err_leg = run_scale(sig_t_R, rho, math.sqrt(2.0))
            except Exception:
                err_leg = float("nan")
            try:
                err_jac = run_scale(sig_t_R, rho, math.sqrt(3.0))
            except Exception:
                err_jac = float("nan")
            try:
                s_opt, e_opt = find_optimal_scale(sig_t_R, rho)
            except Exception:
                s_opt, e_opt = float("nan"), float("nan")

            opts[(sig_t_R, rho)] = (s_opt, e_opt, err_f4)

            print(f"{sig_t_R:>8.1f} {rho:>6.2f} {err_f4:>10.4f} "
                   f"{err_leg:>10.4f} {err_jac:>10.4f} "
                   f"{s_opt:>10.4f} {e_opt:>10.4f}")

    print("\n" + "=" * 78)
    print("Pattern analysis of scale_opt:")
    print("=" * 78)
    print(f"{'σ_t·R':>8} {'ρ':>6} {'scale_opt':>10} "
          f"{'τ(1-ρ)':>8} {'√(3·(?))':>10}")
    for (sig_t_R, rho), (s_opt, e_opt, e_f4) in opts.items():
        tau = sig_t_R
        print(f"{sig_t_R:>8.1f} {rho:>6.2f} {s_opt:>10.4f} "
               f"{tau*(1-rho):>8.2f} {3.0 * (1 - rho):>10.4f}")

    # Is there a principled formula? Let's check:
    #   scale²_opt(τ,ρ) - ???
    print("\n  Check formula: scale² = 2 + C(τ,ρ)")
    for (sig_t_R, rho), (s_opt, e_opt, e_f4) in opts.items():
        if not math.isnan(s_opt):
            excess = s_opt ** 2 - 2.0  # over Legendre baseline
            print(f"    σ_t·R={sig_t_R}, ρ={rho}: scale²-2 = {excess:.4f}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
