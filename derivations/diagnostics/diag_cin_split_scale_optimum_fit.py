"""Focused diagnostic: find scale_opt(τ, ρ) for the constant inner basis and fit.

Built fast. Uses scipy Brent optimization to find scale_opt to 3 digits at
each (σ_t·R, r_0/R) point. Tests simple formulas:
  scale²_opt = 2 + g(τ, ρ)

Hypotheses:
  H1: g = τ·ρ² / (1 + τ·ρ²)  — saturates between 0 (thin) and 1 (thick)
  H2: g = 1/(1 + exp(-τ·ρ)) - 1/2  — sigmoid in τρ
  H3: g = 1 - exp(-τ(1-ρ))  — tracks mean-chord attenuation
  H4: empirical fit (scipy.optimize.curve_fit)

Goal: find a formula that gives scale²_opt accurate to ~1% across the
(σ_t·R, r_0/R) parameter space — if yes, the adaptive closure is
principled and universal.
"""
from __future__ import annotations

import math
import sys

import numpy as np
from scipy.optimize import brent

sys.path.insert(0, '/workspaces/ORPHEUS/derivations/diagnostics')

from diag_cin_aware_split_basis_keff import run_scalar_f4
from diag_cin_split_source_decomposition import make_constant_basis, run_custom_basis

K_INF = 1.5


def err_for_scale(sig_t_R, rho, scale, sig_t=1.0):
    """k_eff error for constant inner-basis with given scale."""
    R = sig_t_R / sig_t
    r_0 = rho * R
    basis = make_constant_basis(scale)
    k = run_custom_basis(r_0, R, sig_t, sig_t / 3.0, sig_t, basis)
    return abs(k - K_INF) / K_INF * 100


def find_scale_opt(sig_t_R, rho, tol=1e-3):
    """Coarse scan + Brent refinement."""
    # Coarse scan over a wide range to locate the basin
    s_grid = np.linspace(0.5, 4.0, 15)
    errs = []
    for s in s_grid:
        try:
            errs.append(err_for_scale(sig_t_R, rho, float(s)))
        except Exception:
            errs.append(1e3)
    idx = int(np.argmin(errs))
    s_best = float(s_grid[idx])
    # Refine around best
    lo = float(s_grid[max(0, idx-1)])
    hi = float(s_grid[min(len(s_grid)-1, idx+1)])

    def f(scale):
        try:
            return err_for_scale(sig_t_R, rho, float(scale))
        except Exception:
            return 1e3

    if lo >= hi:  # edge case
        return s_best, errs[idx]
    try:
        s_opt = brent(f, brack=(lo, s_best, hi), tol=tol, maxiter=20)
        e_opt = f(s_opt)
        return float(s_opt), float(e_opt)
    except Exception:
        return s_best, errs[idx]


def main():
    print("=" * 80)
    print("Scale-optimum fit: scale²_opt(σ_t·R, r_0/R) for constant inner basis")
    print("=" * 80)

    # Strategic grid — focus on σ_t·R ≥ 5 where the closure works.
    # Skip σ_t·R ≤ 2.5 (catastrophic for rank-1 splitbasis, see E3.1).
    points = [
        (5.0,  0.1), (5.0,  0.3), (5.0,  0.5), (5.0,  0.7),
        (10.0, 0.1), (10.0, 0.3), (10.0, 0.5), (10.0, 0.7),
        (20.0, 0.1), (20.0, 0.3), (20.0, 0.5), (20.0, 0.7),
        (50.0, 0.3), (50.0, 0.5), (50.0, 0.7),
    ]

    rows = []
    print(f"\n{'σ_t·R':>8} {'ρ':>6} {'s_opt':>10} {'s²_opt':>10} {'s²-2':>8} {'err_opt':>12} {'F.4 err':>12}")
    print("-" * 80)
    for sig_t_R, rho in points:
        s_opt, e_opt = find_scale_opt(sig_t_R, rho)
        try:
            k_f4 = run_scalar_f4(rho * sig_t_R, sig_t_R, 1.0, 1.0 / 3.0, 1.0)
            e_f4 = abs(k_f4 - K_INF) / K_INF * 100
        except Exception:
            e_f4 = float("nan")
        s_sq = s_opt ** 2
        excess = s_sq - 2.0
        print(f"{sig_t_R:>8.1f} {rho:>6.2f} {s_opt:>10.4f} {s_sq:>10.4f} "
              f"{excess:>8.4f} {e_opt:>12.4f}% {e_f4:>12.4f}%")
        rows.append((sig_t_R, rho, s_opt, s_sq, excess, e_opt, e_f4))

    print("\n" + "=" * 80)
    print("Hypothesis fits (g := s²_opt - 2):")
    print("=" * 80)
    print(f"{'σ_t·R':>8} {'ρ':>6} {'g_meas':>10} "
          f"{'τρ²':>10} {'τρ/(1+τρ)':>12} {'1-exp(-τ(1-ρ))':>16}")
    print("-" * 80)
    for sig_t_R, rho, s_opt, s_sq, excess, e_opt, e_f4 in rows:
        h1 = sig_t_R * rho**2 / (1 + sig_t_R * rho**2)
        h2 = sig_t_R * rho / (1 + sig_t_R * rho)
        h3 = 1.0 - math.exp(-sig_t_R * (1 - rho))
        print(f"{sig_t_R:>8.1f} {rho:>6.2f} {excess:>10.4f} "
              f"{h1:>10.4f} {h2:>12.4f} {h3:>16.4f}")

    # Quick least-squares fit: g(τ, ρ) = a·τρ + b·ρ² + c·τρ²/(1+τρ²) + d
    print("\n" + "=" * 80)
    print("Least-squares fit of g(τ, ρ):")
    print("=" * 80)
    X = []
    y = []
    for sig_t_R, rho, s_opt, s_sq, excess, e_opt, e_f4 in rows:
        if math.isnan(excess):
            continue
        tau = sig_t_R
        X.append([tau * rho, rho**2, tau * rho**2 / (1 + tau * rho**2), 1.0])
        y.append(excess)
    X = np.array(X)
    y = np.array(y)
    coeffs, *_ = np.linalg.lstsq(X, y, rcond=None)
    print(f"  g(τ, ρ) ≈ {coeffs[0]:+.4f}·τρ {coeffs[1]:+.4f}·ρ² "
          f"{coeffs[2]:+.4f}·τρ²/(1+τρ²) {coeffs[3]:+.4f}")
    # RMSE
    y_pred = X @ coeffs
    rmse = math.sqrt(np.mean((y - y_pred) ** 2))
    print(f"  RMSE: {rmse:.4f}   (rel to mean g: {rmse/abs(np.mean(y)):.1%})")


if __name__ == "__main__":
    main()
