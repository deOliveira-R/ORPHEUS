"""Test the candidate formula scale²_opt = (1 + 6ρ) / (3ρ) for constant inner basis.

Derived from E3.1 α-scan observation that α_opt · ρ ≈ 0.3 across σ_t·R ≥ 5.
Converting Jacobi weight α back to constant scale: scale² = α + 2.
So scale²_opt = 1/(3ρ) + 2 = (1 + 6ρ) / (3ρ).

TEST: at each (σ_t·R, ρ) point, compare:
  1. err at scale = formula-predicted (1 + 6ρ)/(3ρ)^0.5
  2. err at scale = √2 (Legendre baseline)
  3. err at scale = √3 (Jacobi-c²)
  4. err at F.4
  5. err at coarsely-scanned empirical optimum

Point count ~12 at σ_t·R ≥ 5. Runtime ~2 minutes.
"""
from __future__ import annotations

import math
import sys

import numpy as np

sys.path.insert(0, '/workspaces/ORPHEUS/derivations/diagnostics')

from diag_cin_aware_split_basis_keff import run_scalar_f4
from diag_cin_split_source_decomposition import make_constant_basis, run_custom_basis

K_INF = 1.5


def err_for_scale(sig_t_R, rho, scale, sig_t=1.0):
    R = sig_t_R / sig_t
    r_0 = rho * R
    basis = make_constant_basis(scale)
    k = run_custom_basis(r_0, R, sig_t, sig_t / 3.0, sig_t, basis)
    return abs(k - K_INF) / K_INF * 100


def empirical_best_scale(sig_t_R, rho, scales):
    best = (float('nan'), float('inf'))
    for s in scales:
        try:
            e = err_for_scale(sig_t_R, rho, s)
            if e < best[1]:
                best = (s, e)
        except Exception:
            pass
    return best


def main():
    print("=" * 90)
    print("Testing scale formula: scale²_opt = (1 + 6ρ) / (3ρ) for constant inner basis")
    print("=" * 90)

    points = [
        (5.0,  0.1), (5.0,  0.3), (5.0,  0.5), (5.0,  0.7),
        (10.0, 0.1), (10.0, 0.3), (10.0, 0.5), (10.0, 0.7),
        (20.0, 0.1), (20.0, 0.3), (20.0, 0.5), (20.0, 0.7),
    ]

    print(f"\n{'σ_t·R':>8} {'ρ':>6} {'sc_form':>8} {'sc_emp':>8} {'err_form':>10} "
          f"{'err_emp':>10} {'err_F.4':>10} {'err_√2':>10} {'err_√3':>10}")
    print("-" * 90)

    # Search grid for empirical optimum
    search_scales = [1.2, 1.4, 1.5, 1.6, 1.7, 1.732, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.5]

    rows = []
    for sig_t_R, rho in points:
        scale_formula = math.sqrt((1 + 6*rho) / (3*rho))
        try:
            err_form = err_for_scale(sig_t_R, rho, scale_formula)
        except Exception:
            err_form = float('nan')
        scale_emp, err_emp = empirical_best_scale(sig_t_R, rho, search_scales)
        try:
            k_f4 = run_scalar_f4(rho * sig_t_R, sig_t_R, 1.0, 1.0/3.0, 1.0)
            err_f4 = abs(k_f4 - K_INF) / K_INF * 100
        except Exception:
            err_f4 = float('nan')
        try:
            err_sqrt2 = err_for_scale(sig_t_R, rho, math.sqrt(2.0))
        except Exception:
            err_sqrt2 = float('nan')
        try:
            err_sqrt3 = err_for_scale(sig_t_R, rho, math.sqrt(3.0))
        except Exception:
            err_sqrt3 = float('nan')
        print(f"{sig_t_R:>8.1f} {rho:>6.2f} {scale_formula:>8.4f} {scale_emp:>8.4f} "
              f"{err_form:>10.4f}% {err_emp:>10.4f}% {err_f4:>10.4f}% "
              f"{err_sqrt2:>10.4f}% {err_sqrt3:>10.4f}%")
        rows.append((sig_t_R, rho, scale_formula, err_form,
                     scale_emp, err_emp, err_f4, err_sqrt2, err_sqrt3))

    # Summary
    print("\n" + "=" * 90)
    print("Wins per point:")
    print("=" * 90)
    print(f"{'σ_t·R':>8} {'ρ':>6} {'formula beats F.4?':>20} {'emp beats F.4?':>18}")
    for r in rows:
        s, rho, sf, ef, se, ee, ef4, e2, e3 = r
        fw = "YES" if ef < ef4 else "no"
        ew = "YES" if ee < ef4 else "no"
        print(f"{s:>8.1f} {rho:>6.2f} {fw:>20} {ew:>18}")


if __name__ == "__main__":
    main()
