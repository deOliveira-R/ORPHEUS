"""Precision check: how low does rank-(1,1,1) with Brent-optimized scale actually go?

E5 showed 0.0000% at display precision. We want to see the TRUE magnitude
of the residual at the optimal scale. Is it really machine-precision, or
just quadrature-floor (~1e-4%)?
"""
from __future__ import annotations

import math
import sys

import numpy as np
from scipy.optimize import minimize_scalar

sys.path.insert(0, '/workspaces/ORPHEUS/derivations/diagnostics')

from diag_cin_aware_split_basis_keff import run_scalar_f4
from diag_cin_split_source_decomposition import make_constant_basis, run_custom_basis

K_INF = 1.5


def err_for_scale(sig_t_R, rho, scale, n_panels=3, p_order=6, n_ang=48):
    R = sig_t_R / 1.0
    r_0 = rho * R
    basis = make_constant_basis(scale)
    k = run_custom_basis(r_0, R, 1.0, 1.0/3.0, 1.0, basis,
                         n_panels=n_panels, p_order=p_order, n_ang=n_ang)
    return abs(k - K_INF) / K_INF * 100, k


def main():
    print("=" * 80)
    print("Precision check — is rank-(1,1,1) with optimal scale really ~0 err?")
    print("=" * 80)

    test_points = [(5.0, 0.3), (10.0, 0.3), (20.0, 0.3),
                    (5.0, 0.5), (10.0, 0.5), (20.0, 0.5),
                    (5.0, 0.7), (10.0, 0.7)]

    print(f"\n{'σ_t·R':>6} {'ρ':>6} {'scale_opt':>10} {'err_ppm':>14} {'err_ppb':>14}")
    print("-" * 80)

    for sig_t_R, rho in test_points:
        def f(s):
            try:
                return err_for_scale(sig_t_R, rho, s)[0]
            except Exception:
                return 1e6

        # Fine search
        s_grid = np.linspace(1.3, 2.5, 25)
        errs = [f(s) for s in s_grid]
        idx = int(np.argmin(errs))
        lo = s_grid[max(0, idx-1)]
        hi = s_grid[min(len(s_grid)-1, idx+1)]

        try:
            res = minimize_scalar(f, bracket=(lo, s_grid[idx], hi),
                                    method='brent',
                                    options={'xtol': 1e-6, 'maxiter': 30})
            s_opt = res.x
            err_opt = res.fun
        except Exception as e:
            s_opt, err_opt = s_grid[idx], errs[idx]

        err_ppm = err_opt * 1e4  # percent → ppm is × 1e4
        err_ppb = err_opt * 1e7  # percent → ppb

        print(f"{sig_t_R:>6.1f} {rho:>6.2f} {s_opt:>10.6f} "
              f"{err_ppm:>14.4f} ppm {err_ppb:>14.4f} ppb")


if __name__ == "__main__":
    main()
