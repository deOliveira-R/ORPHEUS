"""Diagnostic E5: rank-(1,1,2) with adaptive (α_0, α_1) scales per mode.

Created by numerics-investigator on 2026-04-22.
If this test catches a real bug, promote to ``tests/cp/`` — this
diagnostic verifies whether the second inner mode helps over rank-(1,1,1)
scale-only calibration.

HYPOTHESIS (from L8):
The scale gauge DOF applies to EACH inner basis function independently.
At rank-(1,1,2) we have TWO inner modes: φ_0 (scale α_0) and φ_1
(scale α_1). Both are tunable. If 2D tuning gives significantly less err
than rank-(1,1,1) 1D tuning, rank-(1,1,2) adaptive scales are the path
to sub-0.01% universal accuracy.

HOWEVER — existing E2 data (rank-(1,1,N) w/ N ≥ 2 plateaus at 0.99%
with Legendre, 0.073% with Jacobi-c²) suggests mode-1 information is
already captured by rank-(1,1,2); the plateau is structural in the
Petrov-Galerkin weighting of both modes jointly.

The question E5 settles: does DOUBLE scale calibration (α_0, α_1) give
an extra win beyond the rank-(1,1,1) SINGLE scale optimum?

Mode shapes are fixed: φ_0 = constant, φ_1 = 2c - 1 (shifted Legendre
on [0,1]). Scale each: φ_0 → α_0; φ_1 → α_1 · (2c-1). Both in the c-weight
span — basis is not orthonormal after scaling, but closure doesn't require
orthonormality (Petrov-Galerkin with custom W, P, G).

Runtime: ~15-30 min at RICH quadrature.
"""
from __future__ import annotations

import math
import sys
import time

import numpy as np
from scipy.optimize import minimize, minimize_scalar

sys.path.insert(0, '/workspaces/ORPHEUS/derivations/diagnostics')

from diag_cin_aware_split_basis_keff import run_scalar_f4
from diag_cin_split_source_decomposition import make_constant_basis, run_custom_basis

K_INF = 1.5

BASE = dict(n_panels=2, p_order=4, n_ang=32)
MED = dict(n_panels=3, p_order=6, n_ang=48)
TRUE_RICH = dict(n_panels=4, p_order=8, n_ang=64)
RICH = BASE  # <-- E5 uses BASE for 2D Nelder-Mead speed; re-verify winners at MED


# ---------------------------------------------------------------------------
# Basis construction
# ---------------------------------------------------------------------------


def make_two_mode_basis(alpha_0, alpha_1):
    """Rank-2 inner basis: φ_0 = α_0 (constant), φ_1 = α_1 * (2c-1).

    Half-range Legendre mode-1 on [0,1] is P_1(2c-1) = 2c - 1.
    """
    def phi0(c):
        if np.isscalar(c):
            return float(alpha_0)
        return np.full_like(c, float(alpha_0), dtype=float)

    def phi1(c, a1=float(alpha_1)):
        if np.isscalar(c):
            return a1 * (2.0 * c - 1.0)
        return a1 * (2.0 * np.asarray(c, dtype=float) - 1.0)

    return [phi0, phi1]


def make_two_mode_basis_cpow(alpha_0, alpha_1, beta_1):
    """Alternative: φ_0 = α_0, φ_1 = α_1 · c^β_1 (freer shape)."""
    def phi0(c):
        if np.isscalar(c):
            return float(alpha_0)
        return np.full_like(c, float(alpha_0), dtype=float)

    def phi1(c, a1=float(alpha_1), b1=float(beta_1)):
        if np.isscalar(c):
            return a1 * (c ** b1)
        return a1 * (np.asarray(c, dtype=float) ** b1)

    return [phi0, phi1]


# ---------------------------------------------------------------------------
# Evaluation
# ---------------------------------------------------------------------------


def run_keff_rank112(sig_t_R, rho, alpha_0, alpha_1, sig_t=1.0, quad=RICH):
    R = sig_t_R / sig_t
    r_0 = rho * R
    sig_s = sig_t / 3.0
    nsf = sig_t
    basis = make_two_mode_basis(alpha_0, alpha_1)
    k = run_custom_basis(r_0, R, sig_t, sig_s, nsf, basis, **quad)
    return k


def err_rank112(sig_t_R, rho, alpha_0, alpha_1, **kwargs):
    try:
        k = run_keff_rank112(sig_t_R, rho, alpha_0, alpha_1, **kwargs)
        return abs(k - K_INF) / K_INF
    except Exception:
        return 1e6


# ---------------------------------------------------------------------------
# rank-(1,1,1) reference (scale-only)
# ---------------------------------------------------------------------------


def run_keff_rank111(sig_t_R, rho, scale, sig_t=1.0, quad=RICH):
    R = sig_t_R / sig_t
    r_0 = rho * R
    sig_s = sig_t / 3.0
    nsf = sig_t
    basis = make_constant_basis(scale)
    k = run_custom_basis(r_0, R, sig_t, sig_s, nsf, basis, **quad)
    return k


def brent_scale_rank111(sig_t_R, rho, quad=RICH):
    def err(s):
        try:
            k = run_keff_rank111(sig_t_R, rho, s, quad=quad)
            return abs(k - K_INF) / K_INF
        except Exception:
            return 1e6
    res = minimize_scalar(err, bracket=(1.0, 1.7, 2.5), method='brent',
                           options={'xtol': 1e-5, 'maxiter': 30})
    return res.x, res.fun


def optimize_rank112_2d(sig_t_R, rho, quad=RICH, x0=None):
    """Nelder-Mead 2D optimum of (α_0, α_1)."""
    if x0 is None:
        x0 = [math.sqrt(3.0), 1.0]

    def err(x):
        a0, a1 = x
        return err_rank112(sig_t_R, rho, a0, a1, quad=quad)

    res = minimize(err, x0, method='Nelder-Mead',
                   options={'xatol': 3e-4, 'fatol': 1e-6, 'maxiter': 80})
    return tuple(res.x), res.fun


# ---------------------------------------------------------------------------
# Scans
# ---------------------------------------------------------------------------


def scan_key_points(quad=RICH):
    """Main scan — compare rank-(1,1,1) 1D optim vs rank-(1,1,2) 2D optim."""
    points = [
        (5.0,  0.1), (5.0,  0.3), (5.0,  0.5), (5.0,  0.7),
        (10.0, 0.1), (10.0, 0.3), (10.0, 0.5), (10.0, 0.7),
        (20.0, 0.1), (20.0, 0.3), (20.0, 0.5), (20.0, 0.7),
        (50.0, 0.3), (50.0, 0.5),
    ]

    print(f"\n{'σ_t·R':>7} {'ρ':>5} "
          f"{'F.4':>9} "
          f"{'sc111':>7} {'err111':>9} "
          f"{'a0':>7} {'a1':>7} {'err112':>9} "
          f"{'ratio':>7}")
    print("-" * 88)

    rows = []
    for sig_t_R, rho in points:
        # F.4
        try:
            k_f4 = run_scalar_f4(rho * sig_t_R, sig_t_R, 1.0, 1.0/3.0, 1.0, **quad)
            err_f4 = abs(k_f4 - K_INF) / K_INF * 100
        except Exception:
            err_f4 = float('nan')

        # rank-(1,1,1) Brent scale
        try:
            sc, err_111_frac = brent_scale_rank111(sig_t_R, rho, quad=quad)
            err_111 = err_111_frac * 100
        except Exception:
            sc = float('nan'); err_111 = float('nan')

        # rank-(1,1,2) 2D optimum — seed at (sc, 1.0) to give it a fair start
        x0 = [sc if not math.isnan(sc) else math.sqrt(3.0), 1.0]
        try:
            (a0_opt, a1_opt), err_112_frac = optimize_rank112_2d(
                sig_t_R, rho, quad=quad, x0=x0,
            )
            err_112 = err_112_frac * 100
        except Exception:
            a0_opt = a1_opt = float('nan'); err_112 = float('nan')

        # Ratio: how much does 2D help over 1D?
        if err_111 > 0 and not math.isnan(err_112):
            ratio = err_111 / max(err_112, 1e-9)
        else:
            ratio = float('nan')

        print(f"{sig_t_R:>7.1f} {rho:>5.2f} "
              f"{err_f4:>8.4f}% "
              f"{sc:>7.4f} {err_111:>8.4f}% "
              f"{a0_opt:>7.4f} {a1_opt:>7.4f} {err_112:>8.4f}% "
              f"{ratio:>7.2f}×")

        rows.append({
            'sig_t_R': sig_t_R, 'rho': rho,
            'err_f4': err_f4, 'err_111': err_111, 'err_112': err_112,
            'sc': sc, 'a0': a0_opt, 'a1': a1_opt,
        })

    return rows


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main():
    print("=" * 98)
    print("E5 — rank-(1,1,2) adaptive scales vs rank-(1,1,1) scale-only")
    print(f"Quadrature: {RICH}")
    print("=" * 98)

    t0 = time.time()
    rows = scan_key_points(RICH)
    dt = time.time() - t0

    print(f"\nTotal runtime: {dt:.1f}s")

    # Summary
    print("\n" + "=" * 78)
    print("Summary — does rank-(1,1,2) 2D tuning help over rank-(1,1,1) 1D?")
    print("=" * 78)

    improves = 0
    degrades = 0
    neutral = 0
    max_improve = 0.0
    max_degrade = 0.0
    for r in rows:
        if math.isnan(r['err_111']) or math.isnan(r['err_112']):
            continue
        if r['err_112'] < 0.95 * r['err_111']:  # > 5% improvement
            improves += 1
            max_improve = max(max_improve, r['err_111'] / max(r['err_112'], 1e-9))
        elif r['err_112'] > 1.05 * r['err_111']:  # > 5% degradation
            degrades += 1
            max_degrade = max(max_degrade, r['err_112'] / max(r['err_111'], 1e-9))
        else:
            neutral += 1

    print(f"\nrank-(1,1,2) improves (>5%) : {improves}/{len(rows)}")
    print(f"rank-(1,1,2) degrades (>5%) : {degrades}/{len(rows)}")
    print(f"rank-(1,1,2) neutral        : {neutral}/{len(rows)}")
    print(f"Max improvement ratio (1D/2D): {max_improve:.2f}×")
    print(f"Max degradation ratio (2D/1D): {max_degrade:.2f}×")

    # Verdict
    print("\n" + "=" * 78)
    if improves >= 0.7 * len(rows) and max_improve > 2.0:
        print("VERDICT: rank-(1,1,2) 2D adaptive HELPS — shippable path forward.")
    elif improves > degrades:
        print("VERDICT: rank-(1,1,2) 2D adaptive MARGINALLY HELPS — not worth the")
        print("         complexity over rank-(1,1,1) scale-only.")
    else:
        print("VERDICT: rank-(1,1,2) 2D adaptive does NOT help over rank-(1,1,1).")
        print("         The plateau is structural, not a scale DOF issue at mode-1.")
    print("=" * 78)

    return 0


# ---------------------------------------------------------------------------
# Pytest harness
# ---------------------------------------------------------------------------


import pytest


@pytest.mark.slow
def test_rank112_2d_at_least_as_good_as_rank111():
    """2D optim can't beat 1D by a trivial margin → at least tie."""
    sig_t_R, rho = 10.0, 0.3
    _, err_111 = brent_scale_rank111(sig_t_R, rho, quad=RICH)
    _, err_112 = optimize_rank112_2d(sig_t_R, rho, quad=RICH)
    # 2D optim at (sc, 0) should reduce to 1D Brent — so 2D err ≤ 1D err
    # always (up to optimization noise).
    assert err_112 <= 1.05 * err_111, (
        f"rank-(1,1,2) 2D err {err_112*100:.4f}% is much worse than "
        f"rank-(1,1,1) 1D err {err_111*100:.4f}%"
    )


if __name__ == "__main__":
    sys.exit(main())
