"""Diagnostic E4: Regime-switched closure (F.4 at thin τ, split+scale at thick τ).

Created by numerics-investigator on 2026-04-22.
If this test catches a real bug, promote to ``tests/cp/`` — this
diagnostic verifies rank-N universal closures over wide (σ_t·R, ρ) range.

HYPOTHESIS (Direction M from research log):
- L10 established: split basis is catastrophic at σ_t·R ≤ 2.5 (3-86%
  err for ALL inner bases).
- E3.1/E3.7 established: split basis with calibrated scale BEATS F.4
  at σ_t·R ≥ 5 by 2-20× for moderate ρ.

Goal: ship a REGIME-SWITCHED closure that:
- Uses F.4 for τ ≤ THRESHOLD_LOW (thin τ) — safe, bounded.
- Uses split-basis rank-(1,1,1) with calibrated scale for τ ≥
  THRESHOLD_HIGH — wins by 2-20× at thick τ.
- Blends in between.

Two scale-calibration modes:
- (a) formula: scale²_opt = (1+6ρ)/(3ρ)    [no τ dep, simple]
- (b) brent:   1D optimize scale at each problem (expensive, optimal)

For a shippable universal closure:
- Regime-switched(a) should beat F.4 (at matched rich quad) at EVERY (τ, ρ).
- If not, regime-switched(b) should beat F.4 at EVERY (τ, ρ).
- Blend region must be SMOOTH — no dip/bump in err(τ) at threshold.

Runtime: ~15-30 min for full scan.
"""
from __future__ import annotations

import math
import sys
import time

import numpy as np
from scipy.optimize import minimize_scalar

sys.path.insert(0, '/workspaces/ORPHEUS/derivations/diagnostics')

from diag_cin_aware_split_basis_keff import run_scalar_f4
from diag_cin_split_source_decomposition import make_constant_basis, run_custom_basis

K_INF = 1.5

# Quadrature controls — MED is fast enough for 32-point scan but faithful
# enough that F.4 floors out. Default for scan is MED. RICH available for
# targeted re-runs. BASE = baseline.
# We export RICH as alias of MED so the main() code can use a single
# variable name; swap to BASE or true-RICH by overriding.
BASE = dict(n_panels=2, p_order=4, n_ang=32)
MED = dict(n_panels=3, p_order=6, n_ang=48)
TRUE_RICH = dict(n_panels=4, p_order=8, n_ang=64)
RICH = MED  # <-- what main() uses; change here to TRUE_RICH for final run

# Regime-switch thresholds
THRESHOLD_LOW = 3.0   # σ_t·R below this → pure F.4
THRESHOLD_HIGH = 5.0  # σ_t·R above this → pure split+scale


# ---------------------------------------------------------------------------
# Core runners — matched quadrature everywhere
# ---------------------------------------------------------------------------


def run_f4_keff(sig_t_R, rho, sig_t=1.0, **quad):
    R = sig_t_R / sig_t
    r_0 = rho * R
    sig_s = sig_t / 3.0
    nsf = sig_t
    return run_scalar_f4(r_0, R, sig_t, sig_s, nsf, **quad)


def run_split_scale_keff(sig_t_R, rho, scale, sig_t=1.0, **quad):
    R = sig_t_R / sig_t
    r_0 = rho * R
    sig_s = sig_t / 3.0
    nsf = sig_t
    basis = make_constant_basis(scale)
    return run_custom_basis(
        r_0, R, sig_t, sig_s, nsf, basis,
        coupling_weight_power=1, **quad,
    )


def scale_from_formula(rho):
    """scale²_opt = (1 + 6ρ) / (3ρ)  — E3.7 conjecture."""
    return math.sqrt((1.0 + 6.0 * rho) / (3.0 * rho))


def scale_brent(sig_t_R, rho, quad, bounds=(1.0, 2.8)):
    """Bounded 1D minimization on scale (avoids divergent Brent expansion).

    Uses bounded method (safer than Brent which extends brackets)."""
    def err(scale):
        try:
            k = run_split_scale_keff(sig_t_R, rho, scale, **quad)
            v = abs(k - K_INF) / K_INF
            if not np.isfinite(v):
                return 1e6
            return v
        except Exception:
            return 1e6
    res = minimize_scalar(
        err, bounds=bounds, method='bounded',
        options={'xatol': 1e-5, 'maxiter': 25},
    )
    return res.x, res.fun


# ---------------------------------------------------------------------------
# Regime-switched closure variants
# ---------------------------------------------------------------------------


def regime_switched_keff(sig_t_R, rho, *, scale_mode='formula', quad=RICH,
                          threshold_low=THRESHOLD_LOW,
                          threshold_high=THRESHOLD_HIGH):
    """Build regime-switched k_eff.

    scale_mode:
      - 'formula' → scale = √((1+6ρ)/(3ρ)) (cheap, no τ dep)
      - 'brent'   → scale = 1D Brent optimum (expensive, optimal)

    Blending scheme: linear interpolation of k_eff between the two
    closures over [threshold_low, threshold_high].
    """
    if sig_t_R <= threshold_low:
        return run_f4_keff(sig_t_R, rho, **quad), 'f4'
    k_f4 = run_f4_keff(sig_t_R, rho, **quad)
    # Get scale
    if scale_mode == 'formula':
        scale = scale_from_formula(rho)
    elif scale_mode == 'brent':
        scale, _ = scale_brent(sig_t_R, rho, quad)
    else:
        raise ValueError(scale_mode)
    k_split = run_split_scale_keff(sig_t_R, rho, scale, **quad)
    if sig_t_R >= threshold_high:
        return k_split, 'split'
    # Blend linearly in σ_t·R
    w = (sig_t_R - threshold_low) / (threshold_high - threshold_low)
    k_blend = (1.0 - w) * k_f4 + w * k_split
    return k_blend, f'blend(w={w:.2f})'


# ---------------------------------------------------------------------------
# Main scan
# ---------------------------------------------------------------------------


def main(quad=None):
    if quad is None:
        quad = MED
    print("=" * 108)
    print("E4 — Regime-switched closure (F.4 at thin τ, split+scale at thick τ)")
    print(f"Quadrature: {quad}")
    print(f"Thresholds: LOW = {THRESHOLD_LOW}, HIGH = {THRESHOLD_HIGH}")
    print("=" * 108)
    globals()['RICH'] = quad   # so regime_switched uses this in nested calls

    sig_t_Rs = [0.5, 1.0, 2.5, 5.0, 10.0, 20.0, 50.0, 100.0]
    rhos = [0.1, 0.3, 0.5, 0.7]

    # Header
    print(f"\n{'σ_t·R':>8} {'ρ':>5} "
          f"{'F.4':>10} "
          f"{'scale_form':>11} {'err_form':>10} "
          f"{'scale_brent':>12} {'err_brent':>10} "
          f"{'RS_form':>10} {'RS_brent':>10}")
    print("-" * 108)

    # Collect rows for post-summary
    rows = []
    t0 = time.time()
    for sig_t_R in sig_t_Rs:
        for rho in rhos:
            # F.4 baseline
            try:
                k_f4 = run_f4_keff(sig_t_R, rho, **RICH)
                err_f4 = abs(k_f4 - K_INF) / K_INF * 100
            except Exception as e:
                k_f4 = float('nan'); err_f4 = float('nan')

            # Split+formula scale
            sc_form = scale_from_formula(rho)
            try:
                k_form = run_split_scale_keff(sig_t_R, rho, sc_form, **RICH)
                err_form = abs(k_form - K_INF) / K_INF * 100
            except Exception:
                k_form = float('nan'); err_form = float('nan')

            # Split+Brent scale (skip at thin — it'll blow up)
            if sig_t_R >= THRESHOLD_LOW:
                try:
                    sc_br, err_br_frac = scale_brent(sig_t_R, rho, RICH)
                    err_br = err_br_frac * 100
                except Exception:
                    sc_br = float('nan'); err_br = float('nan')
            else:
                sc_br = float('nan'); err_br = float('nan')

            # Regime-switched: formula-scale variant
            try:
                k_rs_f, _ = regime_switched_keff(sig_t_R, rho, scale_mode='formula', quad=RICH)
                err_rs_f = abs(k_rs_f - K_INF) / K_INF * 100
            except Exception:
                err_rs_f = float('nan')

            # Regime-switched: brent-scale variant
            try:
                k_rs_b, _ = regime_switched_keff(sig_t_R, rho, scale_mode='brent', quad=RICH)
                err_rs_b = abs(k_rs_b - K_INF) / K_INF * 100
            except Exception:
                err_rs_b = float('nan')

            print(f"{sig_t_R:>8.2f} {rho:>5.2f} "
                  f"{err_f4:>9.4f}% "
                  f"{sc_form:>11.4f} {err_form:>9.4f}% "
                  f"{sc_br:>12.4f} {err_br:>9.4f}% "
                  f"{err_rs_f:>9.4f}% {err_rs_b:>9.4f}%")

            rows.append({
                'sig_t_R': sig_t_R, 'rho': rho,
                'err_f4': err_f4, 'err_form': err_form, 'err_br': err_br,
                'err_rs_f': err_rs_f, 'err_rs_b': err_rs_b,
                'scale_form': sc_form, 'scale_br': sc_br,
            })

    dt = time.time() - t0
    print(f"\nTotal runtime: {dt:.1f}s")

    # ------------------------------------------------------------------
    # Summary: does regime-switched beat F.4 universally?
    # ------------------------------------------------------------------
    def _count_wins(variant_key):
        wins = 0
        comparable = 0
        max_ratio_lose = 0.0
        max_ratio_win = 0.0
        worst_lose = None
        for r in rows:
            if math.isnan(r['err_f4']) or math.isnan(r[variant_key]):
                continue
            comparable += 1
            if r[variant_key] < r['err_f4']:
                wins += 1
                ratio = r['err_f4'] / max(r[variant_key], 1e-9)
                if ratio > max_ratio_win:
                    max_ratio_win = ratio
            else:
                ratio = r[variant_key] / max(r['err_f4'], 1e-9)
                if ratio > max_ratio_lose:
                    max_ratio_lose = ratio
                    worst_lose = r
        return wins, comparable, max_ratio_win, max_ratio_lose, worst_lose

    print("\n" + "=" * 78)
    print("Summary — universal-closure test")
    print("=" * 78)

    for key, label in [('err_rs_f', 'Regime-switched(formula)'),
                        ('err_rs_b', 'Regime-switched(brent)')]:
        wins, tot, rwin, rlose, worst = _count_wins(key)
        universal = "UNIVERSAL ✓" if wins == tot else f"not universal ({tot - wins} losses)"
        print(f"\n{label}: {wins}/{tot} wins — {universal}")
        print(f"  best win ratio: {rwin:.1f}×")
        print(f"  worst lose ratio: {rlose:.1f}×"
              + (f"  at σ_t·R={worst['sig_t_R']}, ρ={worst['rho']}"
                 if worst else ''))

    # Monotonicity check at the transition
    print("\n" + "=" * 78)
    print("Transition smoothness check at σ_t·R ∈ [2.5, 3.5, 4.0, 4.5, 5.0, 6.0] for ρ=0.3")
    print("=" * 78)
    for stR in [2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0]:
        try:
            k_rs, regime = regime_switched_keff(stR, 0.3, scale_mode='formula', quad=RICH)
            err = abs(k_rs - K_INF) / K_INF * 100
            print(f"  σ_t·R={stR:>5.2f}: err={err:.4f}% [{regime}]")
        except Exception as e:
            print(f"  σ_t·R={stR:>5.2f}: FAIL {e}")

    return 0


# ---------------------------------------------------------------------------
# Pytest harness (subset, runs in <5 min at reduced scan)
# ---------------------------------------------------------------------------


import pytest


@pytest.mark.slow
def test_regime_switched_beats_f4_thick_core():
    """At thick τ and moderate ρ, regime-switched (brent) must beat F.4."""
    sig_t_R, rho = 10.0, 0.3
    k_f4 = run_f4_keff(sig_t_R, rho, **RICH)
    err_f4 = abs(k_f4 - K_INF) / K_INF
    k_rs, _ = regime_switched_keff(sig_t_R, rho, scale_mode='brent', quad=RICH)
    err_rs = abs(k_rs - K_INF) / K_INF
    assert err_rs < err_f4, (
        f"regime-switched (brent) did NOT beat F.4 at σ_t·R={sig_t_R}, ρ={rho}: "
        f"err_rs={err_rs*100:.4f}% vs err_f4={err_f4*100:.4f}%"
    )


@pytest.mark.slow
def test_regime_switched_safe_thin():
    """At thin τ, regime-switched must fall back to F.4 (same k_eff)."""
    sig_t_R, rho = 1.0, 0.3
    k_f4 = run_f4_keff(sig_t_R, rho, **RICH)
    k_rs, regime = regime_switched_keff(sig_t_R, rho, scale_mode='formula', quad=RICH)
    assert regime == 'f4'
    assert abs(k_rs - k_f4) < 1e-10, f"k_rs={k_rs} != k_f4={k_f4}"


if __name__ == "__main__":
    sys.exit(main())
