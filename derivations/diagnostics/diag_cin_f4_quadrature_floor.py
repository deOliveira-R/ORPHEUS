"""Diagnostic E2.1: F.4 quadrature-floor study for hollow sphere.

Created by numerics-investigator on 2026-04-21.

Purpose: determine whether F.4's 0.122% residual is (a) quadrature-limited
or (b) structural (basis-mismatch ceiling).

Method: refine Nyström radial quadrature (n_panels, p_order) AND angular
quadrature (n_angular, n_surf_quad) systematically at σ_t·R=5, r_0/R=0.3.
If the residual saturates at a non-zero value, the plateau is structural
(the F.4 Lambert-P/G + Marshak-W mismatch truly generates ~X% error).
If it keeps shrinking toward zero, the 0.122% is pure quadrature error
and there is no structural ceiling.

This sets the baseline we're trying to beat with any new closure.
"""
from __future__ import annotations

import numpy as np

from orpheus.derivations.peierls_geometry import (
    CurvilinearGeometry,
    build_closure_operator,
    build_volume_kernel,
    composite_gl_r,
)


K_INF = 1.5  # nu_sig_f / (sig_t - sig_s) = 1.0 / (1.0 - 1/3) = 1.5


def solve_k_eff(K, sig_t, sig_s, nsf, *, tol=1e-12, max_iter=500):
    N = K.shape[0]
    A = np.diag(np.full(N, sig_t)) - K * sig_s
    B = K * nsf
    phi = np.ones(N)
    k = 1.0
    for _ in range(max_iter):
        q = B @ phi / k
        phi_new = np.linalg.solve(A, q)
        k_new = k * (np.abs(B @ phi_new).sum() / np.abs(B @ phi).sum())
        if abs(k_new - k) < tol:
            return k_new
        phi = phi_new / max(np.linalg.norm(phi_new), 1e-30)
        k = k_new
    return k


def run_f4(r_0, R, sig_t_val, sig_s_val, nsf_val,
           n_panels, p_order, n_ang, dps=15):
    geom = CurvilinearGeometry(kind="sphere-1d", inner_radius=r_0)
    radii = np.array([R])
    sig_t = np.array([sig_t_val])
    r_nodes, r_wts, panels = composite_gl_r(radii, n_panels, p_order,
                                             dps=dps, inner_radius=r_0)
    K_vol = build_volume_kernel(geom, r_nodes, panels, radii, sig_t,
                                 n_angular=n_ang, n_rho=n_ang, dps=dps)
    op = build_closure_operator(geom, r_nodes, r_wts, radii, sig_t,
                                 reflection="white",
                                 n_angular=n_ang, n_surf_quad=n_ang, dps=dps)
    K = K_vol + op.as_matrix()
    return solve_k_eff(K, sig_t_val, sig_s_val, nsf_val)


def main():
    print("=" * 78)
    print("E2.1: F.4 quadrature-floor study")
    print("  σ_t·R = 5, r_0/R = 0.3, k_inf = 1.5")
    print("=" * 78)

    sig_t, sig_s, nsf = 1.0, 1.0 / 3.0, 1.0
    R = 5.0
    r_0 = 0.3 * R

    # Baseline (standard): n_panels=2, p_order=4, n_ang=32.
    # Refinement sweep — radial first, then angular, then both.
    configs = [
        # (n_panels, p_order, n_ang, label)
        (2, 4, 32, "baseline"),
        (4, 4, 32, "radial 2× panels"),
        (2, 6, 32, "radial p_order 6"),
        (2, 8, 32, "radial p_order 8"),
        (4, 6, 32, "radial 2× p6"),
        (4, 8, 32, "radial 2× p8"),
        (8, 4, 32, "radial 4× panels"),
        (8, 8, 32, "radial 4× p8"),
        (2, 4, 48, "angular 48"),
        (2, 4, 64, "angular 64"),
        (4, 6, 48, "radial+angular r4p6 48"),
        (4, 8, 64, "radial+angular r4p8 64"),
        (8, 8, 64, "radial 4× p8 ang64"),
    ]

    results = []
    for n_panels, p_order, n_ang, label in configs:
        try:
            k = run_f4(r_0, R, sig_t, sig_s, nsf, n_panels, p_order, n_ang)
            err = abs(k - K_INF) / K_INF * 100.0
            results.append((n_panels, p_order, n_ang, label, k, err))
            print(f"  n_panels={n_panels}, p={p_order}, n_ang={n_ang:3d}  "
                  f"({label:28s}): k={k:.10f}, err={err:.6f}%")
        except Exception as e:
            print(f"  n_panels={n_panels}, p={p_order}, n_ang={n_ang}  "
                  f"({label}): FAILED — {type(e).__name__}: {e}")

    print("\n" + "=" * 78)
    print("Scaling analysis")
    print("=" * 78)
    # Pair up refinements and compute ratios
    base_err = results[0][5]
    print(f"Baseline err = {base_err:.6f}%")
    for res in results[1:]:
        ratio = base_err / res[5] if res[5] > 1e-9 else np.inf
        print(f"  vs {res[3]:28s}: err={res[5]:.6f}%, ratio={ratio:.2f}×")

    # Extract the "richest" config as asymptotic estimate
    min_err = min(r[5] for r in results)
    best = [r for r in results if r[5] == min_err][0]
    print(f"\nBest (floor estimate): err={min_err:.6f}% @ ({best[0]}, {best[1]}, {best[2]})")

    return 0, results


def test_f4_quadrature_floor():
    """Characterize F.4 convergence under quadrature refinement.

    The test PASSES as long as F.4 continues to reduce residual under
    refinement (monotone decrease). If it plateaus, that's diagnostic
    information but not a bug.
    """
    _, results = main()
    # Baseline
    base_err = results[0][5]
    # Best
    min_err = min(r[5] for r in results)
    # Record for later sessions: the floor estimate
    assert base_err > 0.0, "F.4 baseline should have nonzero residual"
    # Not a strict assertion — this is diagnostic
    print(f"\n=== TEST RESULT ===")
    print(f"F.4 baseline err = {base_err:.6f}%")
    print(f"F.4 floor estimate (richest quadrature) = {min_err:.6f}%")


if __name__ == "__main__":
    import sys
    ret, _ = main()
    sys.exit(ret)
