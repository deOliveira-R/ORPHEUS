"""Diagnostic: Issue #131 Probe D — P_esc_outer quadrature for multi-region slab.

Created by numerics-investigator on 2026-04-23.

Hypothesis: `compute_P_esc_outer` multi-region slab branch (line ~1619)
integrates on µ ∈ (-1, 1) with NO split at µ=0, then discards negative
µ inside the loop. This wastes half the quadrature nodes and — worse —
places node density symmetrically around µ=0 rather than where the
integrand `E₂`-like mass lives. Compare the two paths at σ_t uniform
(where the multi-region path SHOULD reduce to the homogeneous one).

Test:
  1. Build a "homogeneous" 2-region slab with identical σ_t in both
     regions.
  2. Call compute_P_esc_outer on it (triggering the multi-region GL
     branch at line 1619).
  3. Build the actual homogeneous slab with 1 region and call the
     closed-form branch at line 1609.
  4. Compare — if they differ by more than 1e-10, it's quadrature error.
  5. Scan n_angular and watch convergence (or non-convergence).
"""
from __future__ import annotations

import numpy as np
import pytest


@pytest.mark.slow
def test_issue131_probe_d_pesc_multiregion_quadrature():
    from orpheus.derivations.peierls_geometry import (
        SLAB_POLAR_1D,
        compute_P_esc_outer,
        compute_P_esc_inner,
        compute_G_bc_outer,
        compute_G_bc_inner,
        composite_gl_r,
    )

    # Same σ_t in both halves ≡ homogeneous slab with this σ_t.
    # Pick σ_t = 2.0 to match the shipped fixture's worst case
    # (Region B thermal group).
    sig_t_val = 2.0
    thicknesses = [0.5, 0.5]
    radii = np.cumsum(thicknesses)
    n_panels, p_order, dps = 2, 3, 20

    r_nodes, r_wts, _panels = composite_gl_r(
        np.asarray(radii, dtype=float), n_panels, p_order, dps=dps,
    )

    sig_t_2reg = np.array([sig_t_val, sig_t_val])
    sig_t_1reg = np.array([sig_t_val])
    radii_1reg = np.array([1.0])

    # Homogeneous: exact closed-form E_2.
    P_hom_outer = compute_P_esc_outer(
        SLAB_POLAR_1D, r_nodes, radii_1reg, sig_t_1reg,
        n_angular=24, dps=dps,
    )
    # Multi-region path, same physical σ_t: should match P_hom exactly.
    P_mr_outer = compute_P_esc_outer(
        SLAB_POLAR_1D, r_nodes, radii, sig_t_2reg,
        n_angular=24, dps=dps,
    )

    err_outer = np.max(np.abs(P_hom_outer - P_mr_outer))
    print(f"\nP_esc_outer homogeneous vs 2-region-same-σ_t (n_angular=24):")
    print(f"  max|err| = {err_outer:.3e}")
    for i, (ph, pm) in enumerate(zip(P_hom_outer, P_mr_outer)):
        print(f"  r_nodes[{i}]={float(r_nodes[i]):.4f}: "
              f"hom={ph:.10e}, 2reg={pm:.10e}, diff={ph - pm:.3e}")

    # Convergence scan
    print("\nConvergence of multi-region P_esc_outer in n_angular:")
    for N in [8, 16, 24, 32, 48, 64, 96, 128, 192, 256]:
        P_mr = compute_P_esc_outer(
            SLAB_POLAR_1D, r_nodes, radii, sig_t_2reg,
            n_angular=N, dps=dps,
        )
        err = np.max(np.abs(P_hom_outer - P_mr))
        print(f"  n_angular={N:4d}: max|err| = {err:.3e}")

    # Same test for inner face
    P_hom_inner = compute_P_esc_inner(
        SLAB_POLAR_1D, r_nodes, radii_1reg, sig_t_1reg,
        n_angular=24, dps=dps,
    )
    P_mr_inner = compute_P_esc_inner(
        SLAB_POLAR_1D, r_nodes, radii, sig_t_2reg,
        n_angular=24, dps=dps,
    )
    err_inner = np.max(np.abs(P_hom_inner - P_mr_inner))
    print(f"\nP_esc_inner homogeneous vs 2-region-same-σ_t (n_angular=24):")
    print(f"  max|err| = {err_inner:.3e}")

    # And G_bc (which uses µ ∈ (0,1) quadrature correctly).
    G_hom_outer = compute_G_bc_outer(
        SLAB_POLAR_1D, r_nodes, radii_1reg, sig_t_1reg,
        n_surf_quad=24, dps=dps,
    )
    G_mr_outer = compute_G_bc_outer(
        SLAB_POLAR_1D, r_nodes, radii, sig_t_2reg,
        n_surf_quad=24, dps=dps,
    )
    err_G_outer = np.max(np.abs(G_hom_outer - G_mr_outer))
    print(f"\nG_bc_outer homogeneous vs 2-region-same-σ_t (n_surf_quad=24):")
    print(f"  max|err| = {err_G_outer:.3e}")

    G_hom_inner = compute_G_bc_inner(
        SLAB_POLAR_1D, r_nodes, radii_1reg, sig_t_1reg,
        n_surf_quad=24, dps=dps,
    )
    G_mr_inner = compute_G_bc_inner(
        SLAB_POLAR_1D, r_nodes, radii, sig_t_2reg,
        n_surf_quad=24, dps=dps,
    )
    err_G_inner = np.max(np.abs(G_hom_inner - G_mr_inner))
    print(f"G_bc_inner homogeneous vs 2-region-same-σ_t (n_surf_quad=24):")
    print(f"  max|err| = {err_G_inner:.3e}")

    # If the multi-region GL path is underconverged we expect
    # err_outer to be large at N=24 and slowly decay, while err_G
    # (which uses µ ∈ (0, 1) directly) is tight.
    assert err_outer < 1e-8, (
        f"P_esc_outer multi-region quadrature error = {err_outer:.3e}. "
        f"This is the bug! Closed-form path (homogeneous) vs GL path "
        f"(multi-region) diverge even when σ_t is uniform."
    )


if __name__ == "__main__":
    test_issue131_probe_d_pesc_multiregion_quadrature()
