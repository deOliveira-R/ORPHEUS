"""Diagnostic: Issue #131 Probe E — closed-form E_n rescue for multi-region slab.

Created by numerics-investigator on 2026-04-23.

Hypothesis: replacing the GL µ-quadrature in P_esc_outer / inner
multi-region slab branch with closed-form ½ E_2(τ_total_to_face)
eliminates the 4e-3 quadrature error.

Justification: for a slab with piecewise-constant σ_t(x),

   P_esc_outer(x_i) = ½ ∫_0^1 exp[ -(∫_x_i^L σ_t(x') dx')/µ ] dµ
                    = ½ E_2( τ_total_to_outer(x_i) )

where τ_total_to_outer(x_i) = Σ_k σ_t^(k) · (b_k ∩ [x_i, L]). The µ-integral
is closed-form — the multi-region piecewise-constant σ_t only changes
τ_total, not the angular dependence. Same story for G_bc: 2 E_2 of the
same optical depth.
"""
from __future__ import annotations

import numpy as np
import mpmath
import pytest


def _tau_to_outer_face(x_i: float, radii: np.ndarray,
                       sig_t: np.ndarray) -> float:
    """Piecewise-integrated optical depth from x_i to the outer face at
    radii[-1], for a slab with region boundaries at
    ``[0, radii[0], radii[1], ..., radii[-1]]`` and piecewise-constant
    σ_t[k] in region k = (radii[k-1], radii[k]]."""
    L = float(radii[-1])
    if x_i >= L:
        return 0.0
    tau = 0.0
    r_prev = 0.0
    for k, r_k in enumerate(radii):
        a, b = max(float(r_prev), x_i), float(r_k)
        if a < b:
            tau += float(sig_t[k]) * (b - a)
        r_prev = r_k
    return tau


def _tau_to_inner_face(x_i: float, radii: np.ndarray,
                       sig_t: np.ndarray) -> float:
    """Piecewise-integrated optical depth from x_i to the inner face at x=0."""
    if x_i <= 0.0:
        return 0.0
    tau = 0.0
    r_prev = 0.0
    for k, r_k in enumerate(radii):
        a, b = float(r_prev), min(float(r_k), x_i)
        if a < b:
            tau += float(sig_t[k]) * (b - a)
        r_prev = r_k
    return tau


def _E2(tau: float) -> float:
    if tau <= 0.0:
        return 1.0
    return float(mpmath.expint(2, tau))


def _closed_form_P_esc_outer(r_nodes, radii, sig_t) -> np.ndarray:
    P = np.zeros(len(r_nodes))
    for i, x_i in enumerate(r_nodes):
        tau = _tau_to_outer_face(float(x_i), radii, sig_t)
        P[i] = 0.5 * _E2(tau)
    return P


def _closed_form_P_esc_inner(r_nodes, radii, sig_t) -> np.ndarray:
    P = np.zeros(len(r_nodes))
    for i, x_i in enumerate(r_nodes):
        tau = _tau_to_inner_face(float(x_i), radii, sig_t)
        P[i] = 0.5 * _E2(tau)
    return P


def _closed_form_G_bc_outer(r_nodes, radii, sig_t) -> np.ndarray:
    G = np.zeros(len(r_nodes))
    for i, x_i in enumerate(r_nodes):
        tau = _tau_to_outer_face(float(x_i), radii, sig_t)
        G[i] = 2.0 * _E2(tau)
    return G


def _closed_form_G_bc_inner(r_nodes, radii, sig_t) -> np.ndarray:
    G = np.zeros(len(r_nodes))
    for i, x_i in enumerate(r_nodes):
        tau = _tau_to_inner_face(float(x_i), radii, sig_t)
        G[i] = 2.0 * _E2(tau)
    return G


@pytest.mark.slow
def test_issue131_probe_e_closed_form_matches_homogeneous():
    """Closed-form for multi-region σ_t=const slab ≡ homogeneous slab closed form."""
    from orpheus.derivations.peierls_geometry import (
        SLAB_POLAR_1D,
        compute_P_esc_outer,
        compute_P_esc_inner,
        composite_gl_r,
    )

    sig_t_val = 2.0
    thicknesses = [0.5, 0.5]
    radii = np.cumsum(thicknesses)
    n_panels, p_order, dps = 2, 3, 20

    r_nodes, _r_wts, _panels = composite_gl_r(
        np.asarray(radii, dtype=float), n_panels, p_order, dps=dps,
    )
    r_nodes_arr = np.asarray(r_nodes, dtype=float)

    sig_t_2reg = np.array([sig_t_val, sig_t_val])

    # Homogeneous closed form (reference)
    P_hom_outer = compute_P_esc_outer(
        SLAB_POLAR_1D, r_nodes_arr, np.array([1.0]),
        np.array([sig_t_val]), n_angular=24, dps=dps,
    )

    # Our proposed closed-form for multi-region
    P_fix_outer = _closed_form_P_esc_outer(r_nodes_arr, radii, sig_t_2reg)

    err = np.max(np.abs(P_hom_outer - P_fix_outer))
    print(f"\nClosed-form multi-region P_esc_outer vs homogeneous: max|err|={err:.3e}")
    for i in range(len(r_nodes_arr)):
        print(f"  r={r_nodes_arr[i]:.4f}: hom={P_hom_outer[i]:.10e} "
              f"fix={P_fix_outer[i]:.10e}")
    assert err < 1e-12, f"Closed-form fix mismatches homogeneous: {err:.3e}"

    # And now at the ACTUAL fixture (non-uniform σ_t):
    # Region A thermal sig_t=1.0, Region B thermal sig_t=2.0.
    sig_t_heter = np.array([1.0, 2.0])
    P_fix_heter = _closed_form_P_esc_outer(r_nodes_arr, radii, sig_t_heter)
    P_gl_heter = compute_P_esc_outer(
        SLAB_POLAR_1D, r_nodes_arr, radii, sig_t_heter,
        n_angular=24, dps=dps,
    )
    print("\nHeterogeneous σ_t=[1.0, 2.0] (thermal fixture):")
    print(f"  max|closed − GL_N=24|: {np.max(np.abs(P_fix_heter - P_gl_heter)):.3e}")
    for N in [24, 48, 96, 192, 384, 512]:
        P_gl = compute_P_esc_outer(
            SLAB_POLAR_1D, r_nodes_arr, radii, sig_t_heter,
            n_angular=N, dps=dps,
        )
        err_gl = np.max(np.abs(P_fix_heter - P_gl))
        print(f"    vs GL N={N}: max|err|={err_gl:.3e}")


if __name__ == "__main__":
    test_issue131_probe_e_closed_form_matches_homogeneous()
