"""Sanity gate: closed-form E_2 P_esc vs shipped rank-1 Mark compute_P_esc.

The SymPy derivation in
``scratch/derivations/peierls_class_b_sphere_bickley_naylor.py`` produced a
closed form for the spherical pointwise escape probability:

    P_esc(r) = (1/(4r·Σ_t)) · [
       e^(-Σ_t·(R-r)) - e^(-Σ_t·(R+r))
       + Σ_t·(R+r)·E_2(Σ_t·(R-r))
       - Σ_t·(R-r)·E_2(Σ_t·(R+r))
    ]

This script verifies the closed form matches the shipped GL-quadrature
``compute_P_esc`` to within quadrature precision. If it does, the
closed form is the multi-region-extensible analytical replacement for
the existing rank-1 Mark primitive (analog of slab E_2 piecewise sum,
Issue #131 template).

**Important caveat** (recorded for next-session work): replacing the
GL-quadrature rank-1 Mark P_esc with the analytical E_2 form does NOT
by itself fix the rank-1 Mark CLOSURE error (~27 % at sphere R=1 MFP).
The closure structure (rank-1 outer product u_n ⊗ v_n with isotropic
re-entry) is what causes the Mark error; the analytical evaluation
just removes the GL-quadrature contribution to that error. The full
fix for Class B verification requires either:

(a) A new white-BC kernel via method-of-images that absorbs the
    surface reflection analytically into the volume kernel itself
    (Davison sphere white-BC kernel — needs literature confirmation
    that closed form exists)
(b) Augmented Nyström with the surface partial current as an extra
    unknown (no Mark approximation at all)
(c) Iterative power iteration on the surface source

This script's CURRENT scope: confirm the analytical P_esc primitive
is correct; defer the full closure question to the literature pull.
"""

from __future__ import annotations

import numpy as np
import mpmath

from orpheus.derivations.peierls_geometry import (
    SPHERE_1D,
    compute_P_esc,
)


def P_esc_closed_form(r_val: float, R_val: float, Sigma_t_val: float, dps: int = 50):
    """Closed-form P_esc for spherical solid cell using E_2."""
    a_val = mpmath.mpf(R_val) - mpmath.mpf(r_val)
    b_val = mpmath.mpf(R_val) + mpmath.mpf(r_val)
    Sa = mpmath.mpf(Sigma_t_val) * a_val
    Sb = mpmath.mpf(Sigma_t_val) * b_val
    with mpmath.workdps(dps):
        bracket = (
            mpmath.exp(-Sa) - mpmath.exp(-Sb)
            + Sb * mpmath.expint(2, Sa)
            - Sa * mpmath.expint(2, Sb)
        )
        return float(bracket / (4 * mpmath.mpf(r_val) * mpmath.mpf(Sigma_t_val)))


def main():
    print("=" * 76)
    print("Sanity gate: analytical E_2 P_esc vs shipped GL-quadrature")
    print("=" * 76)

    R_val = 1.0
    Sigma_t_val = 1.0
    radii = np.array([R_val])
    sig_t = np.array([Sigma_t_val])

    # Sample some interior radii
    r_test = np.array([0.01, 0.1, 0.3, 0.5, 0.7, 0.9, 0.95])

    for n_angular in [16, 32, 64, 128]:
        P_shipped = compute_P_esc(SPHERE_1D, r_test, radii, sig_t,
                                   n_angular=n_angular, dps=25)
        print(f"\n  n_angular = {n_angular}:")
        print(f"  {'r':>6} {'P_closed':>16} {'P_shipped':>16} {'rel_diff':>12}")
        for i, r_val in enumerate(r_test):
            p_c = P_esc_closed_form(r_val, R_val, Sigma_t_val)
            p_s = float(P_shipped[i])
            rel = abs(p_c - p_s) / max(abs(p_c), 1e-30)
            print(f"  {r_val:>6.3f} {p_c:>16.12f} {p_s:>16.12f} {rel:>12.3e}")

    print()
    print("=" * 76)
    print("CONCLUSION")
    print("=" * 76)
    print()
    print("  If the closed-form values match the shipped values to better than")
    print("  the GL-quadrature error budget at each n_angular, the analytical")
    print("  E_2 form is verified as a correct replacement for the rank-1")
    print("  Mark P_esc primitive.")
    print()
    print("  Implementing this in production requires:")
    print("    1. Multi-region extension (Σ_t piecewise → piecewise τ sums)")
    print("    2. New compute_P_esc_E2 / compute_G_bc_E2 primitives")
    print("    3. boundary='white_E2_mark' wiring (still rank-1 Mark closure")
    print("       structure, just with analytical primitives)")
    print()
    print("  This does NOT improve the rank-1 Mark closure error (~27% at")
    print("  R=1 MFP). For Class B verification at <1% tolerance, the")
    print("  Davison white-BC kernel via method-of-images is the deeper")
    print("  path — pending literature researcher confirmation.")


if __name__ == "__main__":
    main()
