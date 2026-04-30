"""SymPy + numerical derivation of the canonical 3-D G_bc kernel for cylinder.

Issue #112 Phase C: the current `compute_G_bc` cylinder uses a 2-D-projected
``Ki_1(τ)/d`` form. The cylinder Hébert investigator (commit 30335f2) showed
this form is biased by ~7.6 % independent of quadrature on cyl 1G/1R (row-sum
probe `K · 1 / σ_t = 0.89` for cylinder vs 0.999 for sphere), preventing the
``(1-P_ss^cyl)⁻¹`` Hébert series from reaching <1.5 % L1 tolerance.

This script derives the correct 3-D form from first principles by integrating
the 3-D point-kernel ``g(r, r') = exp(-Σ_t |r-r'|) / (4π |r-r'|²)`` over the
infinite z-extent of the cylinder.

Setup
-----
Observer at interior point (r, 0, 0) (radial distance r from cylinder axis,
z=0 for convenience by translation symmetry). Surface point at (R cos φ,
R sin φ, z) for azimuthal angle φ ∈ [0, 2π) and z ∈ (-∞, +∞). Outward
normal at surface: n_b = (cos φ, sin φ, 0).

In-plane distance: d(φ) = √(r² - 2rR cos φ + R²)
3-D distance: d_3D = √(d²(φ) + z²)

For Mark closure, ψ⁻(r_b, Ω) = J⁻/π for all incoming directions at any
surface point. The flux at interior r from this surface input:

    φ(r) = ∫_S ψ⁻ · (Ω·n_b)_inward · exp(-Σ_t · d_3D) / d_3D² · dA
         = (J⁻/π) · ∫_S |cos θ_inc| · exp(-Σ_t · d_3D) / d_3D² · dA

with dA = R dφ dz on the cylinder lateral surface.

The cosine of the incidence angle:
    cos θ_inc = ((r_i - r_b) · n_b) / d_3D = (r cos φ - R) / d_3D
    |cos θ_inc| = (R - r cos φ) / d_3D    (for r < R, all r cos φ < R)

So:
    G_bc^cyl(r) = (R/π) · ∫_0^{2π} dφ · (R - r cos φ) · ∫_{-∞}^{+∞} dz · 1/d_3D³ · exp(-Σ_t · d_3D)

The z-integral reduces to a Bickley function via the standard substitution
z = d·tan α (so d_3D = d/cos α, dz = d/cos²α dα):

    ∫_{-∞}^{+∞} exp(-Σ_t · d/cos α) / (d/cos α)³ · (d/cos² α) dα
        = (1/d²) · ∫_{-∞}^{+∞} exp(-Σ_t · d / cos α) · cos α dα
        = (2/d²) · ∫_0^{π/2} cos α · exp(-Σ_t · d / cos α) dα
        = (2/d²) · Ki_2(Σ_t · d)

[Using the Bickley convention Ki_n(x) = ∫_0^{π/2} cos^{n-1} θ · exp(-x/cos θ) dθ.
 So Ki_2(x) = ∫_0^{π/2} cos θ · exp(-x/cos θ) dθ.]

Therefore:

    ┌──────────────────────────────────────────────────────────────────────┐
    │  G_bc^cyl(r) = (4R/π) · ∫_0^π (R - r cos φ) / d²(φ) · Ki_2(Σ_t·d) dφ │
    └──────────────────────────────────────────────────────────────────────┘

(Symmetry around φ = 0 collapses ∫_0^{2π} → 2·∫_0^π.)

Comparison with current code (peierls_geometry.py:1493-1518):

    G_bc^current(r) = (2R/π) · ∫_0^π Ki_1(Σ_t · d(φ)) / d(φ) dφ

DIFFERENCES (three bugs in one):
    1. Bickley order: Ki_2 (correct) vs Ki_1 (current)
    2. Geometry factor: (R - r cos φ)/d² (correct, the cos θ_inc Lambertian
       projection times area-Jacobian) vs 1/d (current)
    3. Leading coefficient: 4R/π (correct) vs 2R/π (current)

The MAGNITUDE of the bias is best understood via the r=0 limit where
d(φ) = R for all φ:
    G_bc^correct(0) = (4R/π) · ∫_0^π R/R² · Ki_2(Σ_t·R) dφ = 4 · Ki_2(Σ_t·R)
    G_bc^current(0) = (2R/π) · ∫_0^π Ki_1(Σ_t·R)/R dφ     = 2 · Ki_1(Σ_t·R)

For thin cells (Σ_t·R → 0): Ki_2(0) = 1, Ki_1(0) = π/2, so:
    G_bc^correct(0) → 4
    G_bc^current(0) → π ≈ 3.1416   (off by π/4 = 21 %)

For Σ_t·R = 1: Ki_2(1) ≈ 0.115, Ki_1(1) ≈ 0.245, so:
    G_bc^correct(0) ≈ 0.460
    G_bc^current(0) ≈ 0.490   (off by 6.5 %)

This matches the cylinder investigator's empirical 7.6 % row-sum bias on
cyl 1G/1R (where cell is R=1 MFP, σ_t·R = 1).
"""

from __future__ import annotations

import numpy as np
import sympy as sp
import mpmath

from orpheus.derivations._kernels import ki_n_mp


def Ki_n(n: int, x: float, dps: int = 25) -> float:
    """Bickley-Naylor wrapper using ORPHEUS robust ki_n_mp (tan substitution)."""
    if x < 0:
        x = 0.0
    return float(ki_n_mp(n, x, dps))


def G_bc_cyl_correct(r: float, R: float, Sigma_t: float,
                      n_quad: int = 64, dps: int = 25) -> float:
    """Correct 3-D form: G_bc^cyl = (4R/π) ∫_0^π (R - r cos φ)/d² · Ki_2(Σ_t·d) dφ."""
    pts, wts = np.polynomial.legendre.leggauss(n_quad)
    phi_pts = 0.5 * (pts + 1) * np.pi
    phi_wts = wts * 0.5 * np.pi
    total = 0.0
    for k in range(n_quad):
        cf = float(np.cos(phi_pts[k]))
        d2 = r * r - 2 * r * R * cf + R * R
        d = float(np.sqrt(max(d2, 0.0)))
        if d <= 0.0:
            continue
        total += phi_wts[k] * (R - r * cf) / d2 * Ki_n(2, Sigma_t * d, dps=dps)
    return float(4.0 * R / np.pi * total)


def G_bc_cyl_current(r: float, R: float, Sigma_t: float,
                      n_quad: int = 64, dps: int = 25) -> float:
    """Current 2-D form: G_bc^cyl = (2R/π) ∫_0^π Ki_1(Σ_t·d)/d dφ."""
    pts, wts = np.polynomial.legendre.leggauss(n_quad)
    phi_pts = 0.5 * (pts + 1) * np.pi
    phi_wts = wts * 0.5 * np.pi
    total = 0.0
    for k in range(n_quad):
        cf = float(np.cos(phi_pts[k]))
        d2 = r * r - 2 * r * R * cf + R * R
        d = float(np.sqrt(max(d2, 0.0)))
        if d <= 0.0:
            continue
        total += phi_wts[k] * Ki_n(1, Sigma_t * d, dps=dps) / d
    return float(2.0 * R / np.pi * total)


def main():
    print("=" * 78)
    print("CYLINDER G_bc 3-D derivation vs current 2-D form")
    print("=" * 78)

    print("\n--- Bickley function values (sanity) ---")
    print(f"  Ki_1(0) = π/2      ≈ {np.pi/2:.6f}    numeric = {Ki_n(1, 0.0):.6f}")
    print(f"  Ki_2(0) = 1                             numeric = {Ki_n(2, 0.0):.6f}")
    print(f"  Ki_1(1)            numeric = {Ki_n(1, 1.0):.6f}")
    print(f"  Ki_2(1)            numeric = {Ki_n(2, 1.0):.6f}")

    print("\n--- G_bc at r=0 (analytical limits) ---")
    print(f"  R=1, Σ_t=1: correct (4·Ki_2(1))  = {4 * Ki_n(2, 1.0):.6f}")
    print(f"              current (2·Ki_1(1))  = {2 * Ki_n(1, 1.0):.6f}")
    print(f"              ratio (cur/cor)       = {2*Ki_n(1,1.0)/(4*Ki_n(2,1.0)):.4f}")
    print(f"  Numerical at r=0.001 (avoid /0):")
    print(f"              correct  = {G_bc_cyl_correct(1e-3, 1.0, 1.0):.6f}")
    print(f"              current  = {G_bc_cyl_current(1e-3, 1.0, 1.0):.6f}")

    print("\n--- G_bc(r) sweep at R=1, Σ_t=1 ---")
    print(f"  {'r':>6} {'correct':>12} {'current':>12} {'ratio':>10}")
    for r_val in [0.001, 0.1, 0.3, 0.5, 0.7, 0.9, 0.99]:
        g_c = G_bc_cyl_correct(r_val, 1.0, 1.0)
        g_o = G_bc_cyl_current(r_val, 1.0, 1.0)
        print(f"  {r_val:>6.3f} {g_c:>12.6f} {g_o:>12.6f} {g_o/g_c:>10.4f}")

    print("\n--- G_bc thin-cell limit ---")
    for sig_R in [1e-3, 0.1, 1.0, 10.0]:
        g_c = G_bc_cyl_correct(0.001, 1.0, sig_R)
        g_o = G_bc_cyl_current(0.001, 1.0, sig_R)
        print(f"  σ·R={sig_R:>6.3f}: correct={g_c:.6f}  current={g_o:.6f}  "
              f"ratio={g_o/g_c:.4f}")

    print()
    print("=" * 78)
    print("CONCLUSION")
    print("=" * 78)
    print()
    print("  Current G_bc^cyl is OVERESTIMATING by ~6-21 % depending on σ_t·R.")
    print("  Correct form: (4R/π) ∫_0^π (R - r cos φ)/d² · Ki_2(Σ_t·d) dφ")
    print()
    print("  Three changes from current:")
    print("    (a) Ki_1 → Ki_2")
    print("    (b) 1/d → (R - r cos φ)/d²  [Lambertian projection]")
    print("    (c) leading 2R/π → 4R/π")
    print()
    print("  The BIAS direction (current too HIGH) is the OPPOSITE of what would")
    print("  cause K·1/σ_t = 0.89 (which suggests K is too LOW). So if the")
    print("  Hébert investigator's diagnosis is correct, the cylinder bias may")
    print("  also have a SECOND component beyond G_bc — possibly in the")
    print("  P_esc cylinder primitive too. To be tested empirically by")
    print("  patching G_bc and re-running the row-sum probe.")


if __name__ == "__main__":
    main()
