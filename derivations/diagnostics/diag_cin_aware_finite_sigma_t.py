"""
Extend the c_in-aware split basis to finite σ_t:
  - Compute W_{oi,s}^{mn}(τ) numerically for τ = σ_t R in a useful range.
  - Compute W_{io}^{mn}(τ) via reciprocity + time-reversal.
  - Verify Sanchez-McCormick reciprocity A_k W_{jk}^{mn} = A_j W_{kj}^{nm}.
  - Compute W_{gg} (grazing self-transmission) and its chord geometry.
  - Compute outer white BC decomposition: how the isotropic re-emission
    (constant distribution on [0,1] in µ) projects onto grazing/steep sub-bases.

Derivation context: see diag_cin_aware_basis_derivation.py (§[1]-[3] verified
symbolically). Key formulas:
  • Outer-steep to inner-via-chord transmission matrix:
      W_{oi,s}^{mn}(τ) = (1/ρ) ∫_0^1 P̃_m(c) P̃_n(c) exp(-τ s(c;ρ)) c dc,
    where s(c;ρ) = sqrt(1 - ρ²(1-c²)) - ρ c  (normalized chord length, τ=σ_t R).
  • Outer-grazing to outer-grazing: the chord stays within the shell with
    impact parameter b = R sqrt(1-µ²) > r_0, touching outer at another point.
    By 1D spherical symmetry, the arrival cosine equals the emission cosine;
    chord length L_gg(µ) = 2 R µ (radially-bounced shallow chord).
  • Inner to outer-steep: by ray reversibility, same chord length,
    mapping inner cosine c_I to outer cosine µ_R = sqrt(1 - ρ²(1-c²)).

All numerical integrations use scipy.integrate.quad for adaptive accuracy.
"""
from __future__ import annotations

import math
import numpy as np
import sympy as sp
from scipy import integrate

# ---------------------------------------------------------------------------
# Symbolic setup of half-range Legendre
# ---------------------------------------------------------------------------

mu_sym = sp.Symbol("mu", positive=True, real=True)


def _half_range_legendre_exprs(n_max):
    """Legendre P̃_n on [0, 1] w.r.t. µ-weight (Sanchez µ-ortho)."""
    basis = []
    for k in range(n_max + 1):
        fk = mu_sym**k
        for pj in basis:
            cjk = sp.integrate(fk * pj * mu_sym, (mu_sym, 0, 1))
            fk = fk - cjk * pj
        norm_sq = sp.integrate(fk * fk * mu_sym, (mu_sym, 0, 1))
        fk = sp.simplify(fk / sp.sqrt(norm_sq))
        basis.append(fk)
    return basis


def half_range_legendre_lambdas(n_max):
    """Return list of Python callables P̃_n(x) for n = 0..n_max."""
    exprs = _half_range_legendre_exprs(n_max)
    return [sp.lambdify(mu_sym, e, modules=["numpy"]) for e in exprs]


# ---------------------------------------------------------------------------
# Core geometric quantities
# ---------------------------------------------------------------------------

def mu_crit(rho: float) -> float:
    return math.sqrt(1.0 - rho**2)


def c_I(mu: float, rho: float) -> float:
    """Outer-to-inner angle remapping. Requires mu >= mu_crit."""
    return math.sqrt(1.0 - (1.0 / rho**2) * (1.0 - mu**2))


def mu_R(c: float, rho: float) -> float:
    """Inverse: inner-to-outer arrival cosine."""
    return math.sqrt(1.0 - rho**2 * (1.0 - c**2))


def chord_oi_normalized(c: float, rho: float) -> float:
    """
    Outer→inner chord length / R as a function of inner cosine c.
    s(c; ρ) = sqrt(1 - ρ²(1-c²)) - ρ c = µ_R(c) - ρ c.
    """
    return mu_R(c, rho) - rho * c


def chord_gg_normalized(mu: float) -> float:
    """
    Outer-grazing self-traversal chord length / R, for µ ∈ [0, µ_crit].
    The ray enters outer at cosine µ (to inward normal), stays in shell,
    exits outer at cosine µ (by chord symmetry in sphere), length 2Rµ.
    """
    return 2.0 * mu


# ---------------------------------------------------------------------------
# Transmission matrix elements: numerical integration
# ---------------------------------------------------------------------------

def W_oi_steep(m: int, n: int, tau: float, rho: float, leg: list) -> float:
    """
    W_{oi,s}^{mn}(τ) = (1/ρ) ∫_0^1 P̃_m(c) P̃_n(c) exp(-τ s(c;ρ)) c dc.
    At τ=0 this is (1/ρ) δ_{mn}.
    """
    def integrand(c):
        s = chord_oi_normalized(c, rho)
        return leg[m](c) * leg[n](c) * math.exp(-tau * s) * c

    val, _ = integrate.quad(integrand, 0.0, 1.0, epsabs=1e-12, epsrel=1e-10)
    return val / rho


def W_io(m: int, n: int, tau: float, rho: float, leg: list) -> float:
    """
    W_{io}^{mn}(τ) = inner-to-outer-steep transmission.

    By time-reversal of the outer→inner ray, a ray leaving inner at cosine c
    arrives at outer-steep at cosine µ_R(c) = sqrt(1 - ρ²(1-c²)), traversing
    the same chord length L = R s(c; ρ).

    The matrix element is:
      W_{io}^{mn}(τ) = ∫_0^1 P̃_m(c) · [P̃_n^s at µ_R(c), i.e., (1/ρ) P̃_n(c)]
                      · exp(-τ s(c;ρ)) · c dc
                   = (1/ρ) ∫_0^1 P̃_m(c) P̃_n(c) exp(-τ s(c;ρ)) c dc
                   = W_{oi,s}^{mn}(τ)

    Careful: this is the transmission integrated over the inner surface basis
    with c-weight. Sanchez-McCormick reciprocity involves surface-area factors.
    """
    # By time reversal, the chord and radiance are same; the formula is symmetric.
    return W_oi_steep(m, n, tau, rho, leg)


def W_gg(m: int, n: int, tau: float, rho: float, leg_g: list) -> float:
    """
    W_{gg}^{mn}(τ) = grazing-to-grazing outer self-transmission.

    On outer surface, grazing µ ∈ [0, µ_crit]. By chord symmetry, a grazing
    ray at emission µ arrives back at outer at cosine µ. Chord length = 2Rµ.

    So W_{gg}^{mn}(τ) = ∫_0^{µ_crit} P̃_m^g(µ) P̃_n^g(µ) exp(-τ · 2µ) µ dµ.
    Basis is µ-orthonormal on [0, µ_crit], so at τ=0 this gives δ_{mn}.

    leg_g is the list of grazing basis functions (SymPy-lambdified for
    µ ∈ [0, µ_crit]).
    """
    muc = mu_crit(rho)

    def integrand(mu):
        return leg_g[m](mu) * leg_g[n](mu) * math.exp(-tau * 2.0 * mu) * mu

    val, _ = integrate.quad(integrand, 0.0, muc, epsabs=1e-12, epsrel=1e-10)
    return val


# ---------------------------------------------------------------------------
# Grazing sub-basis on [0, µ_crit]
# ---------------------------------------------------------------------------

def _grazing_exprs(n_max, rho_val):
    """Symbolic Gram-Schmidt of {µ^k} on [0, µ_crit] (numerical rho)."""
    muc = math.sqrt(1.0 - rho_val**2)
    mu = sp.Symbol("mu", positive=True, real=True)
    basis = []
    for k in range(n_max + 1):
        fk = mu**k
        for pj in basis:
            cjk = sp.integrate(fk * pj * mu, (mu, 0, sp.Float(muc)))
            fk = fk - cjk * pj
        norm_sq = sp.integrate(fk * fk * mu, (mu, 0, sp.Float(muc)))
        fk = sp.simplify(fk / sp.sqrt(norm_sq))
        basis.append(fk)
    return basis, mu


def grazing_basis_lambdas(n_max, rho):
    exprs, mu = _grazing_exprs(n_max, rho)
    return [sp.lambdify(mu, e, modules=["numpy"]) for e in exprs]


# ---------------------------------------------------------------------------
# Sanchez-McCormick reciprocity check
# ---------------------------------------------------------------------------

def verify_reciprocity(n_max, tau, rho, leg):
    """
    Sanchez-McCormick:  A_k W_{jk}^{mn} = A_j W_{kj}^{nm}
    with A_outer = 4π R², A_inner = 4π r_0².

    For j=outer-steep (surface area 4π R²), k=inner (surface area 4π r_0²):
      A_inner · W_{oi}^{mn} = A_outer · W_{io}^{nm}
      r_0² W_{oi,s}^{mn} = R² W_{io}^{nm}
      W_{io}^{nm} = (r_0/R)² W_{oi,s}^{mn} = ρ² W_{oi,s}^{mn}

    But we derived W_{io}^{nm} = W_{oi,s}^{mn} (same formula!). This
    indicates a CONVENTION issue:
      - The "matrix element" W as we're computing it treats the basis
        functions symmetrically but does not include surface area factors.
      - Sanchez-McCormick W includes area normalization.

    To satisfy reciprocity we need to rescale one of them. This tells us
    the "physical" W_io^{nm} = ρ² · W_{oi,s}^{mn}, i.e., small.

    This is a BASIS NORMALIZATION point we need to work out carefully.
    """
    print(f"\n  Reciprocity check at τ={tau}, ρ={rho}:")
    print(f"    SanchezMcC: W_{{io}}^{{nm}} = ρ² W_{{oi,s}}^{{mn}} (area factor ρ²={rho**2:.4f})")
    # As-computed W_{io} and W_{oi,s} are equal (by construction in time-reversal).
    # For physical reciprocity we'd apply ρ² factor to one of them.
    for i in range(min(2, n_max + 1)):
        for j in range(min(2, n_max + 1)):
            w_oi = W_oi_steep(i, j, tau, rho, leg)
            w_io = W_io(j, i, tau, rho, leg)  # note transposed indices for reciprocity
            ratio_direct = w_io / w_oi if abs(w_oi) > 1e-14 else float("nan")
            print(f"    W_{{oi,s}}^{{{i},{j}}}={w_oi:+.6e}, W_{{io}}^{{{j},{i}}}={w_io:+.6e}, "
                  f"ratio={ratio_direct:+.6f} (unscaled)")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    print("\n=== c_in-aware split basis: finite-σ_t extension ===\n")

    rho = 0.3  # r_0 / R
    N_max = 3
    tau_list = [0.0, 0.5, 1.0, 2.0, 5.0, 10.0]

    print(f"  ρ = r_0/R = {rho}")
    print(f"  µ_crit = sqrt(1-ρ²) = {mu_crit(rho):.6f}")
    print(f"  Steep cone half-angle ≈ arccos(µ_crit) ≈ {math.degrees(math.acos(mu_crit(rho))):.2f}°")
    print(f"  Grazing fraction of outer hemisphere (solid angle): "
          f"(1 - µ_crit²) = ρ² = {rho**2:.4f} = {rho**2 * 100:.1f}%")
    print(f"  Steep fraction: 1 - ρ² = {1 - rho**2:.4f} = {(1 - rho**2) * 100:.1f}%")
    print(f"  (Note: 'grazing' covers a LARGER solid angle cos-range but smaller "
          f"area in angular current measure due to µ-weight)")

    leg = half_range_legendre_lambdas(N_max)  # on [0, 1]
    leg_g = grazing_basis_lambdas(N_max, rho)  # on [0, µ_crit]

    print(f"\n[1] W_{{oi,s}}^{{mn}}(τ) — should be diagonal at τ=0, near-diagonal elsewhere")
    print(f"    Diagonal value at τ=0 should be 1/ρ = {1/rho:.4f}\n")
    for tau in tau_list:
        W = np.zeros((N_max + 1, N_max + 1))
        for i in range(N_max + 1):
            for j in range(N_max + 1):
                W[i, j] = W_oi_steep(i, j, tau, rho, leg)
        print(f"    τ = {tau:6.2f}  (diag-dominance: off_L∞/diag_L∞ = "
              f"{(np.max(np.abs(W - np.diag(np.diag(W))))/max(np.max(np.abs(np.diag(W))), 1e-30)):.3e})")
        for i in range(N_max + 1):
            row = "      " + "  ".join(f"{W[i,j]:+.4e}" for j in range(N_max + 1))
            print(row)
        print()

    print(f"[2] W_{{gg}}^{{mn}}(τ) — grazing self-transmission, chord=2Rµ")
    for tau in [0.0, 1.0, 5.0]:
        W = np.zeros((N_max + 1, N_max + 1))
        for i in range(N_max + 1):
            for j in range(N_max + 1):
                W[i, j] = W_gg(i, j, tau, rho, leg_g)
        print(f"    τ = {tau:6.2f}")
        for i in range(N_max + 1):
            row = "      " + "  ".join(f"{W[i,j]:+.4e}" for j in range(N_max + 1))
            print(row)
        print()

    print(f"[3] Reciprocity check (convention diagnostic)")
    verify_reciprocity(N_max, 1.0, rho, leg)
    verify_reciprocity(N_max, 5.0, rho, leg)


if __name__ == "__main__":
    main()
