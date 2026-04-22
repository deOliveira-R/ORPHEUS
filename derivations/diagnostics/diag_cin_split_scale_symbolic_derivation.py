"""Symbolic derivation of the 1/ρ geometric factor in scale²_opt = 2 + 1/(3ρ).

Hypothesis: at σ_t → 0, the geometric probability integral for
"volume source → inner surface partial current" has a 1/ρ leading
asymptote that, combined with the Eddington factor 1/3, gives the
empirical scale formula.

Specifically, compute J^+_inner / A_inner / Q_vol:
  = (1/A_inner) · ∫_{shell} dr · 4π r² · P_esc_inner(r)
  = (1/A_inner) · 2π r_0³ · I(ρ)
  where I(ρ) = ∫_ρ^1 [1 - √(1-u²)]/u^4 du  (u = r_0/r)

For small ρ, I(ρ) ≈ (1-ρ²)^{3/2}/(3ρ³) - ... which might give 1/ρ³?

Let me compute I(ρ) exactly via sympy.
"""
from __future__ import annotations

import math
import sympy as sp


def main():
    print("=" * 80)
    print("Symbolic derivation: J^+_inner / A_inner / Q_vol at σ_t=0 in hollow sphere")
    print("=" * 80)

    rho, u = sp.symbols("rho u", positive=True, real=True)

    # Integrand: (1 - sqrt(1 - u²)) / u^4
    integrand = (1 - sp.sqrt(1 - u**2)) / u**4

    # Evaluate the integral from ρ to 1
    I_rho = sp.integrate(integrand, (u, rho, 1))
    I_rho_simplified = sp.simplify(I_rho)
    print(f"\n  I(ρ) = ∫_ρ^1 (1 - √(1-u²))/u^4 du")
    print(f"       = {I_rho_simplified}")
    print(f"  (expanded): {sp.expand(I_rho_simplified)}")

    # Limit as ρ → 0
    limit_small = sp.limit(I_rho, rho, 0, '+')
    print(f"\n  Lim I(ρ) as ρ → 0: {limit_small}")

    # Leading behavior
    series_at_0 = sp.series(I_rho_simplified, rho, 0, 4)
    print(f"  Series at ρ=0: {series_at_0}")

    # J^+_inner / Q at inner surface:
    #   J^+_inner = Q · ∫_{shell} (1 - √(1-(r_0/r)²))/2 · 4π r² dr
    #             = Q · 2π r_0³ · I(ρ)  (from u = r_0/r substitution)
    # A_inner = 4π r_0²
    # J^+_inner / A_inner = Q · r_0/2 · I(ρ) = Q · R·ρ/2 · I(ρ)
    print(f"\n  J^+_inner / A_inner / Q = R · ρ · I(ρ) / 2")
    print(f"    = R · {sp.simplify(rho * I_rho_simplified / 2)}")

    result = sp.simplify(rho * I_rho_simplified / 2)
    series_result = sp.series(result, rho, 0, 4)
    print(f"  Series at ρ=0: {series_result}")

    # Compare to: is it ≈ R/4 for small ρ? Or R/(something·ρ)?
    print(f"\n  At ρ=0 limit: {sp.limit(result, rho, 0, '+')}")
    print(f"  At ρ=1 limit: {sp.limit(result, rho, 1, '-')}")

    # Now: P_0_inner (the constant-basis projection) integral:
    # P_0_inner(r) = ∫_0^1 (exit-prob at c) · c dc
    # This is more complex — needs the angular distribution of rays reaching inner.

    # At σ_t=0, a ray from r in direction with cosine µ_r (to inward normal) reaches
    # inner if |µ_r| ≥ cos(θ_view(r)). Reaches with inner cosine c_r = some value.
    # Actually in free flight, the inner-arriving cosine is determined by geometry:
    # by impact-parameter conservation: r · sin θ_r = r_0 · sin θ_{r_0}
    # cos θ_{r_0} = √(1 - (r/r_0)² sin²θ_r) = √(1 - (r/r_0)²(1 - µ_r²))
    # Where µ_r is the radial-outward cosine at r.

    # Wait, this gets complex because r varies over the shell, and each r has
    # different inner-arriving cosine distribution.

    # For a source at radius r emitting in direction (θ_r, φ_r) with Ω·r̂ > 0
    # (outward from center), ray goes outward and exits at outer.
    # For Ω·r̂ < 0 (inward), ray goes to inner if grazing angle tight enough.
    # At the inner surface, cosine to inward-normal (= -r̂_{r_0}):
    # c_{r_0} = Ω · (-r̂_{r_0})

    # Skip the full derivation — too much for a quick sympy probe.
    # The qualitative observation: I(ρ) in sympy should reveal the 1/ρ-type behavior.


if __name__ == "__main__":
    main()
