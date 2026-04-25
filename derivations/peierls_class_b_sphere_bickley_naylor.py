"""SymPy derivation: analytical white-BC closure for solid Class B sphere.

Created 2026-04-25 for Issue #101 / #132 resolution. The rank-N Marshak
closure has been falsified for Class B (this session) and Class A
hollow (L21 close-out). The canonical alternative is **analytical
evaluation** of the surface response operator — historically written
in terms of Bickley-Naylor functions Ki_n, but careful: those are the
canonical form for CYLINDER (where polar-angle integration over the
cylinder axis introduces the Ki_n integral). For SPHERE, the natural
closed form may instead be in terms of EXPONENTIAL INTEGRALS E_n.

This script derives the spherical P_esc(r) and G_bc(r) from first
principles via SymPy and identifies which canonical functions appear.
The result feeds into:

- ``orpheus/derivations/peierls_geometry.py`` — new
  ``compute_P_esc_bickley`` / ``compute_G_bc_bickley`` primitives
- ``boundary="white_bickley"`` closure routing through
  ``build_closure_operator``
- ``tests/derivations/test_peierls_class_b_bickley.py`` — pin
  Class B sphere k_eff vs ``cp_sphere._build_case`` analytical
  k_inf at < 1e-4 tolerance, MR×MG-capable

Multi-region extension follows the slab Issue #131 template:
piecewise sum of the closed-form integrand across material breakpoints
along each ray.

Mathematical setup
------------------
Sphere of radius R, observer at radius r ∈ [0, R], uniform Σ_t (single-
region first; multi-region in §4 below). The pointwise escape
probability for an isotropic source at r is

.. math::

    P_{\\rm esc}(r) = \\frac{1}{4\\pi}\\int_{4\\pi} e^{-\\Sigma_t \\rho(\\Omega)} d\\Omega
                    = \\frac{1}{2}\\int_0^\\pi e^{-\\Sigma_t \\rho(\\theta)} \\sin\\theta\\,d\\theta

where ρ(θ) is the chord length from r in direction θ to the surface,
solved from |r·r̂ + ρ·Ω̂|² = R²:

.. math::

    \\rho(\\theta) = -r\\cos\\theta + \\sqrt{R^2 - r^2\\sin^2\\theta}

ρ(θ) ranges from R-r (forward, θ=0) to R+r (backward, θ=π) as θ varies.
"""

from __future__ import annotations

import sympy as sp


# ═══════════════════════════════════════════════════════════════════════
# Symbols
# ═══════════════════════════════════════════════════════════════════════

theta, u, t, rho, R, r, Sigma_t = sp.symbols(
    "theta u t rho R r Sigma_t",
    positive=True, real=True,
)


# ═══════════════════════════════════════════════════════════════════════
# Step 1 — Direct angular integration
# ═══════════════════════════════════════════════════════════════════════

print("=" * 76)
print("Step 1: Direct angular form of P_esc(r)")
print("=" * 76)

rho_of_theta = -r * sp.cos(theta) + sp.sqrt(R**2 - r**2 * sp.sin(theta)**2)
print(f"\n   ρ(θ) = {rho_of_theta}")

P_esc_angular = sp.Rational(1, 2) * sp.Integral(
    sp.exp(-Sigma_t * rho_of_theta) * sp.sin(theta),
    (theta, 0, sp.pi),
)
print()
print(f"   P_esc(r) = {P_esc_angular}")


# ═══════════════════════════════════════════════════════════════════════
# Step 2 — Substitution u = cos θ
# ═══════════════════════════════════════════════════════════════════════

print()
print("=" * 76)
print("Step 2: Substitute u = cos θ → simplifies sin θ dθ → -du")
print("=" * 76)

rho_of_u = -r * u + sp.sqrt(R**2 - r**2 * (1 - u**2))
rho_of_u_simplified = sp.simplify(rho_of_u)
print(f"\n   ρ(u) = {rho_of_u_simplified}")

P_esc_u = sp.Rational(1, 2) * sp.Integral(
    sp.exp(-Sigma_t * rho_of_u),
    (u, -1, 1),
)
print()
print(f"   P_esc(r) = {P_esc_u}")
print()
print("   At u=+1 (forward): ρ =", sp.simplify(rho_of_u.subs(u, 1)))
print("   At u=-1 (backward): ρ =", sp.simplify(rho_of_u.subs(u, -1)))
print("   At u=0 (perpendicular): ρ =", sp.simplify(rho_of_u.subs(u, 0)))


# ═══════════════════════════════════════════════════════════════════════
# Step 3 — Substitute t = ρ(u) (chord length as the integration variable)
# ═══════════════════════════════════════════════════════════════════════
#
# From ρ(u) = -ru + sqrt(R²-r²+r²u²), we have:
#   ρ + ru = sqrt(R²-r²+r²u²)
#   ρ² + 2ρru + r²u² = R²-r²+r²u²
#   ρ² + 2ρru = R²-r²
#   u = (R²-r²-ρ²)/(2ρr)
#
# du/dρ = -(R²-r²+ρ²)/(2ρ²r)  (algebraic differentiation below)
#
# As u goes from 1 to -1, ρ goes from R-r to R+r monotonically.

print()
print("=" * 76)
print("Step 3: Substitute t = ρ as the integration variable")
print("=" * 76)

u_of_t = (R**2 - r**2 - t**2) / (2 * t * r)
print(f"\n   u(t) = {u_of_t}")
print(f"   u(R-r) = {sp.simplify(u_of_t.subs(t, R - r))}     (should be +1)")
print(f"   u(R+r) = {sp.simplify(u_of_t.subs(t, R + r))}     (should be -1)")

du_dt = sp.diff(u_of_t, t)
du_dt_simplified = sp.simplify(du_dt)
print(f"\n   du/dt = {du_dt_simplified}")
print()
print("   ⇒ |du/dt| = (R²-r²+t²)/(2 r t²)  for t ∈ [R-r, R+r]")
print()

# Set up the t-substituted integral. The direction reversal (u: 1→-1
# as t: R-r→R+r) cancels with the negative du/dt, so the integral
# limits and integrand both come out positive.
integrand_t = sp.exp(-Sigma_t * t) * (R**2 - r**2 + t**2) / (2 * r * t**2)
P_esc_t_integrand = sp.Rational(1, 2) * integrand_t
P_esc_t = sp.Integral(P_esc_t_integrand, (t, R - r, R + r))
print(f"   P_esc(r) = (1/(4r)) · ∫_(R-r)^(R+r) e^(-Σ_t·t) · (R²-r²+t²)/t² dt")
print()
print(f"   = (1/(4r)) · [∫ e^(-Σ_t·t) dt + (R²-r²) · ∫ e^(-Σ_t·t)/t² dt]")


# ═══════════════════════════════════════════════════════════════════════
# Step 4 — Reduce to E_n form
# ═══════════════════════════════════════════════════════════════════════
#
# Two integrals to evaluate:
#
#   I_1 = ∫_a^b e^(-Σ_t t) dt = (1/Σ_t) [e^(-Σ_t a) - e^(-Σ_t b)]
#
#   I_2 = ∫_a^b e^(-Σ_t t) / t² dt
#       — integration by parts: u = e^(-Σ_t t), dv = dt/t²
#         I_2 = [-e^(-Σ_t t)/t]_a^b - Σ_t · ∫_a^b e^(-Σ_t t)/t dt
#             = e^(-Σ_t a)/a - e^(-Σ_t b)/b - Σ_t · (E_1(Σ_t a) - E_1(Σ_t b))
#
# where a = R-r, b = R+r. Note R²-r² = (R-r)(R+r) = a·b.

print()
print("=" * 76)
print("Step 4: Closed-form reduction via E_n")
print("=" * 76)

a, b = sp.symbols("a b", positive=True)  # a = R-r, b = R+r

I_1 = (sp.exp(-Sigma_t * a) - sp.exp(-Sigma_t * b)) / Sigma_t
print(f"\n   I_1 = ∫_a^b e^(-Σ_t t) dt = {I_1}")

# E_1 in SymPy via integration:
# Use Ei or expint. mpmath uses E_n(x) = ∫_1^∞ e^(-xt)/t^n dt.
# SymPy: sympy.expint(n, x) = E_n(x).
E1_a = sp.expint(1, Sigma_t * a)
E1_b = sp.expint(1, Sigma_t * b)
I_2 = (
    sp.exp(-Sigma_t * a) / a - sp.exp(-Sigma_t * b) / b
    - Sigma_t * (E1_a - E1_b)
)
print(f"   I_2 = ∫_a^b e^(-Σ_t t)/t² dt")
print(f"       = e^(-Σ_t·a)/a - e^(-Σ_t·b)/b - Σ_t·[E_1(Σ_t·a) - E_1(Σ_t·b)]")

# Substitute back:  R² - r² = a·b,  4r = 2(b-a)
P_esc_closed = (I_1 + a * b * I_2) / (2 * (b - a))
P_esc_closed_simplified = sp.simplify(P_esc_closed)
print()
print(f"   P_esc(r) = (1/(2(b-a))) · [I_1 + a·b·I_2]")
print(f"            (with a = R-r, b = R+r, so 2(b-a) = 4r)")
print()
print(f"   Expanded:")
print(f"     = (1/(4r·Σ_t)) · [(e^(-Σ_t·a) - e^(-Σ_t·b))")
print(f"                         + (a·b·Σ_t/a)·e^(-Σ_t·a)")
print(f"                         - (a·b·Σ_t/b)·e^(-Σ_t·b)")
print(f"                         - a·b·Σ_t²·(E_1(Σ_t·a) - E_1(Σ_t·b))]")
print(f"")
print(f"     = (1/(4r·Σ_t)) · [(1 + Σ_t·b)·e^(-Σ_t·a) - (1 + Σ_t·a)·e^(-Σ_t·b)")
print(f"                         - Σ_t²·a·b·(E_1(Σ_t·a) - E_1(Σ_t·b))]")


# ═══════════════════════════════════════════════════════════════════════
# Step 5 — Recognise the E_2 / E_3 form
# ═══════════════════════════════════════════════════════════════════════
#
# Use E_n recursion: E_n(x) = (e^(-x) - x·E_{n-1}(x)) / (n-1) for n ≥ 2.
# Equivalently: x·E_{n-1}(x) = e^(-x) - (n-1)·E_n(x).
#
# So Σ_t·a·E_1(Σ_t·a) = e^(-Σ_t·a) - E_2(Σ_t·a).
# Substituting:
#   Σ_t·a·b·E_1(Σ_t·a) = b·(e^(-Σ_t·a) - E_2(Σ_t·a))
#
# and similarly for b. Hence:
#
#   Σ_t²·a·b·(E_1(Σ_t·a) - E_1(Σ_t·b))
#     = Σ_t·b·(e^(-Σ_t·a) - E_2(Σ_t·a)) - Σ_t·a·(e^(-Σ_t·b) - E_2(Σ_t·b))
#     = Σ_t·b·e^(-Σ_t·a) - Σ_t·a·e^(-Σ_t·b) - Σ_t·b·E_2(Σ_t·a) + Σ_t·a·E_2(Σ_t·b)
#
# Substitute into the bracket:
#   [(1 + Σ_t·b)·e^(-Σ_t·a) - (1 + Σ_t·a)·e^(-Σ_t·b)
#    - Σ_t·b·e^(-Σ_t·a) + Σ_t·a·e^(-Σ_t·b)
#    + Σ_t·b·E_2(Σ_t·a) - Σ_t·a·E_2(Σ_t·b)]
#
#   = e^(-Σ_t·a) - e^(-Σ_t·b)
#       + Σ_t·b·E_2(Σ_t·a) - Σ_t·a·E_2(Σ_t·b)
#
# Hence:
#
#   P_esc(r) = (1/(4r·Σ_t)) · [
#                 e^(-Σ_t·(R-r)) - e^(-Σ_t·(R+r))
#                 + Σ_t·(R+r)·E_2(Σ_t·(R-r))
#                 - Σ_t·(R-r)·E_2(Σ_t·(R+r))
#             ]
#
# This is the canonical closed form for the spherical white-BC escape
# probability — uses E_2 (exponential integral), not Ki_n. The Ki_n
# form is for cylinder; the literature researcher's reference to "Ki_2"
# in the original analysis was a recall error.

print()
print("=" * 76)
print("Step 5: Algebraic simplification → E_2 closed form")
print("=" * 76)
print()
print("   Using E_n recursion x·E_{n-1}(x) = e^(-x) - (n-1)·E_n(x):")
print()
print("   ┌───────────────────────────────────────────────────────────────┐")
print("   │                                                               │")
print("   │   P_esc(r) = (1/(4r·Σ_t)) · [                                 │")
print("   │     e^(-Σ_t·(R-r)) - e^(-Σ_t·(R+r))                           │")
print("   │     + Σ_t·(R+r)·E_2(Σ_t·(R-r))                                │")
print("   │     - Σ_t·(R-r)·E_2(Σ_t·(R+r))                                │")
print("   │   ]                                                           │")
print("   │                                                               │")
print("   └───────────────────────────────────────────────────────────────┘")
print()
print("   Canonical form: spherical pointwise escape probability uses E_2,")
print("   NOT Ki_n. The Bickley-Naylor functions are CYLINDER-specific")
print("   (where polar-angle integration over the axis introduces them).")


# ═══════════════════════════════════════════════════════════════════════
# Step 6 — Numerical verification
# ═══════════════════════════════════════════════════════════════════════

print()
print("=" * 76)
print("Step 6: Numerical verification — closed form vs direct quadrature")
print("=" * 76)

import mpmath
import scipy.integrate as si


def P_esc_closed_form(r_val: float, R_val: float, Sigma_t_val: float, dps: int = 50):
    """Closed-form P_esc using E_2."""
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


def P_esc_direct(r_val: float, R_val: float, Sigma_t_val: float):
    """Direct quadrature via scipy: (1/2) ∫_0^π e^(-Σ_t ρ(θ)) sin θ dθ."""
    def integrand(th):
        rho_th = -r_val * mpmath.cos(th) + mpmath.sqrt(
            R_val**2 - r_val**2 * mpmath.sin(th)**2,
        )
        return float(mpmath.exp(-Sigma_t_val * float(rho_th)) * mpmath.sin(th))
    val, _ = si.quad(integrand, 0, mpmath.pi, limit=200, epsabs=1e-15, epsrel=1e-13)
    return 0.5 * val


print()
print(f"  R=1.0, Σ_t=1.0:")
print(f"  {'r':>6} {'P_esc closed':>18} {'P_esc direct':>18} {'rel diff':>12}")
print("  " + "─" * 60)
for r_val in [0.001, 0.1, 0.3, 0.5, 0.7, 0.9, 0.99]:
    p_closed = P_esc_closed_form(r_val, 1.0, 1.0)
    p_direct = P_esc_direct(r_val, 1.0, 1.0)
    rel = abs(p_closed - p_direct) / max(abs(p_closed), 1e-30)
    print(
        f"  {r_val:>6.3f} {p_closed:>18.12f} {p_direct:>18.12f} {rel:>12.3e}"
    )

print()
print(f"  R=1.0, Σ_t=5.0 (thicker cell):")
print(f"  {'r':>6} {'P_esc closed':>18} {'P_esc direct':>18} {'rel diff':>12}")
print("  " + "─" * 60)
for r_val in [0.001, 0.1, 0.5, 0.9]:
    p_closed = P_esc_closed_form(r_val, 1.0, 5.0)
    p_direct = P_esc_direct(r_val, 1.0, 5.0)
    rel = abs(p_closed - p_direct) / max(abs(p_closed), 1e-30)
    print(
        f"  {r_val:>6.3f} {p_closed:>18.12f} {p_direct:>18.12f} {rel:>12.3e}"
    )

print()
print(f"  Limit check r → 0 (point at center):")
print(f"  Expected: P_esc(0) = exp(-Σ_t · R)")
for Sigma_t_val in [0.5, 1.0, 5.0]:
    expected = float(mpmath.exp(-Sigma_t_val * 1.0))
    actual = P_esc_closed_form(1e-8, 1.0, Sigma_t_val)
    print(f"    Σ_t={Sigma_t_val}: expected={expected:.10f}  closed-form(r=1e-8)={actual:.10f}")


# ═══════════════════════════════════════════════════════════════════════
# Step 7 — G_bc(r) derivation
# ═══════════════════════════════════════════════════════════════════════
#
# For uniform isotropic INWARD partial current J⁻ on the sphere surface,
# the angular flux at the surface in the inward hemisphere is:
#
#   ψ⁻(r_b, Ω) = J⁻ / π    for Ω·n_b < 0
#
# (since J⁻ = ∫_{Ω·n<0} |Ω·n| ψ⁻ dΩ = π·ψ⁻ for isotropic distribution).
#
# The scalar flux at interior r from this surface source — using observer-
# centred angular integration:
#
#   φ(r) = ∫_{4π} ψ⁻(r_b(Ω), Ω) e^(-τ(r,Ω)) dΩ
#
# where r_b(Ω) is the surface point reached by going from r in
# direction Ω (the same chord-length geometry as P_esc but now reading
# the angular flux on the surface). For the unit J⁻ = 1:
#
#   φ(r) = (1/π) · ∫_{4π} e^(-Σ_t ρ(Ω)) dΩ   (over directions where ψ⁻ ≠ 0)
#        = (1/π) · 4π · P_esc(r)              (the angular integral is 4π · P_esc)
#        = 4 · P_esc(r)
#
# Wait — that's not right. Let me redo. ψ⁻(r_b, Ω) is nonzero only for
# Ω·n_b < 0 (incoming). Going from interior r BACK along Ω̂ to r_b means
# Ω̂ at r_b points INTO the sphere, so Ω·n_b < 0 ✓. The integration is
# over directions FROM r (over 4π), and the ray reaches r_b on the
# surface. The angular flux at r_b in direction +Ω (pointing into the
# sphere from r_b's perspective) is J⁻/π.
#
# So φ(r) = ∫_{4π} (J⁻/π) · e^(-Σ_t ρ(Ω)) dΩ = (J⁻/π) · 4π · 2 · P_esc(r)
#        = 8 · J⁻ · P_esc(r)? No that's not right either.
#
# Actually, let me think more carefully:
#   φ(r) = ∫_{4π} ψ(r, Ω) dΩ   [scalar flux = integral of angular flux]
#        = ∫_{4π} ψ⁻(r_b(Ω), Ω) · e^(-τ(r,r_b)) dΩ
#
# Each direction Ω from r corresponds to ONE surface point r_b(Ω) and
# ONE incoming angular flux value. The Ω·n_b condition at r_b: the ray
# arrives at r_b traveling outward (Ω̂ aligned with positive radial
# direction at r_b for all θ if r > 0), so Ω̂·n_b > 0 at r_b. But we
# want INCOMING current at r_b, which has Ω̂·n_b < 0 — opposite sign.
#
# Resolution: for a neutron at interior r in direction Ω, the
# "history" came from r_b in the OPPOSITE direction -Ω. So we read
# ψ⁻(r_b, -Ω) which has -Ω·n_b < 0 ⇒ Ω·n_b > 0. Wait that's still
# outgoing. I'm getting confused with sign conventions.

print()
print("=" * 76)
print("Step 7: G_bc(r) for uniform isotropic surface source — TBD")
print("=" * 76)
print()
print("   Derivation deferred pending literature-researcher confirmation")
print("   of the conventional sign convention. Should be a similar")
print("   E_2 closed form by symmetry. See follow-up commit.")
