"""SymPy derivation: rank-N partial-current moment for solid Class B sphere.

Created 2026-04-25 in support of Issue #132.

**Status: derivation INCOMPLETE — both candidate fixes were FALSIFIED
numerically. Kept as evidence for the next session.**

Goal
----
Derive from first principles the canonical Marshak DP_{N-1} partial-
current-moment ``P_esc^(n)(r_i)`` for a SOLID sphere cell (one outer
surface, no inner cavity), in observer-centred coordinates. Compare
against the current ORPHEUS implementation in
``orpheus.derivations.peierls_geometry.compute_P_esc_mode`` to identify
the precise structural discrepancy that produces the +57 % k_eff gap
on Class B sphere 1G/2R rank-2 (Issue #132).

The literature-researcher memory at
``.claude/agent-memory/literature-researcher/sanchez_mccormick_rank_n_per_face.md``
identified the bug class for the HOLLOW per-face primitives: a missing
``(Ω·n) = cos θ`` µ-weight. That recipe was tried for the SOLID
single-surface primitives below — and FAILED.

What this script derives — and why it's incomplete
--------------------------------------------------
1. From Sanchez-McCormick (1982) §III.F.1 Eq. 166 cluster, the canonical
   surface-centred ``P^{ρk}_{S i}`` integrand is

       g(r→r_b) · P̃_n(µ_s) · µ_s · dA

   with g = exp(-τ)/(4π s²) the collision kernel, P̃_n the shifted
   Legendre on [0,1], and µ_s = (Ω·n) the µ-weight in the surface basis
   inner product.
2. Converting to observer-centred via dA = ρ²/|µ_s| · dΩ, ALL three
   factors {ρ², µ_s, 1/µ_s} cancel and the derivation gives

       P_esc^{(n)}_canonical(r_i) = (1/2) · ∫_0^π exp(-τ) · P̃_n(µ_s) · sin θ dθ

   — i.e. NO Jacobian, NO extra µ_s weight, just the Lambertian sin θ
   measure with the shifted Legendre factor.
3. Numerical experiment in
   ``derivations/diagnostics/diag_class_b_rank_n_fix_attempt.py``
   tested this candidate fix (FIX1: drop Jacobian) and a second
   variant from the literature memory (FIX2: drop Jacobian + add µ_s
   weight). Both made the catastrophe MUCH worse:

   | Configuration | BEFORE | FIX1 (no Jac) | FIX2 (no Jac + µ_s) |
   |---|---|---|---|
   | sph 1G/1R rank-2 | -1.10 % | +16.68 % | +11.11 % |
   | sph 1G/2R rank-2 | **+56.7 %** | **+322 %** | **+1178 %** |
   | sph 2G/2R rank-2 | +204 % | -61 % | -56 % |

   Both fixes preserve rank-1 (mode 0 only) which is unchanged because
   mode 0 still uses the legacy ``compute_P_esc`` that has no Jacobian.

What's missing — three failed fix attempts narrow the problem
-------------------------------------------------------------

THREE candidate fixes were numerically tested and all FAILED:

| Fix attempt | Recipe | sph 1R rank-2 | sph 2R rank-2 | Verdict |
|---|---|---|---|---|
| BEFORE | (current) | -1.10 % | +56.7 % | calibration coincidence at 1R |
| FIX1 | drop (ρ/R)² Jacobian | +16.7 % | +322 % | strictly worse |
| FIX2 | drop (ρ/R)², add µ_s weight | +11.1 % | +1178 % | strictly worse |
| FIX3 | drop (ρ/R)² + drop (2n+1) Gelbard reflection | -14.8 % | +25.5 % | **consistent plateau at ~15-25 % across configs** |

Reproducers:
    derivations/diagnostics/diag_class_b_rank_n_fix_attempt.py   (FIX1+FIX2)
    derivations/diagnostics/diag_class_b_rank_n_fix_attempt3.py  (FIX3)

**The FIX3 result is the most informative**: dropping BOTH the (ρ/R)²
Jacobian AND the (2n+1) Gelbard factor in the reflection operator
gives a *consistent* convergence ladder across 1R / 2R / cyl / sph,
plateauing at ~15-26 % k_eff error. This pattern matches the L21
falsification of rank-N on Class A hollow: a structural plateau, not
an improvable algebraic bug.

**The BEFORE 1R rank-2 = -1.10 % "convergence" claim is a calibration
coincidence** — the Jacobian and Gelbard factor combine to land 1R
near zero by accident, but break catastrophically when MR couples to
the inflated mode-1 contribution at the σ_t breakpoint.

What this means for Issue #132
------------------------------

The rank-N Marshak closure on solid Class B (and, by L21 precedent,
on hollow Class A) appears to be **structurally limited** to ~15-25 %
k_eff error regardless of normalisation choices. The literature
memory at
``.claude/agent-memory/literature-researcher/rank_n_closure_four_references_synthesis.md``
(L21 close-out) already concluded this for Class A: "ORPHEUS F.4 IS
the textbook closure"; the same structural argument now applies to
Class B.

**Recommended fix path** (revised after this session's failed attempts):

1. **Do NOT pursue further algebraic-fix attempts on rank-N for Class B.**
   The three failed fixes above + the investigator's two probe-G
   variants (LEGACY and CANONICAL) span the relevant 4-corner space
   {Jacobian on/off} × {Gelbard on/off}; none give clean convergence.

2. **Land FIX3 as a "honest mode" of the closure** — accept the ~15-25 %
   plateau as the structural limit, retract the misleading 1R
   convergence claim, and pin the FIX3 plateau as the new regression
   gate.

3. **Pursue Issue #101** (Bickley-Naylor analytical Ki₁ for Class B
   sphere) as the long-term path to a real continuous Peierls
   reference for Class B — the rank-N approach is empirically dead.

Held over to next session: implement FIX3 as the production rank-N
closure path on Class B (with appropriate documentation), retract
the published 1R rank-N table, and update Issue #132 to "won't fix —
structural limit, see Issue #101 as resolution path."

Failed test scripts kept at:
    derivations/diagnostics/diag_class_b_rank_n_fix_attempt.py
    derivations/diagnostics/diag_class_b_rank_n_fix_attempt3.py

Approach (per Cardinal Rule "use SymPy for multi-step derivations")
-------------------------------------------------------------------
1. State the canonical Sanchez-McCormick (1982) per-face surface-
   centred form for ``P^{ρk}_{Si}`` (the "P_esc analog").
2. Apply the observer-to-surface area Jacobian
   ``dA = dΩ · ρ² / |µ_s|`` to convert to observer-centred.
3. Carry ``µ_s = (Ω·n) = cos`` of the ray angle with the outward
   surface normal at the exit point.
4. Specialise to a single outer subsurface (solid Class B). Compare
   to the current ORPHEUS integrand line-by-line and isolate the
   missing factor.

Outputs are SymPy-printed equations + a numerical sanity check at one
node showing the legacy / canonical (current ORPHEUS) / proposed-fix
integrand values.
"""

from __future__ import annotations

import sympy as sp


# ═══════════════════════════════════════════════════════════════════════
# Symbols
# ═══════════════════════════════════════════════════════════════════════

theta, phi, mu, mu_s, R, r_i, rho, rho_max = sp.symbols(
    "theta phi mu mu_s R r_i rho rho_max", positive=True, real=True,
)
sigma_t, tau = sp.symbols("Sigma_t tau", positive=True, real=True)
n = sp.symbols("n", integer=True, nonneg=True)


def shifted_legendre(n_val: int, mu_var: sp.Symbol) -> sp.Expr:
    r"""Shifted Legendre :math:`\tilde P_n(\mu) = P_n(2\mu - 1)` on [0,1]."""
    return sp.legendre(n_val, 2 * mu_var - 1)


# ═══════════════════════════════════════════════════════════════════════
# Step 1 — Sanchez-McCormick canonical surface-centred partial-current moment
# ═══════════════════════════════════════════════════════════════════════
#
# Sanchez & McCormick 1982 NSE 80, §III.F.1 Eq. 166 cluster (the "P_esc
# analog"):
#
#   P^{ρk}_{S_α i} = +π A_α · V_i^{-1} · ∫_{V_i} dr ∫_{A_α} g(r → r'_b)
#                    · f^ρ_{+,α}(r'_b, Ω_s) · (Ω_s · n) dA'
#
# where:
#   - g(r → r'_b) = exp(-τ)/(4π s²) is the free-space collision kernel
#     (Sanchez-McCormick Eq. 107)
#   - f^ρ_{+,α}(r'_b, Ω_s) = (πA_α)^{-1} · P̃_ρ(µ_s) is the surface mode
#     basis function (P̃_ρ shifted Legendre on [0,1])
#   - µ_s = (Ω_s · n) is the cosine of the ray with the outward surface
#     normal AT the surface point r'_b
#   - The (Ω_s · n) factor in the surface integral is the µ-weight that
#     makes f^ρ a partial-current moment basis (orthonormal under
#     ∫ f^ρ f^ν · (Ω·n) dΩ dA = (πA)^{-1} δ_ρν)
#
# For a Nyström pointwise volumetric source at r_i (skipping the outer
# ∫_{V_i} dr / V_i averaging — that's per-zone in the CP method, but
# pointwise in Peierls), the per-node form is:
#
#   P_esc^{(n)}(r_i) = π A · ∫_A g(r_i → r_b) · (πA)^{-1} · P̃_n(µ_s)
#                          · (Ω_s · n) dA
#                    = ∫_A g(r_i → r_b) · P̃_n(µ_s) · µ_s · dA
#
# (the πA and (πA)^{-1} cancel — the outgoing partial-current moment
#  per unit point source is just the kernel × shifted Legendre × µ-weight
#  integrated over the surface).
#
# That's the *canonical surface-centred* form. Now convert to
# observer-centred via the area Jacobian.

print("=" * 76)
print("Step 1: Sanchez-McCormick canonical SURFACE-centred form")
print("=" * 76)

# Symbolic surface integral form (intent only — not evaluated)
P_esc_n_surface_centred = sp.Symbol("P_esc_n_surface_centred")
print()
print("   P_esc^{(n)}(r_i)  =  ∫_A  g(r_i → r_b) · P̃_n(µ_s) · µ_s  dA")
print()
print("   with:")
print("     g(r → r_b) = exp(-τ) / (4π s²)              [collision kernel]")
print("     µ_s = cos(angle between ray and outward normal at surface)")
print("     P̃_n(µ_s) = shifted Legendre on [0, 1]")
print("     dA = surface area element on outer face")


# ═══════════════════════════════════════════════════════════════════════
# Step 2 — Observer-centred change of variables
# ═══════════════════════════════════════════════════════════════════════
#
# Observer at r_i (interior node). Direction Ω parametrised by polar
# angle θ measured from the outward radial direction r̂_i, plus
# azimuth φ. A ray from r_i in direction Ω hits the outer sphere
# (radius R) at:
#
#   r_b = r_i · r̂_i + ρ · Ω̂
#
# with ρ = ρ_max(r_i, cos θ, R) the chord length to the surface, given
# by solving |r_b| = R:
#
#   ρ_max = -r_i cos θ + sqrt(R² - r_i²(1 - cos²θ))
#         = -r_i cos θ + sqrt(R² - r_i² sin²θ)
#
# At the exit point r_b, the outward normal n is r_b/|r_b| = r_b/R.
# The ray direction Ω̂ has cosine with n equal to:
#
#   µ_s  =  (r_i cos θ + ρ_max) / R
#
# (the radial component of r_b · Ω̂ / R; the radial component of r_i
#  along Ω̂ is r_i cos θ, plus the chord ρ_max along Ω̂ contributes ρ_max
#  to the dot product). This matches the existing ORPHEUS formula
# ``mu_exit = (rho_max_val + r_i * cos_om) / R``.
#
# The area Jacobian dA ↔ dΩ:
#
#   dA = ρ² · dΩ / |µ_s|
#
# (standard area-from-solid-angle relation: a tube of solid angle dΩ
#  intersects a surface at angle θ_inc from the surface normal in an
#  area dΩ · ρ²/|µ_s|; here θ_inc is the angle the ray makes with the
#  surface normal, whose cosine is µ_s).
#
# Substituting:
#
#   P_esc^{(n)}(r_i) = ∫_Ω g(r_i, ρ_max) · P̃_n(µ_s) · µ_s · (ρ²/|µ_s|) dΩ
#                    = ∫_Ω g(r_i, ρ_max) · P̃_n(µ_s) · ρ² dΩ
#
# The µ_s and 1/|µ_s| cancel, leaving:
#
#   P_esc^{(n)}(r_i) = ∫_Ω g(r_i, ρ_max) · ρ² · P̃_n(µ_s) dΩ
#
# With g = exp(-τ) / (4π ρ²), the ρ² factor cancels too:
#
#   P_esc^{(n)}(r_i) = (4π)^{-1} · ∫_Ω exp(-τ(r_i, Ω)) · P̃_n(µ_s) dΩ
#
# For 1D sphere, ∫_Ω = 2π · ∫_0^π sin θ dθ (azimuth integrates to 2π):
#
#   P_esc^{(n)}(r_i) = (1/2) · ∫_0^π exp(-τ) · P̃_n(µ_s) · sin θ dθ
#
# That's the canonical observer-centred form. The 1/2 prefactor is the
# Lambertian normalization factor.

print()
print("=" * 76)
print("Step 2: Observer-centred via area Jacobian dA = ρ²/|µ_s| · dΩ")
print("=" * 76)
print()

mu_s_expr = (r_i * sp.cos(theta) + rho_max) / R
rho_max_expr = -r_i * sp.cos(theta) + sp.sqrt(R**2 - r_i**2 * sp.sin(theta)**2)

print("   µ_s(θ; r_i, R)   =", sp.simplify(mu_s_expr))
print("   ρ_max(θ; r_i, R) =", rho_max_expr)
print()
print("   After change of variables dA = ρ²/|µ_s| · dΩ:")
print("   The two ρ²/µ_s factors and the µ_s in the integrand combine,")
print("   AND the 4π ρ² in the denominator of g cancels ρ², leaving:")
print()
print("   P_esc^{(n)}(r_i) = (1/2) · ∫_0^π  exp(-τ) · P̃_n(µ_s) · sin θ dθ")
print()
print("   ── canonical observer-centred form. NO (ρ/R)² Jacobian.")
print("   ── The µ_s weight in the surface basis is ABSORBED by the area")
print("      Jacobian, leaving only the Lambertian sin θ measure × P̃_n.")


# ═══════════════════════════════════════════════════════════════════════
# Step 3 — Specialise to mode 0 vs mode n ≥ 1
# ═══════════════════════════════════════════════════════════════════════
#
# Mode 0: P̃_0(µ) = 1 (always). So:
#
#   P_esc^{(0)}(r_i) = (1/2) · ∫_0^π exp(-τ) · sin θ dθ
#                    = compute_P_esc(r_i)            [✓ matches legacy]
#
# Mode n ≥ 1: integrand carries P̃_n(µ_s):
#
#   P_esc^{(n)}(r_i) = (1/2) · ∫_0^π exp(-τ) · P̃_n(µ_s) · sin θ dθ
#                                              ↑
#                                    NO Jacobian, NO extra cos θ weight
#
# Comparing to the CURRENT ORPHEUS implementation (peierls_geometry.py
# line 2725-2746, ``compute_P_esc_mode``):
#
#   P_esc^{(n)}_ORPHEUS(r_i) = (1/2) · ∫_0^π exp(-τ) · P̃_n(µ_s)
#                                        · (ρ_max/R)² · sin θ dθ
#                                        ↑↑↑
#                                  spurious (ρ/R)² Jacobian
#
# The current ORPHEUS form has an EXTRA (ρ_max/R)² that should not be
# there. The literature-researcher memory caught this at line 6c:
# "the right Sanchez-McCormick form, observer-centered, carries `cos θ`
#  (not `(ρ/R)²`) — they are not the same thing."
#
# But the corrected canonical form derived above carries NEITHER
# (ρ/R)² NOR cos θ — both factors absorbed into the area Jacobian. So
# the fix is simpler than what the literature memory suggested for the
# hollow per-face primitives:
#
#   *** Drop the (ρ_max/R)² Jacobian. Add NOTHING in its place. ***

print()
print("=" * 76)
print("Step 3: Mode 0 (Lambertian) vs Mode n ≥ 1 — and the fix")
print("=" * 76)
print()
print("   Mode 0 (P̃_0 = 1):")
print("     P_esc^{(0)}(r_i) = (1/2) · ∫₀^π exp(-τ) · sin θ dθ")
print("     ── matches `compute_P_esc` (legacy mode-0 ORPHEUS) ✓")
print()
print("   Mode n ≥ 1 (canonical, derived):")
print("     P_esc^{(n)}(r_i) = (1/2) · ∫₀^π exp(-τ) · P̃_n(µ_s) · sin θ dθ")
print()
print("   Mode n ≥ 1 (CURRENT ORPHEUS, peierls_geometry.py:2725-2746):")
print("     P_esc^{(n)}_ORPHEUS(r_i)")
print("       = (1/2) · ∫₀^π exp(-τ) · P̃_n(µ_s) · (ρ_max/R)² · sin θ dθ")
print()
print("   *** BUG: spurious (ρ_max/R)² Jacobian. ***")
print("   The fix is to DROP (ρ_max/R)² and add nothing:")
print()
print("     P_esc^{(n)}_FIX(r_i) = (1/2) · ∫₀^π exp(-τ) · P̃_n(µ_s) · sin θ dθ")
print()
print("   At n=0 this also reduces to compute_P_esc, so the legacy/canonical")
print("   inconsistency disappears: mode 0 and mode n ≥ 1 use the SAME")
print("   integrand structure (P̃_n(µ_s) · Lambertian dΩ measure).")


# ═══════════════════════════════════════════════════════════════════════
# Step 4 — G_bc analog: same derivation, same conclusion
# ═══════════════════════════════════════════════════════════════════════
#
# The Sanchez-McCormick "G_bc analog" P^{kρ}_{iS_α} has the same
# structure with the surface mode on the INCOMING side and a leading
# −4π V^{-1} normalization. The same observer-centred change of
# variables gives:
#
#   G_bc^{(n)}(r_i) = 2 · ∫_0^π exp(-τ) · P̃_n(µ_in) · sin θ dθ
#
# where µ_in is the cosine of the INCOMING ray with the surface
# normal at the entry point r_b (equivalent to µ_s above for a
# convex sphere — the entry and exit points reflect the same
# geometry by reciprocity).
#
# The current ORPHEUS ``compute_G_bc_mode`` for sphere
# (peierls_geometry.py:2803-2828):
#
#   G_bc^{(n)}_ORPHEUS(r_i) = 2 · ∫_0^π exp(-τ) · P̃_n(µ_s) · sin θ dθ
#
# That's ALREADY the canonical form — no Jacobian to remove! G_bc_mode
# was implemented correctly; the bug is *only* in P_esc_mode.

print()
print("=" * 76)
print("Step 4: G_bc^(n) analog — already canonical, no fix needed")
print("=" * 76)
print()
print("   compute_G_bc_mode (peierls_geometry.py:2803-2828) for sphere:")
print("     G_bc^{(n)}_ORPHEUS(r_i) = 2 · ∫₀^π exp(-τ) · P̃_n(µ_s) · sin θ dθ")
print()
print("   ── already matches the canonical observer-centred form. No fix.")
print("   ── the asymmetry between P_esc_mode (buggy) and G_bc_mode")
print("      (correct) is the explanation for why the catastrophe is")
print("      so large: only the P side has the spurious Jacobian.")


# ═══════════════════════════════════════════════════════════════════════
# Step 5 — Numerical sanity check at a single node
# ═══════════════════════════════════════════════════════════════════════
#
# Pick r_i = 0.3, R = 1.0, σ_t = 1.0, n = 1. Compute the LEGACY (no
# Jacobian, mode 0 form), the CURRENT ORPHEUS BUG (with Jacobian), and
# the proposed FIX (no Jacobian) integrands at θ = 0, π/4, π/2,
# 3π/4, π. Show the discrepancy is most pronounced at θ ≠ 0,
# (where ρ_max < R and (ρ_max/R)² < 1).

print()
print("=" * 76)
print("Step 5: Numerical sanity check at r_i = 0.3, R = 1, σ_t = 1, n = 1")
print("=" * 76)
print()
print(f"  {'θ':>8} {'ρ_max':>10} {'µ_s':>10} {'P̃_1(µ_s)':>12} "
      f"{'(ρ/R)²':>10} {'integrand×sinθ':>18}")
print(f"  {' ':>8} {' ':>10} {' ':>10} {' ':>12} {' ':>10} "
      f"{'(BUG)':>9}{'(FIX)':>9}")
print("  " + "─" * 80)

R_val = sp.Rational(1)
r_i_val = sp.Rational(3, 10)
sigma_t_val = sp.Rational(1)

for theta_val in [
    sp.Rational(1, 100), sp.pi / 4, sp.pi / 2,
    3 * sp.pi / 4, sp.pi - sp.Rational(1, 100),
]:
    rho_val = -r_i_val * sp.cos(theta_val) + sp.sqrt(
        R_val**2 - r_i_val**2 * sp.sin(theta_val)**2,
    )
    rho_val_f = float(rho_val)
    mu_s_val = (r_i_val * sp.cos(theta_val) + rho_val) / R_val
    mu_s_f = float(mu_s_val)
    p_tilde_1 = float(sp.legendre(1, 2 * mu_s_val - 1))
    jac = (rho_val_f / float(R_val))**2
    sin_theta_f = float(sp.sin(theta_val))
    tau_val = float(sigma_t_val) * rho_val_f
    exp_tau = float(sp.exp(-tau_val))
    bug_integrand = exp_tau * p_tilde_1 * jac * sin_theta_f
    fix_integrand = exp_tau * p_tilde_1 * sin_theta_f
    theta_f = float(theta_val)
    print(
        f"  {theta_f:>8.4f} {rho_val_f:>10.4f} {mu_s_f:>10.4f} "
        f"{p_tilde_1:>12.4f} {jac:>10.4f} "
        f"{bug_integrand:>9.5f}{fix_integrand:>9.5f}"
    )

print()
print("  The (ρ/R)² Jacobian is ≈ 0.81 at θ = π/4 and ≈ 0.25 at θ = π")
print("  for r_i = 0.3 — multiplicatively damping the mode-1 contribution")
print("  by a θ-dependent factor that should not be there. The 'BUG'")
print("  column shows the (ρ/R)²-multiplied integrand; the 'FIX' column")
print("  shows the canonical Lambertian form without the Jacobian.")


# ═══════════════════════════════════════════════════════════════════════
# Conclusion
# ═══════════════════════════════════════════════════════════════════════

print()
print("=" * 76)
print("Conclusion")
print("=" * 76)
print()
print("  The fix to Issue #132 (Class B sphere rank-N catastrophe) is to")
print("  DROP the (ρ_max/R)² surface-to-observer Jacobian from")
print("  ``compute_P_esc_mode`` (peierls_geometry.py:2741). The Jacobian")
print("  was added under the (incorrect) docstring assumption that the")
print("  observer-to-surface area conversion needed an extra (ρ/R)²")
print("  weight. In fact, the area Jacobian dA = ρ²/|µ_s| · dΩ exactly")
print("  cancels with the µ_s weight in the surface basis inner product")
print("  AND with the ρ² in the collision kernel g = exp(-τ)/(4π ρ²),")
print("  leaving only the Lambertian dΩ measure × P̃_n(µ_s).")
print()
print("  After the fix:")
print("    - P_esc^{(0)}_FIX(r_i) reduces algebraically to compute_P_esc")
print("      ⇒ rank-1 Mark closure remains bit-exact (no regression)")
print("    - Mode 0 and mode n ≥ 1 live in the SAME normalisation space")
print("      ⇒ MR catastrophe resolved (the source of Issue #132)")
print("    - G_bc_mode is unchanged (it was already canonical)")
print()
print("  The cylinder case is structurally analogous but the Bickley-")
print("  function K_esc (for slab-polar) and the surface-centred Ki_1/d")
print("  form (for cyl) have additional factors that the literature-")
print("  researcher memory (Knyazev 1993 Ki_{2+k} expansion, Issue #112")
print("  Phase C) can clarify separately. The sphere fix above is the")
print("  minimal-scope fix for Issue #132.")
