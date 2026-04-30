---
name: Phase 5a closeout — continuous-µ specular sphere reference + Phase 5+ blockers
description: Phase 5a shipped Sanchez 1986 Eq. (A6) reference implementation `compute_K_bc_specular_continuous_mu_sphere`, SymPy 4/4 verifications, literature memo. Production wiring deferred to Phase 5+ — blocked by diagonal-singularity treatment + Sanchez↔ORPHEUS K_ij Jacobian conversion.
type: project
---

# Phase 5a closeout — continuous-µ specular sphere (2026-04-28)

Branch `feature/peierls-specular-bc` (Phase 4 + Phase 5a together).

## What shipped

1. **`derivations/peierls_specular_continuous_mu.py`** — SymPy verification script with 4 checks, all PASS:
   - V1: multi-bounce-factor identity `µ·T(µ) → 1/(2a)` at µ→0 (L'Hôpital limit).
   - V2: M1 sketch ↔ Sanchez Eq. (A6) equivalence. **Important finding**: cross-domain-attacker's M1 sketch had `µ` in the WRONG place (numerator of MB factor). Corrected form drops µ from the MB factor; `F_out` carries the `µ_*^{-1}` Jacobian.
   - V3: vacuum reduction `α → 0` zeros the BC kernel cleanly.
   - V4: diagonal-singularity analysis. The Sanchez Eq. (A6) integrand at the diagonal has a non-integrable `1/µ²` singularity at the surface diagonal (ρ = a) and an integrable `1/µ` at interior diagonals. **This is the Phase 5 production blocker** — Sanchez does not specify a numerical method.

2. **Reference implementation** at `compute_K_bc_specular_continuous_mu_sphere` (peierls_geometry.py) — implements Sanchez Eq. (A6) verbatim with α=1, β=0, ω₁=0. Per-(i, j) integration; uses dimensionless cosines `µ_-` and `µ_*` correctly per the paper p. 342. Working in optical units (ρ = σr, a = σR).

3. **Sanchez 1986 literature memo** at `.claude/agent-memory/literature-researcher/phase5_sanchez_1986_sphere_specular.md` — full extraction of Eqs. (A1)–(A7), notation map to ORPHEUS, citation graph (Pomraning-Siewert 1982, Sanchez-McCormick 1982). Sanchez's BC parametrisation supports both specular (`α`) and diffuse (`β`, `χ(µ)`) reflection; ORPHEUS uses `α=1, β=0`.

4. **Dispatch wiring**: `closure="specular_continuous_mu"` registered in `_build_full_K_per_group` but raises `NotImplementedError` with a clear "Phase 5 research prototype, production wiring pending" message. The matrix-Galerkin form `closure="specular_multibounce"` remains the production path.

## What Phase 5a did NOT ship

- **Production closure wiring**. The reference implementation `compute_K_bc_specular_continuous_mu_sphere` produces Sanchez Eq. (A6) directly, but the result is in Sanchez's normalisation conventions (`g_h(ρ' → ρ)` is the Green's function with `4π ρ'²` surface-area baked in via Sanchez Eq. (2)). ORPHEUS's `K_ij` discretisation uses a different convention: `Σ_t·φ_i = Σ_j K_ij · q_j` with explicit `rv = 4π r²` radial-volume weights and `r_wts` quadrature weights. **Bridging the two requires a Jacobian conversion that's not documented in Sanchez 1986** — Phase 5+ work to derive empirically (rank-1 cross-check against `closure="white_hebert"`).

- **Diagonal-singularity treatment**. The Sanchez kernel has `1/µ²` non-integrable singularity at the surface diagonal. ORPHEUS's existing `build_volume_kernel` uses adaptive quadrature for analogous (but weaker, log) singularities in `E_1`. Phase 5+ scope: adaptive µ-quadrature, singularity subtraction, or change-of-variables.

- **Tests**. Without production wiring there's nothing to gate. Phase 5+ tests would target rank-1 ≡ white_hebert algebraic identity, MC ground-truth k_inf match, no-overshoot-at-high-Q regression.

- **Slab companion**. Hébert §3.8.3 has the slab continuous-µ form. Slab matrix form already converges, so this is verification-quality. Punted to Phase 5+.

- **Cylinder**. Knyazev `Ki_(2+k)` does NOT compose with `f_mb(α, θ_p)` cleanly — cross-domain-attacker analysis. Open as a separate Phase 5++ research direction.

## Phase 5+ blockers (open work)

1. **Sanchez ↔ ORPHEUS K_ij Jacobian conversion** — derive the conversion factor `α(r_i, r_j)` such that `K_ij = α · g_h(r_j → r_i)` matches ORPHEUS's discrete Nyström convention. Approach: empirical rank-1 cross-check (Sanchez at single-µ-node should reduce to Hebert white BC; the ratio reveals α).

2. **Diagonal-singularity treatment** — the integrand at ρ' = ρ has a `1/µ²` (surface) or `1/µ` (interior) singularity. Implementation options:
   - Adaptive Gauss-Kronrod quadrature with explicit subdivision near µ=0
   - Singularity subtraction (analog of ORPHEUS's existing E_1 subtraction in build_volume_kernel)
   - Change-of-variables (e.g., `µ = u²` to absorb `1/µ`)
   - Gauss-Jacobi rule weighted by `µ`
   The choice affects production speed; Sanchez gives no guidance.

3. **Multi-region sphere** — Sanchez Eq. (A6) cosh closed forms require homogeneous σ. For multi-region, the spherical-from-slab reduction (Sanchez Eq. 5) breaks. Two paths:
   - Outer-shell-only with Eq. (A7) coupling for inner regions
   - Per-µ matrix inversion in `T(µ)` (each bounce traverses a different annulus mix)

4. **Slab Phase 5** — Hébert §3.8.3 form with multi-region cleanly via `τ_total/µ`. Low effort; ship as verification cross-check against Phase 4 slab MB.

5. **Cylinder Phase 5++** — needs new derivation for the joint (α, θ_p) integration; Knyazev does not compose with f_mb.

## What this Phase 5a result MEANS for users

- **Production users today**: keep using `closure="specular_multibounce"` for sphere/cyl/slab. The Phase 4 matrix-Galerkin form is shipped, tested, and passes 9 regression tests. For sphere/cyl thin cells, stay at N ∈ {1,2,3} (UserWarning at N≥4).

- **Research users**: `compute_K_bc_specular_continuous_mu_sphere` is callable directly as a research-grade reference. The kernel matches Sanchez Eq. (A6) verbatim. Use for:
  - Cross-checking Phase 4 matrix-form K_bc at low N (rank-1 algebraic identity)
  - Studying the structural difference between matrix-Galerkin projection and continuous-µ kernel
  - Pomraning-Siewert 1982 cross-check truth source for ω₁=0 case

- **Sphinx**: Phase 4 docstrings already point at Phase 5 as the "proper fix" for the matrix-Galerkin pathology. Phase 5a doesn't change the user-facing closure list (specular_continuous_mu raises NotImplementedError until production wiring lands).

## Files touched (Phase 5a only)

- `orpheus/derivations/peierls_geometry.py` (+~280 LoC):
  - `_chord_tau_mu_sphere` helper (multi-region antipodal chord τ(µ))
  - `compute_K_bc_specular_continuous_mu_sphere` reference implementation
  - `closure="specular_continuous_mu"` dispatch (raises NotImplementedError)
  - Updated dispatch help message
- `derivations/peierls_specular_continuous_mu.py` (NEW, ~330 LoC) — SymPy verification with 4 checks
- `.claude/agent-memory/literature-researcher/phase5_sanchez_1986_sphere_specular.md` (NEW)
- `.claude/agent-memory/literature-researcher/MEMORY.md` (index entry)
- `.claude/agent-memory/numerics-investigator/specular_continuous_mu_phase5a_closeout.md` (THIS FILE)
- `.claude/agent-memory/numerics-investigator/MEMORY.md` (index entry)
- `.claude/plans/flickering-wibbling-hennessy.md` (the approved plan that drove Phase 5a)

## Lessons learned

- **Cross-domain-attacker's M1 sketch had a µ-weight bug**. Without the SymPy V2 verification we'd have implemented a kernel wrong-by-a-factor-of-µ — exactly the "plausible substitution error" the QA agent guards against. The SymPy script is a permanent regression for any future Phase 5 reformulation.
- **Diagonal singularity in the Sanchez kernel was not flagged by the literature pull**. The agent's extract correctly transcribed Eq. (A6) but didn't analyse the µ→0 behaviour. SymPy V4 surfaced it.
- **Closed-form Green's functions don't trivially port to discrete Nyström**. Sanchez 1986 is a theoretical reduction; numerical implementation requires its own quadrature design. Future literature pulls for analogous theoretical papers should explicitly ask "what's the implementation gap?".
- **Plan estimates were optimistic by ~1 session**. Best-case "2 sessions for sphere" assumed a clean ORPHEUS-mappable kernel; the Jacobian conversion gap was not visible until production wiring was attempted.
