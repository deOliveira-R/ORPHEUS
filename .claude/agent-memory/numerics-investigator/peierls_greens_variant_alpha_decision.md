---
name: Plan 2 B2 — Variant α architectural decision
description: Variant α (integral-operator power iteration on Sanchez 1986 Eq. A6) chosen as Plan 2 prototype path. Variants β/γ rejected with reasoning. Variant δ (Garcia stable P_N) deferred to B5 cross-check. Critical structural note on what makes Variant α genuinely different from Phase 5 retreat.
type: project
---

# Plan 2 B2 — Variant α architectural decision

**Date**: 2026-05-01 (Claude Opus 4.7, 1M context).
**Status**: APPROVED — proceed to B3 (SymPy derivation).
**Predecessors**: B1 literature pull
(`.claude/agent-memory/literature-researcher/peierls_greens_function_lit.md`),
Phase 5 retreat
(`.claude/agent-memory/numerics-investigator/_archive/specular_continuous_mu_phase5_retreat.md`).

## TL;DR

Adopt **Variant α** — phase-space (r, µ) power iteration on the
**angle-resolved** Green's function with bounces summed analytically
(Sanchez 1986 Eq. A5 specular leg, α=1). Variant α evaluates kernels
along **bouncing characteristics** (a sequence of 1-D arc integrals),
not by Nyström sampling of the angle-integrated kernel
`g_α(ρ' → ρ)` that killed Phase 5. The B5 cross-check anchor is the
existing rank-1 algebraic identity T_00^sphere = P_ss^sphere
(specular_multibounce ≡ white_hebert at rank-1).

## Decision

**Path**: Variant α (integral-operator power iteration along
bouncing characteristics).

**Rationale**:

1. **Reuses ORPHEUS primitives**. Sanchez Eq. (A5) `t_h(r', Ω' → r,
   Ω)` for α=1 reduces to a sum over bounces of trajectory integrals
   along a chord of length 2Rµ_surf per bounce. The geometric series
   closes to T(µ_surf) = 1/(1 − e^{−2Σ_t R µ_surf}). All ingredients
   (chord lengths, optical depths, T(µ)) already implemented in
   `compute_T_specular_sphere` and the geometry primitives.
2. **No closure needed at the operator level**. BC absorbed into the
   kernel via Sanchez Eq. (A1) `t = t̄ + t_h`. No K_bc separation,
   no rank-N gating, no `(1−P_ss)^{−1}` Hébert correction at the
   operator level.
3. **Sidesteps the Phase 5 hypersingularity**. The angle-integrated
   kernel `g_α(ρ' → ρ)` is hypersingular at ρ = ρ' = R (1/µ²
   non-integrable spike). Variant α never assembles `g_α` — instead
   it evaluates the angle-RESOLVED Green's function `t̃(r' → r, µ)`
   along characteristics, which is pointwise finite (just an
   exponential along a chord). Singularity emerges only at the final
   step `φ(r) = ∫ ψ(r, µ) 2π dµ`, which we handle with µ-weighted
   Gauss-Jacobi quadrature (or u² substitution).

**Variant β rejected**: surface-IE / BEM-style reformulation with
volume + surface unknowns coupled. Substantially more invasive; no
literature precedent for sphere-specular k_eigenvalue; rank-1
verification anchor unclear.

**Variant γ rejected**: Case-Zweifel singular eigenfunction expansion
is dominated by the modern Garcia 2020/2021 stable-P_N
re-summation (Variant δ). Requires Chandrasekhar-H projection for
continuous-spectrum component — heavy machinery for a prototype.

**Variant δ deferred to B5**: Garcia 2020/2021 JCP stable P_N is the
state-of-the-art numerical reference in the published literature for
homogeneous-sphere isotropic-scattering with reflective BC. Use as
external truth source for B5 cross-verification once user pulls
paywalled PDFs. NOT a primary path because it inherits a P_N
truncation closure (analogous to ORPHEUS's rank-N gating).

## What makes Variant α genuinely different from Phase 5

The Phase 5 retreat documented three failure modes:

1. **Nyström sampling of K_ij = g_α(ρ_j → ρ_i)** — diverges at the
   surface diagonal i = j corresponding to ρ = R.
2. **Galerkin double-integration** `A_ij = ∫∫ L_i L_j g_α dρ' dρ` —
   diverges via log(Q_µ) inside the µ-integration of g_α at coincident
   (ρ, ρ') Lagrange-Gauss collocation pairs.
3. **Bounce-resolved M2 expansion** sampled at K_max=0 — exposed the
   diagonal singularity persists EVEN without the multi-bounce factor
   T(µ).

All three failures share the same root: assembling the
**angle-integrated** kernel `g_α(ρ' → ρ)` produces a function with
a non-integrable spike at the surface diagonal, and that spike is
hit by any quadrature-on-quadrature scheme.

**Variant α never assembles `g_α`.** The iteration variable is the
angle-resolved ψ(r_i, µ_q) on a 2-D phase-space grid. For each grid
point (r_i, µ_q):

```
ψ^(n+1)(r_i, µ_q) = (Σ_s + νΣ_f/k) · ∫_traj e^{−Σ_t s} φ^(n)(r(s)) ds
                  + e^{−Σ_t L_first} · T(µ_surf(r_i, µ_q))
                    · ∫_period e^{−Σ_t s} φ^(n)(r(s)) ds
```

The trajectory integrals are 1-D arc integrals along a single
characteristic — each integrand is bounded (exponential of negative
times a bounded scalar flux). The geometric T(µ_surf) is a finite
prefactor for each (r_i, µ_q) with µ_surf > 0.

Scalar flux extraction at the end:
```
φ^(n+1)(r_i) = 2π ∫_{−1}^{1} ψ^(n+1)(r_i, µ) dµ
```

The µ-integration encounters the µ → 0 singularity of T(µ) ~
1/(2Σ_t R µ), which is integrable (logarithmic) but requires care
in quadrature. **Resolution**: use Gauss-Jacobi quadrature weighted
by µ on [0, 1], OR substitute u = √µ which absorbs the µ measure.

## Critical novelty disclosure

**No published reference uses integral-operator power iteration on
the angle-resolved Green's function for sphere homogeneous specular
k-eigenvalue.** Sanchez 1986 stops at the analytical reduction;
Pomraning-Siewert 1982 stops at the analytical reduction; Lewis-Miller
1984 doesn't cover this approach; Garcia 2020/2021 uses spectral P_N.

This is **genuinely original**. The B3 SymPy step is therefore
**load-bearing** — if the operator action on a polynomial trial
function turns out to be divergent, Variant α inherits Phase 5's
structural pathology despite the reformulation.

**Risk**: low but non-zero. The functional-analysis intuition
(integration over `dr'·dΩ'` on the trial function smooths the
diagonal singularity to a measure-zero set) is standard, but ORPHEUS
has no literature reference where this exact transformation has been
carried out for sphere specular.

## Plan B3 SymPy verifications (revised from B1 memo)

The literature memo identified six SymPy tasks. We refine them to
**three load-bearing verifications**:

- **V_α1**: rank-1 isotropic-trial eigenvalue identity. For closed
  homogeneous sphere with specular BC, the rank-1 isotropic mode
  satisfies (1 − ω₀) · const = 0, i.e., k = k_inf = νΣ_f/Σ_a. SymPy
  proves the operator action on constant trial gives ω₀ · const.
- **V_α2**: T_00^sphere = P_ss^sphere algebraic identity. SymPy proves
  the rank-1 reduction of `compute_T_specular_sphere` equals
  `compute_P_ss_sphere` for homogeneous sphere. Cross-check anchor
  for B5 (rank-1 Variant α ≡ rank-1 specular_multibounce ≡ white_hebert).
- **V_α3**: vacuum-BC reduction. SymPy proves Sanchez Eq. (A6)
  leading factor 2α → 0 at α=0, so g_h vanishes and g_α = ḡ_2 (vacuum
  kernel from Sanchez Eq. 5). Sanity check that Variant α reduces
  cleanly to the existing vacuum sphere reference.

Other five items from B1 memo (operator on polynomial trial, spectral
radius / Perron-Frobenius, k_eff Rayleigh quotient, multi-region
piecewise τ(µ)):

- Operator action on polynomial trial → re-cast as V_α1 (constant
  trial is the simplest polynomial; if rank-1 eigenmode works, polynomial
  trials of higher degree are a B4 prototype concern, not a B3 SymPy
  concern).
- Spectral radius / Perron-Frobenius → trivial for positive kernel
  (every factor in t̃ is positive on the integration domain). No
  SymPy verification needed; documented in the V_α1 docstring.
- k_eff Rayleigh quotient → standard fission-source iteration; no
  Variant-α-specific algebra to verify.
- Multi-region → out of B3 scope.

## File index

- This memo: `.claude/agent-memory/numerics-investigator/peierls_greens_variant_alpha_decision.md`
- B1 literature: `.claude/agent-memory/literature-researcher/peierls_greens_function_lit.md`
- Phase 5 retreat: `.claude/agent-memory/numerics-investigator/_archive/specular_continuous_mu_phase5_retreat.md`
- Sanchez 1986 memo: `.claude/agent-memory/literature-researcher/phase5_sanchez_1986_sphere_specular.md`
- Plan 2: `.claude/plans/peierls-greens-function-approach.md`
- Existing kernel-level SymPy: `orpheus/derivations/continuous/peierls/origins/specular/continuous_mu.py` (V1–V4: kernel form identities)
- B3 SymPy target: `orpheus/derivations/continuous/peierls/origins/specular/greens_function.py` (V_α1–V_α3: operator-level identities)
- B3 test gate: `tests/derivations/test_peierls_greens_function_symbolic.py`
