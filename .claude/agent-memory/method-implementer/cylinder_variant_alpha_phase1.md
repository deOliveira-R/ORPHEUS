---
name: Cylinder Variant α Phase-1 closeout
description: Phase-1 standalone cylinder Variant α Green's function reference — files, commits, accuracy floors, and Phase-2 unification benchmark
type: project
---

# Cylinder Variant α Phase-1 closeout (2026-05-02)

Phase 1 of `peierls-greens-cylinder-and-2bc.md` SHIPPED. Cylinder
Variant α Green's function reference (1G + multi-group homogeneous,
α ∈ [0, 1] specular) lands as a parallel, standalone implementation
mirroring sphere structure end-to-end.

**Branch**: `feature/peierls-greens-cylinder` off `main`.

## Deliverables (all green)

### Files created

- **Branch-1 SymPy**:
  `orpheus/derivations/continuous/peierls/origins/specular/greens_function_cylinder.py`
  (~395 LoC) — V_α1_cyl + V_α1_cyl.geometry + V_α2_cyl + V_α3_cyl
  derivation functions.
- **Branch-2 production**:
  `orpheus/derivations/continuous/peierls/greens_function_cylinder.py`
  (~554 LoC) — `solve_greens_function_cylinder` (1G),
  `solve_greens_function_cylinder_mg` (MG), `_apply_operator_cylinder`
  helper, `_first_leg_2d_chord`, `_impact_parameter`,
  `_bounce_period_2d_chord`, `_scalar_flux_from_psi`,
  `CylinderGreensResult`, `CylinderGreensMGResult` dataclasses.
- **Symbolic test gate**:
  `tests/derivations/test_peierls_greens_function_cylinder_symbolic.py`
  (~210 LoC, 9 foundation-tagged tests).
- **Numerical L1 gate**:
  `tests/derivations/test_peierls_greens_function_cylinder_solver.py`
  (~300 LoC, 7 L1-tagged tests).

### Files modified

- `orpheus/derivations/continuous/peierls/origins/specular/__init__.py`
  — re-export the four cylinder `derive_*()` functions.
- `docs/theory/peierls_greens.rst` — appended cylinder stub
  (~210 lines, 6 equation labels). Sphinx build clean.

### Commits

| SHA      | Message |
| -------- | ------- |
| `159d90c` | feat(derivations): cylinder Variant α — Branch-1 SymPy proofs |
| `d1aac41` | feat(derivations): cylinder Variant α — Branch-2 1G+MG solver + L1 tests |
| `c36a922` | docs(peierls): cylinder Variant α stub on peierls_greens.rst |

## Test results — full peierls greens suite

```
69 passed in 163.64s (0:02:43)
```

- Sphere V_α tests: 53 (V_α1/2/3 SymPy, V_α1.numerical, B5 xverif,
  vacuum, PS-1982, MG, MR, Garcia 2021) — green, no regression.
- Cylinder V_α tests: 16 (9 symbolic + 7 numerical L1) — green.

## Accuracy floors captured (Phase-2 unification benchmark)

### α=1 closed-cylinder k_inf exactness (V_α1_cyl numerical)

fuel-A-like XS (R=5, σ_t=0.5, σ_s=0.38, νσ_f=0.025), τ_R = 2.5,
k_inf = 0.20833333...:

| (n_r, n_μ, n_φ, n_t) | k_eff | rel_err to k_inf | iter |
|---|---|---|---|
| (8, 8, 16, 24) | 0.2083333333333330 | 1.87e-15 | 1 |
| (12, 10, 24, 32) | 0.2083333333333338 | 2.00e-15 | 1 |
| (16, 12, 32, 48) | 0.2083333333333344 | 4.93e-15 | 1 |
| (20, 16, 48, 64) | 0.2083333333333335 | 6.66e-16 | 1 |

**Variant α k_eff = k_inf to machine precision (1e-15) at every
quadrature order, in 1 iteration.** This is V_α1_cyl algebraic
identity reproduced numerically.

### α=0 vacuum k_eff convergence-floor (the non-trivial case)

Same fuel-A-like XS, α=0:

| (n_r, n_μ, n_φ, n_t) | k_eff | iter | time |
|---|---|---|---|
| (8, 8, 16, 24) | 0.1204312553 | 51 | 0.9 s |
| (12, 10, 24, 32) | 0.1204531627 | 51 | 2.1 s |
| (16, 12, 32, 48) | 0.1204496559 | 51 | 4.5 s |
| (20, 16, 48, 64) | 0.1204515822 | 51 | 11.8 s |
| (24, 20, 64, 96) | 0.1204515921 | 51 | 28.6 s |

Self-convergence (rel diff to next-finer):

- (8,8,16,24) → (12,10,24,32): 1.82e-4
- (12,10,24,32) → (16,12,32,48): 2.91e-5
- (16,12,32,48) → (20,16,48,64): 1.60e-5
- (20,16,48,64) → (24,20,64,96): **8.24e-8** ← Phase-2 benchmark

The convergence is non-monotone at intermediate orders (typical for
3D Gauss-Legendre on a problem with grazing-ray loci). Above (n_r,
n_μ, n_φ, n_t) ≥ (20, 16, 48, 64) the answer stabilizes at
**k_eff ≈ 0.12045159** to ~8e-8 self-consistency.

### MG accuracy

2G asymmetric scattering, α=1:
- k_eff = k_inf = 1.0 to 1.5e-11 relative in 48 iterations.
- Each group's scalar flux uniform to ~1e-15 std/mean.

MG G=1 reduction matches 1G solver to 1e-9 at α=0.

## Cross-check disclaimer (Issue #129)

Cylinder Variant α and cylinder Nyström vacuum disagree at the same
XS by ~1.4e-4 relative:

- Variant α (3D angle-resolved): k_eff = 0.1204515921 (n=24,20,64,96)
- Nyström vacuum (axial pre-integrated, n=48): k_eff = 0.1204656335

This is **not a bug** — Issue #129 documents the planar-limit
discrepancy of axially pre-integrated cylinder Bickley-Naylor. Both
formulations are correct per their respective models. The
load-bearing L1 cross-check for cylinder Variant α is the
**k_inf-exactness invariant at α=1** (V_α1_cyl algebraic ground,
structurally independent of any other ORPHEUS solver), where
Variant α matches to 1e-15.

## Structural surprises and gotchas

### What was structurally identical to sphere

- The bounce-sum closure algebra is **literally identical**: surface
  fixed-point ψ_surf = α B / (1 − α e^{−Σ_t L_period}). V_α1_cyl
  inherits V_α1's algebraic cancellation pattern (both L_0 and
  L_period drop out) without any cylinder-specific logic.
- The first-leg + B + closed-form-T architecture transfers without
  modification.
- Source profile cubic-spline interpolation works identically.
- Power iteration + Rayleigh-quotient k update works identically.

### What is genuinely cylinder-specific

- **3D phase-space** (r, μ_axial, φ_az) vs sphere's 2D (r, μ).
  Cylinder needs both axial-cosine AND in-plane azimuth — these
  cannot be decoupled because (a) impact parameter b = r·|sin(φ_az)|
  conserves only via φ_az, (b) 3D path scaling 1/√(1−μ²) requires
  μ_axial. This is the load-bearing addition.
- **Two conserved invariants** (b AND μ_axial) instead of sphere's
  one (μ_surf). The bounce-period chord depends on BOTH, and the
  closed-form geometric series is over a 2D parameter family.
- **3D path parametrization** — integrate in 2D arclength s_2D and
  scale attenuation by 1/sin(θ_axis) to recover 3D Σ_t·s_3D. This
  is a deliberate convention that keeps the chord geometry simple
  while preserving 3D physics.
- **Two grazing-ray loci** instead of sphere's one:
  - Axial grazing μ_axial → ±1 (analog of sphere's μ → 0).
  - Tangential in-plane φ_az → ±π/2 (impact parameter → R, no
    bounces — handled by `_bounce_period_2d_chord` returning 0).
  Gauss-Legendre on full ranges avoids endpoints and these are
  numerically benign.
- **Prefactor convention**: cylinder uses 1/π (azimuthal-folding)
  in the existing Nyström primitives. Variant α uses the standard
  1/(4π) for the isotropic source per steradian and then
  integrates over the FULL φ_az ∈ [0, 2π], which is mathematically
  consistent (not a 1/π convention internally).

### What did NOT manifest (but might in future phases)

- No grazing-ray pathology surfaced in numerical practice. The
  open Gauss-Legendre quadrature on (0, R) × [-1, 1] × [0, 2π]
  avoids the endpoints; iteration converges cleanly. This may need
  reconsideration for Phase 1b multi-region cases where impact-
  parameter sweeps cross interior shell boundaries.
- No SymPy-choke encountered. All four `derive_*()` functions
  closed in milliseconds. The cylinder algebra is at the same
  complexity tier as sphere V_α1/2/3 — algebraic, not analytical
  closed-form for V_α2_cyl (same Ki_3 special function as
  Nyström).

## Bugs caught during development

**None.** The Branch-2 implementation landed clean on first run —
V_α1_cyl gave k_inf to 1e-15 in 1 iteration, MG gave k_inf to 1e-11,
α=0 vacuum gave physics-plausible peaked-at-center flux. No L0
errors to log.

This clean landing reflects the **algebra-of-record discipline**
working as intended: Branch-1 SymPy proofs verified the algebra
BEFORE Branch-2 implementation, so the production code only had to
operationalize what was already known correct algebraically.

## Closing recommendation on Phase-2 unification

**The cylinder math felt structurally identical to sphere.** The
bounce-sum closure is literally word-for-word the same algebra; the
F + B + T architecture transfers without modification; the only
difference is the chord formulas and the third angular variable.

The Phase-2 unification (operator-theoretic refactor under
`CurvilinearGeometry`) is therefore **structurally well-supported** —
the proposed S = α·R_chord scattering operator + T = (I − S)^{−1}
resolvent abstraction maps cleanly onto both sphere (rank-1) and
cylinder (rank-1, parametrised by 2D conserved-invariant family).

**Caveat** (per `feedback_unify_after_two_instances.md`): Phase 2
must preserve the captured accuracy floor (~8e-8 at vacuum, 1e-15
at closed-cylinder). If unification regresses cylinder accuracy
beyond ~1e-6 relative, HALT — the abstraction has hidden a relevant
discretization choice. The accuracy floor in this memo is the
quantitative gate.

## GitHub issue follow-ups

- Issue #101 (cylinder Peierls reference, chord-based Ki_1 Nyström)
  — **comment that Variant α offers a parallel reference path**.
  Not closing — separate decision.
- Phase 1b (deferred to subsequent dispatch):
  - Multi-region cylinder k-eigenvalue (mirror of
    `solve_greens_function_sphere_mr`).
  - Multi-region cylinder fixed-source (mirror of
    `solve_greens_function_sphere_mr_fixed_source`).
  - These are in scope for cylinder unification but were
    explicitly out of scope for this dispatch.

## Pointers

- Plan: `.claude/plans/peierls-greens-cylinder-and-2bc.md`
- Sphere reference closeout (template):
  `.claude/agent-memory/numerics-investigator/peierls_greens_phase1_closeout.md`
- Sphere SymPy reference: `orpheus/derivations/continuous/peierls/origins/specular/greens_function.py`
- Cylinder SymPy reference: `orpheus/derivations/continuous/peierls/origins/specular/greens_function_cylinder.py`
- Sphere production: `orpheus/derivations/continuous/peierls/greens_function.py`
- Cylinder production: `orpheus/derivations/continuous/peierls/greens_function_cylinder.py`
- Sphinx page: `docs/theory/peierls_greens.rst` (cylinder section starts at `:label: peierls-greens-cylinder`)
