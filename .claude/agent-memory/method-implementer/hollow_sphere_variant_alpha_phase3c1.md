---
name: Hollow Sphere Variant α Phase 3C-1 closeout
description: Phase-3C-1 hollow sphere Variant α Green's function reference solver shipped clean — rank-2 BIE block resolvent on through-rays + rank-1 outer-only closure with impact-parameter phase-space partition. First curvilinear 2-surface instance.
type: project
---

# Hollow Sphere Variant α — Phase 3C-1 closeout (2026-05-02)

Branch: `feature/peierls-greens-cylinder`. Predecessor commit:
`9647081 docs(peierls): note ERR-035 fix in slab Variant α theory`.

## Files created/modified

**Created:**
- `orpheus/derivations/continuous/peierls_greens_function/origins/specular/greens_function_hollow_sphere.py`
  (Branch-1 SymPy V_α1_hollow_sph / V_α2_hollow_sph / V_α3_hollow_sph)
- `orpheus/derivations/continuous/peierls_greens_function/greens_function_hollow_sphere.py`
  (Branch-2 production: `solve_greens_function_hollow_sphere`,
  `solve_greens_function_hollow_sphere_mg`)
- `tests/derivations/test_peierls_greens_function_hollow_sphere_symbolic.py`
  (18 foundation-tagged SymPy gates)
- `tests/derivations/test_peierls_greens_function_hollow_sphere_solver.py`
  (11 L1-tagged solver gates)
- `.claude/agent-memory/method-implementer/hollow_sphere_variant_alpha_phase3c1.md`
  (this file)

**Modified:**
- `orpheus/derivations/continuous/peierls_greens_function/origins/specular/__init__.py`
  (export the three new `derive_*` functions)
- `docs/theory/peierls_greens.rst`
  (Phase-3C-1 Sphinx stub with 4 labelled equations + archivist TODO)
- `docs/verification/matrix.rst` (Nexus auto-regenerated on build)

## Commits

1. `c0486c7 feat(peierls): Phase 3C-1 hollow sphere SymPy V_α derivations + foundation gates`
2. `47b745c feat(peierls): Phase 3C-1 hollow sphere Variant α reference solver (1G + MG)`
3. `d2a9b79 docs(peierls): hollow sphere Variant α stub (Phase 3C-1)`

## Test results

**Acceptance gate (full)**: 145/145 passed in 211 s (3.5 min).
- Pre-Phase-3C-1: 116 passed.
- Phase-3C-1 additions: 18 symbolic + 11 solver = 29 new tests.
- Net: 116 + 29 = 145, all pass.

Sphinx -W build: clean (exit 0). 4 pre-existing SN MMS soft-warnings
unrelated to this phase persist.

## Accuracy floors

### V_α1 closed-shell exactness (load-bearing structural composability)

`α_in = α_out = 1`, fuel-A-like XS, R_in=1.0, R_out=3.0:
- Converged in **1 iteration** from constant initial guess.
- `k_eff = k_inf` to **rel err = 4e-16** (machine precision).
- φ spread (std/mean) = **1e-16** (uniform flux).

This proves the impact-parameter partition (b > R_in outer-only AND
b ≤ R_in through-ray) composes cleanly with the rank-2 closure.

### R_in → 0 solid-sphere limit agreement

`α_in = α_out = 0`, vacuum-vacuum, R_in = 1e-3·R_out:
- (n_r=24, n_mu=24, n_traj_quad=48): rel diff vs solid sphere = **9.2e-10**.
- (n_r=32, n_mu=32, n_traj_quad=64): rel diff = **1.9e-9**.

**Plan target was 1e-3; achieved ≤ 1e-9** — six orders of magnitude
better. The phase-space partition reduces correctly: as R_in shrinks
toward zero, the through-ray subset (b ≤ R_in) approaches measure-
zero and the outer-only branch dominates, recovering solid-sphere
behaviour bit-equivalently up to the inner-cavity perturbation.

### Asymmetric flux-shape sanity

- Reflective inner / vacuum outer (`α_in=1, α_out=0`):
  k_eff = 0.0567 < k_inf = 0.208. φ peaks near inner wall;
  inner/outer ratio = 3.18.
- Vacuum inner / reflective outer (`α_in=0, α_out=1`):
  k_eff = 0.174 < k_inf. φ monotone increasing from inner (cavity
  absorber) to outer (reflective). The cavity-absorber physical
  interpretation holds.
- Symmetric reflective `α=0.5`: k_eff = 0.089 < k_inf, well-behaved.

### Multi-group asymmetric scattering at closed shell

2G with asymmetric Σ_s, `α_in=α_out=1`:
- `k_eff = 1.0000000000151812e+00`, rel err vs k_inf = **1.5e-11**.
- Per-group φ spread = 1.6e-16, 3.6e-16 (machine-precision uniform).

ERR-002 transpose drift detector: PASSED.

### Convergence floor (research-grade)

Reflective inner / vacuum outer at quadrature ladder
(16,16,32) → (24,24,48) → (32,32,64) → (48,48,96):
finest pair consistency = `< 1e-3`. Captures convergence behaviour;
the precise asymptote is the eigenvalue of the rank-2 prototype at
this BC corner.

## Sphere/cylinder/slab regression status

**Bit-equal accuracy floors preserved**:
- Sphere closed BC k_eff = 2.0833333333333340e-01 (was 2.083...340e-01).
- Cylinder closed BC k_eff = 2.0833333333333423e-01 (was identical).
- Slab closed BC k_eff = 2.0833333333333467e-01 (was identical).
- Slab asym closed BC k_eff = 2.0833333333333467e-01 (was identical).

All within ≤ 7e-16 absolute floor. Zero regression on existing 116
tests (all pass at original tolerances).

## Honest verdict on the rank-2 + impact-parameter abstraction

### Did the partition split (b > R_in vs b ≤ R_in) compose cleanly with the rank-2 closure?

**YES.** V_α1 closed-shell at machine precision in 1 iteration is
the structural composability check, and it passed. SymPy proofs
verified algebraically (V_α1_hollow_sph proves both branches
independently produce ψ = q/Σ_t at α_in=α_out=1).

The V_α1 algebra is the SAME for both branches: rank-1 outer-only
collapses via `(1-x)/(1-x) = 1` for the bounce-period chord
(identical to V_α1 sphere); rank-2 through-ray collapses via
`(1-e^{-τ})(1+e^{-τ})/(1-e^{-2τ}) = 1` (identical to V_α1
slab-asym). The two branches share the SAME constant ψ = q/Σ_t at
the boundary b = R_in, so phase-space continuity is preserved.

### Is there a regularization concern at b → R_in (tangent rays)?

**NO concrete pathology observed at production grids.** At b → R_in:
- Through-ray branch: τ_step = Σ_t · (sqrt(R_out²-b²) - sqrt(R_in²-b²))
  → sqrt(R_out²-R_in²) − 0 = finite, nonzero, well-conditioned.
- Outer-only branch: same boundary, same chord (sqrt(R_out²-R_in²)
  is the half-period), same `T = 1/(1 - α e^{-2τ})` — finite.
- The two branches **converge to the same numerical value at the
  boundary** because the phase-space subset of measure zero (single
  tangent ray) doesn't break either closure.

Pathological corner: `α_out = 1, b → R_out` would give `τ_period →
0`, `T → 1/(1-α) → ∞`. This is the same grazing-ray issue the solid
sphere has at α=1, μ→0. Mitigated by Gauss-Legendre interior
quadrature avoiding the exact endpoint b = R_out.

### Will the abstraction extend to annulus (Phase 3C-2) without modification?

**Mostly yes; minor adaptation required for cylinder's 3D angular
space.** The structural elements that transfer cleanly:

1. Impact-parameter phase-space partition by inner-radius threshold
   (b > R_in vs b ≤ R_in). Same logic.
2. Rank-2 monodromy on through-rays with single-transit shell-
   traversal optical depth — same algebraic form, just with
   cylinder-specific chord arithmetic
   (`tau_step = Σ_t · sqrt((R_out² - b²) - (R_in² - b²)) / sqrt(1 -
   μ_axial²)` or similar — to be confirmed in Phase 3C-2 derivation).
3. Outer-only rank-1 closure for b > R_in — reuses cylinder solid
   solver primitives.
4. Variant α core primitives (`compute_resolvent_T_rank2`,
   `apply_variant_alpha_closure_rank2`) — geometry-agnostic at the
   operator-symbol level. NO modification needed to the core.

The new piece for annulus: cylinder's angular discretisation is
(μ_axial, φ_az) instead of single μ, so the impact parameter
formula generalises to `b² = r² · (1 − μ_axial² − cos²(φ_az − φ_r))`
or similar. This is a chord-algebra change, not a closure-algebra
change. The rank-2 + impact-parameter abstraction is robust across
this generalisation.

**Architectural verdict**: the rank-2 framework lifted from slab to
hollow sphere with NO modification to `variant_alpha_core.py`. The
core primitives are truly geometry-agnostic. Phase 3C-2 annulus
will require only:
- A new module `greens_function_annulus.py` mirroring
  `greens_function_hollow_sphere.py` structure.
- Cylinder-specific chord arithmetic for impact parameter and
  shell-traversal length.
- New SymPy module `origins/specular/greens_function_annulus.py`.

## L0 bugs caught during development

**NONE.** The hollow sphere prototype shipped clean on first try
across all 29 new tests. The earlier Phase 3A → 3B work caught
ERR-034 (trajectory missing μ factor) and ERR-035 (Phase-3A heuristic
closure mismatch); those FIXES carried forward to Phase 3C-1 because
the hollow-sphere `_apply_operator_hollow_sphere` mirrors the
correct trajectory parametrisation `r_traj_sq = r² - 2r·μ·s + s²`
from sphere/slab-asymmetric.

**No new ERR-NNN entries logged.**

## Cross-domain frame analysis

The hollow sphere is the **first curvilinear 2-surface validation**
of the rank-2 BIE block resolvent frame predicted by
`.claude/agent-memory/cross-domain-attacker/variant_alpha_2surface_bie_frame.md`.

Frame match: STRONG. The slab-asymmetric template lifted to hollow
sphere with only chord-algebra changes; the closure machinery
(`apply_variant_alpha_closure_rank2`, `compute_resolvent_T_rank2`)
is **byte-equal-shared** between slab and hollow-sphere prototypes.

The two new geometric features beyond slab — phase-space partition
by impact parameter, and curvilinear shell-traversal chord — both
compose cleanly with the rank-2 closure. The frame's structural
prediction (rank-1 → rank-2 generalises 1-surface → 2-surface
geometry, regardless of curvature) holds across slab AND hollow
sphere. Phase 3C-2 annulus will be the third instance.

## Lessons for future implementations

1. **Foundation-tagged SymPy gates first.** Writing the SymPy
   derivation BEFORE the production solver caught the algebra
   composability concern (does outer-only AND through-ray
   independently produce q/Σ_t at α=1?) symbolically before any
   numerical work. Both branches passed V_α1 symbolically,
   guaranteeing the numerical V_α1 closed-shell test would also
   pass — and it did, machine precision in 1 iter.

2. **Sign-of-μ × b-vs-R_in case analysis matters.** The first-leg
   backward chord depends on FOUR cases:
   - μ > 0, b > R_in: outer-only, chord to outer surface.
   - μ > 0, b ≤ R_in: through-ray, chord to inner surface.
   - μ < 0, b > R_in: outer-only, chord to outer surface (different
     branch of the discriminant root).
   - μ < 0, b ≤ R_in: through-ray, chord to outer surface.
   Worked through symbolically (sign of `r·μ ± sqrt(disc)`) before
   coding — saved the pattern of "chord goes negative if you pick
   the wrong root" debugging cascade.

3. **The B_in / B_out parametrisation needed care.** The shell-
   traversal chord at conserved b is the SAME geometric chord
   traversed in two directions. `B_out` (inner→outer) parametrises
   r²(u) = b² + (sqrt_disc_in + u)² (outward); `B_in` (outer→inner)
   parametrises r²(u) = b² + (sqrt_disc_out − u)² (inward). The
   slab-asymmetric mapping `α_left ↔ α_in, α_right ↔ α_out,
   B_LR ↔ B_out, B_RL ↔ B_in` was the cleanest and worked first try.

4. **CubicSpline domain extension to the shell**. The radial nodes
   live on `[R_in, R_out]` (open quadrature avoids endpoints which
   are surfaces). The CubicSpline `extrapolate=True` allows
   evaluation slightly outside `[R_in, R_out]` for trajectory points
   near the surface; clipping to `[R_in², R_out²]` on `r_traj_sq`
   keeps physical interpretation safe.

## Open issues / deferred work

- **Multi-region hollow sphere**: deferred. The single-region
  prototype is sufficient to validate the rank-2 framework on
  curvilinear 2-surface geometry. Multi-region hollow can follow
  the sphere `_apply_operator_mr` segmentation pattern when needed.
- **Phase 3C-2 annulus**: next priority per plan. Will mirror this
  structure with cylinder-specific chord arithmetic.
- **Issue #132 Class B catastrophe at hollow sphere**: deferred.
  The Phase-3 work has so far focused on the rank-2 framework
  validation; Issue #132 reproducer for hollow sphere can be added
  as a follow-on test once the multi-region-hollow capability
  exists.

## Provenance

- Plan: `.claude/plans/peierls-greens-cylinder-and-2bc.md` Phase 3C-1.
- Cross-domain frame: `.claude/agent-memory/cross-domain-attacker/variant_alpha_2surface_bie_frame.md`.
- Phase 3B closeout (predecessor): `slab_asymmetric_variant_alpha_phase3b.md`.
- ERR-035 fix (predecessor): `err035_phase3a_delegation_fix.md`.
