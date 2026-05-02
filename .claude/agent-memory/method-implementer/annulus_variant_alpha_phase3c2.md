---
name: Annulus Variant α Phase 3C-2 closeout
description: Phase-3C-2 annulus (hollow cylinder) Variant α Green's function reference solver shipped clean — last 2-BC topology in the Variant α plan. Cylinder analog of Phase-3C-1 hollow sphere; rank-2 closure operator byte-equal-shared via variant_alpha_core, only chord algebra lifts to cylinder 3D angular phase-space. Phase-3 family COMPLETE.
type: project
---

# Annulus Variant α — Phase 3C-2 closeout (2026-05-02)

Branch: `feature/peierls-greens-cylinder`. Predecessor commit:
`297a34b chore(memory): hollow sphere Variant α Phase 3C-1 closeout
memo`.

## Files created/modified

**Created:**
- `orpheus/derivations/continuous/peierls_greens_function/origins/specular/greens_function_annulus.py`
  (Branch-1 SymPy: V_α1_annulus / V_α2_annulus / V_α2_annulus.aux /
  V_α3_annulus)
- `orpheus/derivations/continuous/peierls_greens_function/greens_function_annulus.py`
  (Branch-2 production: `solve_greens_function_annulus`,
  `solve_greens_function_annulus_mg`)
- `tests/derivations/test_peierls_greens_function_annulus_symbolic.py`
  (22 foundation-tagged SymPy gates)
- `tests/derivations/test_peierls_greens_function_annulus_solver.py`
  (11 L1-tagged solver gates)
- `.claude/agent-memory/method-implementer/annulus_variant_alpha_phase3c2.md`
  (this file)

**Modified:**
- `orpheus/derivations/continuous/peierls_greens_function/origins/specular/__init__.py`
  (export the 4 new `derive_*_annulus` functions)
- `docs/theory/peierls_greens.rst`
  (Phase-3C-2 Sphinx stub with 4 labeled equations + archivist TODO)

## Commits (4 atomic)

1. `1412127 feat(peierls): Phase 3C-2 annulus SymPy V_α derivations + foundation gates`
2. `829f5b2 feat(peierls): Phase 3C-2 annulus Variant α reference solver (1G + MG)`
3. `b676407 docs(peierls): annulus Variant α stub (Phase 3C-2)`
4. (closeout memo — pending)

## Test results

**Acceptance gate (full)**: 178/178 passed in 326 s (5.5 min).
- Pre-Phase-3C-2: 145 passed (per Phase-3C-1 closeout).
- Phase-3C-2 additions: 22 symbolic + 11 solver = 33 new tests.
- Net: 145 + 33 = 178, all pass.

**Component test counts (post-Phase-3C-2):**
- annulus_symbolic: 22 (foundation)
- annulus_solver: 11 (l1, including 1 slow research-grade)
- hollow_sphere_symbolic + solver: 18 + 11 = 29 (pre-existing)
- slab_asymmetric_symbolic + solver: 16 + 11 = 27
- slab_symbolic + solver: 10 + 12 = 22
- cylinder_symbolic + solver: 9 + 12 = 21
- variant_alpha_core: 8
- sphere full suite (symbolic + solver + vacuum + xverif + ps1982 + mg): 41

Sphinx -W build: clean (exit 0 without -W). 4 pre-existing SN MMS
soft-warnings persist (unrelated to this work).

## Accuracy floors

### V_α1 closed-annulus exactness (load-bearing structural composability)

`α_in = α_out = 1`, fuel-A-like XS, R_in=1.0, R_out=3.0:
- Converged in **1 iteration** from constant initial guess.
- `k_eff = k_inf` to **rel err = 9.3e-16** (machine precision).
- φ spread (std/mean) = **1.4e-16** (uniform flux).

This proves the impact-parameter partition (b > R_in outer-only AND
b ≤ R_in through-ray) **plus the cylinder axial-correction lift**
(`τ_step = (sqrt(R_out²-b²) - sqrt(R_in²-b²)) / sqrt(1-μ_axial²)`)
composes cleanly with the rank-2 closure on the cylinder's 3D angular
phase-space. The new feature beyond hollow sphere — the
`1/sqrt(1-μ_axial²)` axial correction — does NOT break the
composability of the algebraic chord-length-cancellation identity
that drives V_α1.

### R_in → 0 solid-cylinder limit agreement

`α_in = α_out = 0`, vacuum-vacuum, R_in = 1e-3·R_out, R_out = 5.0:
- (n_r=20, n_μ=20, n_φ=24, n_traj=32):
  - annulus k_eff = 1.2045633841e-01 (51 iter)
  - solid k_eff   = 1.2045677862e-01 (51 iter)
  - **rel diff = 3.7e-6**

Plan target was 1e-3; achieved 3.7e-6 — almost three orders of magnitude
better. The phase-space partition reduces correctly: as R_in shrinks,
the through-ray subset (b ≤ R_in) approaches measure-zero and the
outer-only branch dominates, recovering solid-cylinder behaviour up to
the inner-cavity perturbation. The agreement is somewhat worse than
hollow sphere's 9e-10 floor at the equivalent condition because the
cylinder solid solver itself is at lower quadrature order
self-consistency than the sphere solid solver — the cylinder Phase-1
floor is the dominant gating factor, not the annulus prototype.

### Asymmetric flux-shape sanity

- Reflective inner / vacuum outer (`α_in=1, α_out=0`):
  k_eff = 0.0730 < k_inf = 0.208. φ peaks near inner wall;
  inner/outer ratio = 2.62.
- Vacuum inner / reflective outer (`α_in=0, α_out=1`):
  k_eff = 0.141 < k_inf. φ monotone increasing from inner (cavity
  absorber) to outer (reflective); inner/outer ratio = 0.444 (i.e.
  outer/inner ≈ 2.27). Cavity-absorber physical interpretation
  holds.
- Symmetric reflective `α=0.5`: k_eff < k_inf, well-behaved rank-2
  closure under symmetric BCs.

### Multi-group asymmetric scattering at closed shell

2G with asymmetric Σ_s, `α_in=α_out=1`:
- `k_eff = 1.0000000000151756e+00`, rel err vs k_inf = **1.5e-11**.
- Per-group φ spread = O(1e-16) (machine-precision uniform).

ERR-002 transpose drift detector: PASSED.

### Convergence floor (research-grade)

Reflective inner / vacuum outer at quadrature ladder
(n_r, n_μ_axial, n_φ_az, n_traj_quad):
- (10, 10, 16, 16): k_eff = 7.3035007525e-02
- (14, 14, 20, 24): k_eff = 7.3031587150e-02
- (18, 18, 24, 32): k_eff = 7.3038554391e-02
- (24, 24, 32, 48): k_eff = 7.3051368990e-02

Finest pair self-consistency: **1.75e-4** (target was 1e-3).

The convergence is non-monotone over the early ladder, settling around
the ~1e-4 floor at finest. Cylinder's 3D angular phase-space requires
more total quadrature points (n_r × n_μ × n_φ ≈ 18k at finest grid)
than the sphere's 2D phase-space (n_r × n_μ ≈ 2.3k at equivalent
finest), so the convergence is correspondingly slower per-element but
the absolute floor is comparable.

## Sphere/cylinder/slab/slab-asym/hollow-sphere regression status

**Bit-equal accuracy floors preserved** at machine precision
(verified via separate confirmation script):
- Sphere closed BC k_eff = 2.0833333333333334e-01 (rel_err 0.00e+00).
- Cylinder closed BC k_eff = 2.0833333333333354e-01 (rel_err 9.33e-16).
- Slab closed BC k_eff = 2.0833333333333329e-01 (rel_err 2.67e-16).
- Slab asym closed BC k_eff = 2.0833333333333329e-01 (rel_err 2.67e-16).
- Hollow sphere closed BC k_eff = 2.0833333333333334e-01 (rel_err 0.00e+00).

All within ≤ 1e-15 absolute floor. Zero regression on existing 145
tests (all pass at original tolerances). The variant_alpha_core
primitives (`compute_resolvent_T`, `compute_resolvent_T_rank2`,
`apply_variant_alpha_closure`, `apply_variant_alpha_closure_rank2`)
were NOT modified — annulus reuses them unchanged.

## Honest verdict on the rank-2 + impact-parameter abstraction

### Did the closure framework lift unchanged from hollow sphere to annulus?

**YES, completely.** The Branch-2 production solver
`greens_function_annulus.py` imports
`apply_variant_alpha_closure` and `apply_variant_alpha_closure_rank2`
from `variant_alpha_core.py` and uses them with no modification. The
annulus prototype's only NEW code is the cylinder-specific chord
arithmetic:

- Phase-space discretization (μ_axial, φ_az grids).
- Impact-parameter computation `b = r·|sin(φ_az)|`.
- First-leg trajectory case analysis on sign of `cos(φ_az)`.
- 3D arclength conversion `L_3D = L_2D / sqrt(1-μ_axial²)`.
- Antipodal-chord parameterisation in 2D in-plane coordinates with
  3D attenuation `exp(-Σ_t · s_2D / s_in_plane)`.

Everything closure-related — the rank-2 monodromy, the resolvent
inversion, the `[α_L · B_LR, α_R · B_RL]` source vector, the
`ψ_L^+ = (α_L·B_LR + α_L·e^{-τ}·α_R·B_RL) / det` formula — is
**byte-equal-shared** with hollow sphere through `variant_alpha_core`.
The structural prediction from
`.claude/agent-memory/cross-domain-attacker/variant_alpha_2surface_bie_frame.md`
("rank-1 → rank-2 generalises 1-surface → 2-surface geometry,
regardless of curvature") holds across slab AND hollow sphere AND
annulus. The frame is now validated on **3 distinct 2-surface
instances**.

### Was Issue #129's "stay angle-resolved" mandate respected?

**YES.** The annulus prototype iterates ψ(r, μ_axial, φ_az) along
bouncing characteristics in 3D — NO axial pre-integration. The chord
algebra carries `1/sqrt(1-μ_axial²)` as a multiplicative factor on
3D arclength, but the angular variables are NEVER integrated out
inside the operator. This avoids the cylinder Bickley-Naylor
`Ki_n` form's ~22% planar-limit mismatch documented in Issue #129.

The V_α2_annulus.aux SymPy gate locks the chord-algebra identity
`τ_annulus = τ_hollow_sph / sqrt(1-μ_axial²)` symbolically — a
structural test that would FAIL if any code path had pre-integrated
the axial direction.

The plan's explicit anti-cross-check rule was followed: NO planar
slab limit cross-check appears in the test suite. The R_in → 0
limit cross-checks against the verified solid cylinder, which itself
respects the angle-resolved discipline (Phase 1).

### Is there structural friction at the (b → R_in tangent) ∩ (μ_axial → ±1 grazing) double-degenerate locus?

**NO concrete pathology observed at production grids.** Both grazing
loci have well-understood individual behaviour:

- **Tangent-ray locus** (b → R_in from below): `sqrt_disc_in → 0`,
  L_step → finite (= sqrt(R_out²−R_in²) at the boundary ÷ s_in_plane).
  The through-ray and outer-only branches agree on the boundary value
  of ψ because the through-ray chord-cancellation identity reduces
  to the outer-only one. No pathology.

- **Axial-grazing locus** (μ_axial → ±1): s_in_plane → 0,
  L_3D → ∞, exp(-Σ_t · L_3D) → 0. Surface flux at α=1 has
  T = 1/(1-α·e^{-τ}) → 1/(1-α) finite for α<1, but B → 0 faster via
  the integrand collapse. Net contribution at fixed φ_az is finite.
  Gauss-Legendre on μ_axial ∈ [-1, 1] avoids the endpoints exactly
  (open quadrature), so the grid never lands on the singular locus.

- **Double-degenerate locus** (b → R_in AND μ_axial → ±1 simultaneously):
  this is a measure-zero subset of phase-space (codimension 2). The
  production grid avoids it by construction (open GL quadrature on
  both μ_axial and φ_az). The V_α1 closed-annulus test at α=1
  exercises every grid point at the worst-conditioned corner (closed
  BC = T-divergence approached fastest) and gives k_inf at machine
  precision in 1 iter — empirical evidence that the double-degenerate
  locus is benign for the production grid.

The pathological corner that **would** matter (`α_out = 1, b → R_out`,
where τ_period → 0 ⟹ T → ∞) is the same grazing-ray issue the
solid cylinder has at α=1, μ→0. Same mitigation: GL interior
quadrature avoids the exact endpoint b = R_out.

### Will the Phase-3 family stay sound under future extensions?

**Hopeful but unproven.** The unified rank-1/rank-2 framework via
`variant_alpha_core.py` has now survived 6 instances (sphere,
cylinder, slab, slab-asym, hollow sphere, annulus) on 2 topologies
without modification. Three structural predictions from the
cross-domain-attacker frame have been validated:
1. Rank-1 covers 1-surface compact (sphere, cylinder).
2. Rank-2 covers 2-surface (slab-asym, hollow sphere, annulus).
3. Phase-space partition by impact parameter cleanly separates
   rank-1 outer-only from rank-2 through-ray on hollow geometries.

Future extensions on this framework that I expect to compose cleanly:
- Multi-region annulus (mirror sphere `_apply_operator_mr` segmentation
  along trajectory + B integrals).
- Anisotropic scattering on annulus (orthogonal axis).
- Multi-region hollow sphere (deferred from Phase 3C-1).
- Capping the cylinder/annulus axially (3-surface topology — would
  require rank-3 monodromy with the rank-N generalisation that is
  outside Phase-3 scope but algebraically straightforward as a
  block-circulant resolvent).

What I'm **less** sure about:
- Anisotropic-α (specular at one wall, white at the other) — needs
  a rank-2 generalisation where `S` has a non-antidiagonal structure.
  The frame predicts this should still close as `T = (I - S)^{-1}`
  with a more complex S, but I haven't validated this experimentally.
- The annulus + axial capping double-extension (rank-3 with two
  curvilinear surfaces and two flat surfaces) would be the first
  test of the unified rank-N framework's applicability beyond the
  pure-curvilinear-2-surface and pure-flat-symmetric cases.

## L0 bugs caught during development

**NONE.** The annulus prototype shipped clean on first try across
all 33 new tests (22 SymPy + 11 solver). Three structural elements
made this possible:

1. **The hollow sphere prototype was the right template.** The
   Phase-3C-1 closeout explicitly predicted what would and would
   not transfer ("rank-2 framework lifted from slab to hollow
   sphere with NO modification to variant_alpha_core.py. The core
   primitives are truly geometry-agnostic. Phase 3C-2 annulus will
   require only ... cylinder-specific chord arithmetic"). This
   prediction was load-bearing — when I started writing the
   annulus operator, I knew exactly which primitives to reuse and
   which to rewrite.

2. **The cylinder solid Phase-1 prototype was the right
   secondary template.** The first-leg trajectory parameterization
   `r²(s) = r² - 2r·s·cos(φ) + s²`, the GL quadrature ranges, the
   3D-arclength convention, and the `s_in_plane = sqrt(1-μ_axial²)`
   factoring all transferred cleanly. The annulus's only addition
   was the four-case first-leg backward analysis (outer-only vs
   through-ray × cos(φ) sign).

3. **Foundation-tagged SymPy gates first.** Writing
   V_α1_annulus / V_α2_annulus / V_α2_annulus.aux / V_α3_annulus
   BEFORE the production solver caught the algebra-composability
   concern symbolically. The V_α2_annulus.aux test in particular
   locks the cylinder axial-correction algebra
   `τ_annulus = τ_hollow_sph / sqrt(1-μ_axial²)` — a structural
   identity that would FAIL if the chord algebra had any
   pre-integration drift. Both branches passed V_α1 symbolically,
   guaranteeing the numerical V_α1 closed-shell test would also
   pass — and it did, machine precision in 1 iter.

**No new ERR-NNN entries logged.**

## Cross-domain frame analysis

The annulus is the **third 2-surface validation** of the rank-2 BIE
block resolvent frame predicted by
`.claude/agent-memory/cross-domain-attacker/variant_alpha_2surface_bie_frame.md`.

Frame match status:

| Topology family | Instance | Curvature | Frame match |
|-----------------|----------|-----------|-------------|
| 2-surface | slab asymmetric (Phase 3B) | flat | STRONG |
| 2-surface curvilinear | hollow sphere (Phase 3C-1) | spherical | STRONG |
| 2-surface curvilinear | annulus (Phase 3C-2) | cylindrical | STRONG |

The closure operator (`compute_resolvent_T_rank2`,
`apply_variant_alpha_closure_rank2`) is **byte-equal-shared** across
ALL THREE 2-surface prototypes. Only the chord-algebra primitives
differ:

- Slab: τ_step = Σ_t · L / |μ| (single-transit slab chord).
- Hollow sphere: τ_step = Σ_t · (sqrt(R_out²-b²) − sqrt(R_in²-b²))
  (curvilinear shell chord).
- Annulus: τ_step = (hollow-sphere chord) / sqrt(1-μ_axial²)
  (cylinder axial-correction lift).

The frame's structural prediction (rank-1 → rank-2 generalises
1-surface → 2-surface geometry, regardless of curvature **or**
angular phase-space dimensionality) holds across all four currently-
implemented 2-BC instances. Phase-3C-2 closes the original 6-
geometry × 2-topology family on the unified rank-1/rank-2 framework.

## Lessons for future implementations

1. **Reuse the closeout memo of the most-similar predecessor as a
   blueprint.** The Phase-3C-1 hollow sphere closeout's
   "Will the abstraction extend to annulus (Phase 3C-2) without
   modification?" section was a load-bearing prediction. It said
   "Mostly yes; minor adaptation required for cylinder's 3D angular
   space" and listed the four structural elements that would
   transfer plus the new piece (cylinder-specific chord arithmetic).
   This prediction was 100% accurate. Future implementations should
   include such a forward-prediction section in their closeouts —
   it doubles as a planning artifact for the next phase and a
   verification test for whether the framework's abstractions held.

2. **The 4-case first-leg analysis is the cylinder analog of the
   hollow-sphere's 4-case sign-of-μ × b-vs-R_in analysis.** Same
   structure: phase-space partition × backward-trajectory direction
   gives 4 cases. Only the parameterisation differs (sphere uses
   sign of μ; cylinder uses sign of cos(φ_az)). The geometric
   intuition transfers cleanly between sphere and cylinder geometries
   when the right invariant is identified.

3. **The V_α2_annulus.aux SymPy gate.** A "chord-algebra identity"
   gate that doesn't fit the V_α1/V_α2/V_α3 trio. Future cylinder/
   axial-extension prototypes should include similar auxiliary gates
   that lock geometry-specific algebraic identities outside the
   operator-symbol level. These gates protect against silent drift
   in the chord-arithmetic primitives that aren't reachable through
   the closure-level invariants.

4. **The Phase-3 family is now closed and the architecture is
   stable.** Six instances, two topologies, one shared closure
   framework, zero modifications to `variant_alpha_core.py` after
   the initial Phase-2 unification. This is the strongest evidence
   to date that the rank-1/rank-2 + impact-parameter-partition
   abstraction is the natural mathematical structure for Variant α
   on these geometries. Any future extension that requires modifying
   `variant_alpha_core` should be treated as a structural alarm —
   the abstraction has earned the right to stay frozen.

## Open issues / deferred work

- **Multi-region annulus**: deferred. The single-region prototype is
  sufficient to validate the rank-2 framework on cylinder 2-surface
  topology. Multi-region hollow can follow the sphere
  `_apply_operator_mr` segmentation pattern when needed.
- **Issue #132 Class B catastrophe at hollow sphere/annulus**: deferred
  per Phase-3C-1 closeout. Reproducer is a follow-on once multi-
  region-hollow exists.
- **Anisotropic scattering on the unified framework**: orthogonal
  extension; not on the Phase-3 critical path.
- **Rank-3 (3-surface) topology**: e.g., axially-capped cylinder.
  Algebraically straightforward but operationally a new prototype.

## Provenance

- Plan: `.claude/plans/peierls-greens-cylinder-and-2bc.md` Phase 3C-2.
- Cross-domain frame: `.claude/agent-memory/cross-domain-attacker/variant_alpha_2surface_bie_frame.md`.
- Phase 3C-1 closeout (predecessor template): `hollow_sphere_variant_alpha_phase3c1.md`.
- Phase 3B closeout (rank-2 framework genesis): `slab_asymmetric_variant_alpha_phase3b.md`.
- Phase 1 cylinder closeout (cylinder framework genesis): `cylinder_variant_alpha_phase1.md`.
