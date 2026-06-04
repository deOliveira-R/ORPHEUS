# Plan: Variant α Green's function — cylinder + 2-BC topologies + unified architecture

**Author**: Claude Opus 4.7, 2026-05-02
**Predecessor**: Plan `peierls-greens-function-approach.md` (sphere Variant α landed end-to-end: 1G, MG, MR, MR-fixed-source).
**Proposed branch**: `feature/peierls-greens-topology-unification` off `main`.

## Frame finding (cross-domain attacker, 2026-05-02)

The Variant-α "geometric bounce sum" T(μ_surf) = 1/(1 − α·exp(−Σ_t·L_period))
is the **rank-1 multiplication-operator case** of a resolvent

```
T = (I − S)^{-1},   S = α · R_chord
```

where `S` is the **boundary-to-boundary scattering operator** on
surface phase-space (BIE frame). The 2-surface topologies (slab,
annulus, hollow sphere) are the **rank-2 block case of the same
resolvent**, with:

- vacuum on a face → zeroing the corresponding row of S
- mixed α_in / α_out → diagonal entries of S
- two-chord-period dynamics → off-diagonal entries

**Cylinder remains rank-1** (1-surface compact); only the chord
algebra ρ_max(r, θ), r'(r, ρ, θ) differs from sphere.

**Pathology predictor**: T is bounded iff `1 ∉ ess_range(S)`. For
α=1 grazing rays μ→0 the period collapses (L_period → 0,
exp(−Σ_t·L_period) → 1, S → 1). Sphere has this; slab is immune;
cylinder/annulus inherit iff a vanishing-chord locus reaches the
boundary.

Memo: `.claude/agent-memory/cross-domain-attacker/variant_alpha_2surface_bie_frame.md`.

## Phasing

**Ordering note (user override, 2026-05-02):** original draft put the
operator-theoretic refactor first. User overrode → "build cylinder
first, test it, then unify and see if the tests still hold." Reasoning:
Nyström unification under `CurvilinearGeometry` was complex and
underperforms; unification can introduce accuracy trade-offs that must
be measured, not assumed. Phases reordered: cylinder Phase 1
(standalone, mirror sphere structure), unification Phase 2 (with
both instances in hand and accuracy benchmarks captured). See
[`feedback_unify_after_two_instances.md`](file:///home/vscode/.claude/projects/-workspaces-ORPHEUS/memory/feedback_unify_after_two_instances.md).

### Phase 1 — Cylinder Variant α (1-surface compact, standalone)

Goal: produce an independent, research-grade cylinder Variant-α
implementation that mirrors the sphere structure. Capture the
achieved accuracy as the BENCHMARK any later unification must
preserve.

(was Phase 2 in original draft — promoted to Phase 1)

### Phase OLD-1 — Operator-theoretic refactor (deferred to Phase 2)

Goal: name the resolvent. Establish the abstraction with **two**
instances in hand (sphere + cylinder), proven equivalent to both
standalone implementations (both test suites stay green at their
captured accuracy).

Tasks (deferred):

Tasks:
- A1. In `greens_function.py`, refactor `T(μ_surf) = 1/(1 − …)` into
  a named primitive `_compute_resolvent(S, …)` with S as the
  boundary-to-boundary scattering operator. Rank-1 sphere case is a
  scalar; the function signature must accept rank-2 future cases.
- A2. Extract the three Variant-α primitives — F (first-leg), B
  (bounce-period), S (scattering operator) — into module-level
  pure functions on top of `CurvilinearGeometry`. Sphere becomes the
  reference implementation.
- A3. Rewrite Branch-1 SymPy proofs as resolvent identities:
  - V_α1: closed-sphere bounce sum → direct geometric expansion
    of `(I − S)^{-1}` for rank-1 S.
  - V_α2: T_00 = P_ss → Schur-complement / matrix-element identity
    of the resolvent (zeroth angular moment).
  - V_α3: α=0 reduction → S = 0 trivial case.
  Files: `origins/specular/greens_function.py` rename functions to
  reflect resolvent semantics; tests in
  `tests/derivations/test_peierls_greens_function_symbolic.py`
  retarget the new identities.
- A4. All existing sphere tests must pass unchanged (1G, MG, MR,
  MR-fixed-source, PS-1982, Garcia 2021, vacuum).
- A5. Sphinx page `peierls_greens.rst`: add §"Operator-theoretic
  reformulation" introducing S, T = (I−S)^{-1}, and the rank-1
  vs rank-2 distinction. Don't move existing math; add the
  unifying section.

Acceptance (now for Phase 2): green tests for BOTH sphere AND
cylinder at captured accuracy floor + clean Sphinx build + V_α1/V_α2/V_α3
re-derivations covering both geometries. If unified architecture
regresses cylinder OR sphere accuracy beyond an agreed tolerance,
HALT — re-evaluate per `feedback_unify_after_two_instances.md`.

### Phase 1 detailed tasks (cylinder standalone)

Goal: prove a second geometry can carry Variant α independently.
Mirror sphere structure as a self-contained module. Only chord
algebra is new mathematically.

Tasks:
- B1. Branch-1 SymPy proofs in `origins/specular/greens_function_cylinder.py`:
  - V_α1_cyl: closed-cylinder bounce sum (chord = 2·ρ_max(r,θ),
    angular variable in (r,μ_r,φ) — KEEP angle-resolved per #129).
  - V_α2_cyl: T_00_cyl = P_ss_cylinder (use existing
    `compute_P_ss_cylinder` from geometry.py).
  - V_α3_cyl: α=0 reduction.
- B2. Branch-2 implementation `solve_greens_function_cylinder` in
  `greens_function.py` (or peer module). Mount on
  `CurvilinearGeometry(kind='cylinder-1d')`. Reuse the F/B/S
  primitives from Phase 1. Do NOT pre-integrate axial direction —
  iterate ψ(r, μ_r, φ_az) along the bouncing characteristic.
- B3. Cross-check: at α=1 perfect-reflection homogeneous, k_eff
  must equal k_inf = νΣ_f/Σ_a EXACTLY (analytic invariant —
  reproduce sphere's Variant-α exactness property).
  At α=0 vacuum, cross-check against existing cylinder vacuum
  reference (`compute_P_esc_cylinder_3d`).
- B4. Tests: `tests/derivations/test_peierls_greens_function_cylinder.py`
  with at least: 1G α=1 k_inf, 1G α=0 vs vacuum reference,
  multi-group, and a heterogeneous 2G case.
- B5. **Issue #129 acknowledgment**: do NOT cross-check against thin
  slab via planar limit (chord ≠ L/|μ| in pre-integrated cylinder
  Ki_n). Cross-check against MC if possible, else against the
  Branch-1 SymPy reference at low quadrature order.

Acceptance: cylinder Variant-α passes its own L0/L1 tests AND the
shared S/T abstraction was re-used (no copy of resolvent code).

### Phase 3 — Slab Variant α (2-surface, rank-2 S, simplest)

Goal: prove the rank-2 block resolvent path with the simplest
geometry. This is where the abstraction earns its keep.

Tasks:
- C1. Branch-1 derivations: rank-2 S for slab with reflective
  α_left, α_right. Two-chord period — the bouncing characteristic
  alternates between left and right faces. S becomes:
  ```
      [ 0                       α_left · exp(−Σ_t·L) ]
  S = [                                              ]
      [ α_right · exp(−Σ_t·L)   0                    ]
  ```
  (off-diagonal because each face's outgoing flux feeds the OTHER
  face's incoming flux). Resolvent T = (I − S)^{-1} closed-form.
- C2. Branch-2 `solve_greens_function_slab`. Mount on existing
  `CurvilinearGeometry(kind='slab')`. Reuse F/B primitives;
  resolvent path uses the rank-2 branch.
- C3. Cross-check against existing slab analytical solutions
  (Pomraning, slab E_n forms, white-BC reference).
- C4. Verify the unification: existing sphere code (rank-1, scalar
  T) and new slab code (rank-2, 2x2 T) both call into a single
  resolvent primitive `_compute_resolvent`. If they don't, the
  abstraction is wrong — fix Phase 1 before continuing.

Acceptance: slab Variant-α passes L0/L1 tests + the resolvent
primitive is shared with sphere (no duplication).

### Phase 4 — Hollow sphere & annulus (2-surface curvilinear, rank-2 S)

Goal: extend to the hardest case — 2-surface curvilinear with
impact-parameter-dependent trajectory partitioning.

Tasks:
- D1. Trajectory classification by impact parameter b:
  - b > r_inner: ray bounces only on outer surface (rank-1 S, like
    the solid case).
  - b ≤ r_inner: ray bounces inner-outer-inner-… (rank-2 S, like
    slab but with curved chord algebra).
  Branch-1 SymPy: derive both cases and the b-dependent piecewise
  resolvent.
- D2. Branch-2: `solve_greens_function_hollow_sphere` and
  `solve_greens_function_annulus`. The Nyström side already has
  hollow primitives (`compute_hollow_{cyl,sph}_transmission`,
  `compute_G_bc_inner`); the Variant-α side just needs the F/B
  primitives extended for piecewise trajectory.
- D3. Cross-check against the existing hollow Nyström rank-N closure
  (where it converges — Issue #110 Phase F).
- D4. **Class B catastrophe sentinel**: Issue #132 said rank-N MR
  with strong scatterer overshoots by 57%. Variant α structurally
  avoids this. Reproduce the Issue #132 case for solid sphere
  (already done) AND for hollow sphere AND for annulus, demonstrate
  Variant α gets the right answer where rank-N closure does not.

Acceptance: hollow & annulus pass L0/L1 tests, the resolvent
primitive carries them with no duplication, AND the Class B
catastrophe is structurally cleared in 2-surface topology.

## Verification matrix at completion

| Geometry | Topology | Reference | L0 | L1 | Status |
|----------|----------|-----------|----|----|--------|
| Sphere | 1-surface compact | PS-1982 + Garcia 2021 | ✓ | ✓ | shipped |
| Cylinder | 1-surface compact | k_inf at α=1 + MC at α<1 | new | new | Phase 2 |
| Slab | 2-surface | Pomraning E_n + white-BC | new | new | Phase 3 |
| Hollow sphere | 2-surface | hollow Nyström rank-N (where converged) | new | new | Phase 4 |
| Annulus | 2-surface | hollow cyl Nyström rank-N (where converged) | new | new | Phase 4 |

## Risks acknowledged

- **#129 planar limit**: cylinder Ki_n has axial pre-integrated;
  cyl ↔ slab cross-check disagrees by 22%. Variant α avoids this by
  staying angle-resolved. Don't cross-check via planar limit.
- **#132 Class B MR catastrophe**: rank-N closure fails on solid
  cyl/sph MR. Variant α structurally avoids; reproduce in 2-surface.
- **Phase 5 Hadamard singularity**: angle-integrated specular kernel
  is hypersingular. Variant α avoids by iterating angle-resolved ψ
  along bouncing characteristics. Maintain that discipline in cyl,
  slab, hollow.
- **#100, #103 rank insufficiency**: Variant α has no closure
  operator → not subject to rank gating. Don't regress into rank-N
  during 2-surface extension.
- **α=1 grazing-ray pathology**: T unbounded when L_period → 0
  (μ→0 in sphere, vanishing-chord locus in annulus). Mitigate by
  quadrature endpoint avoidance + ε-regularization at the surface
  (already in sphere code).

## Issues this plan touches / closes

- #110 Phase F (hollow-core): partially advanced (hollow Variant α
  closes the verification ground for Nyström hollow).
- #111 Phase G (slab Peierls in unified framework): the operator
  refactor IS the unification axis Phase 1 establishes.
- #129 planar-limit physics: NOT closed — but the documented
  "don't planar-cross-check" discipline becomes the project policy.
- #132 Class B catastrophe: structurally cleared in 2-surface
  topologies once Phase 4 ships.
- #100, #103 rank-N closure plateau: Variant α offers a
  rank-independent path; this plan does not retire the closures
  but provides a verified alternative.
- #101 cylinder Peierls reference Ki_1 Nyström: superseded by
  cylinder Variant α (Phase 2).
