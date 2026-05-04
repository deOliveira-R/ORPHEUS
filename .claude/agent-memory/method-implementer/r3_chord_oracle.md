---
name: R3 ChordOracle Protocol + per-geometry implementations
description: 2026-05-04 R3 hindsight refactor on `refactor/r3-chord-oracle`; extracts the per-geometry chord-arithmetic from 5 trajectory_resolvent geometry modules into a unified ChordOracle Protocol with 6 concrete frozen-dataclass implementations; bit-equality preserved at IEEE-754 exact-bit level; -762 LOC net in production code.
type: project
---

# R3 — `ChordOracle` Protocol + per-geometry implementations (closeout)

**Branch**: `refactor/r3-chord-oracle` (off `main`)
**Date**: 2026-05-04
**Plan**: `.claude/plans/trajectory_resolvent_hindsight_refactor.md` §B1
**Hindsight memo**:
`.claude/agent-memory/method-implementer/trajectory_resolvent_hindsight_review.md`
**Frame**:
`.claude/agent-memory/cross-domain-attacker/trajectory_resolvent_foreign_frames.md`
(BaseAtlas in the fiber-bundle picture)

## What R3 ships

Three atomic commits:

1. **`3ffbf24`** — `feat(trajectory_resolvent): add R3 ChordOracle
   Protocol + 6 concrete oracles`. Pure addition. Defines
   `ChordOracle` (runtime-checkable Protocol) + six concrete frozen
   dataclasses: `SphereChordOracle`, `MultiRegionSphereChordOracle`,
   `CylinderChordOracle`, `SlabAsymmetricChordOracle`,
   `HollowSphereChordOracle`, `AnnulusChordOracle`. 21 foundation
   tests pin np.int64-view bit-equality between every oracle's
   `apply_operator` and the legacy inlined `_apply_operator_*`.
2. **`5fa83ed`** — `refactor(trajectory_resolvent): rewire
   _apply_operator_* facades to ChordOracle`. Replaces every legacy
   inlined body with a thin facade that constructs the matching
   oracle and delegates to its `apply_operator`. Bit-equality
   preserved end-to-end: 12 canonical k_eff configurations verified
   bit-equal at IEEE-754 via `float.hex(k_eff)`. Net: **-762 LOC** in
   production code (-868 deletions, +106 additions).
3. **`a8a631f`** — `docs(trajectory_resolvent): document ChordOracle
   as the base atlas + export from package`. Adds the *"The chord
   oracle — base atlas of the fiber bundle"* subsection to
   `docs/theory/trajectory_resolvent.rst` and exports the oracles +
   Protocol from the package's `__init__.py`.

## Bit-equality verification

The R3 refactor is a pure relocation of code: every FP operation
runs in the same order. IEEE-754 reproducibility then guarantees
identical numerical output. Three layers of evidence:

1. **Foundation tests at the operator level** (21 tests at
   `tests/derivations/test_trajectory_resolvent_chord_oracle.py`):
   for each oracle and three alpha configurations, the
   `oracle.apply_operator(source, σ_t, n_traj_quad)` output is
   compared bit-for-bit (via `view(np.int64)`) against the legacy
   `_apply_operator_*` body. Every (oracle, alpha) pair passes.
2. **k_eff snapshot at the solver level** (12 canonical
   configurations via `/tmp/baseline_keff.py` snapshot): sphere /
   cylinder / slab / slab_asym / hollow_sphere / annulus × 1G/MG ×
   closed/vacuum/asymmetric. Every k_eff is byte-equal to the pre-
   refactor baseline at `float.hex` precision (IEEE-754 exact-bit).
3. **End-to-end test suites**:
   - `tests/derivations/test_peierls_greens_function_*` (8 files,
     93 tests) + `test_trajectory_resolvent_billiard.py` (15 tests) +
     `test_trajectory_resolvent_power_iterate.py` (6 tests) +
     new `test_trajectory_resolvent_chord_oracle.py` (21 tests) →
     **135 pass** (114 baseline + 21 new).
   - `tests/cross_method/` (43 tests) + 4 xverif/Garcia/MG/MR test
     files → **126 pass** (zero regression).

## The chord-arithmetic concept made first-class

Before R3, chord arithmetic was scattered:

- `greens_function.py` — sphere first-leg + bounce-period (lines
  184-188); MR segmentation (lines 717-794).
- `greens_function_cylinder.py` — `_first_leg_2d_chord` (lines 165-
  192), `_impact_parameter` (lines 194-200), `_bounce_period_2d_chord`
  (lines 202-212), `_apply_operator_cylinder` 3D-lift body (lines
  220-350).
- `greens_function_slab_asymmetric.py` — sign-of-μ branching, dual
  single-transit B integrals (lines 240-302).
- `greens_function_hollow_sphere.py` — b-partition routing, L_step
  shell-traversal, sign-of-μ for inner/outer first arrival (lines
  220-410).
- `greens_function_annulus.py` — cylinder + b-partition + 3D lift
  (lines 286-444).

After R3, the chord arithmetic lives in ONE module (`chord_oracle.py`)
with six concrete implementations + helpers (`_trajectory_segments_oracle`,
`_chord_segments_oracle`, `_region_at_radius_oracle` for MR sphere).
Each implementation is roughly 50-200 lines of pure arithmetic with
the closure call at the end. The legacy `_apply_operator_*` functions
collapse to ~10-line facades.

## Why R3 is "just rewire" and bit-equality is preserved by construction

The oracles were extracted *verbatim* from the legacy `_apply_operator_*`
functions: the body of each oracle's `apply_operator` is byte-equivalent
to the body of the function it replaces. The legacy function then
becomes a facade that constructs the oracle and calls `.apply_operator()`
— a one-step indirection that does not change any FP arithmetic.

The legacy private helpers (`_first_leg_2d_chord`,
`_bounce_period_2d_chord`, `_impact_parameter`, `_trajectory_segments`,
`_chord_segments`, `_region_at_radius`) stay in their original modules
for back-compat:

- `_region_at_radius` is still consumed by
  `solve_greens_function_sphere_mr*` for the `region_at_node`
  labelling — load-bearing.
- The cylinder helpers are no longer called inside the geometry
  module but kept for any external code (notebooks, ad-hoc scripts)
  that may import them. Their convention (single underscore) is
  "private but accessible".

## Architectural placement (Cardinal Rule 2)

Per `feedback_unify_after_two_instances.md`, the R3 oracles stay in
`trajectory_resolvent/` rather than promoting to `derivations/common/`
because today only the Variant α family consumes the abstraction. The
plan's §R4 schedules the promotion when:

- The cross-method regression net (test-architect deliverable, in
  flight on a parallel branch) lands and bedds in.
- A second consumer materializes — most likely
  `peierls_nystrom` (production CP variant), which already has
  similar primitives expressed differently
  (`chord_half_lengths` in `derivations/common/kernels.py` is the
  proto-ChordOracle for concentric annuli/shells).
- Optionally, a third consumer (MoC reference solver).

The pattern from prior phases held — *"build cylinder first, test
it, then unify and see if the tests still hold"* — applied
recursively to abstractions: the oracles are validated inside
trajectory_resolvent for one merge cycle before they cross the
common-layer boundary.

## Cross-domain frame validation (BaseAtlas)

The cross-domain-attacker memo
(`trajectory_resolvent_foreign_frames.md` §"Frame 4: Fiber bundle /
G-structure") named the chord oracle as the **BaseAtlas** in the
fiber-bundle frame:

- *Total space*: the Variant α operator family.
- *BaseAtlas* (R3): per-chart spatial parametrisation, mapping the
  abstract orbit space M/G into ray-traversal data on M.
- *AngularFiber* (R2 partially, R5 deferred): residual symmetry
  group's parametrisation — μ for sphere; (μ_axial, φ_az) for
  cylinder; signed μ for slab.
- *ClosurePrimitives* (Phase 2): rank-1/rank-2 resolvent +
  Variant α closure invariants in `variant_alpha_core.py`.

R3 institutionalises the BaseAtlas layer as a Protocol; the second
"frame validation" instance (after R2's Billiard math-heart) is now
load-bearing for any future structural unification (Frame 1: Tensor
networks/MPO; Frame 6: PIMC verification).

## What surfaced during the refactor

- **Cylinder helper preservation**. `_first_leg_2d_chord`,
  `_impact_parameter`, `_bounce_period_2d_chord` are no longer
  called by the post-rewire `_apply_operator_cylinder` facade.
  Kept as-is because their convention (single underscore) makes
  them potentially consumed externally; deleting them would be
  a breaking API change for any notebook that imports them. If
  the abstraction promotes to `derivations/common/` (R4),
  consider deprecating these in favour of the
  `CylinderChordOracle.compute_chord_data(...)` family of
  inspector methods (currently not exposed; the oracle's API is
  the single `apply_operator` call).
- **MR oracle ignores `sigma_t` argument**. The
  `MultiRegionSphereChordOracle` carries `sigma_t_per_region` as a
  baked-in attribute and ignores the `sigma_t` keyword of the
  Protocol's `apply_operator`. This is API-uniform but slightly
  surprising — documented in the oracle's docstring. A cleaner
  design would split the Protocol into single-region and
  multi-region forms, but that creates more types than it eliminates.
- **The slab-symmetric solver auto-piggybacks on `SlabAsymmetricChordOracle`**.
  `solve_greens_function_slab` delegates to `solve_greens_function_slab_asymmetric`
  with `α_left = α_right = α` (the ERR-035 fix from a prior closeout).
  Post-R3 this is unchanged: the slab-asymmetric oracle is the
  only slab oracle, and the symmetric solver is a thin BC adapter.

## Foundation tests catalog

The 21 R3 foundation tests pin:

- 3 sphere oracle bit-equality (alpha ∈ {1.0, 0.5, 0.0})
- 3 MR-sphere oracle bit-equality
- 3 cylinder oracle bit-equality
- 3 slab-asym oracle bit-equality (alpha pairs (1,1), (0.5, 0.8), (0,0))
- 3 hollow-sphere oracle bit-equality
- 3 annulus oracle bit-equality
- 1 frozen-dataclass invariant
- 1 Protocol conformance for all 6 oracle types
- 1 sphere oracle Protocol conformance

Total wall time: **0.46 s**. The tests run inside the foundation tier
(no equation labels needed; `@pytest.mark.foundation`).

## What's NOT in R3 (deferred per plan)

- **R4**: promotion of `ChordOracle` to
  `orpheus/derivations/common/chord_oracle.py`. Gated on (a) one
  merge cycle of stability + (b) at least one external consumer
  (peierls_nystrom or MoC).
- **R5 (B4)**: parametrise `extract_scalar_flux` as a single
  function over the angular axes. Low priority; opportunistic.
- **R5 (B5)**: `mesh.volume_integral` abstraction for fission rate.
- **R5 (A4)**: MG-as-tensor-product factoring.
- **R5 (A5)**: unified `solve(geometry, ...)` facade.
- **B6** (`MultiRegionChordOracle` decorator over a base oracle —
  shared structure between cylinder MR + slab MR + sphere MR). The
  current `MultiRegionSphereChordOracle` is a sibling of
  `SphereChordOracle`, not a decorator. When cylinder MR / slab MR
  ship (Issues #129, #145), the decorator pattern is the natural
  unification.

## Branch hygiene

Three commits, all on `refactor/r3-chord-oracle`. No parallel-agent
trampling: `git branch --show-current` checked before each commit;
`git status` checked after each commit. Two parallel agents (test-
architect on cross-method TransportSolver Protocol + a third on
trajectory_resolvent test net) are touching adjacent territory but
not the same files.

## Status

- **Production code**: 5 facades rewired; 1 new module
  (`chord_oracle.py`, 1023 LOC).
- **Tests**: 21 new foundation tests. 0 regression in 240+ existing
  tests. Bit-equality at IEEE-754 exact-bit level for 12 canonical
  k_eff configurations.
- **Sphinx**: `-W` build clean.
- **Plan**: §R3 complete. §R4 (`ChordOracle` promotion to
  `derivations/common/`) deferred per `feedback_unify_after_two_instances.md`.
- **Memory**: this entry.

Three commits ready to merge to main when the user requests
(commit hashes `3ffbf24`, `5fa83ed`, `a8a631f` on
`refactor/r3-chord-oracle`).
