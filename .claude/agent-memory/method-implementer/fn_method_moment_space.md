---
name: F_N method MomentSpace class — 2nd math-heart instance
description: Closeout memo for the MomentSpace class on feat/fn-method-moment-space — 2nd concrete instance of the math-heart pattern (1st is Billiard in trajectory_resolvent), shared CriticalSolution/FluxSolution result types, math-rich docstring + theory section.
type: project
---

# F_N method MomentSpace — 2nd math-heart instance

## Branch & commits

- **Branch**: `feat/fn-method-moment-space` (off main, 3 commits ahead of origin/main)
- **Commits**:
  1. `4cb92d1` — `feat(derivations/common): add CriticalSolution + FluxSolution shared types`
  2. `4c8c68b` — `feat(fn_method): add MomentSpace class — 2nd math-heart instance`
  3. `9eadc64` — `docs(fn_method): add Mathematical structure of the F_N moment space section`
- **Date**: 2026-05-04
- **Branch hygiene**: NOT pushed — user merges manually per brief

## What shipped

### 1. Shared result types (commit 1)

`orpheus/derivations/common/solution_types.py` — the FIRST cross-method
shared vocabulary for continuous reference solvers. Two frozen
dataclasses:

- `CriticalSolution(eigenvalue, eigenvalue_kind, parameter_value,
  parameter_kind, converged, metadata)` — result of any
  critical-configuration root-find.
- `FluxSolution(spatial_nodes, scalar_flux, angular_flux,
  angular_nodes, eigenvalue, eigenvalue_kind, spatial_units,
  metadata)` — result of any flux-reconstruction call.

Both populate the cross-method gates (`MomentSpace.solve_critical()`,
`Billiard.solve_critical()`) so a downstream consumer can hold a
`CriticalSolution` without knowing which pillar produced it.

The `metadata: dict[str, Any]` field is open-ended on purpose — every
solver carries its own diagnostic vocabulary, and forcing a closed
schema would be premature standardisation. Cross-method protocol
adapters document expected keys per-adapter.

### 2. MomentSpace class (commit 2)

`orpheus/derivations/continuous/fn_method/moment_space.py` —
the F_N method's natural mathematical home. Class composition:

- **Geometry**: slab / sphere / infinite (cylinder rejected — out of
  pillar per Westfall-Metcalf 1972).
- **Materials**: production-protocol `dict[int, Mixture]`.
- **F_N order**: default 9 (typical Grandjean-Siewert operating point).
- **Flux-reconstruction strategy**: `atkinson_nystrom` (default,
  ERR-036 fix), `legacy_gl` (diagnostic), `none` (solve-only).

Methods:

- `from_problem(materials, geometry, ...)` — recommended factory.
- `solve_critical(...)` → `CriticalSolution` — dispatches to
  slab / sphere / infinite-medium branches based on `GeometrySpec`.
- `reconstruct_flux(n_panels=256, z_eval=None)` → `FluxSolution`.
- `solve_kinf()` → `float` — convenience for k_inf-only access.
- `c` property — returns `(σ_s + νσ_f)/σ_t` for 1G mixtures.

The class is a THIN FACADE over the existing function-level API
(`solve_fn_slab_bare_critical` / `solve_fn_sphere_bare_critical` /
`compute_kinf_*`). Bit-equality with the function-level API
preserved (verified by `float.hex()`).

### 3. Theory section (commit 3)

New section "Mathematical structure of the F_N moment space" in
`docs/theory/fn_method.rst` — ~537 lines of graduate-textbook
narrative covering:

- L²(domain × sphere) Hilbert-space framing
- Half-range decomposition at boundaries
- Galerkin projection: minimise residual orthogonal to basis
- Galerkin orthogonality + system equation
- Collocation closure (slab Grandjean-Siewert vs sphere
  Siewert-Thomas Chebyshev-interior grid; geometry-sign s = ±1
  duality)
- The critical condition det M = 0
- Truncation error analysis (smoothness-dependent O(N⁻ᵖ))
- Multi-region extension via block transfer matrices
- Connection to flux reconstruction + Atkinson product-Nyström
  (ERR-036)
- Relationship to Billiard via the Sanchez-Chandrasekhar
  three-meanings taxonomy
- The architectural payoff (production-protocol inputs, shared
  result types, math-rich docs, bit-equality)

This is the user's explicit "every time I read it I learn"
requirement materialised.

## Test gate

`tests/derivations/test_fn_method_moment_space.py` — 14
foundation-tagged tests pinning:

- Production-protocol input acceptance (`from_problem`)
- Cylinder rejection (out of pillar)
- Unknown flux strategy rejection
- `c` property formula
- **Bit-equality with function-level API** for slab, sphere,
  1G k_inf, 2G k_inf with upscatter (URRb-2-0-IN params)
- `CriticalSolution` metadata fields (n_modes, determinant_residual,
  raw_result)
- `FluxSolution` returned by `reconstruct_flux` (slab atkinson +
  sphere KLL paths)
- `flux_reconstruction='none'` raises in `reconstruct_flux`
- Sub-critical input rejection

All 14 tests pass in 0.79s. Total test suite: 192 baseline + 14 new
= 206 passing.

## What was NOT implemented (deferred scope)

- **Multi-group F_N spatial extension**: `MomentSpace.solve_critical`
  for slab/sphere is currently 1G-only. The Siewert-Thomas 1986 2G
  machinery is shipped at the function level but not wired through
  the class facade — `solve_critical` raises `NotImplementedError`
  for 2G+ slab/sphere. Multi-group infinite medium IS supported.
- **Reflected-slab F_N path**: `solve_fn_slab_reflected_critical`
  (Neshat-Maiorino 1980) is not yet a `MomentSpace` factory variant.
  The class focuses on bare-critical configurations; reflected
  problems use the existing function API.
- **Migrating cross_method/adapters.py to use MomentSpace**: per the
  brief's "follow-up commit if cleaner" guidance, the `FNSlabAdapter` /
  `FNSphereAdapter` adapters in `tests/cross_method/adapters.py` still
  call the function-level API directly. Migrating them would buy
  documentation locality (the cross-method Sphinx sees `MomentSpace`
  rather than 6 separate function names) but no behavioural change.
  Left as a future cleanup.

## Cross-coordination with parallel agent (R2-Billiard)

The parallel agent was building `Billiard` in
`orpheus/derivations/continuous/trajectory_resolvent/` simultaneously
on `refactor/r2-billiard-class`. The shared `CriticalSolution` /
`FluxSolution` types in `derivations/common/solution_types.py` are
the coordination point — I arrived first, defined them in
`derivations/common/`, and the parallel agent will import them.

**Branch hygiene incident** (× 3): the parallel agent's checkout
trampled mine three times during the session. Each time their
`git checkout` switched HEAD on the shared working directory. The
recovery sequence:

1. Verify my branch's tip via `git log feat/fn-method-moment-space`
   (commits are SAFE in git regardless of current checkout — only
   the working-tree contents are trampled).
2. Stash any cross-pollination from the parallel agent
   (`docs/theory/trajectory_resolvent.rst`,
   `orpheus/.../trajectory_resolvent/__init__.py`,
   `orpheus/.../trajectory_resolvent/billiard.py`,
   `tests/derivations/test_trajectory_resolvent_billiard.py`).
3. `git checkout feat/fn-method-moment-space`.
4. Confirm my files are back (`ls
   orpheus/derivations/continuous/fn_method/moment_space.py`).

The previous closeout memo
(`.claude/agent-memory/method-implementer/rename_meshtemplate_to_geometryspec.md`)
documented this exact incident pattern. **Lesson reinforced** —
parallel agents touching adjacent territory must verify
`git branch --show-current` before each commit. The commits never
got contaminated; the working-tree state did.

## Sphinx -W status

Build succeeded. Two warnings encountered during development (now
fixed):

- `:doc:`/skills/vv-principles`` — `/skills/` is not a Sphinx
  document tree. Convention: use literal `.claude/skills/...`
  references (no Sphinx directive). Fixed in commit 3.
- `:doc:`/skills/algebra-of-record`` — same. Fixed in commit 3.

## Architectural payoff

This is the **2nd concrete instance of the math-heart pattern** across
the project. The 1st is `Billiard` in trajectory_resolvent (parallel
work). Per the project's `feedback_unify_after_two_instances.md` memory,
the unifying Protocol over math-heart classes can ONLY be designed
after ≥2 working instances exist. `MomentSpace` is what makes
`Billiard` no longer a one-off — the Protocol design is the next
dispatch (test-architect or another method-implementer for design,
once both instances are bedded in).

The shared `CriticalSolution` / `FluxSolution` are the **eager
unification** (small, stable dataclasses; cost of premature
unification is trivial). The behavioural Protocol over the math-heart
classes themselves is the **deferred unification** — gated on
empirical observation of variation patterns once 2 instances exist.

## Manifest (algebra-of-record discipline check)

- [x] Branch-1 SymPy module: pre-existing across
      `orpheus/derivations/continuous/fn_method/origins/*` —
      ungutted, ungrown by this work.
- [x] Foundation-tagged test gate at
      `tests/derivations/test_fn_method_moment_space.py`
      (14 tests, all passing).
- [x] Branch-2 production code at
      `orpheus/derivations/continuous/fn_method/moment_space.py`
      (the class facade) + the existing function-level solvers
      it wraps.
- [x] L1 cross-check tests pre-existing at
      `tests/derivations/test_fn_la13511_*.py` — unchanged by this
      work; the class facade preserves them via bit-equality.
- [x] Sphinx section with `:label:` cross-refs +
      `:mod:` cross-references at `docs/theory/fn_method.rst`
      (the new "Mathematical structure" section). NOT a stub —
      full rich narrative shipped (the user's "every time I read
      it I learn" requirement was the load-bearing deliverable).
- [x] Closeout memo (this file).
- [ ] DISPATCH_REQUEST to archivist — **NOT EMITTED** per user's
      `feedback_subagent_autonomy.md` guidance to dispatch only when
      asked. The theory section already IS the rich narrative; no
      stub-to-narrative expansion is owed.
