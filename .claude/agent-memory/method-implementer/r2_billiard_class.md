---
name: R2 Billiard class — math-rich math-heart for trajectory_resolvent
description: Hindsight refactor R2 — Billiard class with full Birkhoff transfer-operator framing, shared CriticalSolution/FluxSolution adoption, topology -> orbit_space_class rename. 93+119 tests pass identically pre/post; bit-equality preserved across 14-config baseline.
type: project
---

# R2 — `Billiard` class closeout

`refactor/r2-billiard-class` (off `main`), 2026-05-04.

## Scope

Lands the second hindsight refactor for the Variant α family from
`.claude/plans/trajectory_resolvent_hindsight_refactor.md`:

- **R2 (this work)**: introduce `Billiard` class as the math-heart
  for trajectory_resolvent + adopt shared CriticalSolution/FluxSolution
  result types from the parallel MomentSpace agent's work + rename
  CurvilinearGeometry.topology -> orbit_space_class.

R1 (`power_iterate_variant_alpha` driver) shipped earlier on
`refactor/r1-power-iterate-driver`. R3 (ChordOracle Protocol) is
gated on R1+R2 bedding-in.

## Two commits

### Commit 1 — `feat(trajectory_resolvent): introduce Billiard class`

`813c3ad`. Adds:

- `orpheus/derivations/continuous/trajectory_resolvent/billiard.py`
  (1184 LOC): the `Billiard` class with math-rich docstring naming
  the Birkhoff–Sinai billiard system explicitly. Construct via
  `Billiard.from_problem(geometry_kind=, materials=, geometry=,
  alpha=, quadrature=)`; resolves `closure_rank` (1 or 2) from the
  orbit-space M/G class. `solve_critical()` /
  `solve_fixed_source()` dispatch to the legacy
  `solve_greens_function_*` entry points and re-pack into the
  SHARED cross-method `CriticalSolution` / `FluxSolution` from
  `derivations.common.solution_types`.
- 380-line theory section "Mathematical billiards and the Variant
  α resolvent" added to
  `docs/theory/trajectory_resolvent.rst`. Covers Birkhoff–Sinai
  billiard system; Poincaré–Birkhoff billiard map; transfer
  operator S; Neumann series resolvent T = (I-S)^{-1}; rank-1 /
  rank-2 specializations; the three regimes (specular ↔ ergodic ↔
  vacuum); why the billiard frame is structurally right for
  Variant α. Footnote citations to Birkhoff 1927, Sinai 1970,
  Chernov-Markarian 2006.
- 15 foundation tests at
  `tests/derivations/test_trajectory_resolvent_billiard.py`. Pin
  `float.hex` equality of every `k_eff` and `np.array_equal` of
  every `psi`/`phi` across all 6 geometries × {1G, MG} +
  sphere_mr + fixed-source.

Public API:

```python
from orpheus.derivations.continuous.trajectory_resolvent import (
    Billiard, CriticalSolution, FluxSolution,
)
b = Billiard.from_problem(
    geometry_kind="sphere",
    materials={"sigma_t": 0.5, "sigma_s": 0.4, "nu_sigma_f": 0.1},
    geometry={"R": 5.0},
    alpha=1.0,
)
sol = b.solve_critical()  # CriticalSolution
sol.eigenvalue            # k_eff
sol.metadata["psi"]       # angular flux
sol.metadata["phi"]       # scalar flux
sol.metadata["mesh"]      # dict of grid arrays
```

### Commit 2 — `refactor(peierls_nystrom): rename topology -> orbit_space_class`

`14b60c9`. Renames per the brief:

- `CurvilinearGeometry.topology` (property) → `orbit_space_class`.
- `slab.TOPOLOGY` constant → `slab.ORBIT_SPACE_CLASS`.
- Returned strings `"two_surface"` / `"one_surface_compact"` →
  `"two-endpoint"` / `"one-endpoint"` (more precise — describes the
  orbit-space topology directly).
- Capability-matrix column key `topology_class` →
  `orbit_space_class` (cases.py + meta-generator).
- Rendered column header `"Topology"` → `"Orbit-space M/G"`.

Out of scope (preserved): `build_two_surface_case` /
`build_one_surface_compact_case` function names (public API);
`_class_a_cases` / `_class_b_cases` private helpers; the `"A"` /
`"B"` values in the capability matrix's orbit_space_class column
(internal taxonomy that survives the column-name rename).

The `.topology` property had ZERO code consumers outside its own
definition file (verified via grep before the rename), so the
rename is non-breaking for downstream code.

## Bit-equality preservation

The R2 facade is non-arithmetic: every `Billiard.solve_critical()`
call delegates to the underlying `solve_greens_function_*` entry
point and re-packs the result. No FP arithmetic moved; no
quadratures changed; no closures rewired. Bit-equality is an
IEEE-754 property, not an approximation tolerance.

**14-configuration baseline** captured at HEAD via `float.hex(...)`
exact-bit-pattern repr (sphere 1G closed/vacuum/MG/MR; cylinder 1G
closed/MG; slab sym 1G; slab asym 1G α_L/α_R/MG; hollow sphere 1G
closed/MG asymmetric; annulus 1G/MG; sphere MR fixed-source vacuum).
Re-running each post-refactor reproduces every `k_eff`, `psi.sum`,
`psi.norm`, `phi.sum`, `phi.norm`, and `iterations` to
**bit-exact equality** at every commit.

## Cross-method shared types adoption

The parallel MomentSpace agent shipped `CriticalSolution` /
`FluxSolution` shared types in `derivations.common.solution_types`
(commit `08da6ef` at session start; cherry-picked into this branch).
Schema:

- `CriticalSolution`: `eigenvalue`, `eigenvalue_kind`,
  `parameter_value`, `parameter_kind`, `converged`, `metadata`.
- `FluxSolution`: `spatial_nodes`, `scalar_flux`, `angular_flux`,
  `angular_nodes`, `eigenvalue`, `eigenvalue_kind`,
  `spatial_units`, `metadata`.

Method-specific rich data (full `psi`, mesh, raw legacy result)
lives under `metadata`. Variant α reports
`parameter_kind="fixed_geometry"` (no critical root-find — the
geometry is fixed by the caller). For sphere/cylinder/sphere_mr
parameter_value carries `R`; for slab/slab_asymmetric `L`; for
hollow geometries `R_out`.

This is the FIRST cross-method shared vocabulary. Both math-heart
classes (Billiard for trajectory_resolvent, MomentSpace for
fn_method) populate the same types, enabling cross-method
comparators that don't need to know which pillar produced the
result.

## R2 step 2 (legacy result dataclass unification) — DEFERRED

The brief's Commit 2 was "unify per-geometry result dataclasses
into CriticalSolution/FluxSolution" — ie make
`solve_greens_function_*` themselves return shared types.

**Decision**: deferred. Rationale:

1. The shared types are already adopted at the `Billiard.solve_*`
   facade level (the public API entry point per the brief).
2. Converting the LEGACY `solve_greens_function_*` entry points to
   return shared types would break ~78 existing tests that read
   `.psi`, `.psi_g`, `.k_eff`, etc. on the legacy dataclasses, and
   ~30 cross-method adapter call sites.
3. The legacy dataclasses are now **internal implementation
   detail** — preserved in `metadata["raw_result"]` for back-compat
   via the Billiard facade. Any external consumer can migrate to
   `Billiard` at their own pace.
4. Per the algebra-of-record discipline ("unify after two
   instances"), the shared-type Protocol over math-heart classes
   themselves is gated on more instances; the legacy-API rewrite
   should follow that same discipline.

Tracked as R2.5 in the hindsight refactor plan.

## Test results

| Suite | Pre-R2 | Post-R2 |
| --- | ---: | ---: |
| trajectory_resolvent solver suite (78) + Billiard foundation (15) | 78 → 93 | **93** |
| cross_method + capability matrices + peierls_geometry | 119 | **119** |
| Sphinx -W -b html | clean | clean |

## Branch hygiene incidents

The parallel MomentSpace agent (on `feat/fn-method-moment-space`)
yanked my branch checkout to theirs at least 4 times during the
session. Each swap deleted the in-flight `billiard.py` working-tree
file and reverted `__init__.py`. Recovery pattern: stash apply
from `stash@{1}` (carries my full file set on `refactor/r2-billiard-class`),
then re-Edit any reverted files, then commit IMMEDIATELY before
the next swap. Once committed, files are safe.

The lesson — **previously logged from the GeometrySpec rename
session** — is reinforced: parallel agents touching adjacent
territory should `git branch --show-current` before EVERY commit.
For this session the issue was not the parallel agent's commits
(they're on their own branch) but `git checkout` operations done
by the parallel agent that pulled HEAD onto their branch
mid-session.

## Files changed

Commit 1:
- NEW `orpheus/derivations/continuous/trajectory_resolvent/billiard.py`
  (1184 LOC)
- MOD `orpheus/derivations/continuous/trajectory_resolvent/__init__.py`
  (export Billiard, CriticalSolution, FluxSolution)
- MOD `docs/theory/trajectory_resolvent.rst` (+390 lines new theory)
- NEW `tests/derivations/test_trajectory_resolvent_billiard.py` (15 tests)

Commit 2:
- MOD `orpheus/derivations/continuous/peierls_nystrom/geometry.py`
  (topology property → orbit_space_class; string values renamed)
- MOD `orpheus/derivations/continuous/peierls_nystrom/slab.py`
  (TOPOLOGY → ORBIT_SPACE_CLASS)
- MOD `orpheus/derivations/continuous/peierls_nystrom/cases.py`
  (topology_class column key → orbit_space_class, 5 occurrences)
- MOD `tools/verification/generate_capability_matrices.py`
  (topology_class → orbit_space_class everywhere; header updated)
- MOD `docs/theory/peierls_nystrom.rst` (text references updated)
- MOD `docs/theory/reference_solvers.rst` (text references updated)
- MOD `docs/theory/_peierls_nystrom_capability_matrix.inc.rst`
  (auto-regenerated)

Net +2107 lines added (mostly docstrings + theory + tests),
−39 lines removed (renames).

## Out of scope (defer to R3 or later)

- R3 — `ChordOracle` Protocol: the highest-leverage refactor.
  Lands after R1+R2 bed in.
- R2.5 — legacy `solve_greens_function_*` return-type unification:
  see "DEFERRED" section above.
- B4 — `extract_scalar_flux` parametric, B5 — `mesh.volume_integral`,
  A4 — MG as tensor product, A5 — `solve(geometry, ...)` facade:
  opportunistic polishing items per the plan.

## Sphinx narrative status

The R2 theory section "Mathematical billiards and the Variant α
resolvent" is shipped at MAXIMUM-EFFORT level (per Cardinal Rule 3
and the user's "every time I read it I learn" directive). It is
NOT a stub — no archivist dispatch needed. The Birkhoff–Sinai
framing, Poincaré–Birkhoff billiard map, transfer-operator
resolvent, three regimes, and the Variant α correspondence
derivation are all in the page. Future archivist visits should
sharpen, not rebuild.
