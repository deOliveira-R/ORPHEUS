---
name: a5-boundary-flux-trace-space-rehome-audit
description: Pre-carve audit for field-role-typing A.5 — re-home BoundaryFlux onto the unified TraceSpace (FaceLayout moves field→space). Consumer inventory, gaps, retirement list, recommended minimal-churn carve.
metadata:
  type: project
---

# A.5 — BoundaryFlux → TraceSpace re-home (pre-carve audit)

Branch `refactor/field-role-typing`, worktree
`/Users/rodrigo/git/nuclear/ORPHEUS/.claude/worktrees/field-role-typing`.
A.5 makes `BoundaryFlux` a pure `values + space` Field whose `space` IS a
`TraceSpace` (built A.2/A.3, wired to mesh as `SNMesh.trace` in A.4). The
field-side `layout: FaceLayout` moves onto the space (`TraceSpace.layout`).

**Why:** completes the View-G boundary-space unification the TraceSpace
docstring (`trace_space.py:19-22`) already anticipates — retires the ad-hoc
`FunctionSpace("sn_boundary_flat")` build.

**How to apply:** drives a bit-identical turn-by-turn main-agent carve (per
[[no_method_implementer_for_surgical_carves]]). The `name` flip
`sn_boundary_flat → sn_trace` is the ONE behavioral identity change.

## Current structure (boundary_flux.py, 280 lines)
- `@dataclass(frozen,eq=False,kw_only)` over `Field`. Extra fields beyond
  Field's `values+space`: `layout: FaceLayout` (:100), `mesh: SNMesh` (:101).
- Ad-hoc `FunctionSpace("sn_boundary_flat",(total_size,))` built at
  `zeros_for_sn_mesh` (:207-210) and `from_face_arrays` (:274-277) — THE
  retirement target.
- `_check_partner` (:121-140): super class+space gate, then mesh `is` check
  (:126-130, keep) + layout `is`/`!=` fallback (:131-140, retire candidate).
- `face_view` (:144-173) + `face_views` (:175-183) read `self.layout.faces`.
- `__post_init__` (:105-117) reads `self.layout.total_size`.
- Algebra is Field-inherited via `replace(self, values=...)` — preserves
  layout/mesh/space; no arithmetic touches layout except `_check_partner`.

## Mapping + gaps
TraceSpace gives shape `(total_size,)` (:312), Euclidean weights=None,
`layout` (:262-264 Optional+compare=False), inherits FunctionSpace `__eq__`
on `(name,shape)`. Gaps:
1. `from_mesh_and_quadrature` needs `Quadrature` — but use the CACHED
   `mesh.trace` instead of rebuilding (it carries layout+shape).
2. `mesh.trace` is `Optional` (None for 2-D cylindrical, no SN sweep,
   geometry.py:411,863-881). NO production ctor path hits a trace-less mesh
   (all construction is inside SN solve). Recommend assert-non-None.
3. `face_view` must read `self.space.layout.faces` post-move.
4. `_check_partner` layout block: `space==space` gate already subsumes it
   once all BoundaryFluxes share the cached `mesh.trace` (same layout
   identity). Verify TestMeshBinding (test_boundary_flux.py:274-286) doesn't
   pin the "structurally-distinct layout raises" path before dropping.

## SNMesh wiring (the A.5 inputs)
- `SNMesh.trace` (geometry.py:863-881): cached `self._trace`, built in
  `_resolve_bcs` (:411-420) via
  `TraceSpace.from_mesh_and_quadrature(mesh, quad, boundary_face_layout)`.
  CANONICAL space source. Cached ⇒ stable layout identity across accesses.
- `SNMesh.boundary_face_layout` (geometry.py:884-944): REBUILDS a fresh
  FaceLayout each call (why `_check_partner` needs the is-then-!= fallback
  today). Sourcing space from cached `mesh.trace` fixes this.
- `zeros_for_sn_mesh` can do: `space=mesh.trace; values=np.zeros(space.shape);
  cls(values=,space=,mesh=)` — no separate `layout=` kwarg.

## Direct `BoundaryFlux(...)` ctor sites (pass space=/layout=): ONLY 2, both tests
- tests/transport/fields/test_boundary_flux.py:366 (space built :364)
- tests/sn/sweep/core/test_phase_c_gates.py:159 (space built :155-157)
Everything else uses factories.

## `.layout`-on-field reads that break if layout moves field→space
(If a read-through `@property layout` is kept, ALL of these survive verbatim —
the minimal-churn carve.)
- test_boundary_flux.py: :175,179,185,194,195,205,206,234,246,247,302-306,328
  + frozen test :433-437 (retire/retarget) + ctor :364-371
- test_native_matvec.py: :122,147,258,276,330,335,367
- test_invertible_operator.py: :577
- _test_helpers.py: :260
- test_l2_boundary_face_view.py: :145
NO production code reads `.layout` off a BoundaryFlux externally — all prod
external access is `.face_view(...)` + `.values` (API-stable). Prod breakage
confined to BoundaryFlux internals + the 2 factory space-builds.

## Production consumers (all stable except internals)
zeros_for_sn_mesh callers: scattering.py:1142, fission.py:313, solver.py:
{203 via zeros_boundary_flux, 668,1230,1399}, geometry.py:{821,858}.
sweep.py face_view sites (impl-internal change only): 471,472,487 (1-D),
840-843 (2-D seed), 934-937 (2-D writeback).
timed_full_field.py to_flat:418 / from_flat:457-458 use .values/replace — safe.

## Retirement list
1. `FunctionSpace("sn_boundary_flat")` builds: boundary_flux.py:207-210,
   274-277 → `space=mesh.trace`.
2. field-side `layout: FaceLayout` (:100) → space.layout (keep read-through
   `@property layout` for consumer stability = minimal churn).
3. `__post_init__` layout.total_size refs (:107,116).
4. face_view/face_views self.layout reads (:168,173,183).
5. `_check_partner` layout block (:131-140) — retire if no test pins it.
6. test ctor sites + `sn_boundary_flat` literals (test_boundary_flux.py:364,
   test_phase_c_gates.py:155-157) + frozen-layout test (test_boundary_flux.py:433).
7. docstring `sn_boundary_flat` mentions (boundary_flux.py:80-81,124).
NOT retired: boundary_face_layout, FaceLayout/FaceSlot, zeros_boundary_flux,
all factory call sites.

## Recommended carve
Keep `@property def layout(self): return self.space.layout`. Source
`space=mesh.trace` (assert non-None) in both factories. `name` flips to
`sn_trace`. Only 2 ctor sites + frozen-layout test + any `sn_boundary_flat`
name assertions need editing. Decide `_check_partner` layout-block drop after
checking TestMeshBinding.
