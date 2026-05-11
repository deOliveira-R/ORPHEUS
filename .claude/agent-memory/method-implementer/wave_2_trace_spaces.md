---
name: Wave 2 — Trace function spaces with per-face inflow/outflow masks
description: NEW orpheus/numerics/trace_space.py — TraceSpace ABC + InflowTraceSpace + OutflowTraceSpace frozen dataclasses subclassing FunctionSpace, carrying per-face directional masks built from mesh face normals × quadrature direction cosines.
type: project
---

# Wave 2 closeout — trace function spaces with per-face inflow/outflow masks

**Branch**: `feature/numerics-trace-spaces` (off `refactor/sn-operator-algebra`)
**Date**: 2026-05-11
**Plan**: `/home/vscode/.claude/plans/transient-giggling-cake.md` Wave 2

## Deliverables

- **NEW** `/workspaces/ORPHEUS/orpheus/numerics/trace_space.py` (~370 lines).
  - `TraceSpace(FunctionSpace, ABC)` — abstract dataclass.
  - `InflowTraceSpace(TraceSpace)` — carries `inflow_mask: NDArray[bool]
    (n_faces, n_ordinates)` + `face_names: tuple[str, ...]`.
  - `OutflowTraceSpace(TraceSpace)` — carries `outflow_mask` +
    `face_names` (symmetric).
  - `from_mesh_and_quadrature(mesh, quadrature, faces=None, ng=1)`
    classmethod on both.
  - `inflow_indices_for_face(face)` / `outflow_indices_for_face(face)`
    methods returning 1-D int index arrays (the API surface Wave 5's
    `SNBoundaryRealizer` will consume).
- **NEW** `/workspaces/ORPHEUS/tests/numerics/test_trace_space.py` —
  15 tests (12 from the brief + 3 bonus); 13 `@pytest.mark.l0`, 2
  `@pytest.mark.l1`.

## Design decisions taken

### Frozen-dataclass inheritance pattern (Option 1 from the brief)

`@dataclass(frozen=True) class TraceSpace(FunctionSpace, ABC)` works
in Python 3.12 without dataclass-machinery surprises. Both
`inflow_mask` and `face_names` are `field(default=None,
repr=False, compare=False)` (resp. `default=()`) so that:

- Equality / hash inherit `FunctionSpace`'s `(name, shape)`-only
  semantics. Two `InflowTraceSpace(N=8, ng=1)` built from different
  meshes still compare equal — Test 2 / `test_inflow_trace_space_equality_independent_of_mask`.
- `InflowTraceSpace` vs `OutflowTraceSpace` distinguishable by
  `name=trace_inflow` vs `trace_outflow` — Test 2 /
  `test_inflow_outflow_distinguishable_by_name`.
- Frozen mutation raising `FrozenInstanceError` still works — Test 1
  / `test_inflow_trace_space_constructible_and_frozen`.

Alternative considered: wrapping `FunctionSpace` instead of
inheriting. Rejected — would force every consumer to access
`.space` indirection. The dataclass inheritance is cleaner once
the `compare=False, repr=False` fields handle the
identity-preservation issue.

### Tangential ordinate handling — strict-inequality predicate

Predicate: inflow iff `sign * mu[axis] < -eps`, outflow iff
`> +eps` with `_TANGENTIAL_EPS = 1e-12`. Tangential ordinates
(|Ω · n_f| ≤ ε, e.g. Lebedev's `(0, 0, ±1)` on a face
perpendicular to z-axis) are excluded from BOTH masks. This is
the correct semantics — they neither flow in nor out, they graze
the face. Test 4 explicitly verifies this for Lebedev 11 on face
"xmin" (where the `(0, 0, ±1)` ordinates have Ω · n = 0).

### Face-normal table — name-based lookup

Two module-level dicts `_MESH1D_FACE_NORMALS` / `_MESH2D_FACE_NORMALS`
map face name → `(axis_index, sign)`. Face names match the BC
field naming convention (`bc_left` / `bc_right` for Mesh1D;
`bc_xmin` / `bc_xmax` / `bc_ymin` / `bc_ymax` for Mesh2D) — face
names dropping the `bc_` prefix: `"left"`, `"right"`, `"xmin"`,
`"xmax"`, `"ymin"`, `"ymax"`.

### Curvilinear deferral message — grep-able marker

The NotImplementedError carries the exact text the user's earlier
directive requested:
```
"InflowTraceSpace.from_mesh_and_quadrature for curvilinear
Mesh1D (coord=SPHERICAL/CYLINDRICAL) is deferred until a
curvilinear Krylov consumer arrives.
See plan Wave 2 (transient-giggling-cake.md)."
```
Same deferral message variant for Mesh2D with CYLINDRICAL coord.
Future consumer can grep for "transient-giggling-cake" or
"curvilinear Krylov consumer".

## Test counts

```
tests/numerics/test_trace_space.py: 15 tests
  13 @pytest.mark.l0
   2 @pytest.mark.l1
```

All 15 pass in 0.45 s.

## Regression suites

- `tests/numerics/` — **447 pass** (was 432 pre-Wave-2; +15 from
  test_trace_space.py).
- `tests/geometry/` — **177 pass**.
- `tests/sn/test_quadrature.py` — **49 pass** (focused regression;
  full `tests/sn/` was killed at 38% after 30+ min — pre-existing
  intrinsic SN-suite heaviness, not a trace_space regression. No
  module under `tests/sn/` imports `trace_space`; pure addition.)

Combined `tests/numerics/ tests/geometry/` — **624 pass, 1 warning,
0 failed, 0 skipped** in 0.89 s.

## Bifurcation discipline (per algebra-of-record)

Pure type construction + mask-building, no SymPy reference needed
— the math is `sign(Ω · n)` predicate evaluation, structurally
trivial. **Branch-2 production only**, no Branch-1 SymPy module.
Documented in the file's module docstring under the "Geometric
convention" section.

## Line ranges in `orpheus/numerics/trace_space.py`

- 1–95: module docstring (geometric convention, design rationale).
- 99–105: imports + `_TANGENTIAL_EPS = 1e-12`.
- 108–127: `_MESH1D_FACE_NORMALS` + `_MESH2D_FACE_NORMALS` lookup
  tables.
- 130–144: `_quadrature_axis()` helper (axis-index → mu_x/mu_y/mu_z
  with safe fallback for 1-D quadratures).
- 147–202: `_build_per_face_mask()` — core predicate construction,
  raises `NotImplementedError` for curvilinear meshes.
- 209–222: `TraceSpace(FunctionSpace, ABC)` — abstract base.
- 225–301: `InflowTraceSpace` — class body, factory,
  `inflow_indices_for_face`.
- 304–375: `OutflowTraceSpace` — class body, factory,
  `outflow_indices_for_face`.

## Augmentation response

Both pieces of lateral context from the main agent's augmentation
were honored:

1. Curvilinear NotImplementedError message embeds the
   "curvilinear Krylov consumer" / "transient-giggling-cake.md"
   marker as requested — grep-able for the future consumer.
2. `test_sn_2region_reflective_case_eigenvalue.py` was excluded as
   directed; the broader `tests/sn/` regression was attempted
   twice and timed out (heavy SN suite is intrinsic, not a Wave-2
   regression). Focused `tests/sn/test_quadrature.py` (49 tests,
   the most likely affected file by quadrature-touching changes)
   passes cleanly. Pure-addition + no `tests/sn/` import of
   `trace_space` means risk of SN-suite regression is structurally
   zero.

## Open items / next-wave hooks

- **Wave 3** `BoundaryTraceLaw.assert_source_lives_on_incoming_trace`
  — will consume `InflowTraceSpace` as the type tag for the source
  `q`.
- **Wave 5** `SNBoundaryRealizer.realize` — will iterate
  `for face in inflow_trace.face_names: indices =
  inflow_trace.inflow_indices_for_face(face); mask =
  IncomingOrdinateMaskTensor(indices, n_ordinates=N)`.
- **Wave 8** `SNMesh._resolve_bcs` — will hold the trace spaces as
  mesh-resident metadata constructed once at construction.
- **Curvilinear support** — when a curvilinear Krylov consumer
  arrives, replace the NotImplementedError with a
  geometry-specific predicate (radial-flow vs angular-flow per
  face on Mesh1D-spherical; r-z faces on Mesh2D-cylindrical).

## No Sphinx narrative this wave

Following the Wave-0 / Wave-1 precedent: no Sphinx narrative
shipped until the realizer (Wave 5) lands the unified theory page
that ties the BC primitives + trace spaces + realizer together.
No archivist dispatch.
