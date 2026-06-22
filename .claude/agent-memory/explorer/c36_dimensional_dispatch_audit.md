---
name: c36-dimensional-dispatch-audit
description: Post-S6 (#222) audit of every surviving is_1d/ndim dispatch site in orpheus/sn for the C3.6 N-D tail — the 3 genuine binary-where-ternary sites, the d=3 admission chain in dependency order, and the quadrature-is-NOT-the-blocker finding.
metadata:
  type: project
---

C3.6 pre-implementation audit (2026-06-11, worktree sn-nd-layout, post-#222
re-layering). Line numbers drift; the SHAPE is durable.

**The pre-S6 "5 not-is_1d matvec gates" item is RESOLVED** — S2/S6 collapsed
them; `StreamingOperator.apply/apply_transpose` route polymorphically through
`loss_representation.loss_action(_transpose)`. What survives is a different,
smaller set.

## The 3 genuine binary-where-ternary (b) sites
1. **`ScanMarch` supports/kernel mismatch** (loss_representation.py): `supports
   = is_1d or is_cartesian` ADMITS Cartesian d=3 (and the docstring claims
   "Cartesian any d"), but `sweep`/`loss_action` dispatch `if is_1d: ... else
   <d=2-hardcoded kernel>` — `_sweep_interior`/`_loss_action_interior` hard-
   unpack `ng, nx, ny = sig_t.shape` and `sx_eff, sy_eff = signs_eff`. At d=3
   default_for picks ScanMarch FIRST and dies with a raw unpack ValueError,
   not a principled deferral. Fix = either the recursive scan(x)∘march(y,z)
   or narrow supports to `is_1d or (is_cartesian and ndim <= 2)`.
2. **`_octant_sweep` truncates labels to 2-tuples** (sweep_schedule.py):
   builds `OctantLabel((sign_x, sign_y))` unconditionally — sign_z dropped
   with no ndim parameter. `_OctantWalk` does `signs[:ndim]` of a 2-tuple →
   2 faces for a 3-axis mesh → IndexError inside the (genuinely d-generic)
   FullFieldWavefront kernels. The schedule is the choke between the d-generic
   quadrature octants and the d-generic walk.
3. **`_OUT_FACE`/`_outgoing_faces` hand-list x/y** (sweep_schedule.py): no z
   entries → at d=3 a reflective z face would NEVER be assigned a G-S reflect
   group → silent wrong fixed point (the only SILENT-wrong site; everything
   else fails loud).

Plus silent-wrong-if-reached: `SNMesh.boundary_face_layout` hand-lists the
2-D 4-face layout under `reduced is None` — a d=3 mesh would get a layout
MISSING zmin/zmax.

## The d=3 admission chain (dependency order; first link blocks all)
1. `legacy_mesh_from_axes` raises for ≥3 axes (axis.py) — no Mesh3D dataclass
   (geometry/mesh.py has only Mesh1D/Mesh2D); `SNMesh.from_axes` round-trips
   through it. THE gate.
2. SNMesh ctor isinstance(Mesh1D)/else-Mesh2D metadata (nx/ny/dx/dy/mat_map).
3. `_setup_cartesian` builds only streaming_x/streaming_y; `streaming(axis)`
   body hand-lists `(streaming_x, streaming_y)[axis]` despite a d-generic
   docstring (IndexError at axis=2).
4. `boundary_face_layout` + `_resolve_bcs` hand-list xmin..ymax (bc fields are
   named attributes, NOT FaceLabel-keyed — the C4/#220 seam; consumers
   `SNBoundaryOperator._face_laws` + `_reflective_faces` are already d-generic
   `getattr(mesh, f"bc_{face}")` over `trace.layout.faces`).
5. Schedule layer (sites 2+3 above). 6. ScanMarch kernels (site 1).

**⭐ Quadrature is NOT a blocker**: `Quadrature` is genuinely 3-cosine
(nodes (N,3), mu_z, `axis_cosines(i)` d-generic, `octants` partitions by
sign of ALL components → 8 full octants for level_symmetric/lebedev).
The sweep_graph.py `_FrontierPlan` docstring's "(no 3-D quadrature yet)" is
MISLEADING — what's missing is the 3-axis MESH, not the quadrature.

## What IS d-generic and d=3-PINNED (don't re-derive)
Graph layer fully pinned: `from_cartesian`/`walk_full`/`walk_windowed`/
`cell_kernel_batch`/`residual_kernel_batch` at d=3 — nd_admission B1–B7
(synthetic shapes, all 8 octants, apply↔solve round-trip) + window≡full
bit-id at shapes (3,2,3)/(4,3,2) BOTH directions incl. sheds. d-generic and
unpinned-at-d=3: `_OctantWalk`, `_sweep_scheduled`, `_SweepEmit`,
FullFieldWavefront's interior kernels (blocked by mesh+schedule). NOTHING
runs any representation's `.sweep` at d=3 (impossible: needs an SNMesh).

## Honest-as-is (leave alone)
`CumprodScan.supports` (chain needs total order); `MovingFrontierWindow`
d==2 supports (explicit reason string); ScanMarch/`_DAGWavefront`
`loss_action_transpose` deferral raises ("multi-D" wording); operator.py's 3
`curvature=="cartesian" and not is_1d` breadcrumbs in `_compute_LpC`/
`_compute_decomposition` (S6.3-re-pointed, say "multi-D") — EXCEPT
`_compute_LpC_transpose`'s message still says "2-D Cartesian adjoint"
(stale wording, behavior fine); solver.py `reduced is None` proxy gates
(`_maybe_window`, G-S selector, 1-D cache) — honest while reduced⟺is_1d
coincidence holds.
