---
name: c4-bc-facelabel-producer-audit
description: "#220 C4 carve dependency audit: SNMesh bc_* named-attr surface → FaceLabel-derived producer. Verdict producer-side-only TRUE (2 generic prod consumers); axis.bc[label.endpoint] EXISTS; the missing piece is the FaceLabel→'xmin' name crosswalk; 1-D bc_ymin/ymax placeholders + bc_left/right aliases are production-orphans."
metadata:
  type: project
---

# C4 carve (#220) dependency audit — SNMesh BC producer surface (2026-06-11, worktree sn-nd-layout)

Fact: the carve IS producer-side-only in production. The only two production
readers of `SNMesh.bc_*` are already `getattr(mesh, f"bc_{face}")`-generic over
`trace.layout.faces`: `SNBoundaryOperator._face_laws`
(`orpheus/sn/boundary_operator.py:121-132`) and
`sweep_schedule._reflective_faces` (`orpheus/sn/sweep_schedule.py:208-219`).
`test_sweep_schedule_nd.py:101-119` already duck-fakes `bc_zmin`/`bc_zmax` via
SimpleNamespace — the f-string getattr surface is the consumer CONTRACT at d=3
today. Everything else reading `.bc_left`/`.bc_right` etc. on an SNMesh is
tests (≈6 files) or docs.

**Why:** issue #220 retires the named bc_xmin/.../bc_ymax + bc_left/bc_right
hybrid in `_resolve_bcs` (`orpheus/sn/geometry.py:394-490`) for a
FaceLabel-derived loop; nd_foundation §2.6 needs it before d=3.

**How to apply / load-bearing findings (verify line numbers via Nexus at pickup):**

1. `axis.bc[label.endpoint]` EXISTS TODAY — `AxisMesh.bc` (axis.py:202-204,
   keys = label_low/label_high = "min"/"max") and `RadialAxisMesh.bc`
   (axis.py:268-270, key = "outer"); `Axis1D` protocol declares it (:141).
   `face_labels()` iterates `axis.endpoints` which ARE the bc-dict keys, so
   `self.axes[label.axis_index].bc[label.endpoint]` is well-defined for every
   label. #220's proposed loop is REAL, not aspirational. The `axes` tuple is
   populated from `axes_from_legacy_mesh` (axis.py:406-493) which reads the
   Mesh1D/Mesh2D BC fields — Mesh fields stay; SNMesh attrs go.
2. THE missing piece: NO FaceLabel→"xmin"-style name crosswalk exists.
   `FaceLabel.__str__` = "face_0_min" (axis.py:92-93). The string world
   ("xmin"/"xmax"/"ymin"/"ymax") is wide: `trace.layout.faces`,
   `_FACE_NORMALS` (trace_space.py:167-172, NO z entries — separate d=3
   blocker), `SNMethodSpace.for_face(face=str)`, operator.py 1-D matvec
   hand-lists, loss_representation 1-D sweep hand-lists. Minimal-churn carve =
   string-keyed `SNMesh.bc` isomorphic to `layout.faces`, with names DERIVED
   from FaceLabel via one function using `AXIS_NAMES` (sweep_graph.py:97):
   `f"{AXIS_NAMES[i]}{'max' if ep=='outer' else ep}"`. The "outer"→"xmax"
   translation is currently hand-coded in TWO places (geometry.py:450, :948).
3. Latent d=3 bug to fix in-carve: `_resolve_one`'s
   `axis = "y" if face in ("ymin","ymax") else "x"` (geometry.py:523) maps
   "zmin"→"x". Replace with the axis from the FaceLabel.
4. Production-ORPHANS (retire outright, no shim): (a) 1-D degenerate
   `bc_ymin`/`bc_ymax` placeholders (geometry.py:463-483, minimal-method-space
   ReflectiveBoundary(axis="y") identity ops) — NO production reader (1-D
   `trace.layout.faces` = xmin/xmax only); pinned only by
   `test_snmesh_realizer_wiring.py` y-placeholder test. (b) the 2-D
   `bc_left`/`bc_right` aliases (geometry.py:489-490). (c) 1-D
   `bc_xmin`/`bc_xmax` aliasing — keep ONE canonical surface.
5. `boundary_face_layout` (geometry.py:900-961) has exactly ONE production
   caller: `_resolve_bcs`:430 → `TraceSpace.from_mesh_and_quadrature`.
   Its hand-listed reduced/curvature branches can become a face_labels loop
   (face shape = `(N, ng, *face_shape(label))` matches all three branches).
   `SNMesh.face_labels`/`face_shape`/`face_outflow_ordinates`/
   `n_unknowns_flat` have ZERO production consumers today (tests-only
   substrate) — this carve makes them load-bearing.
6. Test migration list (retirement = test migration):
   `test_boundary_conditions.py:67-129`,
   `test_snmesh_realizer_wiring.py` (many, incl. None-asserts :266-267,:308-309),
   `test_operator_block_role.py:186-190` (hand-lists 4 attr names),
   `test_phase_c_gates.py:666-672` (`bc_right.apply` monkeypatch),
   `test_unified_matvec_{cylinder,sphere}.py` (`sn_mesh.bc_right.apply`),
   plus the two already-generic getattr tests
   (`test_sn_boundary_operator.py:151,207`, `test_sweep_schedule.py:75`).
   Docs to sweep: boundary_conditions.rst (:959,:2211-2267),
   discrete_ordinates.rst (:159-210,:5377-5710), api/geometry.rst:67-88,
   index_convention.rst:759.
7. Pole invariant is already structural: `RadialAxisMesh.endpoints=("outer",)`
   ⇒ sphere face_labels = `(FaceLabel(0,"outer"),)` ⇒ the loop produces no
   pole entry; trace/layout agree ("xmax" only). Defaulting semantics to
   preserve: `None → BC("reflective")`.

Related: [[c36-dimensional-dispatch-audit]], [[sn-phantom-axis-rank-change-audit]].
