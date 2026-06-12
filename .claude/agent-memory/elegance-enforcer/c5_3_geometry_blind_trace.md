---
name: c5-3-geometry-blind-trace
description: C5.3 (#225) geometry-blind trace-space review — AXIS_NAMES relocate-down + 6-face derive + dead gate-only mesh param retire; PASS-with-nits (one forward CONCERN, one out-of-scope docs-rot)
metadata:
  type: project
---

C5.3 (#225, commit `e66025b` + docs follow-up `a9a0faf`, branch
worktree-sn-nd-layout) — the trace layer goes geometry-blind, the d=3 face
table un-blocks. Reviewed PASS-with-nits. Sibling to [[c5-1-axis-primary-snmesh]]
+ [[c5-2-phantom-shim-retirement]]. The C5 BOUNDARY/trace tail before 3-D
admission.

**The carve (4 moves, all sound):**
1. `AXIS_NAMES` RELOCATED DOWN `orpheus/sn/axis.py → orpheus/numerics/face_layout.py`
   (one literal def, face_layout.py:68; verified ONLY def). `sn.axis` re-exports
   via `from orpheus.numerics.face_layout import AXIS_NAMES  # noqa: E402` — ONE
   object, 3 import surfaces (trace_space DIRECT; sn.axis re-export consumed by
   sweep_schedule/geometry/loss_representation). NO twin-import hazard: module-level
   immutable tuple, never re-bound. face_layout imports ONLY stdlib+numpy (verified
   — genuine bottom of dep graph, no sn-ward, no cycle). RIGHT home: FaceLayout is
   the `"{axis}{min|max}"` string-world keeper; AXIS_NAMES is the axis half of that
   crosswalk — Pattern-7 definition-site colocation. PARTIALLY resolves C5.1 NIT 2
   (axis-name convention now single-sourced at numerics; the AxisCoord↔CoordSystem
   triplication is SEPARATE, still bounded).
2. `_FACE_NORMALS` 4-entry hand-list → comprehension over AXIS_NAMES
   `{f"{name}{suffix}":(axis,sign) for axis,name in enumerate(AXIS_NAMES) for
   suffix,sign in (("min",-1),("max",+1))}`. Produces the 4 OLD entries
   byte-identical + zmin/zmax. The old hand-list SILENTLY LACKED z (the pre-C5.3
   d=3 blocker — twin/transcription collapse). sha256 byte-id structurally
   guaranteed: `_build_omega_dot_n` iterates `faces` (layout) not `_FACE_NORMALS`,
   dict insertion order irrelevant to output; 4 existing (axis,sign) unchanged.
3. `_build_omega_dot_n(mesh,quad,faces) → (quad,faces)` — mesh param was GATE-ONLY
   + both gates DEAD: (a) the curvilinear-Mesh2D `NotImplementedError` UNREACHABLE
   (VERIFIED: `axes_from_legacy_mesh` axis.py:558-564 raises `NotImplementedError
   match="non-Cartesian"` for 2-D non-Cartesian DURING SNMesh `__init__` legacy
   adapter → no SNMesh with cylindrical Mesh2D can EXIST; `from_axes` native never
   builds Mesh2D); (b) the isinstance-TypeError carried no data. `from_mesh_and_
   quadrature → from_quadrature_and_layout` rename + bare-ctor error msg migrated
   (aggressive retirement, ZERO surviving refs to the renamed `TraceSpace` method).
4. `SNMesh._resolve_bcs` trace build now UNCONDITIONAL (gate excluded only the
   unconstructible); `SNMesh.trace` typed/doc'd ALWAYS-non-None. `SNMethodSpace.
   for_face` mesh param → `Optional=None` metadata; call site (geometry.py:562)
   still feeds `self.mesh` while the legacy adapter lives — fine.

**Test migration FAITHFUL + complete (44 green, ran -O):** `test_2d_cylindrical_
raises` (pinned the TRACE-level gate) RETIRED WITH its gate; refusal-pin RELOCATED
to the construction surface `test_2d_cylindrical_mesh_refused_at_construction`
(`match="non-Cartesian"` = the axes_from_legacy_mesh msg) — textbook retirement=
test-migration, invariant survives where the dispatch now lives. Unknown-face probe
"zmin"→"wmin" (z now KNOWN). NEW C5-G8..G11 (foundation-tier): G9 6-entry table
vs HAND-TRANSCRIBED expected dict (CORRECT mirror-not-import — deriving expected
from AXIS_NAMES would be vacuous self-compare; the hand-mirror catches a
comprehension bug); G10 zmax==+mu_z / zmin==−mu_z exact; G11 all-six distinct-axis
rows under level_symmetric (Mode-2 swap detector) WITH non-vacuity guard
(`pytest.fail` if cosine arrays coincide); G8 trace-builds-for-every-constructible.

**`a9a0faf` follow-up = my C5.2 CONCERN 1 CLOSED:** the 6 `mesh.nx, mesh.ny`
literal-attr docstring refs (angular_flux/scalar_flux/angular_residual/scalar_
residual/angular_source_sink/scalar_source_sink) → `*mesh.spatial_shape`.
DISCRIMINATION HELD EXACTLY (per C5.2 memory rule): literal-attr derivation prose
FIXED; bare `(ng,nx,ny)` VALUES layout-placeholder one line above LEFT untouched.

**NITS (neither blocks):**
1. (CONCERN, forward-looking) `_quadrature_axis` (PRE-EXISTING, predates C5.3)
   `getattr(quad,"mu_z",zeros)` zero-fallback: a z-face request against a 2-D
   quad silently returns Ω·n=0 (all-tangential) rather than failing. Cannot misfire
   TODAY (boundary_face_layout names z only for genuine 3-D mesh ⇒ 3-D quad). C5.3
   is the commit that newly makes z reachable through `_FACE_NORMALS`, so flag for
   the 3-D admission step: a layout/quad rank MISMATCH should fail loud, not zero-
   pad. Bug habitat = a future partial-3-D wiring that names z on a 2-D quad gets
   silent tangential instead of an error.
2. (out-of-scope, ARCHIVIST) `docs/theory/{boundary_conditions,operator_algebra,
   discrete_ordinates}.rst` carry ~6 stale `:meth:`InflowTraceSpace.from_mesh_and_
   quadrature`` refs. `InflowTraceSpace` is a LONG-RETIRED class (survives only in
   trace_space.py:17 history-mention) — these predate even the TraceSpace
   unification, NOT refs to the symbol C5.3 renamed. C5.3 deepens the rot slightly
   (the method name now exists on NO live class). Pre-existing docs-rot cluster;
   not a finding against this commit's prod code.
3. (cosmetic) `a9a0faf` left scalar_residual.py:78 / scalar_source_sink.py:97 as
   long single physical lines (`shape == (mesh.ng, *mesh.spatial_shape)``. Use...`).
   Pure wrap-width, no bug.

Comment honesty GOOD: the `min`=−axis/`max`=+axis convention comment (trace_space.py
:172) matches the comprehension data AND every face-name producer (loss_repr:335/344,
sweep_schedule:203). The "unreachable" / "carried no data" claims VERIFIED against
axes_from_legacy_mesh, not just asserted.
