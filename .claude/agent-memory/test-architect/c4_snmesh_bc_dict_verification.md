---
name: c4-snmesh-bc-dict-verification
description: C4 carve V&V plan (#220, N-D campaign) — make SNMesh boundary inventory dim-generic. NEW face_name(FaceLabel)->str crosswalk ("min"→xmin/"max"→xmax/"outer"→xmax/"min" on axis2→zmin) single-sourced near FaceLabel; SNMesh.bc becomes a string-face-keyed dict (keys == boundary_face_layout.faces, both derived from face_labels via crosswalk); _resolve_bcs collapses the 1D-curv/1D-slab/2D isinstance split into ONE loop over face_labels; _resolve_one takes the FaceLabel (kills the :523 "y if face in ymin/ymax else x" latent-d3 bug); boundary_face_layout's 3-branch body → ONE from_named_shapes loop; RETIRE bc_xmin/xmax/ymin/ymax/left/right instance attrs + the degenerate 1D y-placeholder block + 2D left/right aliases. Two prod consumers migrate getattr(mesh,f"bc_{face}")→mesh.bc[face] (boundary_operator._face_laws:130, sweep_schedule._reflective_faces:218). Gate design: bit-id via existing solver suites (BC realization is object construction, realizer plumbing unchanged) + NEW structural pins for the dict surface + Mode-9 MIXED-BC observability + crosswalk L0 + d=3 synthetic admission. Ordered by 2-commit landing.
metadata:
  type: project
---

C4 verification plan — **make the SNMesh boundary inventory dimension-generic**
(#220, the N-D-layout campaign). Worktree `worktree-sn-nd-layout`. This is the
proactive `test-architect` gate (operator-algebra carve crossing the
geometry↔BC-resolution↔sweep-schedule boundary). Extends
[[c3-6-honest-d-dispatch-verification]] (same #220, same AXIS_NAMES SoT, same
Mode-8 bare-assert hazard, same synthetic-d3 idiom), and
[[c3-dim-generic-sweep-dag-verification]] (the x↔y-swap moat, the from_named_shapes
ordering). The DAG-side (sweep_graph/sweep_schedule) is C3.6; THIS is the
mesh-side BC inventory.

## ⭐ Claim-layer declaration (vv §1.5) — gate this first

Every C4 gate is a **flux-shape / structural** claim. There is NO new eigenvalue
claim and NO new MMS. The carve is **object construction** (which dict key holds
which realized `_BoundBoundaryOperator`) — the realizer PLUMBING is UNCHANGED
(`_resolve_one` still builds the same `SNMethodSpace.for_face` → `SNBoundaryRealizer`
→ `_BoundBoundaryOperator` chain; only HOW the face is named/looped changes). So:

- **Pillar for the dict-surface pins** = closed-form structural (the face-name
  set and per-face shape are derivable by hand from the axis endpoints — first
  principles, structurally independent of the crosswalk module's own derivation,
  vv L11). The hand-oracle: `(axis_index, endpoint) → AXIS_NAMES[axis_index] +
  {min→min, max→max, outer→max}`.
- **Pillar for the value bit-identity** = inheritance from the EXISTING verified
  solver suites (closed-form `Q/Σ_t`, `k_inf`, the affine sha256 goldens). The
  carve introduces ZERO numerical change; the value gates are "stays green /
  stays byte-identical".
- The eigenvalue safety rides transitively through `test_keff_2d` + the affine
  goldens staying green — NOT a new eigenvalue claim.

NO golden REGENERATION is expected (the carve is byte-identical by construction;
the layout order is preserved — verified below). If any golden drifts, STOP —
it's a bug, not a refactor.

## ⭐ The byte-identity proof (VERIFIED this session — load-bearing)

The `boundary_face_layout` rewrite (carve item 4) is byte-identical IFF the
crosswalk preserves the existing face ORDER (offset assignment is iteration-order
in `FaceLayout.from_named_shapes`). Verified live (`SNMesh` probe):

```
2D face_labels order:  [(0,'min'),(0,'max'),(1,'min'),(1,'max')]
crosswalk →            ['xmin','xmax','ymin','ymax']   == current 2D layout order ✓
2D face_shape((0,'min')) = (ny,)=(2,)  → (N,ng,2)=(N,ng,ny)  byte-id ✓
2D face_shape((1,'min')) = (nx,)=(3,)  → (N,ng,3)=(N,ng,nx)  byte-id ✓
sphere face_labels = [(0,'outer')] → ['xmax']  (single entry) == current ✓
  face_shape((0,'outer')) = ()  → (N,ng)  byte-id ✓
slab face_labels = [(0,'min'),(0,'max')] → ['xmin','xmax']  == current ✓
```

So `FaceLayout.from_named_shapes([(face_name(l), (N, ng, *self.face_shape(l)))
for l in self.face_labels])` reproduces the current layout BYTE-FOR-BYTE
(faces, shapes, offsets, total_size). `_OUTWARD_ENDPOINTS` (axis.py:347) already
includes `"outer"`, so the crosswalk's `outer→max` is consistent with the
outflow convention.

## ⭐ The :523 latent d=3 bug — the carve's correctness WHY (DEMONSTRATED reasoning)

`_resolve_one` (geometry.py:522-524) currently does:
```python
axis = "y" if face in ("ymin", "ymax") else "x"   # ← latent d=3 bug
law = ReflectiveBoundary(axis=axis, albedo=1.0)
```
A `"zmin"`/`"zmax"` face (d=3) falls into the `else` and gets `axis="x"` — a
reflective z-face would build the **x-axis reflection permutation**
(`quad.reflection_index("x")`) instead of z. This is the angular-redistribution
analog of the C3.6 silent-z-face-shed hole: a WRONG reflection partner at the
z-boundary → wrong inflow seed → wrong fixed point (silent; global balance
telescopes). The carve passes the `FaceLabel` so `axis = AXIS_NAMES[label.axis_index]`
— d-generic, correct at d=3 by construction. THIS is the gate-justifying surface
(G-cross + F2 below). NO ERR entry: d=3 is unconstructible today (Mesh3D is C5),
so no production bug shipped — closed-by-construction before d=3 exists.

## ⭐ Blast-radius audit (VERIFIED this session — narrows the consumer set)

LIVE production consumers of `sn_mesh.bc_<name>` ATTRIBUTE ACCESS (the carve's
retirement target): EXACTLY TWO, both the task-named:
- `boundary_operator.py:130` `getattr(self.sn_mesh, f"bc_{face}")` → `mesh.bc[face]`
- `sweep_schedule.py:218` `getattr(sn_mesh, f"bc_{face}") == "reflective"` → `mesh.bc[face] == "reflective"`
  (NB: relies on `_BoundBoundaryOperator.__eq__(str)` returning `self.kind==other`,
  `_bound_compat.py:128-131` — PRESERVE; the dict values are STILL
  `_BoundBoundaryOperator`, so `== "reflective"` still works.)

EVERY OTHER `.bc_left`/`.bc_xmin`/etc read in production is on the **INPUT
geometry mesh** (`Mesh1D.bc_left`/`Mesh2D.bc_xmin` — the user-facing input
convention), NOT on `SNMesh`: `sn/solver.py:72,75` (`mesh.bc_left`,
`getattr(mesh,"bc_xmin")`), `mc/solver.py:176`, `moc/geometry.py:302`,
`cp/solver.py:196` — these are `Mesh1D.bc_right` reads, UNTOUCHED by C4.
`operator.py:1321` ("StreamingOperator reads `sn_mesh.bc_*` directly") is a STALE
DOCSTRING — `operator.py` has ZERO live `sn_mesh.bc_*` attribute access (grepped);
the BC consumption is via `SNBoundaryOperator`/`sweep_schedule`. The docstring
wants a wording sweep but is not a code consumer.
`pole_angular_closure.py:1197` is a DOCSTRING (`if sn_mesh.bc_left is None`
describing a historical branch) — not live.

⟹ The retirement migration is **2 production lines + the test files**. The
docstrings (operator.py:1321, geometry.py:145, _bound_compat.py:87, vacuum.py:56)
are a wording-sweep review item, not code.

## 1. THE GATE LIST — ordered by 2-commit landing

The carve lands as TWO commits:
- **Commit (1)** = `face_name` crosswalk + `boundary_face_layout` loop rewrite
  (the LAYOUT-DERIVATION arm; numerically inert — the layout is byte-id).
- **Commit (2)** = `bc`-dict resolution loop + `_resolve_one(FaceLabel)` +
  consumer migration + RETIREMENT of the legacy attrs/placeholder/aliases.

### Commit (1) — crosswalk + layout derivation

| Gate | Level | File (new/ext) | Marker | Pins |
| ---- | ----- | -------------- | ------ | ---- |
| **G1.1** ⭐ crosswalk L0 exhaustive table (NEW) | foundation | NEW `tests/sn/primitives/test_face_name_crosswalk.py` | `foundation` | `face_name(FaceLabel(0,'min'))=="xmin"`, `(0,'max')="xmax"`, `(1,'min')="ymin"`, `(1,'max')="ymax"`, `(2,'min')="zmin"`, `(2,'max')="zmax"`, **`(0,'outer')="xmax"`** (curvilinear), `(2,'outer')="zmax"` (synthetic). d∈{1,2,3}×{min,max}+outer. PURE function, no mesh — d=3 admission is FREE here. |
| **G1.2** ⭐ crosswalk fail-loud negative (NEW) | foundation | SAME file | `foundation` | `face_name(FaceLabel(0,'bogus'))` RAISES (ValueError/KeyError) — the "fail loud on non-canonical endpoint" requirement. vv L11 §11: a validator needs BOTH a positive AND a negative test. |
| **G1.3** ⭐ layout byte-identity slab/sphere/cyl/2D (EXT — the EXISTING L0 suite) | foundation | EXT `tests/sn/primitives/test_boundary_face_layout.py` (UNCHANGED — it already pins faces/shapes/offsets/total_size) | `foundation` (file `pytestmark`) | the 6 existing tests (`test_slab_layout_has_xmin_and_xmax`, `test_sphere_layout_has_only_xmax`, `test_2d_layout_has_four_faces`, `test_2d_layout_excludes_interior_cells`, `test_layout_is_idempotent_property`) STAY GREEN after the loop rewrite — byte-id by construction. THE primary layout gate. ⚠ Mode-8: these are BARE asserts → run this file WITHOUT `-O` once in the (1) acceptance. |
| **G1.4** layout OFFSET order pin (NEW, sharpens G1.3 for the x↔y arm) | foundation | EXT same L0 suite OR crosswalk file | `foundation` | `list(m.boundary_face_layout.faces) == ["xmin","xmax","ymin","ymax"]` for 2D AND `faces["ymin"].offset == 2*N*ng*ny` (the y-block starts after BOTH x-faces). `np.testing`/`pytest.fail` (NOT bare assert). Pins the crosswalk did NOT reorder (a swap would still total-size-match but mis-offset — `from_named_shapes` would catch a gap, but an x↔y swap of EQUAL-shape faces would NOT; here nx≠ny so shapes differ ⟹ a swap trips `FaceLayout.__post_init__`, but pin it explicitly anyway). |
| **G1.5** A2D-1 + affine sha256 UNTOUCHED (regression) | l2/l3 | EXT `test_affine_carve_bit_identity.py` + `test_streaming_operator.py::TestT4dApply2DCartesianSourceHashPin` | (existing `l2`/`l3`) | the converged-flux sha256 goldens + the streaming-apply source-hash stay BYTE-identical — commit (1) is layout-derivation only, no kernel/no `_apply_2d_cartesian` source touch. **NO regeneration.** |

### Commit (2) — bc-dict + _resolve_one(FaceLabel) + consumer migration + retirement

| Gate | Level | File (new/ext) | Marker | Pins |
| ---- | ----- | -------------- | ------ | ---- |
| **G2.1** ⭐ `set(mesh.bc) == set(boundary_face_layout.faces)` (NEW) | foundation | NEW `tests/sn/primitives/test_snmesh_bc_dict.py` | `foundation` | THE acceptance invariant: slab → `{xmin,xmax}` (2), 2D → `{xmin,xmax,ymin,ymax}` (4), sphere → `{xmax}` (1), cylinder → `{xmax}` (1). And `set(mesh.bc) == set(mesh.boundary_face_layout.faces)` AND `== {face_name(l) for l in mesh.face_labels}` — all three derived from `face_labels` ⟹ structurally coherent. `np.testing`/`pytest.fail`. |
| **G2.2** ⭐ MIXED-BC per-face identity OBSERVABLE (NEW, Mode-9) | foundation | SAME file | `foundation` | the degenerate-blindness closer. 1-D slab `bc_left=reflective, bc_right=vacuum` → `mesh.bc["xmin"]=="reflective"` AND `mesh.bc["xmax"]=="vacuum"` (NOT both-reflective — a swap of the loop's `label.endpoint` would put vacuum on xmin). 2-D `bc_xmin=reflective,bc_xmax=vacuum,bc_ymin=vacuum,bc_ymax=reflective` → each of the 4 keys carries its DECLARED kind (ASYMMETRIC across BOTH axes ⟹ catches x↔y swap AND min↔max swap in the loop). `== "reflective"`/`"vacuum"` via `_BoundBoundaryOperator.__eq__`. |
| **G2.3** ⭐ the :523-fix observable — y-axis law from a y-face (NEW) | foundation | SAME file | `foundation` | the d=3-latent-bug catcher, OBSERVED at d=2. A 2-D `bc_ymin=reflective` realized op must carry the **y-reflection permutation** (`quad.reflection_index("y")`), NOT x. OBSERVE via `_angular_factor(mesh.bc["ymin"].inner).perm` (the `PermutationOperator.perm`, the existing `test_snmesh_realizer_wiring._angular_factor` helper) `== quad.reflection_index("y")` AND `!= quad.reflection_index("x")` (they DIFFER under a 2-D `level_symmetric` quad — assert the inequality so the pin can't pass vacuously if x==y perm). This pins `_resolve_one` derived `axis="y"` from `FaceLabel(1,*)`, not from the retired string-membership test. |
| **G2.4** ⭐ pole invariant — NO inner/pole entry (NEW) | foundation | SAME file | `foundation` | sphere AND cylinder `mesh.bc` has EXACTLY `{xmax}`; `"xmin" not in mesh.bc` AND `"face_0_outer" not in mesh.bc` (the raw FaceLabel str must NOT leak as a key). Preserves the RadialAxisMesh.endpoints=("outer",) ⟹ one entry structural invariant. |
| **G2.5** retired-attr AttributeError (NEW negative) | foundation | SAME file | `foundation` | `with pytest.raises(AttributeError): mesh.bc_left` (and `.bc_xmin`, `.bc_ymin` on a slab, `.bc_right`). vv §e: one explicit pin that the retirement is REAL (no shim survives). `pytest.raises` (fires under -O — it's a context manager, not a bare assert). Worth it — cheap, and a silent surviving attr is exactly the "deprecation shim outlives its cycle" smell. |
| **G2.6** key-miss raises KeyError (NEW negative) | foundation | SAME file | `foundation` | `with pytest.raises(KeyError): slab_mesh.bc["ymin"]` (1-D has no y-face) AND `sphere_mesh.bc["xmin"]` (no pole face). Pins `mesh.bc` is a plain dict (fail-loud on a non-existent face), NOT a defaultdict masking a missing face. |
| **G2.7** ⭐ MIGRATE `test_boundary_conditions.py` literal asserts (EXT) | foundation/l1 | EXT `tests/sn/operators/test_boundary_conditions.py:60-129` | (existing) | `sn.bc_left == "reflective"` → `sn.bc["xmin"] == "reflective"` (×all rows); `test_2d_mesh_resolution` `sn.bc_xmin==…` → `sn.bc["xmin"]==…`. Convert the touched BARE asserts → `assert … , msg`? NO — they're in test modules ⟹ pytest's rewriter keeps them LIVE under -O (L26). But the carve TOUCHES them, so per [[c3-6...]] §F1 discipline, leave as `assert` ONLY if confirmed live; the string-compare ones ARE fine (rewriter handles them). Keep the existing `test_mixed_bcs`/`test_2d_mesh_resolution` — they ARE the Mode-9 anchors, just re-keyed. |
| **G2.8** ⭐ MIGRATE `test_snmesh_realizer_wiring.py` — y-placeholder + curvilinear None (EXT) | foundation/l1 | EXT `tests/sn/operators/test_snmesh_realizer_wiring.py` | (existing) | (a) the y-placeholder test (`:200-230`, `isinstance(sn.bc_ymin, _BoundBoundaryOperator)` etc) DIES with the placeholder block → REPLACE with `test_1d_bc_has_no_y_entries`: `set(sn.bc)=={"xmin","xmax"}` AND `"ymin" not in sn.bc` AND `"ymax" not in sn.bc` (the placeholder's removal is the POINT — 1-D y is no longer a fake face). (b) curvilinear `sn.bc_left is None`/`sn.bc_xmin is None` (`:266-267,308-309`) → `"xmin" not in sn.bc`. The realized-op assertions on the OUTER face (`bc_right.inner` → IncomingOrdinateMaskTensor for vacuum / PermutationOperator for reflective) MIGRATE to `sn.bc["xmax"].inner`. |
| **G2.9** MIGRATE `test_operator_block_role.py:186-190` hand-list (EXT) | foundation | EXT `tests/sn/operators/test_operator_block_role.py` | (existing) | the `for face in ("bc_left","bc_right","bc_xmin","bc_xmax"): getattr(sn,face)` 1-D-slab block-role loop → `for face in ("xmin","xmax"): bc = sn.bc[face]` (drop the left/right duplicates — they were aliases of xmin/xmax; the block-role claim is per UNIQUE face). |
| **G2.10** MIGRATE `test_phase_c_gates.py:666,672` monkeypatch (EXT) | foundation/l1 | EXT `tests/sn/sweep/core/test_phase_c_gates.py` | (existing) | `sn_mesh.bc_right.apply` → `sn_mesh.bc["xmax"].apply` (the sphere outer face); `patch.object(sn_mesh.bc_right, "apply", …)` → `patch.object(sn_mesh.bc["xmax"], "apply", …)`. |
| **G2.11** MIGRATE curvilinear matvec patches (EXT) | l1/l2 | EXT `test_unified_matvec_cylinder.py:79,131` + `test_unified_matvec_sphere.py:87` | (existing) | `sn_mesh.bc_right.apply` → `sn_mesh.bc["xmax"].apply`. ⚠ cyl-matvec is the standing #206 xfail batch — confirm these are SETUP lines (not the xfailed assertions); the re-key must not un-xfail #206. |
| **G2.12** MIGRATE getattr mirrors (EXT) | foundation | EXT `test_sn_boundary_operator.py:151,207` + `test_sweep_schedule.py:75` | (existing) | `getattr(sn, f"bc_{face}")` → `sn.bc[face]` (these ALREADY loop over `sn.trace.layout.faces` / `sn_mesh.bc` candidate faces — the cleanest migration; the loop var `face` is already the string name). |
| **G2.13** ⭐ end-to-end value bit-identity (the inheritance gate) | l1/l2 | EXT solver suites (drive/confirm green) | (existing `l1`/`l2`) | the MINIMAL SUFFICIENT numerical set (see §3) stays green / byte-identical: `test_affine_carve_bit_identity.py` (sha256, NO regen), `test_dd_regression.py` (SAFETY×conv_tol), `test_keff_2d.py` (eigenvalue safety), the fixed-source `Q/Σ_t` + `k_inf` closed-form suites. |
| **G2.14** d=3 synthetic bc-dict admission (NEW, optional but cheap) | foundation | SAME `test_snmesh_bc_dict.py` OR crosswalk file | `foundation` | the `_resolve_bcs` LOOP is d-generic by construction; PROVE the crosswalk+loop COMBINATORICS at d=3 WITHOUT Mesh3D via the PURE pieces: `[face_name(FaceLabel(a,ep)) for a,ep in [(0,'min'),(0,'max'),(1,'min'),(1,'max'),(2,'min'),(2,'max')]] == ["xmin","xmax","ymin","ymax","zmin","zmax"]`. (The full `_resolve_bcs` at d=3 needs axes with `.bc[endpoint]`; a 3-axis `AxisMesh` tuple duck-type drives the LOOP but `_resolve_one` calls `SNMethodSpace.for_face` which needs a real mesh → defer the FULL d=3 resolution to C5/Mesh3D. State this scope boundary.) |

## 2. THE :523-FIX OBSERVABLE — how to OBSERVE the y-axis law (deliverable item b)

The task asks: "a FaceLabel(1,*) must produce a y-axis law — how to observe?
law's axis attr". VERIFIED this session: the realized `ReflectiveBoundary` does
NOT retain its `axis` string on the `_BoundBoundaryOperator` (the `axis` is
consumed INSIDE `SNBoundaryRealizer.realize` at `boundary_realizer.py:185`
`perm = quad.reflection_index(law.axis)` → baked into a `PermutationOperator`).
So the observable is the **PERMUTATION**, not an `axis` attribute:

```python
from <wiring test> import _angular_factor   # extracts the PermutationOperator
ymin_perm = _angular_factor(mesh.bc["ymin"].inner).perm
np.testing.assert_array_equal(ymin_perm, quad.reflection_index("y"))
# AND prove it's NOT the x-perm (else the pin passes vacuously):
assert not np.array_equal(quad.reflection_index("x"), quad.reflection_index("y"))
```

Under a 2-D `level_symmetric(sn_order=4)` quad, `reflection_index("x")` ≠
`reflection_index("y")` (x flips μ_x, y flips μ_y — different ordinate
partners). This is the SHARPEST observable for the :523 axis-dispatch fix: a
y-face that built the x-perm (the retired `face in ("ymin","ymax")` logic, OR a
d=3 `else→x` fallthrough) FAILS the `assert_array_equal`. This is G2.3 above;
it is the d=2 PROXY for the d=3 latent bug (F2 closes the d=3 direction
structurally).

## 3. THE MINIMAL SUFFICIENT NUMERICAL SET (deliverable item a)

BC realization is OBJECT CONSTRUCTION (the realizer plumbing is unchanged). The
value-side bit-identity is INHERITED from the existing suites — the question is
which configs make BC HANDLING OBSERVABLE (a swap/wrong-realization would change
the converged flux). The minimal set that makes each BC-handling path observable:

| Config | Suite | What makes BC observable |
| ------ | ----- | ------------------------ |
| **1-D slab MIXED** (reflective\|vacuum) | `test_affine_carve_bit_identity.py::si_slab_2g_het` + a fixed-source `Q/Σ_t` slab | a reflective↔vacuum swap on xmin/xmax changes the flux profile (vacuum bleeds, reflective doesn't). 2G het ⟹ ≥2G (H1). |
| **2-D Cartesian MIXED** (asymmetric x/y BCs) | `test_affine_carve_bit_identity.py::si_2d_p1_aniso_het` + `test_keff_2d.py` (a mixed-BC case) | x↔y swap in the dict loop changes the converged 2-D flux; the affine sha256 is sub-ULP sharp. THE x↔y-swap catcher on the VALUE side. |
| **sphere/cyl OUTER** (vacuum + reflective) | `test_phase_c_gates.py` + `test_l1_standoff_slab_cylinder.py` + the sphere LS4 pole regression | the single outer-face realization + pole-no-BC invariant; a leaked pole-face entry would change the inner seed. |

⟹ MINIMAL SUFFICIENT VALUE SET = `test_affine_carve_bit_identity.py` (sha256,
the sharpest — slab-2G + 2D-aniso-het, both mixed-realistic) +
`test_dd_regression.py` (the broad SAFETY×conv_tol snapshot, non-square 2-D
cases) + `test_keff_2d.py` (eigenvalue safety) + the curvilinear pole regression
(`test_l1_standoff_slab_cylinder.py` / sphere LS4). The affine sha256 + DD
snapshot are the BYTE-identity goldens that gate "no numerical change". **NO
regeneration** (byte-id by construction; if they drift, it's a bug).

## 4. THE MIGRATION MAP — what each breaking pin BECOMES (deliverable item d)

| File:line | Current | Becomes |
| --------- | ------- | ------- |
| `test_boundary_conditions.py:67-129` | `sn.bc_left == "reflective"`, `sn.bc_xmin == "reflective"`, … | `sn.bc["xmin"] == "reflective"`, … (string-compare via `__eq__`, rewriter keeps live under -O). Keep `test_mixed_bcs`/`test_2d_mesh_resolution` — re-key only; they ARE the Mode-9 anchors. |
| `test_snmesh_realizer_wiring.py:200-230` (y-placeholder) | `isinstance(sn.bc_ymin, _BoundBoundaryOperator)` + perm-identity asserts | DELETE the placeholder pin; REPLACE with `test_1d_bc_has_no_y_entries`: `"ymin" not in sn.bc and "ymax" not in sn.bc and set(sn.bc)=={"xmin","xmax"}`. The placeholder's DEATH is the carve's POINT. |
| `test_snmesh_realizer_wiring.py:266-267,308-309` (curvilinear None) | `sn.bc_left is None`, `sn.bc_xmin is None` | `"xmin" not in sn.bc` (the pole has no dict entry — was None, now absent). The outer-face `bc_right.inner` realized-op asserts → `sn.bc["xmax"].inner`. |
| `test_operator_block_role.py:186-190` | `for face in ("bc_left","bc_right","bc_xmin","bc_xmax"): getattr(sn,face)` | `for face in ("xmin","xmax"): sn.bc[face]` (drop left/right alias duplicates — block-role is per-unique-face). |
| `test_phase_c_gates.py:666,672` | `sn_mesh.bc_right.apply`, `patch.object(sn_mesh.bc_right,…)` | `sn_mesh.bc["xmax"].apply`, `patch.object(sn_mesh.bc["xmax"],…)`. |
| `test_unified_matvec_cylinder.py:79,131` + `test_unified_matvec_sphere.py:87` | `sn_mesh.bc_right.apply` | `sn_mesh.bc["xmax"].apply` (SETUP lines; must NOT un-xfail #206 cyl-matvec). |
| `test_sn_boundary_operator.py:151,207` + `test_sweep_schedule.py:75` | `getattr(sn, f"bc_{face}")` | `sn.bc[face]` (loop var already the string name — cleanest). |
| (hypothetical `test_sweep_schedule_nd.py` d=3 duck-fake) | `SimpleNamespace(**{f"bc_{face}": …})` getattr-attrs | `SimpleNamespace(bc={face: … for face in faces})` — the d=3 duck contract becomes a `bc` DICT entry, not getattr-attrs. (FILE DOES NOT EXIST at HEAD — this is a C3.6-plan artifact; if it lands before C4, migrate; else N/A.) |

**Post-migration grep gate (audit deliverable):** `grep -rn '\.bc_left\b\|\.bc_right\b\|\.bc_xmin\b\|\.bc_xmax\b\|\.bc_ymin\b\|\.bc_ymax\b' tests/sn/` returns ZERO `SNMesh.bc_*` ATTRIBUTE reads (only `Mesh1D`/`Mesh2D` INPUT-mesh reads survive — those are the geometry input convention, untouched). And `grep -rn 'getattr([^,]*, *f"bc_{' tests/ orpheus/` returns ZERO (the two prod mirrors migrated; test mirrors migrated). Confirms full retirement ([[retirement-means-test-migration]] + aggressive-retirement).

## 5. FAILURE MODES THE GATES DO NOT CATCH — + cheapest closer (item c/e)

| # | Foreseen miss | Why blind | Cheapest closer |
| - | ------------- | --------- | --------------- |
| **F1** ⭐ Mode-8 inert layout/bc asserts | `test_boundary_face_layout.py` + `test_boundary_conditions.py` use BARE `assert` (45-pass baseline emits the `-O` ignore warning, CONFIRMED this session). A WRONG byte-id layout could pass vacuously under -O. The string-COMPARE asserts (`sn.bc["xmin"]=="reflective"`) ARE kept live by pytest's rewriter (L26) — but only IN TEST MODULES. | The NEW pins (G1.4/G2.1-2.6/G2.14) MUST use `np.testing.assert_*`/`pytest.fail`/`pytest.raises` (fire regardless). Run `test_boundary_face_layout.py` + `test_boundary_conditions.py` ONCE WITHOUT `-O` in each commit's acceptance to confirm the touched legacy asserts are live (they ARE, via rewriter — but confirm, don't assume). |
| **F2** ⭐ x↔y↔z swap in the crosswalk INVISIBLE on symmetric all-reflective | `face_name` indexes `AXIS_NAMES[axis_index]`. A transposed `AXIS_NAMES=("y","x","z")` produces `ymin` where `xmin` belongs — INVISIBLE on a SQUARE all-reflective box (a swapped name still maps to a reflective face) AND on G2.1 (set-equality is order-blind). | G2.2 (MIXED-BC, ASYMMETRIC across both axes, nx≠ny) + G2.3 (the y-PERM observable) + G1.4 (offset order with nx≠ny) close it. The asymmetry is LOAD-BEARING — a square symmetric box masks the swap (vv §H2 corner). |
| **F3** crosswalk `outer→max` but layout expects `xmax` only ON AXIS 0 | a curvilinear axis is ALWAYS axis 0 (radial). If a future code path put a radial axis at index ≥1, `face_name(FaceLabel(1,'outer'))="ymax"` — correct, but unexercised. | G1.1 includes `(2,'outer')="zmax"` (synthetic) to pin the `outer→max` rule is axis-INDEPENDENT (not hardcoded to axis 0). Free. |
| **F4** ⭐ the d=3 BC-resolution VALUE is NEVER tested (only the crosswalk/loop STRUCTURE) | no Mesh3D exists; `_resolve_one` needs a real `SNMethodSpace.for_face`. G2.14 pins the crosswalk COMBINATORICS at d=3, NOT the realized d=3 z-reflection. The :523-fix CORRECTNESS at d=3 rides the d=2 G2.3 proxy + F2 structural. | STATE the boundary in G2.14's docstring (vv Mode-7 declare-what-is-NOT-activated): "d=3 BC-dict pins the crosswalk+loop STRUCTURE; the realized d=3 z-face reflection VALUE is C5/Mesh3D, gated by an end-to-end FP-invariance + pole/standoff at that phase." Honest analog of the C3.6 ScanMarch d=3 deferral. No new gate — a scope boundary. |
| **F5** retired-attr returns a NEW unrelated attr instead of AttributeError | if the retirement leaves a `@property bc_left` shim (Pattern: deprecation shim outlives cycle), G2.5's `pytest.raises(AttributeError)` FAILS — which is the POINT (it's a NEGATIVE gate). But if someone adds `__getattr__` that returns None for unknown `bc_*`, `mesh.bc_left` returns None silently. | G2.5 already covers it (`pytest.raises(AttributeError)` FAILS if a `__getattr__`-None shim exists ⟹ caught). Plus the grep gate (§4) catches a surviving `@property bc_left` in production. |

## 6. VV TAGGING + DISCIPLINE

- **Markers ([[feedback-vv-tagging]]):** every NEW C4 pin carries NO equation
  `:label:` → `@pytest.mark.foundation` (the crosswalk is a pure string map; the
  bc-dict is a data-structure invariant; the layout is a factory output). NEVER
  `verifies(...)` on these. The MIGRATED solver-suite rows KEEP their existing
  `l1`/`l2`/`l3` + `verifies`/`catches` tags (re-key only, no marker change).
- **Mode-8 (`-O`-strip) — the LIVE hazard (CONFIRMED this session: the 45-pass
  affected-suite run emitted the bare-assert `-O` ignore warning).** EVERY new
  C4 pin uses `np.testing.assert_*`/`pytest.fail`/`pytest.raises`. The migrated
  STRING-compare asserts (`sn.bc["xmin"]=="reflective"`) are kept live by
  pytest's rewriter (L26) — but run the two legacy-bare-assert files WITHOUT -O
  once per commit to confirm.
- **Self-improvement (fires BEFORE delivery):** the :523 `axis="y" if face in
  (...) else "x"` latent d=3 bug is a NEW INSTANCE of the same silent-wrong
  reflective-BC class as C3.6's z-face-shed (vv Mode 9 / Sig 1, ALREADY tabled)
  — but a DISTINCT mechanism: wrong-reflection-PARTNER (z-face builds x-perm) vs
  C3.6's face-never-sheds. This is NOT a new failure-MODE row for vv-principles
  (it's Mode 9). NO ERR entry: d=3 unconstructible today ⟹ no production bug
  shipped ⟹ closed-by-construction before d=3 exists (no "log every caught bug"
  trigger). If a future Mesh3D build ever ships the wrong z-perm, THAT becomes a
  new ERR. NOT tagging G2.3/F2 with `catches(...)` — there is no ERR to catch
  (the d=2 y-perm path is CORRECT today via the string-membership test; the carve
  preserves d=2 correctness while making d=3 correct-by-construction).

## 7. PRE-READS (file:line, VERIFIED this session)

- `orpheus/sn/axis.py:66-93` (`FaceLabel(axis_index, endpoint)`; `__str__` =
  `face_{i}_{ep}` — the raw str that must NOT leak as a bc key), `:100-141`
  (`Axis1D` protocol — `endpoints`, `bc[endpoint]`), `:148-204` (`AxisMesh.bc`
  = `{label_low: bc_low, label_high: bc_high}`, endpoints `(min,max)`),
  `:207-270` (`RadialAxisMesh.bc` = `{label_outer: bc_outer}`, endpoints
  `(outer,)` — the pole invariant source), `:282-308` (`face_labels`/`face_shape`
  dim-agnostic ground truth), `:347` (`_OUTWARD_ENDPOINTS={"max","outer"}` —
  the `outer→max` consistency anchor).
- `orpheus/sn/sweep_graph.py:97` (`AXIS_NAMES=("x","y","z")` — the SoT the
  crosswalk reads; import direction TBD by main agent).
- `orpheus/sn/geometry.py:394-490` (`_resolve_bcs` — the 1D-curv/1D-slab/2D
  isinstance split to collapse; `:463-483` the y-placeholder block to RETIRE;
  `:489-490` the 2D left/right aliases), `:492-535` (`_resolve_one`; **`:522-524`
  the `axis="y" if face in (...) else "x"` latent-d3 dispatch to kill**),
  `:711-735` (`face_labels`/`face_shape` SNMesh accessors), `:900-961`
  (`boundary_face_layout` — the 3-branch body to collapse to ONE
  `from_named_shapes` loop; `:943` `reduced is not None` slab/curv discriminator).
- `orpheus/sn/boundary_operator.py:120-132` (`_face_laws` — prod consumer #1,
  `getattr(sn_mesh, f"bc_{face}")`).
- `orpheus/sn/sweep_schedule.py:208-219` (`_reflective_faces` — prod consumer #2,
  `getattr(sn_mesh, f"bc_{face}") == "reflective"`; relies on
  `_BoundBoundaryOperator.__eq__`).
- `orpheus/geometry/boundary/_bound_compat.py:72-134` (`_BoundBoundaryOperator`;
  `:106` `self.kind`; `:128-131` `__eq__(str)→kind==other` — PRESERVE; dict
  values stay `_BoundBoundaryOperator`).
- `orpheus/sn/boundary_realizer.py:185` (`perm = quad.reflection_index(law.axis)`
  — where the :523 `axis` is consumed → baked into a `PermutationOperator`; the
  observable for G2.3 is the PERM, not an `axis` attr).
- `orpheus/numerics/face_layout.py:108-200` (`FaceLayout`; `:133-155`
  `__post_init__` offset/gap validation; `:157-200` `from_named_shapes` —
  iteration-order = offset assignment, the byte-id ordering anchor).
- Test files to migrate: `test_boundary_conditions.py:60-129`,
  `test_snmesh_realizer_wiring.py:200-230,260-312`,
  `test_operator_block_role.py:180-190`, `test_phase_c_gates.py:660-678`,
  `test_unified_matvec_cylinder.py:79,131`, `test_unified_matvec_sphere.py:87`,
  `test_sn_boundary_operator.py:151,207`, `test_sweep_schedule.py:75`.
- EXISTING bit-id anchors (stay green/byte-id, NO regen):
  `tests/sn/primitives/test_boundary_face_layout.py` (the 6 L0 layout pins),
  `tests/sn/solve/test_affine_carve_bit_identity.py` (sha256),
  `tests/sn/regression/test_dd_regression.py` (SAFETY×conv_tol),
  `tests/sn/operators/test_streaming_operator.py::TestT4dApply2DCartesianSourceHashPin`
  (A2D-1), `test_keff_2d.py` (eigenvalue safety), `test_method_space.py`,
  sphere LS4 pole regression.
- Baseline this session (`-O`): `test_boundary_conditions.py +
  test_snmesh_realizer_wiring.py + test_boundary_face_layout.py +
  test_operator_block_role.py` = **45 passed** (clean; the `-O` bare-assert
  warning IS emitted ⟹ Mode-8 live).
