---
name: issue-251-trace-space-widening
description: "#251 (Leg B of #247) DEEP file:line map of the TraceSpace structure + boundary-trace lifecycle, for widening the 2-D Cartesian LD boundary trace to carry the 2^{d-1} transverse face-slope moment. The trace is purely SCALAR-per-face today (NO moment axis, NO face_moment_tail hook); the cleanest widening lever is geometry.py's boundary_face_layout per-face slot shape (face_shape ⊗ 2^{d-1}). L17 inflow/outflow crosswalks + L20 construct/consume audit + SOLVE-vs-MATVEC site list. Branch refactor/sn-foundation-cleanup HEAD d9396a2."
metadata:
  type: project
---

# #251 — boundary-trace widening for the 2^{d-1} transverse face-slope (Leg B of #247)

**Tree**: main checkout, branch `refactor/sn-foundation-cleanup`, HEAD `d9396a2`
(Leg A of #247 just landed). NOT a worktree. `$CLAUDE_ENVIRONMENT` = host
(`.venv/bin/python`). All line numbers from Read/Grep of the working tree.
**Goal**: Leg A (#247) widened the BULK external slope-source; Leg B (#251) widens
the BOUNDARY TRACE so a moment-resolved manufactured inflow can carry the
`2^{d-1}` transverse (along-face) slope, and the sweep outflow can STORE those
`2^{d-1}` moments instead of collapsing to slot 0. Today the trace is
**scalar-per-face** end-to-end.

---

## 0. THE ONE STRUCTURAL FACT THAT SHAPES EVERYTHING (read first)

**The TraceSpace does NOT own the per-face shape — `FaceLayout` does, and the
`FaceLayout` is built by `SNMesh.boundary_face_layout` from `face_shape(label)`.**

The chain (storage):

```
SNMesh.boundary_face_layout (geometry.py:1004, body 1044-1050)
  → FaceLayout.from_named_shapes([(face_name, (N, ng, *face_shape(label))) ...])
       face_shape(label)  = geometry.py:809 → axis.py:349-360 = the codim-1 cell counts
  → FaceSlot.shape = (N, ng, *face_shape)          # face_layout.py:97, validated flat_size==prod(shape)
  → FaceLayout.total_size = Σ prod(slot.shape)     # face_layout.py:148
TraceSpace.from_quadrature_and_layout(quad, layout)  (geometry.py:516-519)
  → TraceSpace.shape = (layout.total_size,)         # trace_space.py:384  (FLAT whole-boundary)
  → omega_dot_n = (n_faces, n_ord)                  # trace_space.py:374 — PER-ORDINATE ONLY, no spatial
  → inner_product_weights = |Ω·n|⊙w_n broadcast over ALL trailing slot axes  # trace_space.py:285
```

So **there is NO `face_moment_tail` hook and NO `SpatialMomentSpace`-like factor on
the trace today.** `moment_layout.face_moment_tail` is NOT imported by
`trace_space.py` nor by `geometry.boundary_face_layout`. The trace is purely
scalar-per-face: each face slot is `(N, ng)` (1-D), `(N, ng, ny)` / `(N, ng, nx)`
(2-D edges). The `omega_dot_n` table and the metric are **per-ordinate only** —
both broadcast cleanly over any trailing axes (the metric explicitly so:
`face_w_axis0 = face_w.reshape((N,) + (1,)*(len(slot.shape)-1))` then
`broadcast_to(slot.shape)` at trace_space.py:285-286). **A moment axis added to the
slot shape rides the metric + omega_dot_n for free** — they classify/weight by
ordinate (axis 0), agnostic to trailing spatial/moment axes.

**⇒ The cleanest place to grow the `2^{d-1}` transverse-moment axis is the per-face
SLOT SHAPE in `geometry.py:boundary_face_layout`** (append the face moment tail to
`(N, ng, *face_shape(label))`). Everything downstream — `TraceSpace.shape`
(derives from `total_size`), the BoundaryFlux/BoundarySourceSink/BoundaryResidual/
BoundaryDisplacement validation (against `total_size` only, _bases.py:516), the
metric, `face_view` — accommodates it automatically because they are all
flat-shape-from-layout (validate only against `total_size`, never a hardcoded
`(N, ng)`).

**Scheme reachability (verified)**: `self.scheme` is set at `geometry.py:271-272`
BEFORE `self._trace` is built at `geometry.py:517`. So
`self.scheme.spatial_basis_per_axis` IS available to `boundary_face_layout`.
`boundary_face_layout` is currently geometry-blind/scheme-blind; widening it means
giving it ONE read of the scheme's `_n_face_moments = per_axis**(ndim-1)` to append
the tail. This is the central DESIGN QUESTION the method-implementer must answer:
does the moment width live on the layout (scheme-aware `boundary_face_layout`), or
does the trace stay scalar and the moments live ONLY in the interior cochain +
manufactured inflow (no persistent storage)? The persistent-storage path needs the
layout widened; a transient-only path (inflow read + outflow used immediately, never
stored on `mesh.trace`) does not — but Leg B's outflow-store requirement (§3) and
the prescribed-inflow producer (§5) both want a persistent moment trace.

**Negative control (DD/Step byte-identical)** — CONFIRMED:
`LinearDiscontinuous.spatial_basis_per_axis == 2` (linear_discontinuous.py:300);
the DD/Step base `spatial_basis_per_axis == 1` (scheme.py:815). For DD/Step:
`_n_face_moments = 1**(ndim-1) = 1` → `face_moment_tail(1) == ()`
(moment_layout.py:73) → NO trailing axis → slot stays `(N, ng, *face_shape)` →
trace byte-identical. For LD-2D: `_n_face_moments = 2**(2-1) = 2`.
`_n_face_moments` lives at `loss_representation.py:311-321`.

---

## 1. `TraceSpace` itself (`orpheus/numerics/spaces/trace_space.py`, read in full)

A frozen-dataclass `FunctionSpace` subclass. **It is whole-boundary FLAT** —
`shape == (layout.total_size,)` — with per-face access derived from the `layout`
slot. It is `compare`-by-`(name, shape)` only (`layout` + `omega_dot_n` are
`compare=False` leaf data, trace_space.py:327-332).

- **Construction**: `from_quadrature_and_layout(quadrature, layout)`
  (trace_space.py:341-388). Builds: `omega_dot_n = _build_omega_dot_n(...)`
  (`(n_faces, n_ord)`, trace_space.py:200-249), the partial-current metric
  `inner_product_weights = _build_trace_metric_weights(...)`
  (`(total_size,)`, trace_space.py:252-288), then `cls(name="sn_trace",
  shape=(total_size,), inner_product_weights=..., layout=layout,
  omega_dot_n=omega_dot_n)` (trace_space.py:382-388). The ONLY production caller
  is `geometry.py:517`.
- **Internal layout**: there is NO trace buffer owned by the space — the space is
  a SHAPE descriptor. The buffer lives on the boundary FIELDS
  (BoundaryFlux/BoundarySourceSink/...). Per-face data shape = the `FaceSlot.shape`
  from the layout (today `(N, ng, *face_shape)`). Per-face VIEW =
  `layout.faces[face].slice_view(flat_buf)` (face_layout.py:113-121,
  `flat_buf[offset:offset+flat_size].reshape(shape)` — memory-shared, no copy).
- **face_moment_tail / SpatialMomentSpace hook?** NONE. The trace has no
  awareness of spatial moments. The `omega_dot_n` and metric are per-ordinate and
  broadcast over trailing axes, so they are moment-axis-READY by construction, but
  nothing today appends a moment axis.
- **Directional selectors** (the inflow/outflow predicate the whole boundary
  lifecycle keys on): `inflow_indices_for_face(face)` (trace_space.py:407-414,
  `flatnonzero(row < -_TANGENTIAL_EPS)`) / `outflow_indices_for_face(face)`
  (trace_space.py:416-423). These index ORDINATES (axis 0) only — a moment axis on
  the slot does not perturb them (the `[inflow, :]`/`[..., target]` indexing on
  axis 0 stays valid for a `(N, ng, ...transverse..., 2^{d-1})` slot).
- **Cleanest place to grow the moment axis**: NOT in `trace_space.py` (it is
  shape-from-layout). The lever is `boundary_face_layout` (the slot shape) — see §0.
  If chosen, `trace_space.py` needs ZERO changes (the metric broadcast + selectors
  already handle trailing axes). This is the elegant outcome: the trace space is
  ALREADY moment-polymorphic; only its layout-supplier is scalar.

---

## 2. The INFLOW READ path

The sweep reads the GIVEN trace inflow at three layers (SOLVE and MATVEC share the
walk-level lift):

| Layer | Site | What it reads | Shape today | After widening |
|---|---|---|---|---|
| Operator seed (SOLVE) | `operator.py:1107-1115` | `boundary_buf.face_view(face)[:] = rhs.boundary.face_view(face)` — copies the inflow seed into a fresh `BoundaryFlux.zeros_on` | scalar `(N,ng,*face_shape)` whole-slot copy | AUTO-widens (whole-slot copy; grows with the layout) |
| Octant inflow read (SOLVE) | `loss_representation.py:712` (`inflow_of=boundary_flux.face_view`) → 651-654 builds the per-octant `inflow` tuple | scalar per face | AUTO-widens if slot grows; the values carry the moment axis |
| Octant inflow read (MATVEC) | `loss_representation.py:806` (`inflow_of=boundary.face_view`) → `_OctantWalk.loss_action` | scalar per face | same |
| Scalar→moment LIFT (BOTH) | `loss_representation.py:357-378` `_inflow_to_moments` (esp. **375-376**) | takes scalar `(*face,)`, zero-fills `(*face, 2^{d-1})`, slot-0=scalar, **transverse slopes ZEROED** | #251 makes this a no-op / identity IF the trace already carries the moment axis (the inflow tuple arrives moment-valued), OR rewrites it to PASS THROUGH a moment-valued inflow instead of zero-fill |
| Walk-level seed (FULL oracle) | `loss_representation.py:1180-1222` `_allocate_full_cochain` (esp. **1214-1220**) | seeds IN-edge from inflow; `n_face_moments==1` → direct, else slot-0=scalar transverse-zero | #251: seed ALL `2^{d-1}` from a moment-valued inflow |

**KEY**: `_inflow_to_moments` (`loss_representation.py:357-378`) is the boundary
twin of Leg A's `_lift_external_source_to_moments`. Today it INGESTS a scalar
inflow (because the trace stores scalars) and zero-fills the transverse slopes. The
#251 widening has two sub-shapes depending on the §0 design choice:
- **If the trace stores moments** (layout widened): the `inflow` tuple arrives
  already `(*face, 2^{d-1})`-valued; `_inflow_to_moments` becomes the IDENTITY (the
  `n == 1` early-return generalizes to "already moment-shaped → pass through"). The
  transverse slope is read straight off `face_view`.
- **If the trace stays scalar** (no persistent storage): `_inflow_to_moments` must
  instead EVALUATE the manufactured transverse moments from a moment-aware inflow
  source — a heavier change, and it loses the outflow-store path (§3). NOT
  recommended; the layout-widening path is cleaner.

The lift is called at exactly two production consumers (MovingFrontierWindow,
the 2-D production strategy): `_sweep_interior` (SOLVE) `loss_representation.py:1033`,
`_loss_action_interior` (MATVEC) `loss_representation.py:1120`. The
FullFieldWavefront oracle calls it via `_allocate_full_cochain` (the `inflow=...`
walk path; the per-direction wrappers at 1290-1295 / 1364-1368).

---

## 3. The OUTFLOW CAPTURE/STORE path (the deeper change)

The sweep CAPTURES the domain-edge outflow as a `2^{d-1}`-moment object (the
interior cochain is already moment-valued during the walk) and then **collapses it
to slot 0** before writing the scalar trace. The collapse sites — what
`capture = tuple(c[..., AVERAGE_MOMENT] for c in capture)` DROPS:

| Path | Strategy | Capture-collapse site | What is dropped |
|---|---|---|---|
| SOLVE | MovingFrontierWindow (prod) | `loss_representation.py:1068-1069` | the `2^{d-1}` transverse outflow moments; only slot 0 (average) survives |
| MATVEC | MovingFrontierWindow (prod) | `loss_representation.py:1140-1141` | same |
| SOLVE | FullFieldWavefront (oracle) | `loss_representation.py:1319-1320` | same |
| MATVEC | FullFieldWavefront (oracle) | `loss_representation.py:1386-1387` | same |

Each collapse is guarded `if n_face_moments > 1:` — so DD/Step (==1) never
collapses (no axis), byte-identical. **"Un-collapsing"** = delete/skip the collapse
so `capture` retains the `2^{d-1}` axis, AND the SHED target must accept it:

- SOLVE shed: `loss_representation.py:707-708`
  `boundary_flux.face_view(face)[oct_idx] = capture_a`. The `boundary_flux` slot
  must carry the moment axis (i.e. the trace/layout must be widened — §0) for this
  write to land the `2^{d-1}` outflow. With a scalar trace today, the collapse is
  MANDATORY (the write would shape-mismatch otherwise).
- MATVEC shed: `loss_representation.py:800-801`
  `streamed[face][oct_idx] = capture_a`, where `streamed` is allocated
  `zeros_like(boundary.face_view(face))` (`loss_representation.py:777-780`) — so it
  inherits the trace slot shape; widen the trace and `streamed` auto-widens.
- MATVEC boundary-residual `B`-block emit: `_OctantWalk.loss_action`
  `loss_representation.py:812-823` — reads `given = boundary.face_view(face)`
  (815), writes `out_boundary.face_view(face)[out_idx] = streamed[face][out_idx] -
  given[out_idx]` (819-820, the outflow defect) and `[in_idx] = given[in_idx]`
  (823, the inflow identity). All whole-slot/ordinate-indexed → auto-widens with the
  layout.

**⇒ The outflow-store half REQUIRES the trace/layout widening (§0)**; the inflow
half could in principle be transient, but storing the outflow moments needs a
moment-shaped slot. This is why §0 recommends the layout-widening path: it serves
BOTH the inflow read and the outflow store with one lever.

---

## 4. `_inflow_to_moments` (`loss_representation.py:357-378`) — the lift detail

Already mapped in §2. The body (375-376) zero-fills `(*face, n)` and writes slot 0
from the scalar. The `n = self._n_face_moments` (370) = `per_axis**(ndim-1)`
(`loss_representation.py:311-321`). The `n == 1` early-return (371-372) is the
DD/Step identity. **#251 carry-a-moment-inflow rewrite**: when the trace already
carries the moment axis, the incoming `inflow` tuple is `(*face, 2^{d-1})`-shaped;
the method should detect that (rank check, or just return inflow when the trailing
axis matches `n`) and pass through. The transverse slopes then come from the
`face_view` read, not from a zero-fill. Single derivation site — change here, both
production consumers (1033/1120) + the oracle follow.

---

## 5. The prescribed-inflow PRODUCER (`boundary_source_sink.py:187-275`)

`BoundarySourceSink.prescribed_inflow(mesh, face_values)` builds the affine-BC
inhomogeneous term `q` (= `q.boundary`). Today:
- `bss = cls.zeros_on(mesh)` (257) — sized to `mesh.trace` (auto-widens).
- per face: `view = bss.face_view(face)` (265, "# (N, ng)"), shape-checks
  `arr.shape != view.shape` (267-272), then `view[inflow, :] = arr[inflow, :]`
  (274), inflow = `trace.inflow_indices_for_face(face)` (273).
- Shape today: `face_values[face]` is `(N, ng)` (the full per-face slot over
  ordinates); only inflow ROWS read.

**#251 widening**: once the slot is `(N, ng, ...transverse..., 2^{d-1})`, the
producer's `face_values` arrays grow the same axis, and:
- the `arr.shape != view.shape` check (267) auto-validates the wider shape (it
  compares against the slot shape, which grew);
- `view[inflow, :] = arr[inflow, :]` (274) — the `[inflow, :]` indexes ORDINATE
  axis 0 then takes everything trailing; for a `(N, ng, *transverse, 2^{d-1})`
  slot the `, :]` (a single trailing colon) is INSUFFICIENT — it must become
  `[inflow] = arr[inflow]` (or `[inflow, ...]`) to span ALL trailing axes. **This
  is a small but REAL edit** (the `, :` assumes exactly one trailing axis). Flag it
  for the method-implementer.
- the manufactured transverse moments are supplied BY THE CALLER (the MMS) in
  `face_values` — the producer just packs them. So the MMS (the L1 gate, Leg B's
  consumer per the test-architect memo §7) must project the manufactured inflow ψ
  onto the transverse face-Legendre basis and pass the moment-resolved
  `face_values`.

Docstring/`(N, ng)` references to update (not bugs, prose): boundary_source_sink.py
lines 76, 198, 233, 247, 265, 271.

---

## 6. The REFLECTIVE-COUPLING concern (`B` operator — the non-obvious Leg-B hazard)

The sibling `B` operator (`SNBoundaryOperator`, `orpheus/sn/boundary_operator.py`)
applies the reflective `R·G` coupling, reading/writing `face_view`. The core is
`_reflect_trace` (boundary_operator.py:160-228, esp. **219-227**):
```
face_in = boundary.face_view(face)         # 220
full = getattr(law, method)(face_in)       # 221  law.apply / law.apply_transpose
target = trace.inflow_indices_for_face(face) (or outflow for transpose)  # 222-226
out_boundary.face_view(face)[target] = full[target]   # 227
```
The realized `law` for reflective is `albedo * PermutationOperator(perm, axis=0)`
(`boundary_realizer.py:185-205`) — it permutes ORDINATES on axis 0 and is
**identity/broadcast over all trailing axes** (boundary_realizer.py:169-174,
190-204: "trailing axes (group, spatial / face) broadcast"). **So the angular
reflection passes a trailing `2^{d-1}` moment axis through unchanged for STORAGE.**

BUT — the deep physics caveat (NOT a storage issue, a CORRECTNESS one the
method-implementer + a numerics/qa review must settle): reflecting across a face
maps an outflow ordinate to an inflow ordinate. The TRANSVERSE (along-face) slope
moment is a spatial coefficient in the face's tangent plane. Under reflection
across that same face, the tangent-plane coordinate is PRESERVED (the reflection
flips only the normal-direction cosine), so the transverse slope should reflect
WITHOUT a sign flip — the permutation-on-axis-0 + moment-axis-passthrough is
PROBABLY correct as-is for reflective. Confirm against the operator-algebra adjoint
(`op.H`) and an MMS reflective-LD case. The `write[target] = full[target]` (227)
indexing is ordinate-axis-only → auto-widens for storage. **No `B` storage change;
flag the transverse-slope-under-reflection sign for verification.** (For periodic /
white the same passthrough reasoning applies; white AVERAGES over angle —
`AngularAverageOperator`, axis 0 — moment axis broadcasts.)

`reflect_into_inflow` (boundary_operator.py:268+, the direct-loop inflow seed) and
`_apply_faces` (230-264) both route through `_reflect_trace` (the single source of
truth, boundary_operator.py:188-190) — no separate widening.

---

## 7. L17 CROSSWALK — boundary trace (producer shape → consumer shape)

### INFLOW (read into the sweep)

| Boundary | Producer shape (today) | Consumer shape (today) | After #251 widening |
|---|---|---|---|
| `mesh.trace` slot (storage) | `(N, ng, *face_shape)` scalar (geometry.py:1047-1050) | flat `(total_size,)` | append `2^{d-1}` tail → `(N, ng, *face_shape, 2^{d-1})` (the §0 lever) |
| prescribed-inflow `q.boundary` producer | `face_values[face]` = `(N, ng)` (boundary_source_sink.py:233) | `view[inflow, :] = arr[inflow, :]` (274) | `(N, ng, *transverse, 2^{d-1})`; **fix the `, :` → full trailing span** |
| operator SOLVE seed | `rhs.boundary.face_view` whole-slot | `boundary_buf.face_view[:] =` (operator.py:1113-1114) | AUTO (whole-slot copy) |
| octant inflow read | `boundary.face_view(face)` scalar | `inflow` tuple (loss_representation.py:651-654 / 806) | AUTO (values carry moment axis) |
| scalar→moment lift | scalar `(*face,)` | `(*face, 2^{d-1})` slot-0=scalar, **transverse ZEROED** (loss_representation.py:375-376) | identity/pass-through (inflow already moment-valued) |
| full-cochain IN-edge seed | scalar | `(*face, 2^{d-1})` slot-0=scalar, **transverse ZEROED** (loss_representation.py:1214-1220) | seed ALL `2^{d-1}` from moment inflow |
| reflective `B` inflow emit | `face_view` whole-slot | `law.apply(face_in)`, ordinate-permute axis 0, moment passthrough (boundary_operator.py:220-227) | AUTO storage; **verify transverse-slope sign under reflection** |

### OUTFLOW (captured from the sweep, stored to the trace)

| Boundary | Producer shape (today) | Consumer shape (today) | After #251 widening |
|---|---|---|---|
| interior cochain capture | `(*face, 2^{d-1})` moment-valued | — | unchanged (already moment) |
| SOLVE capture→trace collapse | `(*face, 2^{d-1})` | `c[..., AVERAGE_MOMENT]` → scalar (loss_representation.py:1068-1069) | DROP the collapse → keep `2^{d-1}` |
| MATVEC capture→trace collapse | `(*face, 2^{d-1})` | `c[..., AVERAGE_MOMENT]` → scalar (loss_representation.py:1140-1141) | DROP the collapse |
| FullField SOLVE/MATVEC collapse (oracle) | `(*face, 2^{d-1})` | slot-0 scalar (loss_representation.py:1319-1320 / 1386-1387) | DROP the collapse |
| SOLVE shed→trace | scalar capture | `boundary_flux.face_view[oct_idx] =` (loss_representation.py:708) | needs moment-shaped slot (§0) to land `2^{d-1}` |
| MATVEC shed→`streamed` | scalar | `streamed[face][oct_idx] =` ; `streamed = zeros_like(face_view)` (loss_representation.py:777-780, 801) | AUTO (inherits widened slot) |
| MATVEC `B`-residual emit | scalar | `out_boundary.face_view[out_idx] = streamed[out_idx] - given[out_idx]` (loss_representation.py:819-823) | AUTO (ordinate-indexed, whole trailing) |

---

## 8. L20 DEPENDENCY AUDIT — construct/consume + backward-compat

### Who CONSTRUCTS a TraceSpace
EXACTLY ONE production site: `geometry.py:517` (`SNMesh.__init__`,
`TraceSpace.from_quadrature_and_layout`). The bare dataclass ctor is guarded against
(trace_space.py:394-399 raises if `omega_dot_n`/`layout` is None). Tests construct
via `mesh.trace` or the factory; none hand-build a moment trace.

### Who CONSUMES `mesh.trace` / `boundary.face_view` (production)
`mesh.trace` consumers (grep, prod): `transport/residuals/boundary_residual.py`,
`transport/source_sinks/boundary_source_sink.py`,
`transport/displacements/boundary_displacement.py`,
`transport/fields/boundary_flux.py`, `transport/fields/_bases.py`,
`derivations/continuous/mms/sn.py`, `sn/solver.py`, `sn/method_space.py`,
`sn/boundary_operator.py`, `sn/loss_representation.py`, `sn/sweep_schedule.py`.
`face_view` consumers (prod): `boundary_source_sink.py`, `fields/_bases.py`,
`fields/boundary_flux.py`, `solver.py`, `operator.py`, `boundary_operator.py`,
`loss_representation.py`.

**Backward-compat verdict (CONFIRMED): growing a `2^{d-1}` tail is byte-identical
for DD/Step.** The boundary FIELDS (BoundaryFlux/SourceSink/Residual/Displacement,
all `BoundaryField` subclasses) validate ONLY against `(layout.total_size,)`
(`_bases.py:516-525`) — never a hardcoded `(N, ng)` — so they accommodate ANY slot
shape the layout dictates. The `face_view` returns `slice_view` reshaped to the
slot shape (face_layout.py:121), so widening the slot widens the view transparently.
DD/Step → `face_moment_tail(1) == ()` → slot unchanged → every buffer / snapshot /
metric byte-identical. The metric (`_build_trace_metric_weights`,
trace_space.py:285-286) and `omega_dot_n` (per-ordinate) broadcast over the new
trailing axis with no change. The ONLY non-auto edits are the source-of-truth
sites (§9).

### Off-path (not affected)
- 2-D cylindrical `Mesh2D` has NO trace (`mesh.trace is None`) — guarded
  everywhere (e.g. boundary_source_sink.py:250-256); no SN sweep; irrelevant.
- 1-D curvilinear LD is REJECTED (`test_ld_curvilinear_scan_rejected`) — LD is
  2-D-Cartesian + 1-D-slab only. 1-D slab faces are `(N, ng)` (no transverse axis;
  `face_shape == ()` → `_n_face_moments = per_axis**0 = 1` → no moment tail even
  for LD-1D). **So Leg B's transverse face-slope is 2-D-ONLY by construction** (a
  1-D slab face is a point, no along-face direction). Confirm: in 1-D
  `_n_face_moments = per_axis**(1-1) = 1` always. The boundary moment widening is
  purely a 2-D (and higher) concern.

---

## 9. "WHAT CHANGES WHERE" for the method-implementer (minimal site list)

**The design choice (§0)**: persistent moment trace (layout widened) is the clean
path — it serves inflow read AND outflow store with one lever and keeps
`trace_space.py` untouched. The sites below assume that path.

**STORAGE (the lever) — 1 site:**
1. `geometry.py:boundary_face_layout` (1044-1050): append the scheme's face moment
   tail to each slot shape: `(N, ng, *face_shape(label), *face_moment_tail(per_axis**(ndim-1)))`.
   Read `per_axis = self.scheme.spatial_basis_per_axis` (available — scheme set at
   271-272 before trace at 517). DD/Step → `()` → byte-identical. **This is the
   single point that makes the trace moment-valued; everything below follows.**

**INFLOW lift — 1 site (BOTH paths share it):**
2. `loss_representation.py:_inflow_to_moments` (357-378, esp. 375-376): pass through
   a moment-valued inflow instead of zero-filling transverse slopes. (Generalize the
   `n == 1` early-return to "already moment-shaped → identity".) Drives both
   `_sweep_interior` (1033, SOLVE) and `_loss_action_interior` (1120, MATVEC) +
   the `_allocate_full_cochain` IN-edge seed (1214-1220, oracle).

**OUTFLOW store — 4 collapse sites to DROP (guarded `if n_face_moments > 1`):**
3. SOLVE prod: `loss_representation.py:1068-1069` (delete the collapse).
4. MATVEC prod: `loss_representation.py:1140-1141`.
5. SOLVE oracle: `loss_representation.py:1319-1320`.
6. MATVEC oracle: `loss_representation.py:1386-1387`.
   (The sheds at 708 SOLVE / 801 MATVEC + the `B`-residual emit 819-823 then land
   the `2^{d-1}` outflow into the now-moment-shaped slot AUTOMATICALLY.)

**PRODUCER — 1 real edit + docstrings:**
7. `boundary_source_sink.py:prescribed_inflow` line 274: `view[inflow, :] =
   arr[inflow, :]` → `view[inflow] = arr[inflow]` (span ALL trailing axes, not just
   one). Shape-check at 267 auto-validates the wider slot. Docstring `(N, ng)` →
   moment-shape at 76/198/233/247/265/271.

**VERIFY (no code change, flag for numerics/qa):**
8. `boundary_operator.py:_reflect_trace` (219-227): the reflective `PermutationOperator`
   passes the moment axis through for storage (boundary_realizer.py:169-204) — but the
   transverse-slope SIGN under reflection across a face must be verified correct
   (probably no flip; the tangent coordinate is preserved). MMS reflective-LD + `op.H`
   adjoint check.

### SOLVE vs MATVEC distinction (for the implementer)
- SOLVE path: operator seed (operator.py:1107-1115) → octant inflow read
  (loss_representation.py:712) → lift (1033) → walk → capture-collapse (1068-1069,
  DROP) → shed to `boundary_flux` (708).
- MATVEC path: octant inflow read (loss_representation.py:806) → lift (1120) → walk
  → capture-collapse (1140-1141, DROP) → shed to `streamed` (801) → `B`-residual
  emit (819-823). The MATVEC has NO external boundary source (the `B`-block residual
  is `streamed − given`); `Q_zero` (770) carries the moment tail for the bulk.
- They SHARE: `_inflow_to_moments` (357), `_n_face_moments` (311),
  `_spatial_moment_tail` (340), `face_moment_tail` (moment_layout.py:63), the
  `omega_dot_n`/metric (trace_space.py, untouched), and the widened layout (§0/site 1).

---

## Key file:line index

- `geometry.py:1044-1050` — `boundary_face_layout` slot shape (THE storage lever).
- `geometry.py:809` / `axis.py:349-360` — `face_shape(label)` (codim-1 cell counts).
- `geometry.py:271-272` vs `:517` — scheme set BEFORE trace (so per_axis reachable).
- `trace_space.py:341-388` — `from_quadrature_and_layout` (UNCHANGED; moment-ready
  metric+selectors broadcast over trailing axes, 285-286 / 407-423).
- `loss_representation.py:311-321` — `_n_face_moments = per_axis**(ndim-1)`.
- `loss_representation.py:357-378` (375-376) — `_inflow_to_moments` (the boundary lift).
- `loss_representation.py:1068-1069 / 1140-1141 / 1319-1320 / 1386-1387` — outflow
  capture-collapse-to-slot-0 (the 4 DROP sites).
- `loss_representation.py:708 / 801 / 819-823` — SOLVE/MATVEC sheds + `B`-residual emit.
- `loss_representation.py:1214-1220` — full-cochain IN-edge seed (oracle).
- `operator.py:1107-1115` — SOLVE operator boundary seed.
- `boundary_source_sink.py:187-275` (274) — `prescribed_inflow` producer (the `, :` edit).
- `boundary_operator.py:160-228` (219-227) — `_reflect_trace` (B; verify transverse sign).
- `boundary_realizer.py:169-205` — reflective = `PermutationOperator(perm, axis=0)`,
  trailing-axis passthrough.
- `moment_layout.py:60 / 63-73` — `AVERAGE_MOMENT`, `face_moment_tail` (`()` for DD).
- `face_layout.py:97 / 113-121 / 148` — `FaceSlot.shape`, `slice_view`, `total_size`.
- `_bases.py:516-525` — boundary-field validation (total_size ONLY → auto-widens).
- `linear_discontinuous.py:300` (per_axis=2) / `scheme.py:815` (per_axis=1) — neg control.

## Nexus staleness note
Graph built pre-#247-Leg-A; `solve_sn_fixed_source`/`_lift` line numbers in the
graph drift. All line numbers above are from Read/Grep of the working tree at HEAD
`d9396a2`. Trusted the file over the graph throughout.
