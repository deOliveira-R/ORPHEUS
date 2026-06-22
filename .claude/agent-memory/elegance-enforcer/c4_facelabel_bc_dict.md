---
name: c4-facelabel-bc-dict
description: C4 carve (#220) — SN boundary inventory dimension-generic; FaceLabel.face_name crosswalk + face_labels-derived bc dict; PASS clean. The N-D campaign's boundary tail.
metadata:
  type: project
---

# C4 — SN boundary inventory dimension-generic (#220, the N-D-campaign BOUNDARY tail)

PASS clean (commit `515f12f` C4.1 + working-tree C4.2, on `worktree-sn-nd-layout`).
Sibling to [[c3_6_honest_d_dispatch]] (the SWEEP tail) and [[nd_generic_carve_tuple_axis]].
The boundary-inventory analog of the d-generic sweep: kill per-geometry hand-lists,
single-source the axis↔face-name crosswalk, realize pole-as-structural-absence.

**HEADLINE = a latent d=3 CORRECTNESS bug killed by construction.** The pre-C4
`_resolve_one` hand-listed `axis = "y" if face in ("ymin","ymax") else "x"` — which
would have silently built the WRONG reflection permutation for a z-face (mapped every
non-y face to "x"). C4 re-keys on `AXIS_NAMES[label.axis_index]`, correct at any d.
This is the C4 equivalent of C3.6's ERR-056 z-face existence hole — a d=3 EXISTENCE/
CORRECTNESS bug that the dimension-generic rewrite dissolves, pinned by a NEW test
(`test_2d_reflective_y_face_builds_y_axis_permutation`, verifies per-axis partner).

## The SSOT chain (genuinely single-sourced — verified)
`face_labels(axes)` → per-label `FaceLabel.face_name` (axis.py:118) →
BOTH `boundary_face_layout` keys AND `SNMesh.bc` dict keys. BC declaration access
`axes[label.axis_index].bc[label.endpoint]` keys on the SAME endpoint string
`face_labels` iterates (`axis.endpoints`), and `AxisMesh.bc`/`RadialAxisMesh.bc` are
keyed on those labels — coextensive BY CONSTRUCTION. The two pre-C4 "outer"→"max"
translation sites (geometry.py:450 + :948) are GONE; `_ENDPOINT_SUFFIX`
(`{"min":"min","max":"max","outer":"max"}`, axis.py:85) is the sole suffix crosswalk.

**NON-DUPLICATE to not mistake for a twin:** `_OUTWARD_ENDPOINTS = frozenset({"max",
"outer"})` (axis.py:401) co-locates "outer" with "max" same as `_ENDPOINT_SUFFIX`, BUT
it is a DIFFERENT crosswalk (endpoint→outflow-SIGN for `face_outflow_ordinates`, a
different semantic layer = which ordinates leave the domain) — encodes the same physical
fact ("solid-radial outer surface behaves like a max face") for two distinct purposes.
NOT a Pattern-2 violation. When auditing C5, do not "unify" these.

## Pattern 4 — pole-as-structural-absence FULLY realized
Pre-C4 curvilinear pole = `bc_left=None`/`bc_xmin=None` sentinel. Post-C4 = NO dict
entry. Grep-verified ZERO `bc[...] is None` checks and ZERO `bc.get(...)` calls anywhere.
Live pole detection = `"xmin" in boundary.layout.faces` (operator.py:403/637/862) —
membership on the trace inventory (face_labels-derived), not a None test. The meaningless
pole-BC is now unspellable.

## The fail-loud ValueError placement (principled, NOT a dodge)
`FaceLabel.face_name` raises ValueError on a non-canonical endpoint (overridden
`AxisMesh.label_low="left"`). RULING: principled intermediate boundary. `AxisMesh`
DELIBERATELY permits label overrides (axis.py:209); `face_name` is where a renamed label
first collides with the `"{axis}{min|max}"` world. VERIFIED the override is a GUARDED-
FUTURE concern not a live bug: the only `label_low="left"`/`label_outer="surface"`
constructions are in `test_axis_primitive.py`; production `axes_from_legacy_mesh`
(axis.py:472-540) never overrides. Making the illegal label unrepresentable at the
`AxisMesh` CONSTRUCTOR is C5's (axis-native constructor) territory — destination already
named in the carve's own docstrings. DISCRIMINATOR for C5 review: if C5 lands the native
constructor, the `face_name` ValueError should become unreachable-by-construction (or
move to construction).

## Dict ergonomics ruling: plain dict + KeyError is CORRECT (no MappingProxyType)
Sibling post-init `SNMesh` attrs (`_streaming_axes`, `_trace`) are PLAIN mutable attrs,
no MappingProxyType on the class. `SNMesh` is construct-once immutable-by-discipline.
Wrapping `bc` alone = inconsistent ceremony (anti-pattern #10). Plain-dict KeyError-on-
miss (no masking `.get` default) is the right fail-loud, pinned at
`test_bc_dict_misses_and_retired_attributes_fail_loud`. GENERAL RULE for this codebase:
do NOT recommend MappingProxyType for a single post-init mesh attr when its siblings are
plain attrs — match the convention, don't invent ceremony.

## Retirement = COMPLETE
Zero surviving `SNMesh.bc_xmin`/`bc_left`/etc in prod OR tests (residual grep hits =
`AxisMesh.bc` the producer + unrelated `bc_kind` in mc/moc). `mesh.bc_left/bc_right` reads
in solver.py:72 / mc / cp / moc / axis.py converter are on the INPUT Mesh1D/Mesh2D (the
declaration), correctly untouched — DISCRIMINATOR: `mesh.bc_*` = input mesh OK; `sn_mesh.
bc_*` / `self.bc_*` on SNMesh = retired. Curvilinear bypass + 1D/2D isinstance split +
y-placeholder block (with `SNMethodSpace.minimal`) deleted outright, no shims. All stale
docstrings swept (operator/vacuum/_bound_compat/boundary_operator/pole_angular_closure).

## Test migration = EXEMPLARY (retirement = test migration done right)
- Crosswalk pins (test_face_name_crosswalk.py): hand-transcribed `(axis,endpoint)→name`
  table with EXPLICIT "must NOT be generated from AXIS_NAMES/_ENDPOINT_SUFFIX" note —
  genuine mirror-not-import.
- `test_bc_inventory_equals_face_layout_across_geometries`: TWO-CLAUSE — (1) `set(sn.bc)`
  vs HAND-transcribed face set per geometry (catches bc+layout drifting TOGETHER), (2)
  `set(sn.bc)==set(boundary_face_layout.faces)` (catches drift APART). The hand-set is the
  non-tautological anchor. This two-clause shape is the right way to pin an SSOT-derived
  pair (a pure cross-key check alone is blind to co-movement).
- y-placeholder TP pin MIGRATION (test_bound_compat.py:238): HISTORY note + the migrated
  pin STRENGTHENS (old = "some PermutationOperator"; new = `assert_array_equal(perm,
  reflection_index("x"))`). Moved onto the slab's REAL xmin/xmax faces.

## Gates
`tests/sn/operators tests/sn/primitives` 770 passed / 4 skip / 1 xfail (incl new pins).
Migrated curvilinear/sweep/bound_compat suites 46 passed / 31 xfailed exit 0 (31 = pre-
existing #212 deselect, not introduced). Curvilinear matvec suites are SLOW (~12 min).
