---
name: c5-1-axis-primary-snmesh
description: C5.1 (#225) axis-primary SNMesh inversion review — two-surface-funnel-into-_init_core PASS-with-nits; mat_map twin-validation + 3-site AxisCoord↔CoordSystem map nits
metadata:
  type: project
---

C5.1 (#225, commit `2dd16c5`, branch worktree-sn-nd-layout) — the axis-primary
inversion keystone before 3-D admission. Reviewed PASS-with-nits. The build:
`SNMesh` (plain class, NOT a dataclass) gets TWO construction surfaces funneling
into ONE body `_init_core` (kw-only): `__init__(mesh,...)` = legacy inbound
adapter (`axes_from_legacy_mesh` ONCE + extract `mesh.mat_ids`/`mat_map`);
`from_axes` = native surface via `cls.__new__(cls)` + `_init_core` (stores axes
VERBATIM, no round-trip). `__new__` is the CORRECT tool here precisely because
SNMesh is a plain class — no dataclass `__init__`/frozen/`__post_init__` to bypass.

**Why this is the right shape (reinforce):** the inversion KILLS a real bug
habitat — the pre-C5.1 axes→mesh→axes round-trip silently reset custom endpoint
labels to min/max/outer. Now custom labels FAIL LOUD (FaceLabel.face_name
ValueError) — Pattern-4 illegal-state surfaced at construction. SSOT for material
assignment is `self.mat_map` set in `_init_core`; the legacy `mesh.mat_ids`/
`.mat_map` field is INBOUND PROVENANCE only (`_validate_materials` reads
`self.mat_map`, geometry.py:595). The `assert isinstance(mesh, Mesh1D)` lines in
the curvilinear match arms are -O-stripped type-narrowing, sound (the 1-axis
curvilinear adapter is always Mesh1D by construction).

**NIT 1 (Pattern-2 twin-validation) — mat_map shape validated TWICE on the
`from_axes` path.** `legacy_mesh_from_axes` (axis.py:623) validates `mat_map.shape
!= spatial_shape(axes)` AND `_init_core` (geometry.py:346) validates the SAME
mat_map against the SAME shape. Pre-C5.1 the validation lived ONLY in
`legacy_mesh_from_axes` (parent geometry.py had NO mat_map validation — just
`self.mat_map = mesh.mat_ids`); C5.1 ADDED the `_init_core` copy. Bug habitat:
the two error messages differ ("legacy_mesh_from_axes: ..." vs "SNMesh: ...");
the day the adapter retires (C5.2–C5.5) the axis.py copy goes with it, and IF the
`_init_core` copy were ever the one removed the `from_axes` shape guard would
vanish. Low-severity NOW because both validate identically and `_init_core` is the
keeper (it's on the SSOT path). Acceptable-for-now; the adapter-side copy
dissolves with the adapter.

**NIT 2 (Pattern-7) — AxisCoord↔CoordSystem mapping now in 3 places.**
(a) `coord_system` axis.py:479 NEW pure primitive (fwd AxisCoord→CoordSystem);
(b) `axes_from_legacy_mesh` axis.py:532-561 (rev CoordSystem→AxisCoord);
(c) `legacy_mesh_from_axes` axis.py:631-654 (fwd, embedded in adapter ctor).
(a) and (c) encode the SAME fwd map. (c) can't just call (a) — it needs per-branch
`assert isinstance` + per-branch BC-field wiring (bc_low/bc_high vs bc_outer). On
the documented retirement path. Defer-acceptable; flag if a 4th site appears.

**NIT 3 (comment honesty, minor over-claim).** `test_d2_..._axis_vs_legacy`
comment says nx≠ny + asymmetric BCs make an axis-swap (Mode-2) detectable — TRUE,
but via spatial_shape/(nx,ny)/widths/streaming/BC, NOT via mat_map (the test's
`np.arange(4*7)%1` mat_map is all-zeros, swap-blind). The asymmetry coverage is
real; the mat_map just doesn't contribute to it.

Tests EXEMPLARY otherwise: np.testing.assert_* only (Mode-8 -O-safe), verbatim
`is`-identity pin, custom-label fail-loud, curvilinear adapter path
(sphere+cylinder parametrized), d≥3 gate names C5.5, coord_system primitive
direct. d≤2 sha256 affine goldens byte-identical (313+925 green).
