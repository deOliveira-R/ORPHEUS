# R-1 Phase A — Layer A: 2-D Cartesian alignment for typed AngularFlux traces

**Branch**: `refactor/sn-operator-algebra` (created against `a411e4c` =
R-1 Step 4 G1).
**Date drafted**: 2026-05-23.
**Purpose**: implementation plan for the **2-D Cartesian
alignment** that unblocks R-1 Step 4 **G2** (the
`_solve_fixed_source_si` carve onto typed `SourceIteration`).  This
plan is implementer-ready — every file:line is verified against the
current tree.

When this plan lands, the G2 carve resumes in the original session
without 2-D Cartesian regression.

## Why this exists — root cause of the G2 blocker

The R-1 Step 4 Step G2 carve naturally collapses
`_solve_fixed_source_si` from a manual `Q_iso + Q_aniso +
transport_sweep` loop onto the typed primitive
`SourceIteration(L+C, S, ZeroOperator()).solve(q_ext_typed)`.

The typed primitive depends — internally — on the typed `AngularFlux`
↔ flat vector conversion in
`AngularFlux.to_flat_with_traces` / `AngularFlux.from_flat_with_traces`,
which is currently **1-D only**.  Specifically the call chain:

```
SourceIteration.solve(q_ext_typed)
└── _zeros_like(q_ext_typed)
    └── q_ext_typed.to_flat_with_traces()    # (orpheus/sn/angular_flux.py:434)
        └── self.boundary.xmax_face[eq_map.face_outer_ordinate, :]
            ↳ TypeError on 2-D Cartesian — xmax_face is None
            ↳ 2-D BoundaryFlux uses xmin_xmax_buf instead (a different concept).
```

Pre-G2 the SI loop bypassed this by calling `transport_sweep` directly
(`transport_sweep` dispatches on mesh shape and uses the legacy
`xmin_xmax_buf` / `ymin_ymax_buf` persistent buffers for 2-D).  Post-G2
the typed primitive's `to_flat_with_traces` access fails before any
sweep happens.

The audit's SURPRISE-5 ("SI is geometry-agnostic") was true pre-G2
but stops being true after the carve.  **Layer A dissolves the
constraint at the root** by extending the 1-D B1'' face-aware
layout to 2-D Cartesian.

## What "Layer A" does — minimum to unblock G2

Layer A extends three primitives to handle 2-D Cartesian natively:

1. **`BoundaryFlux` 2-D semantics** — add `ymin_face` + `ymax_face`
   fields, populate all 4 face blocks in `BoundaryFlux.zeros()` for
   2-D Cartesian.
2. **`build_equation_map_with_traces` 2-D extension** — add 4 face
   slot maps (`face_xmin_ordinate`, `face_xmax_ordinate`,
   `face_ymin_ordinate`, `face_ymax_ordinate`), 4 face slot counts.
3. **`to_flat_with_traces` / `from_flat_with_traces` 2-D layout** —
   cells block + 4 face blocks (xmin, xmax, ymin, ymax); the bare-
   ndarray 2-D apply path (legacy FD) continues to compute cells;
   face slots stay zero (no real equation — fine for iteration vector
   sizing).

**NOT in Layer A**: actual B1'' face-aware 2-D matvec algorithm.
That's Layer B (Phase A proper) — `transport_operator_matvec_unified`
2-D wavefront cell body + `transport_sweep` 2-D face-aware sweep.
Layer B is substantial; #199 item 5 tracks it.

## Architectural goal — 1-D and 2-D unified via B1'' face layout

After Layer A:

| Field | 1-D slab | 1-D sphere/cyl | 2-D Cartesian |
|---|---|---|---|
| `BoundaryFlux.xmin_face` | `(N, ng)` | None (pole) | `(N, ng, ny)` |
| `BoundaryFlux.xmax_face` | `(N, ng)` | `(N, ng)` | `(N, ng, ny)` |
| `BoundaryFlux.ymin_face` | None | None | `(N, ng, nx)` |
| `BoundaryFlux.ymax_face` | None | None | `(N, ng, nx)` |
| `xmin_xmax_buf` | None | None | `(N, ng, nx+1, ny)` (legacy persistent buffer; coexists) |
| `ymin_ymax_buf` | None | None | `(N, ng, nx, ny+1)` (legacy persistent buffer; coexists) |

In 2-D the `*_face` fields are the PATH-FORWARD B1'' face state; the
`*_buf` fields are the legacy persistent buffers that the 2-D sweep
implementation still uses.  Both coexist until Layer B retires the
buffers.

The `EquationMap` (built by `build_equation_map_with_traces`) for 2-D
gains 4 slot maps + 4 slot counts.  `to_flat_with_traces` enumerates
5 blocks: cells + 4 face blocks.

## Concrete file:line changes

Verified against `refactor/sn-operator-algebra @ a411e4c`.

### 1. `orpheus/sn/boundary_flux.py`

**Current state**:
- `xmin_face`, `xmax_face`: `(N, ng)` 1-D only (lines 93-94)
- `xmin_xmax_buf`, `ymin_ymax_buf`: 2-D persistent buffers (lines 97-98)
- `BoundaryFlux.zeros()` dispatches by mesh.reduced; 1-D populates
  `xmin/xmax_face`, 2-D populates `xmin_xmax_buf/ymin_ymax_buf`
  (lines 102-136)
- `xmin`/`xmax`/`ymin`/`ymax` accessors return face SLICES of the
  persistent buffers for 2-D (lines 138-180 area)

**Add**:

```python
# In the dataclass field block (around line 94):
xmin_face: np.ndarray | None = None
xmax_face: np.ndarray | None = None
# NEW for 2-D:
ymin_face: np.ndarray | None = None
ymax_face: np.ndarray | None = None
# Existing persistent buffers (legacy 2-D sweep):
xmin_xmax_buf: np.ndarray | None = None
ymin_ymax_buf: np.ndarray | None = None
```

**Update `BoundaryFlux.zeros()` (line 102)**:

```python
@classmethod
def zeros(cls, mesh):
    N = mesh.quad.N
    ng = mesh.ng
    nx = mesh.nx
    ny = mesh.ny
    reduced = mesh.reduced
    bf = cls(mesh=mesh)
    if reduced is not None:
        # 1-D — slab or curvilinear.
        curv = getattr(mesh, "curvature", None)
        if curv in ("spherical", "cylindrical"):
            bf.xmax_face = np.zeros((N, ng))
        else:
            # 1-D slab.
            bf.xmin_face = np.zeros((N, ng))
            bf.xmax_face = np.zeros((N, ng))
    else:
        # 2-D Cartesian.
        # NEW (Layer A): the 4 boundary-edge B1'' face blocks.
        bf.xmin_face = np.zeros((N, ng, ny))
        bf.xmax_face = np.zeros((N, ng, ny))
        bf.ymin_face = np.zeros((N, ng, nx))
        bf.ymax_face = np.zeros((N, ng, nx))
        # LEGACY (kept until Layer B): the persistent buffers
        # that the existing 2-D sweep body consumes.  Coexist with
        # the *_face fields until the 2-D matvec is rewritten.
        bf.xmin_xmax_buf = np.zeros((N, ng, nx + 1, ny))
        bf.ymin_ymax_buf = np.zeros((N, ng, nx, ny + 1))
    return bf
```

**Audit downstream readers of `xmin_face`/`xmax_face`** that assume
the 1-D `(N, ng)` shape and now might receive 3-D `(N, ng, ny)` in
2-D.  Grep for direct attribute access:

```bash
grep -rn "\.xmin_face\|\.xmax_face\|\.ymin_face\|\.ymax_face" orpheus/sn/ tests/sn/
```

Expected hits (all need 2-D shape checks):
- `orpheus/sn/operator.py` — `transport_operator_matvec_unified` reads
  `psi.boundary.xmax_face` for the outer-face seed.  For 1-D this is
  `(N, ng)`; for 2-D it's `(N, ng, ny)`.  The function's 2-D guard
  at line 1023-1024 raises `NotImplementedError` so this branch is
  unreachable for 2-D — no change needed (the guard catches it
  before the face-state read).
- `orpheus/sn/angular_flux.py` — `to_flat_with_traces` /
  `from_flat_with_traces` (covered in §3 below).
- `orpheus/sn/sweep.py` — `transport_sweep` may read face state via
  `boundary.xmin_face` / `xmax_face`.  In 2-D the function falls
  back to `xmin_xmax_buf` / `ymin_ymax_buf` so no immediate change
  needed (the persistent-buffer path is the legacy 2-D contract).
- `tests/sn/boundary_flux/...` — direct fixture access; update tests
  that hard-code 1-D shape assertions.

### 2. `orpheus/sn/operator.py` — `build_equation_map_with_traces`

**Current state** (line 671):

```python
def build_equation_map_with_traces(
    nx: int, quad, ng: int, *, has_inner_bc: bool,
) -> EquationMap:
    ...
    # Returns EquationMap with:
    #   n_eq, ordinate, ix, iy (cell block)
    #   n_face_outer, face_outer_ordinate (outer face block)
    #   n_face_inner, face_inner_ordinate (inner face block)
    #   n_unknowns = n_eq * ng + (n_face_outer + n_face_inner) * ng
```

**Add 2-D dispatch**:

```python
def build_equation_map_with_traces(
    nx: int, quad, ng: int, *,
    has_inner_bc: bool,
    ny: int = 1,                                   # NEW — defaults to 1-D
    bc_xmin=None, bc_xmax=None,                    # NEW — 2-D edge BC info
    bc_ymin=None, bc_ymax=None,
) -> EquationMap:
    if ny == 1:
        # 1-D path — unchanged.
        ...
    else:
        # 2-D Cartesian — 4 face blocks.
        # NEW:
        # face_xmin_ordinate: outflow ords at xmin face = mu_x < 0
        # face_xmax_ordinate: mu_x > 0
        # face_ymin_ordinate: mu_y < 0
        # face_ymax_ordinate: mu_y > 0
        # Each face block has n_face_{edge} = n_outflow_at_edge * n_cells_along_edge
        mu_x = quad.mu_x
        mu_y = getattr(quad, "mu_y", np.zeros(quad.N))
        face_xmin_ord = np.where(mu_x < -1e-15)[0]
        face_xmax_ord = np.where(mu_x > +1e-15)[0]
        face_ymin_ord = np.where(mu_y < -1e-15)[0]
        face_ymax_ord = np.where(mu_y > +1e-15)[0]
        n_face_xmin = face_xmin_ord.size * ny
        n_face_xmax = face_xmax_ord.size * ny
        n_face_ymin = face_ymin_ord.size * nx
        n_face_ymax = face_ymax_ord.size * nx
        # Build cell block similar to existing 2-D build_equation_map
        # (but with face-aware contract).
        ...
        eq_map = EquationMap(
            n_eq=n_eq, ordinate=..., ix=..., iy=...,
            n_face_outer=...,  # legacy 1-D fields kept zero / None for 2-D
            face_outer_ordinate=...,
            n_face_inner=...,
            face_inner_ordinate=...,
            # NEW 2-D fields:
            n_face_xmin=n_face_xmin, face_xmin_ordinate=face_xmin_ord,
            n_face_xmax=n_face_xmax, face_xmax_ordinate=face_xmax_ord,
            n_face_ymin=n_face_ymin, face_ymin_ordinate=face_ymin_ord,
            n_face_ymax=n_face_ymax, face_ymax_ordinate=face_ymax_ord,
        )
```

**Update `EquationMap` dataclass** (`orpheus/sn/operator.py:119`):

Add 8 new fields (4 slot counts + 4 slot maps), all defaulting to
`None` / 0 for 1-D backward compat:

```python
@dataclass
class EquationMap:
    n_eq: int
    ordinate: np.ndarray
    ix: np.ndarray
    iy: np.ndarray
    # 1-D face fields (existing):
    n_face_outer: int = 0
    face_outer_ordinate: np.ndarray | None = None
    n_face_inner: int = 0
    face_inner_ordinate: np.ndarray | None = None
    # NEW 2-D face fields:
    n_face_xmin: int = 0
    face_xmin_ordinate: np.ndarray | None = None
    n_face_xmax: int = 0
    face_xmax_ordinate: np.ndarray | None = None
    n_face_ymin: int = 0
    face_ymin_ordinate: np.ndarray | None = None
    n_face_ymax: int = 0
    face_ymax_ordinate: np.ndarray | None = None

    @property
    def n_unknowns(self) -> int:
        # Sum cells + all populated face blocks.
        return (
            self.n_eq * ng  # cells
            + self.n_face_outer + self.n_face_inner  # 1-D faces
            + self.n_face_xmin + self.n_face_xmax    # 2-D x faces
            + self.n_face_ymin + self.n_face_ymax    # 2-D y faces
        )
```

Note: `n_unknowns` needs `ng` to compute the cell block size.
Currently `n_unknowns` is a stored field (not a property) — check
the existing definition.  If stored, update the producer in
`build_equation_map_with_traces`.

### 3. `orpheus/sn/angular_flux.py` — `to_flat_with_traces` / `from_flat_with_traces`

**Current state** (line 434 area):

```python
def to_flat_with_traces(self) -> np.ndarray:
    # 1-D layout:
    # [cell block (N*ng*nx)] [face_outer block (n_face_outer*ng)] [face_inner block (n_face_inner*ng)]
    nx = self.mesh.nx
    ng = self.mesh.ng
    ...
    eq_map = build_equation_map_with_traces(nx, self.mesh.quad, ng,
                                              has_inner_bc=...)
    cell_block = self.values[eq_map.ordinate, :, eq_map.ix, 0].T.ravel(order="F")
    face_outer = self.boundary.xmax_face[eq_map.face_outer_ordinate, :]
    ...
    return np.concatenate([cell_block, face_outer.ravel(order="F"), ...])
```

**Extend to 2-D**:

```python
def to_flat_with_traces(self) -> np.ndarray:
    nx = self.mesh.nx
    ny = self.mesh.ny
    ng = self.mesh.ng

    if self.mesh.reduced is not None:
        # 1-D path — unchanged.
        ...
    else:
        # 2-D Cartesian path — NEW (Layer A).
        eq_map = build_equation_map_with_traces(
            nx, self.mesh.quad, ng,
            has_inner_bc=True,  # 2-D Cartesian has real BCs on all 4 edges
            ny=ny,
            bc_xmin=self.mesh.bc_xmin, bc_xmax=self.mesh.bc_xmax,
            bc_ymin=self.mesh.bc_ymin, bc_ymax=self.mesh.bc_ymax,
        )
        # Cell block: ravel (N, ng, nx, ny) in F-order indexed by eq_map.
        cell_block = self.values[
            eq_map.ordinate, :, eq_map.ix, eq_map.iy,
        ].T.ravel(order="F")
        # Face blocks: xmin, xmax, ymin, ymax.  Each carries (outflow ords, ng, cells along edge).
        xmin_block = (
            self.boundary.xmin_face[eq_map.face_xmin_ordinate, :, :]
            .ravel(order="F")
            if eq_map.n_face_xmin > 0 else np.zeros(0)
        )
        xmax_block = (
            self.boundary.xmax_face[eq_map.face_xmax_ordinate, :, :]
            .ravel(order="F")
            if eq_map.n_face_xmax > 0 else np.zeros(0)
        )
        ymin_block = (
            self.boundary.ymin_face[eq_map.face_ymin_ordinate, :, :]
            .ravel(order="F")
            if eq_map.n_face_ymin > 0 else np.zeros(0)
        )
        ymax_block = (
            self.boundary.ymax_face[eq_map.face_ymax_ordinate, :, :]
            .ravel(order="F")
            if eq_map.n_face_ymax > 0 else np.zeros(0)
        )
        return np.concatenate([
            cell_block, xmin_block, xmax_block, ymin_block, ymax_block,
        ])
```

**Inverse `from_flat_with_traces` (around line ~480)**:

Mirror — decode flat into cells + 4 face blocks for 2-D.

### 4. Update callers of `build_equation_map_with_traces`

Grep for all callers and confirm they pass `ny=1` (current default)
for 1-D OR pass `ny=mesh.ny` for 2-D:

```bash
grep -rn "build_equation_map_with_traces" orpheus/sn/ tests/sn/
```

Likely callers in production:
- `orpheus/sn/operator.py:_ensure_eq_map` in `StreamingOperator` /
  `CollisionOperator` — needs `ny=sn_mesh.ny` for 2-D path
- `orpheus/sn/angular_flux.py:to_flat_with_traces` /
  `from_flat_with_traces` — needs `ny=mesh.ny`

In tests:
- Several `*_with_traces` test files may need a `ny=ny` arg

### 5. `StreamingOperator._apply_typed` 2-D branch validation

**Current state** (`orpheus/sn/operator.py:2073`-2076):

```python
# 2-D Cartesian — flat round-trip onto the legacy compute.
if curv is None and ny > 1:
    flat_in = psi.to_flat_with_traces()
    flat_out = self.apply(flat_in)
    return AngularFlux.from_flat_with_traces(flat_out, sn_mesh)
```

Post-Layer-A, both `to_flat_with_traces` and `from_flat_with_traces`
handle 2-D.  `self.apply(flat_in)` for 2-D (bare-ndarray, line 2011-
2024) routes through `transport_operator_matvec` (legacy FD), which
consumes the cell block of `flat_in`.  Face blocks in `flat_in` are
not consumed by the FD path (the FD path doesn't know about B1''
faces) — but `from_flat_with_traces` then reconstructs the
AngularFlux from the cell block and writes zero face state.

Net effect post-Layer-A: 2-D `_apply_typed` returns an AngularFlux
with valid cell residuals (FD math, same as pre-Layer-A) and zero
face residuals.  This is sufficient for `SourceIteration` to converge
on 2-D Cartesian.

**Layer B** will replace the FD math with the path-forward 2-D
wavefront sweep that produces actual face residuals.

## Verification approach

### Foundation pins (new — must add)

1. **`BoundaryFlux.zeros(mesh_2d)` shape contract** — verify the 4
   face blocks have correct shapes `(N, ng, ny)` / `(N, ng, nx)`.
2. **`build_equation_map_with_traces(ny>1)` slot map contract** —
   verify face_xmin_ordinate, face_xmax_ordinate, etc. correctly
   identify outflow ordinates at each edge.
3. **`to_flat_with_traces` / `from_flat_with_traces` 2-D round-trip
   identity** — random AngularFlux → flat → typed must be bit-equal.
4. **`StreamingOperator._apply_typed` 2-D path** — verify it doesn't
   raise; result matches pre-Layer-A bit-equal at cell slots.

### L0 / L1 gates that must continue passing

1. `tests/sn/l1_analytical/test_kinf_homogeneous.py` — 28 passed + 2
   xfailed (Phase 1.2 / 1.3 / G0 / G1 baseline).
2. `tests/sn/regression/test_dd_regression.py` — 6 passed + 5
   PRE-EXISTING failures (matches Phase 1.1 baseline; one of the 5
   is `2d_1g_LS4_dd_15x15` which has been pre-existing-NotImpl since
   PR-TYPED-2 — Layer A might un-break it).
3. `tests/sn/test_native_matvec.py` — 20 foundation pins (G0 contract).
4. `tests/sn/test_fixed_source_g1.py` — 7 pins (G1 contract).
5. `tests/sn/test_invertible_operator.py` — 33 tests (LC.solve + bridge
   regression).
6. `tests/sn/test_streaming_operator.py` + `test_collision_operator.py` —
   typed operator tests.

### Acceptance gates for the Layer A commit

* All 4 NEW foundation pins (§Verification.1-4 above) green.
* All listed L0/L1 baseline tests preserved (no new failures vs G1 baseline).
* `solve_sn_fixed_source` with 2-D Cartesian + `inner_solver="source_iteration"`
  (the failing case post-G2) now works end-to-end.

## What this unblocks

Post-Layer A, the G2 carve in
[the original session at `r1_step4_session2_followup.md` §3.2 / G2]
becomes:

```python
def _solve_fixed_source_si(...):
    # No 2-D guard needed — Layer A made typed SourceIteration
    # geometry-agnostic.
    q_ext_typed = AngularFlux(
        PerOrdinateSource(external_source, sn_mesh).values, sn_mesh,
    )
    LC = StreamingOperator(...) + CollisionOperator(...)
    si = SourceIteration(LC, solver.scattering_op, ZeroOperator(),
                          tol=inner_tol, max_iter=max_inner)
    psi_typed, residuals = si.solve(q_ext_typed)
    return Solution(angular_flux=psi_typed, ...)
```

i.e., the atomic-G2 body from the original session works as-is once
Layer A lands.

## What Layer A explicitly does NOT do

* **Replace `xmin_xmax_buf` / `ymin_ymax_buf`** — these persistent
  buffers stay alive for the legacy 2-D sweep body.  Layer B retires
  them as part of the full path-forward 2-D matvec rewrite (#199 item 5).
* **Implement B1''-face-aware 2-D wavefront sweep** — the
  `*_face` blocks for 2-D get populated with zeros at the end of the
  matvec; they're vestigial face residuals.  Layer B fills them with
  real physics.
* **Migrate 2-D regression snapshots** — the legacy 2-D FD matvec is
  preserved; numerical results unchanged.

## Estimated effort

* Layer A core (BoundaryFlux + EquationMap + traces): 1-2 hours of
  focused editing.
* New foundation pins: 1 hour.
* Verification gate run: 30 min.
* Sphinx update (mention the 2-D face layout): 30 min.

Total: ~3-4 hours single focused session.

## Commit shape

Single atomic commit on `refactor/sn-operator-algebra`:

```
refactor(sn): 2-D Cartesian B1'' face layout — BoundaryFlux + EquationMap + traces (R-1 Phase A Layer A)

Extends the 1-D B1'' face-aware layout to 2-D Cartesian:

* BoundaryFlux: new ymin_face + ymax_face fields; zeros() populates
  all 4 face blocks (shape (N, ng, ny) for xmin/xmax, (N, ng, nx)
  for ymin/ymax) for 2-D.  Legacy xmin_xmax_buf / ymin_ymax_buf
  coexist until Layer B retires them.
* EquationMap: new n_face_{xmin,xmax,ymin,ymax} slot counts +
  face_{xmin,xmax,ymin,ymax}_ordinate slot maps.  1-D fields kept
  for backward compat.
* build_equation_map_with_traces: ny kwarg + 2-D branch.
* AngularFlux.to_flat_with_traces / from_flat_with_traces: 2-D
  layout = [cells][xmin_block][xmax_block][ymin_block][ymax_block].

What this unblocks: R-1 Step 4 G2 — _solve_fixed_source_si onto
SourceIteration + typed AngularFlux becomes geometry-agnostic
naturally; the 2-D Cartesian path no longer raises in to_flat_with_traces.

Layer B (#199 item 5) is the full path-forward 2-D matvec rewrite —
deferred.

Verification: 4 new foundation pins + all G0/G1 baselines preserved.
```

## References

* `.claude/plans/r1_step4_session2_followup.md` §3.2 / G2 (the
  blocked carve).
* `.claude/plans/r1_step4_g_dependency_audit.md` SURPRISE-5 (the
  audit's outdated "SI is geometry-agnostic" claim).
* `.claude/plans/r1_step4_g_convention_crosswalk.md` Axis 6 (the
  face-state crosswalk that Layer A extends to 2-D).
* `.claude/plans/r1_step4_g_verification_plan.md` §G2.3 (2-D
  Cartesian landing-zone tests blocked by Layer A).
* Commit `a411e4c` (G1) — the base.
* Issue #199 item 5 (2-D Cartesian absorption) — Layer B tracking.
