---
name: issue-247-moment-source-consumption-path
description: "#247 file:line map of the CURRENT external-source + boundary-trace consumption path — where the slope/transverse moments are zeroed, where a moment-valued Q̂ must land to join the EXISTING Σ_s·φ̂ slope-source matvec, the L17 producer→consumer crosswalk, and the L20 caller audit proving a typed-union widening is backward-compatible."
metadata:
  type: project
---

# #247 — moment-resolved external source + boundary trace: consumption-path map

**Goal of #247** (from `orpheus/derivations/continuous/mms/sn.py:1142-1159`): the
production `solve_sn_fixed_source` consumes a **flat-in-moment** external source
(slope rows `Q̂` zeroed) and a **scalar** boundary trace (transverse face-slope
zeroed). The manufactured slope-source `Q̂` and the transverse-face-moment inflow
ARE derived in the MMS (Branch 1 is slope-source-READY); only the production
CONSUMPTION path is deferred. #247 = wire the consumption so a moment-resolved
`Q̂` / boundary trace LANDS in the SAME matvec/sweep path the scattering slope
source `Σ_s·φ̂` (the "S3" moment matvec) already travels — it joins an EXISTING
path, NOT a new one.

**Tree read**: main checkout, branch `refactor/sn-foundation-cleanup`, HEAD
`dd4f542`. NOT a worktree. Nexus graph built at `f6475e5` (stale by the #241/#238
commits) — all line numbers below are from Read/Grep of the working tree (the
graph's `solve_sn_fixed_source@1932` / `_lift@1896` match the file; the brief's
~1876/~1884 estimates were the drift).

---

## 1. External-source entry — `orpheus/sn/solver.py`

### 1a. `solve_sn_fixed_source` (def @ **1932**; signature 1932-1945)

```python
def solve_sn_fixed_source(
    materials, mesh, quadrature,
    external_source: "np.ndarray | TimedFullField",   # <- the typed union ALREADY exists
    boundary_condition="vacuum", scattering_order=0,
    max_inner=1000, inner_tol=1e-12, inner_solver=None,
    inner_schedule="gauss_seidel", mat_map=None, scheme=None,
) -> Solution:
```

- `external_source` is ALREADY a typed union: `np.ndarray (N,ng,*spatial)`
  bulk-only/vacuum, OR a `TimedFullField` composite (bulk ⊕ prescribed-inflow
  boundary). Docstring 1964-1988.
- The only normalisation site is `_build_fixed_source_rhs(external_source, sn_mesh)`
  at **2084** (Cardinal Rule 2 — "one construction point; shape validation lives
  inside the helper"). Both SI (`_solve_fixed_source_si`@2087) and Krylov
  (`_solve_fixed_source_krylov`@2096) paths consume what it returns.

### 1b. `_build_fixed_source_rhs` (def @ **1810**) — the validation + lift

- **The flat-shape validation** (lines **1875-1880**), against
  `expected = (N, ng, *sn_mesh.spatial_shape)` (built @ 1847):
  ```python
  # Issue #196 PR-INDEX-5: bulk source principled (N, ng, *spatial).
  if bulk_values.shape != expected:
      raise ValueError(
          f"fixed-source bulk shape {bulk_values.shape} does not match "
          f"(N, ng, *spatial) = {expected}")
  ```
  **This is the gate #247 must widen.** A moment-resolved
  `(N,ng,*spatial,2^d)` source FAILS here today (rank too high → `!= expected`).
  For `TimedFullField` input the bulk is `external_source.bulk.values` (1850); the
  boundary is validated separately by `total_size` (1851-1860) and re-homed via
  `BoundarySourceSink.from_mesh(..., sn_mesh)` (1860). A bare non-array object →
  `TypeError` (1862-1872).
- The **moment-lift call + the "Q̂≠0 ansatz is S4" docstring note** (lines
  **1881-1887**):
  ```python
  # The external source carries the trailing 2^d spatial-moment axis at a
  # multi-moment closure (#240 D5b-S3) so it composes with the moment-carrying
  # scattering source ``S.apply(ψ)`` in the SI rhs ``q_ext + S.apply(ψ)``.  The
  # EXTERNAL source is flat-in-moment (Q̂ = 0 — the slope rows are zero; the
  # strengthened Q̂≠0 ansatz is S4): lift onto slot 0 (average), rest zero.
  # DD/Step (per_axis == 1) → no lift, byte-identical.
  bulk_values, per_axis = _lift_external_source_to_moments(bulk_values, sn_mesh)
  ```
- Returns the composite (1888-1893): `TimedFullField(bulk=AngularSourceSink.from_mesh(
  bulk_values, sn_mesh, spatial_moments=per_axis), boundary=boundary)`. The
  `spatial_moments=per_axis` sets the bulk's `SpatialMomentSpace` factor — this is
  what makes `q_ext.space == S.apply(ψ).space` so the typed `+` is legal (see §3).

### 1c. `_lift_external_source_to_moments` (def @ **1896**) — WHERE the slope rows are zeroed

```python
def _lift_external_source_to_moments(bulk_values, sn_mesh) -> tuple[np.ndarray, int]:
    from orpheus.sn.spatial._ubld import AVERAGE_MOMENT, face_moment_tail
    per_axis = sn_mesh.scheme.spatial_basis_per_axis                  # 1907
    tail = face_moment_tail(per_axis ** sn_mesh.ndim)                 # 1908 — cell tail per_axis^d
    if tail == ():                                                    # 1909 DD/Step → byte-identical
        return bulk_values, per_axis
    lifted = np.zeros((*bulk_values.shape, *tail), dtype=bulk_values.dtype)  # 1911 ALL-ZERO buffer
    lifted[..., AVERAGE_MOMENT] = bulk_values                         # 1912 average ← flat; SLOPES STAY 0
    return lifted, per_axis                                           # 1913
```

**The slope-zeroing IS line 1911-1912**: a fully-zero `(*flat, 2^d)` buffer, with
only slot `AVERAGE_MOMENT (=0)` written from the flat input — every other moment
(the `Q̂` slope rows) is left zero. This is the single point #247 changes: the
moment-resolved branch must instead copy ALL `2^d` rows from a moment-valued input
(or evaluate a callable producing them) rather than zero-fill + slot-0.

`AVERAGE_MOMENT` / `face_moment_tail` are re-exported from `orpheus.sn.spatial._ubld`
but the canonical home is `orpheus/numerics/moment_layout.py` (#245 layering fix):
- `AVERAGE_MOMENT = 0` (`moment_layout.py:55`) — slot-0 cell/face average index.
- `face_moment_tail(n)` (`moment_layout.py:58-68`): `() if n==1 else (n,)` — "append
  a trailing moment axis iff >1 moment" (DD/Step → `()` keeps buffers byte-identical).

**Callers of `_lift_external_source_to_moments`** (Nexus + grep): EXACTLY ONE
production caller — `_build_fixed_source_rhs` (`solver.py:1887`). (The eigenvalue
path does NOT call it — see §5. The docstring's "shared by the fixed-source and
eigenvalue paths" refers to the moment-lift CONCEPT being single-sourced for any
future eigenvalue external-source hook, not a current call.)

---

## 2. Boundary path — `_inflow_to_moments` + the `trace` producer

### 2a. The `trace` producer — `mesh.trace` is a `TraceSpace`, NOT in geometry.py per se

`SNMesh.trace` (`orpheus/sn/geometry.py:938` property) returns `self._trace`, built
once in the SNMesh constructor at **`geometry.py:516-519`**:
```python
from orpheus.numerics.spaces.trace_space import TraceSpace
self._trace = TraceSpace.from_quadrature_and_layout(
    self.quad, self.boundary_face_layout,
)
```
The trace is **geometry-blind** (quadrature + face names only; comment 515). It
carries **one SCALAR value per face per ordinate per group** — there is NO
`2^{d-1}` transverse-face-moment axis on the trace today. The boundary inflow the
sweep reads (`boundary.face_view(face)` — e.g. `loss_representation.py:806`,
`operator.py:1108-1115`) is therefore scalar-per-face.

### 2b. `_inflow_to_moments` (def @ `orpheus/sn/loss_representation.py` **357-378**) — WHERE the transverse face-slope is zeroed

On `_LossRepresentation` (base). Width source = `_n_face_moments` (property @ 311-321):
`per_axis ** (ndim - 1)` = `2^{d-1}` per face for LD, `1` for DD/Step.
```python
def _inflow_to_moments(self, inflow: tuple[np.ndarray, ...]) -> tuple[np.ndarray, ...]:
    n = self._n_face_moments                              # 370  = per_axis^(d-1)
    if n == 1:                                            # 371  DD/Step → identity
        return inflow
    widened = []
    for face in inflow:                                   # 374
        buf = np.zeros((*face.shape, n))                 # 375  ALL-ZERO 2^{d-1} moment buffer
        buf[..., AVERAGE_MOMENT] = face                  # 376  average ← scalar; SLOPES STAY 0
        widened.append(buf)                              # 377
    return tuple(widened)                                # 378
```
**The transverse face-slope zeroing IS line 375-376**: a fully-zero
`(*face, 2^{d-1})` buffer with only slot 0 written from the scalar inflow. The
docstring (367-368) flags this explicitly: *"the transverse slopes zero ... The
non-vanishing domain-inflow moment trace is the #240 D5b-S4 boundary widening."*
#247's boundary half changes line 375-376 to carry a moment-valued inflow.

**Consumers of `_inflow_to_moments`** (both in `MovingFrontierWindow`, the
production 2-D strategy): `_sweep_interior` (SOLVE) @ **1033**, and
`_loss_action_interior` (APPLY/matvec) @ **1120**. After the walk, the domain-edge
capture is reduced back to slot 0 for the scalar boundary trace at **1068-1069**
(solve) and **1136-1139 / 1320 / 1387** (apply) — `capture = tuple(c[...,
AVERAGE_MOMENT] for c in capture)`. So even with a moment-valued inflow, the
OUTflow trace currently collapses to the average; #247's boundary half also needs
the trace to carry the `2^{d-1}` outflow moments (a `TraceSpace` widening — the
deeper change).

---

## 3. Downstream landing — where a moment-valued Q̂ must LAND

The external source and the scattering slope source MUST agree on the SAME
moment-carrying buffer, because they are summed by typed `+` in the SI rhs.

### 3a. The SI rhs — `q_ext + S.apply(ψ)` (`orpheus/numerics/iteration.py:503-505`)
```python
rhs = q_ext                                    # 503  the moment-lifted composite
for g in self.gains:                           # 504  gains = (S, B)
    rhs = rhs + g.apply(psi)                   # 505  TimedFullField.__add__
```
`q_ext` is the `q_ext_composite` from `_build_fixed_source_rhs` (passed at
`solver.py:2202-2204` via `si.solve(q_ext_composite, ...)`). The Krylov matvec is
the sibling `(L+C).apply − Σ gains.apply` (`iteration.py:578`).

### 3b. The typed `+` REQUIRES matching moment width (the load-bearing constraint)
`TimedFullField.__add__` (`timed_full_field.py:282`) → `Field.__add__`
(`numerics/field.py:262`) → `_check_partner` (`field.py:235-260`): rejects unless
`self.space == other.space` (line 256-260). The bulk space carries the
`SpatialMomentSpace` factor (set via `spatial_moments=per_axis`), so `q_ext` MUST
be lifted to the same `2^d` width as `S.apply(ψ)` or the `+` raises `ValueError`.
**This is WHY `_lift_external_source_to_moments` exists** — the lift is mandatory
for the typed sum, not optional. #247's widened lift keeps this invariant (still
produces a `2^d`-wide bulk) but fills the slope rows.

### 3c. The scattering slope source `S_full` ALREADY populates `(N,ng,*spatial,2^d)`
`ScatteringOperator.apply` (`orpheus/sn/scattering.py`) → `_assemble_per_ordinate_source`
(@ **496-544**). The iso accumulator is sized to the driving flux's moment width
at **532-535**:
```python
spatial_moments = _spatial_moments_of(phi)                    # 532 reads φ̂'s width OFF its space
iso = ScalarSourceSink.zeros_on(mesh, spatial_moments=spatial_moments)  # 533-535
iso = self.add_iso_source(iso, phi)   # 536  Σ_s0·φ̂ into ALL moments (the slope source)
iso = self.add_n2n_source(iso, phi)   # 537
```
The anisotropic `ℓ≥1` half is `_aniso_source_from_moment_values` (@ 434-494,
`R·Λ_{ℓ≥1}·φ̂`). So the scattering source IS a genuine moment-valued
`(N,ng,*spatial,2^d)` buffer whenever `per_axis>1` — **the external `Q̂` joins this
existing buffer**, it does not create a new path.

### 3d. The sweep / matvec landing — `(L+C).solve` and `loss_action`
- SOLVE: `InvertibleOperator.solve` (`operator.py:904`) sweeps via
  `self.loss_representation.sweep(rhs.bulk.values, self.sigma, boundary_buf, ...)`
  at **operator.py:1137-1143**. `rhs.bulk.values` IS the moment-lifted source
  (comment 1121: "by producer contract"). The sweep output carries the trailing
  `2^d` axis (1144-1157).
- MATVEC: `StreamingOperator.apply` (`operator.py:355`) → `loss_representation.loss_action`
  (`operator.py:414`). The matvec twin shares the SAME representation instance (L21).

### 3e. The rank discriminator at the cell — `_CellSolve.cell` (`orpheus/sn/sweep_graph.py`)
This is the seam the brief flags. `_CellSolve.cell` (@ **894-965**), specifically
lines **916-936**, documents EXACTLY how a threaded-slope external source interacts
with `is_moment_valued_by_rank`:
```python
moment_valued = self.scheme.is_multi_moment                  # 916
# ... comment 917-930: the cell SOLVE source can arrive EITHER moment-lifted
# (production solve_sn_fixed_source runs _lift_external_source_to_moments, so
# Q_cells is (N_oct, ng, n_diag, 2^d) — its trailing axis IS the moment axis,
# slopes carry the scattering source Σ_s·φ̂ that needs global→sweep re-signing)
# OR flat (the lower-level sweep API passes (N_oct, ng, n_diag) that _ubld_system
# lifts onto slot 0 itself). Discriminate by RANK ... is_moment_valued_by_rank.
source_is_moment = is_moment_valued_by_rank(Q_cells, reaction_xs)   # 931
psi_avg, psi_out = self.scheme.cell_kernel_batch(
    psi_in=psi_in, s_axes=s_axes, reaction_xs=reaction_xs,
    Q_cells=_reframe(Q_cells, self.moment_frame_signs,
                     is_moment_valued=source_is_moment))     # 934-936 global→sweep re-sign
psi_avg = _reframe(psi_avg, self.moment_frame_signs,
                   is_moment_valued=moment_valued)            # 940-942 sweep→global
```
- `is_moment_valued_by_rank(array, reference)` (`moment_layout.py:71-90`): `True`
  iff `array.ndim > reference.ndim + 1`. A genuine moment buffer
  `(N,ng,*spatial,2^d)` sits at `reference.ndim+2`; a flat `(N,ng,*spatial)` at
  `reference.ndim+1`. RANK-based (NOT trailing-size) so a coincidental `n_diag==2^d`
  anti-diagonal can't mis-fire (#246).
- **Interaction with #247**: a threaded-slope external `Q̂` arrives at the cell as
  a moment-valued `Q_cells` (rank `reference.ndim+2`) EXACTLY like the scattering
  source does today. `is_moment_valued_by_rank` returns `True` → `_reframe` applies
  the `2^d` octant involution (`octant_moment_frame_signs`) to the slope rows
  (global→sweep), as it ALREADY does for the scattering slope source. So #247's
  moment-valued `Q̂` is re-signed CORRECTLY for free — the rank discriminator does
  not distinguish external-slope from scattering-slope (both are genuine `2^d`
  moment sources). **No new branch at the cell.** The only thing the discriminator
  cares about is rank; today the external source is `2^d`-wide but its slope rows
  are zero (frame-invariant), so the reframe is a no-op on them; #247 makes them
  non-zero and the SAME reframe handles them. The matvec twin `_CellResidual.cell`
  (@ 983+, see 1008/1021/1026) hard-codes `is_moment_valued=False` for its
  `Q_cells` (it uses `Q_zero`, line operands at `loss_representation.py:770`) — the
  matvec has NO external source (b−Ax handles it at the scipy level), so #247 does
  not touch `_CellResidual`'s source path.

### 3f. The pure-z degenerate broadcast — `_moment_broadcast_sigma`
`loss_representation.py:494-515` (called @ 698 solve, 789 matvec): `sig[..., None]
if is_moment_valued_by_rank(moment_valued, sig) else sig`. The pure-z balance
`Q/σ_t` (solve) and `σ_t·ψ̄` (matvec) broadcast σ over the moment axis — a
moment-valued external `Q̂` rides this unchanged (it is `2^d`-wide already).

### 3g. The lower-level `_ubld_system` flat-source lift (the "OR flat" arm)
The strategy-level `sweep(Q, ...)` API can receive a FLAT `(N,ng,*spatial)` Q and
lifts it internally onto slot 0 (the `_ubld_system` / `ScanMarch` arms): see
`loss_representation.py:2737-2746` (`lifted[..., AVERAGE_MOMENT] = Q_per_ord`) and
the `AVERAGE_MOMENT` capture-reductions at 1069/1320/1387/1220/2927. This is the
SAME slot-0/slope-zero policy as `_lift_external_source_to_moments` but at a deeper
layer. **#247 caveat**: if a moment-valued source is threaded all the way from
`solve_sn_fixed_source`, the production path (already moment-lifted at 1887) passes
it through `solve`→`loss_representation.sweep(rhs.bulk.values,...)` where
`rhs.bulk.values` is ALREADY `2^d`-wide — so the deeper `_ubld_system` flat-lift is
bypassed (it only fires for operator-free callers passing a flat Q). The two lift
sites must stay consistent on the slope-fill policy (today both zero; #247 changes
only the top one, `_lift_external_source_to_moments`).

**Production sweep strategies** (which cell path runs): 2-D Cartesian =
`MovingFrontierWindow` (`loss_representation.py:970`, drives `_CellSolve`/`_CellResidual`
via `walk_windowed`); 1-D slab/curvilinear = `ScanMarch` (`@1396`, a parallel-prefix
scan, NOT a wavefront — its own moment handling at 2737-2746/2927). LD is currently
2-D-and-1-D-slab; curvilinear LD is rejected (`test_ld_curvilinear_scan_rejected`).

---

## 4. The L17 crosswalk (producer shape → consumer shape at each subsystem boundary)

### External source

| Boundary | Producer shape | Consumer shape | Slope/moment status |
|---|---|---|---|
| Caller → `solve_sn_fixed_source` (`solver.py:1936`) | flat `(N,ng,*spatial)` ndarray OR `TimedFullField` (bulk flat) | typed union | TODAY flat-only; #247 widens to accept `(N,ng,*spatial,2^d)` or a callable producing it |
| `_build_fixed_source_rhs` validation (`solver.py:1875-1880`) | `bulk_values.shape` | `expected=(N,ng,*spatial)` | **HARD-FAILS a `2^d`-rank source today** — the gate to widen |
| `_lift_external_source_to_moments` (`solver.py:1896-1913`) | flat `(N,ng,*spatial)` | `(N,ng,*spatial,2^d)`, slot-0=flat, **slopes ZEROED (1911-1912)** | Q̂=0 ("S4 ansatz"); #247 fills slopes |
| `AngularSourceSink.from_mesh(..., spatial_moments=per_axis)` (`solver.py:1889-1890`) | `(N,ng,*spatial,2^d)` values | typed bulk w/ `SpatialMomentSpace` factor | width must match `S.apply(ψ)` for typed `+` (`field.py:256`) |
| SI rhs `q_ext + S.apply(ψ)` (`iteration.py:503-505`) | both `(N,ng,*spatial,2^d)` | summed `TimedFullField` | scattering ALREADY `2^d` (`scattering.py:532-537`); external joins it |
| `(L+C).solve` → `loss_representation.sweep(rhs.bulk.values,...)` (`operator.py:1137`) | `(N,ng,*spatial,2^d)` | sweep `Q` | passed through to cells |
| cell `_CellSolve.cell` Q_cells (`sweep_graph.py:931-936`) | `(N_oct,ng,n_diag,2^d)` | rank-discriminated → `_reframe` global→sweep | moment-valued ⇒ slopes re-signed (SAME as scattering) |

### Boundary trace

| Boundary | Producer shape | Consumer shape | Slope/moment status |
|---|---|---|---|
| `SNMesh.trace` = `TraceSpace` (`geometry.py:516-519`, prop `:938`) | scalar per face/ordinate/group (NO moment axis) | `boundary.face_view(face)` | scalar-per-face TODAY; #247 needs `TraceSpace` to carry `2^{d-1}` |
| sweep inflow read (`loss_representation.py:806`, `operator.py:1108-1115`) | scalar per face | `inflow` tuple | scalar |
| `_inflow_to_moments` (`loss_representation.py:357-378`) | scalar per face | `(*face,2^{d-1})`, slot-0=scalar, **transverse slopes ZEROED (375-376)** | "#240 D5b-S4 widening"; #247 fills transverse moments |
| outflow capture reduction (`loss_representation.py:1068-1069 / 1136-1139 / 1320 / 1387`) | `(*face,2^{d-1})` | `c[...,AVERAGE_MOMENT]` → scalar | collapses to average; #247 needs the trace to STORE the `2^{d-1}` outflow |
| prescribed-inflow producer `BoundarySourceSink.prescribed_inflow` (`boundary_source_sink.py:188`) | per-face inflow | `q.boundary` (affine BC term) | the non-vacuum inflow channel; today scalar-per-face |

---

## 5. Dependency audit (L20) — typed-union widening is backward-compatible

### `_lift_external_source_to_moments` — 1 production caller, private
| Caller | file:line | passes |
|---|---|---|
| `_build_fixed_source_rhs` | `orpheus/sn/solver.py:1887` | flat `(N,ng,*spatial)` (the ONLY caller; private helper) |

No test or external caller calls `_lift` directly (grep: only the def, the
`_build` call, the `sweep_graph.py:918` comment, and `mms/sn.py:1144` + a
`test_mms_ld_2d.py:320` comment reference it by name). **Widening its slope-fill is
internal** — change the body, keep `(lifted, per_axis)` return contract.

### `solve_sn_fixed_source` — NO production callers (L1-verification-only entry)
70 callers via Nexus; ALL are tests / diagnostics / the regression snapshot
generator — confirming the docstring: *"This entry point exists for L1
verification via MMS, not for engineering problems"* (`solver.py:2032-2036`). Each
passes a **flat** `np.ndarray` or a vacuum `TimedFullField` today. Representative
buckets (file:line — all FLAT/vacuum, the existing path):
- `tests/sn/verification/mms/test_mms*.py` (slab/2d/aniso/curvilinear/het/LD): the
  MMS gates — flat per-ordinate `Q` from `mms/sn.py`. The LD-2-D gates
  (`test_mms_ld_2d.py:46/92/193/232/323/386`) are the ones that WOULD exercise a
  moment `Q̂` once #247 lands (today flat; comment `test_mms_ld_2d.py:320` notes the
  deferral).
- `tests/sn/solve/test_fixed_source_*.py`, `test_si_single_primitive_contract.py:98`,
  `test_fixed_source_2d_equivalence.py:58/116`, `test_d3_admission.py:129/171`,
  `test_2d_anisotropic_windowing.py:101`, `test_affine_carve_bit_identity.py:169` —
  all flat.
- `tests/sn/verification/analytical/test_mms_prescribed_inflow.py:74` — the
  non-vacuum `TimedFullField` boundary path (still scalar-per-face inflow).
- `tests/sn/regression/_generate_snapshots.py:575` — the baseline generator (flat).
- `derivations/diagnostics/*` and `scratch/derivations/diagnostics/*` — flat.

**Backward-compat verdict (CONFIRMED)**: a typed-union widening — accept the
existing flat `(N,ng,*spatial)` OR a new moment-resolved `(N,ng,*spatial,2^d)`
ndarray (or a callable producing the slope rows) — leaves the flat path UNTOUCHED:
- `external_source: np.ndarray | TimedFullField` already exists; the new shape is a
  third recognised ndarray rank, discriminated by `bulk_values.ndim` against
  `expected` at `solver.py:1876`. DD/Step (`per_axis==1`, `face_moment_tail()==()`)
  rejects a `2^d`-axis input naturally (there is no moment axis) — so the new shape
  is only meaningful for LD meshes, and the flat path stays byte-identical (the
  `_lift` early-return at 1909-1910 is the negative control).
- The cell-level rank discriminator (`sweep_graph.py:931`) ALREADY handles a
  moment-valued source identically to the scattering source — no new cell branch.
- The typed `+` invariant (`field.py:256`) is preserved (both operands stay
  `2^d`-wide).

### Note on the eigenvalue path
`solve_sn` / the eigenvalue inner build the FISSION source separately (the
within-group `F` / `solver.py:630/653` use `spatial_moments=...` on their own
accumulators) and do NOT route through `_build_fixed_source_rhs` /
`_lift_external_source_to_moments`. #247 is scoped to the FIXED-SOURCE entry; an
eigenvalue external-source hook does not exist today (the `solve_sn` external hook
is the engineering TODO at `solver.py:2036`).

---

## Key file:line index for the method-implementer

- `solver.py:1875-1880` — flat-shape validation gate to WIDEN.
- `solver.py:1896-1913` (esp. **1911-1912**) — external-source slope-zeroing to FILL.
- `loss_representation.py:357-378` (esp. **375-376**) — boundary transverse-slope zeroing to FILL.
- `geometry.py:516-519` / prop `:938` — `TraceSpace` producer (scalar-per-face; the trace-widening site for the boundary half).
- `loss_representation.py:1068-1069 / 1136-1139` — outflow-capture slot-0 collapse (trace-storage half).
- `scattering.py:496-544` (esp. **532-537**) — the EXISTING `2^d` scattering slope source the external Q̂ joins.
- `iteration.py:503-505` — the `q_ext + S.apply(ψ)` typed sum.
- `numerics/field.py:256` — the space-equality invariant the lift must preserve.
- `sweep_graph.py:916-936` — the rank-discriminated cell re-frame (handles moment Q̂ for free).
- `moment_layout.py:55 / 58-68 / 71-90` — `AVERAGE_MOMENT`, `face_moment_tail`, `is_moment_valued_by_rank`.
