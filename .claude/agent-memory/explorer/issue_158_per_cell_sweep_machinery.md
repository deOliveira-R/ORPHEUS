---
name: issue-158-per-cell-sweep-machinery
description: #158 Increment A — exact slab per-cell SOLVE-walk machinery for OneDimPerCellWalk (dag_walk(ordinate_idx=n, mu_level_idx=None) yields slab CellVisits in chain order; face_view inflow/outflow slots; the psi_in=outgoing_spatial_flux chain; (2,ng) moment source slicing; supports/registry gap). Verified against working tree 2026-06-14.
metadata:
  type: project
---

# #158 Increment A — per-cell SOLVE-walk machinery (slab/Cartesian 1-D)

Verified against the working tree 2026-06-14 (post-#206, `feature/sn-space-angle-tier2`).
SHAPE is durable; the `file:line` below drift — re-confirm with Nexus at pickup.
Companion: [[tier2-space-angle-live-map]] (the registry/seam map).

## The decisive structural fact (answers A directly)

**The SLAB per-cell visit producer is `SNMesh.dag_walk(ordinate_idx=n, mu_level_idx=None)`**
and it is ALREADY the proven recipe — the sweep CACHE BUILDER uses exactly this call to
build slab chains. `SNMesh.dag_walk` (`geometry.py:1052`) → for `CoordSystem.CARTESIAN` →
`_iter_cartesian_visits(ordinate_idx)` (`geometry.py:1289`):

```python
mu_n = float(self.quad.mu_x[ordinate_idx])
cell_indices = range(self.nx) if mu_n >= 0 else range(self.nx - 1, -1, -1)
for i in cell_indices:
    st = self.reduced.streaming_terms(cell_idx=i, direction_idx=ordinate_idx)
    yield CellVisit(cell_idx=i, streaming_terms=st, face_area_downstream=1.0)
```

- For SLAB pass `mu_level_idx=None` (slab/sphere have NO levels; the cylinder is the only
  geometry that requires `mu_level_idx`, `geometry.py:1144`).
- `ordinate_idx=n` IS the global ordinate index for slab/sphere (`geometry.py:1077-1082`).
- Visits come in **chain/sweep order, upstream→downstream**: `range(nx)` for μ≥0,
  `range(nx-1,-1,-1)` for μ<0 (`geometry.py:1303-1304`).
- Each visit carries populated `streaming_terms` (from `reduced.streaming_terms`) and
  `face_area_downstream=1.0` for slab (`CellVisit`, `cell_update.py:181`; the 1.0 neutral
  curvature documented `cell_update.py:222`).
- **PROOF this is the right call**: `SweepCoefficientCache._build` (`sweep_cache.py:243-267`)
  does `level_visits_iter = [(None, n, n) for n in range(N)]` for CARTESIAN, then
  `dag_walk(ordinate_idx=ordinate_idx, mu_level_idx=mu_level_idx)` and reads
  `v.streaming_terms.*` + `v.face_area_downstream` per visit. So `dag_walk` for slab is NOT
  awkward — it is the canonical slab-visit primitive. **Do NOT hand-build CellVisit from
  `geom.chain_idx`** — `dag_walk` is cleaner and already chain-ordered.

## A. Slab cell walk — the call

```python
visits = list(self.mesh.dag_walk(ordinate_idx=n, mu_level_idx=None))   # chain order, slab
for visit in visits:
    i = visit.cell_idx
    ... cell_update.update(visit=visit, total_xs=..., source=..., upstream_state=...)
```

`face_area_downstream=1.0` (slab) ⟹ DD's `update` (`diamond.py:166`) emits
`outgoing_spatial_flux = 2ψ̄ − ψ_in` (the chain link). `is_degenerate` is NEVER True for
slab/sphere — it is set only for cylindrical pure-azimuthal ordinates
(`sweep_cache.py:305`, gated on `face_area_downstream==0.0 and abs_mu<1e-15`). So a slab
walk is purely the non-degenerate chain — there is NO degenerate branch to handle.

## B. BC inflow + outflow persist (slab) — face_view API

`boundary_flux.face_view("xmin")` / `face_view("xmax")` (`_bases.py:402`) return a writable
view (memory-shared with `.values`) of shape **`(N, ng)`** for slab — indexed by **GLOBAL
ordinate index** on axis 0. The proven slab read/write pattern (`loss_representation.py`):

- Setup (`:2475-2478`):
  ```python
  xmin_face = boundary_flux.face_view("xmin")   # (N, ng)
  xmax_face = boundary_flux.face_view("xmax")   # (N, ng)
  inflow_left  = xmin_face   # incoming-ord slots = seeded inflow (forward μ≥0 enters xmin)
  inflow_right = xmax_face   # incoming-ord slots = seeded inflow (backward μ<0 enters xmax)
  ```
- Read inflow seed BEFORE the walk (`:2537-2540`): forward (μ≥0) ordinate `n` seeds from
  `inflow_left[n]` (= `xmin_face[n]`); backward (μ<0) seeds from `inflow_right[n]`
  (= `xmax_face[n]`). Shape `(ng,)` per ordinate.
- Persist outflow AFTER the walk (`:2589-2592`): forward chain's last face → `xmax_face[n]`;
  backward chain's last face → `xmin_face[n]`. Shape `(ng,)`.
- **Bare-sweep contract** (`:2456-2467`): reading the inflow ords (before) and writing the
  outflow ords (after) touch DISJOINT ordinate sets, so aliasing the face view is safe. The
  sweep does NOT re-apply the reflective `R·G` — the caller supplies the seeded inflow and
  reflects after (the sibling `−B`).

For a per-cell walk, the slot mapping is identical: which boundary face is INFLOW vs OUTFLOW
is decided per-ordinate by `sign(μ_n)` — μ≥0: inflow=xmin[n], outflow→xmax[n]; μ<0:
inflow=xmax[n], outflow→xmin[n].

## C. The chain across cells — psi_in threading

For a **non-degenerate slab** ordinate, after each cell you MUST set
`psi_in = result.outgoing_spatial_flux` (shape `(ng,)`). The degenerate-cyl loop
(`loss_representation.py:2687-2705`) does NOT do this — it has NO spatial flow
(`face_area_downstream=0.0` → DD returns `outgoing_spatial_flux=None`, `diamond.py:165-167`),
and reads `result.cell_average_flux` + `result.outgoing_angular_state` only. **The slab walk
is the OPPOSITE**: it has spatial flow, no angular state.

- `psi_in` for **cell 0** (the first visited cell) = the BC inflow face value, `(ng,)` —
  for μ≥0 that is `xmin_face[n]`; for μ<0 that is `xmax_face[n]`.
- The matvec twin confirms the chain shape (`_sweep_direction`, `:2032`):
  `psi_face_in = sn_mesh.cell_update.outgoing_face_from_average(psi_cell, psi_face_in)`. The
  SOLVE-direction analogue is `psi_in = result.outgoing_spatial_flux` (DD: `2ψ̄−ψ_in`).
- The LAST visit's `outgoing_spatial_flux` is the outflow to persist on the downstream
  boundary face (B above).

DD `CellResult` for slab: `cell_average_flux=ψ̄ (ng,)`, `outgoing_spatial_flux=2ψ̄−ψ_in (ng,)`,
`outgoing_angular_state=None` (slab `angular_upstream` is None → `diamond.py:174-179`).

## D. Source shape generality — slicing on the cell axis

`update(...)` takes `source` of shape `(ng,)` for DD; LD wants `(2, ng)` (the moment source).
Cell-axis = the LAST axis of the per-ordinate array. Slicing `[..., i]`:

- `(ng, nx)[..., i]  → (ng,)`     ✓  (DD per-ordinate source)
- `(2, ng, nx)[..., i] → (2, ng)` ✓  (LD per-ordinate moment source)

`...` (ellipsis) handles both ranks cleanly — `QV[n][..., i]` works for `(N, ng, nx)` and
`(N, 2, ng, nx)` source layouts identically. The proven slab `update` call slices
`QV_full[:, i]` (`loss_representation.py:2698`) where `QV_full = QV_per_ord[n]` is `(ng, nx)`;
the rank-agnostic generalization is `QV_per_ord[n][..., i]`.

**V-multiply broadcast** (`loss_representation.py:2447`): `QV_per_ord = Q_per_ord * V[None,None,:]`
for `(N, ng, nx)` with `V` of shape `(nx,)`. For an LD `(N, 2, ng, nx)` layout the multiply is
`Q * V[None, None, None, :]` — V indexes the trailing cell axis. CLEANEST rank-agnostic form:
`Q * V[(...,) + (slice(None),)]`... simpler: broadcast V against the trailing axis with
`V.reshape((1,) * (Q.ndim - 1) + (nx,))` OR just `Q * V` IF V is shaped to the trailing axis.
Since `V` is `(nx,)` and the cell axis is ALWAYS last, `Q * V` broadcasts natively (numpy
right-aligns) for BOTH `(N, ng, nx)` and `(N, 2, ng, nx)` — V `(nx,)` aligns to the trailing
nx axis in both. **`Q * V` (no `[None,...]` prefix) is the rank-agnostic V-multiply.**
(The existing `V[None,None,:]` is the explicit 3-rank form; drop the prefix for rank-agnostic.)

## E. supports / registry placement

Confirmed shape:
```python
@classmethod
def supports(cls, mesh: "SNMesh") -> Compatibility:
    return Compatibility(
        mesh.is_1d and not mesh.cell_update.is_affine_scannable,
        "requires a 1-D mesh with a non-affine-scannable cell-update scheme",
    )
```
- `mesh.is_1d` (`geometry.py:707`) = `ndim == 1` (genuine, not phantom ny=1).
- `mesh.cell_update.is_affine_scannable` (`cell_update.py:553` default False on
  `CellUpdateBase`; `diamond.py:127` True on DD). LD sets it False in Increment A.
- `Compatibility(ok: bool, reason: str)` frozen dataclass (`loss_representation.py:163`).

Registry: `LOSS_REPRESENTATIONS` tuple (`loss_representation.py:1460`) currently
`(CumprodScan, ScanMarch, MovingFrontierWindow, FullFieldWavefront)`. `default_for`
(`:1468`) returns the FIRST whose `supports().ok`. **Slot `OneDimPerCellWalk` AFTER
`CumprodScan`** (CumprodScan claims the affine-scannable 1-D fast path first; OneDimPerCellWalk
catches the non-affine 1-D remainder — the two are mutually exclusive on
`is_affine_scannable`, so order between them is immaterial for correctness, but after-CumprodScan
reads as "the slow 1-D fallback"). Recommended tuple:
`(CumprodScan, OneDimPerCellWalk, ScanMarch, MovingFrontierWindow, FullFieldWavefront)`.

**The gap it fills (confirmed)**: `transport_sweep` (`:1653`) dispatches via
`default_for(sn_mesh).sweep(...)`. A 1-D non-affine mesh is rejected by EVERY current
strategy — CumprodScan (`:701` needs `is_affine_scannable`), ScanMarch (`:1217` needs
`is_affine_scannable`), `_DAGWavefront` family (`:789` needs `is_cartesian AND ndim==2`) —
so `default_for` raises `IncompatibleRepresentation` (`:1494`). OneDimPerCellWalk is exactly
the missing leaf.

## F. gate-2 forcing (DD through the per-cell walk)

DD's `update` works on the SAME per-cell walk — gate-2 forces
`OneDimPerCellWalk(dd_mesh).sweep(...)` DIRECTLY (bypassing `supports`, which would reject DD
since DD `is_affine_scannable=True`). DD's `update` (`diamond.py:138`) takes `source` of shape
`(ng,)`, returns `cell_average_flux (ng,)` + `outgoing_spatial_flux=2ψ̄−ψ_in (ng,)`. The walk
chains `psi_in = result.outgoing_spatial_flux` (C above), so DD-via-PerCell reproduces
CumprodScan's answer.

**CONSTRUCTION GUARD CONFIRMED — gate-2 MUST route around it.** `_LossRepresentation.__post_init__`
(`loss_representation.py:311-318`) DOES run `type(self).supports(self.mesh)` and raises
`IncompatibleRepresentation` if not `ok`. Since DD has `is_affine_scannable=True`,
`OneDimPerCellWalk.supports(dd_mesh).ok` is **False** → `OneDimPerCellWalk(dd_mesh)` RAISES at
construction. gate-2 therefore CANNOT do `OneDimPerCellWalk(dd_mesh).sweep(...)` directly. Options
(pick the cleanest at implementation):
  1. Refactor `sweep` so the walk body is a module-level/free function (or a `_OneDimPerCellWalk`
     inner holder mirroring `_OneDimScanWalk`) that gate-2 calls WITHOUT the guarded wrapper —
     this is the structurally-clean answer and mirrors how `CumprodScan.sweep` delegates to the
     guard-free `_OneDimScanWalk(self.mesh).sweep(...)`. **Recommended**: put the walk in a
     guard-free inner so both the production `OneDimPerCellWalk.sweep` AND gate-2 share it.
  2. Construct via `object.__new__(OneDimPerCellWalk)` + manual `object.__setattr__("mesh", ...)`
     to skip `__post_init__` (hacky — avoid).
The cleanest factoring (option 1) also keeps the matvec twin able to reuse the body — prefer it.

## Recommended skeleton — OneDimPerCellWalk.sweep (slab)

```python
class OneDimPerCellWalk(_LossRepresentation):
    r"""1-D per-cell DAG SOLVE walk — the non-affine-scannable 1-D fallback.

    Walks each ordinate's cell chain calling ``cell_update.update`` per cell
    (the reference per-cell path the affine scan optimizes). Slab/Cartesian
    first; the home for schemes (LD, SC) whose ψ̄ is NOT an affine face-pair
    function so CumprodScan's scan cannot apply.
    """

    @classmethod
    def supports(cls, mesh: "SNMesh") -> Compatibility:
        return Compatibility(
            mesh.is_1d and not mesh.cell_update.is_affine_scannable,
            "requires a 1-D mesh with a non-affine-scannable cell-update scheme",
        )

    def sweep(self, Q, sig_t, boundary_flux, *, initial_guess=None,
              moment_projection=None):
        if moment_projection is not None:
            raise ValueError(
                "OneDimPerCellWalk.sweep: moment output is 2-D Cartesian only."
            )
        mesh = self.mesh
        quad = mesh.quad
        N, nx = quad.N, mesh.nx
        ng = Q.shape[1]                          # (N, ng, nx)
        weights, mu = quad.weights, quad.mu_x
        cell_update = mesh.cell_update
        V = mesh.volumes                          # (nx,)

        QV = Q * V                                # rank-agnostic V-multiply (cell axis = last)
        sig_t_gx = sig_t                          # (ng, nx)

        angular_flux = np.zeros((N, ng, nx))
        scalar_flux = np.zeros((ng, nx))

        xmin = boundary_flux.face_view("xmin")    # (N, ng)
        xmax = boundary_flux.face_view("xmax")    # (N, ng)

        for n in range(N):
            mu_n = mu[n]
            w_n = weights[n]
            forward = mu_n >= 0
            psi_in = (xmin[n] if forward else xmax[n]).copy()   # (ng,) BC inflow seed

            QV_n = QV[n]                           # (ng, nx)  [or (2, ng, nx) for LD]
            last_out = None
            for visit in mesh.dag_walk(ordinate_idx=n, mu_level_idx=None):
                i = visit.cell_idx
                result = cell_update.update(
                    visit=visit,
                    total_xs=sig_t_gx[:, i],            # (ng,)
                    source=QV_n[..., i],               # (ng,) DD / (2, ng) LD
                    upstream_state=UpstreamState(
                        spatial_upstream=psi_in,        # (ng,)
                        angular_upstream=None,          # slab: no angular thread
                    ),
                )
                psi_avg = result.cell_average_flux      # (ng,)
                angular_flux[n, :, i] = psi_avg
                scalar_flux[:, i] += w_n * psi_avg
                psi_in = result.outgoing_spatial_flux   # (ng,) chain link
                last_out = psi_in
            # persist outflow on the downstream face
            if forward:
                xmax[n] = last_out
            else:
                xmin[n] = last_out

        return angular_flux, scalar_flux           # (N, ng, nx), (ng, nx)
```

### Adapt notes
- `UpstreamState`, `CellVisit`, `CellResult` import from `orpheus.sn.spatial.cell_update`.
- `angular_flux[n, :, i] = psi_avg` records the cell-AVERAGE (what DD/the scan store, NOT the
  face flux) — matches the scan's `angular_flux[ords] = psi_avg_cell_order` (`:2579`).
- For LD `angular_flux` per cell is `(ng,)` from `result.cell_average_flux` — CONFIRMED: LD's
  `update` (`linear_discontinuous.py:303-327`) returns `cell_average_flux=ψ̄ (ng,)`,
  `outgoing_spatial_flux=ψ̄+ψ̂ (ng,)`, `outgoing_angular_state=None`. The slope ψ̂ is recoverable
  from the (average, outflow) pair — NO `CellResult` field change, so LD's `CellResult` is the
  SAME rank-`(ng,)` shape as DD's. The angular_flux buffer stays `(N, ng, nx)` and the walk
  skeleton above works for LD UNCHANGED (only the per-cell `source` slice differs: `(2, ng)`).
  LD `update` calls `_require_slab(upstream_state)` (`:319`) — it ALREADY enforces slab-only,
  matching Increment A's slab-first scope.
- `loss_action` / `loss_action_transpose` (matvec twin) for OneDimPerCellWalk: the apply walk
  `_OneDimScanWalk._apply_walk` (`:1873`) currently HARDCODES the WDD
  `cell_balance_for_streaming` denom inline — a non-DD scheme needs its `residual` routed there
  instead. That is the matvec lockstep deliverable (live map §2 "matvec lockstep task"); for
  Increment A the SOLVE walk is the priority, but `loss_action` must at minimum route the
  scheme's `residual` for the Krylov path. If Increment A is SOLVE-only (SI fixed-source),
  `loss_action` can raise NotImplementedError with a clear message until the residual twin lands.

## Anti-surprises
- slab `is_degenerate` is ALL False — no degenerate branch in the slab walk.
- slab `outgoing_angular_state` is None — do NOT thread an angular accumulator (that is the
  curvilinear M-M path only).
- the scan path stores `psi_avg` (cell average) into `angular_flux`, NOT the face flux — match it.
- `Q * V` broadcasts rank-agnostically because the cell axis is ALWAYS trailing; the explicit
  `V[None,None,:]` is the 3-rank-only form.
- gate-2 must bypass `supports` (DD is affine-scannable so `supports` rejects it) — verify the
  base has no construction guard that blocks direct instantiation (the docstring claims a
  `__post_init__` guard but the base body did not show one — RE-CHECK at implementation).
