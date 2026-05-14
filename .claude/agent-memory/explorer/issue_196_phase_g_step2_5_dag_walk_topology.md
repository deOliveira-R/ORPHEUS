---
name: issue-196-phase-g-step2-5-dag-walk-topology
description: Topology / iteration-side characterisation of the Phase G "fold over DAG visits" pattern. Inventories iter_cell_visits / iter_cells_by_direction, the CellVisit + UpstreamState dataclasses, and the 2-D wavefront sweep. Concludes the fold shape is already implemented for 1-D; the 2-D path uses a vectorised batched companion (NOT a per-visit fold) by design.
metadata:
  type: project
---

# DAG-walk topology — what already exists for the Phase G fold (Issue #196 Step 2.5)

Companion to `issue_196_phase_g_step2_5_dd_polymorphism.md` (the per-cell-strategy side).
This memo answers the iteration / topology side: does `dag_walk(mesh, direction) -> Iterator[CellVisit]` already exist, does `CellVisit` carry what a fold needs, and what's the 2-D state of the codebase?

Files inspected (in full where relevant): `orpheus/sn/geometry.py`, `orpheus/sn/sweep.py`, `orpheus/sn/sweep_graph.py`, `orpheus/sn/spatial/cell_update.py`, `orpheus/sn/operator.py` (lines 787, 820, 1027, 1074 — the matvec call sites).

There is NO `Mesh1D` class on the SN side — `SNMesh` is the canonical mesh. The `Mesh1D` / `Mesh2D` names in `orpheus/geometry/` are upstream of `SNMesh` and feed it as input; the SN sweep / iteration API is anchored on `SNMesh`.

---

## Q1 — Cell-iteration primitive

**Both primitives already exist on `SNMesh`** (no `Mesh1D`/`Mesh2D` separate API). The closest match to the proposed `dag_walk(mesh, direction)` is **`SNMesh.iter_cells_by_direction(direction_sign, mu_level_idx=None)`** at `orpheus/sn/geometry.py:586-711`. Its older sibling `SNMesh.iter_cell_visits(ordinate_idx, mu_level_idx=None)` at `orpheus/sn/geometry.py:425-528` is the per-ordinate-indexed primitive that all curvilinear 1-D sweep call sites currently consume.

What each yields:

| Method | Signature | Yields | Notes |
|---|---|---|---|
| `iter_cell_visits` | `(ordinate_idx, mu_level_idx=None)` | `CellVisit` | Per-ORDINATE. Slab/sphere read `mu_x[ordinate_idx]` to pick sweep sign; cylindrical needs the within-level azimuthal index + level. |
| `iter_cells_by_direction` | `(direction_sign, mu_level_idx=None)` | `CellVisit` | Per-SIGN-CLASS. Picks a representative non-degenerate ordinate and delegates to `iter_cell_visits` — single source of truth via `iter_cell_visits` (line 695-697). |

Both yield the same `CellVisit` packet defined at `cell_update.py:278-347`. Geometry dispatch happens inside `iter_cell_visits` (`geometry.py:507-528`):

- `CARTESIAN` → `_iter_cartesian_visits` (`geometry.py:530-553`): `range(nx)` forward for `μ ≥ 0`, reversed for `μ < 0`. `face_area_downstream=None` (slab DD doesn't read it).
- `SPHERICAL` → `_iter_spherical_visits` (`geometry.py:555-584`): same range pattern; `face_area_downstream` resolved to outer for outward, inner for inward.
- `CYLINDRICAL` → `_iter_cylindrical_visits` (`geometry.py:713-776`): per-level; degenerate `|η|<1e-15` iterates forward with `face_area_downstream=None`.

**The three current 1-D iteration sites** (the three `_sweep_1d_*` functions):

| Sweep | File:line | Uses `iter_cell_visits`? |
|---|---|---|
| `_sweep_1d_cumprod` (slab) | `sweep.py:229-364` | **NO** — drives `for n in range(n_half)` and applies the cumprod recurrence vectorised across the whole row (`sweep.py:337-362`). The DD math is inlined here for bit-identity with the pre-Wave-D path; it never enters the per-cell strategy. |
| `_sweep_1d_spherical` | `sweep.py:397-661` | **YES** at `sweep.py:625` — `for visit in sn_mesh.iter_cell_visits(ordinate_idx=n)`. Calls `cell_update.update(visit, ...)` at line 631. |
| `_sweep_1d_cylindrical` | `sweep.py:668-880` | **YES** at `sweep.py:839-841` — `for visit in sn_mesh.iter_cell_visits(ordinate_idx=m_local, mu_level_idx=p)`. Calls `cell_update.update(...)` at line 847. |

The `SNStreamingOperator` matvec path also folds the same generator: `operator.py:787` (`for visit in outward_visits`), `:820` (inward), `:1027` and `:1074` (cylindrical), where the visits come from `iter_cells_by_direction` rather than the per-ordinate variant. This is the Phase G "apply uses the same generator the solve uses" idiom already in production for sphere and cylinder.

So the proposed `dag_walk(mesh, direction)` is essentially `iter_cells_by_direction` — already implemented, already consumed by both the SI sweep AND the matvec, with cylindrical-level support and degenerate-axis handling. The only branch that does NOT yet go through it is the slab `_sweep_1d_cumprod` (bit-identity legacy).

---

## Q2 — `CellVisit` + `UpstreamState` field inventory and fold-compatibility

`CellVisit` (frozen, slotted, `cell_update.py:278-347`):

| Field | Type | Static-per-cell? |
|---|---|---|
| `cell_idx` | `int` | Yes — geometric. |
| `streaming_terms` | `StreamingTerms` | Yes — pure geometry (cell volume, face areas, α, ΔA/w, τ_mm, abs_mu, signed mu). |
| `face_area_downstream` | `float | None` | Yes — sweep-direction-resolved BUT determined once when the visit is emitted, never mutated. |

`UpstreamState` (frozen, slotted, `cell_update.py:354-376`):

| Field | Type | Dynamic across fold? |
|---|---|---|
| `spatial_upstream` | `np.ndarray (ng,)` | YES — the spatial face flux flowing into this cell. |
| `angular_upstream` | `np.ndarray (ng,) | None` | YES — the upstream half-angle state for curvilinear; `None` for slab. |

`CellResult` (frozen, slotted, `cell_update.py:383-411`):

| Field | Type | Becomes the next-cell upstream? |
|---|---|---|
| `cell_average_flux` | `np.ndarray (ng,)` | No — output written to ψ_avg buffer. |
| `outgoing_spatial_flux` | `np.ndarray (ng,) | None` | YES — becomes the next visit's `spatial_upstream` (sweep.py:646, 865). |
| `outgoing_angular_state` | `np.ndarray (ng,) | None` | YES — becomes the next ordinate's `angular_upstream[i]` (sweep.py:639, 855). |

**Fold-compatibility verdict: CLEAN MATCH.** The separation is exactly what a fold needs:

- `CellVisit` is **immutable** (`@dataclass(frozen=True, slots=True)`) — the per-visit packet emitted by the iterator does NOT carry mutable state. There is one foundation test (`tests.sn.spatial.test_cell_update_protocol.TestCellVisitPacket.test_cell_visit_is_frozen`) that pins this.
- `UpstreamState` carries the **dynamic** scalar state that needs threading through the fold step: spatial face flux + (optionally) angular half-angle state. Note this is per-cell upstream, not "one upstream blob for the whole sweep" — `angular_upstream` is indexed by `i` in the calling sites (e.g. `sweep.py:629`).
- The fold step shape `upstream → (visit, upstream) → result` is realised today as:
  ```
  upstream = UpstreamState(spatial_upstream=psi_face,
                           angular_upstream=psi_angle[visit.cell_idx])
  result = cell_update.update(visit, total_xs[i], source[i], upstream)
  psi_face = result.outgoing_spatial_flux           # threads spatial
  psi_angle[visit.cell_idx] = result.outgoing_angular_state  # writes angular
  ```
  (sphere: `sweep.py:625-646`; cylinder: `sweep.py:839-865`).

**What's missing for a pure fold:** the angular half-angle state is **per-cell** (`psi_angle[i]`), so the natural fold accumulator is not the scalar `UpstreamState` but `(psi_face, psi_angle[:])` — i.e. one scalar (spatial) and one array (angular) per step. The sweep currently realises this as `spatial_upstream = psi_face` (threaded as a scalar across visits) plus `angular_upstream = psi_angle[visit.cell_idx]` (read/written by index into a persistent per-cell array). This is NOT a clean Haskell-style fold — the angular state lives outside the accumulator. For the Step 1 sibling `cell_update.residual` method, this is fine: residual takes a probe `cell_avg` and computes the operator action; no state is threaded between residual calls. The fold framing the user proposed is true for `(L+C).solve(q)` IF you accept that `psi_angle[:]` is a persistent array that lives alongside, AND for `L.apply` (no folding required at all — each cell's residual is independent given the input ψ_avg field).

No `CellVisit` mutation is required to make the fold work. The state taxonomy is already correct: static-per-cell → `CellVisit`; dynamic-across-fold → `UpstreamState`.

---

## Q3 — 2-D Cartesian sweep — present, but uses a DIFFERENT shape

**Present.** `_sweep_2d_wavefront` lives at `orpheus/sn/sweep.py:897-1065`. It is the production 2-D path, dispatched from `transport_sweep` for Cartesian 2-D meshes. The infrastructure:

- **Per-octant DAG**: `orpheus/sn/sweep_graph.py:128-260` defines `SweepDependencyGraph`, built once at mesh construction via `SweepDependencyGraph.from_cartesian_2d(nx, ny, label)` (one per `(sign_x, sign_y) ∈ {-1,+1}²`). The graph stores `levels: tuple[(ii, jj), ...]` — for each topological level `k`, the cell indices on the anti-diagonal `local_i + local_j = k`. Built at `sweep_graph.py:241-251` for `k in range(nx + ny - 1)`.
- **Iteration shape**: NOT per-cell. `SweepDependencyGraph.apply` (`sweep_graph.py:264-352`) walks levels with `for ii, jj in self.levels:` and dispatches the WHOLE LEVEL to `cell_update.update_batch(slice_args)` (line 337). `slice_args` is a `SweepCellSlice` packet (`cell_update.py:172-271`) carrying `(N_oct, n_diag, ng)`-shaped arrays — the per-cell DD math vectorises over BOTH the ordinate axis (all ordinates in the current octant) AND the anti-diagonal axis (all cells on this level) simultaneously.
- **Cell update strategy**: same `CellUpdate` Protocol, but uses the optional `update_batch(slice_args)` method (`cell_update.py:665-697`). Default `update_batch` raises `NotImplementedError`; `DiamondDifference` overrides it. Step / LD / EC do not yet. Per the `_sweep_2d_wavefront` module docstring `sweep.py:883-895`: "Parameterising anti-diagonal-vectorised cell updates by the `CellUpdate` Protocol while preserving 1-ULP equality is the open architectural problem for Wave C-extension's LD / EC / Step rollout."

**Key consequence for the Phase G fold framing:** the per-cell `Iterator[CellVisit]` fold pattern that works for 1-D **is intentionally not how 2-D operates today**. 2-D walks levels (not cells) and vectorises across `(N_oct, n_diag)`. So the proposed `dag_walk(mesh, direction) -> Iterator[CellVisit]` does NOT generalise unchanged to 2-D — the 2-D analog yields **levels**, not visits, and the fold body consumes a slice rather than a single packet.

`SNMesh.iter_cell_visits` explicitly rejects 2-D Cartesian meshes (`geometry.py:497-503`): "2-D Cartesian wavefront sweeps use anti-diagonal scheduling, not per-cell visits." The docstring at `geometry.py:489-495` is candid about this: "2-D Cartesian wavefront scheduling is intentionally not encapsulated here — its anti-diagonal vectorisation operates on cell slices, not per-cell visits."

What WOULD need to slot in for a unified `dag_walk`: the level-walking iterator already exists conceptually — `SweepDependencyGraph.levels: tuple[(ii, jj), ...]`. A future `dag_walk(mesh, direction) -> Iterator[CellVisitSlice]` could yield `SweepCellSlice` packets (which already exist!) from this generator, with the slot dimension `n_diag` collapsing to 1 in the 1-D case. But this is a forward-looking unification; nothing in the codebase currently exposes a per-level iterator with a 1-D ⇄ 2-D-uniform signature.

---

## Concrete sketch — what `dag_walk(mesh, direction)` looks like AFTER consolidating the three 1-D iteration sites

The 1-D consolidation is essentially **adopting `iter_cells_by_direction` everywhere a direction-sign suffices**:

```python
# Geometry-aware (on SNMesh; ALREADY EXISTS at geometry.py:586):
def iter_cells_by_direction(self, direction_sign: int,
                             mu_level_idx: int | None = None
                             ) -> Iterator[CellVisit]:
    ...  # picks a representative non-degenerate ordinate, delegates
         # to iter_cell_visits — single source of truth for cell order

# Geometry-blind (lives in transport_sweep / SNStreamingOperator):
def sweep_one_direction(mesh, direction_sign, *,
                        cell_update, source, total_xs, bc, ...,
                        mu_level_idx=None):
    psi_face = bc.inflow(direction_sign)
    for visit in mesh.iter_cells_by_direction(direction_sign,
                                              mu_level_idx=mu_level_idx):
        upstream = UpstreamState(
            spatial_upstream=psi_face,
            angular_upstream=psi_angle[visit.cell_idx],  # None for slab
        )
        result = cell_update.update(visit, total_xs[visit.cell_idx],
                                    source[visit.cell_idx], upstream)
        psi_avg[visit.cell_idx] = result.cell_average_flux
        if result.outgoing_spatial_flux is not None:
            psi_face = result.outgoing_spatial_flux
        if result.outgoing_angular_state is not None:
            psi_angle[visit.cell_idx] = result.outgoing_angular_state
```

The slab `_sweep_1d_cumprod` is the only 1-D site that doesn't already match this shape (it operates on whole rows via cumprod for bit-identity); collapsing it carries the same bit-identity caveat as the DD polymorphism memo §5.1.

**2-D slot.** The fold-over-visits shape DOES NOT extend to 2-D as-is. The 2-D analog is fold-over-LEVELS where each step calls `cell_update.update_batch(slice_args)`. A future `dag_walk` that yields `SweepCellSlice` with `n_diag=1` in 1-D would unify the two — but the user explicitly said NOT to design that now, so this is a forward note only.

---

## Gotchas / mismatches

1. **Carlson half-angle seed and M-M angular redistribution are SEPARATE PASSES that run BEFORE the per-ordinate sweep — confirmed.** Both curvilinear sweeps build their per-cell `psi_angle` array via `carlson_inward_sweep_from_source(...)` BEFORE the `for n in range(N)` ordinate loop (sphere: `sweep.py:557-568`; cylinder: `sweep.py:789-798` for each μ-level). The per-ordinate sweep then reads/writes `psi_angle[i]` cell-by-cell inside its fold over visits. So the Step 2 canonical-closure framing (seed before, fold during) is consistent with the actual code. Within one ordinate, the cell update reads the upstream half-angle from `psi_angle[i]` and writes the downstream back into `psi_angle[i]` for the NEXT ordinate to consume — this is the M-M angular DAG and it is per-ordinate, not per-cell, hence "separate pass" relative to the per-cell spatial fold.

2. **`iter_cell_visits` vs `iter_cells_by_direction` co-existence.** Both are public, both yield the same packet type, but they take different arguments. The SI sweep call sites use `iter_cell_visits(ordinate_idx=n)` because they need to bind to a specific ordinate (`mu_x[n]`, `weights[n]`, `Q_aniso[n]`). The matvec call sites use `iter_cells_by_direction(±1)` because the sweep-frame apply operates on all outward / inward ordinates simultaneously. Both work because `iter_cells_by_direction` internally delegates to `iter_cell_visits` (`geometry.py:695-711`) — single source of truth for cell ordering. A Phase G unification could collapse these to one method (`iter_cells_by_direction` with an optional `ordinate_idx` for source/weight resolution at the call site), but that's downstream of the four-operator algebra work.

3. **`angular_upstream` is per-cell, not per-fold-step.** The fold framing in the user's brief has `upstream` as a single accumulator threaded through the visits. That's clean for the spatial face flux (slab + non-degenerate curvilinear), but the angular half-angle state is a persistent **array indexed by `cell_idx`** that lives outside the fold. The sweep reads `psi_angle[visit.cell_idx]` to populate `UpstreamState.angular_upstream` and writes back the result. This is principled — the M-M redistribution couples ordinate-to-ordinate at each cell independently, not cell-to-cell along the spatial DAG — but it means a pure Haskell-style `foldl` does NOT capture the angular state. The cleaner framing is "fold over visits, with a side-effect on `psi_angle[i]`" (which is what the sphere and cylinder sweeps do today).

4. **Slab `_sweep_1d_cumprod` does NOT fold over `iter_cell_visits`.** It operates on whole-row cumprod / cumsum vectorised ops (`sweep.py:367-390`). The bit-identity contract (documented at `diamond.py:173-194` and re-explained in the DD polymorphism memo §5.1) is the reason. Collapsing this into the unified fold inherits the bit-identity caveat from that memo — `np.array_equal` → `np.allclose(rtol=1e-13)` re-baselining is needed for the slab tests. The cumprod IS algebraically a fold of the same per-cell DD update; the project chose to keep the legacy operation order for migration safety.

5. **2-D wavefront uses a TOTALLY different generator: `SweepDependencyGraph.levels`**, not `iter_cell_visits`. The level packet is `(ii, jj)` arrays (cell indices on the anti-diagonal); the consumer (`SweepCellSlice` + `cell_update.update_batch`) is shape-`(N_oct, n_diag, ng)`. So the per-cell `for visit in mesh.iter_cells_by_direction(...)` fold pattern that the brief proposes is intentionally absent in 2-D and `iter_cell_visits` explicitly raises on 2-D meshes (`geometry.py:497-503`). The shipped 2-D fold IS still a forward-substitution on (L+C) — just batched along the anti-diagonal — but the iterator shape and the cell-update signature differ.

6. **`face_area_downstream` is `None`-coded for "no spatial flow" today**, NOT `0.0`. Slab carries `None` because slab DD doesn't read face areas at all (chord_length plays the dual role); the cylindrical pure-azimuthal degenerate case carries `None` because the cell has no radial face flow on that ordinate. The DD polymorphism memo §5.5 proposes flipping to `0.0` so the strategy reads ONE number unconditionally. That is a CellVisit field change orthogonal to the fold-shape question; today's iterator emits `None` and downstream consumers check `is None`.
