---
name: issue-196-phase-g-replan-sn-structure
description: Structural map of the SN module for the Phase G replan of Issue #196. Confirms SNMesh IS the discrete-ordinates phase space, surfaces the half-shipped DiscreteOrdinatesPhaseSpace/StreamingOperator/CollisionOperator from the prior Phase G attempt, and lists the anti-recommendations the Step 2/3 agent must internalise to avoid re-creating types that already exist.
metadata:
  type: project
---

## TL;DR — the one-line answers

- **SNMesh IS the SN-augmented discrete-ordinates phase space.** It
  is the layer that fuses (Mesh1D / Mesh2D) + AngularQuadrature +
  precomputed streaming / connection coefficients + cell-update
  strategy + pole-angular closure + resolved boundary operators into
  one mathematically-honest object. `n_groups` is the only thing it
  does **not** already carry — and that single missing scalar is the
  entire surface area of the type the previous Phase G attempt
  reinvented as `DiscreteOrdinatesPhaseSpace`.

- **A `DiscreteOrdinatesPhaseSpace` already exists in `orpheus/sn/streaming.py`.**
  It is a 3-field frozen dataclass `(sn_mesh, quadrature, n_groups)`.
  Two of the three fields are properties already on SNMesh
  (`SNMesh.mesh`, `SNMesh.quad`). It is **not** a Cardinal Rule 2
  victim of duplication accidentally — it was a deliberate Step 2
  shim, but it is exactly the wrapper the user rejected. Step 2 of
  the replan must NOT keep this type as designed; SNMesh + ng (a
  scalar) is the same information.

- **The Phase G Step 2 work landed half-shipped on
  `refactor/sn-operator-algebra`.** `orpheus/sn/streaming.py` carries
  `StreamingOperator`, `CollisionOperator`,
  `solve_within_group_sn_curvilinear`, `build_sn_operators` and is
  wired into the curvilinear SI sweep
  (`orpheus/sn/sweep.py:255-296`). It has a **broken TYPE_CHECKING
  import** (`from .mesh import SNMesh` — there is no `orpheus/sn/mesh.py`;
  SNMesh lives in `orpheus/sn/geometry.py`). The body of
  `StreamingOperator.apply` still delegates to the procedural
  `transport_operator_matvec_*` functions verbatim — so the
  reorganisation is currently only the **outer shell**, the inner
  procedural code that Manifestation #7 traces to was never replaced.

---

## 1. SNMesh — full surface area

File: `orpheus/sn/geometry.py:51-868`

### 1.1 Constructor

```python
class SNMesh:
    BOUNDARY_OPERATOR_REGISTRY: ClassVar[dict[str, type[BoundaryTraceLaw]]] = {
        "vacuum": VacuumInflow,
        "reflective": ReflectiveBoundary,
    }

    def __init__(
        self,
        mesh: Mesh1D | Mesh2D,
        quadrature: AngularQuadrature,
        cell_update: CellUpdate | None = None,
        pole_angular_closure: PoleAngularClosure | None = None,
    ) -> None:
```
Source: `geometry.py:125-191`.

| Parameter              | Type                      | Default                                          | Physical / architectural meaning                                                                                          |
| ---------------------- | ------------------------- | ------------------------------------------------ | ------------------------------------------------------------------------------------------------------------------------- |
| `mesh`                 | `Mesh1D \| Mesh2D`        | required                                         | Base spatial mesh — coordinate system, widths, material map, declared face BCs (`bc_left`, `bc_xmin`, ...).               |
| `quadrature`           | `AngularQuadrature`       | required                                         | Angular discretisation Ω. Carries `mu_x / mu_y / mu_z`, weights, `reflection_index(axis)`, `level_indices` (cylindrical). |
| `cell_update`          | `CellUpdate \| None`      | `DiamondDifference()` (geometry.py:142-144)      | Per-cell spatial closure strategy (DD; future LD / EC / Step). Consumed by `transport_sweep`, NOT by the matvec apply.    |
| `pole_angular_closure` | `PoleAngularClosure \| None` | `MorelMontryAngularSweep()` (geometry.py:188-191) | Angular redistribution strategy for curvilinear (sphere / cylinder). Consumed by both apply matvec AND solve sweep.       |

### 1.2 Attributes set in `__init__`

All inline in `__init__` between `geometry.py:132-191`:

| Attribute                | Source line   | Carries                                                                                | Consumed by                                                                  |
| ------------------------ | ------------- | -------------------------------------------------------------------------------------- | ---------------------------------------------------------------------------- |
| `self.mesh`              | `:132`        | Base Mesh1D / Mesh2D handle.                                                           | Everywhere (call sites read `sn_mesh.mesh.volume_measure`, `mat_ids`, ...).  |
| `self.quad`              | `:133`        | The angular quadrature.                                                                | Every operator (`L.apply`, `S.apply`, sweep, BC realiser).                   |
| `self.cell_update`       | `:142`        | DD / future LD / EC / Step strategy.                                                   | `transport_sweep` (solve path) ONLY.                                         |
| `self.pole_angular_closure` | `:188`     | M-M / Bailey / legacy-τ strategy.                                                      | Both `transport_operator_matvec_*` (matvec apply) and `transport_sweep`.     |

Spatial scalars (set in the coord-system dispatch at `:194-207`):

| Attribute       | Shape             | Meaning                                                          |
| --------------- | ----------------- | ---------------------------------------------------------------- |
| `self.nx`       | int               | Cells along x (1-D) or x-axis (2-D).                             |
| `self.ny`       | int               | Cells along y (2-D); `1` on 1-D meshes.                          |
| `self.dx`       | `(nx,)`           | Radial / x widths.                                               |
| `self.dy`       | `(ny,)`           | y widths (or `np.array([1.0])` on 1-D).                          |
| `self.mat_map`  | `(nx, ny)`        | Per-cell material id.                                            |
| `self._volumes` | `(nx, ny)`        | Per-cell volume (exposed as `self.volumes` property).            |

Streaming-stencil dispatch (`:226-257`):

| Attribute                                       | Coord  | Source                                                                  | Meaning                                                            |
| ----------------------------------------------- | ------ | ----------------------------------------------------------------------- | ------------------------------------------------------------------ |
| `self.reduced`                                  | all 1D | `slab_streaming` / `spherical_streaming` / `cylindrical_streaming` (Wave B Issue #6) | Canonical `ReducedStreamingOperator` accessor with `streaming_terms(cell_idx, direction_idx, mu_level_idx)` packets — face areas, ΔA/w, α_half, τ_mm. |
| `self.reduced`                                  | 2-D Cartesian | `None`                                                              | 2-D Cartesian uses inlined `streaming_x` / `streaming_y` instead.  |
| `self.curvature`                                | sphere | `"spherical"`                                                          | Dispatch key consumed by `SNStreamingOperator._ensure_eq_map` and `_build_rhs_*`. |
| `self.curvature`                                | cyl    | `"cylindrical"`                                                        | Same dispatch key.                                                 |
| `self.curvature`                                | slab/2D| `None`                                                                 | Same.                                                              |
| `self.streaming_x`, `self.streaming_y`          | 2-D Cart | `_setup_cartesian` (`:780-803`)                                       | `(N, nx)` and `(N, ny)` precomputed `2|μ|/Δ` for DD denominators.  |
| `self.sweep_graphs`                             | 2-D Cart | `_setup_cartesian` (`:817-823`)                                       | Per-octant `SweepDependencyGraph` (4 octants). Not built for curvilinear (1-D wavefronts are trivial). |

Boundary resolution (`:259-407`, then attributes on `_resolve_bcs`):

| Attribute       | Set by                  | Carries                                                          |
| --------------- | ----------------------- | ---------------------------------------------------------------- |
| `self._inflow_trace`, `self._outflow_trace` | `_resolve_bcs:290-307` | Per-mesh-per-quad `InflowTraceSpace` / `OutflowTraceSpace`. Built once. |
| `self.bc_left`, `self.bc_right`           | `_resolve_one(...)`     | `_BoundBoundaryOperator` wrapping a 1-arg `LinearOperator` (per realiser dispatch). 1-D meshes only. |
| `self.bc_xmin`, `self.bc_xmax`, `self.bc_ymin`, `self.bc_ymax` | `_resolve_one(...)` or 1-D-bootstrap aliases | Same shape — but 1-D meshes alias `bc_xmin`↔`bc_left`, `bc_xmax`↔`bc_right`, and the y-faces are realised `ReflectiveBoundary(axis="y", α=1.0)` placeholders (degenerate; identity on GL1D quadratures). |

### 1.3 `@property` surface

| Property             | Source line       | Returns                                                                                                       |
| -------------------- | ----------------- | ------------------------------------------------------------------------------------------------------------- |
| `volumes`            | `:411-414`        | `(nx, ny)` cell volumes.                                                                                      |
| `is_1d`              | `:416-419`        | `True` iff `ny == 1`.                                                                                         |
| `face_areas`         | `:839-849`        | **Deprecated** — proxies `self.reduced.face_areas` with `DeprecationWarning`.                                  |
| `delta_A`            | `:851-861`        | **Deprecated** — proxies `self.reduced.delta_A` with `DeprecationWarning`.                                    |

(Note: `geometry.py:863-868` records that six more Wave-D-Round-1.1
property shims — `alpha_half`, `redist_dAw`, `tau_mm`, and the
cylindrical per-level analogues — were RETIRED in Wave E Round 2.
Callers read `self.reduced.<name>` directly. Cleanup candidate: the
two remaining `DeprecationWarning` shims (`face_areas`, `delta_A`) are
also unused in production paths — see §6.)

### 1.4 Sweep DAG traversal methods

| Method                                              | Source                | Yields                                                       | Consumed by                                                                                |
| --------------------------------------------------- | --------------------- | ------------------------------------------------------------ | ------------------------------------------------------------------------------------------ |
| `iter_cell_visits(ordinate_idx, mu_level_idx=None)` | `:425-528`            | `CellVisit(cell_idx, streaming_terms, face_area_downstream)` | The unified sweep (`sweep.py`) and curvilinear apply matvecs.                              |
| `iter_cells_by_direction(direction_sign, mu_level_idx=None)` | `:586-711`   | Same                                                         | The Phase C sweep-frame apply matvec (`operator.py:751-765`, `:1006-1071`).                |
| `_iter_cartesian_visits(ordinate_idx)`              | `:530-553`            | Slab visits, no face area packed.                            | Internal helper for `iter_cell_visits`.                                                     |
| `_iter_spherical_visits(ordinate_idx)`              | `:555-584`            | Spherical visits with downstream face area.                  | Internal.                                                                                   |
| `_iter_cylindrical_visits(ordinate_idx, mu_level_idx)` | `:713-776`         | Per-level visits with pure-azimuthal degenerate branch.      | Internal.                                                                                   |

### 1.5 BOUNDARY_OPERATOR_REGISTRY + `_resolve_one`

| Method            | Source        | Behaviour                                                                                                                                          |
| ----------------- | ------------- | -------------------------------------------------------------------------------------------------------------------------------------------------- |
| `_resolve_bcs(mesh)` | `:264-356` | Builds the inflow / outflow trace spaces; dispatches to `_resolve_one` per face; aliases bc_xmin/left and produces degenerate y-face shims on 1-D. |
| `_resolve_one(bc, face)` | `:365-407`| Looks up `bc.kind` in `BOUNDARY_OPERATOR_REGISTRY`, constructs the trace-law instance, builds an `SNMethodSpace.for_face(...)`, hands both to `SNBoundaryRealizer().realize`, wraps the resulting 1-arg `LinearOperator` in `_BoundBoundaryOperator` with the `kind` tag. |
| `_face_names_for_mesh(mesh)` | `:358-363` | Returns `("left", "right")` for Mesh1D, `("xmin", "xmax", "ymin", "ymax")` for Mesh2D. |

The registry is intentionally minimal (only `"vacuum"` and
`"reflective"` map — `_resolve_one` raises on anything else). 5 other
BC kinds (`white`, `periodic`, `albedo`, `prescribed_inflow`, `mixed`)
are handled by `SNBoundaryRealizer` but not yet exposed at the SN
mesh level (geometry.py:117-123 comment explains the gating: they
need SN-sweep-side wiring for cycles / periodic / etc.).

### 1.6 One-paragraph summary — what SNMesh IS

`SNMesh` is **the SN-augmented mesh = the discrete-ordinates phase
space**. It composes the spatial mesh (`mesh`), the angular
discretisation (`quad`), the precomputed streaming connection
coefficients (`reduced`, `streaming_x/y`, `sweep_graphs`), the
boundary-condition realisations (`bc_left/right/xmin/...`), and the
two solver strategies (`cell_update`, `pole_angular_closure`) into
one typed object. The continuous transport phase-space
:math:`V = X \otimes \Omega \otimes G` corresponds to (SNMesh, ng) —
where SNMesh carries X⊗Ω plus all auxiliaries, and `ng` is the
single missing scalar an operator constructor needs to know how
many groups the cross-section tensors carry. There is no
"missing" `DiscreteOrdinatesPhaseSpace` type: SNMesh **is** it.

---

## 2. Four-operator architecture surface today

### 2.1 The user-rejected wrapper that already exists

**`orpheus.sn.streaming.DiscreteOrdinatesPhaseSpace`**
(`streaming.py:120-141`)

```python
@dataclass(frozen=True, slots=True)
class DiscreteOrdinatesPhaseSpace:
    sn_mesh: "SNMesh"
    quadrature: "AngularQuadrature"
    n_groups: int

    @property
    def geometry(self) -> str | None:
        return getattr(self.sn_mesh, "curvature", None)
```

It is a 3-field NamedTuple-like wrapper around `(sn_mesh, quadrature,
n_groups)`. **`quadrature` duplicates `sn_mesh.quad`**;
**`geometry` duplicates `sn_mesh.curvature`**. Only `n_groups` is
new information, and it is a single int. This is exactly the wrapper
the user rejected ("SNMesh IS the DiscreteOrdinatesPhaseSpace").

### 2.2 The four operators today

| Symbol | Class                                                | File:line                  | LinearOperator?       | Constructor signature                                          | Apply takes               | Apply returns          |
| ------ | ---------------------------------------------------- | -------------------------- | ----------------------- | -------------------------------------------------------------- | ------------------------- | ---------------------- |
| L (legacy) | `SNStreamingOperator(LinearOperatorMixin)`       | `operator.py:1114-1473`    | Yes (`{apply, solve, apply_transpose}`) | `(sn_mesh: SNMesh, sig_t: np.ndarray)` (dataclass)             | packed `(n_unknowns,)`    | packed `(n_unknowns,)` |
| L (Phase G) | `StreamingOperator(LinearOperatorMixin)`       | `streaming.py:250-408`     | Yes (`{apply}`)        | `(V: DiscreteOrdinatesPhaseSpace, boundary: BoundaryOperator)` (dataclass) | packed `psi, *, sigma_t=` | packed                 |
| C (Phase G) | `CollisionOperator(LinearOperatorMixin)`       | `streaming.py:416-512`     | Yes (`{apply, solve, apply_transpose}`) | `(V: DiscreteOrdinatesPhaseSpace, sigma_t: np.ndarray)` (dataclass) | packed or `(N,nx,ny,ng)` | same shape             |
| S      | `ScatteringOperator(LinearOperatorMixin)`            | `scattering.py:202-524`    | Yes (`{apply}`)        | dataclass: `n_ordinates, nx, ny, ng, scattering_order, sig_s, sig2, sig_s0, Y, weights, cells_by_mat` plus `LegendreMomentScattering` for the Pℓ block | various                  | various                |
| F      | `FissionOperator(LinearOperatorMixin)`               | `fission.py:82-165`        | Yes (`{apply}`)        | dataclass: `(chi, sig_p)` shape `(nx, ny, ng)` each            | `(nx, ny, ng)`            | `(nx, ny, ng)`         |

Construction call sites:
- `SNSolver.__init__` builds the (L, S, F) triple at `solver.py:240-262`.
  This is the canonical operator-triple construction site. Note: it
  uses **`SNStreamingOperator`** (the legacy form), not the Phase G
  `StreamingOperator`. The two are parallel right now.
- `sweep.py:255-296`'s `_curvilinear_sweep_via_canonical_operator`
  uses **`StreamingOperator`** (the Phase G form), via
  `build_sn_operators(V, materials)` (streaming.py:697-740). This is
  the only call site that consumes the Phase G operators in
  production today.

### 2.3 CollisionOperator — does it exist as a standalone unit?

YES — `streaming.py:416-512` ships `CollisionOperator`. It is
separate from `SNStreamingOperator.apply`, which still applies σ_t
internally via `transport_operator_matvec_*` (the σ_t term is the
last line of each `lhs[:, ks] = streaming + redistribution +
collision` assembly in `operator.py:799, 832, 1039, 1086`).

So **C** does not yet have a single canonical home: the Phase G
`CollisionOperator` exists but `L = SNStreamingOperator` still
folds it inside `L`. The current production view is **(L+C) as one
indivisible operator** (`SNStreamingOperator`); the Phase G view is
**L + C as two operators composed**.

### 2.4 The `transport_operator_matvec_*` procedural family

| Function                                       | File:line          | Apply takes                                                | Apply returns                  | Called by                                                                  |
| ---------------------------------------------- | ------------------ | ---------------------------------------------------------- | ------------------------------ | -------------------------------------------------------------------------- |
| `transport_operator_matvec` (Cartesian)        | `operator.py:412-458` | `(solution, eq_map, quad, sig_t, nx, ny, ng, dx, dy, bc_xmin=, bc_xmax=, bc_ymin=, bc_ymax=)` | packed `(n_unknowns,)`         | `SNStreamingOperator.apply` (`operator.py:1353-1360`) AND `StreamingOperator.apply` (`streaming.py:386-393`) |
| `transport_operator_matvec_spherical`          | `operator.py:571-838` | `(solution, eq_map, quad, sig_t, nx, ng, face_areas, volumes, alpha_half, redist_dAw, tau_mm, *, sn_mesh, bc_outer, pole_angular_closure)` | packed | Same two call sites (`SNStreamingOperator.apply` + `StreamingOperator.apply`). |
| `transport_operator_matvec_cylindrical`        | `operator.py:851-1107` | analogous, per-level arrays                            | packed                         | Same two call sites.                                                       |

Layout: ψ comes in as a 1-D packed vector with `eq_map.n_unknowns =
n_eq * ng` entries; `n_eq` is the number of `(ordinate, cell)` pairs
that are **unknowns** (i.e., excluding incoming-at-reflective-boundary
slots that are determined by the BC). The interior body uses
`solution_to_angular_flux*` to scatter the packed vector into a
`(ng, N, nx, ny)` cell-centre angular flux, then runs WDD propagation
+ pole-angular-closure + collision, then re-packs.

---

## 3. The two production sweep / matvec paths

### 3.1 Sweep paths (solve via `transport_sweep`)

Entry point: `orpheus.sn.sweep.transport_sweep(Q, sig_t, sn_mesh,
psi_bc, Q_aniso=None)` (`sweep.py:112-153`).

Dispatch:
- `sn_mesh.reduced.requires_upstream_angular_state` False → `_cartesian_sweep`.
- `... True` → `_curvilinear_sweep` (sweep.py:202-252).
- `sn_mesh.reduced is None` (2-D Cartesian) → fallthrough to Cartesian.

| Path                                                  | File:line        | What it does                                                                                                                                                            |
| ----------------------------------------------------- | ---------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `_cartesian_sweep`                                    | `sweep.py:160-195` | Gates on `(GL1D + DD + isotropic)` → 1-D cumprod fast path. Otherwise → 2-D wavefront.                                                                                  |
| `_sweep_1d_cumprod`                                   | `sweep.py:303-...` | Bit-identical cumprod path: vectorised DD recurrence with `streaming_x` denominator, calls `sn_mesh.bc_left.apply` / `sn_mesh.bc_right.apply` for face reflection / vacuum. |
| `_sweep_2d_wavefront`                                 | (further in sweep.py) | Anti-diagonal scheduling; per-cell uses `cell_update.update(visit, sig_t_cell, source, upstream)` for DD math.                                                          |
| `_curvilinear_sweep`                                  | `sweep.py:202-252` | **NEW Phase G**: dispatches to `_curvilinear_sweep_via_canonical_operator`. The historical `_sweep_1d_spherical` / `_sweep_1d_cylindrical` procedural sweeps were replaced. |
| `_curvilinear_sweep_via_canonical_operator`           | `sweep.py:255-296` | Builds `V = DiscreteOrdinatesPhaseSpace(sn_mesh, quad, ng)`, `L = StreamingOperator(V, boundary=sn_mesh.bc_right)`, `C = CollisionOperator(V, sigma_t=sig_t)`, then calls `solve_within_group_sn_curvilinear(L, C, Q, ...)`. |

### 3.2 Matvec apply paths (Krylov)

Entry points:
- `SNStreamingOperator.apply(psi)` (`operator.py:1287-1360`) —
  dispatches by `sn_mesh.curvature` to one of the three
  `transport_operator_matvec_*`.
- `StreamingOperator.apply(psi, *, sigma_t)` (`streaming.py:310-393`)
  — dispatches by `self.V.geometry` (a property returning
  `sn_mesh.curvature`) to the same three procedural matvecs.

Both reach the same procedural matvecs. The Phase G shell does not
yet hold its own implementation — it delegates.

### 3.3 Top-level fixed-source dispatch

In `solver.py`:

| Function                                       | File:line              | Role                                                                                                                                                                                            |
| ---------------------------------------------- | ---------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `solve_sn_fixed_source(...)`                   | `:924-1064`            | Public entry. Auto-flips curvilinear to `"krylov"`; Cartesian default `"source_iteration"`. Builds `SNMesh`, `SNSolver`, dispatches to `_solve_fixed_source_si` or `_solve_fixed_source_krylov`. |
| `_solve_fixed_source_si(solver, sn_mesh, ext_src, ...)` | `:1067-1129`  | Cartesian default. Bit-identical Wave A-D source iteration loop wrapping `transport_sweep`.                                                                                                     |
| `_solve_fixed_source_krylov(solver, sn_mesh, ext_src, ...)` | `:1132-1289` | Curvilinear-default. Outer source iteration; inner GMRES on `solver.L.apply` (`SNStreamingOperator`); decode via `solution_to_angular_flux_{spherical,cylindrical,...}`.                       |
| `solve_sn(materials, mesh, quadrature, ...)`   | `:824-921`             | Eigenvalue. Builds SNMesh + SNSolver, calls `power_iteration(solver, max_iter=...)`. Final sweep for angular flux. Returns `SNResult`.                                                          |
| `SNSolver.solve_fixed_source(...)`             | `:281-293`             | Within-eigenvalue inner solve. Dispatches to `_solve_source_iteration` or `_solve_krylov` based on `inner_solver` field set at `__init__`.                                                       |
| `SNSolver._solve_krylov(...)`                  | `:424-551`             | GMRES on `solver.L.apply` (`SNStreamingOperator`) with `_make_sweep_preconditioner`. Builds packed RHS via `_build_rhs_cartesian` / `_build_rhs_spherical` / `_build_rhs_cylindrical`.            |
| `SNSolver._solve_source_iteration(...)`        | `:374-420`             | Sweep-driven within-group iteration via `transport_sweep`.                                                                                                                                       |

`solve_sn` returns `SNResult(keff, keff_history, angular_flux,
scalar_flux, geometry, quadrature, eg, elapsed_seconds)`
(`solver.py:100-115`).

---

## 4. SNMesh ↔ BCs relationship

### 4.1 The realised contract

After `SNMesh.__init__`, every face has a `_BoundBoundaryOperator`
shim wrapping a 1-arg `LinearOperator` produced by
`SNBoundaryRealizer.realize(law, method_space)`:

```python
sn_mesh.bc_xmin.apply(psi_face_out)  # → psi_face_in
```

The shim carries a `.kind` tag for legacy string comparisons
(`bc.kind == "vacuum"`). For specular reflective the realised op is a
`PermutationOperator(quad.reflection_index("x"))`; for vacuum it is
`IncomingOrdinateMaskTensor(inflow_indices, n_ordinates, axis=0)`;
for white it is `AngularAverageOperator.from_quadrature(...)`. See
`boundary_realizer.py:122-194` for the full dispatch table.

### 4.2 Are BCs first-class `LinearOperator`s today?

**YES.** The realised op IS a `LinearOperator` (specifically one of the
Wave-0/Wave-1 primitives in `orpheus.numerics.operator`). It exposes
`.apply(x)`. The Wave-9 single-arg migration completed: every SN call
site reads `bc.apply(psi)` with the quadrature ALREADY BOUND inside
the realisation. `_BoundBoundaryOperator` is a thin wrapper that
preserves the `LinearOperator` contract while also carrying the
`.kind` string for legacy comparisons.

### 4.3 Production read sites for `bc_*`

| Read site                                                  | Behaviour                                                                                                                                            |
| ---------------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------- |
| `_sweep_1d_cumprod` (`sweep.py:386-388`)                   | `bc_left_obj = sn_mesh.bc_left; bc_right_obj = sn_mesh.bc_right` → `bc.apply(psi_face_full)`.                                                        |
| `transport_operator_matvec_spherical` (`operator.py:713, 815`) | `outer_inflow_estimate = bc_outer.apply(fi[:, :, -1, 0].T)` (pre-sweep BC trace), then `inflow_full = bc_outer.apply(outflow_at_boundary.T)` (Phase C boundary trace law). |
| `transport_operator_matvec_cylindrical` (`operator.py:954, 1048`) | Analogous — first for Carlson context, then for inward sweep seed.                                                                                  |
| `transport_operator_matvec` Cartesian (`operator.py:308-335`) | Per-face `bc_xmin.apply(outgoing_xmin.T)` → fill incoming slots.                                                                                    |
| `SNStreamingOperator.apply` (`operator.py:1335-1359`)      | Passes `sn_mesh.bc_xmin/xmax/ymin/ymax` or `sn_mesh.bc_right` through to the matvec function as kwargs.                                              |
| `StreamingOperator` (Phase G) construction (`sweep.py:288`) | `L = StreamingOperator(V=V, boundary=sn_mesh.bc_right)` — picks JUST `bc_right` for curvilinear, where it represents the outer face.                |

### 4.4 The contract — does `sn_mesh.bc_outer.apply(psi_face)` produce the inflow trace?

YES, post-Wave-9. The 1-arg apply contract is: takes the outflow face
flux (full ordinate vector at the face), returns the inflow face
flux (same shape, with the inflow-ordinate slots populated per the
law). Vacuum returns zero for inflow ordinates and preserves outflow
ordinates (per Wave-5 semantic correction at
`boundary_realizer.py:18-25`). Specular returns `out[ref_x[n]]` at
the partner ordinate. White returns the angular average.

---

## 5. Cross-references to relevant tests

### 5.1 Regression snapshots (the algebra-of-record gate)

Location: `tests/sn/regression/snapshots/`

Snapshot inventory (from `tests/sn/regression/snapshots/*.npz`):
- `2d_1g_LS4_dd_15x15.npz`
- 6 × `2d_octant_equivalence_*.npz`
- 3 × `cyl_*.npz`
- Plus spherical snapshots (further in directory listing).

Driver: `tests/sn/regression/test_dd_regression.py:32-46` —
re-runs each case, asserts `np.array_equal` against the frozen
`.npz`.

### 5.2 Variant α cross-checks (Gate 4.2 / Phase E sentinel)

File: `tests/sn/test_phase_c_crosscheck.py`

Notable cases:
- Lines 151-291: SN snapshot k_eff vs `trajectory_resolvent` Variant α
  reference for the 5 P0 snapshots that admit a Variant α equivalent.
- Snapshot #3 (`sphere_2g_p1_aniso`) skipped — Variant α is isotropic.
- Bare-sphere / bare-cylinder MR Variant α helpers at `:535-552`.
- Currently with a 2% drift band per the docstring `:390` — Phase F
  is the targeted closure.

### 5.3 Streaming operator core tests

File: `tests/sn/test_snstreamingoperator.py` — bit-identity vs legacy
matvec, reciprocity round-off, linearity probes, capability checks.

Key tests:
- `test_apply_slab_bit_identical_to_legacy` (`:175-210`)
- `test_apply_spherical_bit_identical_to_legacy` (`:210-260`)
- `test_apply_cylindrical_bit_identical_to_legacy` (`:260-300`)
- `test_apply_2d_cartesian_bit_identical_to_legacy` (`:300-352`)
- `test_solve_slab_bit_identical_to_transport_sweep` (`:352-374`)
- `test_solve_spherical_bit_identical_to_transport_sweep` (`:374-390`)
- `test_solve_cylindrical_bit_identical_to_transport_sweep` (`:390-422`)
- `test_reciprocity_round_off` (`:423-466`)
- `test_apply_is_linear` (`:518-565`)

### 5.4 Phase C gates (sphere/cyl L0 + flat-flux invariants)

File: `tests/sn/test_phase_c_gates.py` — Gate 1.1 (flat flux on
homogeneous reflective sphere), Gate 1.6 (verifies decorators), etc.

### 5.5 SI-vs-Krylov manifestation

File: `tests/sn/test_sweep_operator_inconsistency.py` — pins the
Krylov-vs-sweep deviation on collisionless reflective slab,
demonstrating Manifestation #7's symptomatic form. `:167-188` tests
that Krylov gives exact flat flux while the sweep deviates.

### 5.6 BC tests

Files:
- `tests/sn/test_boundary_conditions.py`
- `tests/sn/test_snmesh_realizer_wiring.py` (10 tests — Wave 8 / C188 wiring)
- `tests/sn/test_method_space.py`
- `tests/sn/test_angular_average_operator.py`

### 5.7 Sweep vs apply consistency

File: `tests/sn/spatial/test_sweep_vs_apply_consistency.py` — the
direct gate for Manifestation #7's resolution.

### 5.8 Phase G Step 2 diagnostics (the prior failed attempt)

Files (5 diagnostic scripts):
- `tests/sn/diagnostics/phase_g_step2_00_baseline.py`
- `tests/sn/diagnostics/phase_g_step2_01_psi_comparison.py`
- `tests/sn/diagnostics/phase_g_step2_02_sncell_residual.py`
- `tests/sn/diagnostics/phase_g_step2_03_closure_audit.py`
- `tests/sn/diagnostics/phase_g_step2_04_fixed_source.py`
- `tests/sn/diagnostics/phase_g_step2_05_homogeneous.py`

These document the prior attempt's evidence chain and live next to
the `numerics-investigator` memo at
`.claude/agent-memory/numerics-investigator/issue_196_phase_g_step2_diagnostic.md`.

---

## 6. Anti-recommendations (the Step-2 agent MUST internalise these)

The Step 2 brief must be written so the executing agent does NOT
do any of the following:

1. **SNMesh IS the SN phase space. Do not create
   `DiscreteOrdinatesPhaseSpace`.** A `(sn_mesh, quadrature,
   n_groups)` wrapper duplicates 2/3 of its fields (`quadrature` is
   `sn_mesh.quad`; `geometry` is `sn_mesh.curvature`). The only new
   scalar — `n_groups` — must be passed as a constructor parameter to
   any operator that needs it (or held on the operator instance
   itself, as `SNStreamingOperator` does today via `sig_t.shape[2]`).
   If the four-operator type **requires** a single bundled "phase
   space" handle for ergonomic reasons, the work is to **add `ng` or
   `n_groups` as an attribute on SNMesh** (one line), NOT to introduce
   a new wrapper. The existing `DiscreteOrdinatesPhaseSpace` in
   `orpheus/sn/streaming.py:120` must be **DELETED** as part of Step
   2's cleanup.

2. **`StreamingOperator(LinearOperator)` already exists at
   `orpheus/sn/streaming.py:250`. EXTEND it; do not re-create it.**
   The class is in the right place architecturally but: (a) its
   `apply` delegates verbatim to the procedural
   `transport_operator_matvec_*` (the very procedural code that
   Manifestation #7 traces to), (b) it has a broken TYPE_CHECKING
   import `from .mesh import SNMesh` (`mesh.py` doesn't exist; SNMesh
   is in `geometry.py`). Step 2's actual mathematical work is to
   replace the `apply` body with primitive `LinearOperator`
   compositions per the
   `.claude/agent-memory/explorer/issue_196_sn_operator_architecture_audit.md`
   plan (WDD as one primitive, M-M as a second, BC trace as a third,
   collision as the fourth, source assembly as the fifth). The
   wrapper class is fine; the wrapper class's BODY is what must be
   rewritten.

3. **`CollisionOperator` already exists at `orpheus/sn/streaming.py:416`.
   EXTEND it; do not re-create it.** It correctly carries `sigma_t` as
   its only material data, advertises `{apply, solve, apply_transpose}`,
   and is self-adjoint. Its packed-vector apply path is currently
   fragile (it imports `build_equation_map_spherical` and calls it
   unconditionally — `streaming.py:499` — which will be wrong for slab
   and 2-D Cartesian). Step 2 must fix the EquationMap dispatch via
   `sn_mesh.curvature`, not by writing a second `CollisionOperator`
   class somewhere else.

4. **`ScatteringOperator` (S) lives in `orpheus/sn/scattering.py:202`
   and is wired through `SNSolver.scattering_op` at
   `solver.py:240-252`. DO NOT create a parallel S.** Use
   `ScatteringOperator.from_solver_data(...)` — that classmethod is
   the canonical construction site. Its constructor reads exactly
   the precomputed structures `SNSolver` already builds; Step 2's
   four-operator algebra should obtain its `S` from there.

5. **`FissionOperator` (F) lives in `orpheus/sn/fission.py:82` and is
   wired through `SNSolver.fission_op` at `solver.py:253-255`. DO NOT
   create a parallel F.** Use `FissionOperator.from_solver_data(chi=,
   sig_p=)`. The :math:`1/k` factor stays at the eigenvalue iteration
   layer, NOT inside F (see `fission.py:35-54` — this is an
   architectural decision documented in the module doc).

6. **`SNStreamingOperator` (legacy L) at `orpheus/sn/operator.py:1114`
   STILL EXISTS and is the active production L for Cartesian and the
   Krylov inner solver across all geometries.** It is currently fused
   `(L+C)` (σ_t inside `transport_operator_matvec_*`). Step 2 must
   decide whether to (a) delete it once the Phase G `StreamingOperator`
   is functional, or (b) reuse its `_ensure_eq_map`,
   `_ensure_dense_matrix`, and `apply_transpose` machinery in the
   Phase G class. Option (b) is materially easier — the EquationMap
   dispatch + dense-probe-transpose are non-trivial — but the
   sweep-frame matvec body (operator.py:412-1107) is exactly the
   procedural code that the four-operator algebra must replace.

7. **`bc_left/right/xmin/xmax/ymin/ymax` are already realised
   `LinearOperator`s. DO NOT re-realise them.** When constructing
   `B` for `L = StreamingOperator(V, boundary=B)`, the imports must
   read:

   ```python
   from orpheus.sn.geometry import SNMesh
   # NOT: from orpheus.sn.mesh — there is no orpheus/sn/mesh.py
   
   B = sn_mesh.bc_right  # for curvilinear; the 1-arg LinearOperator
   ```

   `streaming.py:111` has the wrong import path
   (`from .mesh import SNMesh`) — this must be corrected to
   `from .geometry import SNMesh` as part of Step 2.

8. **`transport_sweep` (`sweep.py:112`) is the production sweep entry
   point. The four-operator algebra's `(L+C).solve` is supposed to
   eventually REPLACE it** (or at least be the curvilinear path's
   sole implementation). Phase G Step 2 already wired the curvilinear
   branch through `_curvilinear_sweep_via_canonical_operator`. Step 2
   of the **replan** must:
   - Verify this wiring is sound (broken `mesh.py` import suggests
     it may currently fail import at test-collection time — check
     before assuming the Phase G code path is exercised by CI).
   - Decide whether to similarly route the Cartesian sweep through
     `(L+C).solve` (and if so, write that primitive — the current
     `solve_within_group_sn_curvilinear` raises if `curv not in
     ("spherical", "cylindrical")`, see `streaming.py:602-606`).

9. **The `cell_update` (`DiamondDifference`) and
   `pole_angular_closure` (`MorelMontryAngularSweep`) strategies live
   ON SNMesh.** They are NOT operator-algebra parameters. Step 2 must
   not move them into L or C or anywhere else; they are properly
   scoped on the SN-augmented mesh as part of the discretisation
   choice. The operator-algebra layer reads them via the SNMesh.

10. **The realised BC's `axis` must match the face being modelled.**
    `_resolve_one` (`geometry.py:393-396`) inspects the face name to
    pick `axis="x"` or `axis="y"` for `ReflectiveBoundary`. Anything
    that constructs a `BoundaryOperator` independently of SNMesh's
    `_resolve_bcs` must replicate this logic — DON'T. Always read
    `sn_mesh.bc_*` for the realised op.

### 6.1 Imports the Step-2 agent MUST use (not reinvent)

```python
# SNMesh and the auxiliary phase-space data
from orpheus.sn.geometry import SNMesh

# The legacy/production L and its companion machinery
from orpheus.sn.operator import (
    SNStreamingOperator,
    build_equation_map,
    build_equation_map_spherical,
    build_equation_map_cylindrical,
    transport_operator_matvec,
    transport_operator_matvec_spherical,
    transport_operator_matvec_cylindrical,
    solution_to_angular_flux,
    solution_to_angular_flux_spherical,
    solution_to_angular_flux_cylindrical,
)

# Scattering (S) and fission (F) — canonical sites
from orpheus.sn.scattering import ScatteringOperator, LegendreMomentScattering
from orpheus.sn.fission import FissionOperator

# Phase G Step 2 partial-ship — to be cleaned / extended (NOT recreated)
from orpheus.sn.streaming import (
    StreamingOperator,
    CollisionOperator,
    # DiscreteOrdinatesPhaseSpace,  # DELETE this — SNMesh IS the phase space
    build_sn_operators,
    solve_within_group_sn_curvilinear,
    build_carlson_context,
)

# Sweep + the unified sweep dispatch (algebra-of-record solve path)
from orpheus.sn.sweep import transport_sweep

# BC realiser (already wired into SNMesh._resolve_bcs)
from orpheus.sn.boundary_realizer import SNBoundaryRealizer
from orpheus.sn.method_space import SNMethodSpace

# Operator-algebra primitives — the building blocks for L's apply body
from orpheus.numerics.operator import (
    LinearOperator,
    LinearOperatorMixin,
    CAP_APPLY, CAP_SOLVE, CAP_APPLY_TRANSPOSE,
    IdentityOperator, ZeroOperator, ScaledOperator,
    PermutationOperator,
    IncomingOrdinateMaskTensor,
    PeriodicWrapOperator,
    OperatorSum, OperatorProduct,
)
```

---

## 7. Structural surprises worth surfacing

1. **`DiscreteOrdinatesPhaseSpace` ALREADY EXISTS in `orpheus/sn/streaming.py`.**
   The previous Phase G attempt did create it. The brief's wording
   "I'm not sure why there is not a DiscreteOrdinatesPhaseSpace" is
   factually inverted — the type exists, and it is exactly the
   redundant wrapper the user objected to. The replan must surface
   this and decide: delete it outright; or repurpose SNMesh as the
   exposed phase-space type, ensuring SNMesh carries an `ng` /
   `n_groups` attribute or the operator constructors take `ng` as a
   parameter.

2. **`orpheus/sn/streaming.py:111` has a broken import path.**
   `TYPE_CHECKING` branch says `from .mesh import SNMesh` —
   `orpheus/sn/mesh.py` does not exist. SNMesh is in
   `orpheus/sn/geometry.py`. The TYPE_CHECKING block being inside an
   `if TYPE_CHECKING:` guard means this only breaks static type
   analysis, not runtime — but it is a smoking gun that the previous
   Phase G work was authored by an agent that did not verify SNMesh's
   actual location.

3. **`SNStreamingOperator` (operator.py) and `StreamingOperator`
   (streaming.py) coexist in production today.** SNSolver wires
   the legacy `SNStreamingOperator` (`solver.py:260` →
   `self.L = SNStreamingOperator(...)`); `_curvilinear_sweep`
   (`sweep.py:285-296`) wires the Phase G `StreamingOperator`
   inside a different code path. They both delegate to the same
   underlying `transport_operator_matvec_*` procedural functions,
   so functionally they behave identically — but the duplication
   itself is a Cardinal Rule 2 violation. The replan must pick ONE
   class name and one class home, then delete the other.

4. **`CollisionOperator.apply` (`streaming.py:493-499`)
   unconditionally calls `build_equation_map_spherical`** —
   regardless of geometry. This is a clear bug in the partial-ship:
   slab and 2-D Cartesian callers will get the wrong EquationMap. No
   tests exist for `CollisionOperator.apply` on slab or 2-D Cartesian
   today (only the curvilinear path goes through
   `solve_within_group_sn_curvilinear` and even there `C` is consumed
   via its `sigma_t` attribute, NOT via its `apply` method — see
   `streaming.py:637-638`). So the bug is dormant; but it must be
   fixed when Step 2 actually consumes `C.apply`.

5. **`solve_within_group_sn_curvilinear` is GMRES, not classical
   Picard.** The function name suggests "sweep" but the body
   (`streaming.py:520-688`) is `scipy.gmres` on
   `L.apply(psi, sigma_t=sigma_t)`. The docstring (`:541-557`)
   explains why: pure Picard with lagged redistribution diverges for
   curvilinear SN (ρ > 1). This is the "principled" Phase G
   alternative — both Krylov and SI paths now converge to the same
   fixed point. The architectural consequence: the four-operator
   algebra's `(L+C).solve` IS GMRES-on-L.apply for curvilinear; it is
   NOT a separate sweep algorithm.

6. **The curvilinear SI path now calls Krylov internally.** Reading
   `sweep.py:202-296`: `_curvilinear_sweep` → `solve_within_group_sn_curvilinear`
   → GMRES. So `inner_solver="source_iteration"` on a sphere/cylinder
   currently calls GMRES under the hood. This is correct (per the
   Phase G design memo) and closes Manifestation #7 by construction —
   but it means **the historical sweep code path** (`_sweep_1d_spherical`,
   `_sweep_1d_cylindrical`) has been REMOVED. The replan must verify
   this is intentional and that the audit memo's call-graph at
   `:60-130` (referencing `_sweep_1d_spherical:397-595`) is
   currently STALE.

7. **The previous Phase G work did not delete the legacy procedural
   matvecs.** `transport_operator_matvec_spherical`,
   `transport_operator_matvec_cylindrical`, and
   `transport_operator_matvec` (Cartesian) all still exist in
   `operator.py`. The `StreamingOperator.apply` method routes through
   them. So the "four-operator algebra" claim — that L's apply is a
   composition of primitives — is **not yet true**; it is still
   delegation to the procedural code. Step 2 of the replan is where
   that delegation must be replaced.

8. **`OctantLabel` / `SweepDependencyGraph` (`sweep_graph.py`) is
   2-D-Cartesian-only.** They are built only in `_setup_cartesian`
   (geometry.py:817-823); for curvilinear they are `None`. Any
   four-operator-algebra design that wants to expose "sweep DAG" as
   a primitive must NOT generalise this across geometries — it is
   SN 2-D-Cartesian-specific by design. (MoC, per the user memory
   `project_moc_structure`, will use fiber bundles + solution
   sheaves, not a shared SweepGraph.)

---

## Memo-output checklist (parent-agent action)

If the user accepts this plan, the parent agent should copy lines
above (from `# Phase G Replan — SN structural backbone` through
section 7 inclusive, replacing the section-1 title with the
explorer-memory frontmatter below) into:

`/Users/rodrigo/git/nuclear/ORPHEUS/.claude/agent-memory/explorer/issue_196_phase_g_replan_sn_structure.md`

with frontmatter:

```markdown
---
name: issue-196-phase-g-replan-sn-structure
description: Structural map of the SN module for Phase G replan of Issue #196. Confirms SNMesh IS the discrete-ordinates phase space, surfaces the half-shipped DiscreteOrdinatesPhaseSpace + StreamingOperator + CollisionOperator from the prior Phase G attempt, and lists the anti-recommendations the Step 2 agent must internalise to avoid re-creating types that already exist.
metadata:
  type: project
---
```

And add to `MEMORY.md`:

```markdown
- [Phase G replan SN structure map](issue_196_phase_g_replan_sn_structure.md) — Structural surface of SNMesh + the four-operator layer for Phase G Step 2 replan
```
