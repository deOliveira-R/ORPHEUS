# R-1 Step 4 Step G — convention crosswalk for `solve_sn_fixed_source`

**Branch**: `refactor/sn-operator-algebra` @ `757695c`
**Date**: 2026-05-22
**Trigger**: session-2 plan §2.2 + `.claude/lessons.md` L17.  Companion
to `.claude/plans/r1_step4_g_dependency_audit.md` (callgraph; what
retires).  This document covers **conventions** (what each value
means at each subsystem boundary).

> **L17 (verbatim):**  "A carve that crosses subsystem boundaries
> (operator algebra ↔ sweep, scalar ↔ per-ordinate, packed ↔ typed)
> MUST start with a written crosswalk table enumerating each
> subsystem's input/internal/output convention.  The crosswalk is the
> architecture; the code is its transcription.  Skipping it costs
> ~3× debug time."

**Scope** — every value that crosses the boundary between any two of:

* the public API (`solve_sn_fixed_source` entry point)
* `_solve_fixed_source_si` (SI inner loop)
* `_solve_fixed_source_krylov` (Krylov inner loop)
* `_build_rhs_*` (RHS packing for Krylov; retires in G3a)
* `transport_sweep` (path-forward sweep primitive)
* `solver.L.apply` (legacy SNStreamingOperator matvec; retires in G3f)
* `KrylovAcceleration` / `SourceIteration` (path-forward iteration primitives)
* `_make_sweep_preconditioner` (packed-vector preconditioner; retires in G3a)
* `solution_to_angular_flux*` (packed→typed decoder; retires in G3d/G3e)

---

## Subsystem inventory (today's tree, `a638798`/`757695c`)

| Subsystem | Role | Lives in | Path-forward fate |
|---|---|---|---|
| `solve_sn_fixed_source` | Public entry; dispatches on `inner_solver` | `solver.py:1127` | Stays; internals swap (G1, G2) |
| `_solve_fixed_source_si` | Cartesian-default SI loop | `solver.py:1278` | Carved onto `SourceIteration` (G2) |
| `_solve_fixed_source_krylov` | Curvilinear-default GMRES loop | `solver.py:1374` | Carved onto `KrylovAcceleration` (G1; 1-D only — 2-D NotImpl) |
| `_build_rhs_{cartesian,spherical,cylindrical}` | Packs in-scatter + (n,2n) + external into packed-vector RHS; applies `/sum_w` at packing site | `solver.py:839 / 938 / 982` | **RETIRES G3a** (dissolves with packed-vector layout) |
| `_make_sweep_preconditioner` | Wraps WDD sweep as scipy `LinearOperator` for GMRES | `solver.py:680` | **RETIRES G3a** |
| `solver.L.apply` (SNStreamingOperator matvec) | Packed `(n_unknowns,)` → packed `(n_unknowns,)` matvec | `operator.py:1493` | **RETIRES G3f** |
| `build_equation_map*` | Build `EquationMap` slot dispatch | `operator.py:218 / 509 / 637` | **RETIRES G3e** (Cart-1-D deferred to Phase A; G3h) |
| `solution_to_angular_flux*` | Packed → typed `(N, ng, nx, ny)` decoder | `operator.py:262 / 535 / 638` | **RETIRES G3e** (Cart-1-D deferred to Phase A; G3h) |
| `transport_sweep` | Typed sweep primitive (`PerOrdinateSource` → `AngularFlux`) | `sweep.py:transport_sweep` | Stays (already path-forward post-Phase-1.1/1.2) |
| `KrylovAcceleration` / `SourceIteration` | Typed iteration primitives | `numerics/iteration.py` | Stays (consumed by G1 / G2) |
| `transport_operator_matvec_unified` (cell body) | Cell-centre matvec; consumes `(N, ng, nx, ny)`; **face state still uses packed slot maps** | `operator.py:904` | Cell body stays; **face-state signature retires G3c** in favour of native `BoundaryFlux` consumption |
| `EquationMap` (class) + `_with_traces` family | Slot-dispatch packed layout | `operator.py:119 + 677..768` | **RETIRES G3d** (whole legacy layout retires per user direction 2026-05-22) |

---

## Convention crosswalk — per-axis tables

The crosswalk has six axes.  Every table reads **left to right** as
the value flows through the call chain.  The **Bridge** rows show
where today's code applies a unit-conversion or layout-conversion
hop; the **Path forward** rows show where the value lands post-Step-G
(the bridge dissolves at the producer per Pattern 7 / L18).

### Axis 1 — `external_source` shape + magnitude

| Site | Today (legacy) | Path forward (post-G) |
|---|---|---|
| **Caller of `solve_sn_fixed_source`** | `np.ndarray`, shape `(N, ng, nx, ny)`, **per-ord density** (R-1 Step 4 A1 — `/sum_w` projected at producer via `PerOrdinateSource.from_isotropic`) | `AngularFlux` (or `PerOrdinateSource`), shape `(N, ng, nx, ny)`, **per-ord density** — same magnitude convention; type tag stronger |
| **`solve_sn_fixed_source`** | Validates shape `== (N, ng, nx, ny)`; passes opaque to `_solve_fixed_source_*` | Validates typed object; passes through |
| **`_solve_fixed_source_si`** | Adds to `iso_per_ord.values + Q_aniso_p1`; wraps in `PerOrdinateSource` for `transport_sweep` | Direct consumer — already typed at API boundary; no wrap step |
| **`_solve_fixed_source_krylov`** | **Packs to legacy layout** via `eq_map.ordinate` / `ix` / `iy` loop (`solver.py:1424-1428`); ravels `(ng, n_eq)` → `(n_unknowns,)` flat | Direct consumer — `KrylovAcceleration.solve(q_ext_typed)` per Step D pattern; **no packing** |
| **Bridge** (today) | The packing loop at `solver.py:1424` is an `AngularFlux → packed-1D` adapter | Dissolves at G1 — `KrylovAcceleration` consumes typed natively |

### Axis 2 — In-scatter + (n,2n) source magnitude

| Site | Today (legacy SI) | Today (legacy Krylov) | Path forward |
|---|---|---|---|
| **Producer (`ScatteringOperator.apply`)** | Returns `AngularFlux`, per-ord density (Pattern 7 producer; post-A1 commit `de8822d`) | Not used in Krylov path — `_build_rhs_*` rebuilds the in-scatter source from scratch | Returns `AngularFlux`, per-ord density — consumed natively by both inner solvers |
| **Krylov `_build_rhs_*`** | n/a | Builds per-equation `qS = sig_s[mid][0].T @ phi_cell / sum_w`; **applies `/sum_w` at packing site** | Dissolves — `KrylovAcceleration(L+C, S, ZeroOp)` reuses producer-side normalised `S.apply` (mirror of Step D's `_solve_krylov` carve) |
| **(n,2n) `q2`** | `solver._add_n2n_source(Q_iso, phi)` adds to `Q_iso`; `PerOrdinateSource.from_isotropic(Q_iso, sn_mesh)` projects (per-ord density at producer factory boundary) | `_build_rhs_*` computes `q2 = 2.0 * (sig2[mid].T @ phi_cell) / sum_w`; **applies `/sum_w` at packing site** | Single-source — `Q_iso` collected as iso scalar, projected via `PerOrdinateSource.from_isotropic` once; `KrylovAcceleration` consumes the typed result |
| **Bridge** (today) | Producer-side for SI (good); consumer-side `/sum_w` AGAIN inside `_build_rhs_*` (bridge duplication — the bug class L18 warned against) | — | Single producer-side normalisation; consumers see per-ord density only |

**Observation**: the Krylov path's `_build_rhs_*` re-applies `/sum_w`
at the consumer side independently of the producer-side fix that A1
landed.  This is invisible today because the legacy `solver.L.apply`
(SNStreamingOperator) consumes the same packed-vector convention.
After G1, the typed `KrylovAcceleration` carve uses producer-side
normalisation only — the consumer-side `/sum_w` in `_build_rhs_*`
dissolves with the function.

### Axis 3 — Operator matvec (`L·ψ`)

| Site | Today (legacy Krylov) | Path forward (post-G1+G3b+G3c) |
|---|---|---|
| **Input ψ shape** | Packed `(n_unknowns,)` ravel-`'F'` with `EquationMap` slot dispatch (cell unknowns + B1'' face-outer + face-inner if slab) | Typed `AngularFlux` — `values: (N, ng, nx, ny)`, `boundary.xmax_face: (N, ng)`, `boundary.xmin_face: (N, ng)` (slab only) |
| **Cell-centre body** | `solver.L.apply` → `transport_operator_matvec_unified(psi_view, ...)` — takes native `(N, ng, nx, ny)` for cells but **packed face sub-arrays** `(n_face_outer, ng)` + `face_outer_ordinate` int slot map | Native `(N, ng, nx, ny)` end-to-end |
| **Face-state input** | `psi_face_outer = psi[face_outer_block_slice]` gathered from packed via `solution_to_angular_flux_with_traces` | `psi.boundary.xmax_face[mask_inflow, :]` where `mask_inflow` derived from `sn_mesh.quad` (e.g. `quad.mu < 0` for outer face) — **no `eq_map` slot lookup** |
| **σ_t·ψ subtraction** | `sigma_packed = self.sigma_t[:, eq_map.ix, eq_map.iy].ravel('F')`; `m_full - sigma_packed * psi` (at cell slots only; face slots carry no collision) | `m_cell - self.sigma_t[None, :, :, :] * psi.values` (native broadcasting) |
| **Output shape** | Packed `(n_unknowns,)` (cell + face slots) | Typed `AngularFlux` (`values` + `boundary` carrying face residual) |
| **Bridge** (today) | Two adapters: `solution_to_angular_flux_with_traces` (packed → typed face slots) at gather; `pack_with_traces` at scatter | Dissolves G3c — the unified matvec gains a native-shape signature; G3d retires the `_with_traces` family entirely |

### Axis 4 — Preconditioner (sweep-as-precond for GMRES)

| Site | Today (legacy Krylov) | Path forward (post-G1) |
|---|---|---|
| **Wrapping** | `solver._make_sweep_preconditioner(eq_map, n_unknowns)` → scipy `LinearOperator` with `matvec=lambda r: decode(transport_sweep(pack(r), ...))`; reads `eq_map` for slot dispatch | `KrylovAcceleration(preconditioner=lambda q: q)` — explicit identity (Step D pattern; per #200 the block-inverse face precond replaces this later) |
| **Input convention** | Packed residual vector `(n_unknowns,)`; the preconditioner internally packs `transport_sweep`'s typed output back to packed via `solution_to_angular_flux*` decoders | Typed `AngularFlux` residual; identity precond means `lambda q: q` returns it unchanged |
| **Sweep dispatch** | The wrapper packs the residual into a `PerOrdinateSource`, calls `transport_sweep`, decodes back to packed | Direct — `transport_sweep` natively consumes typed `PerOrdinateSource` and emits typed `AngularFlux` (path-forward primitive; no decode needed) |
| **Bridge** (today) | Triple adapter: residual packed→typed for the sweep, sweep output typed→packed for GMRES, eq_map slot lookup for both directions | Dissolves G3a (`_make_sweep_preconditioner` deletes) — sweep-as-preconditioner becomes a one-line `preconditioner=lambda q: LC.solve(q)` reusing the InvertibleOperator's `.solve` (already a pure function post-Phase-1.2) when issue #200 ships the block-inverse face precond |

### Axis 5 — `phi` (scalar flux) shape

| Site | Today | Path forward |
|---|---|---|
| **SI loop** | `phi: (ng, nx, ny)`; updated each iterate via `transport_sweep` return | Same shape; same producer (`transport_sweep` returns native shape) — path-forward already |
| **Krylov outer** | `phi: (ng, nx, ny)`; decoded each outer from packed `solution` via `solution_to_angular_flux* → _scalar_flux_from_angular` | `phi = psi_typed.integrate_angular().values` — `(ng, nx, ny)`; native reduction on `AngularFlux` |
| **Solution boundary** | `ScalarFlux(phi, sn_mesh)` — typed at the return | Same — `Solution.scalar_flux` already typed |

### Axis 6 — Boundary face state (B1'' face-aware)

| Site | Today | Path forward |
|---|---|---|
| **Outer face inflow** | `psi_face_outer: (n_face_outer, ng)`; populated by reading packed slot `face_outer_block_slice`; `face_outer_ordinate: (n_face_outer,)` int array maps slot → ordinate index | `boundary.xmax_face: (N, ng)`; inflow ordinates picked via `quad.mu < 0` (slab) / `quad.mu_x < 0` at outer radius (curvilinear) — **no slot map** |
| **Inner face inflow (slab)** | `psi_face_inner: (n_face_inner, ng)`; same shape & dispatch as outer | `boundary.xmin_face: (N, ng)`; inflow ordinates picked via `quad.mu > 0` |
| **Inflow-ordinate mask source** | Precomputed `EquationMap.face_outer_ordinate` / `face_inner_ordinate` int arrays | Derived on-demand from `sn_mesh.quad`'s direction-sign properties |
| **Bridge** (today) | `eq_map` lookup at every gather/scatter; `solution_to_angular_flux_with_traces` adapter at decode | Dissolves G3b (`_ensure_eq_map` removed from typed operators) and G3d (`_with_traces` family retires) |

---

## Bridge sites in today's tree — what dissolves and when

Each bridge below is a load-bearing site for the L17/L18 lesson —
the **producer-side fix** removes the bridge entirely, not just
moves it elsewhere.  Ordered by retirement step:

1. **`solver.py:1424-1428`** — `external_source: (N, ng, nx, ny) → ext_packed: (ng, n_eq)` via `eq_map.ordinate[k] / ix[k] / iy[k]` loop.  Dissolves **G1**: `KrylovAcceleration.solve(q_ext_typed)` consumes typed directly.
2. **`_build_rhs_*` consumer-side `/sum_w`** (`solver.py:908, 925, 930, 966, 972`) — duplicates the producer-side A1 normalisation at the packing site.  Dissolves **G1+G3a**: the typed Krylov path uses producer-side `S.apply` only.
3. **`solver._make_sweep_preconditioner.matvec`** — packs residual, runs sweep, decodes; triple adapter for one operation.  Dissolves **G3a**: `KrylovAcceleration(preconditioner=lambda q: q)` is the production choice today; issue #200 plans a block-inverse face precond that consumes typed natively.
4. **`solver.L.apply` (SNStreamingOperator matvec) packed-1D in/out** — wraps `transport_operator_matvec_unified` with eq_map slot dispatch.  Dissolves **G3f**: `SNStreamingOperator` retires; `KrylovAcceleration` calls the typed `(L+C).apply` directly via the operator-algebra Protocol.
5. **`transport_operator_matvec_unified` packed-face-slot signature** — `psi_face_outer: (n_face_outer, ng)` + `face_outer_ordinate: (n_face_outer,)` int slot map.  Dissolves **G3c**: new native signature `face_outer: (N, ng)` (sliced internally by `quad.mu < 0`).
6. **`StreamingOperator._ensure_eq_map`, `CollisionOperator._ensure_eq_map`, `_eq_map` cache, `n_unknowns`** — eq_map machinery inside the typed operators (sign that the typed operators still wrap legacy compute).  Dissolves **G3b** alongside (5).
7. **`solution_to_angular_flux*` decoders** — packed → typed.  Dissolves **G3e**.
8. **`AngularFlux.to_flat_with_traces` / `from_flat_with_traces`** — typed ↔ legacy bridges; no consumer left to feed.  Dissolves **G3g**.

---

## Path-forward inner-solver shapes (canonical, post-Step G)

### `_solve_fixed_source_krylov` (post-G1+G3a)

```
external_source (AngularFlux)         ──┐
        │  q_ext_typed = external_source │
        ▼                                 │
KrylovAcceleration(                       │  L,C,S,F all consume/emit AngularFlux
    L + C,                                │  natively (post-G3b); transport_sweep
    S,                                    │  is the production sweep primitive;
    ZeroOperator(),                       │  preconditioner=lambda q: q (until #200).
    preconditioner=lambda q: q,           │
)                                         │
        │                                 │
        ▼                                 │
psi_typed: AngularFlux  ──────────────────┘
        │
        ▼
phi = psi_typed.integrate_angular().values     # (ng, nx, ny)
```

No `eq_map`, no `pack/unpack`, no `_build_rhs_*`.  2-D Cartesian
raises `NotImplementedError` from `KrylovAcceleration.solve` (mirroring
`_solve_krylov`; Phase A absorbs).

### `_solve_fixed_source_si` (post-G2)

```
external_source (AngularFlux)
        │
        ▼
SourceIteration(
    L + C,                            # InvertibleOperator (sweep .solve)
    S,                                # ScatteringOperator (typed apply)
    ZeroOperator(),                   # no fission
).solve(q_ext_typed)                 # iterates until ‖φ−φ_prev‖/‖φ‖ < tol
        │
        ▼
psi_typed: AngularFlux
```

Geometry-agnostic — `transport_sweep` (used inside
`InvertibleOperator.solve`) dispatches on `sn_mesh.curvature` /
`reduced` and handles slab + sphere + cylinder + 2-D Cartesian
natively (SURPRISE-5 in the dependency audit).

### Public API stays unchanged

```
solve_sn_fixed_source(
    materials, mesh, quadrature,
    external_source: np.ndarray,         # (N, ng, nx, ny) per-ord density
    ...
) -> Solution
```

The signature accepts `np.ndarray` today; under Step G it can EITHER
(a) keep the bare-array signature and wrap into `AngularFlux`
internally at the entry (preserving public API stability), OR
(b) tighten to `AngularFlux` over a deprecation cycle.  **Decision
deferred to G1 design** — preserve bare-ndarray entry by default
since `solve_sn_fixed_source` is L1-MMS-facing and external callers
exist.

---

## Path-forward convention contract (one-page reference)

The full canonical shapes + magnitudes the path-forward solver
enforces.  Each row IS the API contract for that field.

| Field | Shape | Magnitude / Convention | Notes |
|---|---|---|---|
| `external_source` (API input) | `(N, ng, nx, ny)` | Per-ord density (`/sum_w` projected at producer) | Caller uses `PerOrdinateSource.from_isotropic` for iso scalar sources |
| `external_source` (internal, typed) | `AngularFlux(values, sn_mesh)` | Per-ord density | Wrapped at entry or directly typed |
| `AngularFlux.values` | `(N, ng, nx, ny)` | Per-ord density (`ψ_n` ordinate-by-ordinate) | The path-forward canonical layout |
| `AngularFlux.boundary.xmax_face` | `(N, ng)` | Per-ord face flux at outer boundary | Inflow ordinates: `quad.mu < 0` (slab) or `quad.mu_x < 0` at outer r (curvilinear) |
| `AngularFlux.boundary.xmin_face` | `(N, ng)` | Per-ord face flux at inner boundary (slab only) | Inflow ordinates: `quad.mu > 0` |
| `ScalarFlux.values` (`phi`) | `(ng, nx, ny)` | Iso scalar flux density | `phi = ψ.integrate_angular().values` |
| `PerOrdinateSource.values` | `(N, ng, nx, ny)` | Per-ord density | Same as `AngularFlux.values` by shape; semantic distinction is "source vs flux" (dimensional sin acknowledged in #205) |
| `ScatteringOperator.apply(AngularFlux)` | input `(N, ng, nx, ny)` → output `(N, ng, nx, ny)` | Per-ord density (`/sum_w` at producer; Pattern 7) | Producer-side normalisation post-A1 |
| `FissionOperator.apply(AngularFlux)` | input `(N, ng, nx, ny)` → output `(N, ng, nx, ny)` | Per-ord density | Same shape as scattering |
| `(L+C).apply(AngularFlux)` | input `(N, ng, nx, ny)` → output `(N, ng, nx, ny)` + `boundary` face residual | Per-ord density for cells; B1'' face residual at outer (+ slab-inner) | Native matvec post-G3c |
| `(L+C).solve(AngularFlux, *, initial_guess=None)` | rhs `(N, ng, nx, ny)` → output `(N, ng, nx, ny)` | Per-ord density | `transport_sweep` underneath; pure function post-Phase-1.2 |
| `KrylovAcceleration(L+C, S, ZeroOp).solve(q_ext_typed)` | `(N, ng, nx, ny)` → `(N, ng, nx, ny)` | Per-ord density | GMRES on `(L+C-S)·ψ = q_ext` |
| `SourceIteration(L+C, S, ZeroOp).solve(q_ext_typed)` | same | same | Picard fixed-point |
| `sn_mesh.quad.weights.sum()` | scalar | `W = sum_w` (quadrature normalisation) | Used ONLY at producer sites (`PerOrdinateSource.from_isotropic`, `ScatteringOperator.apply` typed branch); never at consumer sites |
| Inflow-ordinate mask for outer face | `(N,)` bool or `(n_inflow_outer,)` int | Derived from `quad`; e.g. `quad.mu < 0` for slab | NO precomputed `EquationMap.face_outer_ordinate` |
| Inflow-ordinate mask for inner face (slab) | `(N,)` bool or `(n_inflow_inner,)` int | Derived from `quad`; e.g. `quad.mu > 0` | Symmetric to outer |

---

## Bug-class containment — what this crosswalk prevents

| Bug class | How it would manifest today (sans crosswalk) | How the crosswalk prevents it |
|---|---|---|
| ERR-049 (convention drift) re-introduction | Carve threads per-ord rhs through `_build_rhs_*` which re-applies `/sum_w`; `k_eff` lands at `k_inf / (W·c_s)` | Axis 2 makes the duplication explicit; producer-side normalisation is the only `/sum_w` in the chain |
| ERR-050 (silent-precond-fallback) re-introduction | Stateful read inside `L.solve` for the sweep precond on GMRES residual vector | Axis 4 shows preconditioner is identity today; sweep-as-precond is principled via Phase 1.2 pure-function `LC.solve` |
| Packed/typed shape swap | `transport_sweep` consumes packed `(n_unknowns,)` instead of typed `PerOrdinateSource`; silent shape mismatch | Axis 1 + path-forward shapes table fix the typed contract; bare-ndarray entry only at API |
| Face-state slot-map drift | `face_outer_ordinate` slot map encodes a different inflow-ordinate set than `quad.mu < 0` derives | Axis 6 removes the slot map; mask is derived from `quad` at call time |
| 2-D Cartesian Krylov regression on path forward | G1 carve silently breaks 2-D Cartesian (no typed Krylov support) | Path-forward block says `KrylovAcceleration` raises `NotImplementedError` on 2-D; SI handles 2-D via `transport_sweep` (SURPRISE-5) |

---

## Cross-references

* `.claude/lessons.md` L17 (convention crosswalk before carve), L18
  (Pattern 7 at producer), L19 (None-default stateful coupling), L20
  (retirement dependency audit), L21 (sweep & matvec share one
  strategy).
* `.claude/skills/coding-elegance/SKILL.md` §"Convention crosswalk
  template" + Pattern 7.
* `.claude/skills/vv-principles/error_catalog.md` ERR-049
  (convention drift, closed by Phase 1.1) + ERR-050 (silent precond
  fallback, closed by Phase 1.2).
* `.claude/plans/r1_step4_g_dependency_audit.md` — companion (the
  callgraph; what retires).  CORRECTION section is the source of
  truth for the legacy vs path-forward classification used here.
* `.claude/plans/r1_step4_session2_followup.md` §2.2 (this artifact)
  + §3-§4 (G1 / G2 / G3 implementation phases).
* GitHub issue #199 (Step G omnibus); #200 (block-inverse face
  preconditioner — replaces today's identity precond); #174
  (`_build_rhs_cartesian` softer refactor — subsumed by G3a); #205
  (cross-method field architecture — picks up the `AngularSource` vs
  `AngularFlux` dimensional split deferred from R-1).
