---
name: issue-206-phase-b-frame-blast-radius
description: #206 Phase B — extracting a shared 1-D-scan SWEEP frame; the _OctantWalk template anatomy, the 1-D free-function inventory (_run_1d_sweep/_sweep_1d_unified/_ensure_*_cache), zero external callers, _SweepEmit modes, the two-stratum cache, and the coupling hazards (Carlson seed / pole_angular_closure / initial_guess / CumprodScan-vs-ScanMarch-1D divergence).
metadata:
  type: project
---

# #206 Phase B — shared 1-D-scan SWEEP frame blast radius

VERIFIED on branch `refactor/sn-cellupdate-seam-slab` (HEAD `f61e1b0`) by Read/Grep,
2026-06-14. ALL code now lives in `orpheus/sn/loss_representation.py` (the historical
`sweep.py` dissolved at S6.4(f); the 1-D body relocated here with the walks). Nexus graph
was stale (`docs/bailey-bib-migration`); these line numbers are file-truth.

## Headline finding (the relocation is CLEAN)

`_run_1d_sweep`, `_sweep_1d_unified`, `_ensure_geom_cache`, `_ensure_coll_cache` have
**ZERO external callers** — entirely internal to `loss_representation.py`. The only two
production call sites are `CumprodScan.sweep` (`:726`) and `ScanMarch.sweep` 1-D branch
(`:1246`), and they invoke `_sweep_1d_unified` IDENTICALLY (same arg order, same kwargs).
Tests / derivations mention the symbols only in **docstrings/comments** (no live imports):
`tests/sn/primitives/test_boundary_face_layout.py:18`, `tests/sn/sweep/core/test_sweep_graph_nd_admission.py:20`,
`tests/sn/operators/test_streaming_operator.py:1095`, `derivations/diagnostics/diag_krylov_iter_breakdown.py:10`,
`derivations/diagnostics/diag_si_cyl_20cell_nan_step5_root_cause.py:136`. So Phase B is a
pure intra-module relocation → bit-identical; no import-surface churn, no test rewiring.

## The `_OctantWalk` template anatomy (the 2-D pattern to mirror)

Block: `loss_representation.py:444-682`. `@dataclass(frozen=True)`, sole field `mesh: SNMesh`.
NO `__init__` (dataclass-generated). Instantiated as `_OctantWalk(sn_mesh)` at the call sites
(`:925`, `:1124`, `:1368`, `:2314`).

- **The shared seam** = `_interior_walk` (`:481-524`): the ONE octant frame. Loops
  `for sweep in sweeps`; reads `sweep.label.signs` (the in-plane projection, already done by
  the schedule's `_octant_sweep` — NO re-truncation here, by design); `if not any(signs)` →
  `pure_z(oct_idx)` (the degenerate volumetric branch); else maps grazing-0→+1 (`signs_eff`),
  reads inflow via `inflow_of(face)[oct_idx]` over `_inflow_faces(signs_eff)`, runs
  `interior(oct_idx, signs_eff, inflow)`, sheds each capture into `_outflow_faces(signs_eff)`.
- **Two public directions, both route through `_interior_walk`:**
  - `sweep_group(group, *, operands: _SolveOperands, emit: _SweepEmit, boundary_flux, interior)`
    (`:526-580`) — SOLVE. Closes over `pure_z` (`emit.pure_z(oct_idx, Q/sig_t)`), `run_interior`
    (calls the supplied `interior(operands, emit, oct_idx, signs_eff, inflow)`), `shed`
    (writes `boundary_flux.face_view(face)[oct_idx]`). inflow_of = `boundary_flux.face_view`.
  - `loss_action(operator, psi: TimedFullField, interior)` (`:582-682`) — APPLY (matvec).
    Builds `_ApplyOperands`, allocates `LpC`, reads inflow from `psi.boundary.face_view`,
    runs ONE Jacobi group, then assembles the O.4b boundary residual + the typed
    `TimedFullField`.
- **HOW it forks** (the anti-bool-flag design): the direction is carried by the injected
  OBJECTS, never a flag —
  1. the **cell-kernel** = the `interior` callable each representation supplies
     (window's `walk_windowed`, scan-march's row-march, oracle's `walk_full`), in its solve
     (`cell_kernel_batch`) or apply (`residual_kernel_batch`) flavor; and
  2. the **emit policy** = `_SweepEmit` (solve) / the apply accumulator + O.4b residual.
- **The tripwire forbidding a bool `is_solve`**: docstring `:466-468` —
  `tests/sn/operators/test_one_octant_walk.py` enforces the kernel/emit-OBJECT shape.

## The 1-D sweep free functions (the relocation targets)

- `_sweep_1d_unified(Q, sig_t, sn_mesh, boundary_flux, *, initial_guess=None)` (`:1719-1779`).
  Thin: `geom = _ensure_geom_cache(sn_mesh)`; `coll = _ensure_coll_cache(sn_mesh, sig_t, geom)`;
  `return _run_1d_sweep(Q, sig_t, sn_mesh, boundary_flux, geom, coll, initial_guess=...)`.
- `_ensure_geom_cache(sn_mesh) -> GeometryCoefficients` (`:1782-1788`): lazy getattr/build/stash
  on `sn_mesh._geom_cache` (the solver normally pre-stashes; ad-hoc callers build lazily).
- `_ensure_coll_cache(sn_mesh, sig_t, geom) -> CollisionCache` (`:1791-1814`): same lazy
  getattr/build/stash on `sn_mesh._coll_cache`; builds `CollisionCache.from_geometry(geom, sig_t, cell_update)`.
- `_run_1d_sweep(Q, sig_t, sn_mesh, boundary_flux, geom, coll, *, initial_guess=None)` (`:1817-2206`).
  The body. Splits at `is_slab = coord is CARTESIAN`:
  - **SLAB joint-batch** (`:1961-2033`): partition ordinates by `mu>=0`; ONE `ordinate_scan`
    per chain direction (2 scans/sweep regardless of N/ng); no M-M thread; reads/writes
    `boundary_flux.face_view("xmin"/"xmax")`; `initial_guess` UNUSED on slab.
  - **CURVILINEAR per-ordinate** (`:2040-2197`): reads `closure = sn_mesh.pole_angular_closure`;
    extracts `_initial_guess_values(initial_guess)` (`:2053`); per-µ-level Carlson seed via
    `closure.psi_half_seed(psi_level, CarlsonSweepContext(...))` (`:2080-2088`); coupled-pole
    mirror capture/consume (ERR-058, `:2066-2119`); degenerate cyl-axis ordinates fall to the
    slow `cell_update.update` per-cell path (`:2122-2146`); reads single `xmax` face (the pole
    is a regularity condition, no BC face).
- `_initial_guess_values(initial_guess) -> ndarray|None` (`:1675-1716`): container-agnostic
  extractor (legacy `AngularFlux.values` vs composite `TimedFullField.bulk.values` via duck-typed
  `.bulk`). Used ONLY by the curvilinear branch. ALSO a relocation companion.

## `_SweepEmit` (the solve-direction emit policy, `:386-441`)

`@dataclass(frozen=True)`. Fields: `weights (N,)`, `angular_flux|None`, `scalar_flux|None`,
`moment_buf|None`, `Y|None`. `__post_init__` guard: EXACTLY one of {angular pair} XOR {moment
pair} wired (illegal half-state unrepresentable). Modes: **angular** (per-octant write +
`Σ w_n ψ_n` accumulate) and **moment** (`φ_ℓ^m += Σ w_n Y ψ_n`; 2-D-Cartesian peak-mem win).
`pure_z(oct_idx, psi_avg)` emits the no-face volumetric balance using `buf[...] +=`. NOTE for
Phase B: the 1-D body does NOT use `_SweepEmit` today — it writes `angular_flux`/`scalar_flux`
ndarrays directly and REFUSES moment mode (`moment_projection is not None` → ValueError at
`CumprodScan.sweep:716` and `ScanMarch.sweep:1239`). A 1-D frame CAN adopt `_SweepEmit` (angular
mode only) for symmetry, but the angular write granularity differs (per-chain joint-batch slab /
per-ordinate curvilinear, not per-octant), so it would need an emit shape that accepts a chain/
ordinate-granular write — OR keep direct ndarray writes (simpler, bit-identical).

## The cache layer (`orpheus/sn/spatial/sweep_cache.py`, the frame must host)

Two frozen `slots=True` dataclasses, separated by MUTATION CADENCE (cross-domain-attacker Smell #16):
- `GeometryCoefficients` (`:117-334`): Stratum 1, geometry × quadrature ONLY, NO ng axis.
  Built once per SNMesh via `from_mesh_and_quad` (iterates `sn_mesh.dag_walk`, slow Python, ONCE).
  Survives σ_t rebinds / BC changes / every iteration. Fields: chain_idx/chain_idx_inv (N,nx),
  abs_mu/c_in/c_out/tau_inv/mm_a_in_coeff/mu_start/is_degenerate (N,), A_down/A_total/dA_w/V
  (N,nx), level_ordinates (curvilinear-only tuple|None).
- `CollisionCache` (`:342-474`): Stratum 2, geometry × σ_t. `from_geometry(geom, sig_t, cell_update)`
  delegates the (a, 1/denom) math to `cell_update.affine_scan_coefficients` (Issue #236 §2 —
  scheme owns math, cache owns storage+lifetime), adds `cumprod_a` (prefix product along cell
  axis 2). Fields: inverse_denom/a_attenuation/cumprod_a (N,ng,nx). Class-var `_build_count`
  instruments the cache-invariance gate (`test_sweep_cache.py`).
- **Mutation cadence / access**: `SNSolver.__init__` (`solver.py:902-914`) builds BOTH when
  `sn_mesh.reduced is not None` (1-D) and STASHES on `sn_mesh._geom_cache`/`_coll_cache`.
  `rebind_cross_sections` (`solver.py:958-963`) rebuilds ONLY `coll_cache` (geom survives).
  The sweep accesses via the lazy `_ensure_*_cache` getattr on the mesh. A Phase-B frame holding
  the cache must preserve this mesh-stash contract (the solver pre-stash + the lazy fallback)
  so it stays bit-identical and the cache-invariance counter stays at 1.

## Coupling hazards for the relocation

1. **pole_angular_closure / Carlson seed**: curvilinear-only. `_run_1d_sweep:2049` reads
   `sn_mesh.pole_angular_closure`; `:2080-2088` builds per-level `CarlsonSweepContext` and calls
   `closure.psi_half_seed`. This is Pattern-2 single-source-of-truth with the matvec
   (`MorelMontryAngularSweep.precompute_psi_state`). The frame must thread `sn_mesh` (it already
   holds it as the sole field — clean) and keep the closure read in the curvilinear arm.
2. **boundary_flux threading**: slab reads/writes BOTH `xmin`+`xmax` face views (`:1916-1919`,
   `:2031-2033`); curvilinear reads/writes ONLY `xmax` (`:1931`, `:2197`). The inflow-read /
   outflow-write touch DISJOINT ordinate slots (bare-sweep, no internal `bc.apply`) — preserve.
3. **initial_guess threading**: slab IGNORES it (no M-M thread). Curvilinear consumes it via
   `_initial_guess_values` for the Carlson seed + pole-cell upstream. The frame's `sweep` entry
   must keep `initial_guess` kwarg passthrough.
4. **CumprodScan vs ScanMarch-1D divergence — NONE in the 1-D body.** Both `CumprodScan.sweep:726`
   and `ScanMarch.sweep:1246` call `_sweep_1d_unified(Q, sig_t, self.mesh, boundary_flux,
   initial_guess=...)` IDENTICALLY. Both raise the SAME moment-mode ValueError before the call.
   The ONLY structural difference is the surrounding `is_1d` guard in ScanMarch (`:1234`) — its
   2-D arm goes to `_sweep_jacobi`; CumprodScan is 1-D-only (`supports` = `mesh.is_1d`). So a
   shared 1-D frame is a drop-in for BOTH `.sweep` bodies' 1-D path with no behavioral fork.

## What the shared 1-D frame must host + how it forks

HOST: (a) the two-stratum cache ensure/stash (`_ensure_geom_cache`/`_ensure_coll_cache` +
the `sn_mesh._{geom,coll}_cache` mesh-stash contract); (b) the slab joint-batch scan body +
the curvilinear per-ordinate M-M/Carlson body (`_run_1d_sweep`); (c) `_initial_guess_values`;
(d) the `boundary_flux` face I/O (slab xmin/xmax, curvilinear xmax-only); (e) the moment-mode
refusal. FORKS: at the geometry branch `is_slab` (slab joint-batch vs curvilinear per-ordinate),
NOT at solve/apply — Phase B wires only SOLVE. Unlike `_OctantWalk` (forks on injected
kernel/emit OBJECTS), the 1-D scan's variation is GEOMETRY (slab vs curvilinear), already an
internal `is_slab` branch, and the future matvec fork (Phase C) is at the cell-kernel
solve(`ordinate_scan`)/apply(`_compute_LpC`) — the 1-D matvec is currently
`operator.M_spatial._compute_LpC` (CumprodScan.loss_action:741, ScanMarch.loss_action:1367),
a DIFFERENT object than the sweep's `ordinate_scan`, so the Phase-C fork point is the
per-ordinate cell recurrence, mirroring `_OctantWalk`'s kernel injection.
