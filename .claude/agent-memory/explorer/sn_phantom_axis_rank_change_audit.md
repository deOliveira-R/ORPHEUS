---
name: sn-phantom-axis-rank-change-audit
description: Operator-machinery impact of dropping the phantom ny=1 axis (1-D → genuine rank-1 (N,ng,nx)); the single-lever _shape_for_mesh + the local/structured operator split + the latent to_flat vs n_unknowns_flat mismatch.
metadata:
  type: project
---

The SN "drop the phantom `ny=1`" carve (worktree `sn-nd-layout`, branch
`worktree-sn-nd-layout`): make every array genuinely `(N, ng, *spatial)` with
spatial rank == ndim. 1-D today is `(N, ng, nx, 1)`; target `(N, ng, nx)`.

**Why:** architecture-forward dim-agnostic SN ([[feedback-architecture-forward-not-legacy-fit]]).
**How to apply:** the SHAPE below is durable; re-confirm `file:line` via Nexus at pickup.

## The single lever (headline)
`orpheus/transport/fields/_bases.py` — `BulkField._shape_for_mesh` returns
`(mesh.quad.N, mesh.ng, mesh.nx, mesh.ny)`; `ScalarField._shape_for_mesh` returns
`(mesh.ng, mesh.nx, mesh.ny)`. These TWO methods are the single source of truth
for EVERY typed bulk field's shape (AngularFlux/ScalarFlux/source-sinks/residuals/
displacements). Change them to `(N, ng, *mesh.spatial_shape)` / `(ng, *mesh.spatial_shape)`
and every typed field becomes rank-d at once. `Field.__post_init__` enforces
`values.shape == space.shape` STRICTLY → every stale rank-2 producer fails LOUD
(a fail-loud carve, not silent). `SNMesh.spatial_shape` ALREADY returns genuine
`(nx,)` for 1-D; `ndim`/`axes` are rank-d-ready. `HarmonicMomentField.from_mesh_and_L`
needs the matching change (its `(mesh.ng, mesh.nx, mesh.ny)`).

## Local vs structured split (CONFIRMED empirically)
- **Spatially-LOCAL (C, S, F, sources):** math is rank-agnostic; they hardcode
  rank-2 only in POSITIONAL INDEXING, not the contraction.
  - `CollisionOperator` (operator.py): `self.sigma[None, :, :, :] * psi.bulk.values`
    — the `[None,:,:,:]` (4 indices) hardcodes rank-2. Fix → `self.sigma[None]`
    (single prepend, ellipsis-broadcast). Same in `.solve` (`q.bulk.values / sigma[None,:,:,:]`).
  - `ScatteringOperator`: einsums in `MomentProjection`/`HarmonicMomentReconstruction`
    (`projection.py`) ARE ellipsis (`"n,nlm,n...->lm..."`) → survive untouched. BUT
    `AngularFlux.integrate_angular` uses `"n,ngij->gij"` (angular_flux.py:131) — hardcoded
    rank-2; fix → `"n,ng...->g..."`. AND `material_xs_field` per-material verbs
    (`apply_p0_in_scatter`/`apply_n2n`/`apply_legendre_scattering_moments`) index
    `Q[:, ix, iy]` with the 2-D `cells_by_material` `(ix,iy)` pair — the einsum after
    is rank-agnostic (`[:, ix, iy]` flattens spatial to one cell list `c`), but the
    gather/scatter is rank-2-positional. Rank-d form: `cells_by_material` yields
    `np.nonzero(mat_map)` of length ndim, or flat-cell indices into a spatial-flattened view.
  - `FissionOperator`: `RankOneOperator(chi, sig_p, axis=0) & IdentityOperator()` —
    docstring SAYS `(N_axis, *spatial)`, reduces axis=0, broadcasts `*spatial`. GENUINELY
    rank-agnostic → survives UNTOUCHED (only needs rank-d chi/sig_p from material_xs_field).
- **Spatially-STRUCTURED (L / sweep):** rank matters; this is where `(...,nx,1)` bites.
  - 1-D matvec `_MSpatialOperatorSum._compute_LpC` / `_compute_decomposition` /
    `_compute_LpC_transpose` (operator.py): phantom axis baked HARD — `np.zeros((ng,N,nx,ny))`
    with ny=1, `psi_view[:,:,0,0]` (pole seed), `[:, global_dir, i, 0]` writes,
    `sigma_t[:,:,0]`, `volumes[:,0]`, `.transpose(1,0,2,3)`, `[..., 0]`.
  - 1-D sweep `_sweep_1d_unified` (sweep.py): the CLEAN pattern — drops y on ENTRY
    (`Q[:,:,:,0]`, `sig_t[:,:,0]`, `volumes[:,0]`), works in genuine rank-1 internally,
    RE-ADDS ny=1 at RETURN (`scalar_flux[:,:,None]`, `angular_flux (N,ng,nx,1)`). The
    pole closure (`pole_angular_closure.py`) ALSO drops y on entry (`[..., 0]`) and works
    rank-1 internally (`_dAw_per_level (nx,N)`, `_V = volumes[:,0]`). So the curvilinear
    pole machinery is ALREADY rank-1 inside; only the operator.py orchestrator carries phantom.
  - 2-D Cartesian (`_apply_2d_cartesian`, `_sweep_2d_wavefront`, DD `cell_kernel_batch`/
    `update_batch`, sweep_graph) STAYS `(N,ng,nx,ny)` — NOT touched by the 1-D rank change.

## The latent to_flat ↔ n_unknowns_flat mismatch (a SEPARATE bug, NOT the phantom axis)
`TimedFullField.to_flat`/`from_flat` is RANK-AGNOSTIC — `bulk.values.ravel()` /
`reshape(template.bulk.values.shape)`. Round-trips bit-identically against its own
template at ANY rank. The production Krylov path uses ONLY to_flat/from_flat.
BUT `axis.n_unknowns_flat` (the theoretical pack-size) counts only OUTFLOW ordinates
per face; `BoundaryFlux` stores ALL N ordinates per face. So `to_flat.size != n_unknowns_flat`
for ALL geometries (slab 28 vs 24, sphere 24 vs 22, cyl 144 vs 132) — the gap is the
boundary block (inflow ords stored but not counted), NOT the phantom bulk axis (the bulk
`N·ng·nx·1 == N·ng·nx` coincides). `n_unknowns_flat` has ZERO production callers (only the
`SNMesh.n_unknowns_flat` property + one primitive test). NOT load-bearing for the rank change.
The 1-D face layout (`geometry.py` `face_layout`) is ALREADY rank-correct `(N, ng)` — faces
carry NO phantom axis. Boundary handling needs NO change.

## Minimal change-set partition
MUST change (rank-2-positional): `_bases.py` 2× `_shape_for_mesh` (LEVER) + `HarmonicMomentField`
shape + `CollisionOperator` `sigma[None,:,:,:]` (apply+solve) + `AngularFlux.integrate_angular`
einsum + `material_xs_field` `_ensure_cell_views` reshape + `cells_by_material`/4 per-material
`[:, ix, iy]` verbs + the WHOLE 1-D matvec triple in operator.py + `_sweep_1d_unified` entry/exit
+ `geometry.__init__` `(nx,ny)` normalization (`volumes.reshape(N,1)`, `mat_map.reshape(N,1)`,
`ny=1`, `is_1d`, `volumes/areas` shape) + solver.py initial-guess/validation/rate-reduction sites.
SURVIVES UNTOUCHED (genuinely rank-agnostic): `FissionOperator`+`RankOneOperator`, all
`MomentProjection`/`Reconstruction` einsums, `ordinate_scan` (scan.py, leading-axis+trailing-lanes),
iteration drivers (`SourceIteration`/`KrylovAcceleration`/`KEigenvalue` — shape-agnostic via
ravellable protocol), `to_flat`/`from_flat`, the 1-D face layout, the 2-D path.

## Surprise/blocker
NO load-bearing `ny=1` broadcast found — the phantom axis is never RELIED upon for a
broadcast that rank-1 would break; it's purely vestigial positional indexing. The one
asymmetry: C/S/F use `[None,:,:,:]`/`[:, ix, iy]` POSITIONAL indexing that must become
`[None]`/ndim-generic — trivial but MANY sites. Highest-risk surface = the 1-D curvilinear
matvec (operator.py, open #206/#195/#196), but the pole closure underneath is already rank-1,
so the risk is contained to the orchestrator's index bookkeeping.
