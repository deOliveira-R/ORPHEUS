---
name: sn-space-angle-discretization-coupling
description: ORPHEUS SN couples spatial+angular discretization ONLY in curvilinear (M-M angular thread injects into the spatial scan); Cartesian is FULLY separable. Two SEPARATE strategy registries exist (spatial CellUpdate / angular PoleAngularClosure) but the curvilinear cell-balance entangles them through the shared denom. Spatial scheme axis is a seam with ONE occupant (DiamondDifference); no Step/LD/EC exist.
metadata:
  type: project
---

## VERDICT (2026-06-13): the claim is CONFIRMED-WITH-NUANCE, geometry-split

User claim: "spatial and angular discretization are currently the same
[coupled] but don't have to be; can they be independently selected/combined?"

**Cartesian (slab/2-D/3-D): ALREADY INDEPENDENT.** No coupling. Spatial DD
and angular GL quadrature are orthogonal. EMPIRICAL: slab-iso MMS O(h²) rate
is bit-identical across N=4/8/16 (4.01/4.00/4.00); space×angle cross-term
|M|/max ≈ 0.000–0.005. The Cartesian wavefront path (`sweep_graph._CellSolve`
→ `DiamondDifference.cell_kernel_batch`) consumes ONLY `s_axes` (per-axis
streaming) + `sigt_cells` — NO τ, NO angular redistribution. Each ordinate's
spatial scan is independent (slab joint-batch = 2 scan calls per sweep
regardless of N, `loss_representation.py:1949-2024`).

**Curvilinear (sphere/cyl): GENUINELY COUPLED — shared CODE, not co-location.**
The M-M angular thread `ψ_{n+1/2}=(ψ̄−(1−τ)ψ_{n−1/2})/τ` couples ordinates
sequentially within a μ-level. The coupling is bidirectional at TWO sites in
`_run_1d_sweep` (`loss_representation.py`):
- L2142-2150: the angular upstream `ang_contrib = (ΔA/w)·c_in·ψ_{n−1/2}` is
  ADDED into the spatial scan's source `b` → angle feeds space.
- L2167-2176: the spatial cell-average `psi_avg_chain` is consumed to emit the
  next ordinate's angular state → space feeds angle.
The shared kernel is `cell_balance_terms` / `cell_balance_for_streaming`
(`orpheus/sn/spatial/cell_balance.py`): the SAME denom carries
`2|μ|·A_down` (spatial WDD) + `(ΔA/w)·c_out` (angular, c_out=α_out/τ) +
`Σ_t·V`. τ (the M-M angular weight) enters the SPATIAL denom that produces
ψ̄ that the SPATIAL closure uses. This is shared algebra, not adjacency.

## ARCHITECTURE: TWO separate strategy registries (seams already half-built)

1. SPATIAL scheme axis: `CellUpdateBase(RegistryMixin, key=...)` Protocol
   (`orpheus/sn/spatial/cell_update.py`). `SNMesh(cell_update=...)` injects it
   (`geometry.py:271`, default `DiamondDifference()`). The cell-update separates
   "spatial closure" (WDD ψ_out=2ψ̄−ψ_in) from "angular closure"
   (M-M) INSIDE its `update()` body via two structural-presence `if`s
   (`face_area_downstream>0`, `angular_upstream is not None`) — NOT geometry
   dispatch. Docstrings name the intended axis: "DD/LD/EC/Step supply the
   kernel pair". **BUT ONLY `DiamondDifference` is implemented** — Step/LD/EC
   exist only as docstring examples + `NotImplementedError` defaults on
   `cell_kernel_batch`/`residual_kernel_batch`. The spatial axis is a
   SEAM-WITH-ONE-OCCUPANT (#158 framing).
2. ANGULAR scheme axis: `PoleAngularClosureBase(RegistryMixin, key=...)`
   Protocol (`orpheus/sn/spatial/pole_angular_closure.py`).
   `SNMesh(pole_angular_closure=...)` injects it (`geometry.py:460`).
   `default_angular_closure_class`: CARTESIAN→`IdentityAngularClosure` (returns
   (0,0); no redistribution), SPHERICAL/CYLINDRICAL→`MorelMontryAngularSweep`.
   This axis has 2 REAL occupants (Identity + M-M) + the swappable
   `psi_half_seed` sub-strategy (`AngularEdgeExtrapolation` default,
   `CarlsonInwardSweep`, `ZeroSeed`). The angular axis is GENUINELY selectable.

So the (space ⊗ angle) PRODUCT exists structurally: two independent injection
points on SNMesh. The Cartesian corner is fully realized (Identity angular ×
DD spatial = no interaction). The curvilinear corner couples the two chosen
strategies through the shared cell-balance denom.

## INDEPENDENT SELECTION difficulty
- Cartesian: EASY (already independent; Identity angular closure already
  decouples by construction).
- New ANGULAR scheme on existing spatial DD: MEDIUM — the `PoleAngularClosure`
  Protocol + registry are built; add a `key="..."` subclass implementing
  `precompute_psi_state`/`cell_contribution`/`angular_adjoint`. Both sweep
  (L2040 `sn_mesh.pole_angular_closure`) and matvec consume it.
- New SPATIAL scheme (LD/EC/Step) on existing M-M angular: MEDIUM-HARD — the
  `CellUpdate` Protocol seam exists, but the FAST curvilinear path
  (`_run_1d_sweep`) inlines DD's `ordinate_scan` (Blelloch closed form) +
  the `0.5*(in+out)` DD closure HARDCODED (L2004, L2167); it does NOT route
  through `cell_update.update` except for the degenerate cyl-axis ordinate
  (L2113). A non-DD spatial scheme would need its own scan kernel + the
  Cartesian wavefront `cell_kernel_batch` override. The Protocol promises
  selectability the fast path doesn't yet honor for curvilinear.

## INTERFERENCE (empirical, the #233 spatial vs #229 angular floors)
Do the two floors ADD or INTERACT? Probes (sphere/cyl/slab aniso MMS,
n_cells × n_ordinates grid):
- SLAB: cross-term ≈ 0 → fully separable (control).
- CYLINDER: error ≈ dominated by azimuthal angular floor (#229); fixed n_phi →
  mesh refinement SATURATES (h-ratios ≈ 1.00); halving n_phi halves error.
  Cross-term small (0.02–0.08) because one term swamps. `E ≈ E_angle(n_phi)`.
- SPHERE: the floor each mesh-refinement saturates at is GATED by N. N=8 →
  spatial plateaus ~4.7e-3 (even rises nx80→160); N=32 → keeps converging to
  3.6e-4. Cross-term LARGE (|M|/max=0.37 at N8→16). This is MIN/GATING
  coupling `E(h,N)≈max(E_space(h),E_angle(N))`-like, NOT additive. They
  INTERFERE: you cannot harvest fine-h benefit without also refining N.
CONCLUSION on "harvest the benefits of both": in Cartesian YES (orthogonal).
In curvilinear the angular floor GATES the spatial accuracy — refining only
one axis stalls; both must advance together. The #233 pole-cell O(h) and #229
azimuthal floor are NOT independent contributors; the angular thread's
interpolation error caps what the spatial scheme can resolve.

Probes: `$CLAUDE_JOB_DIR(84fd66f8)/tmp/diag_sep_{space_angle,cyl,slab,slab_iso}.py`.
Relates: #158 (spatial-scheme axis), #6 (LD as ANGULAR — note the #158/#6
naming conflates LD-spatial vs LD-angular), #229 (azimuthal floor), #233
(pole-cell O(h)), ERR-026 family (curvilinear closure).
