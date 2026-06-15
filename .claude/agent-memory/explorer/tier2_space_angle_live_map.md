---
name: tier2-space-angle-live-map
description: Post-#206 LIVE map of the SN (space⊗angle) discretization machinery for Tier 2 (#158 spatial / #235 angular); registries, dispatch seams, #219 collision, verification hooks, smallest-first-step per tier.
metadata:
  type: project
---

# Tier 2 (space ⊗ angle) LIVE map — verified against main@cba6d2f (graph) / working tree post-#206

Built 2026-06-14 (VERY THOROUGH explorer pass). The `.claude/plans/sn_space_angle_discretization_plan.md`
§2/Tier line refs (`_run_1d_sweep:2004/2113/2167`, `_compute_LpC`) are STALE —
#206 relocated everything onto `_OneDimScanWalk`. This map is the post-#206 truth.

**Why:** ground the Tier 2a (#158 non-DD spatial CellUpdate) + Tier 2b (#235 cyl
2-D (η,φ) angular closure) implementation plans against the actual current tree.
**How to apply:** read this before planning either tier; the SHAPE is stable, the
exact line numbers below drift — re-confirm with Nexus `context`/`query` at pickup.

## 1. SPATIAL registry — `orpheus/sn/spatial/cell_update.py`

- **Protocol** `CellUpdate` (`:326`, runtime_checkable) + **concrete ABC**
  `CellUpdateBase(RegistryMixin, ABC)` (`:509`). Self-registration via `key=`.
- **Three class-level traits:** `is_linear`, `is_positivity_preserving`,
  `is_affine_scannable` (`:374-376` Protocol; `:551-559` ABC; default
  `is_affine_scannable=False` is the OPT-OUT gate — the seam discriminator).
- **Methods a subclass declares:**
  - `@abstractmethod update(visit, total_xs, source, upstream_state) -> CellResult` (`:565`) — the SOLVE direction. ALWAYS required.
  - `@abstractmethod residual(cell_avg, visit, total_xs, source, upstream_state) -> np.ndarray` (`:575`) — the APPLY/matvec twin. ALWAYS required (lockstep with `update`).
  - `cell_kernel_batch(...)` (`:594`) — batched SOLVE for the DAG wavefront walks (Cartesian). Default RAISES NotImplementedError. Only needed to join `FullFieldWavefront`/`MovingFrontierWindow`.
  - `residual_kernel_batch(...)` (`:628`) — batched APPLY twin. Default RAISES.
  - **Scan-family triple (gated on `is_affine_scannable`):** `affine_scan_coefficients(...)` (`:672`, returns `(a_attenuation, inverse_denom)`), `cell_average_from_faces(face_in, face_out)` (`:710`, DD: ½(in+out)), `outgoing_face_from_average(cell_avg, face_in)` (`:730`, DD: 2ψ̄−in). ALL default RAISE NotImplementedError.
- **Occupants:** EXACTLY ONE — `DiamondDifference(CellUpdateBase, key="diamond_difference")`
  at `orpheus/sn/spatial/diamond.py:95`. `is_affine_scannable=True` (`:127`),
  implements all 7 methods (`update:138`, `residual:189`, `cell_kernel_batch:283`,
  `residual_kernel_batch:341`, `affine_scan_coefficients:379`,
  `cell_average_from_faces:467`, `outgoing_face_from_average:482`).
- **Step / LD / ExponentialCharacteristic: DO NOT EXIST as code.** Confirmed by
  grep across all of `orpheus/`: the only "Step"/"LinearDiscontinuous" hits are a
  DOCSTRING example (`cell_update.py:538` shows `class Step(...)`) and a docstring
  mention (`geometry.py:270`). No `NotImplementedError` stub class, no registry entry.
- **What a NEW non-DD occupant must provide (minimum):** `is_linear`,
  `is_positivity_preserving`, `is_affine_scannable` (set `False` for LD/SC — they
  couple ≥2 face moments), `update`, `residual`. That alone gives a constructible
  scheme honored by the per-cell `cell_update.update` path. To run on Cartesian 2-D
  wavefronts it ALSO needs `cell_kernel_batch`/`residual_kernel_batch`. It does NOT
  need the scan triple (those are the DD fast-path only).

## 2. The per-cell dispatch seam — `orpheus/sn/loss_representation.py` (`_OneDimScanWalk`)

`_OneDimScanWalk` (`:1723`). The 1-D sweep is `_run` (`:2380`); the 1-D matvec is
`_apply_walk` (`:1873`). (Public `sweep` `:1740` → `_run`; `loss_action` `:1802`
+ `loss_action_decomposed` `:1826` + `loss_action_transpose` `:2164` → `_apply_walk`.)

### Sweep `_run` — TWO spatial paths:
- **Fast affine-scan path (DD only):** slab joint-batch (`:2520-2592`) and the
  curvilinear non-degenerate per-ordinate scan (`:2707-2756`). Both call
  `ordinate_scan(a_atten, b, psi_in)` (the Blelloch closed form,
  `orpheus/sn/spatial/scan.py:85`) then close the diamond via
  `self.mesh.cell_update.cell_average_from_faces(...)` (`:2570`, `:2735`). This is
  the path that REQUIRES `is_affine_scannable` (it reads the `a_attenuation` cache
  + the scan triple).
- **Slow per-cell path (the seam to generalize):** the degenerate cyl-axis ordinate,
  `:2681-2705`. Branch `if geom.is_degenerate[global_n]:`, builds `visits =
  list(self.mesh.dag_walk(ordinate_idx=…, mu_level_idx=level))`, then
  `for visit in visits: result = cell_update.update(visit=…, total_xs=sig_t_p[:,i],
  source=QV_full[:,i], upstream_state=UpstreamState(spatial_upstream=psi_in,
  angular_upstream=psi_angle[:,i]))`. **This is the ONLY production site that calls
  the polymorphic `cell_update.update` per cell.** A non-affine-scannable scheme
  would route through a GENERALIZED version of exactly this loop — same
  `dag_walk` + `cell_update.update` shape, but gated on `not is_affine_scannable`
  (or on a scheme capability) instead of on `is_degenerate`.

### Matvec `_apply_walk` — already PER-CELL (asymmetry to note):
- The apply walk is genuinely a per-cell loop (`_sweep_direction.for i in cell_indices`,
  `:2003-2035`): `cell_balance_for_streaming(...)` → `m_full = (denom·ψ̄ −
  numer_upstream)/V` → `psi_face_in = cell_update.outgoing_face_from_average(psi_cell,
  psi_face_in)` (`:2032`). It does NOT use `ordinate_scan` (no fast scan on apply).
  So the matvec twin already reads the seam via `outgoing_face_from_average`; a new
  scheme's `residual` would slot here. BUT it currently hardcodes the WDD
  `cell_balance_for_streaming` denom inline — a non-DD scheme needs its `residual`
  routed here, not the shared streaming-balance helper. This is the matvec lockstep
  task for #158.

### The supports gates (where a non-affine-scannable scheme is REJECTED):
- `CumprodScan.supports` (`:701`): `mesh.is_1d AND mesh.cell_update.is_affine_scannable`.
- `ScanMarch.supports` (`:1217`): `(mesh.is_1d OR (is_cartesian AND ndim==2)) AND mesh.cell_update.is_affine_scannable`.
- `_DAGWavefront.supports` (`:789`, base of MovingFrontierWindow/FullFieldWavefront):
  `is_cartesian AND ndim==2` — does NOT gate on affine-scannable (it uses
  `cell_kernel_batch`, not the scan triple). So a non-affine Cartesian-2D scheme
  with `cell_kernel_batch` would ALREADY select through the wavefront family.
- `default_for` (`:1468`) iterates `LOSS_REPRESENTATIONS = (CumprodScan, ScanMarch,
  MovingFrontierWindow, FullFieldWavefront)` (`:1460`), first `supports().ok`.
  **The hole for #158:** a 1-D non-affine-scannable scheme matches NO strategy →
  `IncompatibleRepresentation`. The seam to open = a 1-D per-cell `cell_update.update`
  strategy (a new `LossRepresentation`, e.g. `PerCellWalk`/`DAGScan-1D`) whose
  `supports = mesh.is_1d AND NOT is_affine_scannable`, OR a generalized branch
  inside `CumprodScan`/`_OneDimScanWalk._run` that falls to the degenerate per-cell
  loop for non-affine schemes. Polymorphic-dispatch (a new strategy) matches the
  sn_sweep_strategy.md "construct general, select narrow" principle better than a
  bool branch.

## 3. §1 axis-typing state

- **NO `SpatialScheme` / `AngularScheme` enum exists.** Confirmed. The two ABCs
  ARE the axis disambiguation: `CellUpdateBase` (spatial) vs `PoleAngularClosureBase`
  (angular). There is NO single `LD`/scheme enum conflating the axes — and the user
  already decided (plan session 2) to DROP the enum idea: a spatial-LD is a
  `CellUpdateBase` subclass (#158), an angular-LD is a `PoleAngularClosureBase`
  subclass (#6). They live in different registries with different `registry` dicts,
  so the same display name "LD" on each axis is already unambiguous by type.
- **What "axis-disambiguate LD" concretely needs:** essentially NOTHING structural —
  the type system already separates them. §1 is reduced to: when LD lands, name it
  `LinearDiscontinuous(CellUpdateBase, key="linear_discontinuous")` for #158 and a
  DISTINCT `LinearDiscontinuousAngular(PoleAngularClosureBase, key=...)` for #6, so
  the two registry keys don't collide. The registries are SEPARATE dicts
  (`CellUpdateBase.registry` `:549` vs `PoleAngularClosureBase.registry` `:324`), so
  even an identical key string would not collide. §1 is therefore the cheapest tier
  and is largely a naming-convention decision, not a refactor.

## 4. ANGULAR registry — `orpheus/sn/spatial/pole_angular_closure.py`

- **Protocol** `PoleAngularClosure` (`:192`) + **concrete ABC**
  `PoleAngularClosureBase(RegistryMixin, ABC)` (`:295`, `registry` dict `:324`,
  `is_linear` trait, `@abstractmethod __call__`).
- **Occupants (genuinely TWO + selectable):**
  - `MorelMontryAngularSweep(PoleAngularClosureBase, key="morel_montry_angular_sweep")` (`:488`) — curvilinear default. Real production methods: `precompute_psi_state(:838)`, `cell_contribution(:888)`, `angular_adjoint(:925)`, `psi_half_seed`/`level_indices`. Carries the Carlson coupled-pole seed.
  - `IdentityAngularClosure(PoleAngularClosureBase, key="identity_angular_closure")` (`:1189`) — Cartesian default (zero redistribution); `precompute_psi_state(:1265)`, `cell_contribution(:1276)`, `angular_adjoint(:1288)`.
- **Injection:** `SNMesh(pole_angular_closure=…)` (`geometry.py:204`). Bound at
  `geometry.py:459-465`: user-supplied verbatim, else
  `default_angular_closure_class(self.coord)(self)` —
  `default_angular_closure_class` (`pole_angular_closure.py:1324`): CARTESIAN→Identity,
  SPHERICAL/CYLINDRICAL→MorelMontry.
- **The angular axis IS LIVE in production** (unlike the spatial seam, which is inert
  in the production default sweep). The curvilinear sweep `_run` reads
  `self.mesh.pole_angular_closure` (`:2608`) and calls `closure.psi_half_seed(...)`
  (`:2647`); the matvec `_apply_walk` reads `pole_angular_closure` (`:1940`) and calls
  `precompute_psi_state(:1971)` + `cell_contribution(:2006)`. So a new angular closure
  is honored end-to-end TODAY.
- **How #235 cyl 2-D (η,φ) closure slots in:** subclass
  `PoleAngularClosureBase(key="...")`, implement `__call__` + the three production
  hooks `precompute_psi_state`/`cell_contribution`/`angular_adjoint` +
  `psi_half_seed`/`level_indices`. The φ-azimuthal structure is already carried by
  the quadrature: `Quadrature.product(n_mu, n_phi)`
  (`orpheus/numerics/quadrature/directional.py:545`) builds η=μ_x (radial) × ξ=μ_y
  (azimuthal); `quad.level_indices` (`:451`, one entry per η-level) and
  `quad.reflection_index` (`:365`) are what the M-M closure sweeps per-level today.
  A 2-D (η,φ) closure replaces the 1-D-per-level φ-march with a genuine 2-D angular
  cell update. Inject via `SNMesh(pole_angular_closure=CylinderEtaPhiClosure(...))`;
  no sweep/matvec change needed (the seam is live).

## 5. The #219 MethodSpace-home question — NAME COLLISION CONFIRMED

- **`SNMethodSpace` ALREADY EXISTS** at `orpheus/sn/method_space.py:72` — a frozen
  dataclass carrying `(quadrature, face, inflow_indices, mesh, trace)`. It is the
  BC-realizer argument (`SNBoundaryRealizer.realize(law, method_space)`), NOT the
  aspirational #219 registry-home. Constructed at `geometry.py:571`
  (`SNMethodSpace.for_face(...)` in `_resolve_bcs`). Re-exported from
  `boundary_realizer`. **This is exactly the collision the plan warned about.**
- **NO `MethodSpace` ABC / registry-home exists.** The #219 "MethodSpace foundational
  layer" is purely aspirational — nothing built. User decision (plan session 2):
  #219 is DEFERRED; the registries stay IN-PLACE on `SNMesh` (`cell_update`,
  `pole_angular_closure`). The registry home IS decided: it's `SNMesh`, not a new
  MethodSpace. If #219 ever materializes, it MUST rename the existing BC-realizer
  struct first (rename hazard).
- **Construction/default sites:** `cell_update` default `geometry.py:271-273`
  (→`DiamondDifference()`); `pole_angular_closure` default `geometry.py:459-465`
  (→coord-dispatched). Both also threaded through `SNMesh.from_axes` (`:850-851`,
  `:912-913`) and `__init__` kwargs (`:203-204`).

## 6. Verification hooks for Tier 2

- **#233 spatial characterization gate:** `tests/sn/verification/mms/test_curvilinear_pole_cell_characterization.py`.
  Four L1 tests, all `catches("ERR-059")`: `test_sphere_global_L2_second_order_dual_reference`
  (`:201`), `test_cylinder_global_L2_second_order` (`:256`),
  `test_sphere_pole_cell_first_order_and_Linf_dominant` (`:297`, pole `orders > 0.8`
  LOWER-bound at `:357` — goes GREEN at higher order when #158 lands),
  `test_cylinder_pole_first_order_vs_volume_average_masked_by_midpoint` (`:373`, va
  `orders > 0.8` at `:424`). The pole-cell tests carry NO `verifies(...)` (they pin a
  LIMITATION, not a correctness claim) — a new higher-order spatial scheme makes
  these PASS at order ~2.0 without changing the assertion.
- **Spatial-convergence MMS oracle (DD second-order):**
  `tests/sn/verification/mms/test_mms_curvilinear.py`:
  `test_sn_spherical_mms_converges_second_order` (`:66`),
  `test_sn_cylindrical_mms_converges_second_order` (`:117`) — `catches("ERR-058")`,
  `orders > 1.9` on ladder [20,40,80,160]. A new spatial scheme is verified it stays
  ≥ O(h²) here (must not regress the global L2).
- **Angular-floor / angular-convergence gate (the #235 / #229 target):**
  `tests/sn/verification/mms/test_curvilinear_aniso_convergence.py`:
  `test_cyl_aniso_floor_scales_with_quadrature` (`:107`, `verifies(
  "sn-mms-cylindrical-aniso-spatial-convergence")`, `catches("ERR-026")`) — the
  cylinder azimuthal floor SCALES with `n_phi` (8→16→32: 1.90e-2→7.37e-3→3.10e-3),
  FLAT in `n_mu`. This is the gate #235 lifts: a 2-D (η,φ) closure should drop the
  floor or restore an O(h²) window. Sphere companion
  `test_w1_aniso_sphere_floor_scales_with_quadrature`.
- **Separability probes (`diag_sep_*`): NOT in the tree.** The plan references
  `$CLAUDE_JOB_DIR/tmp/diag_sep_*` as PROMOTION CANDIDATES (Tier-3 §5 gate) — they
  live in job-tmp, never committed. `derivations/diagnostics/` has no `diag_sep_*`.
  Tier 3 must author/promote them fresh.

## SMALLEST FIRST STEP per tier + the one architectural decision each needs

### Tier 2a — #158 (non-DD SPATIAL CellUpdate)
- **Smallest first step:** add a new `CellUpdateBase` occupant in
  `orpheus/sn/spatial/` (e.g. `step.py`: `class StepCharacteristic(CellUpdateBase,
  key="step")` with `is_affine_scannable=False`, implementing `update` + `residual`
  only). Verify in ISOLATION against the per-cell contract (the `update`/`residual`
  round-trip, `cell_update.py:446` invariant) — NO sweep wiring yet. This is a pure
  additive unit, testable without touching `loss_representation.py`.
- **Seam it touches (when wired):** the 1-D per-cell dispatch at
  `loss_representation.py:2681-2705` (generalize the `is_degenerate` per-cell
  `cell_update.update` loop to fire for non-affine schemes) + a new
  `LossRepresentation.supports` admitting `is_1d AND NOT is_affine_scannable` +
  `default_for` ordering. The matvec twin: route the scheme's `residual` into
  `_apply_walk._sweep_direction` (`:2003-2035`) instead of the hardcoded WDD
  `cell_balance_for_streaming`.
- **Architectural decision for the user:** does the non-affine spatial path become a
  NEW `LossRepresentation` strategy (polymorphic, e.g. `PerCellDAGScan` — matches
  sn_sweep_strategy.md "construct general, select narrow"), OR a generalized branch
  inside `CumprodScan`/`_OneDimScanWalk._run`? (Recommended: new strategy — keeps DD's
  fast scan untouched, makes the slow per-cell path first-class and selectable.)

### Tier 2b — #235 (cylinder 2-D (η,φ) ANGULAR closure)
- **Smallest first step:** add a new `PoleAngularClosureBase` occupant in
  `orpheus/sn/spatial/pole_angular_closure.py` (e.g. `class CylinderEtaPhiClosure(
  PoleAngularClosureBase, key="cyl_eta_phi")`) implementing `__call__` +
  `precompute_psi_state` + `cell_contribution` + `angular_adjoint` +
  `psi_half_seed`/`level_indices`. Because the angular seam is ALREADY LIVE in
  production, injecting it via `SNMesh(pole_angular_closure=...)` makes it run
  end-to-end immediately — start by reproducing the M-M result, then improve the
  φ-closure.
- **Seam it touches:** ZERO sweep/matvec wiring needed — the curvilinear `_run`
  (`:2608`, `:2647`) and `_apply_walk` (`:1940`, `:1971`, `:2006`) already read
  `pole_angular_closure`. The φ-azimuthal data is in
  `Quadrature.product`/`level_indices`/`reflection_index`. The only structural
  question is whether the closure needs richer per-(η,φ) connection coefficients than
  the current per-level `cell_contribution` interface exposes.
- **Architectural decision for the user:** does the 2-D (η,φ) closure FIT the existing
  `PoleAngularClosure.__call__` + `cell_contribution(psi_state, i, p, within_positions)`
  per-level interface (φ-march inside one η-level), OR does a genuine 2-D angular cell
  update need a WIDER contract (a 2-D angular stencil, not a per-level 1-D recurrence)?
  This decides whether #235 is an occupant-only add or also an interface widening.

## Cross-cutting anti-surprises
- Spatial seam is INERT in the production DEFAULT sweep (DD inlines its fast scan;
  only the degenerate-cyl ordinate calls `cell_update.update`). Angular seam is LIVE.
  (Plan session-2 empirical finding, still true post-#206.)
- The matvec (`_apply_walk`) is already a per-cell loop (no fast scan on apply), so
  the apply-side already reads `outgoing_face_from_average` — a non-DD scheme's
  `residual` is the apply-side lockstep deliverable.
- 1-D is a parallel-prefix SCAN, not a wavefront DAG — do NOT force a 1-D non-DD
  scheme onto the Cartesian `_DAGWavefront` family (that's `is_cartesian AND ndim==2`
  only). Its home is a 1-D per-cell walk.
