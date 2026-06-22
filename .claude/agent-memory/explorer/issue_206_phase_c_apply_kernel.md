---
name: issue-206-phase-c-apply-kernel
description: Phase C design — relocate the 1-D (L+C)ψ matvec off _MSpatialOperatorSum onto _OneDimScanWalk as the apply-kernel; cell-kernel/emit fork verdict, _compute_decomposition twin collapse, transpose+Carlson+degenerate handling, ordered relocation sequence.
metadata:
  type: project
---

# Issue #206 Phase C — 1-D matvec apply-kernel relocation onto `_OneDimScanWalk`

Design audit for the surgical carve that moves the 1-D `(L+C)ψ` matvec OFF the
operator (`orpheus/sn/operator.py` `_MSpatialOperatorSum._compute_LpC` /
`_compute_LpC_transpose` / `_compute_decomposition`) INTO the shared 1-D sweep
frame (`orpheus/sn/loss_representation.py` `_OneDimScanWalk`), mirroring the 2-D
`_OctantWalk` carve (S6.3b, commit `3a79ab3`). Line anchors verified at the
post-Phase-B HEAD (graph stale @ `7f4d871`; read source directly).

**Why:** Cardinal Rule 2 — `_compute_LpC` (matvec) and `_OneDimScanWalk.sweep`
(sweep) are TWO directions of the SAME `(L+C)` operator currently living in
SEPARATE objects. They must be ONE walk forking only at cell-kernel + emit, so
sweep≡matvec is a CODE FACT not a coincidence (L21).

**How to apply:** This is the implementation strategy for the highest-risk
Phase-C carve. Main agent implements inline turn-by-turn (see
[[feedback-no-method-implementer-for-surgical-carves]]). Bit-identity on the
closure + nULP on the denom IS achievable (the denom algebra is already
single-sourced).

---

## A. Anatomy of the three operator methods

### `_compute_LpC(psi) → TimedFullField` (operator.py:343-595)
Full `(L+C)ψ` single-emission matvec. Production hot path for all 1-D geometries.

- **Guards**: `curvature ∈ {spherical, cylindrical, cartesian}` else ValueError
  (376-377); `cartesian and not is_1d` → **NotImplementedError** (378-385, the
  2-D dispatch — 2-D routes through the wavefront family, NOT here); slab
  requires `xmin` face populated (407-411).
- **Setup**: `pole_angular_closure` resolved (387-389, defaults to
  `MorelMontryAngularSweep` for curvilinear); `precompute_psi_state` builds the
  per-level Carlson half-grids (417-420). `face_outer = boundary.face_view("xmax")`,
  `face_inner = boundary.face_view("xmin")` (slab only).
- **Core**: nested `_sweep_direction(sign, psi_face_in_init)` (422-473):
  - loops levels → `dag_walk_cell_indices(direction_sign, mu_level_idx=p)` (the
    cell chain order, matches `geom.chain_idx`).
  - per cell: `cell_contribution` → `(angular_denom_term, angular_numer_upstream)`;
    `cell_balance_for_streaming(...)` → `(denom, numer_upstream)`;
    **`m_full = (denom·psi_cell − numer_upstream)/V[i]`** (467 — the apply output);
    advances `psi_face_in = cell_update.outgoing_face_from_average(psi_cell, psi_face_in)`
    (469 — the **α=−1, β=2ψ̄ face advance already seamed in Phase A**).
- **Seed logic**: backward `_sweep_direction(-1, face_outer)` first (482);
  curvilinear pole seed = `outflow_at_inner.T[reflection_index("x")]` (the Carlson
  ψ(0,+μ)=ψ(0,−μ) coupled-pole, ERR-058/#195, 484-495); slab seed = `face_inner`
  (497-501); then forward `_sweep_direction(+1, pole_face_seed)` (502).
- **Degenerate cylinder branch** (504-544): `|μ_x| < eps` ordinates — see §A.5.
- **Emit**: `m_cell = out_g_first.swapaxes(0,1)` (bulk); `m_boundary` = the O.4a.2
  trace block (565-588): OUTFLOW slots = `streamed − given` (self-consistency
  defect); INFLOW slots = `given` (identity). NO BC reflection (that's sibling −B).

### `_compute_LpC_transpose(phi) → TimedFullField` (operator.py:597-772)
Euclidean transpose `(L+C)ᵀφ` — reverse-substitution sweep (O.2b, #208).

- **Guards**: `cartesian and not is_1d` → NotImplementedError (636-641).
- Reverses the boundary writeback (666-684), reverses both DD sweeps per level
  with **`for i in reversed(cells)`** (714), reverses the `psi_face_in=2ψ−in`
  advance (731-733), reverses `m=(denom·ψ − ...)/V` (734-737), accumulates into
  `numer_bar` per level.
- **Seed adjoint**: curvilinear routes the seed cotangent through `mirror[gd]`
  into the −1 reversal (the Carlson seed adjoint, 739-751).
- **Angular factor**: `closure.angular_adjoint(numer_bar, sigma_t)` (754-756) —
  the SECOND triangular factor, zero for slab's IdentityAngularClosure.
- Coefficients recomputed via the SAME `cell_balance_for_streaming` (Pattern 2,
  no twin algebra). Verified bit-for-bit vs dense-probe oracle
  `derivations/diagnostics/diag_p42_adjoint_oracle.py`.

### `_compute_decomposition(psi) → (M_spat, M_ang)` (operator.py:774-1088)
**BYTE-TWIN of `_compute_LpC`** — same `_sweep_direction` (907-966), same seed
logic (968-988), same degenerate branch (990-1036). The ONLY difference is
**dual emission**: per cell it ALSO computes
`m_ang = (angular_denom_term·psi_cell − angular_numer_upstream)/V` (954-957) and
`m_spat = m_full − m_ang` (958), writing both into separate buffers. M_ang.boundary
is zero (BulkOperator per MA-Q4); M_spat carries the same trace block as `_compute_LpC`.

### Consumer / caller map (grep — Nexus `callers` returns empty for these private methods)
- `_compute_LpC`: **loss_representation.py:741** (`CumprodScan.loss_action`),
  **:1367** (`ScanMarch.loss_action`, 1-D branch). 2 production sites.
- `_compute_LpC_transpose`: **loss_representation.py:755** (`CumprodScan.loss_action_transpose`),
  **:1435** (`ScanMarch.loss_action_transpose`, 1-D branch). 2 production sites.
- `_compute_decomposition`: **THREE** production consumers:
  1. **operator.py:1158** — `_MSpatialOperatorSum.apply` (`m_spat, _ = ...`).
  2. **operator.py:1276** — `AngularRedistributionOperator.apply` (`_, m_ang = self.m_spatial._compute_decomposition(psi)`).
  3. **operator.py:269** — `_SpatialSweepDirection.apply` standalone slow path
     (builds a transient `_MSpatialOperatorSum`, reads BOTH, masks by direction sign).
  Plus tests: `test_streaming_operator_decomposition.py:333`.

### The operator glue is ALREADY `loss_action(ψ) − σ_t·ψ`
`StreamingOperator.apply` (operator.py:1427-1506): line 1499
`lpc = self.loss_representation.loss_action(self, psi)`; line 1500
`out_bulk = lpc.bulk.values − self.sigma_t[None]·psi.bulk.values`. Mirror for
`apply_transpose` (1508-1556, `loss_action_transpose`). **So Phase C does NOT
touch the operator glue — it only re-points `loss_action`/`loss_action_transpose`
to own the matvec via `_OneDimScanWalk` instead of delegating to `_compute_LpC`.**

---

## B. The current `_OneDimScanWalk` frame (loss_representation.py:1720-2217)

`@dataclass`-style frame holding `mesh: SNMesh` (the 1-D analogue of `_OctantWalk`).
- `sweep(Q, sig_t, boundary_flux, *, initial_guess)` (1737) → `_ensure_geom_cache`
  + `_ensure_coll_cache` + `_run`.
- `_ensure_geom_cache` (1799) → `GeometryCoefficients.from_mesh_and_quad`, stashed
  on `mesh._geom_cache`.
- `_ensure_coll_cache(sig_t, geom)` (1807) → `CollisionCache.from_geometry(geom,
  sig_t, cell_update)`, stashed on `mesh._coll_cache`. Cache holds `inverse_denom`,
  `a_attenuation`, `cumprod_a` — the **same `(a, 1/denom)` from
  `DiamondDifference.affine_scan_coefficients`** the matvec denom must fold onto.
- `_run(...)` (1832-2217): the geometry-blind sweep body.
  - **SLAB joint-batch** (1972-2044): partition ordinates by sign, ONE
    `ordinate_scan` per direction, `cell_average_from_faces` closure, direct
    `angular_flux[ords]=...` + `scalar_flux += einsum`, persist outflow to
    `xmax_face[ords]`/`xmin_face[ords]`.
  - **CURVILINEAR per-ordinate** (2046-2208): per-level Carlson seed via
    `closure.psi_half_seed(psi_level, CarlsonSweepContext(...))` (2099); per
    ordinate scan with M-M angular thread; degenerate ordinate → slow
    `cell_update.update` per-cell path (2132-2157); pole-outflow capture for the
    mirror coupled-pole seed (2177-2181).

**The apply-kernel seam:** the frame docstring (1730-1732) explicitly reserves
the seam: *"The MATVEC (apply direction) attaches in Phase C as the per-ordinate
apply-kernel (the α=−1, β=2ψ̄ scan), mirroring `_OctantWalk`'s cell-kernel
injection."* Phase C adds a `loss_action(operator, psi)` +
`loss_action_transpose(operator, phi)` method to this frame (or a
`run_apply`/`run_apply_transpose`), forking from `_run` at the cell-kernel + emit.

---

## C. The unification VERDICT — shared frame, fork at cell-kernel + emit (PREFERRED)

**They ARE two directions of ONE walk.** The 2-D family already PROVES this: the
2-D apply-kernel `ScanMarch._loss_action_interior` (loss_representation.py:1372-1428)
and sweep-kernel `ScanMarch._sweep_interior` (:1260-1335) share `_x_scan_faces`
(scan.py:251), which the docstring (scan.py:264-269) says explicitly serves BOTH
the SOLVE coefficients (`α=2s/D−1, β=2(Q+...)/D`) AND the APPLY coefficients
(`α=−1, β=2ψ̄`). The 1-D apply-kernel is the `s_y=0` degeneration of that 2-D
apply-kernel — drop the y-march.

### The two cell-kernels in ONE notation

Let a 1-D cell `i` on a chain have incoming face `in_x`, outgoing face `out_x`,
cell-average `ψ̄`, denom `D = 2|μ|A_down + ang_denom + Σ_t·V`, upstream numerator
`U = |μ|A_total·in_x + ang_numer`.

- **SWEEP (solve) cell-kernel** — unknown is ψ̄, given Q:
  `ψ̄ = (Q·V + U)/D` then `out_x = 2ψ̄ − in_x` (the `ordinate_scan` closed form
  with `a = 2|μ|A_total/D − 1`, `b = 2(QV+ang)/D`). Consumes
  `affine_scan_coefficients` (a, 1/denom) from the σ_t cache. → emits ψ̄ + scalar.

- **APPLY (matvec) cell-kernel** — given ψ̄ (the probe), output is the residual:
  reconstruct `out_x = 2ψ̄ − in_x` (= `outgoing_face_from_average(ψ̄, in_x)`,
  α=−1/β=2ψ̄ scan), then **`m_full = (D·ψ̄ − U)/V`** at zero source. This is
  EXACTLY `_compute_LpC`'s line 467-469. → emits m_full bulk + trace residual.

**Same shape, fork only at the kernel object.** Both consume the SAME `D` (the
matvec's `cell_balance_for_streaming` denom == the scan cache's
`affine_scan_coefficients` denom — cell_balance.py:247 vs diamond.py:460, byte-for-byte
the same `2|μ|A_down + (dA_w·c_out=ang_denom) + Σ_t·V`). Both consume the SAME
`outgoing_face_from_average` for the face advance. The fork is:
- sweep kernel: `D, U` known → divide for ψ̄, then reconstruct out-face;
- apply kernel: ψ̄ known → reconstruct out-face (face scan), then `D·ψ̄ − U` for m_full.

This is the elegant `_OctantWalk` pattern (fork at kernel object, NEVER a bool
`is_solve` flag). The denom/face primitives are already single-sourced, so:

**Bit-identity verdict:**
- The relocated 1-D apply IS `_compute_LpC` moved verbatim into the frame (same
  `cell_balance_for_streaming`, same `cell_contribution`, same chain order via
  `dag_walk_cell_indices`, same `outgoing_face_from_average`). → **bit-identical**
  to the current operator path is achievable IF the relocated kernel keeps the
  exact same op order (the safe relocation: move the body, don't re-derive it).
- A FULLER unification (re-expressing the apply via `affine_scan_coefficients` /
  `_x_scan_faces` like the 2-D kernel) would be **principled-equivalent (nULP),
  NOT bit-identical** — the cache's `1/denom` (one reciprocal) vs the matvec's
  `denom·ψ̄ − U` then `/V` reassociates the FP. RECOMMENDATION: **relocate
  verbatim FIRST (bit-identical), then optionally fold onto the shared cache as a
  SEPARATE principled-equiv step** with the fuller-view oracle pinning it. Do not
  combine relocation + re-expression in one step (two failure modes at once).

### The emit fork (§2)
Sweep emits angular flux + scalar (direct ndarray writes — Phase-B note "DO NOT
retrofit `_SweepEmit`" because the slab joint-batch + curvilinear per-ordinate
write at different granularities). The apply emits `(L+C)ψ` as bulk + boundary.
**1-D apply emit policy: direct ndarray writes** mirroring `_compute_LpC` — write
`out_g_first[:, global_dir, i] = m_full` per cell, assemble `m_boundary` trace
block at the end. Do NOT introduce `_SweepEmit` for 1-D apply (same reasoning as
the sweep). The 2-D apply uses an `LpC` accumulator + `streamed` dict (the
`_OctantWalk.loss_action` frame, :630-637) — the 1-D apply is its scalar-march
analogue: an `out_g_first` accumulator + per-face outflow capture.

---

## D. Concrete ordered relocation sequence

Each step is independently testable and bit-identical (or principled-equiv where
flagged). The operator glue (`StreamingOperator.apply` / `apply_transpose`) is
UNTOUCHED throughout — only `loss_action`/`loss_action_transpose` re-point.

**Step 0 (proactive):** dispatch `test-architect` (the carve crosses the
operator↔loss_representation boundary — a PROACTIVE trigger). The verification
plan: which tests pin the legacy `_compute_LpC` path (test_native_matvec.py
foundation pins, test_loss_action_convention.py, test_streaming_operator_decomposition.py,
test_g_adjoint_reciprocity.py, test_phase_c_gates.py Gate 1.3), the
bit-identity golden, the negative gate (2-D still raises).

**Step 1 — slab forward `loss_action`.** Add `_OneDimScanWalk.loss_action(operator,
psi)` carrying the slab arm of `_compute_LpC` (the slab path has zero angular
contribution — `IdentityAngularClosure.cell_contribution` returns zeros, so no
Carlson/M-M). Re-point `CumprodScan.loss_action` (741) + `ScanMarch.loss_action`
1-D branch (1367) to it. `_compute_LpC` stays (still used by `_compute_decomposition`
consumers). Slab test: test_loss_action_convention.py + native_matvec pins
bit-identical.

**Step 2 — sphere + cylinder forward `loss_action`** (the curvilinear arm:
`precompute_psi_state` + Carlson pole seed + M-M `cell_contribution` + the
degenerate branch §A.5). After this, the relocated `loss_action` covers ALL 1-D
geometries. `_compute_LpC` is now reachable ONLY via `_compute_decomposition`.
Curvilinear matvec tests (test_unified_matvec_cylinder.py,
test_coupled_pole_mu_level_invariant.py) bit-identical.

**Step 3 — transpose.** Add `_OneDimScanWalk.loss_action_transpose(operator, phi)`
carrying `_compute_LpC_transpose` (reverse sweep + `closure.angular_adjoint` for
curvilinear). Re-point CumprodScan:755 + ScanMarch:1435. KEEP routing the angular
factor through `closure.angular_adjoint` and the Carlson seed adjoint through the
mirror permutation (do NOT re-inline). test_g_adjoint_reciprocity.py +
test_phase_c_gates.py Gate 1.3 + the diag_p42 oracle bit-identical.

**Step 4 — collapse the `_compute_decomposition` twin (see §C decision below).**
Re-express M_ang via `cell_contribution` reusing the relocated coefficients
(`m_ang = (angular_denom_term·ψ − angular_numer_upstream)/V`); M_spat =
loss_action − M_ang. Keep the three consumers green (operator.py:269, 1158, 1276).

**Step 5 — retire `_compute_LpC` / `_compute_LpC_transpose` from the operator**
once Steps 1-4 leave them with zero production callers. Migrate any direct-call
tests to the relocated path ([[feedback-retirement-means-test-migration]]).

**Final operator-side surface:** `StreamingOperator.apply`/`apply_transpose`
(the `loss_action ± σ_t·ψ` glue) stay. `_MSpatialOperatorSum` keeps the
per-direction `_SpatialSweepDirection` summands + the decomposition machinery IF
the split stays on the operator (see §C twin decision). `_compute_LpC` /
`_compute_LpC_transpose` GONE (bodies live in `_OneDimScanWalk`). The 2-D
NotImplementedError guards move with the matvec into the frame (the
representation `supports`/`is_1d` already gates 2-D away from `_OneDimScanWalk`).

---

## C-twin. The `_compute_decomposition` collapse decision (the hard call)

`_compute_decomposition` is a byte-twin of `_compute_LpC` that ADDS per-cell
`m_ang = (angular_denom_term·ψ − angular_numer_upstream)/V`. The `m_ang` term is
**purely from `cell_contribution`** — it needs NO spatial scan, just the per-cell
angular contributions + ψ_cell + V. So:

- **`m_ang` is cheaply re-expressible** by a single pass over cells calling
  `cell_contribution` (or by capturing the per-cell `(ang_denom, ang_numer)` the
  relocated `loss_action` already computes). It is ZERO for slab/Cartesian.
- **M_spat = (L+C)ψ − M_ang** = `loss_action(ψ).bulk − m_ang` (the additive
  cell-balance algebra, operator.py:798-807). Bit-exact modulo the per-cell FP
  subtraction (~ULP drift on slab, already documented at operator.py:836).

**RECOMMENDATION:** Re-express `_compute_decomposition` to call the relocated
`loss_action` for the FULL `(L+C)ψ`, then compute `m_ang` via a thin
`cell_contribution` pass and `m_spat = full − m_ang`. This collapses the byte-twin
(the second copy of `_sweep_direction` + seed + degenerate logic) WITHOUT breaking
the split. The decomposition machinery STAYS on the operator (it is operator
algebra — `M_spatial`/`M_angular_redist` are operator leaves), but its WALK
delegates to the relocated `loss_action` (single source).

**Trade-off / alternative:** if the thin `m_ang` pass proves awkward (it needs
the per-level `cell_contribution` re-evaluation, which duplicates the closure
setup), the SAFER fallback is: leave `_compute_decomposition` on the operator
as-is, retire ONLY `_compute_LpC`/`_compute_LpC_transpose` (the FULL-matvec
twins), and have `_compute_decomposition` ALSO delegate its full-`(L+C)` arm to
`loss_action` while keeping its own dual-emission. Either way the THREE
`_compute_decomposition` consumers (269, 1158, 1276) stay green because they call
`_compute_decomposition` (whose signature is unchanged). **The split does NOT
need to move off the operator — only the full-matvec twin must.** This is the
cleaner minimal carve: relocate the full matvec, leave the decomposition's
dual-emission on the operator but source its full-`(L+C)` from the relocated walk.

---

## A.5. The degenerate cylinder ordinate branch (the "hard part")

`_compute_LpC` lines 504-544 (twin at 990-1036). For `|μ_x| < eps` ordinates
(cylinder pure-azimuthal, no radial streaming → NO downstream face, NO scan/slab
analogue):
- Build `global_deg`, `deg_level`, `deg_within` mapping each degenerate global
  ordinate to its level + within-level position (508-517).
- Per cell `i ∈ range(nx)` (NOT a chain — every cell independent): per degenerate
  ordinate call `cell_contribution` (one-element mask, 524-528); then
  `cell_balance_for_streaming(abs_mu=abs_mu_deg, A_downstream=0.0, ...,
  psi_face_in=zero_face, ...)` (533-542) — **A_downstream=0.0 and zero inflow**
  (the geometric truth: no radial face); `m_full = (denom·psi_cell − numer)/V[i]`
  (543), write `out_g_first[:, global_deg, i]`.

**Preservation in the relocated kernel:** the sweep frame ALREADY handles
degenerate ordinates in the SOLVE direction via the slow `cell_update.update`
per-cell path (loss_representation.py:2132-2157, gated `geom.is_degenerate[global_n]`).
The relocated APPLY kernel must carry this branch verbatim — it is `A_down=0`,
`psi_in=0`, per-cell (no chain). Use `geom.is_degenerate` (the same geometric
discriminator the sweep uses) NOT the matvec's inline `|μ_x| < eps` recompute, to
single-source the degeneracy test. This is the one place the apply-kernel cannot
reuse the face-scan (no face to advance); it is a per-cell volumetric balance,
same as the sweep's degenerate arm.

---

## E. DO-NOT list (this carve)

1. **Do NOT re-inline the Carlson coupled-pole seed or the M-M angular closure.**
   Keep routing through `closure.psi_half_seed` / `precompute_psi_state` /
   `cell_contribution` / `angular_adjoint` (ERR-058, #195 — re-inlining the
   pole-cell-centre read was the O(h)-wrong bug that every flat-flux gate missed).
2. **The relocated matvec returns `(L+C)ψ`, NOT `Lψ`.** The `−σ_t·ψ` (Resolution-A
   `−C`) subtraction stays the operator's ONLY glue (operator.py:1500/1550). Do
   not fold it into the kernel.
3. **Keep `cell_balance_for_streaming` AND `affine_scan_coefficients` as the
   single denom source.** They produce byte-identical denoms; the carve must NOT
   spawn a third denom expression. (Optional Step 4b folds the matvec onto the
   cache — separate principled-equiv step, fuller-view oracle pinned.)
4. **Keep `cell_balance_for_streaming` itself** — multi-consumer
   (`DiamondDifference.residual` at diamond.py:267 ALSO uses it at n_mask=1).
   Same for `outgoing_face_from_average`/`cell_average_from_faces` (shared by
   sweep, 2-D march, residual).
5. **Preserve the 2-D NotImplementedError guard** — the relocated `_OneDimScanWalk`
   must never accept a 2-D mesh; the representation `supports`/`is_1d` already
   routes 2-D to the wavefront family. The negative gate (test_native_matvec.py
   pin #7, test_bc_extraction_matvec.py:345-400) must stay green.
6. **Preserve the O.4a.2 boundary trace block semantics** (OUTFLOW=defect,
   INFLOW=identity; NO BC reflection — that's sibling −B). The relocated apply's
   `m_boundary` is the `_OctantWalk.loss_action` trace block (:664-675) verbatim
   in 1-D.
7. **Keep the three `_compute_decomposition` consumers green** (operator.py:269,
   1158, 1276) — they call `_compute_decomposition` by name; its signature must
   not change. The recommended carve sources its full-`(L+C)` arm from the
   relocated walk but keeps its dual-emission + boundary block on the operator.
8. **Relocate verbatim before re-expressing.** Bit-identity first (move the
   body); fold onto the scan cache only as a clearly-separated principled-equiv
   follow-up with an oracle ([[feedback-aggressive-retirement]] fuller-view
   exception applies if you keep `_compute_LpC` as a transient oracle).
