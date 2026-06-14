# #206 Phase C (matvec relocation) — blast-radius / dependency map

> Persisted 2026-06-14 (session 3) from the explorer dispatch during #206 planning. Line numbers
> verified against `refactor/sn-cellupdate-seam-slab` @ `eab05ab` (they have drifted from both the
> issue text and earlier memory — corrected throughout). The 2-D analogue (`3a79ab3`, S6.3b) is the
> template; the 1-D matvec is the last `(L+C)` action still living ON the operator. Companion:
> `.claude/plans/issue_206_carve.md` (the carve plan) + the test-architect spec
> `.claude/agent-memory/test-architect/issue_206_cellupdate_seam_verification.md`.

## 1. `_MSpatialOperatorSum` anatomy

**What it is** (`orpheus/sn/operator.py:304`): `StreamingOperator.M_spatial` — an `OperatorSum`
subclass of two `_SpatialSweepDirection(±1)` summands. Constructed in `__init__` (:334) from
`(sn_mesh, sigma_t)`. It is the operator-algebra home of the 1-D `(L+C)` matvec body (slab + sphere
+ cylinder); 2-D Cartesian is explicitly `NotImplementedError`-guarded (:378, :634, :855) and routes
through `loss_representation` instead.

Properties / surface (`cached_property` on `StreamingOperator`):
- `M_spatial` (operator.py:1580) → builds `_MSpatialOperatorSum(sn_mesh, sigma_t)`
- `M_angular_redist` (operator.py:1598) → `ZeroOperator` (Cartesian) / `AngularRedistributionOperator` (curvilinear, holds a ref to `M_spatial`)
- `loss_representation` (operator.py:1556) → `default_for(sn_mesh)` → `CumprodScan` (1-D) / `ScanMarch` (2-D)

The three methods:
| Method | file:line | Returns | Walk |
|---|---|---|---|
| `_compute_LpC` | operator.py:343 | one `TimedFullField` = full `(L+C)ψ` (bulk + boundary block) | single-emission bidirectional sweep |
| `_compute_decomposition` | operator.py:772 | `(M_spat, M_ang)` tuple | dual-emission, SAME bidirectional sweep |
| `_compute_LpC_transpose` | operator.py:595 | `(L+C)ᵀφ` | reverse-substitution sweep |

**`_compute_LpC` and `_compute_decomposition` are byte-level twins that do NOT share a walk** — the
walk is copy-pasted. Both contain an identical nested `_sweep_direction` closure (`_compute_LpC`
:422–471 vs `_compute_decomposition` :905–962), identical `-1`/`+1` driver + Carlson seed (480–500 vs
968–984), identical degenerate-ordinate branch (503–542 vs 986–1032), identical boundary block
(563–586 vs 1046–1069). ONLY difference: `_compute_decomposition` splits each cell's `m_full` into
`m_ang`/`m_spat` (956) + writes two buffers; `_compute_LpC` writes `m_full` (468). `apply` (:269 /
:1154) recombines the pair. **This is itself a Pattern-2 duplication inside the operator** — #206's
relocation should collapse it.

**The algorithm `_compute_LpC` runs** (density-form `(L+C)ψ`):
- `_sweep_direction(sign, seed)` (422): per μ-level, mask in-direction ordinates, walk
  `dag_walk_cell_indices`; per cell call `pole_angular_closure.cell_contribution` →
  `(angular_denom_term, angular_numer_upstream)`, then `cell_balance_for_streaming` (cell_balance.py:120)
  → `(denom, numer_upstream)`; emit `m_full = (denom·ψ_cell − numer_upstream)/V` (467); advance the WDD
  face `psi_face_in = 2·ψ_cell − psi_face_in` (469). **`denom` here is the #206 duplication site.**
- Driver (Wave-O bare-sweep): backward sweep first, seeded from the GIVEN outer inflow trace
  `face_outer` (480) — the reflective re-apply keystone is DELETED.
- Curvilinear Morel–Montry angular thread: `precompute_psi_state` (417) + per-cell `cell_contribution`;
  slab `IdentityAngularClosure` returns zeros.
- Carlson coupled-pole seed (ERR-058 / #195): forward `+1` seed = backward sweep's pole-face outflow at
  the mirror ordinate `outflow_at_inner.T[quad.reflection_index("x")]` (493); slab reads `face_inner` (499).
- Boundary block (563–586): outflow slots = self-consistency defect `streamed − face_outer`; inflow
  slots = identity `ψ.inflow`. Reflective `−B` is a sibling operator, NOT here.

`_compute_LpC_transpose` (595): reverse-mode adjoint — reversed cell traversal (712 `reversed(cells)`),
boundary block swapped, angular factor reversed via `closure.angular_adjoint` (752), Carlson-seed
adjoint into the mirror ordinate (745). Recomputes `denom` through the SAME `cell_balance_for_streaming`
(718) — a THIRD copy of the denom recurrence.

## 2. Consumer map (L20 three-surface)

### `_compute_LpC` (operator.py:343)
- Production: `CumprodScan.loss_action` (loss_representation.py:738); `ScanMarch.loss_action` 1-D branch (:1362).
- Test (direct): **NONE** — no test calls it by path.
- Test (indirect): `test_loss_action_convention.py` (4 tests, via CumprodScan, :60); `test_bc_extraction_matvec.py` (O.4a.2 carve gate, via `apply`, :12).
- `mcp__nexus__impact`: total_affected = 11, every code node internal. Production reach = the 2 `loss_action` delegations.

### `_compute_LpC_transpose` (operator.py:595)
- Production: `CumprodScan.loss_action_transpose` (:752); `ScanMarch.loss_action_transpose` 1-D branch (:1428).
- Test (direct): NONE.
- Test (indirect): `test_g_adjoint_reciprocity.py` (sphere/cyl), Gate 1.3 (via `apply_transpose`).

### `_compute_decomposition` (operator.py:772) — the ENTANGLED one
- Production: `_MSpatialOperatorSum.apply` (`M_spatial.apply`, :1154); `_SpatialSweepDirection.apply` (transient orchestrator, :269); `AngularRedistributionOperator.apply` (`M_angular_redist.apply`, :1272).
- Test (direct): `test_streaming_operator_decomposition.py::test_subtractive_L_differs_from_matvec_at_zero_sigma_t` (`orch_zero._compute_decomposition`, :333).
- Test (indirect): `_test_helpers._LC_matvec` (:297/:322, route through public `(L+C).apply`); `test_coupled_pole_mu_level_invariant.py` (:22).

**Key L20 finding:** `_compute_LpC`/`_compute_LpC_transpose` have a clean, narrow production surface
(2 `loss_representation` delegations each). `_compute_decomposition` has THREE production consumers
(the `M_spatial`/`M_angular_redist` standalone algebra) — Phase C touches `_compute_LpC` cheaply but
must NOT casually delete `_compute_decomposition`.

## 3. The "single source of truth for denom/a/b" duplication (#206 core)

Duplication is NOW asymmetric — the sweep side already landed its half (the `affine_scan_coefficients`
seam), the matvec has not:
| Path | denom/a/b by | file:line |
|---|---|---|
| Sweep (`_run_1d_sweep`) | `CollisionCache.from_geometry` → `DiamondDifference.affine_scan_coefficients` | sweep_cache.py:457 → diamond.py |
| Matvec (`_compute_LpC`) | `cell_balance_for_streaming` (inline) | operator.py:457 → cell_balance.py:120 |
| Matvec-T (`_compute_LpC_transpose`) | `cell_balance_for_streaming` (inline) | operator.py:718 |

Both compute IDENTICAL `denom = 2|μ|·A_down + dA_w·c_out + σ_t·V` (the closure's `angular_denom_term`
= `(ΔA/w)·c_out`). The WDD denominator is at **three sites** → #206 wants ONE; the sweep already
delegates to `DiamondDifference`, so **Phase C makes the matvec delegate to the SAME seam**.

The matvec/sweep duplication is two-layered: (1) the coefficient recurrence (`denom`/`a`) — sweep done,
matvec pending; (2) the diamond face closure `2ψ̄−in` — BOTH still inline (the seam
`outgoing_face_from_average` exists but neither consumes it). The matvec face-advance
`psi_face_in = 2.0*psi_cell − psi_face_in` (operator.py:469, :960) is `outgoing_face_from_average` un-routed.

## 4. The S6.3b template (`3a79ab3`)

S6.3b moved the `(L+C)` walk OFF the operator INTO `LossRepresentation`. Before: matvec lived on
`StreamingOperator` as one private method per representation (`_apply_1d`, `_apply_2d_cartesian`,
`_apply_2d_cartesian_scanmarch`, `_apply_full_field`, `_apply_1d_transpose`) — 456 lines. After:
- each `LossRepresentation.loss_action(operator, psi)` OWNS the `(L+C)ψ` walk + RETURNS `(L+C)ψ` (NOT `Lψ`).
- `operator.apply = loss_action(self, psi) − σ_t·ψ.bulk` — the Resolution-A `L=(L+C)−C` subtraction is the
  ONLY remaining operator glue, applied ONCE (was 5×). `apply_transpose` mirrors it.
- the 5 `_apply_*` deleted; byte-identical; pinned by `test_loss_action_convention.py` (flat-reflective
  `L·ψ=0 ⟹ loss_action=σ_t·ψ` proves it returns FULL `(L+C)`).

**Shared-frame shape — `_OctantWalk`** (loss_representation.py:445), the 2-D template: ONE traversal
frame (`_interior_walk`, :481); sweep + matvec **fork ONLY at two injected objects** — the cell kernel
(`cell_kernel_batch` solve / `residual_kernel_batch` apply) + the emit policy (`_SweepEmit` sweep /
`(L+C)ψ` bulk + O.4b boundary matvec). **NEVER a boolean `is_solve` flag** (tripwire :466,
`test_one_octant_walk.py`).

**For 1-D the analogous shared frame does NOT exist yet.** `_run_1d_sweep` (scan via `ordinate_scan`)
and `_compute_LpC` (hand-rolled cell loop) are SEPARATE walks. BUT the matvec's per-cell forward loop
`for i in cell_indices: ... psi_face_in = 2ψ−psi_face_in` (operator.py:449–469) IS a first-order
recurrence — **the 1-D matvec CAN be expressed as a `_x_scan_faces`-style scan** (how
`ScanMarch._loss_action_interior` reconstructs faces with α=−1, β=2ψ̄, scan.py:246). The clean 1-D
relocation: matvec consumes `affine_scan_coefficients` + reconstructs faces via `outgoing_face_from_average`
/ an `_x_scan_faces`-style closure — making the matvec's `cell_balance_for_streaming` call redundant
with the sweep cache.

## 5. Coupling hazards (what breaks / must move together)

| Hazard | Site | Risk |
|---|---|---|
| Transient orchestrator | operator.py:268–269 (`_SpatialSweepDirection.apply` builds throwaway `_MSpatialOperatorSum`, calls `_compute_decomposition`) | Standalone per-direction slow path; needs `m_spat+m_ang` recombine. If `_compute_decomposition` retires, rebuild from `_compute_LpC` + angular-only, or retire with consumers. |
| `M_spatial.apply` / `M_angular_redist.apply` | operator.py:1154, 1272 | Depend on `_compute_decomposition`'s DUAL emission + `id(ψ)` cache. The split is FINER than `(L+C)`. **Decision:** keep `_compute_decomposition` as split-emitter OR re-express as `(L+C) − angular` reusing the relocated coefficients (collapses the byte-twin). |
| Curvilinear angular adjoint | `_compute_LpC_transpose:752` `closure.angular_adjoint(...)` | The M-M SECOND triangular factor. Pinned by `test_g_adjoint_reciprocity` (-O-firing since S6.3a). A forward-scan reframe does not auto-give the reverse scan — the transpose may relocate as a sibling. |
| Carlson coupled-pole seed (ERR-058 / #195) | `_compute_LpC:493`, `_compute_decomposition:976`, `_compute_LpC_transpose:745` | `pole_face_seed = outflow_at_inner.T[reflection_index("x")]`. ALREADY single-sourced through the closure (`_run_1d_sweep` uses the SAME `psi_half_seed`). Relocation MUST keep routing through `pole_angular_closure`, not re-inline. Pinned by `test_coupled_pole_mu_level_invariant.py`. |
| `materialize_inverse_cache` | operator.py:1086 | Builds `CollisionCache.from_geometry` from M_spatial σ_t — the existing dual-view bridge; #206 makes it literal. Opportunity, not break. |
| 2-D `NotImplementedError` guards | operator.py:378, 634, 855 | All three raise on 2-D (2-D matvec lives in ScanMarch/_OctantWalk). Relocation must preserve geometry dispatch; `ScanMarch.loss_action` `is_1d` branch (:1359) is the seam. |
| `cell_balance_for_streaming` multi-consumer | cell_balance.py:120 | Also consumed by `DiamondDifference.residual` (per-ordinate `n_mask=1`, the L0/L1 cell-balance reference). #206 must NOT delete it — only stop the MATVEC calling it. |

**Bit-identity gate:** the matvec's `(denom·ψ−numer)/V` vs the sweep's `b=2·QV·inv_denom` scan are
algebraically equal but **not bit-identical** across the V-multiply/÷ → expect principled-equivalence
(per `vv-principles`), not byte-identity, unless the matvec is rewritten to the exact density-form
grouping. The face-closure routing (`2ψ̄−in`) IS bit-identical (same algebra).

## Gaps / notes
- The pre-existing #206 numerics red (`test_unified_cylinder_matches_hand_reference`, ~18% cylinder-matvec)
  lives in the curvilinear path — confirm whether Phase C's 1-D-scan reframe touches the cylinder matvec
  (shares `_compute_LpC` via the curvilinear branch). The cylinder degenerate-ordinate branch (operator.py:503)
  has NO slab-scan analogue — relocating it to a scan frame is the hard part.
- `_compute_LpC` and `_compute_decomposition` are themselves a Pattern-2 duplication (copy-pasted walk) —
  cleanest #206 outcome collapses both.
- Issue #206 text is STALE on locations (`operator.py:904` / `sweep.py:371` / `transport_operator_matvec_unified`
  gone) — live duplication is `cell_balance.py:120` (matvec) vs `diamond.py` `affine_scan_coefficients` (sweep).

**Key files:** `orpheus/sn/operator.py` (343, 595, 772, 1154, 1272, 1556–1631), `orpheus/sn/loss_representation.py`
(445, 582, 727–752, 1359–1433, 1810–2008), `orpheus/sn/spatial/diamond.py`, `orpheus/sn/spatial/cell_balance.py` (120),
`orpheus/sn/spatial/sweep_cache.py`, `tests/sn/_test_helpers.py` (297, 322), `tests/sn/operators/test_bc_extraction_matvec.py`,
`tests/sn/operators/test_streaming_operator_decomposition.py` (333), `tests/sn/operators/test_loss_action_convention.py`,
`tests/sn/sweep/curvilinear/test_coupled_pole_mu_level_invariant.py`.
