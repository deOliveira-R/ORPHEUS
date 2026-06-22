---
name: issue-206-cellupdate-seam-verification
description: Verification spec for #206 3-phase SN operator-algebra carve (Phase A closure-route-through-cell_update, Phase B extract 1-D-scan frame, Phase C move 1-D matvec off operator). Per-phase bit-id/principled-equiv gate tables, L14 four-leg standoff coverage map, L17 operator↔loss_representation crosswalk, Mode-9 anisotropic hazard, pre-existing reds to route around.
metadata:
  type: project
---

# Issue #206 cell_update-seam carve — verification spec (test-architect)

Branch `refactor/sn-cellupdate-seam-slab`. The carve routes the shared 1-D scan
core's closures through the `cell_update` scheme (Phase A), extracts the shared
1-D-scan sweep frame (Phase B), then moves the 1-D matvec off the operator
(Phase C) — the 1-D analogue of S6.3b (`3a79ab3`, which moved the 2-D matvec
off the operator). **Why:** the matvec (`operator.py` `_compute_LpC` at
`:469` does inline WDD `psi_face_in = 2·psi_cell − psi_face_in` + `cell_balance_for_streaming`)
and the sweep (`loss_representation.py` `_run_1d_sweep` `:2004`/`:2167` routes
through `ordinate_scan` + `cell_update`) are TWIN PATHS of the same discrete
operator that DON'T share the closure math — Phase A makes them share
`cell_update.outgoing_face_from_average`/`cell_average_from_faces`/`affine_scan_coefficients`
(Pattern 2). **How to apply:** gate Phase A+B as PURE relocations (bit-identical);
Phase C as the four-leg standoff. NO golden regeneration unless an FP-reduction
tree genuinely changes (then narrow per `vv-principles` bit-id-vs-principled).

## The existing gate inventory (verified live in tree, this session)

⭐ **The `--capture-baseline` random-ψ bit-identity gate ALREADY EXISTS** —
`tests/sn/operators/test_bc_extraction_matvec.py` `TestVacuumMatvecBitIdentity`
(`test_vacuum_bulk_bit_identical_1d` slab/sphere/cyl + 2-D-cartesian +
`test_vacuum_boundary_slot_bit_identical_zero_input_1d`). Mechanism: fixed-seed
random ψ → `(L+C).apply` → committed `.npz` under `tests/sn/_data/bc_extraction_baseline/`
→ `np.testing.assert_array_equal` (rtol=0). The `--capture-baseline` flag is
declared in the ROOT `tests/conftest.py` (the gotcha — `pytest_addoption` only
fires there); flag present = WRITE+skip, absent = READ+assert. **This is the
exact gate the prompt asked whether existed: it does. Re-use it, do not rebuild.**

End-to-end CONVERGED-flux sha256 bit-id (≥2G het, `-O`-safe via `raise
AssertionError`, foundation): `tests/sn/solve/test_affine_carve_bit_identity.py`
(slab-2g-het SI + 2-D-2g-aniso-het SI/Krylov). The slab-2g-het row is the
end-to-end PhaseA/B bit-id anchor.

Regression-snapshot harness with the SoT tolerance helper:
`tests/sn/regression/test_dd_regression.py` → `_regression_assert.assert_regression`
(`kind="iterative"` → `SAFETY(10)×conv_tol` read off RUN_CONFIG; `kind="direct"`
→ `nulp(reduction_depth)`; `DriftWarning` tripwire escalatable to strict bit-id
via `-W error::DriftWarning`; `-O`-safe, no bare assert). For PhaseA/B run
`pytest -W error::DriftWarning` to assert STRICT bit-identity (the snapshots
that already pre-drift — the `2d_2g_p1_aniso_dd_8x4_het_si` ~6920 ULP — are
2-D and untouched by the 1-D carve, so the 1-D rows stay strict-clean).

Dual-view cache contract (Pattern 2 anchor):
`tests/sn/sweep/core/test_sweep_cache.py::test_cache_populator_matches_cell_balance_terms`
(cache `(a, 1/denom)` ≡ `cell_balance_terms` at rtol=1e-14). MUST stay green —
it pins that the cache/scan path and the matvec's `cell_balance_for_streaming`
derive the same `(a, denom)`. **This test is the Phase-A invariant**: it already
proves the two algebra surfaces agree; Phase A makes them the SAME surface.

## ⚠ Mode-8 hazard sweep (do BEFORE landing — the gate-runtime-mode decision)

The canonical ORPHEUS invocation is `python -O -m pytest` (memory:
default-test-mode-is-optimize). Bare `assert` is stripped under `-O`. These
carve-relevant gates use BARE assert and are INERT under `-O`:

- `tests/sn/operators/test_streaming_operator.py` — pervasive bare `assert`
  (capability, decomposition-invariant `assert isinstance`, `test_boundary_carries_face_residual` `assert face_max > 1e-12`). The `assert_array_equal`/`assert_allclose` rows DO fire under `-O`; the bare-`assert` rows do not.
- `tests/sn/operators/test_bc_extraction_matvec.py::test_flat_flux_per_ordinate_balance_no_pole_spike` (`assert per_cell[0] <= 2.0*median` — the curvilinear pole-spike L0 diagnostic) — **bare assert, INERT under -O. This is the missing-ΔA/w Mode-3 catcher; if Phase A mis-routes the curvilinear seed this is the gate, but it cannot fire under -O.** Either run this gate WITHOUT `-O` or rewrite to `pytest.fail`/`np.testing.*` at the carve commit (Mode-8 fix-at-touch).
- `test_streaming_operator_decomposition.py` — `assert rel_bulk < 1e-14` etc.

GREEN/SAFE (function-call asserts, fire under -O): `test_g_adjoint_reciprocity.py`
(all `pytest.fail`), `_regression_assert` (all `np.testing.*`/`raise`),
`test_affine_carve_bit_identity.py` (`raise AssertionError`),
`test_bc_extraction_matvec` snapshot rows (`assert_array_equal`),
`test_sweep_vs_apply_consistency` (`assert_allclose`).

**Decision per gate:** the bit-identity / value gates (snapshot, sha256,
regression) all fire under -O. The structural/spike bare-assert gates must run
WITHOUT -O for the carve PR (or be Mode-8-migrated at touch). Stage a one-line
`pytest ... -p no:cacheprovider` run WITHOUT -O on the bare-assert files as a
carve gate.

---

## PHASE A — closure routing (slab-first). Gate table.

Claim layer: **flux-shape + (downstream) eigenvalue**, both INHERITED — the
converged value MUST NOT move. Pillar: bit-identity-by-inheritance (the existing
verified references stay the structurally-independent ground). The carve is a
single-source-of-truth move; principled-equivalence at WORST (if the closure
reorders the FP reduction tree), strict bit-id at BEST.

| Gate | Test (existing) | Tol | What it pins | -O? |
|------|-----------------|-----|--------------|-----|
| A1 sweep snapshot | `test_dd_regression` (1-D slab rows) under `-W error::DriftWarning` | strict bit-id | converged slab flux/keff unmoved | ✓ |
| A2 sweep sha256 | `test_affine_carve_bit_identity::si_slab_2g_het` | sha256 | end-to-end slab-2g-het converged bytes | ✓ |
| A3 matvec random-ψ | `test_bc_extraction_matvec::test_vacuum_bulk_bit_identical_1d[SLB]` | array_equal | matvec bulk on random ψ unmoved | ✓ |
| A4 cache≡balance | `test_sweep_cache::test_cache_populator_matches_cell_balance_terms` | 1e-14 | scan `(a,1/denom)` ≡ `cell_balance_terms` | ✗(allclose ok) |
| A5 sweep≡matvec | `test_loss_action_convention::test_apply_equals_loss_action_minus_collision` (≥2G) | — | `apply.bulk == loss_action − C.apply` | ⚠ check |
| A6 pole-spike | `test_bc_extraction_matvec::test_flat_flux_per_ordinate_balance_no_pole_spike[SPH/CYL]` | 2×median | missing ΔA/w in routed seed (Mode-3) | ✗ bare-assert (run w/o -O) |
| A7 per-ord flat-flux | `test_quadrature.py::test_per_ordinate_flat_flux_consistency` (catches ERR-006/007) | — | per-ordinate streaming+redist=0 | ⚠ check |

**New-gate recommendation A-NEW (the SHARP Phase-A gate, not yet existing):**
a `--capture-baseline` single-MATVEC AND single-SWEEP bit-id gate on a
**heterogeneous, non-flat random ψ** for slab/sphere/cyl ≥2G — `array_equal`
rtol=0 on BOTH `(L+C).apply(ψ).bulk` and ONE `_sweep_1d_unified(Q,σ_t,...)`
output, comparing pre-carve `.npz` (captured at `eab05ab` BEFORE Phase A) to
post-carve. The existing `test_vacuum_bulk_bit_identical_1d` covers the matvec
on random ψ but at VACUUM bc with the (L+C) public surface; it does NOT pin a
single bare `_sweep_1d_unified` call. Add a `tests/sn/sweep/core/` row that
captures the raw sweep output (angular_flux + scalar_flux) on a fixed-seed
random Q + heterogeneous σ_t (so the curvilinear redistribution is ACTIVE — vv
§H2: flat ψ NULLS redistribution → Phase A's curvilinear closure routing would
pass a flat-ψ gate while a real bug hides). Heterogeneous + non-flat is
mandatory because the closure route touches `dA_w·c_out` (curvature
redistribution denom term) which is dead on flat ψ. Use the SoT
`assert_regression(kind="direct", reduction_depth=nx)` so the gate is
ULP-principled if Phase A reorders the reduction, strict otherwise.

## PHASE B — extract shared 1-D-scan frame (PURE RELOCATION). REFINED 2026-06-14 @ HEAD f61e1b0.

Phase A LANDED + committed (`60998e4` A1 closure single-source + `eb2d556`
ScanMarch-2D + `f61e1b0` closeout). The old "A2" matvec-denom single-source
was FOLDED INTO Phase C (user decision); Phase B is SWEEP-SIDE relocation
ONLY (matvec stays on operator until C). The state moved on from the spec
above: the production frame is now the `SweepStrategy` family in
`loss_representation.py` (`CumprodScan:690` / `ScanMarch:1169` /
`FullFieldWavefront` / `MovingFrontierWindow`), NOT free `_sweep_1d_*` only.
But the relocation TARGETS still exist as MODULE-LEVEL free helpers:
`_sweep_1d_unified:1719`, `_ensure_geom_cache:1782`, `_ensure_coll_cache:1791`,
`_run_1d_sweep:1817` — called by BOTH `CumprodScan.sweep:726` AND
`ScanMarch.sweep:1246` (the `is_1d` branch). The 2-D template is `_OctantWalk:445`
(forks at cell-kernel + `_SweepEmit` OBJECTS, NEVER a bool `is_solve`; the
anti-degradation tripwire is `tests/sn/operators/test_one_octant_walk.py`).
Phase B builds the 1-D analogue and relocates the 4 free helpers into it,
shared by CumprodScan + ScanMarch-1D (NOT folded into CumprodScan).

Claim layer: NONE new — pure relocation. The CELL-KERNEL (`ordinate_scan`:
slab joint-batch / curvilinear per-ordinate) + EMIT (`angular_flux` write +
`scalar_flux += w_n·psi_avg`) are the fork; relocation must NOT reorder FP.

⭐⭐ **CRITICAL FINDING (changes the gate design) — `-W error::DriftWarning`
is INERT inside `tests/sn/regression/`.** `tests/sn/regression/conftest.py:25-31`
runs `warnings.simplefilter("always", DriftWarning)` + appends
`always::...DriftWarning` to ini filterwarnings at `pytest_configure`, which
WINS over the command-line `-W error::...` — VERIFIED LIVE: the drifting
`cyl_1g_homogeneous_product_dd_n20` row (scalar 297893 ULP) PASSES under
`-W error::...DriftWarning` (1 passed, 3 warnings). The "strict bit-id via
`-W error`" recipe in the carve plan + the gate table below is a NO-OP for the
regression suite. The genuine strict escalation works ONLY where there is NO
such conftest override — i.e. `tests/sn/sweep/core/` (no DriftWarning conftest;
`-W error::...` DOES escalate there). ⟹ **the A-NEW gate
`test_affine_carve_baseline.py` is the real strict-bit-id Phase-B gate; the
regression suite is a TOLERANCE gate that surfaces (not escalates) drift.**

⭐⭐ **CRITICAL FINDING — the 1-D curvilinear regression rows ALREADY drift at
the pre-Phase-B HEAD (f61e1b0).** Phase A's curvilinear closure routing
reordered the n20 curvilinear reduction tree (within tolerance, NOT
regenerated). Pre-Phase-B "before" baseline (the L12 paste-back, VERIFIED LIVE
2026-06-14):
- `cyl_1g_homogeneous_LS4_dd_n20`: k_eff 1 ULP, scalar **272005 ULP** / 3.56e-11
- `cyl_1g_homogeneous_product_dd_n20`: k_eff 3 ULP, scalar **297893 ULP** / 3.90e-11
- `sphere_2g_homogeneous_dd_n20`: k_eff **3032 ULP**, scalar 30580 ULP / 6.40e-12
- `sphere_2g_p1_aniso_dd_n20`: scalar 278 ULP / 5.20e-14
- ALL SLAB 1-D rows: STRICT-CLEAN (zero DriftWarning).
- 2-D rows (`2d_2g_LS4`/`2d_2g_p1_aniso`) drift too but are UNTOUCHED by the
  1-D carve. The prior spec's "1-D rows stay strict-clean" is WRONG for
  curvilinear post-Phase-A — true for SLAB only.
**Phase B (pure relocation) MUST REPRODUCE THIS EXACT DRIFT PROFILE byte-for-byte
in the warnings summary** (same rows, same ULP counts), NOT achieve zero drift.
A relocation that CHANGES the ULP counts (up OR down) reordered the reduction —
a Phase-B bug. The gate is "warnings summary identical to the before-baseline",
NOT "no warnings".

⭐ The A-NEW gate PASSES strict under `-W error::...DriftWarning` at HEAD (6/6,
VERIFIED) — its small configs (N=4, nx=8, slab/sphere/cyl) are bit-identical to
the pre-carve `be4a57b` (=`eab05ab`-tree) baseline. So Phase A's curvilinear
drift is config-specific (shows at n20 product/LS4, NOT at the A-NEW configs).
The A-NEW gate is therefore the SHARP strict gate; the regression drift is the
tolerance-gate residual to be REPRODUCED.

| Gate | Test (path) | Strict mechanism | What it pins | -O? |
|------|-------------|------------------|--------------|-----|
| B1-strict | `tests/sn/sweep/core/test_affine_carve_baseline.py` (6) under `-W "error::tests.sn.regression._regression_assert.DriftWarning"` | **REAL** (no conftest override in sweep/core) | slab/sphere/cyl single-sweep+matvec bit-id to `be4a57b` baseline | ✓ |
| B1-suite | `tests/sn/sweep/core` (443p/1s/4xf) under same `-W` | **REAL** | every sweep/core invariant strict-bit-id | ✓ |
| B1-regr | `tests/sn/regression/test_dd_regression.py -k "not 2d"` | TOLERANCE (DriftWarning surfaced not escalated — see finding) | 1-D drift profile UNCHANGED vs before-baseline (compare warnings summary) | ✓ |
| B2 import-surface | NEW (recommend, see below) | function-call asserts | relocated symbols resolve from new home AND old call sites still bind | ✓ |
| B3 cumprod≡wavefront | `tests/sn/sweep/core/test_wavefront_cumprod_equivalence.py` (EXISTS) | nulp(128) + closed-form k_inf anchor | CumprodScan(d=1) ≡ FullFieldWavefront(d=1) through the real strategy API | ✓ |
| B4 dispatch pins | `tests/sn/sweep/core/test_unified_sweep_dispatch.py` (12) | pytest.fail | default_for selection + transport_sweep delegation unchanged | ✓ |
| B5 solve bit-id | `tests/sn/solve/test_affine_carve_bit_identity.py` (slab-2g-het sha256) | `raise AssertionError` | end-to-end slab-2g-het converged bytes unmoved | ✓ |

**B2 import-surface gate (the one NEW gate Phase B needs).** Add a class to the
EXISTING `tests/sn/sweep/core/test_unified_sweep_dispatch.py` (it already pins
the strategy import surface + dispatch — natural home) OR a new
`tests/sn/sweep/core/test_one_dim_walk.py` mirroring `test_one_octant_walk.py`
(the `_OctantWalk` anti-degradation tripwire). Assert: (1) the relocated frame
class resolves from `loss_representation` (its new home); (2) `CumprodScan.sweep`
+ `ScanMarch(is_1d).sweep` still produce the documented `(angular,scalar)`
2-tuple via the frame (a thin spy that the frame's interior was hit, mirroring
`test_unified_sweep_dispatch::test_delegates_to_selected_strategy`); (3) the
fork is OBJECTS not a bool — assert `_run_1d_sweep`/`_sweep_1d_unified` no longer
carry an `is_solve`/`is_apply` bool param (signature inspection), mirroring
`test_one_octant_walk`'s anti-bool tripwire. Use `pytest.fail`/`np.testing.*`
(-O-safe). NO golden regeneration.

**B3 cumprod≡wavefront — EXISTS + still pins it.** `test_wavefront_cumprod_equivalence.py`:
`test_cumprod_1d_equals_full_field_spine[vacuum/reflective]` (nulp=128, ≥2G het
non-flat + BC inflow — Mode-9 off the degenerate box) + the closed-form anchor
`test_cumprod_path_hits_analytical_kinf` (k_inf=1.875 transfer-matrix,
structurally-independent ground). GREEN at HEAD (inside the 443 sweep/core
pass). No proposal needed — re-use.

**B4 NEW-gate (frame-level oracle) — NOT needed.** The A-NEW gate
(`test_affine_carve_baseline.py`) already captures the RAW single
`transport_sweep` output (angular+scalar) on het non-flat slab/sphere/cyl and
asserts it bit-identical to the pre-carve baseline. Since `transport_sweep`
routes through `CumprodScan.sweep`→the frame, A-NEW IS the "solve via the new
frame ≡ pre-relocation `_run_1d_sweep` output" oracle the prompt asks about — a
relocation that perturbs the frame would move the A-NEW snapshot and `-W error`
would fire (strict, no conftest override). No additional frame oracle required.

## PHASE C — move 1-D matvec off operator (THE RISKIEST). Four-leg standoff (L14).

`_compute_LpC`/`_compute_LpC_transpose`/`_compute_decomposition` relocate off
`_MSpatialOperatorSum` into the shared 1-D-scan frame / loss_representation.
**Entanglement (must preserve all):** `_compute_decomposition` consumed by the
transient orchestrator (`operator.py:269` `_SpatialSweepDirection.apply`), by
`StreamingOperator.apply` path (`:1154` `_MSpatialOperatorSum.apply` → `m_spat,_`),
by `AngularRedistributionOperator.apply` (`:1272` `_,m_ang`), and by
`StreamingOperator.apply_transpose` (`:1504` → `_compute_LpC_transpose`).
`M_angular_redist` standalone (`:1599`) rides the same decomposition. Shared by
CumprodScan + ScanMarch-1D (both `loss_action` → `operator.M_spatial._compute_LpC`).

### The four-leg coverage map (slab + sphere + cylinder)

- **Leg 1 (sweep ≡ structurally-independent reference):**
  `tests/sn/verification/analytical/test_phase_c_crosscheck.py` — SN 1-D solve
  vs **trajectory_resolvent Variant α Green's-function** (semi-analytical pillar)
  + `kinf_homogeneous` (closed-form pillar). Covers spherical kinf 2G recovery +
  the 5 P0 MR snapshots that admit a Variant α reference. THIS is the leg-1
  ground (NOT another ORPHEUS solver — structurally independent per
  `vv-principles` §1). MUST stay green.
- **Leg 2 (matvec ≡ reference):** `test_bc_extraction_matvec.py`
  `TestStreamingEquilibriumValue::test_flat_flux_per_ordinate_balance_no_pole_spike`
  (Q/Σ_t closed-form value anchor) + the dense-probe matvec reference described
  in that file's docstring (the matvec bulk on random ψ snapshotted at a
  pre-carve HEAD). For the adjoint leaf: `_compute_LpC_transpose` ≡ dense-probe
  transpose oracle `derivations/diagnostics/diag_p42_adjoint_oracle.py`.
- **Leg 3 (sweep ≡ matvec — the TWIN):**
  `test_loss_action_convention.py` (apply ≡ loss_action − C.apply, ≥2G slab
  reflective) + `test_sweep_vs_apply_consistency.py`
  (`test_solve_sn_si_vs_krylov_consistency_homogeneous_sphere` — SI sweep ≡
  Krylov matvec keff) + `test_streaming_operator_decomposition.py`
  (`(L+C).apply ≡ M(ψ;σ_t)` bit-exact, all geometries). The TWIN leg is the
  whole reason for the carve (Phase A unifies the closure; Phase C unifies the
  home) — leg 3 is the regression that the unification held.
- **Leg 4 (under refinement):** `test_phase_c_crosscheck.py` flux-shape rows +
  any mesh-refinement convergence in the 1-D suite. Convergence-RATE necessary
  not sufficient (vv §5) — leg 1 supplies the converged-VALUE ground.

### `_compute_LpC_transpose` (the adjoint) — the curvilinear-angular Mode-9 leg.

`test_g_adjoint_reciprocity.py::test_g_adjoint_reciprocity_full_block`
[slab/sphere/cyl/slab_2g/sphere_2g] — `⟨Aψ,φ⟩_G = ⟨ψ,A.Hφ⟩_G` for
`A=L+C−B`, random NON-flat ψ/φ on bulk+trace, `pytest.fail` (`-O`-SAFE). This
walks `_compute_LpC_transpose` via `A.H.apply`. PAIRED with the L11 negative
control `test_wrong_trace_metric_breaks_reciprocity` (dropping |Ω·n| MUST break
reciprocity — proves §2 not metric-blind, per ERR-051). MUST stay green through
Phase C — it is the only `-O`-firing gate that exercises the curvilinear angular
adjoint (`closure.angular_adjoint`). **This is the load-bearing Phase-C adjoint
gate.**

### Tests that MOVE (L20 test-migration list — reach into `_compute_*` BY PATH)

Most tests reach the matvec through PUBLIC surfaces (`L.apply`,
`(L+C).apply`, `rep.loss_action(L,psi)`) — they survive Phase C transparently.
The genuine path reach-ins to migrate:

- `tests/sn/operators/test_streaming_operator_decomposition.py:328` — COMMENT
  references `_MSpatialOperatorSum._compute_decomposition`; the test BODY uses
  `L.apply`. Update the comment to the new home; body unchanged.
- `tests/sn/_test_helpers.py:297,330` — `_M_matvec`/`_LC_matvec` shims mention
  `_compute_decomposition` in DOCSTRINGS only; bodies route through
  `(L+C).apply`. Update docstrings to new home; bodies unchanged.
- `tests/sn/sweep/curvilinear/test_coupled_pole_mu_level_invariant.py:22-23` —
  docstring names `_compute_LpC`/`_compute_decomposition` as the twins; body
  uses public surface. Update docstring.

**NO test reaches `operator.M_spatial._compute_*(...)` by actual call** in the
1-D path — the matvec internals were already privatised behind `apply`/
`apply_transpose`/`loss_action` (Wave-T T.5 retirement). The L20 migration is
therefore DOCSTRING-ONLY plus the import-surface relink in
`loss_representation.loss_action`/`loss_action_transpose` (`:738`/`:752`) which
call `operator.M_spatial._compute_LpC[_transpose]` — those CALL SITES move with
the methods (intra-production, not a test migration).

### Tests that MUST stay green (do NOT move)

All leg-1..4 tests above + `test_streaming_operator.py` capability/decomposition
+ `test_bc_extraction_matvec` snapshots + `test_affine_carve_bit_identity`.

## L17 crosswalk — operator ↔ loss_representation boundary

| Convention axis | operator.py (`_compute_LpC`) | loss_representation (`_run_1d_sweep`) | Bridge / drift risk |
|---|---|---|---|
| `(L+C)ψ` vs `Lψ` (Resolution A) | returns FULL `(L+C)ψ`; `StreamingOperator.apply` subtracts `σ_t⊙ψ` → `L` | `loss_action` returns `(L+C)ψ`; operator does `−C` | **KEEP**: both return (L+C); the `−C` glue lives at the operator. Phase C must not let the relocated matvec start returning `L` (would double-subtract). Pin: `test_loss_action_convention`. |
| per-ordinate vs iso | matvec is per-ordinate `(ng,N,nx)` throughout | sweep is per-ordinate (slab joint-batched over K; curv per-ord) | both per-ordinate; no iso fold in 1-D. SAFE. |
| WDD face closure | INLINE `psi_face_in = 2·psi_cell − psi_face_in` (`:469`) | `ordinate_scan` + `0.5*(in+out)` (`:2004`/`:2167`) | **THE Phase-A bridge**: route BOTH through `cell_update.outgoing_face_from_average`(=2ψ̄−in)/`cell_average_from_faces`(=½(in+out)). Drift risk = the matvec computes ψ_cell directly (`(denom·ψ−numer)/V`) NOT via the scan recurrence — the closure-route must preserve that the matvec's `2·psi_cell−psi_face_in` equals `outgoing_face_from_average(psi_cell, psi_face_in)` exactly (it does, same algebra). |
| boundary block / O.4b defect | `m_boundary` outflow = `streamed − stored`, inflow = `ψ.inflow` (`:563`) | sweep persists `xmax_face`/`xmin_face` outflow | the boundary-residual convention (outflow defect KEPT, no BC reflection — moves to −B sibling) is matvec-only; the sweep does not emit a residual. Phase C must carry `m_boundary` UNCHANGED. Pin: `test_bc_extraction_matvec` `TestLFullOutflowDefectKept`/`TestVacuumBoundaryDefectKept`. |
| curvilinear angular adjoint | `_compute_LpC_transpose` delegates to `closure.angular_adjoint` (`:613`) | sweep has the FORWARD M-M thread `tau_inv·ψ̄ − mm_a_in·ψ_a_in` (`:2172`) | transpose is matvec-only (sweep is forward). Phase C keeps the reverse factor. Pin: `test_g_adjoint_reciprocity` sphere/cyl. |

**Flag:** the Carlson coupled-pole seed (ERR-058) crosses the boundary —
matvec reads `outflow_at_inner.T[reflection_index]` (`:493`), sweep reads
`pole_outflow[mirror[global_n]]` (`:2110`). Both consume the SAME mirror
permutation + the inward ordinate's pole-face outflow (Pattern 2 — already
unified). Phase C must keep this; the regression is invisible to flat ψ (the
O(h)-wrong pole-cell-centre read was exact on flat ψ — vv §H2). The catcher is
`test_coupled_pole_mu_level_invariant` + the sphere `test_phase_c_crosscheck`
flux-shape rows on NON-flat profiles.

## vv Mode-9 hazard — anisotropic / heterogeneous / non-flat configs ONLY

The 1-D matvec touches the curvilinear M-M angular closure + the Carlson
coupled-pole seed. A wrong angular relocation is INVISIBLE on the degenerate
flat-flux box (vv §H2: flat ψ nulls redistribution; per-ordinate balance
telescopes by construction — vv §H3). Gates that MUST run on stressing configs:

- **Heterogeneous σ_t** (drives `dA_w·c_out` redistribution out of cancellation):
  `test_g_adjoint_reciprocity` uses heterogeneous random σ_t (`_sig_t_heterogeneous`);
  `test_phase_c_crosscheck` MR snapshots are heterogeneous. A-NEW gate MUST use
  heterogeneous σ_t (not uniform).
- **Anisotropic (P1) scattering** (group-coupling + ℓ≥1 moments):
  `test_affine_carve_bit_identity::si_slab_2g_het` is P1; sphere needs an
  aniso companion — `test_curvilinear_aniso_convergence` (if 1-D-reachable) or
  the curvilinear MMS rows. ≥2G mandatory (the 1-group-degeneracy rule — 1G eigenvalue
  degenerate).
- **NON-flat ψ, per-ordinate (L27), NOT weight-summed:** the adjoint
  reciprocity inner is per-ordinate; the pole-spike gate is per-cell-per-ordinate
  L2 (NOT scalar-flux). NEVER assert only on weight-summed scalar flux — a
  μ_x/μ_y swap or wrong mirror index can preserve Σw·ψ while corrupting
  per-ordinate ψ (Mode-2). The decomposition invariant
  `(L+C)ψ ≡ M_spat·ψ + M_ang·ψ` (`test_streaming_operator.py`
  `TestT4cAlgebraDecompositionInvariantCurvilinear`) is per-ordinate-bulk — keep it.

## Pre-existing reds — route AROUND (do NOT conflate with carve regressions)

- `test_streaming_operator.py::test_sphere_1g_apply_bit_identical` /
  `test_sphere_2g_apply_bit_identical` — STALE post-ERR-058 (#195). They use
  `assert_allclose(rtol=1e-13)` against the legacy `pre_t4_snapshots.npz` bundle
  whose sphere arms predate the ERR-058 closure-seed fix (the bundle was
  re-captured 2026-06-12 for slab in the rank-d layout but the sphere arms
  drift at ~1e-13). Confirmed pre-existing. **EXCLUDE from the Phase-C gate
  set** — `-k "not (sphere_1g_apply_bit_identical or sphere_2g_apply_bit_identical)"`.
  (The cylinder arms + the slab arms are the live bit-id pins; sphere bit-id is
  covered by `test_phase_c_crosscheck` value + `test_g_adjoint_reciprocity`
  structure instead.)
- `#212 test_keff_slab::test_heterogeneous_absolute_keff` — `continuous_get`
  hang. **DESELECT** (`--deselect tests/.../test_keff_slab.py::test_heterogeneous_absolute_keff`).
  Orthogonal to the carve.

## Carve gate-run recipe (REFINED 2026-06-14 — Phase-B specific)

⚠ `DriftWarning` MUST use the QUALIFIED path
`error::tests.sn.regression._regression_assert.DriftWarning` — the bare
`error::DriftWarning` raises `AttributeError`. AND ⚠ it is INERT inside
`tests/sn/regression/` (conftest forces `always`); it only ESCALATES under
`tests/sn/sweep/core/` + `tests/sn/solve/`.

1. Phase-B strict bit-id (the REAL gate — sweep/core, where -W escalates):
   `.venv/bin/python -O -m pytest tests/sn/sweep/core tests/sn/solve/test_affine_carve_bit_identity.py -W "error::tests.sn.regression._regression_assert.DriftWarning" -p no:cacheprovider`
   (HEAD baseline: 443p/1s/4xf sweep/core + slab-2g-het sha256 GREEN.)
2. Phase-B drift-profile reproduce (regression, TOLERANCE, 1-D only):
   `.venv/bin/python -O -m pytest tests/sn/regression/test_dd_regression.py -k "not 2d" -p no:cacheprovider`
   then DIFF the DriftWarning summary vs the before-baseline (the 4 curvilinear
   rows + their ULP counts above). Identical summary = pure relocation. Any
   row/ULP change = relocation reordered the reduction (bug).
3. Mode-8 bare-assert gates (NO -O — `test_phase_c_gates.py` Gate 1.2/1.3 use
   bare `assert np.array_equal`/`assert rel<1e-12`; `test_streaming_operator.py`
   pervasive bare assert): `.venv/bin/python -m pytest tests/sn/sweep/core/test_phase_c_gates.py tests/sn/operators/test_streaming_operator.py -k "not (sphere_1g_apply or sphere_2g_apply)" -p no:cacheprovider`
4. `--deselect tests/sn/eigenvalue/test_keff_slab.py::test_heterogeneous_absolute_keff` (#212 hang).
   A-NEW baseline ALREADY captured (`be4a57b`, pre-Phase-A) + committed.

## Self-improvement note (no new failure mode)

This carve introduces NO failure mode absent from the `vv-principles` table —
it is the canonical Pattern-2 twin-path unification (the same shape as Phase F /
S6.3b). The Mode-8 hazard on `test_flat_flux_per_ordinate_balance_no_pole_spike`
(bare-assert L0 pole-spike gate inert under -O) is a candidate for a
Mode-8-fix-at-touch but is already catalogued — no skill-table append needed.
