# Phase 2 — O.2a: the honest `L+C−S−F−B` driver — L17 convention crosswalk + verification plan

**Author:** test-architect (proactive-trigger dispatch, 2026-06-05)
**Worktree:** `/Users/rodrigo/git/nuclear/ORPHEUS/.claude/worktrees/field-role-typing`
**Branch:** `refactor/field-role-typing` (HEAD `216f4f6`; Phase 1 DONE).
**Authoritative plan:** `.claude/plans/sn_development_sequence.md` §3 Phase 2.
**This file is the L17 artefact — it MUST be written and reviewed BEFORE any Phase 2 production code.**

> Environment discipline (L22): every Read/Grep/Bash MUST target the worktree above.
> The MAIN checkout `/Users/rodrigo/git/nuclear/ORPHEUS` is an OLDER branch — reading it silently reads the WRONG source.
> Test invocation: `PYTHONPATH=<worktree> /Users/rodrigo/git/nuclear/ORPHEUS/.venv/bin/python -O -m pytest ...` (default `-O`; sentinel runs WITHOUT `-O`). Redirect to a file then `echo exit=$?`; never `| tail`. Use `gtimeout`. **NO `continuous_get` in the hot path** (#212 hang) — deselect `test_heterogeneous_*` / `test_keff_slab` / `test_mms_heterogeneous` if they appear in a run.

---

## 0. Re-confirmed anchors (the plan's anchors HAD DRIFTED — these are verified via `inspect.getsourcelines` under worktree PYTHONPATH)

| Plan's stated anchor | ACTUAL location (verified 2026-06-05) |
|---|---|
| `ScatteringOperator` in `orpheus/numerics/operator.py` | **WRONG** — it is `orpheus/sn/scattering.py:423` (`class ScatteringOperator(LinearOperatorMixin)`). `numerics/operator.py` has NO `ScatteringOperator`. |
| `SNStreamingOperator` | renamed → `StreamingOperator` at `orpheus/sn/operator.py:1077`. |
| `_compute_LpC` / `_compute_decomposition` on `SNStreamingOperator` | both are methods on **`_MSpatialOperatorSum`** (`operator.py:336` and `operator.py:578`). `_MSpatialOperatorSum` = `StreamingOperator.M_spatial`, subclass of `OperatorSum` at `operator.py:297`. |
| `_apply_2d_cartesian` | `StreamingOperator._apply_2d_cartesian` at `operator.py:1336` (routes 2-D via `dag_walk` / `graph.residual`). |
| `_scattering_with_boundary_op` (S+B fold) | property on `SNSolver` at `solver.py:407`; body is `self.scattering_op + SNBoundaryOperator(self.sn_mesh)` (`:448`). |
| `_within_group_triple` / `_within_group_krylov` (Phase 1 SSoT) | module functions `solver.py:114` / `solver.py:155`. |
| `_reflect_outflow_into_inflow` (whole-trace −B) | `solver.py:1006`. |
| `OperatorSum` domain-compat check (the tripwire mechanism) | `OperatorSum.__init__` `numerics/operator.py:682`; the check is at the `IncompatibleOperatorComposition` raise inside it. |
| `SNBoundaryOperator` | `orpheus/sn/boundary_operator.py:68`; `domain`/`codomain` = `self.sn_mesh.trace` (`:133`/`:137`); `block_role = BlockRole.BOUNDARY`. |
| `B.reflect_into_inflow(boundary)` (sub-step 4 target) | **DOES NOT EXIST YET** — sub-step 4 creates it. |
| #196 canary `test_phase_e_trajectory_resolvent_flux_shape_crosscheck` | `tests/sn/verification/analytical/test_phase_c_crosscheck.py:619`. **NOTE: this is NOT the right RED-first gate — see §5.** |

### Verified tripwire state (sub-step 1)

```
ScatteringOperator.domain  -> None        (scattering.py:420 fget returns None)
ScatteringOperator.codomain-> None        (scattering.py:424)
ScatteringOperator.block_role = BlockRole.BULK
SNBoundaryOperator.domain  -> sn_mesh.trace (trace space)
SNBoundaryOperator.codomain-> sn_mesh.trace
SNBoundaryOperator.block_role = BlockRole.BOUNDARY
```

`OperatorSum.__init__` skips the domain-compat check because `S.domain is None`:
```python
if (a_dom is not None and b_dom is not None and a_dom != b_dom):
    raise IncompatibleOperatorComposition(...)   # numerics/operator.py
```
Sub-step 1 sets `S.domain = <bulk space>`. Then `S.domain != B.domain` (bulk ≠ trace) and `S + B` throws `IncompatibleOperatorComposition` at construction — exactly the forcing function the docstring at `solver.py:428` documents.

---

## 1. THE L17 CONVENTION CROSSWALK (one row per value crossing a subsystem boundary)

Bridge column applies L18: **fix the convention at the producer (definition site), not at every consumer.**

| # | Value crossing a boundary | Subsystem producer → consumer | Input convention | Internal convention | Output convention | BRIDGE (which site converts; L18 verdict) |
|---|---|---|---|---|---|---|
| C1 | **`S` codomain typing** | `ScatteringOperator` (bulk leaf) → `OperatorSum(S, B)` | TODAY: `S.domain = None` (untyped, predates function-space tagging) | `BlockRole.BULK`; reads/writes bulk `AngularFlux` only | bulk `TimedFullField` | **PRODUCER (`ScatteringOperator.domain` fget, `scattering.py:420`).** Sub-step 1 makes it return the bulk FunctionSpace. Once typed, `S + B` is illegal-by-construction (bulk ⊕ trace mismatch) — illegal state unrepresentable (elegance Pattern 4). |
| C2 | **`B` codomain typing** | `SNBoundaryOperator` (trace leaf) → driver / direct-sum carrier | `B.domain = B.codomain = sn_mesh.trace` (ALREADY typed) | `BlockRole.BOUNDARY`; reads/writes the `.boundary` trace block | trace `BoundaryFlux` | **PRODUCER (already correct).** No bridge change — `B` is the already-typed reference against which `S` typing is checked. |
| C3 | **`S+B` folded source vs honest `−S` / `−B` summands** | driver rhs (`_within_group_triple` arg 2) → SI sweep seed / Krylov matvec | TODAY: `(S+B).apply(ψ)` returns one `TimedFullField` whose `.bulk` = scatter source AND `.boundary` = `B·ψ.outflow` (folded) | within-group rhs: `bulk` read as cell source, `boundary` read as inflow seed | honest direct-sum block `L+C−S−F−B` acting on `(bulk ⊕ trace)` carrier | **BRIDGE = the FOLD ITSELF, which sub-step 2 RETIRES.** Honest composition puts `−S` (bulk) and `−B` (trace) as separate summands of the block operator; `OperatorSum.apply` distributes `S.apply(ψ) + B.apply(ψ)` over disjoint blocks. Producer-side: the block operator owns the disjoint-block sum. |
| C4 | **`−B` delivery: driver-fold route vs whole-trace-helper route (R3)** | `B·ψ.outflow` → inflow seed of the next sweep | route A (driver): `B·ψ.outflow` rides in `rhs.boundary` (the fold), sweep reads inflow slots; route B (helper): `_reflect_outflow_into_inflow` writes `ψ.inflow = B·ψ.outflow` IN-PLACE on the buffer | both call the IDENTICAL `SNBoundaryOperator` (`B` single-sourced); only the plumbing differs | one inflow-seed delivery via the composed block | **PRODUCER (`SNBoundaryOperator`).** Both routes already single-source `B`; sub-step 2 collapses the DRIVER fold route, sub-step 4 collapses the HELPER route into trace-only `B.reflect_into_inflow`. NOTE (L18 nuance): Phase 1 R1 ALREADY routed the direct fixed-source SI loop through the fold; the ONLY surviving `_reflect_outflow_into_inflow` call site is `solver.py:1178` (the eigenvalue final reconstruction sweep) — verified by grep. |
| C5 | **`/W` normalisation (per-ordinate vs iso scalar)** | `ScatteringOperator.apply` → SI sweep / Krylov | source already carries `(P0 in-scatter + (n,2n))/W + Pℓ` divided by `W = Σw` | per-ordinate angular flux, `/W` applied | normalised bulk source, no consumer-side rescale | **PRODUCER (`ScatteringOperator.apply`, R-1 Step 4 A1) — already correct.** Documented at `_within_group_triple` docstring (`solver.py:139`). Phase 2 MUST NOT introduce a consumer-side rescale when composing the honest block. Guards ERR-004 / ERR-025 family. |
| C6 | **Sign convention `(L + C − S − F − B)`** | block operator algebra → matvec/sweep | leaves emit `+L`, `+C`, `+S`, `+F`, `+B` actions | the loss operator subtracts S, F, B: matvec `= (L+C).apply − (S+B).apply − F.apply` | one composed `L+C−S−F−B` | **PRODUCER (the block-operator composition expression).** Today the sign lives in `KrylovAcceleration`/`SourceIteration` (they SUBTRACT the `S` arg). Sub-step 2 moves the sign into the visible `−` of the composed expression (`A_loss = L + C - S - F - B`). Guards failure-mode #1 (sign flip): a flipped `−` is visible in the diff. |
| C7 | **1-D `(L+C)` solve-sweep convention vs 2-D `graph.residual` convention (R4 = #206)** | `_compute_LpC` (1-D matvec) / `_compute_decomposition` (1-D dual-emit) / `InvertibleOperator.solve` (sweep) / `_apply_2d_cartesian` (2-D matvec) | 1-D matvec computes `m_full = (denom·ψ̄ − numer_upstream)/V` in DENSITY units via `cell_balance_for_streaming`; sweep consumes `CollisionCache.from_geometry` precomputed `inv_denom` (#206); 2-D routes via `dag_walk` + `graph.residual` | the WDD cell-balance `denom = 2|μ|·A_down + dA_w·c_out + σ_t·V` is the shared math | one geometry-agnostic `(L+C)` action via `graph.residual` shape extended to 1-D | **BRIDGE = the R4 UNIFICATION (sub-step 3).** Today the denom is computed in TWO utilities (matvec `cell_balance_for_streaming` `spatial/cell_balance.py` vs sweep `CollisionCache.from_geometry` `spatial/sweep_cache.py`) — #206. Producer: ONE `(L+C)` action. The 2-D `_apply_2d_cartesian` already routes through `graph.residual`; 1-D `_compute_LpC` already partially uses `dag_walk` (`operator.py:443`). |
| C8 | **packed vs typed layout** | block operator ↔ `SourceIteration` / `KrylovAcceleration` | D-I.3d (2026-05-29) RETIRED the bare-ndarray packed-vector contract; operators consume/emit `TimedFullField` ONLY (`StreamingOperator.apply` raises `TypeError` on bare ndarray, `operator.py:1232`) | typed `TimedFullField` (`bulk: AngularFlux` ⊕ `boundary: BoundaryFlux`) end-to-end | typed `TimedFullField` | **PRODUCER — already collapsed (D-I.3d / D-J).** Phase 2's composed block MUST keep the typed carrier; the `(bulk ⊕ trace)` direct-sum is the natural home (illegal-state-unrepresentable). NO new packed↔typed bridge. This was the load-bearing R-1 Step 4 mistake — already fixed; Phase 2 MUST NOT regress it. |
| C9 | **preconditioner identity vs block-inverse (#200)** | `_within_group_krylov` → GMRES | TODAY: `preconditioner=lambda q: q` (explicit identity, `solver.py:177`) | GMRES iterates on the un-preconditioned operator-algebra matvec | block-inverse / sweep-as-`(L+C)⁻¹` preconditioner | **PRODUCER (`_within_group_krylov`, the ONE Phase-1 SSoT site).** Sub-step 5 replaces the lambda. The preconditioner must respond to BOTH cell AND face blocks of the GMRES residual (the trace residual `(WDD-propagated face) − ψ_face`) — per #200. This is a CONVERGENCE-RATE change, NOT a fixed-point change. |

### Crosswalk summary (the bridges that move)

- **C1** producer-side typing flip = sub-step 1 (the tripwire).
- **C3 + C6** the fold retires; sign moves into the visible composed expression = sub-step 2.
- **C4** `−B` route collapse: driver fold gone (sub-step 2), helper → trace-only verb (sub-step 4).
- **C7** the 3-impl `(L+C)` unification = sub-step 3 (#206 / #196).
- **C9** preconditioner = sub-step 5 (#200).
- **C2, C5, C8** are ALREADY producer-correct — Phase 2 MUST NOT regress them (they are negative constraints).

---

## 1.5 Claim layer + pillar gate (vv-principles §Hierarchical claim taxonomy — MUST pass before any test row is written)

| Sub-step | Primary claim layer | Reference pillar | Pillar proves the layer? | Structural-independence terminus |
|---|---|---|---|---|
| 1 (S domain) | foundation (software invariant: typed composition is illegal) | n/a (negative test) | yes (raises-on-illegal + must-not-raise-on-legal, per anti-pattern #11) | the type system (FunctionSpace identity) |
| 2 (compose `L+C−S−F−B`) | flux-shape + eigenvalue | bit-identity inheritance (vacuum) + closed-form k_∞ (homogeneous reflective) | yes | `k_∞ = λ_max(A⁻¹F)` (closed-form, NOT another SN solver) |
| 3 (R4 unify `(L+C)`) | flux-shape (Cartesian: principled-equiv; curvilinear: CORRECTNESS change) | MMS (convergence-order/flux-shape) + closed-form k_∞ (eigenvalue) + trajectory-resolvent semi-analytical (curvilinear flux-shape) | MMS proves order/shape NOT eigenvalue → eigenvalue leg uses k_∞ + trajectory-resolvent | trajectory-resolvent Variant α (`power_iterate_variant_alpha`, ZERO `orpheus.sn` imports — verified) + closed-form k_∞ |
| 4 (`B.reflect_into_inflow`) | foundation (trace-only verb equivalence) | bit-identity inheritance | yes (the new verb must reproduce the old helper's buffer write bit-exactly) | the old `_reflect_outflow_into_inflow` output (regression inheritance) |
| 5 (#200 preconditioner) | convergence-rate (iterations-to-converge); fixed-point UNCHANGED | closed-form k_∞ + SI≡Krylov fixed-point | yes | `k_∞` + the un-preconditioned Krylov fixed point (same value, fewer iters) |

**GATE RESULT: PASS.** No row pairs an eigenvalue claim with an MMS-only reference. Sub-step 3's eigenvalue leg is explicitly carried by k_∞ (closed-form) + trajectory-resolvent (semi-analytical), NOT MMS. MMS carries only the convergence-order/flux-shape leg. The chain terminates in structurally-independent grounds (closed-form k_∞ + the import-disjoint trajectory-resolvent frame).

---

## 2. TEST INVENTORY BY ROLE

### 2(a) EXISTING tests that pin legacy behaviour — stay-green vs re-baseline

| Test (file::name) | Disposition | WHY |
|---|---|---|
| `tests/sn/operators/test_scattering_operator.py::*` | **stay-green, then EXTEND** | Pins `ScatteringOperator.apply` semantics. Sub-step 1 adds a `domain`; the apply value MUST NOT change (typing is metadata, not arithmetic). Add a NEW positive+negative pair (§2b N1). |
| `tests/sn/operators/test_sn_boundary_operator.py::*` | **stay-green** | `B` typing unchanged; `B.reflect_into_inflow` (sub-step 4) is NEW surface, add N4 here. |
| `tests/sn/operators/test_streaming_operator_decomposition.py::TestLCDecomposition::test_bit_exact_uniform_sigma_t` | **stay-green, possibly narrow** | Pins `(L+C).apply(ψ) ≡ M(ψ;σ_t)` bit-exact across geometries (`np.testing.assert_array_equal`). R4 (sub-step 3) re-wires the `(L+C)` action through `graph.residual`. If the reduction tree changes, this is the bit-identity tripwire — apply the 3-criteria gate; if it must relax, narrow to `assert_array_almost_equal_nulp(nulp=reduction_depth)` WITH the documented justification. **Verified GREEN today** (`/tmp/p2_core.txt`: 100 passed). |
| `tests/sn/operators/test_streaming_operator_decomposition.py::test_subtractive_L_differs_from_matvec_at_zero_sigma_t` | **stay-green** | Reaches into `_MSpatialOperatorSum._compute_decomposition` directly (`operator.py:327` in test). R4 keeps `_compute_decomposition` as the dual-emit introspection arm; if it is renamed/merged, MIGRATE this test onto the unified action (retirement = test migration). |
| `tests/sn/operators/test_invertible_operator.py::*` (100 tests) | **stay-green** | Pins `InvertibleOperator.solve` (the sweep `(L+C)⁻¹`). Verified GREEN today. R4 unifies the sweep's `(L+C)` with the matvec's; the `.solve` VALUE must not change (the sweep fixed point is the reference). |
| `tests/sn/operators/test_native_matvec.py::*` | **stay-green** | Pins `_compute_LpC` 1-D matvec output. Verified GREEN today. |
| `tests/sn/operators/test_solver_components.py::TestTransportSweep::test_matches_saved_reference` | **already xfail (pre-existing snapshot drift)** | `sweep_ref_2g.npy` predates Wave-T reduction-tree refactors; drifted past `rtol=1e-14`. ORTHOGONAL to Phase 2 — do NOT touch. |
| `tests/sn/solve/test_si_single_primitive_contract.py::*` | **stay-green** | Phase 1 R1 structural spy (ONE `SourceIteration` primitive). Sub-step 2 changes the operator the primitive consumes (composed block vs triple) but NOT the primitive count. |
| `tests/sn/solve/test_fixed_source_2d_equivalence.py::*` | **stay-green** | Phase 1 NEW-1: closed-form `Q/Σt` + 2-D SI≡Krylov twin. Sub-step 2's honest composition MUST keep SI≡Krylov on 2-D. |
| `tests/sn/sweep/core/test_sweep_vs_apply_consistency.py::test_solve_sn_si_vs_krylov_consistency_homogeneous_sphere` | **stay-green (WEAK form)** | Homogeneous reflective sphere, k_eff only (`< 1e-6`). Degenerate (flat eigenmode). Verified GREEN today. The test's own docstring §9 flags the need for a STRONGER heterogeneous flux-shape pin — that is NEW-3 (§2b). |
| `tests/sn/verification/analytical/test_l1_standoff_slab_cylinder.py::test_cylinder_l1_sweep_vs_krylov_twin_path` | **RE-BASELINE: REMOVE the xfail-strict marker** | **CRITICAL FINDING.** This is `@pytest.mark.xfail(strict=True)` but ran **XPASS(strict) → FAILED** today (`/tmp/p2_twin.txt`, exit=1, `[XPASS(strict)]`). The cylinder twin-path ALREADY agrees within `rtol=1e-5` — the B1''/D-K work closed it; the xfail is STALE. Phase 2 sub-step 3 MUST remove the marker. This IS the #196 eigenvalue-twin RED-first gate — see §5. |
| `tests/sn/verification/analytical/test_l1_standoff_slab_cylinder.py::test_cylinder_l1_refinement_both_paths[20/40/80]` | **RE-BASELINE: re-validate, likely remove xfail-strict** | Same stale xfail-strict family; nx=40 twin assertion is identical to the leg above (which XPASSes). Re-run during sub-step 3; remove the marker if it XPASSes at all 3 meshes. (Leg is slow ~80s/nx; full re-validation is a Phase-2 deliverable, not a precondition.) |
| `tests/sn/verification/analytical/test_l1_standoff_slab_cylinder.py::test_cylinder_l1_sweep_vs_trajectory_resolvent` | **stay-green** | Leg 2 of the L14 standoff (sweep ≡ trajectory-resolvent, `rtol=3e-2`). The structurally-independent reference leg. |
| `tests/sn/verification/mms/test_curvilinear_aniso_convergence.py::test_sn_{spherical,cylindrical}_aniso_mms_spatial_convergence_phase_c` | **stay xfail — do NOT expect Phase 2 to flip these** | `@pytest.mark.catches("ERR-026")`, xfail(strict=False). Verified xfailed today with orders `spherical [0.68, −0.08, −0.04]`, `cyl [0.45, 0.10, 0.03]` — STALLED, not O(h²). These use `solve_sn_fixed_source` (auto-flips to **krylov** for curvilinear). Their stall is the **angular-redistribution closure** (Mode-7 ansatz activates `(1−μ²)B/r`) gated on the pole-face spatial closure (#195-adjacent), which the `solve_sn_fixed_source` docstring (`solver.py:1283`) says is a SEPARATE follow-up. **R4 (twin collapse) will NOT make these flip** — see §3 sub-step 3 and §5 the false-flip warning. |
| `tests/sn/verification/mms/test_mms.py::test_sn_1d_slab_mms_converges_second_order` | **stay-green** | Slab O(h²) (`orders > 1.9`). Verified GREEN today (`/tmp/p2_slab2d.txt`). Slab has no angular redistribution → unaffected by curvilinear closure; the canary that R4 must not break on Cartesian. |
| `tests/sn/verification/mms/test_mms_2d.py::test_sn_2d_cartesian_mms_converges_second_order` + `test_sn_2d_cartesian_2g_heterogeneous_mms_converges_second_order` | **stay-green** | 2-D Cartesian O(h²) (1G + 2G het). Verified GREEN today. R4 extends the 2-D `graph.residual` shape to 1-D — these pin the 2-D shape stays correct. |
| `tests/sn/verification/analytical/test_phase_c_crosscheck.py::test_phase_e_trajectory_resolvent_flux_shape_crosscheck` | **stay-green (NOT a #196 flip gate)** | L1 slow, `tol_per_cell` up to 1.2e-1 (12%). Loads SI-generated SN snapshots, compares to trajectory-resolvent at 12%. Too loose to flip on the twin collapse. The plan named this as the #196 canary — **it is NOT the right gate** (§5 explains; the L14 twin-path leg is). |

### 2(b) NEW tests that catch the new convention

| ID | Test | Level | Groups | Geometry | Catches | Crosswalk row |
|---|---|---|---|---|---|---|
| **N1** | `test_scattering_domain_makes_fold_illegal` — positive: `ScatteringOperator(...).domain is <bulk space>` (must equal `B`'s domain class? NO — bulk ≠ trace); negative: `S + SNBoundaryOperator(mesh)` raises `IncompatibleOperatorComposition`. AND a must-NOT-raise positive: `S + CollisionOperator(...)` (both bulk) still composes. | foundation | n/a | slab+sphere | C1; sub-step 1 tripwire; anti-pattern #11 (positive+negative) | C1 |
| **N2** | `test_honest_block_recovers_kinf_homogeneous_reflective_2g` — compose `A_loss = L+C−S−F−B`, solve eigenvalue, assert `k == k_∞ = λ_max(A⁻¹F)` (closed-form, mixture A 2g) at `rtol=1e-10`. **≥2 GROUPS (mandatory).** Both SI and Krylov inner. | L1 | 2g | sphere (reflective) + slab | C3, C6; sub-step 2 honest composition keeps the eigenvalue correct; closed-form pillar (NOT another SN solver) | C3, C6 |
| **N3** | `test_si_vs_krylov_flux_shape_heterogeneous_cylinder_2g` — the STRONG form the existing weak test asks for in its §9. 3-region 2G ABA cylinder; assert SI and Krylov produce the SAME flux SHAPE (per-cell L∞-normalised, `rtol` = twin tolerance, NOT the reference tolerance) AND the SAME k_eff. **This is the #196 flux-shape RED-first gate (§5).** | L1 | 2g | cylinder (reflective, heterogeneous) | C7; sub-step 3 dissolves the SI/Krylov twin by construction | C7 |
| **N4** | `test_reflect_into_inflow_bit_identical_to_legacy_helper` — call `B.reflect_into_inflow(boundary)` and the old `_reflect_outflow_into_inflow(boundary, mesh)` on the SAME outflow trace; assert the inflow slots are `np.array_equal`. Vacuum (B=0 → zeros) + reflective + albedo. | foundation | n/a | slab + sphere | C4; sub-step 4 trace-only verb regression-inherits the helper | C4 |
| **N5** | `test_vacuum_composition_bit_identical_pre_phase2` — snapshot the `solve_sn` k_eff + scalar_flux on a vacuum slab + vacuum 2-D Cartesian at HEAD `216f4f6` BEFORE sub-step 2; after each sub-step assert `np.array_equal` (vacuum: B=0, no reflective inflow → the composition change is a pure re-wiring, bit-identity inheritance). | L1 | 2g | slab + 2-D Cartesian (vacuum) | C3, C6, C7 (vacuum leg of the bit-identity ladder) | C3,C6,C7 |
| **N6** | `test_preconditioner_reduces_iteration_count` — `solve_sn_fixed_source(...).history.n_inner` with the #200 block-inverse preconditioner ≤ `0.75 ×` the un-preconditioned baseline (mixture B reflective, c→1). Fixed point UNCHANGED (k_∞ + SI≡Krylov guards). | L1 | 2g | slab + sphere (reflective) | C9; sub-step 5 is a rate change not a value change | C9 |
| **N7** | `test_LpC_unified_action_one_primitive` — STRUCTURAL spy (foundation): assert that after R4 the matvec hot path (`StreamingOperator.apply`), the dual-emit introspection, and `InvertibleOperator.solve` all delegate to ONE `(L+C)` action object (count the call sites of the shared `graph.residual` closure; no second WDD-denom utility). Mirrors Phase 1 NEW-3. | foundation | n/a | slab+sphere+cyl+2D | C7; the SSoT structural proof for #206 | C7 |

### 2(c) L1 MMS gates that MUST continue passing (the four geometries)

| Geometry | Gate | Current status (verified 2026-06-05) | Phase 2 expectation |
|---|---|---|---|
| 1-D slab | `test_mms.py::test_sn_1d_slab_mms_converges_second_order` | GREEN (orders > 1.9) | **STAY GREEN** — slab has no angular redistribution; R4 must not regress it. |
| 2-D Cartesian | `test_mms_2d.py::test_sn_2d_cartesian_mms_converges_second_order` (+ 2G het) | GREEN | **STAY GREEN** — R4 extends the 2-D `graph.residual` shape; the 2-D MMS is the reference for the shape R4 lifts to 1-D. |
| sphere | `test_curvilinear_aniso_convergence.py::..._spherical_aniso_..._phase_c` | xfail (orders `[0.68,−0.08,−0.04]`, STALLED) | **STAY XFAIL.** Gated on the angular-redistribution / pole-face closure (#195-adjacent), NOT the twin. R4 does NOT flip it. (False-flip warning §5.) |
| cylinder | `test_curvilinear_aniso_convergence.py::..._cylindrical_aniso_..._phase_c` | xfail (orders `[0.45,0.10,0.03]`, STALLED) | **STAY XFAIL** (same rationale). |
| sphere/cyl (isotropic companion) | `test_mms_curvilinear.py::*` | (re-confirm green during Phase 2; isotropic ansatz nulls the redistribution term) | STAY GREEN — the Mode-7 PAIR. If isotropic passes and anisotropic stalls, the bug is in the angular-redistribution path (the pair IS the diagnostic). |

### 2(d) The #196 canary that flips green

The plan names `test_phase_e_trajectory_resolvent_flux_shape_crosscheck` as the canary. **THIS IS WRONG (§5).** The honest #196 RED→GREEN gates are:

- **EXISTING, stale-xfail (REMOVE marker in sub-step 3):** `test_l1_standoff_slab_cylinder.py::test_cylinder_l1_sweep_vs_krylov_twin_path` (XPASS today; flips by un-marking) + `::test_cylinder_l1_refinement_both_paths[20/40/80]`.
- **NEW (write in sub-step 3):** N3 (flux-shape SI≡Krylov on heterogeneous cylinder) — the strong form the weak test's §9 calls for.

---

## 3. PER-SUB-STEP GATE + bit-identity-vs-principled-equivalence determination

For each, the 3-criteria gate (vv-principles §Bit-identity vs principled-equivalence): (1) named principled intermediate, (2) structurally-independent reference, (3) FP-non-associativity dimensionally explainable.

### Sub-step 1 — Give `ScatteringOperator` a domain (the tripwire)
- **Gate:** N1 (positive `S + C` composes; negative `S + B` raises `IncompatibleOperatorComposition`; positive `S.apply` value unchanged). All existing `test_scattering_operator.py` + `test_sn_boundary_operator.py` stay green.
- **Determination:** **BIT-IDENTITY** on `S.apply` (typing is metadata, the arithmetic is untouched — `np.array_equal`). The composition behaviour is a NEW negative test (raises), not a value change. Criteria trivially satisfied (no FP change).
- **Failure-mode coverage:** anti-pattern #11 (positive + negative both required).

### Sub-step 2 — Compose `L+C−S−F−B`; retire the fold + collapse the driver `−B` route (R3)
- **Gate:** N5 (vacuum bit-identity, slab + 2-D) + N2 (k_∞ closed-form, 2G reflective sphere+slab, both SI and Krylov) + existing `test_fixed_source_2d_equivalence` (2-D SI≡Krylov) stays green + slab/2-D MMS stay green.
- **Determination:**
  - **Vacuum (N5): BIT-IDENTITY** (`np.array_equal`). B=0, no reflective inflow → the fold-vs-honest-composition is a pure re-wiring of `+0`; the reduction tree is unchanged. Bit-identity inheritance from the pre-Phase-2 snapshot.
  - **Reflective (N2): PRINCIPLED-EQUIVALENCE.** The honest composition distributes `S.apply(ψ) + B.apply(ψ)` over disjoint blocks instead of folding into one `(S+B).apply`. (1) named intermediates: the per-block bulk source `S·ψ` and trace inflow `B·ψ.outflow` are both reactor-physics quantities (scatter source; reflected current). (2) structurally-independent reference: `k_∞ = λ_max(A⁻¹F)` closed-form, NOT another SN solver. (3) drift bounded by `iteration_count × ULP` (the disjoint-block sum reorders nothing within a block; only the inter-block add order can differ, depth 1). If reflective drift exceeds `rtol=1e-10` vs k_∞, investigate (algorithmic change, not FP).
- **Failure-mode coverage:** #1 sign flip (the visible `−` in `A_loss = L+C−S−F−B`); #6 convention drift (C5 `/W` must stay producer-side).

### Sub-step 3 — Unify the 3 `(L+C)` impls (R4 = #206; fixes #196)
- **Gate:** N7 (structural one-action spy) + `test_streaming_operator_decomposition` bit-exact stays green (or narrows per 3-criteria) + `test_invertible_operator` / `test_native_matvec` stay green + N3 (SI≡Krylov flux-shape heterogeneous cylinder) + the L14 standoff legs (§4) + REMOVE the stale xfail on the twin-path legs.
- **Determination — SPLIT BY GEOMETRY:**
  - **Cartesian (slab + 2-D): PRINCIPLED-EQUIVALENCE (or bit-identity).** R4 routes the 1-D `(L+C)` through the SAME `graph.residual` closure the 2-D path already uses. The FP reduction tree changes (one WDD-denom utility instead of two: `cell_balance_for_streaming` ∥ `CollisionCache.from_geometry` collapse). (1) the named intermediate is the WDD cell-balance `denom = 2|μ|·A_down + dA_w·c_out + σ_t·V` — a single principled quantity (#206). (2) reference: closed-form k_∞ + slab/2-D MMS O(h²). (3) drift bounded by `reduction_depth × ULP` per cell. If `test_streaming_operator_decomposition::test_bit_exact_uniform_sigma_t` (`assert_array_equal`) breaks, apply the 3-criteria and narrow to `assert_array_almost_equal_nulp` WITH documented justification; preserve bit-identity elsewhere.
  - **Curvilinear (sphere + cylinder): CORRECTNESS CHANGE.** This is NOT FP drift — it is the #196 fix. Today SI (WDD sweep, O(h^1.3)) and Krylov (apply-matvec) converge to DIFFERENT fixed points (the twin). R4 makes them ONE action → SI≡Krylov by construction (the asymmetry dissolves). **Structurally-independent reference for the correctness change:** (i) `k_∞ = λ_max(A⁻¹F)` closed-form (eigenvalue leg — MMS does NOT prove eigenvalues); (ii) trajectory-resolvent Variant α (`power_iterate_variant_alpha`, verified ZERO `orpheus.sn` imports — semi-analytical, structurally-independent frame: bouncing-characteristic Green's function vs WDD sweep); (iii) L0 per-ordinate flat-flux streaming-equilibrium (`tests/sn/sweep/curvilinear/test_streaming_equilibrium_curvilinear.py`, the second independent reference). The L14 four-leg standoff (§4) is the instrument.
- **CRITICAL SCOPE FLAG (false-flip warning, §5):** R4 dissolves the SI-vs-Krylov ASYMMETRY (the twin = #196), but does NOT by itself achieve curvilinear MMS O(h²). The anisotropic MMS rate tests are STALLED on a SEPARATE defect (the pole-face / angular-redistribution closure, #195-adjacent, explicitly out of Phase 2 scope per `solver.py:1283`). Do NOT mark those MMS tests as the #196 success criterion.
- **Failure-mode coverage:** #1, #3 (missing ΔA/w — caught by the per-ordinate flat-flux L0); #5 index drift (non-uniform mesh in the refinement leg); Signature 1 (curvilinear sweep divergence).

### Sub-step 4 — Trace-only `B.reflect_into_inflow(boundary)`; retire the zero-bulk-probe shim
- **Gate:** N4 (`B.reflect_into_inflow` bit-identical to `_reflect_outflow_into_inflow`, vacuum + reflective + albedo) + the eigenvalue final-reconstruction-sweep path (solver.py:1178) re-validated (the ONLY surviving caller, verified by grep).
- **Determination:** **BIT-IDENTITY** (`np.array_equal` on the inflow slots). The new verb is a trace-only refactor of the same `R·G` reflection; it must reproduce the helper's buffer write exactly. The zero-bulk probe shim is retired (it constructed a zero `AngularFlux` only to read the boundary block — the trace-only verb skips the bulk allocation).
- **Failure-mode coverage:** #2 variable swap (face inflow/outflow indices — caught by per-face N4); anti-pattern #11 (vacuum positive: must produce zeros).

### Sub-step 5 — #200 real preconditioner
- **Gate:** N6 (`n_inner` ≤ 0.75 × un-preconditioned baseline) + k_∞ guard (fixed point unchanged) + SI≡Krylov fixed-point guard.
- **Determination:** **PRINCIPLED-EQUIVALENCE on the VALUE, rate CHANGE on iterations.** A preconditioner changes the Krylov iteration path, NOT the fixed point. (1) the fixed point IS the named principled quantity (the discrete `(L+C−S−F−B)⁻¹q`). (2) reference: k_∞ closed-form + the un-preconditioned Krylov fixed point (same value). (3) the converged value drift is `iteration_count × condition_number × ULP` — and `iteration_count` DROPS, so drift if anything shrinks. The MEASURAND of success is `n_inner` (iterations), not the value.
- **API GAP (flag for sub-step 5 PR):** eigenvalue `solve_sn` returns `n_inner=None`; `solve_sn_fixed_source(...).history.n_inner` is the measurable seam. N6 uses the fixed-source path. (Cross-ref the SI Gauss-Seidel rate-recovery memo — same measurement gap.)
- **Failure-mode coverage:** convergence-rate regression (the #200 silent-fallback bug is gone post-Phase-1.2; N6 is the re-enable gate).

---

## 4. THE L14 FOUR-LEGGED STANDOFF FOR R4 (curvilinear `(L+C)` unification)

The instrument ALREADY EXISTS: `tests/sn/verification/analytical/test_l1_standoff_slab_cylinder.py`. Phase 2 re-validates and un-gates it. The four legs (cylinder, 3-region 2G ABA — heterogeneous, multi-group, mesh-refined; mandatory per vv-principles H1/H2):

| Leg | Assertion | Reference | Status today | Phase 2 |
|---|---|---|---|---|
| (i) Krylov ≡ independent reference | `\|k_krylov − k_ref\| / k_ref < 3e-2` | trajectory-resolvent Variant α (`solve_greens_function_cylinder_mr`; structurally independent — ZERO `orpheus.sn` imports) | part of refinement leg | re-validate |
| (ii) SI ≡ independent reference | `\|k_sweep − k_ref\| / k_ref < 3e-2` | same trajectory-resolvent | GREEN (`test_cylinder_l1_sweep_vs_trajectory_resolvent`) | stay-green |
| (iii) SI ≡ Krylov (the twin) | `\|k_sweep − k_krylov\| / k_sweep < 1e-5` | each other (twin) — valid ONLY because legs (i)+(ii) anchor to the independent ref | **XPASS today (stale xfail)** | **REMOVE the xfail marker** |
| (iv) all three under refinement | legs (i)(ii)(iii) at `nx ∈ {20, 40, 80}` | trajectory-resolvent + twin | xfail-strict (re-validate) | re-validate; remove marker if XPASS at all 3 |

**Mesh-refinement sequence:** `nx ∈ {20, 40, 80}` (the existing parametrization in `test_cylinder_l1_refinement_both_paths`). The "right rate to right limit": both algorithms converge to the trajectory-resolvent k_ref within its quadrature budget (3e-2), AND agree with each other at 1e-5 (the twin) — the latter is the #196 fix, the former is the structural-independence anchor.

**References (spelled out):**
- **Curvilinear semi-analytical:** `orpheus.derivations.continuous.trajectory_resolvent` Variant α — `power_iterate_variant_alpha` / `solve_greens_function_cylinder_mr`. Verified structurally independent (bouncing-characteristic Green's function; ZERO imports of `orpheus.sn`).
- **Closed-form eigenvalue (the layer MMS cannot reach):** `k_∞ = λ_max(A⁻¹F)`, `A = diag(Σ_t) − SigS^T`, `F = χ⊗νΣ_f` (homogeneous reflective leg, mixture A 2G) — used by N2.
- **L0 second independent reference:** per-ordinate flat-flux streaming-equilibrium (`tests/sn/sweep/curvilinear/test_streaming_equilibrium_curvilinear.py`; Gate 1.1) — proves streaming + redistribution = 0 per ordinate (NOT just summed — H3).

**Why this is L14-complete:** legs (i)+(ii) terminate in a structurally-independent ground (trajectory-resolvent, import-disjoint). Leg (iii) (SI≡Krylov) is the twin, which is NEVER sufficient alone (two SN paths agreeing is cross-implementation agreement) — it is valid ONLY because (i)+(ii) anchor both paths to the independent reference. Leg (iv) is "right rate to right limit" (convergence-order is necessary, never sufficient — the limit is the trajectory-resolvent value, not self-referential).

---

## 5. RED-FIRST GATE PROPOSAL for the #196 fix

### The plan's named canary is the WRONG gate

`test_phase_e_trajectory_resolvent_flux_shape_crosscheck` (`test_phase_c_crosscheck.py:619`) is L1/slow at `tol_per_cell` up to **1.2e-1 (12%)**, loads SI-generated SN snapshots, compares to the trajectory-resolvent at 12%. It is too loose to flip on the twin collapse and it does not isolate the SI-vs-Krylov asymmetry. **Do NOT use it as the #196 RED→GREEN gate.**

### The honest #196 RED-first gate (two parts)

**#196 IS the SI-vs-Krylov asymmetry** (GitHub #196 title: "residual O(h) SI-vs-Krylov WDD spatial-closure asymmetry"). The RED-first gate measures that asymmetry, NOT an absolute O(h²) rate.

**Part A (EXISTING, un-mark in sub-step 3):**
- Gate name: `test_cylinder_l1_sweep_vs_krylov_twin_path` (`test_l1_standoff_slab_cylinder.py:211`).
- Assertion: `abs(k_sweep − k_krylov) / k_sweep < 1e-5` on the 3-region 2G ABA cylinder at nx=40.
- Tolerance: `_CYL_TWIN_RTOL = 1.0e-5`.
- RED-first status: it is `@pytest.mark.xfail(strict=True)` — was RED (twin diverged ~4e-3 pre-D-K). **VERIFIED TODAY: XPASS(strict) → FAILS the suite (`/tmp/p2_twin.txt`, exit=1).** The twin already agrees; the marker is STALE. Sub-step 3 removes the marker → the test goes GREEN.
- This is the eigenvalue leg of #196.

**Part B (NEW, write in sub-step 3) = N3:**
- Gate name: `test_si_vs_krylov_flux_shape_heterogeneous_cylinder_2g` (NEW, in `test_l1_standoff_slab_cylinder.py` or a sibling).
- Assertion: on the 3-region 2G ABA cylinder (heterogeneous, multi-group — H1/H2 satisfied), solve once with `inner_solver="source_iteration"` and once with `"krylov"`; per-group L∞-normalise the scalar flux; assert per-cell `max|φ_SI_norm − φ_kr_norm| < 1e-4` (flux SHAPE agreement, the strong form the weak test's §9 asks for) AND `abs(k_SI − k_kr)/k_SI < 1e-5`.
- Tolerance: flux-shape `1e-4` (well below the O(h) ~`4e-3` twin divergence #196 documents at nx=40; above the iterative-solver convergence floor). k_eff `1e-5` (matches `_CYL_TWIN_RTOL`).
- RED-first status: write it during sub-step 3 WITH `@pytest.mark.xfail(strict=True)` if it is RED before the R4 unification lands (the apply-matvec and sweep are still two `(L+C)` actions → flux shapes differ); flip to un-marked GREEN after R4 (one action → identical shape by construction). If a quick pre-R4 run shows the flux shape ALREADY agrees (as the eigenvalue already does), write it un-marked and treat it as the regression pin that R4 must keep green.
- This is the flux-shape leg of #196 — the one the homogeneous weak test (`test_solve_sn_si_vs_krylov_consistency_homogeneous_sphere`) is degenerate against (flat eigenmode).

### What is NOT a #196 gate (false-flip warning)

The anisotropic curvilinear MMS rate tests (`test_sn_{spherical,cylindrical}_aniso_mms_spatial_convergence_phase_c`, xfail, orders STALLED) are gated on the angular-redistribution / pole-face spatial closure (#195-adjacent), a SEPARATE defect explicitly out of Phase 2 scope (`solver.py:1283`: switching curvilinear to krylov "regresses the rate to O(h^1)" because of the FD outer-face cell-center-as-face approximation). **R4 dissolves the twin but does NOT fix the rate.** If a Phase 2 review claims "#196 fixed → these MMS tests should flip," that is a false-flip — the MMS stall is #195, not #196. Keep them xfail.

---

## 6. GATE SEQUENCE (the order to land + verify Phase 2)

1. **Capture pre-Phase-2 baselines** (HEAD `216f4f6`): N5 vacuum snapshot (slab + 2-D); the operator-algebra-core green set (`test_invertible_operator`, `test_native_matvec`, `test_streaming_operator_decomposition`, `test_si_single_primitive_contract`, `test_fixed_source_2d_equivalence` — verified 100 passed today); slab + 2-D MMS green; the cylinder twin-path XPASS state.
2. **Sub-step 1** (S domain) → gate N1 + existing scattering/boundary tests green. BIT-IDENTITY on `S.apply`.
3. **Sub-step 2** (compose `L+C−S−F−B`, retire fold + driver `−B`) → gate N5 vacuum bit-identical + N2 k_∞ + `test_fixed_source_2d_equivalence` green + slab/2-D MMS green.
4. **Sub-step 3** (R4 unify `(L+C)`, #206/#196) → gate N7 structural spy + decomposition bit-exact (or narrowed) + N3 flux-shape twin + L14 standoff (§4) + REMOVE stale xfail on twin-path legs. **The #196 RED→GREEN gate fires here.**
5. **Sub-step 4** (`B.reflect_into_inflow`) → gate N4 bit-identical + reconstruction-sweep path (solver.py:1178) re-validated.
6. **Sub-step 5** (#200 preconditioner) → gate N6 `n_inner` ≤ 0.75× baseline + k_∞ + SI≡Krylov fixed-point.
7. **Full re-validation:** the L14 refinement leg `nx ∈ {20,40,80}` (slow ~80s/nx); the curvilinear MMS pair STAYS xfail (false-flip check); sentinel set without `-O`.

---

## 7. Self-improvement triggers (fired during this plan)

- **No NEW failure mode** introduced — the failure modes covered (sign flip #1, missing factor #3, variable swap #2, index drift #5, convention drift #6, Signature 1 curvilinear divergence, anti-pattern #11 positive+negative) are all in the vv-principles tables. No skill append required.
- **Test-design correction logged:** the plan's named #196 canary (`test_phase_e_trajectory_resolvent_flux_shape_crosscheck`, 12% tol) is too loose AND conflates #196 (twin) with #195 (pole-face rate). The honest gates are the L14 twin-path legs (eigenvalue) + N3 (flux-shape). This is captured in §5 and in agent memory (the Phase-2 crosswalk memory entry) for future sessions.

---

## Appendix A — verbatim baseline pytest stdout (lesson L12)

### A.1 — curvilinear aniso MMS + homogeneous SI-vs-Krylov (`/tmp/p2_baseline.txt`, exit=0)
```
tests/sn/verification/mms/test_curvilinear_aniso_convergence.py spherical_aniso errors: [0.08362901 0.05210107 0.05498337 0.05636904]
spherical_aniso orders: [ 0.6826905  -0.0776825  -0.03590766]
xcyl_aniso errors: [0.09406382 0.06891607 0.06425937 0.06288598]
cyl_aniso orders: [0.4487995  0.1009337  0.03116841]
x
tests/sn/sweep/core/test_sweep_vs_apply_consistency.py .
=================== 1 passed, 2 xfailed, 1 warning in 38.46s ===================
```
Reading: curvilinear MMS orders STALLED (not O(h²)) — gated on #195 pole-face closure, NOT #196. Homogeneous SI-vs-Krylov passes (degenerate flat eigenmode). 2 xfailed = the curvilinear MMS pair stays xfail.

### A.2 — cylinder twin-path leg (`/tmp/p2_twin.txt`, exit=1 — XPASS(strict))
```
tests/sn/verification/analytical/test_l1_standoff_slab_cylinder.py F
=================================== FAILURES ===================================
__________________ test_cylinder_l1_sweep_vs_krylov_twin_path __________________
[XPASS(strict)] PR-TYPED-6.5 Phase 5 — cylinder twin-path divergence (rel ≈ 4e-3 at nx=40) is the L14 manifestation-#6 signature the B1'' face-state architecture fixes. ...
======================== 1 failed, 1 warning in 53.20s =========================
```
Reading: the xfail-strict twin-path leg now XPASSes — the cylinder SI≡Krylov twin ALREADY agrees within 1e-5. The marker is STALE; Phase 2 sub-step 3 removes it (the #196 eigenvalue gate is already met at nx=40; flux-shape N3 + refinement leg complete the proof).

### A.3 — operator-algebra core green set (`/tmp/p2_core.txt`, exit=0)
```
100 passed, 1 skipped, 1 xfailed, 12 warnings in 11.98s
```
(`test_solver_components::TestTransportSweep::test_matches_saved_reference` is the 1 xfail — pre-existing snapshot drift, ORTHOGONAL to Phase 2.)

### A.4 — slab + 2-D Cartesian MMS O(h²) (`/tmp/p2_slab2d.txt`, exit=0)
```
2 passed, 1 warning in 7.78s
```
(Slab + 2-D Cartesian MMS converge O(h²) — the must-stay-green Cartesian legs.)

### A.5 — refinement leg `test_cylinder_l1_refinement_both_paths[20/40/80]`
```
exit=124
```
`gtimeout` killed the run at the 590s budget (3 nx × ~80s sweep + ~80s krylov + uncached ~30s trajectory-resolvent reference per nx). NOT completed within the bounded session window. The per-nx assertion is identical to A.2 (the same SI≡Krylov twin; at nx=40 A.2 XPASSes). Phase 2 sub-step 3 MUST re-run this leg to completion (budget ≥ 10 min, or pre-warm the reference cache) and remove the xfail marker if it XPASSes at all three meshes (deliverable, not a precondition for landing sub-steps 1–2).
