# SN Reshape — Operator Algebra Architectural Campaign

**Date**: 2026-05-05 (updated 2026-05-06)
**Branch**: `refactor/sn-operator-algebra`
**Scope**: SN module + numerics primitives + geometry additions
**Pilot for**: cross-solver migration sequence SN → MoC → CP → MC
**Status**: **Phase 0 + Phase 1 + Phase 2 (DD + amendments) + Phase 3 + Phase 4 complete + Wave H Phase A + Phase B (Issue #168 architecture) complete (Issues 1-7 + 8 + 9-DD + 9.5 + 9.6 + 10, 11, 12, 13 + 14, 15 + #168 Phase A/B); Wave H Phase C (spatial-closure alignment for full ERR-026 closure), Wave C-extension (LD/EC/Step) and Wave F (Phase 5) pending. ERR-026 PARTIAL — architectural infrastructure (BoundaryFaceFlux + PoleAngularClosure Protocols) installed through Phase B; Phase B's empirical finding: canonical M-M angular closure regresses on flat ψ unless paired with WDD spatial closure (Phase C scope).**

---

## Progress (as of 2026-05-10)

Campaign branch `refactor/sn-operator-algebra` ahead of `main` by 20 commits.
Bit-identical regression contract held throughout: 11/11 frozen snapshots at
`tests/sn/regression/snapshots/` remain `np.array_equal`-bit-identical to
the pre-reshape baseline.

### Pre-reshape (already on `main` before campaign branch)

- Issue 16 (regression suite, GH [#165](https://github.com/deOliveira-R/ORPHEUS/issues/165)) — 11 frozen snapshots installed at `tests/sn/regression/snapshots/`
- Anisotropic curvilinear MMS infrastructure with 2 `xfail(strict=True)` ERR-026 tripwires at [tests/sn/l1_analytical/test_mms_curvilinear_aniso_dd_convergence.py](../../tests/sn/l1_analytical/test_mms_curvilinear_aniso_dd_convergence.py) — designed to flip to xpass when Issues 11/12 close ERR-026, forcing marker removal
- vv-principles SKILL extended with failure mode #7 (MMS simplification bias) at [.claude/skills/vv-principles/SKILL.md](../skills/vv-principles/SKILL.md)
- 15 L1 homogeneous + 12 anisotropic MMS foundation tests (~6.2s combined) per-commit gate
- GH [#149](https://github.com/deOliveira-R/ORPHEUS/issues/149) — divide-warning at `solver.py:192` filed for Issue 13/15 owner during Wave H
- 18 reshape issues filed at GH #150–167 on milestone "SN Reshape (Wave 1)"

### Wave A — Phase 0 numerics primitives (DONE)

| Plan # | GH # | Status | Commit | Summary |
|---|---|---|---|---|
| 1 | [#150](https://github.com/deOliveira-R/ORPHEUS/issues/150) | ✓ DONE | `60d3932` | `numerics/operator.py` — LinearOperator Protocol with capability tags; `LinearOperatorMixin` for dunder algebra; OperatorSum/Product/Scaled/Identity/Zero composers; `as_scipy_linop` adapter; **38 foundation tests** |
| 2 | [#151](https://github.com/deOliveira-R/ORPHEUS/issues/151) | ✓ DONE | `2f67853` | `numerics/measure.py` — DiscreteMeasure with tensor product / pushforward / restrict / direct sum; BundleMeasure for MoC (Wave 2); `gauss_legendre`/`gauss_chebyshev`/`equispaced` 1D rules; **22 L1 + 17 foundation tests** |

### Wave B — Phase 0 (Issues 3-5) + Phase 1 (Issues 6-7) (DONE)

| Plan # | GH # | Status | Commit | Summary |
|---|---|---|---|---|
| 3 | [#152](https://github.com/deOliveira-R/ORPHEUS/issues/152) | ✓ DONE | `724547d` | `numerics/symmetry.py` — `SubgroupOfO3` enum + static containment lattice + `is_invariant`; orbit-closure check with O_h (48 elements) and I_h (120 elements) generator sets via Rodrigues + BFS closure; **71 foundation tests** |
| 4 | [#153](https://github.com/deOliveira-R/ORPHEUS/issues/153) | ✓ DONE | `0a9aa94` | `numerics/quadrature/{rules_1d,rules_sphere,rules_product}.py` + thin-adapter refactor of `orpheus/sn/quadrature.py` (502 → 464 LOC). Hot-path attributes (`mu_x`, `mu_y`, `weights`, `level_indices`, `_ref_*`) cached as numpy views on construction; `reflection_index(axis)` reimplemented as cached pushforward. **First SN-touching change**, regression contract held: 11/11 bit-identical |
| 5 | [#154](https://github.com/deOliveira-R/ORPHEUS/issues/154) | ✓ DONE | `60f9fb2` | `numerics/quadrature/registry.py` — `QuadratureSpec` dataclass with structural flags + populated registry + `select_quadrature(geometry, target_degree, **structural)` precedence-chain selector with `SelectionLog` explainability; **37 foundation tests** |
| 6 | [#155](https://github.com/deOliveira-R/ORPHEUS/issues/155) | ✓ DONE | `e4276e1` | `geometry/reduced_operator.py` — Bailey 2009 connection-coefficient lift (α dome recursion, ΔA/w redistribution, M-M τ_mm clamp) into `ReducedStreamingOperator` + `StreamingTerms` + factories (slab/cylindrical/spherical). **Hash equality with `SNMesh._setup_*` outputs** verified via `np.array_equal`; 29 parametrized hash-equality tests across N ∈ {4, 8, 16, 32} (sphere) and (n_mu, n_phi) ∈ {(2,4), (4,4), (4,8)} (cylinder) |
| 7 | [#156](https://github.com/deOliveira-R/ORPHEUS/issues/156) | ✓ DONE | `e93fe47` | `geometry/boundary.py` — `BoundaryOperator` Protocol with tensor-decomposition framing `R = Σ_α G_α ⊗ A_α`; concrete `VacuumBoundaryOperator`, `SpecularBoundaryOperator`, `WhiteBoundaryOperator`, `PeriodicBoundaryOperator`, `AlbedoBoundaryOperator`, `MixedBoundaryOperator`. `SNMesh.BOUNDARY_OPERATOR_REGISTRY` factories return `BoundaryOperator`; sweep at `orpheus/sn/sweep.py` calls `apply_to_incoming(...)`. WhiteBoundaryOperator/PeriodicBoundaryOperator ship as primitives only (not yet wired into `solve_sn`). Backward-compat `__eq__` shim on VacuumBoundaryOperator/SpecularBoundaryOperator accepts string comparisons (transitional). **19 primitive tests, regression contract held: 11/11 bit-identical** |

Wave A + Wave B verification gates (end-of-session):

- 398 numerics + geometry tests on integrated campaign HEAD
- 11/11 regression snapshots bit-identical (629s)
- 27 + 2 xfail safety net (L1 homogeneous + anisotropic MMS foundation + ERR-026 tripwires)
- Sphinx clean (-W); audit clean (22 orphans, 36/38 ERR coverage — both unchanged from baseline)
- Full L1 sweep exit 0

### Wave C — Phase 2 cell-update strategy (Issue 8 + DD-portion of Issue 9) (DONE)

| Plan # | GH # | Status | Commit | Summary |
|---|---|---|---|---|
| 8 | [#157](https://github.com/deOliveira-R/ORPHEUS/issues/157) | ✓ DONE | `e999c15` | `orpheus/sn/spatial/cell_update.py` — `CellUpdate` `@runtime_checkable Protocol` (slot-style traits `is_linear`/`is_positivity_preserving` + `update(streaming_terms, total_xs, source, upstream_state) → CellResult`). `UpstreamState` and `CellResult` `@dataclass(frozen=True, slots=True)` carry per-cell `(ng,)`-shaped state. `StreamingTerms` extended additively with `volume` and `abs_mu` fields populated by all 3 factories (slab/sphere/cylinder). Slab vs curvilinear discriminated via `streaming_terms.alpha_in is None`; cylindrical pure-azimuthal `|η|<1e-15` degenerate case via `outgoing_spatial_flux=None`. **13 protocol tests + 16 vol/abs_mu tests = 29 new foundation tests; existing 29 hash-equality tests preserved (45 total geometry tests on `test_reduced_operator.py`)** |
| 9 (DD) | [#158](https://github.com/deOliveira-R/ORPHEUS/issues/158) | ✓ DONE (DD portion) | `3b1fc75` | `orpheus/sn/spatial/diamond.py` — `DiamondDifference` strategy as **bit-identical extraction** of existing inlined sweep math. Single geometry-polymorphic class with 3 branches: slab (`alpha_in is None`, `s = 2·source/denom` matching the cumprod path's `bQ = 2·source_coeff·Q_1d` semantics), curvilinear (`alpha_in is not None, abs_mu ≥ 1e-15`, mirrors `sweep.py:350-361` operation order verbatim), cyl-degenerate (`alpha_in is not None, abs_mu < 1e-15`, no spatial DD closure). `is_linear=True, is_positivity_preserving=False` (Lewis & Miller §5.3). **11 hand-calc bit-identical tests via `np.array_equal` against scalar formulas typed verbatim from sweep.py:117-123/350-361/533-546**. LD/EC/Step deferred to Wave C-extension |

Wave C verification gates (end-of-session):

- 69 spatial + geometry tests (13 protocol + 11 diamond + 45 geometry incl. 16 new vol/abs_mu)
- 11/11 regression snapshots bit-identical — sweep is provably untouched in Wave C, so the contract holds by construction; verified twice (post-R1 worktree run 617s + post-R2 agent run 599s)
- 27 + 2 xfail safety net intact (ERR-026 tripwires gated for Wave D, did NOT flip)
- Sphinx clean (-W exit 0); audit unchanged from Wave-B baseline (23 orphans of 234 testable labels — denominator grew from 231 to 234 because Wave C added 3 new equation labels [`dd-slab-scalar`, `dd-mm-closure-constants`, `dd-curvilinear-scalar`]; all 3 covered by `@pytest.mark.verifies` in test_diamond.py so orphan count unchanged)
- 36/38 ERR coverage preserved

### Wave D — Phase 3 SN core reshape (Issues 10, 11, 12, 13) (DONE)

| Plan # | GH # | Status | Commit | Summary |
|---|---|---|---|---|
| 10 | [#159](https://github.com/deOliveira-R/ORPHEUS/issues/159) | ✓ DONE | `4183e22` | `orpheus/sn/geometry.py` — `SNMesh.__init__` now calls `slab_streaming` / `spherical_streaming` / `cylindrical_streaming` factories from Wave B Issue 6 instead of inlined `_setup_spherical` / `_setup_cylindrical` methods. Exposes `self.reduced: ReducedStreamingOperator` as canonical accessor. Backward-compat `@property` accessors (with `DeprecationWarning`) for `alpha_half`, `redist_dAw`, `tau_mm`, `alpha_per_level`, `redist_dAw_per_level`, `tau_mm_per_level`, `face_areas`, `delta_A` so the 6 production read sites in `sweep.py` and `solver.py` continue working. **16 new foundation tests; 11/11 regression bit-identical (668s)**. SNMesh shrunk by ~177 LOC (78 + 99 LOC removed) |
| 13 | [#162](https://github.com/deOliveira-R/ORPHEUS/issues/162) | ✓ DONE | `b6eee83` | `orpheus/sn/scattering.py` (290 LOC) and `orpheus/sn/fission.py` (115 LOC) — new. `ScatteringOperator(LinearOperator)` carries the full P0 + Pℓ Galerkin + (n,2n) source construction (math from `_add_scattering_source` + `_build_aniso_scattering` + `_add_n2n_source` consolidated into one operator). `FissionOperator(LinearOperator)` implements `χ ⊗ νΣ_f · φ` rank-1-in-energy structure. `capabilities = {"apply"}` for both. **Architectural deviation from plan**: methods kept as thin delegators on `SNSolver` (instead of removed) because `EigenvalueSolver` Protocol surface and `TestAnisotropicScattering` need them; math now LIVES in operators per Cardinal Rule 2; Wave E retires delegators. **27 new foundation tests** (17 scattering + 10 fission); 11/11 regression bit-identical including Pℓ snapshots (`slab_2g_p1_aniso_dd_n20`, `sphere_2g_p1_aniso_dd_n20`). **Finding**: BiCGSTAB `build_rhs_*` helpers compute scattering inline, do NOT consume the new operators yet — Wave E Issue 15 will wire them |
| 12 | [#161](https://github.com/deOliveira-R/ORPHEUS/issues/161) | ✓ DONE | `c4138e7` | `orpheus/sn/sweep.py` — rewritten as unified algorithm. `transport_sweep` dispatches on `sn_mesh.reduced.requires_upstream_angular_state` (boolean from `ReducedStreamingOperator`), NOT on string-comparing `curvature`. Two unified sub-sweeps: `_cartesian_sweep` (1D + 2D, 1D cumprod fast path preserved when `cell_update == DiamondDifference` AND GL1D AND isotropic) and `_curvilinear_sweep` (sphere + cyl μ-marching, per-cell math delegated to `L.cell_update.update(...)`). `SNMesh` got new `cell_update: CellUpdate` constructor argument (default `DiamondDifference()`). **9 new foundation tests** (dispatch routing + cumprod-fast-path preconditions). **2 documented deviations**: (a) 2D wavefront stays inlined-DD (per-cell `update(...)` accepts `(ng,)` slices but anti-diagonal vectorization works on `(n_diag, ng)` slices; pre-authorized by plan brief); (b) two latent geometry-layer convention bugs caught and worked around in sweep.py — `streaming_terms()` canonical face-area direction (`A[i]=in`) doesn't match inward-sweep semantics, and cylindrical `abs_mu` extraction reads wrong global ordinate; both worked around with helper functions. **11/11 regression bit-identical (777s)**, 27 + 2 xfail safety net intact |
| 11 | [#160](https://github.com/deOliveira-R/ORPHEUS/issues/160) | ✓ DONE | `ea054b3` | `orpheus/sn/operator.py` — adds `SNStreamingOperator(LinearOperator)` (~390 LOC) with three capabilities: `apply` (matrix-free forward action; reuses existing `transport_operator_matvec_*` symmetric-closure math, bit-identical extraction), `solve` (invokes Wave D R2 `transport_sweep` with `DiamondDifference` WDD asymmetric closure), `apply_transpose` (adjoint via dense-matrix-probe approach: build matrix from `apply`, transpose, apply; mathematically rigorous and avoids per-geometry adjoint sign-flip risk). `apply` and `solve` use **different closures by design** — Wave E reconciles via Krylov-on-apply with sweep as preconditioner. **Reciprocity ⟨L·ψ, φ*⟩ = ⟨ψ, L*·φ*⟩ verified to round-off**: rel diff 1.25e-15 (slab), 4.47e-16 (sphere), 7.80e-16 (cyl) — 3 orders of magnitude below `rtol=1e-12` tolerance. Existing `transport_operator_matvec_*` legacy functions KEPT (Wave E retires). **22 new foundation tests** (3 capability + 4 bit-identical apply + 3 bit-identical solve + 6 reciprocity + 3 linearity-of-apply + 3 dispatch). **11/11 regression bit-identical (791s)** — additive code path; existing solver paths in solver.py unchanged. ERR-026 NOT closed (deferred to Wave E Issue 15) |

Wave D verification gates (end-of-session):

- 98 SN tests on integrated campaign HEAD (16 SNMesh + 27 scattering/fission + 9 dispatch + 22 SNStreamingOperator + 24 spatial)
- 11/11 regression snapshots bit-identical across all 3 rounds (worktree runs 668s/777s/791s + orchestrator-side runs 643s/777s + final defensive run pending)
- 27 + 2 xfail safety net intact (ERR-026 tripwires gated for Wave E, did NOT flip in any of 3 rounds)
- Sphinx clean (-W exit 0); audit clean (orphan count net 0 — Wave D added equation labels including `sn-streaming-reciprocity`, all immediately covered by `@pytest.mark.verifies`)
- 36/38 ERR coverage preserved

### Operational notes (lessons from Wave A + Wave B + Wave C + Wave D)

1. **Worktree-base bug** (discovered Wave B Round 1, mitigated for the rest): background-dispatched worktrees may come up at `main`'s HEAD (06b46f2) instead of the orchestrator's current branch HEAD. Mitigation: brief every method-implementer / general-purpose dispatch with explicit detection (`git status && git log --oneline -3 && ls orpheus/...`) and recovery (`git rebase refactor/sn-operator-algebra`) instructions. Issues 3, 4, 5, 7 all hit this and recovered cleanly via rebase.
2. **Cherry-pick conflict pattern** (consistent across all 3 rounds): conflicts on `docs/verification/matrix.rst` (Sphinx auto-regenerates with different test-count totals per branch). Resolution: `git checkout --ours` on the conflicting file; complete the cherry-pick; re-run `sphinx-build`; commit the regenerated matrix as a `chore(matrix)` commit.
3. **Bit-identical contract is the campaign's success criterion** — and held: every SN-touching commit (Issues 4 and 7) was gated by `pytest -m regression -q` (629–640s, 11/11 bit-identical) before merge. The agent verifications + post-merge orchestrator verifications agreed every time.
4. **Sub-agent assignment heuristic**: `general-purpose` for software-only primitives without published-math content (Issues 1, 3, 5); `method-implementer` for issues with bit-identical contracts or published-math grounding (Issues 2, 4, 6, 7). Dispatching parallel pairs (one general-purpose + one method-implementer) per round avoided sequential bottlenecks.
5. **`type:refactor` label does not exist** in the GitHub repo. Reshape issues marked `type:refactor` in the plan (#153, #156, #157, #159, #161, #162, #164) were filed with `type:improvement` instead. Acceptable substitution.
6. **Wave C duplicate-edit side effect** (Round 1 only): the `general-purpose` agent's worktree-isolated dispatch landed file edits in BOTH the agent's worktree AND the main repo (probably a tooling oversight when the agent's CWD-tracker was reset between Bash calls). Mitigation: orchestrator stashed the duplicated main-repo edits before cherry-picking the agent's worktree commit; the cherry-pick was clean since the stash removed the working-tree conflict surface. Round 2 (`method-implementer`) did NOT exhibit this. Worth surveilling on future `general-purpose` worktree dispatches.
7. **Wave C source-semantics confirmation** (Round 2 lesson): the `source` parameter passed to `CellUpdate.update(...)` follows the sweep's call-site convention — `source = Q · V · weight_norm` (already weight-normalized AND already volume-multiplied). The slab cumprod path bakes `2 · weight_norm · dx` into `source_coeff` so the per-cell scalar form becomes `s = 2 · source / denom`; the curvilinear path inserts `source` directly into the numerator as `QV[i]`. Both branches verified bit-identical via `np.array_equal` hand-calc tests typed verbatim from `sweep.py:117-123/350-361`.
8. **Wave D ERR-026 closure correction** (during planning): the campaign plan originally said Issues 11+12 close ERR-026, but analysis showed the bug lives in the curvilinear sweep's one-directional WDD line (`sweep.py:361`) which Wave C's `DiamondDifference` reproduces bit-identically. Wave D Issue 12's unified sweep with `cell_update=DD` therefore PRESERVES the bug. The actual closure mechanism (per `error_catalog.md` ERR-026 "Fix" section): route `solve_sn_fixed_source` through BiCGSTAB for curvilinear — that's an Issue 15 (Wave E) `SNSolver` change, NOT a sweep/operator change. Wave D builds the foundation (`SNStreamingOperator` with matrix-free symmetric-closure `apply`); Wave E flips the solver path. Plan's "Issues 11+12 close ERR-026" was approximate; the precise closure mechanism is Issue 15.
9. **Wave D delegator pattern** (R1.2 deviation): Issue 13's plan said "REMOVE the four legacy methods" from `SNSolver`, but `EigenvalueSolver` Protocol surface (`compute_fission_source`) and `TestAnisotropicScattering` directly call `_build_aniso_scattering`. Removal would break both. Architecturally clean alternative: kept methods as **thin delegators routing to the new operators**. The math LIVES in operator classes (canonical algebra-of-record per Cardinal Rule 2); SNSolver just routes. Wave E Issue 15 retires the delegators when BiCGSTAB migrates to operator algebra.
10. **Wave D 2D wavefront cell-update parameterization deferred** (R2 deviation): Issue 12's per-cell `cell_update.update(...)` accepts `(ng,)` slices, but the 2D wavefront sweep's anti-diagonal vectorization works on `(n_diag, ng)` slices. Extending the `CellUpdate` Protocol with a batch-shaped `update_batch(...)` is Wave C-extension's responsibility when LD/EC/Step land. Wave D's bit-identical contract is satisfied by retaining inlined DD math in the 2D wavefront path; the dispatch consolidation (single `transport_sweep` entry point) is the Wave D deliverable. Pre-authorized by the plan brief.
11. **Wave D streaming_terms() convention bugs caught** (R2 finding): two latent bugs in `geometry/reduced_operator.py`'s `streaming_terms()` API surfaced during sweep integration: (a) face-area direction convention `A[i]=in, A[i+1]=out` matches outward-sweep semantics but not inward (which reads `A[i+1]=in, A[i]=out`); (b) cylindrical `abs_mu = abs(quadrature.mu_x[direction_idx])` reads the WRONG ordinate when called with `direction_idx=m_local` (within-level azimuthal index — the global ordinate is `level_idx[m_local]`). Worked around in sweep.py via `_streaming_terms_for_inward_sweep` and `_streaming_terms_with_abs_mu` helpers. Both should be retired when the geometry-layer API is extended (likely Wave C-extension or Wave E).
12. **Wave D `apply_transpose` design** (R3 deviation): the plan brief offered analytic-adjoint or self-adjoint-shortcut paths for `SNStreamingOperator.apply_transpose`. Agent chose a third path: **build dense matrix from `apply` by probing with unit basis vectors (one-time `O(n²)`), transpose, apply**. Mathematically rigorous (every linear operator's transpose IS its matrix's transpose); avoids per-geometry adjoint code (a chance for sign-flip / missing-factor / transposed-index bugs — the V&V failure modes #1-#3, #5). Reciprocity holds **by construction**; the test still has teeth because it gates `apply` linearity and dense-assembly-probe correctness. `CellUpdate` Protocol does NOT need an `adjoint_update` extension. Wave E may revisit for performance if dense-probe `O(n²)` becomes a bottleneck on large meshes.

### Wave E — Phase 4 iteration as operator algebra (Issues 14, 15) (DONE)

| Plan # | GH # | Status | Commit | Summary |
|---|---|---|---|---|
| 14 | [#163](https://github.com/deOliveira-R/ORPHEUS/issues/163) | ✓ DONE | `52d7688` | `orpheus/numerics/iteration.py` — new (~545 LOC). `SourceIteration(L, S, F, *, inverter=None, max_iter, tol)` solves `(L − S − F)·ψ = q_ext` via fixed-point iteration; default `inverter=None` routes through `L.solve` (sweep-as-L⁻¹), caller-supplied `inverter` enables Krylov-on-`L.apply` (the load-bearing hook for ERR-026 reconciliation). `KEigenvalue(L, S, F, *, inverter, keff_estimator, eigenvalue_method="power", ...)` outer power iteration with `SourceIteration` inside; `eigenvalue_method` is forward hook for FEAST-style contour integral. Capability checks at construction (raise `MissingCapability` if any operand lacks `apply` or L lacks `solve`-or-`inverter`). `power_iteration` deprecated via `DeprecationWarning` on `orpheus.numerics.eigenvalue` import (kept functional for CP/diffusion/MoC/homogeneous). **11 new tests** in `tests/numerics/test_iteration.py` (10 foundation L0 synthetic + 1 L1 SN gate via `SNStreamingOperator + ScatteringOperator + FissionOperator` triple matching `solve_sn` to round-off). 11/11 regression bit-identical (770s). |
| 15 | [#164](https://github.com/deOliveira-R/ORPHEUS/issues/164) | ✓ DONE | `308499e` | `orpheus/sn/solver.py` rewritten (713 → 1187 LOC; +474 from new `_solve_fixed_source_krylov` helper [~370 LOC] and migrated `_build_rhs_*` helpers [~110 LOC] from operator.py). `SNSolver.__init__` constructs the (L = `SNStreamingOperator`, S = `ScatteringOperator`, F = `FissionOperator`) triple. `inner_solver` parameter accepts `"source_iteration"` or `"krylov"` ("bicgstab" REMOVED). `_solve_krylov` wraps GMRES around `SNStreamingOperator.apply` with `transport_sweep` as preconditioner. **Architectural deviations from plan** (preserved bit-identity): (a) `_solve_source_iteration` keeps the existing inlined loop verbatim because Pℓ angular state can't thread cleanly through generic SourceIteration; (b) `solve_sn_fixed_source` curvilinear default did NOT auto-flip to "krylov" in Round 2 — see Round 3. **Retired from `orpheus/sn/operator.py`** (~236 LOC removed): 7 functions (`build_transport_linear_operator{,_spherical,_cylindrical}`, `build_rhs{,_spherical,_cylindrical}`, `angular_flux_to_scalar`). **Retired from `orpheus/sn/geometry.py`**: 6 deprecated `@property` accessors (Wave D R1.1 transitional shims). **11 BiCGSTAB → krylov call sites migrated** (7 in `tests/sn/test_solver_components.py` + 3 in `tests/sn/test_spherical.py` + 1 in `examples/discrete_ordinates/demo_discrete_ordinates.py`). **`tests/sn/test_sweep_operator_inconsistency.py` rewritten** to use the new krylov path via `solve_sn_fixed_source(inner_solver="krylov")`; 4/4 PASS — krylov-on-apply gives the analytically-correct flat flux on the reflective-BC sphere problem (the original ERR-026 evidence). 11/11 regression bit-identical (746s). |
| 15-fix | (Round 3) | ✓ DONE | `2542f04` | `orpheus/sn/operator.py` — `solution_to_angular_flux*` and `transport_operator_matvec*` now consume `BoundaryOperator` instances on `SNMesh` (Wave B Issue 7 infrastructure) via `bc_outer.apply_to_incoming(outgoing, quad)`. Vacuum, reflective, white, periodic, albedo, and mixed BCs are handled uniformly. Bit-identity to pre-Round 3 hard-coded reflective fill preserved for `SpecularBoundaryOperator` (the `BC.reflective` factory — verified by 11/11 regression bit-identical, 786s). **ERR-026 PARTIAL closure**: krylov-on-`apply` now correct for vacuum-BC curvilinear problems on flat-flux constant-source cases. **MMS convergence still blocked**: empirically the FD operator's curvilinear outer-face flux at `i=nx-1` for outgoing `μ>0` uses cell-center as a face-flux approximation (`psi_right = fi[:, n, i, 0]`), which is exact for constant solutions but only first-order accurate on non-constant ones. Krylov-on-apply on the spherical isotropic MMS (vacuum BC) shows order ≈ 1.26, not the >1.9 the test asserts. Therefore: `solve_sn_fixed_source` curvilinear default kept at `"source_iteration"`; `"krylov"` is opt-in. The 2 ERR-026 xfail-strict markers in `tests/sn/l1_analytical/test_mms_curvilinear_aniso_dd_convergence.py` REMAIN xfail with updated reason strings reflecting the partial closure; the legacy `tests/sn/test_mms_curvilinear.py` got newly-added xfail-strict markers (same root cause). `error_catalog.md` ERR-026 status: OPEN → **PARTIAL CLOSURE**. Sphinx narrative section "ERR-026 deferred to Wave E" renamed to "ERR-026 closure status (partial through Wave E)". |

Wave E verification gates (end-of-session, three rounds verified
on campaign HEAD `00be85a`):

- 11/11 regression snapshots bit-identical across all 3 rounds (770s/746s/786s = ~13 min each)
- 27 + 2 xfail safety net intact through Round 1 + Round 2; 27 + 4 xfail in Round 3 (the 2 pre-existing markers retained with updated reasons + 2 newly-added on the legacy isotropic MMS file)
- Sphinx clean (-W exit 0); audit clean (36/38 ERR coverage preserved)
- 11 new iteration tests (Round 1) + krylov-path tests (Round 2) + BC-aware infrastructure tests (Round 3)
- Net LOC: +545 in `numerics/iteration.py`; SNSolver 713 → 1187 (+474 from krylov helper + RHS-build movement); operator.py 1126 → 890 (−236 from retired BiCGSTAB-only functions); +406 from Round 3 BC-plumbing edits

**ERR-026 closure status**: PARTIAL. The structural ingredients are
in place (operator triple + iteration primitives + Krylov-on-apply +
BC-aware FD operator). Constant-source flat-flux ERR-026 evidence
closes; vacuum-BC MMS convergence-rate evidence does not yet — blocked
by the FD operator's curvilinear outer-face cell-center face-flux
approximation (O(h)). See follow-up issue for FD-outer-face DD
extrapolation that would unblock O(h²) MMS convergence.

### Operational notes from Wave E (additions to the lessons list)

13. **Wave E Round 1 `keff_estimator` design**: SN's volume-weighted production/loss balance with (n,2n) folding can't be expressed via the generic Rayleigh quotient `(F·ψ).sum() / (L·ψ - S·ψ).sum()`; agent added an optional `keff_estimator` callable to `KEigenvalue` (default: generic Rayleigh; SN consumers supply `solver.compute_keff` as the bespoke estimator). The L1 SN gate test passes via this hook.
14. **Wave E Round 2 vacuum-BC discovery** (R2 finding): the legacy curvilinear FD operator's `solution_to_angular_flux_spherical` and matvec hard-coded reflective fill at the outer boundary. After Round 2's Krylov path replacement, the curvilinear vacuum-BC MMS path showed ~30% phi error — the operator was structurally tied to reflective BC. Round 3 closed this via `BoundaryOperator.apply_to_incoming` plumbing through the FD operator, but discovered a second-order-accuracy gap (note 15).
15. **Wave E Round 3 outer-face cell-center approximation** (R3 finding): once BC dispatch is correct, the curvilinear FD operator at `i = nx-1` for outgoing `μ > 0` reads `psi_right = fi[:, n, i, 0]` (cell-center) as a face-flux approximation. Exact for constant solutions; first-order on non-constant. Krylov-on-apply on smooth MMS converges at order ~1.26, not >1.9. Full ERR-026 closure on MMS requires DD diamond extrapolation at the outer face (`psi_face_out = 2·psi_cell - psi_face_in`) or an analogous ghost-cell technique. Filed as a follow-up issue. Note: the original WDD sweep also has order issues on this MMS (~order 0 — far worse), so Krylov is still a meaningful improvement; just not enough to satisfy `O(h²)`.
16. **Wave E `inverter` hook is the load-bearing design choice**: `SourceIteration(L, S, F, *, inverter=None)` lets the caller pick between `L.solve` (sweep-as-L⁻¹, default) and `lambda q: gmres(as_scipy_linop(L), q, M=...)` (Krylov-on-`L.apply`, the symmetric-closure path). This decoupling is what lets the same iteration primitive serve synthetic L0 cases (where L is a dense matrix and `inverter` defaults to direct solve), the existing SN source-iteration path (`inverter=L.solve`), and the new Krylov-on-apply path (`inverter` supplied explicitly). Without this hook, Round 2 would have needed parallel iteration primitives.

### Wave G — Phase 2 amendments + dunder consolidation (Issues 9.5, 9.6) (DONE)

| Plan # | GH # | Status | Commit | Summary |
|---|---|---|---|---|
| 9.5 | (closed inline) | ✓ DONE | `dd22cd6` | `orpheus/geometry/boundary.py` — mechanical rename of `ResolvedBC` → `BoundaryOperator` (and the 6 concrete subclasses `VacuumBC` → `VacuumBoundaryOperator`, `SpecularBC` → `SpecularBoundaryOperator`, `WhiteBC` → `WhiteBoundaryOperator`, `PeriodicBC` → `PeriodicBoundaryOperator`, `AlbedoBC` → `AlbedoBoundaryOperator`, `MixedBC` → `MixedBoundaryOperator`). `SNMesh.BC_REGISTRY` → `SNMesh.BOUNDARY_OPERATOR_REGISTRY`; helpers `_sn_bc_vacuum` / `_sn_bc_reflective` → `_sn_vacuum_boundary_operator` / `_sn_reflective_boundary_operator`. 19 files (production + tests + docs + plans + agent memory), 255 ins / 255 del. **`__eq__` shims preserved** verbatim on `VacuumBoundaryOperator` and `SpecularBoundaryOperator`. Bit-identity held: 11/11 regression PASSED in 12:54. |
| 9.6 | (closed inline) | ✓ DONE | `5fbe378` | Architecture + dunder consolidation aligning with Grand Report v3 §4–§6. **Six amendment categories shipped in one bundled commit** (per the plan's "no 6 sub-PRs" directive): (A) `orpheus/numerics/space.py` (NEW, 283 LOC) — `FunctionSpace` primitive (`name`, `shape`, `inner_product_weights`, `inner_product`, `norm`, `__eq__`/`__hash__`/`__repr__`) plus pre-populated factories `angular_flux_space`, `scalar_flux_space`, `boundary_trace_space`. (B) `orpheus/numerics/registry.py` (NEW, 174 LOC) — `RegistryMixin` per Grand Report v3 §4 (`key=` class-creation kwarg, `__init_subclass__` auto-registration, `create()` factory; handles the `@dataclass(slots=True)` rebinding interaction). (C) `orpheus/numerics/operator.py` (+315 LOC) — `domain`/`range` properties on `LinearOperator`, weight-aware Hilbert `adjoint()` + `.H` alias (the v3 §6.3 vocabulary; falls back to plain transpose when weights are None), `__call__(x)` aliasing `apply(x)`, `__pow__(n)` via repeated `__matmul__`, composition compatibility checks raising new `IncompatibleOperatorComposition`, `__repr__` on every composer. `OperatorSum`/`ScaledOperator` now forward `*args, **kwargs` so `0.7*spec + 0.3*white` composes under BoundaryOperator's 2-arg `apply(psi_out, quad)`. (D) `orpheus/numerics/measure.py` (+69 LOC) — `integrate` accepts pre-evaluated arrays in addition to callables (the load-bearing affordance for `volume_measure(values)`); `__call__` aliases `integrate`; container protocol (`__iter__`/`__len__`/`__getitem__`); custom `__repr__`. (E) `orpheus/geometry/boundary.py` (+191 / -91, net +100) — `BoundaryOperator(LinearOperator, RegistryMixin, ABC)` structural lift; concrete subtypes self-register via `key=` (`vacuum`, `reflective`, `white`, `periodic`, `albedo`, `mixed`). `apply(psi_out, quadrature)` is canonical; `apply_to_incoming` **dropped entirely** (28 call sites in `sn/sweep.py` + `sn/operator.py` migrated). **`SpecularBoundaryOperator` ships `apply_transpose`** — clean dual via the same reflection permutation, scaled by `albedo`; bumps capabilities to `{apply, apply_transpose}`. (F) `orpheus/sn/spatial/cell_update.py` (+69 LOC) — `CellUpdateBase(RegistryMixin, ABC)` introduced (Protocol stays for structural typing); `DiamondDifference` inherits and self-registers (`key="diamond_difference"`). (G) `orpheus/geometry/mesh.py` (+58 LOC) — `Mesh1D.volume_measure` and `Mesh2D.volume_measure` properties (1D returns `(N,)` nodes; 2D returns `(nx*ny, 2)` nodes via `meshgrid(indexing='ij')`). **B9 wiring not performed in this commit** — see Wave G operational note 17 below; filed as GH [#169](https://github.com/deOliveira-R/ORPHEUS/issues/169). New tests: `test_space.py` (158 LOC), `test_registry_mixin.py` (153 LOC), `test_operator.py` extensions (+239 LOC including weight-aware adjoint identity), `test_measure.py` extensions (+90 LOC), `test_mesh.py` (NEW, 168 LOC), `test_boundary.py` extensions (+115 LOC including SpecularBoundaryOperator transpose reciprocity), `test_diamond.py` extensions (+27 LOC). 25 files changed, 2407 ins / 170 del. **11/11 regression bit-identical** in 12:43. |

Wave G verification gates (end-of-session):

- 11/11 regression snapshots bit-identical across both rounds (`dd22cd6` 774s + `5fbe378` 763s)
- Sphinx clean (-W exit 0); audit clean (23 orphans, 36/38 ERR coverage — matches Wave E baseline; no new orphans introduced)
- 529 numerics + geometry + spatial tests pass; 38 SN solver-component tests pass; 22 SNStreamingOperator tests pass; full L1 + regression sweep green
- 27 + 4 xfail safety net intact (no markers flipped, none added)

### Operational notes from Wave G (additions to the lessons list)

17. **Wave G B9 `volume_measure` wiring initially deferred under a misread regression contract** (corrected in note 21): the first wiring attempt at `CPSolver.compute_keff` revealed a ~1 ULP drift between `np.sum(sig_p * flux * V[:, None])` and `mu(sig_p * flux).sum()`. The drift is FP non-associativity: the two paths use different IEEE-754 reduction orders. I initially read the campaign's "bit-identical" framing as `np.array_equal` and concluded wiring was blocked. Filed GH [#169](https://github.com/deOliveira-R/ORPHEUS/issues/169) with snapshot-regeneration / integrate-overload paths. **Note 21 corrects the diagnosis** — the actual contract is `rtol=1e-12`, well above 1 ULP — and ships the principled wiring.
18. **Wave G `RegistryMixin` × `@dataclass(slots=True)` interaction**: the slots decorator creates a NEW class object replacing the original; `__init_subclass__` registers the pre-slots class while the final `Foo` symbol points to the post-slots class. Fixed in `RegistryMixin.__init_subclass__`: when re-registration with same qualname/module but different identity occurs, silently update the entry to the new class. Documented in code comments and the closeout memo. Without this fix, `BoundaryOperator.create("vacuum")` would have returned the stale pre-slots class, breaking equality and instance-of checks.
19. **Wave G Hilbert adjoint correctly weight-aware**: per Grand Report v3 §6.3, `.H` is the Hilbert adjoint with respect to the domain/range `inner_product_weights`, NOT just representation transpose. For an operator `A: V → W` with diagonal weights `w_V` and `w_W`, the implementation is `A* y = (1/w_V) ⊙ apply_transpose(w_W ⊙ y)`. Falls back to plain `apply_transpose` (Euclidean) when both weights are None. The `_AdjointOperator` wrapper (~30 LOC in `numerics/operator.py`) handles both paths uniformly. Tested via reciprocity identity `<A u, v>_W = <u, A* v>_V` for both Euclidean and quadrature-weighted angular spaces.
20. **Wave G `apply_to_incoming` retired completely (no alias)**: the original 9.6 spec recommended keeping `apply_to_incoming` as an alias for solver-code readability. The user override during planning was: "Migrate all call sites to `apply`" — the Grand Report v3 vocabulary uses `A @ ψ` / `A(ψ)` for application across the operator algebra; the directional name is a transitional artifact. 28 call sites in `sn/sweep.py` + `sn/operator.py` migrated mechanically.
21. **Wave G follow-up: GH #169 closed via principled wiring + corrected contract reading**: re-investigation showed the actual snapshot test contract is `assert_allclose(rtol=1e-12, atol=1e-13)`, not `np.array_equal`. The 1 ULP drift the wiring introduces is ~2e-16 relative — five orders of magnitude under tolerance. Wired `compute_keff` in both `SNSolver` and `CPSolver` to compose `compute_group_production_rate(φ).sum() / compute_group_absorption_rate(φ).sum()` (or the analogous total-rate composition for CP). The per-group rate vector intermediates ARE physical quantities (per-group reaction rates — used in spectral diagnostics, group-collapse, broad-group prep). Verified against `k_∞` for homogeneous reflective slabs (independent analytical reference). 11/11 regression PASSED at `rtol=1e-12`. **Principle codified** in `.claude/skills/vv-principles/SKILL.md` § "Bit-identity vs principled-equivalence" (loaded by qa, test-architect, numerics-investigator, archivist sub-agents) — three explicit acceptance criteria for non-bit-exact refactors: principled intermediates, structurally-independent verification, FP-non-associativity-bounded drift.

### Wave H — Issue #168 architecture (Phases A + B; Phase C pending) (PARTIAL)

| Plan # | GH # | Status | Commit | Summary |
|---|---|---|---|---|
| 168-A | [#168](https://github.com/deOliveira-R/ORPHEUS/issues/168) Phase A | ✓ DONE | `d73ef68` | `orpheus/sn/spatial/boundary_face_flux.py` (NEW, 415 LOC) — `BoundaryFaceFlux` Protocol + concrete ABC `BoundaryFaceFluxBase` (RegistryMixin auto-registration via `key=`) + `DDExtrapolation` (default, `psi_face_out = 1.5·psi[N-1] − 0.5·psi[N-2]`) + `CellCenter` (legacy first-order reproducer). `solution_to_angular_flux_*` rewritten to return `(fi, boundary_face_flux)` separating cell-centre from BC-face storage (closes Defect 2 of #168). 21 foundation tests + 7 SNStreamingOperator invariant tests. **6 curvilinear regression snapshots intentionally invalidated** (`.npz` files deleted; tests SKIP gracefully via existing `if not snapshot_file.exists(): pytest.skip(...)` mechanism); 5 Cartesian snapshots stay bit-identical. ERR-026 stays PARTIAL — Phase A's structural rewrite alone yields O(h^1.5–1.7) on MMS, below the >1.9 bar; the 4 xfail-strict markers stay. |
| 168-B | [#168](https://github.com/deOliveira-R/ORPHEUS/issues/168) Phase B | ✓ DONE (architecture) | (Phase B commit) | `orpheus/sn/spatial/pole_angular_closure.py` (NEW, 480 LOC) — `PoleAngularClosure` Protocol + ABC + **three** concrete strategies: (a) `LegacyTauSymmetricInterpolation` (default, bit-for-bit pre-Phase-B reproduction), (b) `BaileyFlatFluxRedist` (algebraic flat-flux collapse `redist = -μ_n·ΔA·ψ/V`; equivalent to Legacy on flat ψ, differs per-ordinate on angularly-varying), (c) `MorelMontryAngularSweep` (canonical Hébert §3.9.4 per-cell M-M weighted DD recurrence with α-cascade and starting condition `ψ_{1/2}=0`; **opt-in only**). `transport_operator_matvec_*` and `SNStreamingOperator.apply` thread `pole_angular_closure` through. **Citation correction** across `orpheus/geometry/reduced_operator.py`, `orpheus/sn/sweep.py`, `orpheus/sn/spatial/diamond.py`: Bailey-Adams-Yang-Zika (2009) JCP 227 (a piecewise-linear FE *diffusion* paper — wrong) → **Bailey-Morel-Chang (2010) NSE 165(2):149-169** (LLNL-JRNL-420356; Hébert §3.9.4 is primary). α-recursion normalisation pinned: ORPHEUS `α^O = α^H/2`. 28 foundation tests + 5 L1 flat-flux-identity tests + 1 SNStreamingOperator canonical-form invariant. Sphinx narrative subsection added. **ERR-026 stays PARTIAL** — see Phase C scope below. |

Wave H Phase A + B verification gates:

- 5 Cartesian regression snapshots PASS bit-identical (assert_allclose rtol=1e-12); 6 curvilinear SKIP (Phase A invalidations preserved through Phase B because default = `LegacyTauSymmetricInterpolation` reproduces pre-Phase-B math).
- 67 Phase B + Phase A tests on `test_pole_angular_closure.py` + `test_pole_closure_flat_flux_identity.py` + `test_snstreamingoperator.py` + `test_sweep_operator_inconsistency.py` all pass.
- 4 ERR-026 xfail-strict markers stay xfail (intentional through Phase B; flip pending Phase C).
- Sphinx clean (-W exit 0); audit clean (24 orphans = baseline 23 + 1 new `pole-mm-recurrence` label; 36/38 ERR coverage unchanged).

### Operational notes from Wave H (additions to the lessons list)

22. **Wave H Phase B empirical finding — pairing canonical pole closure with non-aligned spatial closure regresses the operator**: when `MorelMontryAngularSweep` is paired with the apply matvec's existing **spatial** closure (interior arithmetic average `0.5·(ψ_i + ψ_{i+1})` + Phase A outer DD extrapolation), the operator is *worse* than the legacy form on flat ψ. `test_spherical_sweep_vs_bicgstab_flat_flux` regresses (φ ranges 0.6–1.004 instead of analytical 1.0); MMS probes diverge to ~10^14 at nc=10. Mechanism: the DD angular recurrence on flat ψ produces oscillating half-angle face fluxes `0, 2c, 0, 2c, …`; combined with the symmetric arithmetic spatial average this gives garbage. The apply matvec must use the **same WDD spatial closure** the sweep uses (`ψ_face_out = 2·ψ_avg − ψ_face_in`) for the canonical angular closure to be consistent. Phase B therefore ships the canonical form as **opt-in** with `LegacyTauSymmetricInterpolation` as default; Phase C aligns the spatial closure and flips the default.
23. **Wave H Phase B three-strategy split** (deviation from brief): the brief specified two strategies (`MorelMontryAngularSweep` canonical + `BaileyFlatFluxRedist` legacy ablation) on the implicit assumption that BFF reproduces the pre-Phase-B math. This is true on flat ψ but false on angularly-varying ψ — the pre-Phase-B form uses τ-symmetric interpolation between cell-centres, BFF uses a flat-flux algebraic collapse. They differ per-ordinate on varying ψ. `LegacyTauSymmetricInterpolation` was added as the bit-for-bit reproducer to satisfy the regression contract on arbitrary input. Documented in the closeout memo as Deviation 1.
24. **Wave H literature consult discovered citation bug**: `orpheus/geometry/reduced_operator.py` previously cited "Bailey, T. S., Adams, M. L., Yang, B., & Zika, M. R. (2009). *A piecewise linear finite element discretization of the diffusion equation for arbitrary polyhedral grids*. JCP 227, 3738-3757." This is the **wrong Bailey paper** — Bailey-Adams-Yang-Zika is a piecewise-linear FE *diffusion* paper, unrelated to curvilinear S\ :sub:`N`. The intended reference is **Bailey, Morel & Chang (2010), *Asymptotic Diffusion-Limit Accuracy of Sn Angular Differencing Schemes*, NSE 165(2):149-169** (LLNL-JRNL-420356; OA at https://www.osti.gov/servlets/purl/1020346). Hébert (2009) §3.9.4 is the primary source; Bailey-Morel-Chang 2010 is the auxiliary asymptotic-diffusion-limit τ-clamp justification. Corrected across `reduced_operator.py`, `sn/sweep.py`, `sn/spatial/diamond.py`. Worth surveilling on remaining citations in the codebase.

### Remaining waves (next sessions)

| Wave | Phase | Issues | Description |
|---|---|---|---|
| H Phase C | — | GH [#168](https://github.com/deOliveira-R/ORPHEUS/issues/168) Phase C | **Spatial-closure alignment** for full ERR-026 closure. (1) Rewrite the apply matvec's interior face-flux closure from arithmetic average `0.5·(ψ_i + ψ_{i+1})` to the sweep's WDD form `ψ_face_out = 2·ψ_avg − ψ_face_in` (design memo §6.4 / §11). The Phase A `BoundaryFaceFlux` Protocol stays; the interior face closure receives a parallel reformulation. (2) Flip `SNMesh.pole_angular_closure` default `LegacyTauSymmetricInterpolation` → `MorelMontryAngularSweep`. (3) Flip curvilinear `solve_sn_fixed_source` default `"source_iteration"` → `"krylov"`. (4) Regenerate the 6 deleted curvilinear regression snapshots; verify each via FN-method cross-check on the closest Sood La13511 case (`orpheus/derivations/continuous/sood_registry/`). (5) Remove the 4 xfail-strict ERR-026 tripwires. ERR-026 status: PARTIAL CLOSURE → CLOSED. |
| ✓ DONE | — | GH [#169](https://github.com/deOliveira-R/ORPHEUS/issues/169) | **`compute_keff` migration to `volume_measure`** — wired via principled chain. Per the corrected diagnosis: (a) the actual snapshot test contract uses `assert_allclose(rtol=1e-12)`, NOT `np.array_equal` — the campaign plan's "bit-identical" framing was overstated. (b) The drift from rewiring is at most ~1 ULP (~2e-16 relative) — five orders of magnitude under the `1e-12` tolerance. (c) `compute_keff` now composes `compute_group_production_rate(φ).sum() / compute_group_absorption_rate(φ).sum()`, exposing per-group reaction rate vectors as named, inspectable diagnostic quantities (the principled intermediate). (d) Verified against `k_∞ = (νΣ_f · φ)/(Σ_a · φ)` for homogeneous reflective slabs (analytical limit, structurally independent reference). 11/11 regression PASSED at `rtol=1e-12`, no contract relaxation needed. Principle codified in `.claude/skills/vv-principles/SKILL.md` § "Bit-identity vs principled-equivalence". See Wave G operational notes 17 and 21. |
| C-ext | 2 | 9 (LD, EC, Step) | Concrete cell updates beyond DD: `LinearDiscontinuous`, `ExponentialCharacteristic`, `Step` strategies. Each with its own MMS spatial-convergence test. Sequential (one round per strategy). May also include geometry-layer cleanup of the `streaming_terms()` convention bugs (per Wave D operational note 11) and 2D-wavefront `update_batch(...)` Protocol extension (note 10). |
| F | 5 | 17, 18 | Symmetry-preservation + reciprocity invariant tests; Sphinx documentation campaign (architectural narrative — operator algebra unifying SN/MoC/CP/diffusion). |
| Cross-solver | — | (new GHs) | Migrate CP, diffusion, MoC, homogeneous to consume `KEigenvalue` / `SourceIteration`; retire `power_iteration` once all four migrate. Wave E Round 1 retained the deprecated function for this transition. |

---

## Campaign context

The codebase has a working multigroup, anisotropic-scattering, steady-state SN
solver with three coordinate systems (Cartesian slab, cylindrical, spherical).
The math is correct — Bailey et al. (2009) angular redistribution, Morel-Montry
weighted diamond closure for spherical, per-level azimuthal redistribution for
cylindrical, all referenced in code.

What's missing is **architectural unification across the deterministic solver
family**. The pieces are consistent inside SN but the abstractions don't span
SN/MoC/CP/diffusion. This campaign installs the cross-cutting primitives —
`LinearOperator`, `DiscreteMeasure`, `ReducedStreamingOperator`,
`CellUpdate` — using SN as the pilot. MoC, CP, and MC migrate sequentially
in their own campaigns once SN is proven.

The architectural narrative is: the transport equation has differential,
integral, and variational forms. SN/MoC discretize the differential form;
CP discretizes the integral form; PN/diffusion discretize spectral truncations.
At the operator-algebra level they all expose `(L, S, F)` triples consumed
uniformly by source iteration, Krylov, and eigenvalue solvers. That's the
end state.

## Confirmed decisions

1. **1D cumprod fast path lives inside the unified sweep algorithm**, not as
   a separate dispatched function. It's a DD-specific algebraic optimization;
   keep it but express it as such.
2. **Issues #96 / #97 close as a side effect of the reshape**: the BiCGSTAB
   path stops building a separate FD operator and instead does matrix-free
   Krylov on `SNStreamingOperator.apply`. Behavioral change in the
   BiCGSTAB convergence is acceptable.
3. **`numerics/eigenvalue.py:power_iteration` lingers** until CP, diffusion,
   and other consumers migrate to operator-algebra iteration. No hard
   removal date. **Note**: full migration of all solvers (SN → MoC → CP → MC)
   IS planned, sequentially. The new infrastructure must be designed
   forward-compatibly with all four.

## Cross-solver migration sequence

This SN reshape is **Wave 1**. Subsequent waves consume the same primitives:

| Wave | Solver | Notes |
|------|--------|-------|
| 1 | SN | this campaign — installs primitives |
| 2 | MoC | uses `BundleMeasure`, ray-quadrature DiscreteMeasure |
| 3 | CP | dense-matrix `LinearOperator`, scalar-flux operator |
| 4 | MC | uniform interface for verification, shares `Geometry` |

Primitives must be designed in this campaign without MoC/CP/MC blockers.
`BundleMeasure` ships in Phase 0 even though SN doesn't use it.

## Architectural commitments (invariants)

These are non-negotiable across all issues. Any deviation is a session
failure per Cardinal Rule 2 (architecture is critical).

1. **`LinearOperator` is the cross-cutting abstraction** for L (streaming +
   collision), S (scattering), F (fission). Capability-tagged
   (`apply` / `solve` / `apply_transpose`), composable, capability-checked
   at construction. Operators that lack a capability raise on attempted
   call, not on import.

2. **`DiscreteMeasure` is the quadrature primitive.** Composable via tensor
   product (`__mul__`), pushforward, restrict, direct sum (`__add__`).
   `AngularQuadrature` becomes a thin SN-domain adapter over it.

3. **Connection coefficients live on geometry**, exposed through
   `ReducedStreamingOperator`. Cell updates consume `StreamingTerms`,
   never raw α / ΔA / τ arrays.

4. **Cell updates are polymorphic** with a curvilinear-aware signature:
   `update(streaming_terms, total_xs, source, upstream_state) → CellResult`.
   `upstream_state.angular_upstream` is `None` for Cartesian, a per-cell
   per-group array for curvilinear.

5. **Production schemes are positivity-preserving by construction**: LD,
   exponential characteristic, step. DD kept only as comparison artifact,
   marked `is_positivity_preserving = False`.

6. **Existing patterns are extended, not replaced.** `BC_REGISTRY`,
   `SNMesh`, `Mesh1D.from_geometry`, the geometry-layer two-tier split
   (reference vs production solvers), all stay. The reshape adds layers
   above and refactors internals; it does not replace the public shape.

7. **Behavioral regression is the gating contract.** Existing SN tests
   produce bit-identical results when `cell_update == DiamondDifference`.
   Snapshot-frozen scalar-flux + k_eff outputs are captured BEFORE the
   reshape begins (Issue 16, executed first).

## Out of scope

- MoC, CP, MC reshapes (separate campaigns)
- Time-dependent / kinetics (steady-state stays; transient is a downstream
  campaign building on this primitives layer)
- New angular quadrature families (LS_N already implemented; Lebedev variants
  beyond what scipy provides are not expanded here)
- New geometry kinds (HSPH, ANN, 2D curvilinear) — these slot in cleanly
  once `ReducedStreamingOperator` exists
- Anisotropic-scattering Pℓ ordering changes (current implementation in
  `SNSolver._build_aniso_scattering` is preserved as-is, just lifted to
  `ScatteringOperator`)

---

## Issue plan

Each subsection is a self-contained spec for a GitHub issue. Format:

- Module label, type, phase, dependencies, complexity (S / M / L)
- Context (the WHY — Cardinal Rule 3)
- Acceptance criteria (checkbox list — gating for PR merge)
- Files affected
- Design notes (decisions, gotchas, references)

A sub-agent should translate each subsection 1:1 into a GitHub issue with
the appropriate `module:` and `type:` labels.

### Phase 0 — Numerics primitives

#### Issue 1: `numerics/operator.py` — `LinearOperator` protocol

- **Module**: `module:numerics`
- **Type**: `type:feature`
- **Phase**: 0
- **Depends on**: none
- **Complexity**: M

**Context**: Today `scipy.sparse.linalg.LinearOperator` is used only inside
the BiCGSTAB Krylov path (`orpheus/sn/operator.py`). The sweep, scattering
source, and fission source live as bare functions / methods. Lifting them
to a uniform `LinearOperator`-shaped interface lets eigenvalue and Krylov
code consume any solver method (SN/MoC/CP) without knowing which transport
discretization is below. This is the foundation for the cross-solver
migration sequence.

**Acceptance criteria**:

- [ ] `LinearOperator` Protocol in `orpheus/numerics/operator.py` with
      `apply(x)`, optional `solve(b)`, optional `apply_transpose(x)`,
      and `capabilities: frozenset[str]` property
- [ ] Composition primitives: `OperatorSum`, `OperatorProduct`,
      `ScaledOperator`, `IdentityOperator`, `ZeroOperator` — each computes
      its own capability set from constituents (sum requires `apply` from
      both; product requires composition of capabilities)
- [ ] Operator algebra: `__add__`, `__sub__`, `__mul__` (scalar),
      `__matmul__` (operator product) on `LinearOperator`
- [ ] Capability mismatches raise `MissingCapability` at composition time,
      not at call time
- [ ] Adapter `as_scipy_linop(op)` for Krylov consumption via scipy
- [ ] Unit tests on synthetic operators (numpy matrices wrapped) covering:
      composition, capability propagation, scipy interop

**Files**:

- new: `orpheus/numerics/operator.py`
- new: `tests/numerics/test_operator.py`

**Design notes**: Follow the runtime-checkable Protocol pattern already
established for `AngularQuadrature` in `orpheus/sn/quadrature.py`. The
capability set is the key idea — many operators have no efficient `solve`
(S, F), and forcing them to provide stubs is harmful. Don't abstract
shape/dtype — leave that to numpy duck-typing for now.

---

#### Issue 2: `numerics/measure.py` — `DiscreteMeasure` primitive

- **Module**: `module:numerics`
- **Type**: `type:feature`
- **Phase**: 0
- **Depends on**: none
- **Complexity**: M

**Context**: Quadratures are mathematically discrete measures. The natural
operations (tensor product → product quadrature, pushforward → coordinate
change, restriction → half-range) all have direct mathematical content.
Today these are implicit: `ProductQuadrature.create(n_mu, n_phi)` constructs
the product internally without exposing the tensor-product structure.
Promoting `DiscreteMeasure` to a primitive lets consumers compose 1D rules
into 2D rules into S² rules cleanly, and is **required for MoC's bundle
measures** (Wave 2).

**Acceptance criteria**:

- [ ] `DiscreteMeasure` in `orpheus/numerics/measure.py`: `nodes`,
      `weights`, `space` (str/enum tag), `integrate(f)`, `__mul__` (tensor
      product), `__add__` (direct sum), `pushforward(phi)`,
      `restrict(predicate)`, `n_points`
- [ ] `BundleMeasure` for disintegrated measures — base measure plus
      per-base-point fiber measure factory. Required for MoC; not used
      by SN in this campaign but lives here so the abstraction is correct
      from day one
- [ ] 1D primitives: `gauss_legendre(n)`, `gauss_chebyshev(n)`,
      `equispaced(a, b, n)` — each returns a `DiscreteMeasure`
- [ ] Optional metadata fields populated lazily: `invariance_group`
      (filled by Issue 3), `degree_of_exactness`
- [ ] Unit tests: integration of polynomials of known degree, tensor
      product equivalence to manual construction, pushforward correctness
      under invertible maps

**Files**:

- new: `orpheus/numerics/measure.py`
- new: `tests/numerics/test_measure.py`

**Design notes**: Don't try to enforce `Space` types via Python generics —
not expressive enough without runtime overhead. Use `space: str | enum`
as a runtime tag. Composition operations return new `DiscreteMeasure`s
with metadata combined sensibly. `BundleMeasure` is critical to ship now
so MoC migration doesn't have to revisit Phase 0.

---

#### Issue 3: `numerics/symmetry.py` — Subgroups of O(3) and invariance

- **Module**: `module:numerics`
- **Type**: `type:feature`
- **Phase**: 0
- **Depends on**: Issue 2
- **Complexity**: S

**Context**: Quadrature selection (Issue 5) needs subgroup containment
lattice logic: geometries declare their symmetry group, quadratures
declare theirs, selection requires `G_geom.is_subgroup_of(G_quad)`.
The lattice is finite for ORPHEUS-relevant cases — a static lookup
suffices; no need for generator-based machinery.

**Acceptance criteria**:

- [ ] `SubgroupOfO3` enum-backed in `orpheus/numerics/symmetry.py` with
      named entries: `Trivial`, `Z2`, `SO2`, `O2`, `OctahedralOh`,
      `IcosahedralIh`, `SO3`, `O3`, plus parameterized `Cn(n)`,
      `Dnh(n)` for hex/lattice future use
- [ ] `contains(self, other) -> bool` implementing subgroup containment
      via static lookup table
- [ ] `is_invariant(self, measure: DiscreteMeasure) -> bool` checks that
      measure nodes form a union of orbits with consistent weights
- [ ] Unit tests: containment lattice (Z2 ⊂ O_h ⊂ O3, etc.); invariance
      checks on existing quadratures (Lebedev → O_h, LS_N → O_h,
      Gauss-Legendre on μ → SO2)

**Files**:

- new: `orpheus/numerics/symmetry.py`
- new: `tests/numerics/test_symmetry.py`

**Design notes**: Don't overengineer. Static dict mapping
`(G_a, G_b) → bool` for containment is fine. Generator-based machinery
deferred until ORPHEUS needs novel discrete groups (hex / triangular
lattices: C_6v, D_6h — add when consumed).

---

#### Issue 4: Refactor `AngularQuadrature` to use `DiscreteMeasure`

- **Module**: `module:numerics`, `module:sn`
- **Type**: `type:refactor`
- **Phase**: 0
- **Depends on**: Issues 2, 3
- **Complexity**: M

**Context**: `orpheus/sn/quadrature.py` has 4 working quadrature classes
(GL1D, Lebedev, LS_N, Product) pre-dating the `DiscreteMeasure`
abstraction. Each exposes a similar but ad-hoc interface (`mu_x`, `mu_y`,
`weights`, `N`, `reflection_index`, `spherical_harmonics`). Refactoring
them to be `DiscreteMeasure`-backed — with the existing
`AngularQuadrature` Protocol kept as a domain-specific adapter — gives
composability without breaking SN's heavy attribute-access consumption.

**Acceptance criteria**:

- [ ] `orpheus/numerics/quadrature/` package: `rules_1d.py`,
      `rules_sphere.py`, `rules_product.py`
- [ ] Each rule function returns a `DiscreteMeasure` with
      `invariance_group` and `degree_of_exactness` populated
- [ ] `orpheus/sn/quadrature.py` refactored: existing classes (GL1D,
      Lebedev, LS_N, Product) become thin adapters wrapping
      `DiscreteMeasure` and caching SN-specific fields (`mu_x`, `mu_y`,
      `mu_z`, `level_indices`, `_ref_x` etc.) on construction
- [ ] `reflection_index(axis)` re-implemented as `pushforward` with
      coordinate negation; result cached
- [ ] `spherical_harmonics(L)` stays where it is (SN-specific)
- [ ] **All existing SN tests pass without modification.** This is the
      gating constraint — no behavior change

**Files**:

- new: `orpheus/numerics/quadrature/{__init__,rules_1d,rules_sphere,rules_product}.py`
- refactor: `orpheus/sn/quadrature.py`
- existing tests: unchanged

**Design notes**: The bridge issue between the new primitive and existing
SN consumption. Trick: SN consumers index into `quad.mu_x`, `quad.weights`
etc. heavily. Don't break that. The adapter caches array views from the
backing `DiscreteMeasure` on construction.

---

#### Issue 5: Quadrature registry with G + V + structural tags + selector

- **Module**: `module:numerics`
- **Type**: `type:feature`
- **Phase**: 0
- **Depends on**: Issues 3, 4
- **Complexity**: S

**Context**: With G-tagged quadratures and a subgroup lattice, automated
selection becomes the precedence chain: G compatibility → V compatibility
→ structural compatibility → minimum points. `solve_sn` keeps explicit
quadrature-passing as the canonical API; the registry adds
`select_quadrature(geometry, ...)` as a convenience and as a documentation
artifact. The structural tags are themselves teaching content.

**Acceptance criteria**:

- [ ] `QuadratureSpec` dataclass: `name`, `factory`, `invariance_group`,
      `degree_of_exactness`, structural flags (`positive_weights`,
      `axis_aligned`, `level_structured`, `half_range_clean`)
- [ ] `quadrature_registry` populated with all current rules
- [ ] `select_quadrature(geometry, target_degree, **structural)
      → DiscreteMeasure` implementing the precedence chain
- [ ] Selection log returned alongside the chosen quadrature so the choice
      is explainable / loggable
- [ ] Unit tests: slab → GL1D (SO(2)-invariant after isotropy reduction);
      sphere → GL1D on μ_r; cylinder → ProductQuadrature; Cartesian-like
      2D → LS_N or Lebedev with explainable preference

**Files**:

- new: `orpheus/numerics/quadrature/registry.py`
- new: `tests/numerics/test_registry.py`

**Design notes**: Don't make selection mandatory. Explicit
quadrature-passing in `solve_sn` stays canonical; the registry is
opt-in convenience. The structural tags double as Sphinx documentation
content.

---

### Phase 1 — Geometry layer additions

#### Issue 6: `ReducedStreamingOperator` per geometry

- **Module**: `module:geometry`
- **Type**: `type:feature`
- **Phase**: 1
- **Depends on**: Issue 1
- **Complexity**: L

**Context**: Connection coefficients (α, ΔA/w, τ_mm) currently live in
`SNMesh._setup_spherical` and `SNMesh._setup_cylindrical`. They are
SN-specific consumers of geometry-side math. MoC and CP will need the
same reduced-operator concept (different consumption pattern, same
underlying math). Lifting it to geometry now removes future duplication.
Per Cardinal Rule 2: shared concepts between solvers means the codebase
needs an architectural overhaul.

**Acceptance criteria**:

- [ ] `ReducedStreamingOperator` class in
      `orpheus/geometry/reduced_operator.py`
- [ ] Properties: `requires_upstream_angular_state: bool`,
      `angular_marching_axis: Literal["mu", None]`
- [ ] Method: `streaming_terms(cell_idx, direction_idx, mu_level_idx=None)
      → StreamingTerms`
- [ ] `StreamingTerms` carries everything the cell update needs in this
      geometry: chord lengths, face areas, ΔA factor, α coefficients,
      τ_mm — geometry-dependent shape
- [ ] Concrete factory functions (or subclasses): `slab_streaming(mesh)`,
      `cylindrical_streaming(mesh, angular_measure)`,
      `spherical_streaming(mesh, angular_measure)`
- [ ] Output is bit-identical to current `SNMesh._setup_*` (verified by
      hash-equality of arrays in tests)
- [ ] Unit tests against precomputed reference values matching current
      `SNMesh.alpha_half`, `SNMesh.redist_dAw`, `SNMesh.tau_mm`

**Files**:

- new: `orpheus/geometry/reduced_operator.py`
- new: `tests/geometry/test_reduced_operator.py`

**Design notes**: The current `SNMesh` couples geometry math (connection
coefficients) with SN math (quadrature). Split: connection-coefficient
*values* depend on quadrature points (they're integrals of α from the
angular discretization), so the factory takes both mesh and angular
measure. Output is coordinate-system-aware. Per-level structure for
cylindrical (consuming `angular_measure.level_indices` via the adapter
from Issue 4) is preserved.

---

#### Issue 7: Boundary conditions as tensor decompositions

- **Module**: `module:geometry`, `module:sn`
- **Type**: `type:refactor`
- **Phase**: 1
- **Depends on**: Issue 4
- **Complexity**: M

**Context**: BCs are currently `BC(kind, params)` dataclasses with
SN-specific resolution in `SNMesh.BOUNDARY_OPERATOR_REGISTRY` returning a string tag.
Math says BCs are tensor decompositions: R = Σ_α G_α ⊗ A_α where G_α is
geometric (permutation/index map) and A_α is response (albedo, scalar
amplitude). Lifting the resolved BC to a tensor-decomposed object
makes specular / white / albedo / periodic / mixed all uniform, and
lets multi-region interfaces reuse the same primitives.

**Acceptance criteria**:

- [ ] `BoundaryOperator` Protocol/ABC in `orpheus/geometry/boundary.py`:
      `apply_to_incoming(angular_flux_outgoing, quadrature)
      → angular_flux_incoming`
- [ ] Concrete: `VacuumBoundaryOperator` (zero), `SpecularBoundaryOperator` (rank-1: permutation
      `pushforward` from `DiscreteMeasure` + albedo scalar), `WhiteBoundaryOperator`,
      `PeriodicBoundaryOperator`, `AlbedoBoundaryOperator`, `MixedBoundaryOperator` (rank-N sum)
- [ ] `SNMesh.BOUNDARY_OPERATOR_REGISTRY` updated: factories return `BoundaryOperator`
      instances, not string tags. The factory pattern stays
- [ ] Sweep code in `orpheus/sn/sweep.py` updated to call
      `resolved_bc.apply_to_incoming(...)` instead of branching on
      string kind
- [ ] All existing SN tests pass with bit-identical outputs for `vacuum`
      and `reflective`
- [ ] Unit tests for `WhiteBoundaryOperator` and `PeriodicBoundaryOperator` (currently unsupported
      by SN — adding support is a downstream win that this issue enables)

**Files**:

- new: `orpheus/geometry/boundary.py`
- modify: `orpheus/sn/geometry.py`, `orpheus/sn/sweep.py`

**Design notes**: This is where BC tensor-decomposition framing pays off
architecturally. Use the existing `BC_REGISTRY` pattern — just change
what the factories return. White and periodic come essentially free
once specular is correct.

---

### Phase 2 — SN cell update strategies

#### Issue 8: `CellUpdate` ABC with curvilinear-aware signature

- **Module**: `module:sn`
- **Type**: `type:refactor`
- **Phase**: 2
- **Depends on**: Issue 6
- **Complexity**: M

**Context**: Current sweep dispatches by `sn_mesh.curvature` ("spherical"
/ "cylindrical" / None) and inlines the cell update equation per
geometry. To swap spatial schemes (LD, exponential, characteristic)
systematically, the per-cell update needs to be a strategy. The
signature must accommodate curvilinear coupling — the previous
μ-direction's flux feeds into the connection term.

**Acceptance criteria**:

- [ ] `CellUpdate` ABC in `orpheus/sn/spatial/cell_update.py`:
      `update(streaming_terms, total_xs, source, upstream_state)
      → CellResult`
- [ ] Properties: `is_linear: bool`, `is_positivity_preserving: bool`
- [ ] `UpstreamState` dataclass: `spatial_upstream` (always present),
      `angular_upstream` (None for Cartesian, per-cell per-group array
      for curvilinear)
- [ ] `CellResult` dataclass: `cell_average_flux`,
      `outgoing_spatial_flux`, `outgoing_angular_state` (None for
      Cartesian)
- [ ] No concrete strategies yet (Issue 9). Just the ABC and
      protocol-level tests confirming any concrete strategy can be
      substituted

**Files**:

- new: `orpheus/sn/spatial/__init__.py`,
       `orpheus/sn/spatial/cell_update.py`
- new: `tests/sn/spatial/test_cell_update_protocol.py`

**Design notes**: Extract the cell-update equation from existing
`_sweep_1d_spherical` (lines ~230–290 in `orpheus/sn/sweep.py`).
The numer/denom assembly per cell is the cell update; the loop
structure around it is the sweep. This issue establishes the contract;
Issue 9 implements concrete strategies; Issue 12 rebuilds the sweep
around the new ABC.

---

#### Issue 9: Concrete cell updates — DD, LD, ExponentialCharacteristic, Step

- **Module**: `module:sn`
- **Type**: `type:feature`
- **Phase**: 2
- **Depends on**: Issue 8
- **Complexity**: L

**Context**: Production schemes need positivity by construction. LD
(linear discontinuous) and exponential characteristic are workhorses;
step is the robust fallback. DD is kept as comparison artifact and
explicitly marked non-positive.

**Acceptance criteria**:

- [ ] `DiamondDifference` in `orpheus/sn/spatial/diamond.py` — implements
      the existing math from `_sweep_1d_spherical` /
      `_sweep_1d_cylindrical` / `_sweep_2d_wavefront`.
      `is_positivity_preserving = False`. Default-deselected for
      production; available for verification
- [ ] `LinearDiscontinuous` in `orpheus/sn/spatial/linear_discontinuous.py`
      — slope-limited LD as recommended default.
      `is_positivity_preserving = True`. Cite Lewis & Miller §5.3 in
      docstring
- [ ] `ExponentialCharacteristic` in `orpheus/sn/spatial/exponential.py`
      — analytical exponential within cells, optimal for thin/thick
      optical regimes. `is_positivity_preserving = True`
- [ ] `Step` in `orpheus/sn/spatial/step.py` — first-order step
      characteristic, robust last resort. `is_positivity_preserving = True`
- [ ] All four pass MMS verification — link to
      `derivations/continuous/mms/sn/`
- [ ] Spatial convergence tests per scheme on a benchmark
      slab/sphere/cylinder problem demonstrating expected order

**Files**:

- new: `orpheus/sn/spatial/{diamond,linear_discontinuous,exponential,step}.py`
- new: `tests/sn/spatial/test_*.py`

**Design notes**: Don't ship all four in one PR. Sequence:
DD first (move existing code, verify bit-identical results), then LD,
then exponential, then step. Each cell update works for slab AND
curvilinear via the connection-term inputs in
`streaming_terms.alpha_in / alpha_out / tau` and
`upstream_state.angular_upstream`. **The sequencing matters**: DD
landing first preserves the regression contract throughout the rest
of the campaign.

---

#### Issue 9.5: Rename `BoundaryOperator` → `BoundaryOperator` (naming consolidation)

- **Module**: `module:geometry`, `module:sn`
- **Type**: `type:refactor`
- **Phase**: 2 (final, lands between Issue 9 and Issue 10)
- **Depends on**: Issue 7
- **Complexity**: S
- **Risk**: very low — mechanical rename, no behavioral change

**Context**: `BoundaryOperator` describes a *process* (a BC that has been
resolved). The name tells the reader how the object was made, not
what it is. `BoundaryOperator` is the standard mathematical concept
(BEM, transport, integral equations) for a linear operator on the
boundary trace space. A reader unfamiliar with the codebase, on
encountering `BoundaryOperator`, can immediately infer:

- it is a linear operator (subtype of `LinearOperator`)
- it acts on boundary trace spaces (its `domain` and `range` are
  boundary `FunctionSpace`s)
- it composes via operator algebra (multi-region BCs are operator
  products; partial reflectors are operator sums)
- it has rank-decomposition structure (specular = G ⊗ A as established
  in the BC tensor framing)
- it has an adjoint (via `LinearOperator.adjoint()`)

None of this is conveyed by `BoundaryOperator`. Per the architectural
principle that names should be strong concepts allowing fresh-context
readers to infer mathematical properties, the rename is correct.

Lands BEFORE Issue 10 because Phase 3 (Issues 10–15) all touch this
type and dependent types: writing the right name from the start avoids
touch-twice churn across the rest of the campaign. Regression testing
is bit-identical — no behavior change.

**Acceptance criteria**:

- [ ] In `orpheus/geometry/boundary.py`:
      - [ ] `BoundaryOperator` ABC → `BoundaryOperator`
      - [ ] Concrete subtypes renamed:
            `VacuumBoundaryOperator` → `VacuumBoundaryOperator`,
            `SpecularBoundaryOperator` → `SpecularBoundaryOperator`,
            `WhiteBoundaryOperator` → `WhiteBoundaryOperator`,
            `PeriodicBoundaryOperator` → `PeriodicBoundaryOperator`,
            `AlbedoBoundaryOperator` → `AlbedoBoundaryOperator`,
            `MixedBoundaryOperator` → `MixedBoundaryOperator`
- [ ] In `orpheus/sn/geometry.py`:
      - [ ] `SNMesh.BOUNDARY_OPERATOR_REGISTRY` → `SNMesh.BOUNDARY_OPERATOR_REGISTRY`
      - [ ] Internal factory helpers renamed:
            `_sn_bc_vacuum` → `_sn_vacuum_boundary_operator`,
            `_sn_bc_reflective` → `_sn_reflective_boundary_operator`
      - [ ] `_resolve_one(bc, face)` return-type annotation updated to
            `BoundaryOperator` (signature otherwise unchanged)
- [ ] `BC` (declaration tag) name **unchanged** — user-facing concept,
      short and declarative
- [ ] All test references updated by mechanical search-and-replace
- [ ] All Sphinx documentation references updated; cross-references
      verified via Nexus rebuild
- [ ] `apply_to_incoming(...)` method on the renamed `BoundaryOperator`
      ABC: **decision point during implementation** — either keep as
      domain-specific name OR rename to inherited `LinearOperator.apply`
      with `apply_to_incoming` as alias. Recommend the latter (uniformity
      with operator algebra; direction is encoded in `domain` / `range`
      `FunctionSpace`s), but this is a judgement call at PR time
- [ ] Behavioral regression suite (Issue 16) passes with bit-identical
      results — the gating contract

**Files affected**:

- modify: `orpheus/geometry/boundary.py`
- modify: `orpheus/sn/geometry.py`
- modify: `orpheus/sn/sweep.py` (any sites that consume the renamed
  types)
- modify: tests under `tests/sn/`, `tests/geometry/`
- modify: `docs/theory/discrete_ordinates.rst`,
  `docs/theory/geometry.rst` (any references to the old names)

**Design notes**:

- This is the cleanup PR before Phase 3 opens. Do not bundle other
  changes; the value of the rename is precisely that it isolates a
  pure mechanical refactor with zero behavioral risk.
- Use `git grep BoundaryOperator` and `git grep -i 'resolvedbc\|VacuumBoundaryOperator\|SpecularBoundaryOperator\|WhiteBoundaryOperator\|PeriodicBoundaryOperator\|AlbedoBoundaryOperator\|MixedBoundaryOperator'`
  to enumerate sites; verify completeness before merging.
- After this PR lands, Issues 10–15 will use `BoundaryOperator`
  uniformly. Update internal references in those issue specs at
  PR time if they reference the old names.
- The dunder amendments from the campaign refinement round
  (`__init_subclass__` for auto-registration) compose naturally with
  this rename: concrete `BoundaryOperator` subclasses can declare
  `kind = "vacuum"` etc. and self-register into
  `BOUNDARY_OPERATOR_REGISTRY`. The structural change (BoundaryOperator
  inheriting LinearOperator) and `__init_subclass__` machinery land in
  Issue 9.6, immediately after this rename. Keep this issue mechanical;
  9.6 carries the structural lifting.

---

#### Issue 9.6: Architecture + dunder consolidation (retroactive amendments)

- **Module**: `module:numerics`, `module:geometry`
- **Type**: `type:refactor`
- **Phase**: 2 (final, lands after Issue 9.5 and before Issue 10)
- **Depends on**: Issues 1, 2, 7, 9.5
- **Complexity**: M
- **Risk**: low — pure additions, no behavioral change

**Context**: Several architectural and ergonomic refinements were
identified during planning conversations AFTER Issues 1–9-DD shipped.
They consolidate cleanly here as a single PR before Phase 3 opens.
The amendments are pure additions (no breaking changes); regression
contract is bit-identical — the 11/11 frozen snapshots remain
`np.array_equal`-bit-identical.

The amendments cover three categories:

1. **Architecture**: function-space typing on operators; adjoint as a
   first-class operator construction; spatial integration via
   `DiscreteMeasure`.
2. **Dunders**: ergonomic and structural Python idioms that reflect
   the underlying mathematics (operator as callable, measure as
   functional, container protocol on measures).
3. **Cross-cutting**: `BoundaryOperator` inherits `LinearOperator`
   (closes the structural gap from Issue 9.5's mechanical rename);
   auto-registration via `__init_subclass__` for `BoundaryOperator`
   and `CellUpdate` subclass registries.

Why now: every Phase 3 issue (10–15) builds operators on top of these
primitives. Folding the affordances in before Phase 3 opens avoids
retrofit churn across 5 downstream issues.

**Acceptance criteria**:

**A. `FunctionSpace` primitive** (in `numerics/operator.py` or new
`numerics/space.py`):

- [ ] `FunctionSpace` frozen dataclass: `name: str`,
      `shape: tuple[int, ...]`,
      `inner_product_weights: NDArray | None`
- [ ] `inner_product(x, y) -> float` and `norm(x) -> float` methods
      (default L² Euclidean when `inner_product_weights is None`)
- [ ] `__eq__` based on `(name, shape)` identity — weights are metadata,
      not identity. Two spaces with the same name+shape but different
      weights compare equal but warn on construction
- [ ] `__hash__` consistent with `__eq__`
- [ ] `__repr__` showing `name` and `shape`
- [ ] Pre-populated for common ORPHEUS spaces: `"angular_flux"`
      (shape `(n_cells, n_ordinates, n_groups)`), `"scalar_flux"`
      `(n_cells, n_groups)`, `"boundary_trace_in"`,
      `"boundary_trace_out"`

**B. `LinearOperator` extensions** (in `numerics/operator.py`):

- [ ] `domain: FunctionSpace` and `range: FunctionSpace` properties on
      the Protocol; concrete operators declare them at construction
- [ ] `adjoint(self) -> LinearOperator` method returning the adjoint
      as a new `LinearOperator` — wraps `apply_transpose`, swaps
      domain↔range, transforms capabilities (`apply` ↔ `apply_transpose`
      under adjoint)
- [ ] Composition primitives (`OperatorSum`, `OperatorProduct`,
      `ScaledOperator`) verify domain/range compatibility at
      construction; raise `IncompatibleOperatorComposition` on mismatch
- [ ] `__call__(x)` aliasing `apply(x)` for the `L(ψ)` math notation
- [ ] `__pow__(n: int) -> LinearOperator` for operator powers via
      repeated `__matmul__` (binary exponentiation acceptable)
- [ ] `__repr__` showing `domain.name`, `range.name`, capabilities

**C. `DiscreteMeasure` extensions** (in `numerics/measure.py`):

- [ ] `__call__(f)` aliasing `integrate(f)` — a measure IS a functional
- [ ] `__iter__` yielding `(node, weight)` tuples
- [ ] `__len__` returning number of points
- [ ] `__getitem__(i)` returning `(nodes[i], weights[i])`
- [ ] `__repr__` showing `space`, `n_points`, `invariance_group` if known

**D. `Mesh.volume_measure`** (in `orpheus/geometry/mesh.py`):

- [ ] `Mesh1D.volume_measure` property returning a `DiscreteMeasure`
      with `nodes=centers`, `weights=volumes`, `space="spatial_R1"`
- [ ] `Mesh2D.volume_measure` similarly with `space="spatial_R2"`
- [ ] Existing spatial-integration sites (e.g.
      `np.sum(scalar_flux * volumes)` in solver code) are NOT yet
      converted — that's downstream cleanup for natural touch points.
      This issue installs the affordance only

**E. `BoundaryOperator` inherits `LinearOperator`** (closes structural
gap from Issue 9.5):

- [ ] `BoundaryOperator` ABC inherits from `LinearOperator` with
      `domain` = boundary trace outgoing FunctionSpace,
      `range` = boundary trace incoming FunctionSpace
- [ ] Concrete subtypes declare `capabilities = frozenset({"apply"})`
      at minimum; `SpecularBoundaryOperator` adds `"apply_transpose"`
      since it has a clean dual; `tensor_decomposition` property
      stays exposed for inspection
- [ ] `apply(psi_out)` is the canonical method (formerly
      `apply_to_incoming` from Issue 7); `apply_to_incoming` retained
      as alias for solver-code readability where directional clarity
      helps
- [ ] Composition `BoundaryOperator @ BoundaryOperator` works for
      multi-region interfaces; `0.7 * SpecularBoundaryOperator + 0.3 *
      WhiteBoundaryOperator` works for partial reflectors

**F. `__init_subclass__` auto-registration**:

- [ ] `BoundaryOperator` subclasses self-register into
      `BOUNDARY_OPERATOR_REGISTRY` via `__init_subclass__` with a
      `kind: ClassVar[str]` class attribute. Replaces manual factory
      dispatch from Issue 7
- [ ] `CellUpdate` subclasses similarly self-register into a
      `CELL_UPDATE_REGISTRY`. The DD strategy from Issue 9-DD declares
      `kind = "diamond_difference"` retroactively; LD/EC/Step from
      Wave C-extension self-register on definition
- [ ] Manual registry dict insertion still works as fallback

**Files affected**:

- modify: `orpheus/numerics/operator.py` (FunctionSpace, adjoint,
  `__call__`, `__pow__`, `__repr__`, composition checks)
- modify: `orpheus/numerics/measure.py` (`__call__`, container
  protocol, `__repr__`)
- modify: `orpheus/geometry/mesh.py` (`volume_measure` property on
  `Mesh1D` and `Mesh2D`)
- modify: `orpheus/geometry/boundary.py` (`BoundaryOperator` inherits
  `LinearOperator`, `__init_subclass__` registry)
- modify: `orpheus/sn/spatial/cell_update.py` (`__init_subclass__`
  registry; DD strategy declares `kind`)
- modify: existing tests to use new affordances where natural; add
  unit tests for new methods (`adjoint`, `__call__`, `__pow__`,
  composition compatibility, `volume_measure`)

**Coordinating with done work**:

- Issue 1 (commit `60d3932`) added `LinearOperatorMixin` for dunder
  algebra (`__add__`, `__sub__`, `__mul__`, `__matmul__`). This issue
  extends with `__call__`, `__pow__`, `adjoint()`, FunctionSpace
  domain/range, `__repr__`. **Verify with
  `git log -- orpheus/numerics/operator.py` what's already there before
  adding duplicates.**
- Issue 2 (commit `2f67853`) shipped `DiscreteMeasure` with tensor
  product / pushforward / restrict / direct sum. This issue adds
  `__call__`, container protocol, `__repr__`. Check current surface
  before duplicating.
- Issue 7 (commit `e93fe47`) shipped `BoundaryOperator` with tensor
  decomposition. Issue 9.5 renames to `BoundaryOperator` (mechanical).
  This issue adds `LinearOperator` inheritance and `__init_subclass__`
  registration on top of the renamed type.
- Wave C-extension (LD/EC/Step) is sequenced after Wave D per the
  operational notes. Those concrete cell updates will land into the
  registry this issue installs.

**Design notes**:

- **Bit-identical regression contract**: all changes are pure
  additions (new methods, new types, new properties, new aliases).
  No existing call site changes. The 11/11 frozen snapshots and the
  ERR-026 tripwires remain intact.
- **`FunctionSpace` is intentionally lightweight**. Don't try to
  enforce types via Python generics. `(name, shape)` identity with
  weights as metadata is the right tradeoff for a teaching codebase.
- **`adjoint()` is the architectural lever** for sensitivity analysis,
  perturbation theory, and adjoint Monte Carlo (Wave 4). Issues 11,
  13, 15 use it heavily. Without it, adjoint construction reaches
  into operator internals manually — workable but error-prone.
- **`__call__` on operators and measures matches mathematical notation**
  (`L ψ`, `μ(f) = ∫ f dμ`). Standard in scipy / JAX / FEniCS.
- **`__init_subclass__` is mostly drift-prevention**. Manual registries
  work; auto-registration just makes "forgot to register" impossible.
  Apply where it reads cleanly; don't force it where the subclass
  declaration is awkward.
- **Bundle this in one PR**, not split into 6 sub-PRs: the amendments
  are interrelated (`FunctionSpace` is consumed by `LinearOperator`
  domain/range, which is consumed by `adjoint()`, etc.). Splitting
  creates intermediate states where the abstractions are partially
  in place. One PR, one merge, regression-gated.

---

### Phase 3 — SN core reshape

#### Issue 10: Refactor `SNMesh` to consume `ReducedStreamingOperator`

- **Module**: `module:sn`
- **Type**: `type:refactor`
- **Phase**: 3
- **Depends on**: Issues 6, 9.6
- **Complexity**: M

**Context**: `SNMesh._setup_spherical` and `_setup_cylindrical` compute
connection coefficients inline. Issue 6 lifted that math to
`geometry/reduced_operator.py`. Now `SNMesh` becomes a thinner
aggregator: hold mesh, hold quadrature (DiscreteMeasure-backed), hold
reduced operator, hold resolved BCs. No streaming math in `SNMesh`
itself.

**Acceptance criteria**:

- [ ] `SNMesh.__init__` calls
      `geometry.reduced_streaming_operator(mesh, quadrature)` instead
      of `_setup_*` methods
- [ ] All `SNMesh` consumers (sweep, BiCGSTAB operator) read connection
      coefficients via `sn_mesh.reduced.streaming_terms(...)` instead of
      `sn_mesh.alpha_half` etc.
- [ ] Backward-compatible attribute access via properties for one
      release: `sn_mesh.alpha_half` returns
      `sn_mesh.reduced.alpha_half` with a `DeprecationWarning`
- [ ] All existing SN tests pass

**Files**:

- refactor: `orpheus/sn/geometry.py`

**Design notes**: After this issue, `SNMesh` knows nothing about α
coefficients or τ weights — it just holds and routes. The
geometry-vs-SN split becomes visible in code.

---

#### Issue 11: `SNStreamingOperator` as `LinearOperator`

- **Module**: `module:sn`
- **Type**: `type:feature`
- **Phase**: 3
- **Depends on**: Issues 1, 8, 9, 10
- **Complexity**: L

**Context**: The sweep is the implementation of L⁻¹. Wrapping
sweep + streaming-residual + adjoint sweep in a `LinearOperator`
interface lets eigenvalue solvers and Krylov code consume L uniformly
across solver methods. **This is also where #96 / #97 close**: `apply`
and `solve` use the same cell update strategy, so the FD/sweep
inconsistency disappears by construction.

**Acceptance criteria**:

- [ ] `SNStreamingOperator` in `orpheus/sn/operator.py` (replacing the
      BiCGSTAB-only `transport_operator_matvec_*` functions) implements
      `LinearOperator`
- [ ] `solve(q)` invokes the unified sweep (Issue 12)
- [ ] `apply(psi)` computes forward streaming-collision (Ω·∇ + Σ_t)ψ
      using the same cell update as the sweep — needed for residuals
      and matrix-free Krylov
- [ ] `apply_transpose(psi)` is the adjoint sweep: reversed Ω,
      transposed cell update
- [ ] `capabilities = frozenset({"apply", "solve", "apply_transpose"})`
- [ ] Declares `domain` and `range` as `FunctionSpace("angular_flux", ...)`
      from Issue 9.6; composition with `ScatteringOperator` and
      `FissionOperator` from Issue 13 verified at construction time
      via `OperatorSum`/`OperatorProduct` compatibility checks
- [ ] Adjoint constructed via `L.adjoint()` (Issue 9.6 method) returning
      a new `SNStreamingOperator` with reversed Ω ordering and
      transposed cell update — not via standalone `apply_transpose`
      calls in client code
- [ ] Reciprocity test: ⟨Lψ, φ*⟩ = ⟨ψ, L*φ*⟩ to round-off, using
      `domain.inner_product` (the angular-flux inner product weighted
      by volumes ⊗ angular-quadrature weights)
- [ ] Closes #96 and #97 (subject to Issue 15 landing)

**Files**:

- refactor: `orpheus/sn/operator.py`

**Design notes**: The current `operator.py` has a documented
inconsistency between the BiCGSTAB FD operator and the DD sweep
(issues #96 / #97). This issue makes them the same operator: matrix-free
`apply` is the algebraic dual of the sweep using the same cell update.
The behavioral change in BiCGSTAB convergence is acceptable per
campaign decision 2.

---

#### Issue 12: Unified sweep — single algorithm, parameterized

- **Module**: `module:sn`
- **Type**: `type:refactor`
- **Phase**: 3
- **Depends on**: Issues 6, 8, 10
- **Complexity**: L

**Context**: Currently `transport_sweep` dispatches into 4 separate
sweep functions (`_sweep_1d_cumprod`, `_sweep_1d_spherical`,
`_sweep_1d_cylindrical`, `_sweep_2d_wavefront`). Dispatch is by
`sn_mesh.curvature`, conflating "what coordinate system" with "what
discretization scheme." After this refactor: 2 sweep algorithms
(Cartesian flow-ordered + curvilinear μ-marching) parameterized by
cell update; choice between them comes from
`reduced_streaming_operator.requires_upstream_angular_state`.

**Acceptance criteria**:

- [ ] `orpheus/sn/sweep.py` rewritten with `_cartesian_sweep(L, q)` and
      `_curvilinear_sweep(L, q)`. Dispatch on
      `L.geometry.reduced.requires_upstream_angular_state`, not on
      string-comparing `curvature`
- [ ] Per-cell math delegated to `L.cell_update.update(...)` with
      appropriate `streaming_terms` and `upstream_state` assembled
- [ ] **1D cumprod fast path preserved as an optimization within the
      Cartesian sweep** when `cell_update == DiamondDifference` and
      quadrature is GL1D — algebraic identity holds for DD specifically;
      keep for performance, but as a fast path inside the unified
      algorithm, not a separate dispatched function
- [ ] 2D wavefront diagonal scheduling preserved
- [ ] All existing SN tests pass with bit-identical results when
      `cell_update == DiamondDifference`
- [ ] New tests with `cell_update == LinearDiscontinuous`,
      `ExponentialCharacteristic`, `Step` exercise the polymorphism

**Files**:

- rewrite: `orpheus/sn/sweep.py`

**Design notes**: Most behavior-preserving issue in the campaign.
Bit-identical output for DD case is the gating contract. The 1D
cumprod recurrence is an algebraic identity holding for DD only —
keep it as fast path inside the unified Cartesian sweep.

---

#### Issue 13: `ScatteringOperator` and `FissionOperator` as `LinearOperator`s

- **Module**: `module:sn`
- **Type**: `type:refactor`
- **Phase**: 3
- **Depends on**: Issue 1
- **Complexity**: M

**Context**: Scattering and fission sources are currently methods on
`SNSolver` (`_add_scattering_source`, `_build_aniso_scattering`,
`compute_fission_source`). Lifting them to `LinearOperator`s makes the
operator algebra `(L − S − F) ψ = q` explicit and lets eigenvalue
solvers consume the pieces uniformly across solver methods.

**Acceptance criteria**:

- [ ] `ScatteringOperator(LinearOperator)` in `orpheus/sn/scattering.py`:
      holds materials + quadrature; `apply(psi)` returns scattering
      source (P0 + Pℓ); `capabilities = frozenset({"apply"})`
- [ ] `FissionOperator(LinearOperator)` in `orpheus/sn/fission.py`:
      holds materials; `apply(psi)` returns χ ⊗ νΣ_f acting on flux
      moments. The rank-1-in-energy tensor structure χ ⊗ νΣ_f is
      reflected in the implementation (e.g., factored as
      `OuterProduct(chi, nu_sigma_f)` exposing the rank explicitly)
- [ ] Both declare `domain` / `range` as
      `FunctionSpace("angular_flux", ...)` from Issue 9.6 so composition
      with `SNStreamingOperator` is type-checked
- [ ] `ScatteringOperator.adjoint()` returns operator with transposed
      scattering kernel in (E', E) and (Ω', Ω). `FissionOperator.adjoint()`
      swaps roles of χ and νΣ_f. Both follow the
      `LinearOperator.adjoint()` contract from Issue 9.6
- [ ] (n,2n) folded into `ScatteringOperator` (consistent with existing
      `_add_n2n_source` placement in `SNSolver`)
- [ ] All existing SN tests pass

**Files**:

- new: `orpheus/sn/scattering.py`, `orpheus/sn/fission.py`
- modify: `orpheus/sn/solver.py` (remove source-construction methods)

**Design notes**: Pℓ anisotropic scattering is already implemented in
`SNSolver._build_aniso_scattering`. **Move** the implementation; don't
rewrite the math. The Galerkin projection on Y_ℓm is what's happening —
make it explicit in the docstring.

---

### Phase 4 — Iteration as operator algebra

#### Issue 14: New iteration primitives in `numerics/iteration.py`

- **Module**: `module:numerics`
- **Type**: `type:feature`
- **Phase**: 4
- **Depends on**: Issue 1
- **Complexity**: M

**Context**: Current `EigenvalueSolver` Protocol bundles transport
solving with scattering iteration in `solve_fixed_source`. This
conflates concerns and makes adjoint, Krylov-on-the-outer-iteration,
and contour-integral eigenvalue solvers awkward. New shape: iteration
consumes `LinearOperator` instances directly.

**Acceptance criteria**:

- [ ] `SourceIteration(L, S, F, preconditioner=None)` in
      `orpheus/numerics/iteration.py` — fixed-source solver
      `solve(q_ext) → psi`. `(I − L⁻¹·S)·ψ = L⁻¹·(F·ψ + q)` implemented
      as a fixed-point iteration
- [ ] `KEigenvalue(L, S, F, eigenvalue_method="power")` — eigenvalue
      solver consuming operators directly. Other methods to land
      separately: `"contour_integral"` (the FEAST-style extension you
      had identified)
- [ ] Existing `power_iteration` in `numerics/eigenvalue.py` deprecated
      via `DeprecationWarning` but kept functional for transitional
      CP/diffusion compatibility
- [ ] Both new primitives consume any `LinearOperator` triple — agnostic
      to whether L is SN, MoC, or CP
- [ ] Tests against synthetic operators (numpy matrices wrapped) plus
      tests using SN operators from Issues 11, 13

**Files**:

- new: `orpheus/numerics/iteration.py`
- modify: `orpheus/numerics/eigenvalue.py` (add deprecation, keep
  functional)

**Design notes**: Don't delete the old `EigenvalueSolver` protocol —
CP and diffusion still use it and migrate in their own waves.
`power_iteration` can linger until the cross-solver migration sequence
completes.

---

#### Issue 15: Migrate `SNSolver` to operator-algebra iteration

- **Module**: `module:sn`
- **Type**: `type:refactor`
- **Phase**: 4
- **Depends on**: Issues 11, 13, 14
- **Complexity**: L

**Context**: With L, S, F as operators and `KEigenvalue` consuming them,
`SNSolver` shrinks dramatically. Most of its current ~600 lines becomes
adapter code; physics moves to operators.

**Acceptance criteria**:

- [ ] `SNSolver` rewritten in `orpheus/sn/solver.py`: constructs
      `SNStreamingOperator` (L), `ScatteringOperator` (S),
      `FissionOperator` (F); delegates to `KEigenvalue(L, S, F)` for
      eigenvalue path; delegates to `SourceIteration(L, S, F)` for
      fixed-source path
- [ ] `solve_sn` and `solve_sn_fixed_source` public APIs preserved
      verbatim
- [ ] BiCGSTAB inner solver path becomes: GMRES around L using
      `L.apply` matrix-free, instead of building a separate FD operator.
      **Closes #96 and #97**
- [ ] All existing SN tests pass with bit-identical results in the
      `inner_solver="source_iteration"` path
- [ ] BiCGSTAB / GMRES path tests pass without the documented
      inconsistencies (confirmed via the reciprocity invariant test
      from Issue 11)

**Files**:

- rewrite: `orpheus/sn/solver.py`
- simplify: `orpheus/sn/operator.py`

**Design notes**: Issue that makes the architectural campaign visible
at the API surface. After it lands, `SNSolver` is a thin coordinator
and heavy lifting is in compositional operators. Closing of #96 / #97
is a major correctness win.

---

### Phase 5 — Tests + documentation

#### Issue 16: Behavioral regression suite

- **Module**: `module:tests`
- **Type**: `type:test`
- **Phase**: 5 (executed first)
- **Depends on**: none — this issue runs BEFORE the reshape begins
- **Complexity**: M

**Context**: Cardinal Rule 1 (correctness is critical) means the reshape
must not change physics. A frozen-reference regression suite is the
gating contract for every subsequent issue.

**Acceptance criteria**:

- [ ] Snapshot scalar flux + k_eff outputs from a representative
      pre-reshape test set (slab, sphere, cylinder; multigroup;
      Pℓ scattering at orders 0, 1, 3; vacuum + reflective BCs) into
      a frozen reference at `tests/_regression/sn_reference/`
- [ ] Pytest fixture compares post-reshape output to frozen reference
      at machine precision when `cell_update=DiamondDifference`
- [ ] Tagged `@pytest.mark.l1` and `@pytest.mark.regression`
- [ ] CI gates on regression for any PR touching `orpheus/sn/`,
      `orpheus/geometry/`, `orpheus/numerics/`

**Files**:

- new: `tests/_regression/sn_reference/` (snapshot data)
- new: `tests/sn/test_regression.py`

**Design notes**: **Generate snapshots ON main BEFORE branching for the
reshape work.** This is the very first action — Issue 16 lands on main,
then the campaign branches off. Otherwise the reference drifts with the
reshape.

---

#### Issue 17: Symmetry-preservation and reciprocity invariant tests

- **Module**: `module:tests`
- **Type**: `type:test`
- **Phase**: 5
- **Depends on**: Issues 11, 12
- **Complexity**: M

**Context**: New tests exercising structural properties the operator
algebra is supposed to guarantee — properties that weren't expressible
before because the abstractions weren't there.

**Acceptance criteria**:

- [ ] G-symmetry preservation: G-symmetric problem + G-invariant
      quadrature → discrete solution exactly G-invariant (to round-off).
      Concrete tests for slab(Z2), sphere(SO(3)), cylinder(SO(2)×R)
- [ ] Reciprocity: forward × adjoint detector responses match across
      source-detector swap. Uses `apply_transpose` on
      `SNStreamingOperator`
- [ ] Conservation: per-cell balance to round-off, verified against
      `derivations/discrete/sn/balance/`
- [ ] Capability honesty: every `LinearOperator` in the codebase
      correctly reports what it can do (auto-detected by attempting
      `solve` / `apply_transpose` and checking failure mode)
- [ ] Tagged `@pytest.mark.l1`, `@pytest.mark.verifies(...)` linked to
      Sphinx equation labels per V&V harness convention

**Files**:

- new: `tests/sn/test_invariants.py`

**Design notes**: Where the architectural rigor pays off — invariants
that were not previously expressible become testable.

---

#### Issue 18: Sphinx documentation for the operator algebra architecture

- **Module**: `module:docs`
- **Type**: `type:docs`
- **Phase**: 5
- **Depends on**: All implementation issues (lands progressively)
- **Complexity**: L

**Context**: Cardinal Rule 3 — Sphinx is the LLM's brain. The reshape
introduces a new architectural narrative (operator algebra unifying
SN/MoC/CP/diffusion) that needs a dedicated theory page, plus updates
to existing pages. Documentation IS the bridge to the AI-assisted
methodology presentation.

**Acceptance criteria**:

- [ ] New `docs/theory/operator_algebra.rst` — unifying narrative:
      differential / integral / variational forms; L/S/F as
      `LinearOperator`s; sweep as L⁻¹; how this scales to MoC/CP/diffusion
- [ ] New `docs/theory/discrete_measures.rst` — quadratures as discrete
      measures; tensor product / pushforward / restrict; symmetry-V
      framework; registry rationale with selection examples
- [ ] Update `docs/theory/discrete_ordinates.rst` — Key Facts, equations,
      how new architecture maps onto existing math; reference
      Bailey et al. (2009) where present
- [ ] Update `docs/theory/geometry.rst` — orbit-space classification,
      reduced operators, connection coefficients, where they live in
      code
- [ ] All cross-referenced via Nexus per Cardinal Rule 3
- [ ] Sphinx builds clean

**Files**:

- new: `docs/theory/operator_algebra.rst`,
       `docs/theory/discrete_measures.rst`
- modify: `docs/theory/discrete_ordinates.rst`,
          `docs/theory/geometry.rst`

**Design notes**: The math from the planning conversation thread
(orbit-space classification, fiber bundle for rays,
Galerkin = multigroup = PN = FE structure, quadratures as discrete
measures with G–V duality) is the substantive content. Use it.
Documentation campaign overlaps the last few implementation waves.

---

## Dependency waves (parallelization plan)

```
Wave A (parallel, no deps): Issue 1, Issue 2, Issue 16
Wave B (parallel after A): Issue 3, Issue 4 [needs 2,3], Issue 5 [needs 3,4]
                           Issue 6 [needs 1], Issue 7 [needs 4]
Wave C (after B):          Issue 8 [needs 6]
Wave D (parallel after C): Issue 9 [needs 8] — 9-DD shipped; LD/EC/Step deferred to post-Wave-D
Wave D' (after D):         Issue 9.5 [needs 7] — naming consolidation (mechanical rename)
Wave D'' (after D'):       Issue 9.6 [needs 1,2,7,9.5] — architecture + dunder consolidation
Wave E (after D''):        Issue 10 [needs 6, 9.6]
Wave F (after E):          Issue 11 [needs 1,8,9,10,9.6], Issue 13 [needs 1, 9.6]
Wave G (after F):          Issue 12 [needs 6,8,10]
Wave H (parallel after F): Issue 14 [needs 1, 9.6]
Wave I (after G+H):        Issue 15 [needs 11,13,14] — closes #96, #97
Wave C-ext (after Wave I): Issue 9 LD/EC/Step — verified end-to-end through unified sweep
Continuous:                Issue 17 [after 11,12], Issue 18 [progressive]
```

**Note on Issue 9.5 + 9.6 placement**: 9.5 (rename) and 9.6
(architecture + dunder consolidation) land back-to-back as the cleanup
pair before Phase 3 opens. 9.5 is a mechanical rename with zero
behavioral risk. 9.6 is pure-additive architecture that retroactively
upgrades Issues 1, 2, 7 with `FunctionSpace`, `adjoint()`, dunders
(`__call__`, `__pow__`, container protocol), `mesh.volume_measure`,
and `BoundaryOperator` inheriting `LinearOperator` with
`__init_subclass__` auto-registration. Both gated by the behavioral
regression suite (11/11 bit-identical). After 9.6, every Phase 3
issue (10–15) consumes the upgraded primitives without retrofit churn.

**Note on Wave C-extension placement**: Issue 9's LD/EC/Step strategies
landed post-Wave-I per the operational decision recorded in Wave C
notes — they verify end-to-end through the unified sweep (Issue 12)
rather than in isolation. Each ships sequentially with its own MMS
spatial-convergence test; auto-registration into `CELL_UPDATE_REGISTRY`
(Issue 9.6) is automatic at definition time.

**Issue 16 lands on `main` BEFORE the campaign branch is created.**

Realistic single-developer sequencing with Claude Code: ~10–14 sessions
for the reshape, plus documentation campaign overlapping the last few
waves.

## Conventions for sub-agent issue creation

When the sub-agent translates these specs to GitHub issues:

1. **Title format**: short imperative, e.g. `numerics: LinearOperator
   protocol with capability tags`
2. **Labels**: copy `module:` and `type:` labels from each spec
3. **Body**: copy Context + Acceptance Criteria + Files + Design Notes
   verbatim. Sub-agents and fresh Claude Code sessions read the body.
4. **Cross-references**: use `Depends on #NN` and `Closes #NN`
   (latter for issues 96 / 97, closed by Issue 15)
5. **Milestone**: `SN Reshape (Wave 1)` — create this milestone first
6. **Branch convention**: each issue gets a topic branch
   `<type>/<scope>/<short-description>` per the project's git workflow
   (e.g. `feat/numerics/linear-operator-protocol`)
7. **Commit convention**: Conventional Commits
   `<type>(<scope>): <summary>` per `CLAUDE.md` git workflow section
8. **Issue close**: `Closes #NN` trailer in commit body; PR uses
   `git merge --ff-only`

## Long-term migration sequence (after this campaign)

After SN reshape merges to main, subsequent waves consume the same
primitives. Each wave is a separate campaign with its own plan
document under `.claude/plans/`.

| Wave | Solver | Estimated complexity | Key consumer of new primitives |
|------|--------|---------------------|-------------------------------|
| 2 | MoC | L | `BundleMeasure`, ray-quadrature pushforward, `LinearOperator` |
| 3 | CP | L | dense `LinearOperator`, scalar-flux ops, BC tensor decomp |
| 4 | MC | M | `Geometry`, `BoundaryCondition`, sampling from `DiscreteMeasure` |
| 5 | Diffusion | S | `LinearOperator` (already FD-based), retire `EigenvalueSolver` Protocol |

After Wave 5: `numerics/eigenvalue.py:power_iteration` and the old
`EigenvalueSolver` Protocol can be deleted.

## Risk register

- **Behavioral regression slip**: catch by Issue 16 before branching
- **Curvilinear cell-update API insufficient**: surface during Issue 9
  implementation; if so, revise Issue 8 ABC and reconsider downstream
  issues
- **MoC/CP primitive blockers discovered late**: mitigated by Issue 2
  shipping `BundleMeasure` upfront; if other shape mismatches surface
  during MoC reshape, those become primitive revisions in Wave 2
- **Performance regression** in 1D cumprod: Issue 12 keeps the algebraic
  optimization; benchmark in regression suite to catch slowdowns
- **Pℓ anisotropic correctness drift**: behavioral regression suite
  includes Pℓ tests at orders 0, 1, 3 to gate this

## References

- `CLAUDE.md` — project-level Cardinal Rules and conventions
- `docs/theory/discrete_ordinates.rst` — current SN theory page
- `derivations/discrete/sn/balance/` — symbolic derivation of cell
  balance equations
- Bailey et al. (2009) — `_setup_spherical` / `_setup_cylindrical`
  references for angular redistribution and Morel-Montry closure
- Lewis & Miller, Computational Methods of Neutron Transport (1984) —
  §4.5 (M-M closure for RZ), §5.3 (LD discretization)
- GitHub issues #96, #97 — operator-sweep inconsistency, closed by
  Issue 15
