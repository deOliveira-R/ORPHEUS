# Phase G replan — SN four-operator unification (Issue #196)

**Tracking issue**: Issue [#196](https://github.com/deOliveira-R/ORPHEUS/issues/196) — ERR-026 manifestation #7 (SI-vs-Krylov O(h) WDD spatial-closure asymmetry on curvilinear SN). The original Phase G plan (`.claude/plans/issue_196_phase_g_four_operator_unification.md`) is **superseded** by this replan.

**Branch**: `refactor/sn-operator-algebra` (continue from Phase F tip `b0cc1b1`; Step 1's mis-abstraction commits `6eeff94/2b73e6e/dda6f28` and the diagnostic commit `be89c4e` to be partially preserved — see Step 0).

---

## Context — why this replan

The original Phase G plan failed in two phases:

**Step 1 (committed as `6eeff94/2b73e6e/dda6f28`)** promoted `DiamondDifference` and `MorelMontryAngularSweep` to `LinearOperator` subclasses (`SNCellOperator`, `AngularRedistribution`). They are **strategies**, not operators. A `CellUpdate` Protocol already exists at `orpheus/sn/spatial/cell_update.py:418`; `PoleAngularClosure` already exists at `orpheus/sn/spatial/pole_angular_closure.py:179`. The wrappers added one new piece of functionality (a per-cell residual computation) and otherwise delegated. The right home for the residual is a new `CellUpdate.residual` method, not a `LinearOperator` wrapper.

**Step 2 (uncommitted, working tree)** routed the curvilinear SI sweep through `scipy.sparse.linalg.gmres`, dismissing the actual closure-math fix as "Picard diverges" without showing the work. This violates Cardinal Rule 1 (correctness over shipping; no lazy solutions). It also created `DiscreteOrdinatesPhaseSpace` as a wrapper around `SNMesh` — Pattern 2 violation: SNMesh IS the SN-augmented phase space.

The three exploration memos commissioned for this replan establish the corrected baseline:

- `/Users/rodrigo/.claude/plans/crystalline-wondering-token-agent-aebacd5082553635d.md` — SNMesh + structural backbone.
- `/Users/rodrigo/.claude/plans/crystalline-wondering-token-agent-a20b4e206a7ec0c78.md` — SN strategy Protocol map (CellUpdate, PsiHalfAngleSeed, PoleAngularClosure, BoundaryRealizer).
- `/Users/rodrigo/.claude/plans/crystalline-wondering-token-agent-a7f5a72ef8522d183.md` — LinearOperator algebra + iteration + S/F surface area inventory.

**Key findings** that reshape the plan:

1. **Almost everything is already built.** `LinearOperator`, `LinearOperatorMixin`, `OperatorSum`, `OperatorProduct`, `ScaledOperator`, `IdentityOperator`, `ZeroOperator`, `_AdjointOperator` (with `.H` propagation), `as_scipy_linop`, `MissingCapability` — all in `orpheus/numerics/operator.py:88-1368`. `ScatteringOperator(LinearOperator)` with `S = R · Λ · M` factoring (imperatively composed inside `apply`) in `orpheus/sn/scattering.py:254-524`. `FissionOperator(LinearOperator)` rank-1 `χ ⊗ νΣ_f` in `orpheus/sn/fission.py:82-165`. `SourceIteration(L, S, F, q_ext)` and `KEigenvalue(L, S, F, *, keff_estimator)` in `orpheus/numerics/iteration.py:164-569` — **exactly the Step 3 vision shape**. `compute_keff`, `compute_group_production_rate`, `compute_group_absorption_rate` (Issue #169) in `orpheus/sn/solver.py:295-358`.

2. **`SNMesh` IS the SN-augmented phase space.** Constructor at `orpheus/sn/geometry.py:125` takes `(mesh, quadrature, cell_update, pole_angular_closure)`. It carries the streaming stencil, the BC realisations (`bc_left/right/xmin/xmax/ymin/ymax` are already realised `LinearOperator`s per Wave 9 / 11), and the cell-update + pole-angular-closure strategies. The only thing it does **not** carry is `n_groups`. Two options for the four-operator constructors: (a) take `n_groups` as a separate scalar parameter, or (b) add `n_groups` to `SNMesh.__init__` (one line). Either is fine; **no wrapper type is needed**.

3. **The canonical closure math is well-characterised.** Phase G Step 2 diagnostic (committed at `be89c4e`) attributes the SI-vs-Krylov drift to three structural differences: (1) M-M recurrence runs inside per-cell update on SI vs as separate pass on Krylov — Krylov form preserves Pomraning 1989 isotropy at the pole, SI form breaks it (22% L0 error, mesh-independent); (2) Carlson seed source-driven vs ψ-driven, equivalent at fixed point; (3) BC trace face-flux (SI form correct) vs cell-centre proxy (Krylov apply has O(h) artefact). The canonical closure is **separate-pass M-M + WDD streaming + face-flux BC trace + Carlson seed (either form)**. Picard iteration on this corrected closure has bounded spectral radius — the previous Step 2 agent's "Picard diverges" claim was about a different (broken) operator.

4. **The 4-operator algebra is one structural cleanup from existing.** The split into `L = pure streaming + Σ_t lifted out`, `C = collision`, `S = scattering`, `F = fission` is partially shipped in `orpheus/sn/streaming.py:251-512` but with two flaws (DiscreteOrdinatesPhaseSpace wrapper, dormant `build_equation_map_spherical` bug in `CollisionOperator.apply`, broken `from .mesh import SNMesh` import). The legacy `SNStreamingOperator` in `orpheus/sn/operator.py:1115-1473` still carries Σ_t internally and is what `SNSolver.__init__` uses today (`solver.py:240-262`).

5. **The L1 SN bridge test already passes.** `tests/numerics/test_iteration.py::test_keigenvalue_matches_solve_sn_2g_slab` (line 328) constructs `(SNStreamingOperator, ScatteringOperator, FissionOperator)` and runs `KEigenvalue` directly against `solve_sn` reference — proves the iteration primitives are ready to consume the SN operator triple. Adapter shims on each operator (lines 402-440) reconcile shape mismatches; these are the "small bridges" Step 4 removes.

This means the **remaining Phase G work** is smaller and sharper than the original plan implied: roll back the wrong abstractions, fix the strategy-layer extension with `CellUpdate.residual`, fix the canonical curvilinear closure (without GMRES shortcuts), wire `SNSolver` to consume `SourceIteration`/`KEigenvalue` directly, and implement `.H` on leaves.

---

## Acceptance criteria for Phase G closure

ALL of the following hold simultaneously after this replan executes:

1. **No twin-path abstractions.** The four user-facing operators are L, C, S, F plus B composed inside L. No `SNCellOperator`/`AngularRedistribution` wrapper types. No `DiscreteOrdinatesPhaseSpace` wrapper. No `SNStreamingOperator`-vs-`StreamingOperator` duplication.

2. **The strategy layer is properly extended.** `CellUpdate` Protocol has a `residual(...)` method; `DiamondDifference` implements both `update` (the solve direction) and `residual` (the apply direction); both consume the shared `cell_balance_terms` helper in `orpheus/sn/spatial/cell_balance.py`.

3. **The curvilinear SI sweep is correct.** Picard iteration on the canonical closure (separate-pass M-M + WDD + face-flux BC + Carlson seed). NO `scipy.sparse.linalg.gmres` anywhere in the SI sweep path. `inner_solver="source_iteration"` uses a sweep; `inner_solver="krylov"` uses GMRES on `(L+C).apply`. The L0 streaming-equilibrium test (promoted from `tests/sn/diagnostics/phase_g_step2_05_homogeneous.py`) passes at machine precision on BOTH inner solvers.

4. **L is pure streaming, C is separate.** `L = StreamingOperator(SNMesh, *, boundary=...)` carries no σ_t. `C = CollisionOperator(SNMesh, sigma_t)`. `(L + C).solve(q)` is the canonical curvilinear sweep via `OperatorSum.solve` fusion or via a wired `solve` method. `SNSolver.L` becomes the pure-streaming version; σ_t lives on `SNSolver.C`. The legacy `SNStreamingOperator` (`orpheus/sn/operator.py:1115`) is either deleted or kept as a backward-compat alias for one merge cycle then removed.

5. **`SNSolver` consumes the algebra.** `_solve_source_iteration` and `_solve_fixed_source_si` collapse to one call site that delegates to `SourceIteration(L, S, F, q_ext).solve(...)`. `_solve_krylov` and `_solve_fixed_source_krylov` collapse to one call site that builds `A_loss = L + C - S` and calls `scipy.sparse.linalg.gmres` via `as_scipy_linop(A_loss)` with `(L+C).solve` as preconditioner. `power_iteration` is replaced by `KEigenvalue(L, S, F, *, keff_estimator=lambda L_,S_,F_,phi: solver.compute_keff(phi))`.

6. **`.H` works.** `solve_sn_adjoint(materials, mesh, quad, response, *, tol)` is implemented in ≤ 8 lines as `(L + C - S).H.solve(response.as_source())` (via GMRES on `as_scipy_linop((L+C-S).H)`). `L.H` reverses streaming direction (= adjoint sweep). `C.H = C`. `S.H` advertises `CAP_APPLY_TRANSPOSE` and propagates `R.H @ Λ.H @ M.H`. `F.H` swaps the rank-1 factors. The 2G slab reciprocity test `<response, ψ_forward> = <ψ_adjoint, source>` passes at `rtol=1e-10`.

7. **Issue #196 manifestation table closes.** Manifestation #7 closes by construction (one canonical closure consumed by both inner solvers). Manifestation #6 fully closes (Phase F's partial closure resolves under correct Step 2). The newly-discovered L0 manifestation closes via the promoted L0 test. Manifestation #5 (L1 magnitude pre-asymptotic) stays in Issue #195 — separate.

8. **Snapshot principled-equivalence.** The 5 Cartesian snapshots stay bit-identical at `rtol=1e-12`. The 6 curvilinear snapshots regenerate under the canonical closure with three-pillar attestation: (a) L0 streaming-equilibrium analytical answer, (b) Pomraning 1989 pole isotropy, (c) Variant α integral form via existing `tests/sn/test_phase_c_crosscheck.py`. Documented per `vv-principles` §"Bit-identity vs principled-equivalence".

9. **Phase E sentinel re-enabled.** `tests/sn/test_phase_c_crosscheck.py::test_phase_e_trajectory_resolvent_flux_shape_crosscheck` is predicted XPASS after Step 2; the `xfail-strict=True` marker is removed.

10. **Sphinx narrative wired.** `docs/theory/operator_algebra.rst` (existing 796-line Phase 0 stub) is extended with the Phase G unification chapter. `verifies(...)` decorators pin tests to the operator-algebra `:label:` references. Sphinx `-W` builds clean.

---

## Step 0 — Rollback the mis-abstractions

**Goal**: clean working tree. Get to a state where the four-operator algebra can be built on top of correctly-extended strategies, with only the actually-salvageable pieces preserved.

### What stays committed (do NOT revert)

- Commits before Step 1: `b0cc1b1` (Phase F tip) and prior. These are sound.
- The diagnostic commit `be89c4e` — diagnostic scripts under `tests/sn/diagnostics/`, the gate memo at `.claude/agent-memory/test-architect/issue_196_phase_g_step1_verification_gates.md`, and the diagnostic memo at `.claude/agent-memory/numerics-investigator/issue_196_phase_g_step2_diagnostic.md`. These ground the replan; preserve them.

### What gets reverted

- `6eeff94 feat(sn): SNCellOperator + AngularRedistribution LinearOperator promotion` — the production code in `orpheus/sn/spatial/operators.py` plus the export in `orpheus/sn/spatial/__init__.py`. **Preserve** the residual-computation math (`_apply_curvilinear_residual`, `_apply_cylindrical_degenerate_residual`) in a scratch buffer — Step 1 ports it into `DiamondDifference.residual`.
- `2b73e6e test(sn): apply-vs-solve invariants for cell + angular ops` — `tests/sn/spatial/test_sncell_operator.py`, `test_angular_redistribution.py`, and the two new appendices to `test_sweep_vs_apply_consistency.py`. **Preserve** the test contracts that exercise underlying invariants (bit-identity vs `DD.update`, M-M algebraic identities, Carlson seed equivalence) — Step 1 ports them to test the strategies, not the wrappers.
- `dda6f28 chore(agent-memory): index Phase G Step 1 closeout` — keep the closeout memo on disk but mark its decisions as superseded by this replan.

### What gets dropped from the working tree

- `orpheus/sn/streaming.py` (untracked) — the half-shipped Step 2 work with `DiscreteOrdinatesPhaseSpace`, the GMRES `solve_within_group_sn_curvilinear`, the broken `from .mesh import SNMesh` import. **Delete entirely.** Step 3 rebuilds `StreamingOperator` and `CollisionOperator` correctly from scratch (in a different module home; see Step 3).
- Working-tree modifications to `orpheus/sn/sweep.py` — `_curvilinear_sweep` currently dispatches to `_curvilinear_sweep_via_canonical_operator` (GMRES). **Revert to the pre-Step-2 sweep dispatch**, which means restoring `_sweep_1d_spherical` and `_sweep_1d_cylindrical` (if they were deleted; verify with `git diff be89c4e -- orpheus/sn/sweep.py`).
- Working-tree modifications to `orpheus/sn/spatial/diamond.py` and `orpheus/sn/spatial/operators.py` — these consume `cell_balance_terms`. **Keep these modifications conceptually** (cell_balance.py is the qa CONCERN-2 fix); but commit them properly in Step 1 alongside the new `CellUpdate.residual` method.
- `orpheus/sn/spatial/cell_balance.py` (untracked, 235 lines) — **keep on disk; commit in Step 1**.

### Process

1. Inspect `git diff be89c4e -- orpheus/sn/sweep.py orpheus/sn/spatial/` to see exactly what the Step 2 agent changed.
2. Create a scratch buffer with the residual-computation math from `orpheus/sn/spatial/operators.py:319-367` for Step 1's port.
3. `git restore` the modified-tracked files. `rm orpheus/sn/streaming.py`.
4. Revert commits `6eeff94`, `2b73e6e`, `dda6f28` via `git revert` (one revert commit per source commit, with clear "Revert" subjects citing this replan). The closeout memo at `.claude/agent-memory/method-implementer/issue_196_phase_g_step1_closeout.md` can be left in place but should gain a top-line note "SUPERSEDED by replan at .claude/plans/crystalline-wondering-token.md".
5. After revert, verify the test suite is GREEN at the Phase F tip + diagnostic commit state: `pytest tests/sn/ -q`.
6. Reapply `cell_balance.py` cleanly in Step 1 (not in Step 0).

### Acceptance

- `git status` shows clean working tree.
- `git log` shows the diagnostic commit at HEAD plus three revert commits below it (or one squashed revert commit if cleaner).
- `pytest tests/sn/ -q` runs to the same green count as immediately after `be89c4e` was committed (i.e., 236 passed / 5 skipped / 1 xfailed for the sweep+regression suite, plus the full SN suite at its prior state).
- `tests/sn/spatial/test_sncell_operator.py` and `tests/sn/spatial/test_angular_redistribution.py` no longer exist (they were committed in `2b73e6e`).
- `orpheus/sn/spatial/operators.py` no longer exists (committed in `6eeff94`).
- `tests/sn/diagnostics/phase_g_step2_*.py` still exist (they're in the diagnostic commit, which we keep).

### Deliverable

Single commit (or up to three revert commits if cleaner): `revert(sn): Phase G Step 1 + Step 2 mis-abstractions (replan)`. Body cites this plan file and the user's correction.

---

## Step 1 — Strategy layer extension: `CellUpdate.residual`

**Goal**: extend the existing `CellUpdate` Protocol with the per-cell residual operation, on the strategy layer where it belongs. No new `LinearOperator` types. Preserve the salvageable math + tests from the reverted Step 1.

### Pre-step: read these files in full

- `orpheus/sn/spatial/cell_update.py` (612 lines) — the `CellUpdate` Protocol at line 418, `CellUpdateBase` ABC at line 510, companion dataclasses `CellVisit`/`UpstreamState`/`CellResult`/`SweepCellSlice`.
- `orpheus/sn/spatial/diamond.py` (645 lines) — `DiamondDifference(CellUpdateBase, key="diamond_difference")` at line 280, with `_update_slab`/`_update_curvilinear`/`_update_cylindrical_degenerate` as `@staticmethod`s.
- `orpheus/sn/spatial/cell_balance.py` (untracked, 235 lines) — the `CellBalanceTerms` dataclass + `cell_balance_terms` / `cell_balance_terms_degenerate` helpers. This is the qa CONCERN-2 fix.
- The reverted commits' residual-computation math (preserved in scratch buffer from Step 0): `_apply_curvilinear_residual` and `_apply_cylindrical_degenerate_residual`.

### Anti-recommendations (the agent MUST NOT do these)

- Do NOT create `SNCellOperator(LinearOperatorMixin)`, `AngularRedistribution(LinearOperatorMixin)`, or any wrapper type. The Step 0 revert removed these for a reason.
- Do NOT add `apply` / `solve` / `capabilities` / `.H` / `.T` to `CellUpdate` or `DiamondDifference`. The per-cell strategy is not a `LinearOperator`. It does not participate in `A + B` or `A @ B` at the cell level.
- Do NOT introduce a separate module for the residual math. It lives next to the `update` math — same module, same class.

### What to deliver

1. **Add `residual` to the `CellUpdate` Protocol** at `orpheus/sn/spatial/cell_update.py:418`:

   ```python
   @runtime_checkable
   class CellUpdate(Protocol):
       is_linear: bool
       is_positivity_preserving: bool

       def update(
           self,
           visit: CellVisit,
           total_xs: np.ndarray,
           source: np.ndarray,
           upstream_state: UpstreamState,
       ) -> CellResult: ...

       def residual(
           self,
           cell_avg: np.ndarray,
           visit: CellVisit,
           total_xs: np.ndarray,
           source: np.ndarray,
           upstream_state: UpstreamState,
       ) -> np.ndarray: ...
   ```

   Docstring on `residual` states: computes the per-cell discrete operator action `L_cell · ψ̄ − q`. At the converged cell average (i.e., when `cell_avg == update(...).cell_average_flux`), the residual is zero to FP rounding. This is the operator-direction companion to `update`'s solve-direction.

2. **Add `residual` as an abstractmethod on `CellUpdateBase`** at `cell_update.py:559`. Mirror the Protocol signature.

3. **Implement `DiamondDifference.residual`** at `orpheus/sn/spatial/diamond.py`. Three branches mirroring the existing three `_update_*` branches:
   - `_residual_slab(cell_avg, st, total_xs, source, upstream_state)` — closed-form `2|μ|·(cell_avg − ψ_in) + chord·Σ_t·cell_avg − source`.
   - `_residual_curvilinear(cell_avg, st, A_downstream, total_xs, source, upstream_state)` — consume `cell_balance_terms(st, A_downstream, total_xs, upstream_state)` and compute `terms.denom · cell_avg − (source · volume + terms.numer_upstream)`.
   - `_residual_cylindrical_degenerate(cell_avg, st, total_xs, source, upstream_state)` — consume `cell_balance_terms_degenerate(st, total_xs, upstream_state)` similarly.

   These three branches are direct ports of `_apply_curvilinear_residual` / `_apply_cylindrical_degenerate_residual` from the scratch buffer.

4. **Commit `cell_balance.py`** as part of the same step.

5. **Wire the Wave-D `_update_curvilinear` and `_update_cylindrical_degenerate` to consume `cell_balance_terms`**. The current `diamond.py` has duplicated algebra; the qa CONCERN-2 fix makes both `update` and `residual` consume the same helper. Operation-order must be preserved so the slab + Cartesian regression snapshots stay bit-identical.

6. **Port the salvageable tests** from the reverted `tests/sn/spatial/test_sncell_operator.py` to test `DiamondDifference.residual` directly (not a wrapper):
   - Apply-vs-solve consistency: `DiamondDifference.residual(cell_avg, ...) == 0` to `rtol=1e-12` when `cell_avg = DiamondDifference.update(...).cell_average_flux`. Cover slab + sphere + cylinder + cylindrical-degenerate. Heterogeneous + multi-group bias.
   - Linearity in `source` (the per-cell `residual` is affine in source; linear in `cell_avg`).
   - Bit-identity vs the previous `update` invariants (preserve them on the existing `update` API — no change).

   Put these in `tests/sn/spatial/test_diamond.py` (existing file) under a new `TestResidual` class. Do NOT create `test_sncell_operator.py`.

### Mechanism criteria (the agent's correctness contract)

- **The new method lives on the strategy**: `DiamondDifference.residual(...)`. Not `SNCellOperator.apply(...)`. Not a free function in a new module.
- **No new types are added.** `CellUpdate` Protocol gains one method. `CellUpdateBase` gains one abstractmethod. `DiamondDifference` gains three private branches.
- **`cell_balance_terms` is the single source of truth** for the curvilinear cell-balance algebra. Both `update` and `residual` consume it.
- **Bit-identity preserved on all 11 regression snapshots** at `rtol=1e-12` (because `update`'s operation order is preserved when refactoring to consume the helper; `residual` is new — does not affect existing solve results).

### Acceptance

- `tests/sn/spatial/test_diamond.py` passes with the new `TestResidual` class added (target: ~50 new tests, parametrized over geometry × n_groups × source kind).
- `tests/sn/regression/` 11/11 snapshots pass bit-identical (`rtol=1e-12`).
- `pytest tests/sn/ -q` regression: green-count delta is + new tests, no losses elsewhere.
- `mypy orpheus/sn/spatial/` (if installed) or `python -m mypy ...` passes — Protocol/abstractmethod added cleanly.

### Deliverable

Two commits:
- `feat(sn): CellUpdate.residual on the strategy layer (Issue #196 Step 1 replan)` — production code in `cell_update.py`, `diamond.py`, `cell_balance.py`, `__init__.py` export.
- `test(sn): DiamondDifference.residual + cell_balance bit-identity (Issue #196 Step 1 replan)` — test additions to `test_diamond.py`.

Closeout memo: `.claude/agent-memory/method-implementer/issue_196_phase_g_step1_replan_closeout.md`. Cite the rollback, list the salvaged math/test contracts, document the bit-identity preservation.

### Sub-agent dispatch

method-implementer with brief shaped per the Brief Template below.

---

## Step 2 — Canonical curvilinear SI sweep math (no GMRES)

**Goal**: fix the actual closure-math bug in `_sweep_1d_spherical` and `_sweep_1d_cylindrical`. The SI sweep, with `inner_solver="source_iteration"`, must run a Picard iteration on the canonical separate-pass M-M + WDD + face-flux BC operator. Issue #196 manifestation #7 closes by construction once both `transport_operator_matvec_*` (Krylov) and `_sweep_1d_*` (SI) consume the **same canonical operator definition**. Picard convergence on the corrected operator is to be empirically verified, not assumed.

This is the load-bearing step. It is also the step the previous agent shortcut. Treat it with extreme care.

### Pre-step: read these files in full

- `.claude/agent-memory/numerics-investigator/issue_196_phase_g_step2_diagnostic.md` — the H2 verdict + empirical numbers + closure attribution.
- `orpheus/sn/sweep.py` (after Step 0 revert) — `_sweep_1d_spherical`, `_sweep_1d_cylindrical`, the SI sweep dispatch.
- `orpheus/sn/operator.py:412-1107` — `transport_operator_matvec`, `transport_operator_matvec_spherical`, `transport_operator_matvec_cylindrical`. This is the apply-matvec form. Sections to identify: M-M angular recurrence as separate pass (per-level), WDD streaming closure, Carlson seed construction, BC trace via `bc_outer.apply(...)`.
- `orpheus/sn/spatial/pole_angular_closure.py:340-458` — `_mm_weighted_angular_recurrence_single_level`, the algorithmic core. Both the apply-matvec and the SI sweep should consume this directly.
- `orpheus/sn/spatial/psi_half_angle_seed.py:363-419` — `carlson_inward_sweep_from_source`, the shared Carlson seed kernel.
- `tests/sn/diagnostics/phase_g_step2_05_homogeneous.py` — the L0 streaming-equilibrium probe. Currently FAILS on SI (22% pole error), passes on Krylov (machine precision). This test is to be promoted to a permanent gate.

### Anti-recommendations (the agent MUST NOT do these)

- **Do NOT call `scipy.sparse.linalg.gmres` anywhere in `_sweep_1d_*`, `_curvilinear_sweep`, `transport_sweep`, or any code path consumed by `inner_solver="source_iteration"`.** This is the Cardinal Rule 1 violation the previous agent committed. The SI sweep is a Picard iteration on a sweep operator. If the Picard iteration on the canonical (corrected) operator does not converge, **STOP** and dispatch literature-researcher / numerics-investigator. Do not pivot to Krylov as a "principled" fallback.
- **Do NOT dismiss "Picard diverges" without empirical evidence.** The previous agent claimed this without showing the work. The correct test: implement the canonical operator, run Picard with `ψ_0 = 0` on the L0 streaming-equilibrium problem at n=40, plot `||ψ_k − ψ_{k-1}||_2 / ||ψ_k||_2` vs iteration index k. If the residual norm decays geometrically (factor < 1 per step), Picard converges. If not, characterise the spectral radius numerically and report back. Do NOT improvise the algorithm.
- **Do NOT delete `transport_operator_matvec_*` or `_sweep_1d_*`.** They become the canonical implementations. The previous agent's working tree appears to have removed `_sweep_1d_spherical/cylindrical` from `sweep.py` — Step 0's revert must restore them. They are then **edited in place** with the canonical closure math, not replaced.
- **Do NOT create a `streaming.py` module.** L's home is to be decided in Step 3; Step 2 only fixes the existing sweep / matvec code in `orpheus/sn/sweep.py` and `orpheus/sn/operator.py`.
- **Do NOT create a `DiscreteOrdinatesPhaseSpace`.** SNMesh IS the phase space.

### The canonical closure (from the diagnostic)

For curvilinear geometry, the canonical per-cell math is:

1. **Source-driven Carlson seed** (or ψ-driven; equivalent at fixed point) for the half-angle flux at μ = −1: `carlson_inward_sweep_from_source(Q_bar, sigma_t, dr, bc_outer_value)` where `Q_bar = ½ · Σ_t · φ_0` (SI) or derived from current ψ iterate (Krylov). Per Hébert §3.9.4 Eqs. (3.432)-(3.435). The shared free function at `psi_half_angle_seed.py:363` is the canonical kernel.

2. **Separate-pass M-M angular redistribution**. Compute `R = MorelMontryAngularSweep(psi_cells, alpha_half, redist_dAw, tau_mm, volume, level_indices, carlson_context)` ONCE on the input ψ, not inside the per-cell update. `R` has shape `(ng, N, nx, ny)`. This is exactly what `_mm_weighted_angular_recurrence_single_level` already implements per level. The apply-matvec already consumes it this way; the SI sweep currently couples it inside the cell update — that's the bug.

3. **WDD streaming + collision** per cell using `_cell_balance_terms` (from Step 1). The per-cell balance is `denom · cell_avg = (source · volume + numer_upstream − R_{n,i,g} · volume)`. Note: the R term enters as a source contribution, not as a per-cell coefficient — this is what makes "separate-pass M-M" work.

4. **Face-flux BC trace**: `outer_inflow = bc_outer.apply(outflow_at_boundary_face)` — call `bc.apply` on the FACE flux, not on the cell-centre. The apply-matvec currently uses a cell-centre proxy (`fi[:, :, -1, 0]`); this is the O(h) BC trace artefact the diagnostic identified. The SI sweep uses face flux correctly. Both must use face flux after Step 2.

### What to deliver

1. **Rewrite `_sweep_1d_spherical`** (`orpheus/sn/sweep.py`) so the SI sweep math IS the canonical math above. Body shape:
   - Compute Carlson seed once per outer iteration: `carlson_seed = carlson_inward_sweep_from_source(Q_bar, sigma_t, dr, bc_outer.apply(psi_face))`.
   - Compute M-M redistribution once per outer iteration as a separate pass over the current `psi_cells`: `R = pole_angular_closure(psi_cells, alpha_half, redist_dAw, tau_mm, volume, carlson_context=...)`.
   - Inside the per-ordinate / per-cell loop, run WDD streaming + collision using `DiamondDifference.update(visit, total_xs, source - R[n,i,g]/V_i, upstream_state)`. The M-M contribution is folded into the per-cell source — it's not a separate term in the cell update body.
   - BC application: face-flux. `bc_outer.apply(psi_face_outflow)` returns `psi_face_inflow` per Wave 9 contract.

2. **Rewrite `_sweep_1d_cylindrical`** analogously, with per-μ-level Carlson context.

3. **Rewrite `transport_operator_matvec_spherical` and `transport_operator_matvec_cylindrical`** (`orpheus/sn/operator.py:571-1107`) to use the same canonical math: same separate-pass M-M, same face-flux BC trace. After this, the apply-matvec path and the SI sweep path execute the same closure — manifestation #7 closes by construction.

4. **Promote the L0 streaming-equilibrium probe** to a permanent gate at `tests/sn/spatial/test_streaming_equilibrium_curvilinear.py`. Decorators: `@pytest.mark.l0 @pytest.mark.catches("ERR-026")` (and a new ERR-NNN entry — see point 6 below). Verifies: SI sweep AND Krylov apply-matvec BOTH give φ = 10.0 everywhere, ψ_n = 5.0 per ordinate, to `rtol=1e-9`. Promoted from `tests/sn/diagnostics/phase_g_step2_05_homogeneous.py`.

5. **Empirical Picard convergence study**. Either as a permanent test or as a one-off characterisation captured in the closeout memo: run SI on the L0 streaming-equilibrium problem and on `sphere_2g_3reg` n=40 at the canonical closure, log `||ψ_k − ψ_{k-1}||_2 / ||ψ_k||_2` per iteration. Document the empirical spectral radius. This is the evidence the previous agent did not produce.

6. **Add ERR-NNN entry to `.claude/skills/vv-principles/error_catalog.md`** for the L0 manifestation:
   - ERR-NNN (next sequential ID; check the current tail of the catalog)
   - Failure mode: #4 (wrong recurrence — M-M angular state propagation INSIDE the per-cell update breaks per-ordinate isotropy at the spherical pole)
   - How it hid: Phase F's diagnostic compared SI-vs-Krylov k_eff (which agreed to 0.286% at n=40) but never tested the L0 streaming-equilibrium limit. The 6 curvilinear snapshots were SI-generated under the bug, so bit-identity preserved the bug.
   - Which test catches it: `tests/sn/spatial/test_streaming_equilibrium_curvilinear.py`.
   - Lesson: in curvilinear SN, the M-M angular redistribution MUST be a separate pass, not inside the per-cell update. Coupling it inside the cell update breaks per-ordinate inversion at the pole.

7. **Update ERR-026 manifestation table** in the same catalog file:
   - #6: PARTIAL CLOSURE → CLOSED (Phase G Step 2 finishes Phase F's partial closure).
   - #7: OPEN → CLOSED (manifestation #7 by construction — same operator on both paths).
   - New L0 entry (manifestation #8): CLOSED via the L0 test.

8. **Regenerate the 6 curvilinear regression snapshots** under the corrected closure:
   - `sphere_2g_homogeneous_dd_n20.npz`, `sphere_2g_3reg_dd_n40.npz`, `sphere_2g_p1_aniso_dd_n20.npz`
   - `cyl_1g_homogeneous_LS4_dd_n20.npz`, `cyl_1g_homogeneous_product_dd_n20.npz`, `cyl_2g_3reg_LS4_dd_n40.npz`

   **Principled-equivalence justification** (per `vv-principles` §"Bit-identity vs principled-equivalence"):
   - Principled new formulation: each intermediate is named (separate-pass `R`, WDD-closure `cell_balance_terms`, face-flux BC trace).
   - Three-pillar structurally-independent reference: (a) L0 streaming-equilibrium analytical answer, (b) Pomraning 1989 pole isotropy, (c) Variant α integral form via `tests/sn/test_phase_c_crosscheck.py::test_phase_d_trajectory_resolvent_crosscheck` at Phase E rtols.

9. **Remove the Phase E sentinel's xfail-strict marker** at `tests/sn/test_phase_c_crosscheck.py::test_phase_e_trajectory_resolvent_flux_shape_crosscheck`. The test is predicted XPASS after Step 2. If it does NOT xpass, the closure fix is incomplete — investigate, don't improvise.

### Mechanism criteria

- **`grep -rn "scipy.sparse.linalg.gmres" orpheus/sn/sweep.py orpheus/sn/operator.py`** returns NOTHING after Step 2. GMRES is only consumed in `_solve_krylov` and `_solve_fixed_source_krylov` (the explicit Krylov inner solvers).
- **The L0 streaming-equilibrium test passes** on both `inner_solver="source_iteration"` AND `inner_solver="krylov"` at `rtol=1e-9` for the homogeneous reflective sphere with mixture B and isotropic Q=1 at n_cells ∈ {20, 40, 80}.
- **Both apply-matvec and SI sweep call `_mm_weighted_angular_recurrence_single_level`** as a separate pass, NOT inside the per-cell update.
- **Both apply-matvec and SI sweep call `bc_outer.apply(...)` on a face flux**, not on a cell-centre proxy.
- **`pytest tests/sn/spatial/test_sweep_vs_apply_consistency.py -v`** — every existing test continues to pass; the previously-xfail twin-path defense test now XPASSES.

### Decision-point checkpoints

If you (the agent) reach for any of the following, STOP and dispatch:

- "Picard appears to diverge" → STOP. Capture the residual-norm history; dispatch numerics-investigator with the residual history and the operator-action characterisation. The previous agent stopped here with "GMRES is fine" — do not.
- "The face-flux BC trace produces NaN at the boundary" → STOP. Likely a sign error or a face-vs-cell-centre confusion. Dispatch the diagnostic memo at `numerics-investigator/issue_196_phase_g_step2_diagnostic.md` outputs A/D for re-reading.
- "The L0 streaming-equilibrium test still shows 22% error on SI after the fix" → STOP. The M-M redistribution is still inside the per-cell update somewhere; trace and remove.
- "The 6 curvilinear snapshots regenerate to values that disagree with Variant α by more than the Phase E tolerance" → STOP. The canonical closure is not yet correct; do not regen the snapshots until the three-pillar attestation holds.

### Acceptance

- L0 streaming-equilibrium test passes at machine precision on both solvers.
- Phase E sentinel XPASSES; marker removed.
- ERR-026 manifestation table closes #6, #7, and a new #8.
- 6 curvilinear snapshots regenerate; 5 Cartesian snapshots stay bit-identical.
- All `tests/sn/` tests pass (with the regenerated snapshots).
- New ERR-NNN entry added; closeout memo cites three-pillar attestation.

### Deliverable

Atomic commit chain:
- `fix(sn): canonical curvilinear closure — separate-pass M-M + face-flux BC trace (Issue #196 Step 2 replan)` — the math fix in `_sweep_1d_*` and `transport_operator_matvec_*`.
- `test(sn): promote L0 streaming-equilibrium probe; Phase E sentinel re-enabled (Issue #196 Step 2 replan)`.
- `chore(sn): regenerate 6 curvilinear regression snapshots under canonical closure (Issue #196 Step 2 replan)`. Body documents the three-pillar attestation.
- `docs(vv): ERR-NNN L0 pole-closure manifestation; ERR-026 #6+#7+#8 CLOSED (Issue #196 Step 2 replan)`.

Closeout memo: `.claude/agent-memory/method-implementer/issue_196_phase_g_step2_replan_closeout.md`. Empirical evidence (Picard convergence history; three-pillar agreement; L0 test passes), files touched, what does NOT close (Step 3-6 still pending).

### Sub-agent dispatch

method-implementer with a brief sized to ONE step. The brief MUST explicitly forbid GMRES in the SI sweep path. If the implementer signals difficulty with Picard convergence, dispatch a numerics-investigator follow-up before proceeding (do NOT switch to GMRES).

---

## Step 3 — L pure-streaming + C separate, built on `SNMesh`

**Goal**: properly construct the four-operator algebra's first two leaves. `L = StreamingOperator(SNMesh, *, boundary)` with NO σ_t. `C = CollisionOperator(SNMesh, sigma_t)`. `(L + C).solve(q)` is the canonical curvilinear sweep (delegating to Step 2's corrected `_sweep_1d_*`). No `DiscreteOrdinatesPhaseSpace`. No duplication with `SNStreamingOperator` (which Step 3 retires or aliases).

### Pre-step: read these files in full

- `orpheus/sn/operator.py:1115-1473` — `SNStreamingOperator(LinearOperatorMixin)`, the existing monolithic L. Note its capability set `{apply, solve, apply_transpose}`, its `apply` body dispatching to `transport_operator_matvec_*`, its `solve` body delegating to `transport_sweep`.
- `orpheus/numerics/operator.py:489-694` — `OperatorSum`, `OperatorProduct`, `ScaledOperator` capability propagation. Note: `(L+C).solve` does NOT propagate automatically from `OperatorSum` — `solve` does not distribute over sums. Either wire it explicitly on a new `FusedStreamingCollisionRepresentation`, or have `OperatorSum.solve` look for a `(L+C)` SN-specific fusion hook.
- `orpheus/sn/geometry.py:51-868` — `SNMesh` full surface. `BOUNDARY_OPERATOR_REGISTRY`, `_resolve_one`, `bc_left/right/xmin/...`.
- `orpheus/numerics/operator.py:1318-1368` — `as_scipy_linop` adapter for the Krylov path.

### Anti-recommendations

- **Do NOT create `DiscreteOrdinatesPhaseSpace`.** SNMesh IS the phase space. The four-operator constructors take `SNMesh` directly. `n_groups` can either be (a) a separate constructor parameter, (b) read from the cross-section data (`sigma_t.shape[-1]`), or (c) added to `SNMesh.__init__` as a new constructor param. Pick (b) if the operator can derive `n_groups` from the data it already takes; (a) if not.
- **Do NOT inherit from `SNStreamingOperator`.** Replace it. The legacy class is the monolithic `(L+C)` and is exactly what Step 3 splits.
- **Do NOT put L and C in a new `orpheus/sn/streaming.py` module.** The previous attempt created that module with `DiscreteOrdinatesPhaseSpace` and a broken `from .mesh import SNMesh` import; the rejected work is gone. Put them in `orpheus/sn/operator.py` next to where `SNStreamingOperator` lived, then delete `SNStreamingOperator` once the migration is verified.
- **Do NOT call GMRES inside `(L+C).solve`.** That's the canonical sweep — delegate to `transport_sweep` (which after Step 2 is the corrected SI sweep).

### What to deliver

1. **New `StreamingOperator(LinearOperatorMixin)`** in `orpheus/sn/operator.py`:
   - Constructor: `StreamingOperator(sn_mesh: SNMesh, *, boundary: LinearOperator | None = None)`. `boundary` defaults to `sn_mesh.bc_right` for curvilinear (or `bc_xmax`/etc. for Cartesian via a small selector helper).
   - `capabilities = frozenset({CAP_APPLY, CAP_APPLY_TRANSPOSE})`. No `CAP_SOLVE` — L alone is not invertible by the WDD sweep (only `L + C` is).
   - `apply(psi, *, sigma_t=None)` — the pure-streaming directional derivative. Routes through the geometry-dispatching `transport_operator_matvec_*` (post-Step-2 corrected) BUT with `sigma_t = 0` if `sigma_t` is None, OR the caller-supplied σ_t (this is how OperatorSum compositions can subtract C's contribution). The cleanest implementation: split `transport_operator_matvec_spherical` into two helpers, one for streaming + redistribution and one for adding collision; `L.apply` calls only the streaming helper.
   - `apply_transpose(psi)` — analogous, with reverse-direction sweep.

2. **New `CollisionOperator(LinearOperatorMixin)`** in the same module:
   - Constructor: `CollisionOperator(sn_mesh: SNMesh, sigma_t: np.ndarray)`. `sigma_t` shape `(nx, ny, ng)` matching `SNMesh` conventions.
   - `capabilities = frozenset({CAP_APPLY, CAP_SOLVE, CAP_APPLY_TRANSPOSE})`. C is self-adjoint (real diagonal); `apply_transpose = apply`.
   - `apply(psi)` — `Σ_t · ψ`. Shape-agnostic: handle both packed `(n_unknowns,)` and structured `(N, nx, ny, ng)` shapes via duck-typing.
   - `solve(q)` — `q / Σ_t` (per-element). Reject σ_t = 0 cells with a clear error.

3. **Wire `(L + C).solve`** so it routes through the canonical sweep. Two options:
   - Option A: add a `FusedStreamingCollisionRepresentation` class that `OperatorSum.__init__` constructs when it detects `(StreamingOperator, CollisionOperator)` operands; its `solve` delegates to `transport_sweep`.
   - Option B: leave `OperatorSum.solve` raising `MissingCapability` for the general case; provide a free function `solve_sn_within_group(L: StreamingOperator, C: CollisionOperator, q, ...) -> psi` that callers invoke explicitly. `SNSolver` calls this; `KEigenvalue`'s `inverter` hook (`iteration.py:209-229`) is given a closure over it.

   **Option B is simpler and more honest** (composition fusion is a special-case optimization that hides where the sweep lives). Pick Option B unless concrete reason for Option A emerges during implementation.

4. **Retire `SNStreamingOperator`**. After `StreamingOperator + CollisionOperator` are functional and consumed by `SNSolver` (Step 4), the legacy class becomes dead code. Delete it. `SNSolver.L = StreamingOperator(sn_mesh)`; `SNSolver.C = CollisionOperator(sn_mesh, sig_t)`.

5. **Extend `build_sn_operators(sn_mesh, materials) → (L, C, S, F)`** in `orpheus/sn/operator.py` (or wherever the existing builder lives; if no canonical builder, create one). S and F are constructed via their existing factories (`ScatteringOperator.from_solver_data`, `FissionOperator.from_solver_data`).

### Mechanism criteria

- **`StreamingOperator` has no `sigma_t` field** and `StreamingOperator.apply` never multiplies by σ_t.
- **`CollisionOperator.apply` works on slab, sphere, cylinder, and 2-D Cartesian** — fix the dormant `build_equation_map_spherical` bug. Geometry dispatch via `sn_mesh.curvature`.
- **`SNSolver.__init__` constructs `(L, C, S, F)` via `build_sn_operators(sn_mesh, materials)`.** No more `SNStreamingOperator` instance.
- **The 11 regression snapshots all stay bit-identical** because Step 3 is a code reorganisation, not a math change — Step 2 already changed the curvilinear math. If a snapshot breaks at Step 3, the split between L and C has introduced an arithmetic-order change; investigate.

### Acceptance

- `StreamingOperator + CollisionOperator` exist; `SNStreamingOperator` is gone.
- All 11 regression snapshots pass bit-identical.
- `tests/numerics/test_iteration.py::test_keigenvalue_matches_solve_sn_2g_slab` passes with the new operator triple (the adapter shims at `test_iteration.py:402-440` should be removable; if not, the shape mismatch is a hint that Step 3 isn't complete).
- `pytest tests/sn/ -q` green.

### Deliverable

Commit chain:
- `feat(sn): StreamingOperator (pure streaming) + CollisionOperator on SNMesh (Issue #196 Step 3 replan)`.
- `refactor(sn): retire SNStreamingOperator; SNSolver consumes (L, C, S, F) via build_sn_operators (Issue #196 Step 3 replan)`.
- `test(sn): direct LinearOperator coverage for L, C (Issue #196 Step 3 replan)`.

Closeout memo: `.claude/agent-memory/method-implementer/issue_196_phase_g_step3_closeout.md`.

### Sub-agent dispatch

method-implementer.

---

## Step 4 — Wire `SourceIteration` and `KEigenvalue` as canonical consumers

**Goal**: collapse the duplicated inner-solver code in `SNSolver`. `_solve_source_iteration` and `_solve_fixed_source_si` become one call site that invokes `SourceIteration(L, S, F, q_ext).solve(...)`. `_solve_krylov` and `_solve_fixed_source_krylov` collapse similarly. `power_iteration` is replaced by `KEigenvalue` (with the SN-aware `keff_estimator`).

### Pre-step: read these files in full

- `orpheus/numerics/iteration.py` — `SourceIteration` (line 164) and `KEigenvalue` (line 360). The signatures are already `(L, S, F)` + `q_ext`. The `inverter` hook (line 209) accepts an arbitrary inner-solve closure.
- `orpheus/sn/solver.py:281-625, 1067-1289` — the four duplicated solver functions.
- `tests/numerics/test_iteration.py:328-460` — `test_keigenvalue_matches_solve_sn_2g_slab`. Shows the bridge with adapter shims; Step 4 removes the shims.

### Anti-recommendations

- **Do NOT write a new `SourceIteration` or `KEigenvalue`.** They exist. Consume them.
- **Do NOT add a new `inner_solver` switch.** The existing string switch at `solver.py:162-170` stays as the dispatch; the bodies of both branches now consume the same iteration primitives, just with different `inverter` closures.
- **Do NOT write a `PreconditionedGMRES` LinearOperator wrapper** before checking the inline use is the only consumer. There is exactly one Krylov inner-solve call site to be collapsed (`_solve_krylov` and `_solve_fixed_source_krylov` become one). Use inline `gmres(as_scipy_linop(A_loss), q, M=as_scipy_linop((L+C).inverse_via_sweep))` — write the wrapper only if a second consumer appears.

### What to deliver

1. **Collapse `_solve_source_iteration` and `_solve_fixed_source_si`** to one helper `_si_solve(L, C, S, F, q_ext, ...)` that builds and runs `SourceIteration(L, S, F, q_ext)` with `inverter=lambda q: solve_sn_within_group(L, C, q, ...)`. The two existing functions become one-line wrappers that compose the right operands and call `_si_solve`.

2. **Collapse `_solve_krylov` and `_solve_fixed_source_krylov`** to one helper `_gmres_solve(L, C, S, F, q_ext, ...)` that constructs `A_loss = L + C - S`, wraps with `as_scipy_linop`, and calls `scipy.sparse.linalg.gmres`. Use `(L + C).solve_via_sweep` (or whatever Step 3 named the canonical sweep) as the left preconditioner. The two existing functions become one-line wrappers.

3. **Replace `power_iteration` in `solve_sn`** with `KEigenvalue(L, S, F, inverter=..., keff_estimator=lambda L_,S_,F_,phi: solver.compute_keff(phi))`. `power_iteration` continues to ship for the cross-solver migration (CP, diffusion, MoC, homogeneous) per its `DeprecationWarning`.

4. **Use `as_scipy_linop`** at the two Krylov sites instead of the bare `ScipyLinearOperator((n,n), matvec=...)` lines. One-line swap per site; free cleanup.

5. **Remove the adapter shims** from `test_keigenvalue_matches_solve_sn_2g_slab` (test_iteration.py:402-440) once the shape mismatches are gone.

### Mechanism criteria

- `_solve_source_iteration` body is ≤ 10 lines (the rest is the helper).
- `_solve_krylov` body is ≤ 10 lines.
- `solve_sn` calls `KEigenvalue(...).solve(...)`, not `power_iteration(...)`.
- All `tests/sn/regression/` and `tests/numerics/test_iteration.py` tests pass.
- `power_iteration` is no longer called from `orpheus/sn/`; the deprecation warning issued at import is for cross-solver migration only.

### Acceptance

- `pytest tests/sn/ tests/numerics/ -q` green.
- Line counts in `solver.py` drop by ~250 lines (the two ~125-line duplicated solver functions become 10-line wrappers).

### Deliverable

Commit chain:
- `refactor(sn): SourceIteration / KEigenvalue consumers; _si_solve and _gmres_solve helpers (Issue #196 Step 4)`.
- `refactor(sn): solve_sn uses KEigenvalue; power_iteration retired from SN paths (Issue #196 Step 4)`.

Closeout memo: `.claude/agent-memory/method-implementer/issue_196_phase_g_step4_closeout.md`.

### Sub-agent dispatch

method-implementer.

---

## Step 5 — BC verification + minor cleanup

**Goal**: confirm BCs are already first-class `LinearOperator`s composed at the trace-edge layer (Wave 9/11 work), and that the BC composition into L's apply / sweep is correct after Steps 1-4. Minor cleanup of dormant deprecation warnings.

### What to deliver

1. **Read `orpheus/sn/boundary_realizer.py:106-200`** and verify `SNBoundaryRealizer.realize` dispatch is complete (vacuum, specular, white, albedo, periodic, prescribed_inflow). Confirm `sn_mesh.bc_*` exposes realised `LinearOperator`s with the 1-arg `apply` contract.

2. **Verify the BC composition is at the trace-edge layer** in the post-Step-2 `_sweep_1d_*` and `transport_operator_matvec_*` paths. Per-ordinate `bc.apply` slicing inside the per-cell loop should no longer exist (Wave 11 mixed removal landed). If any per-ordinate slicing remains, factor it out.

3. **Confirm Wave 11 `bc_xmin/xmax/ymin/ymax` cleanup is complete**. Any deprecation warnings issued by SNMesh's BC accessors that should now be removed.

### Mechanism criteria

- `grep -rn "for n in range" orpheus/sn/sweep.py orpheus/sn/operator.py | grep "bc"` returns nothing — no per-ordinate BC slicing inside the sweep body.
- `pytest tests/sn/test_boundary_conditions.py tests/sn/test_snmesh_realizer_wiring.py` green.

### Acceptance

- BC tests green; no per-ordinate `bc.apply` calls in the sweep / matvec bodies.
- One-paragraph confirmation in the closeout memo that the Wave-9/11 BC architecture is intact after Phase G.

### Deliverable

If cleanup is needed, single commit `refactor(sn): clean up BC composition at trace-edge layer (Issue #196 Step 5)`. If not, the closeout memo confirms no change needed.

### Sub-agent dispatch

qa with a focused "BC architectural audit" brief.

---

## Step 6 — `.H` on leaves + `solve_sn_adjoint`

**Goal**: the elegance acceptance criterion. `solve_sn_adjoint` writes itself in 6 lines because `.H` propagates through `OperatorSum` to the leaves. Implement `.H` on each leaf.

### Pre-step: read these files in full

- `orpheus/numerics/operator.py:415-486` — `_AdjointOperator`. Its `apply` does `(1/w_V) · A.apply_transpose(w_W · y)`. This means a leaf advertising `CAP_APPLY_TRANSPOSE` gets `.H` for free.
- `orpheus/sn/scattering.py:317-319` — `ScatteringOperator.capabilities = frozenset({CAP_APPLY})`. No `apply_transpose` yet. Step 6 adds it.
- `orpheus/sn/fission.py:120-122` — same for `FissionOperator`. Step 6 adds it.
- `orpheus/sn/operator.py:1442-1472` — `SNStreamingOperator.apply_transpose` (was implemented by dense probing; will be retired with `SNStreamingOperator`). Step 6 ensures the new `StreamingOperator` has an analytic transpose (reverse-direction sweep).

### Anti-recommendations

- **Do NOT implement `apply_transpose` as a dense-probe of `apply`.** That's O(n²) and was a placeholder in the legacy `SNStreamingOperator`. The analytic transpose for `StreamingOperator` is the reverse-direction sweep: swap `iter_cells_by_direction(+1)` for `(-1)` and vice versa. For `CollisionOperator` it's identity (self-adjoint). For `ScatteringOperator` it's `R.H @ Λ.H @ M.H` (each piece already has `apply_transpose` per the `GalerkinProjection` Protocol).
- **Do NOT add `apply_transpose` to operators that don't have a clear mathematical transpose**. The "harmful stub" anti-pattern (raises `NotImplementedError`) is documented in `docs/theory/operator_algebra.rst:149-156`. If a leaf isn't ready, omit `CAP_APPLY_TRANSPOSE` and `S.H` correctly raises `MissingCapability` at composition.

### What to deliver

1. **`StreamingOperator.apply_transpose`** — reverse-direction sweep. Add `CAP_APPLY_TRANSPOSE` to capabilities.

2. **`CollisionOperator` already advertises `CAP_APPLY_TRANSPOSE`** with `apply_transpose = apply` (self-adjoint, diagonal). Verify after Step 3.

3. **`ScatteringOperator.apply_transpose`** — `R.H @ Λ.H @ M.H` composed. Implement and advertise `CAP_APPLY_TRANSPOSE`. The Galerkin projection's `M.H` and `R.H` already exist (`projection.py`); `Λ.H` is the transpose of the per-ℓ block-diagonal scattering matrix.

4. **`FissionOperator.apply_transpose`** — for `F = χ ⊗ (νΣ_f)^T`, `F^T = (νΣ_f) ⊗ χ^T`. Swap the rank-1 factors.

5. **`solve_sn_adjoint(materials, mesh, quadrature, response, *, tol=1e-8)`** in `orpheus/sn/solver.py` (or a new `orpheus/sn/adjoint.py` if appropriate). Body shape:

   ```python
   def solve_sn_adjoint(materials, mesh, quadrature, response, *, tol=1e-8):
       sn_mesh = SNMesh(mesh, quadrature)
       L, C, S, F = build_sn_operators(sn_mesh, materials)
       A_loss = L + C - S
       q_dagger = response.as_source()
       A_loss_H_scipy = as_scipy_linop(A_loss.H, shape=...)
       psi_dagger, _ = scipy.sparse.linalg.gmres(A_loss_H_scipy, q_dagger, rtol=tol)
       return psi_dagger
   ```

   ≤ 8 lines. The `.H` propagates through `OperatorSum` to the leaves; each leaf's `apply_transpose` does the rest.

6. **Reciprocity tests** in a new `tests/sn/test_adjoint.py`:
   - 2G slab reciprocity: `<response, ψ_forward> = <ψ_adjoint, source>` at `rtol=1e-10` for a homogeneous reflective slab.
   - 2G heterogeneous slab reciprocity.
   - Sphere reciprocity (curvilinear; verifies `L.H` against the forward sweep on a self-adjoint problem like isotropic scattering).
   - For a self-adjoint problem (homogeneous, P0, vacuum BC), the adjoint flux should equal the forward flux at the same eigenvalue. Test this as a sanity check.

### Mechanism criteria

- `(L + C - S).H` is constructed via `OperatorSum.H` and `_AdjointOperator` — no manual transpose composition.
- Each leaf operator advertises `CAP_APPLY_TRANSPOSE` only if its `apply_transpose` is the **analytic** transpose, not a dense probe.
- `solve_sn_adjoint` body is ≤ 8 lines.

### Acceptance

- `tests/sn/test_adjoint.py` ≥ 4 cases passing.
- The "adjoint reads like math" sphinx narrative renders.

### Deliverable

Commit chain:
- `feat(sn): .H on StreamingOperator, ScatteringOperator, FissionOperator leaves (Issue #196 Step 6)`.
- `feat(sn): solve_sn_adjoint via A_loss.H.solve (Issue #196 Step 6)`.
- `test(sn): adjoint reciprocity + self-adjoint sanity (Issue #196 Step 6)`.

Closeout memo: `.claude/agent-memory/method-implementer/issue_196_phase_g_step6_closeout.md`.

### Sub-agent dispatch

test-architect for the adjoint reciprocity test design (pre-step), then method-implementer for the leaf `.H` implementations and `solve_sn_adjoint`.

---

## Step S — Sphinx narrative (parallel from Step 2 onward)

**Goal**: extend `docs/theory/operator_algebra.rst` with the Phase G unification chapter. Wire `verifies(...)` decorators on tests to existing `:label:` references in the operator-algebra Sphinx page. Add a `docs/theory/sn_adjoint.rst` (or extend `discrete_ordinates.rst`) for Step 6's adjoint theory.

### Anti-recommendations

- **Do NOT rewrite the existing 796-line `operator_algebra.rst` stub.** It already documents the algebra, capability semantics, composition closures, harmonic factoring. Extend it with a Phase G unification section near the end; cross-reference Phase F's narrative.
- **Do NOT create equation labels that duplicate existing ones.** Use the existing labels (`operator-fixed-source`, `operator-eigenvalue`, `operator-apply`, `wdd-closure`, `mm-weights`, `hebert-3-432`, ...).

### What to deliver

1. **Phase G unification section in `docs/theory/operator_algebra.rst`** — narrative on the four-operator algebra application to SN: how L, C, S, F compose, how the SI sweep and Krylov apply consume the same operator, the L0 streaming-equilibrium gate, the `.H` propagation.

2. **Wire `verifies(...)` decorators** on all new Phase G tests against the existing `:label:` references. The Phase G replan agent must NOT create new labels — use the ones in `operator_algebra.rst` and `discrete_ordinates.rst`.

3. **Sphinx `-W` builds clean** at each step's commit. If a Sphinx warning is introduced, fix it before the commit.

### Deliverable

One or more commits over the course of the replan, owned by the Archivist agent dispatched in parallel from Step 2 onward. Final commit `docs(sn): Phase G unification narrative; verifies decorators wired (Issue #196)`.

### Sub-agent dispatch

archivist.

---

## Brief template for sub-agent dispatch

Every sub-agent brief for this replan MUST follow this template. The previous attempts failed because briefs were too long, too bundled, and lacked mechanism criteria. Use this template verbatim; fill the per-step body.

```
## Background (mandatory reads, in order)

- This plan: `.claude/plans/crystalline-wondering-token.md` §"Step N — ..."
- Predecessor closeout: `.claude/agent-memory/.../issue_196_phase_g_step{N-1}_closeout.md`
- The three exploration memos: <paths>
- Specific source files: <list with paths AND line ranges>

## Anti-recommendations (the agent MUST NOT do these)

- Do NOT create <existing-type>. <Existing-type> at <path:line> is the existing
  one. Extend it.
- Do NOT <specific failure mode that occurred in the prior attempt>. <Why>.
- Do NOT <another specific>.

(At least 5 anti-recommendations per brief, each naming the existing type with
file:line and stating why the alternative is wrong.)

## Mechanism criteria (the correctness contract)

- The new <method/type/file> lives at <path>. Not at <wrong-path>.
- <grep-able invariant> returns <expected> after the work is done.
- <existing-test> passes at <tolerance>.

(Mechanism criteria, not outcome criteria. "The L0 test passes" is outcome.
 "The SI sweep calls _mm_weighted_angular_recurrence_single_level as a separate
 pass, NOT inside the per-cell update" is mechanism.)

## Scope hard limits

IN scope:
- <bullet>
- <bullet>

OUT of scope (deferred):
- <bullet> → Step <M> handles this.
- <bullet> → Tracked in Issue <N>.

## Decision-point checkpoints

If you reach for any of the following, STOP and dispatch / ask:
- "<rationalisation prose>" → STOP. <action>.
- "<observed anomaly>" → STOP. Dispatch <which agent> with <what payload>.

## Test pin

- This step must turn <test path> from <state> to <state>.
- This step must keep <test paths> at <state> with <tolerance>.

## Commit scope

This step ships <N> commits:
- `<type>(<scope>): <subject> (Issue #196 Step <K>)` — <body>
- ...

## Closeout memo

Path: `.claude/agent-memory/.../issue_196_phase_g_step<K>_closeout.md`.

Required content:
- Empirical evidence (numbers, not prose)
- Files touched
- What this does NOT close
- Next step pointer
```

The previous briefs averaged 200 lines and bundled 7 deliverables per step. The corrected briefs should average 60-100 lines and ship one focused deliverable per dispatch.

---

## Verification (Phase G closure gate)

After Step 6 lands, run the following sequence. Any failure means Phase G has NOT closed; investigate before declaring done.

```bash
# Step 0 acceptance
git log --oneline | head -10  # No SNCellOperator commits; revert commits visible

# Step 1 acceptance
.venv/bin/python -m pytest tests/sn/spatial/test_diamond.py -v

# Step 2 acceptance — the load-bearing test
.venv/bin/python -m pytest tests/sn/spatial/test_streaming_equilibrium_curvilinear.py -v
.venv/bin/python -m pytest tests/sn/test_phase_c_crosscheck.py::test_phase_e_trajectory_resolvent_flux_shape_crosscheck -v  # XPASS

# Step 3-4 acceptance
.venv/bin/python -m pytest tests/numerics/test_iteration.py -v
.venv/bin/python -m pytest tests/sn/regression/ -v

# Step 6 acceptance
.venv/bin/python -m pytest tests/sn/test_adjoint.py -v

# Phase G end-to-end
.venv/bin/python -m pytest tests/sn/ tests/numerics/ -v
.venv/bin/python -m sphinx -W docs docs/_build/html

# ERR-026 manifestation table closure
grep -A3 "ERR-026" .claude/skills/vv-principles/error_catalog.md | head -40
# Should show: #5 OPEN (in #195), #6 CLOSED, #7 CLOSED, #8 CLOSED

# Issue #196 closes
gh issue view 196  # Status: OPEN, closes on merge per commit trailer
```

---

## Files this replan will modify (canonical list)

**Strategy layer (Step 1)**:
- `orpheus/sn/spatial/cell_update.py` — add `residual` to Protocol + ABC.
- `orpheus/sn/spatial/diamond.py` — implement `residual` in three branches.
- `orpheus/sn/spatial/cell_balance.py` — commit the existing untracked helper.
- `orpheus/sn/spatial/__init__.py` — export `cell_balance_terms` if needed.
- `tests/sn/spatial/test_diamond.py` — port salvageable invariants as a `TestResidual` class.

**Canonical closure (Step 2)**:
- `orpheus/sn/sweep.py` — rewrite `_sweep_1d_spherical`, `_sweep_1d_cylindrical`.
- `orpheus/sn/operator.py` — rewrite `transport_operator_matvec_spherical`, `transport_operator_matvec_cylindrical`.
- `tests/sn/spatial/test_streaming_equilibrium_curvilinear.py` — promoted L0 gate (new file).
- `tests/sn/test_phase_c_crosscheck.py` — remove `xfail-strict=True` from Phase E sentinel.
- `tests/sn/regression/snapshots/sphere_2g_*.npz`, `cyl_*.npz` — regenerate (6 files).
- `.claude/skills/vv-principles/error_catalog.md` — ERR-NNN entry + ERR-026 table update.

**L, C, build_sn_operators (Step 3)**:
- `orpheus/sn/operator.py` — add `StreamingOperator`, `CollisionOperator`; delete `SNStreamingOperator`.
- `orpheus/sn/solver.py:240-262` — `SNSolver.__init__` uses `build_sn_operators(sn_mesh, materials)`.
- `tests/sn/test_snstreamingoperator.py` — rename or adapt to test the new `StreamingOperator + CollisionOperator` triple.

**Iteration consumers (Step 4)**:
- `orpheus/sn/solver.py:281-1289` — collapse the four duplicated solver functions to two helpers.
- `tests/numerics/test_iteration.py:402-440` — remove adapter shims.

**Adjoint (Step 6)**:
- `orpheus/sn/operator.py` — `StreamingOperator.apply_transpose`.
- `orpheus/sn/scattering.py` — `ScatteringOperator.apply_transpose` + `CAP_APPLY_TRANSPOSE`.
- `orpheus/sn/fission.py` — `FissionOperator.apply_transpose` + `CAP_APPLY_TRANSPOSE`.
- `orpheus/sn/solver.py` or new `orpheus/sn/adjoint.py` — `solve_sn_adjoint`.
- `tests/sn/test_adjoint.py` — reciprocity tests (new file).

**Sphinx (Step S)**:
- `docs/theory/operator_algebra.rst` — Phase G unification section.
- `docs/theory/sn_adjoint.rst` (new) or `docs/theory/discrete_ordinates.rst` extension.

---

## Out of scope (tracked separately)

- **Issue #195** — L1 curvilinear MMS magnitude at n=160. ERR-026 manifestation #5. Separate investigation; not closed by this replan.
- **Issue #197** — NewType-based unit/quantity discipline (Option C). Adjacent improvement; can land independently.
- **MoC / CP / diffusion four-operator extension** — out of scope per `feedback_unify_after_two_instances.md`. SN is instance 1. Issue #11 (GMRES for SN) and Issue #161 (unified sweep) merge into this replan; their close conditions are achieved by Step 4.
- **Phase D / Phase E narrative cleanup** — already in `docs/theory/discrete_ordinates.rst`. Step S adds Phase G; does not rewrite Phase D-F.
- **Plan note from `[[feedback-check-literature-folder-first]]`** — if Step 2 dispatches literature-researcher, the agent must `ls scratch/literature/` before concluding paywall.

---

## Final note for the agent picking this up

This plan is your **decision instrument**. Sections that look like checklists (anti-recommendations, mechanism criteria, decision-point checkpoints) are the load-bearing parts — not the narrative prose. If a sub-agent's brief is being drafted and ANY anti-recommendation from the relevant section is missing, the brief is incomplete.

Three meta-rules learned from the failed attempts:

1. **Name existing types in the brief.** "SNMesh IS the SN phase space" beats "we need a phase-space type". The agent only knows what's in its brief plus what it explores; if you don't name SNMesh, it invents.

2. **State mechanism criteria, not outcome criteria.** "The L0 test passes" is testable by routing through Krylov. "The SI sweep calls Picard, not GMRES" forces the right architecture.

3. **One deliverable per dispatch.** A brief that asks for 7 things gets 5 done well and 2 done sloppily. Split.
