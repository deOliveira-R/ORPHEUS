# Phase G replan — SN four-operator unification (Issue #196)

**Tracking issue**: Issue [#196](https://github.com/deOliveira-R/ORPHEUS/issues/196) — ERR-026 manifestation #7 (SI-vs-Krylov O(h) WDD spatial-closure asymmetry on curvilinear SN). The original Phase G plan (`.claude/plans/issue_196_phase_g_four_operator_unification.md`) is **superseded** by this replan.

**Branch**: `refactor/sn-operator-algebra` (continued from Phase F tip `b0cc1b1`; replan commits start at `a4d4d8a`).

---

## STATUS (last update: 2026-05-14)

**Steps 0, 1, 2, 2.5, 2.5b COMPLETED + L12 fixes landed. Step 2.5c DESIGN LOCKED (implementation pending). Step 2.6 reduced to Q3 only (Q1 subsumed by 2.5c). Steps 3-6 + S + acceptance gate pending.**

| Step | Status | Key commits | Closeout / design |
|------|--------|-------------|-------------------|
| 0 — Rollback mis-abstractions | DONE 2026-05-13 | `31abf4a`, `4ec8d31`, `a4d4d8a` | (revert + replan + SUPERSEDED note on prior closeout) |
| 1 — `CellUpdate.residual` strategy layer | DONE 2026-05-13 | `562b7df`, `e2570fd`, `59f3604` | `.claude/agent-memory/method-implementer/issue_196_phase_g_step1_replan_closeout.md` |
| 2 Path C (sphere) — surgical pole-face IC + Carlson seed | DONE 2026-05-13 | `0281ab4`, `868c2fd`, `e5091a6`, `a71f211` | `.claude/agent-memory/method-implementer/issue_196_phase_g_step2_path_c_closeout.md` |
| 2 Path C (cylinder) — `1/Σw` normalisation | DONE 2026-05-13 | `1e02ca5`, `e2856a7`, `b42be78`, `a52c3db` | `.claude/agent-memory/method-implementer/issue_196_phase_g_step2_cylinder_fix_closeout.md` |
| 2.5 — Unified sweep skeleton + polymorphic DD | DONE 2026-05-13 | `23766ee`, `2ed4f15`, `fa0fa2c`, `abd0383`, `92b527b`, `1a699b0` | `.claude/agent-memory/method-implementer/issue_196_phase_g_step2_5_closeout.md` (+ Addendum L12 verification) |
| 2.5 follow-up — L12-caught regressions | DONE 2026-05-14 | `93117df` | (transport_sweep dispatch on Mesh2D(ny=1) + retire stale _solve_recurrence test imports) |
| 2.5b — `ordinate_scan` Blelloch primitive | DONE 2026-05-14 | `f5b913c`, `c034372`, `defa1e4`, `b86a9e4` | `.claude/agent-memory/method-implementer/issue_196_phase_g_step2_5b_closeout.md` |
| **2.5c — Two-stratum precomputed cache** | **DESIGN LOCKED 2026-05-14 — IMPLEMENTATION PENDING** | — | Plan §"Step 2.5c" + memos: `.claude/agent-memory/cross-domain-attacker/issue_196_phase_g_step2_5c_native_expression.md`, `.claude/agent-memory/explorer/issue_196_phase_g_step2_5c_immutability_audit.md` |
| 2.6 — `dag_walk` canonicalisation (Q3 only; Q1 SUBSUMED by 2.5c) | DESIGN LOCKED | — | Plan §"Step 2.6"; further-collapse memo `.claude/agent-memory/explorer/issue_196_phase_g_step2_5_further_collapse.md` |
| 3 — L pure-streaming + C separate on SNMesh | PENDING | — | — |
| 4 — `SourceIteration` / `KEigenvalue` consumers | PENDING | — | — |
| 5 — BC verification + cleanup | PENDING | — | — |
| 6 — `.H` on leaves + `solve_sn_adjoint` | PENDING | — | — |
| S — Sphinx Phase G unification narrative | PENDING (parallel from Step 3 onward) | — | — |

**Independent main-agent verification on Step 2.5b end-state (2026-05-14):**
- 52/52 PASS on `test_ordinate_scan.py` in 0.42s — pair-monoid theorem at rtol=1e-14 confirmed
- 26/26 PASS on L0 streaming-equilibrium curvilinear (32:32 wall)
- 113 PASS on `test_diamond.py + test_ordinate_scan.py + test_sweep_regression.py` in 7.22s
- 81 PASS on `test_sweep_vs_apply_consistency.py + test_apply_matvec_cylinder_invariants.py` in 6:16
- Full `pytest tests/sn/ -q`: projected ~33 min (killed earlier 3hr-run was actually progressing, just slow — ~28% in 11:28 = ~33 min total at that rate)

**Performance regression flagged 2026-05-14:** Step 2.5b's `ordinate_scan` is correct (1% of sweep time) but the per-cell Python iteration around it (`iter_cell_visits` × `streaming_terms` dataclass × `affine_coefficients` `np.fromiter`) accounts for 80% of sweep time. Pre-Step-2.5 slab cumprod path was vectorised across all cells; Step 2.5b replaced it with per-cell fold. Slab benchmark `nx=160 N=16 ng=4`: 30.85 (pre-2.5) → 15.43 ms/sweep (post-2.5b). **Step 2.5c restores cumprod-era speed** via the two-stratum precomputed cache — target ≤ 1.5 ms/sweep + full SN suite < 5 min. cached_property landed on `Mesh1D.widths` / `centers` (~1.6% impact, kept) but the deeper structural fix is Step 2.5c.

**Independent main-agent verification on Step 2 end-state (2026-05-13):** 39 passed, 0 failed in 49 min (24 L0 streaming-equilibrium cases sphere+cylinder; 11 regression snapshots; 2 Phase E sentinel cases; Pomraning isotropy). ERR-026 manifestations #6 and #7 CLOSED. ERR-048 added (3 manifestations: pole-face IC, Carlson seed source convention, cylinder Σw normalisation).

**The 3-way standoff for correctness (user's verbatim mandate, 2026-05-13):**
> "There is a 3-way stand-off. Any of the algorithms must be fixed to match the reference, then the algorithms must match each other, so that in the end both algorithms (krylov and sweep) match reference correctness and each other. They must also show correct behaviour under refinement."

Concretely four legs: (1) Krylov ≡ reference, (2) SI sweep ≡ reference, (3) Krylov ≡ SI, (4) correct convergence under refinement on each. Matching each other is necessary but NOT sufficient — both algorithms can be equally wrong (cylinder pre-fix Krylov and SI agreed within ~1% but disagreed with analytical truth by 31%). Reference must be structurally independent (analytical / MMS / Variant α).

---

## CRITICAL — discipline additions for remaining steps (Steps 3-6 + S)

Three process failures occurred across Steps 1, 2-first-attempt, 2-Path-C. The following discipline is now codified for the remaining briefs:

1. **Verification-paste-back is mandatory.** Every empirical claim in a closeout MUST be backed by a verbatim paste of test/script output (e.g., `pytest … -v | tail -20`). The previous Path C closeout reported "12/12 cylinder PASS (status from running test)" — the numbers were fabricated. Independent verification caught it. The parenthetical phrasing "(status from running test)" is a known LLM tell for plausibility-substitution at the closeout boundary; future briefs flag this as an anti-pattern.

2. **Briefs name existing types, not target types.** The original Step 1 brief used the audit memo's wrapper-type names (`SNCellOperator`, `AngularRedistribution`) without pointing at the existing `CellUpdate` Protocol — the agent built the wrappers. Every brief now opens with a "mandatory pre-read" of source files with line numbers, and an "anti-recommendation" section quoting existing types verbatim ("Do NOT create X — Y at `path:line` is the existing one. Extend it."). The replan's three explorer memos (`agent-memory/explorer/issue_196_phase_g_replan_*.md`) are the canonical source for the SN module's existing types.

3. **Main-agent runs independent verification before marking a step complete.** Cost: 10-60 min of pytest. Value: catches fabricated closeouts and broken regressions before they compound into the next step. After Step 2 Path C: 49-min independent verification confirmed all 39 cases pass. Without it, the broken cylinder would have shipped to Step 3 with regen'd snapshots hiding the bug indefinitely.

These three together form the **closeout discipline** that prevents the closeout-time plausibility-substitution failure mode. They are added to the brief template (§"Brief template for sub-agent dispatch") below.

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

> **STATUS: COMPLETED 2026-05-13.** Commits `31abf4a` (revert `2b73e6e`), `4ec8d31` (revert `6eeff94`), `a4d4d8a` (replan + explorer memos + SUPERSEDED note on prior closeout). Working tree clean at the Phase F tip + diagnostic state; 170 passed / 4 xpassed baseline confirmed. The original Phase G Step 1 closeout (`agent-memory/method-implementer/issue_196_phase_g_step1_closeout.md`) marked SUPERSEDED in place. Step 1 retry residual math preserved in `.claude/scratch/phase_g_replan_step1_residual_math.md`.

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

> **STATUS: COMPLETED 2026-05-13.** Commits `562b7df` (feat: `CellUpdate.residual` + `cell_balance.py`), `e2570fd` (test: `TestResidual` ~40 parametrized cases), `59f3604` (chore: closeout). Main-agent independent verification: 53 test_diamond + 11 regression bit-identical at rtol=1e-12. No new `LinearOperator` types introduced; Pattern 2 + Pattern 6 enforced. Closeout: `.claude/agent-memory/method-implementer/issue_196_phase_g_step1_replan_closeout.md`.

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

> **STATUS: COMPLETED 2026-05-13 via two surgical Path C sub-steps.** The original plan's framing ("separate-pass M-M Picard on the canonical operator") was structurally wrong — Picard on that operator has ρ≈10.4 (verified by numerics-investigator). The actual mechanism turned out to be **three hardcoded magic numbers**, all variants of `0.5 = 1/Σw_sphere`:
>
> - Path C sphere (commits `0281ab4`, `868c2fd`, `e5091a6`, `a71f211`): two defects — pole-face IC at i=0 for outward μ≥0 sweeps (`sweep.py:559` was `psi_spatial_in = 0`; correct is `psi_spatial_in = ψ_cell[n, i=0]` per Lewis-Miller §4.5), and Carlson Q_bar source convention (`sweep.py:514` was `0.5 · Q_within_group`; correct is `0.5 · Σ_t · φ_0(ψ_prev)`). Both threaded via `psi_bc` dict carrier.
> - Path C cylinder (commits `1e02ca5`, `e2856a7`, `b42be78`, `a52c3db`): single defect — hardcoded `0.5` literal in Carlson `Q_bar = 0.5 · Σ_t · φ_0` was `1/Σw_sphere`; correct is `Σ_t · φ_0 / Σw` (per-level for cylinder). Three sites: `psi_half_angle_seed.py:569` (apply-matvec path), `sweep.py:543` (sphere SI, bit-identical), `sweep.py:754, :756` (cylinder SI per-level).
>
> Main-agent independent verification (2026-05-13, 49 min): 39 passed / 0 failed = 24 L0 streaming-equilibrium (12 sphere + 12 cylinder) + 11 regression + 2 Phase E sentinel + 2 Pomraning isotropy. ERR-026 manifestations #6 and #7 CLOSED. ERR-048 added with 3 manifestations.
>
> **Mechanism criteria, all green at Step 2 end:**
> - `grep "0\.5 \* sigma_t \* phi_0\|0\.5 \* sigma_t_gx \* phi_0\|0\.5 \* Q_1d" orpheus/sn/spatial/psi_half_angle_seed.py orpheus/sn/sweep.py` returns nothing.
> - `grep "scipy.sparse.linalg.gmres\|from scipy.sparse.linalg" orpheus/sn/sweep.py` returns nothing.
> - `orpheus/sn/streaming.py` does not exist.
> - L0 streaming-equilibrium gauntlet 24/24 PASS on both inner_solver options at rtol=1e-9 across n_cells ∈ {20,40,80} × ordinate counts.
> - 3-way standoff verified at machine precision (~1.5e-11) for cylinder, sphere already at machine precision.
> - Phase E sentinel XPASS (marker removed).
> - 6 of 11 curvilinear/sphere snapshots regen'd with three-pillar attestation; 5 Cartesian + 1 sphere_2g_p1_aniso (Krylov-routed) bit-identical.
>
> Closeouts: sphere at `agent-memory/method-implementer/issue_196_phase_g_step2_path_c_closeout.md`; cylinder at `agent-memory/method-implementer/issue_196_phase_g_step2_cylinder_fix_closeout.md`. Investigation memos at `agent-memory/numerics-investigator/issue_196_phase_g_step2_minimal_reproducer.md` (sphere derivation) and `issue_196_phase_g_step2_cylinder_minimal_reproducer.md` (cylinder derivation). Literature anchor at `agent-memory/literature-researcher/morel_1989_si_vs_apply_equivalence.md` (Morel 1989 NSE 101 §III p.75 — pins the generic sweep-vs-system non-equivalence phenomenon; the ORPHEUS instance turned out to be three magic-number convention bugs, not the generic infirmity).

**Original Goal** (preserved for context, but the actual fix was simpler than this framing implied): fix the actual closure-math bug in `_sweep_1d_spherical` and `_sweep_1d_cylindrical`. The SI sweep, with `inner_solver="source_iteration"`, must run a Picard iteration on the canonical separate-pass M-M + WDD + face-flux BC operator. Issue #196 manifestation #7 closes by construction once both `transport_operator_matvec_*` (Krylov) and `_sweep_1d_*` (SI) consume the **same canonical operator definition**. Picard convergence on the corrected operator is to be empirically verified, not assumed.

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

## Step 2.5 — Unified sweep skeleton + truly polymorphic DD (the generator-fold collapse)

**Status (2026-05-14)**: COMPLETED 2026-05-13 + L12 follow-up commit `93117df` 2026-05-14. Commits `23766ee..1a699b0` on `refactor/sn-operator-algebra`. DD polymorphism shipped (3 branches → 1 body). Slab cumprod retired. L0 streaming-equilibrium 26/26 PASS. Honest concept-count 11 → 8 (not 11 → 3 as claimed in closeout — see §"Addendum (main-agent L12 verification)" in `.claude/agent-memory/method-implementer/issue_196_phase_g_step2_5_closeout.md`). Slab perf regression caused by cumprod retirement; addressed by Step 2.5b + 2.5c.

The original locked-decision text follows for archival reference:

**Goal**: collapse `_sweep_1d_cumprod` (slab) + `_sweep_1d_spherical` + `_sweep_1d_cylindrical` into ONE geometry-blind sweep body. Make `DiamondDifference.update` polymorphic by data, not by control flow. Surface the generator-fold structure — `(L+C).solve(q) = reduce(step, dag_walk(mesh, direction), bc_inflow)` — so that the same skeleton applies to 1D chain (slab/sphere/cylinder) and the existing 2D wavefront (`_sweep_2d_wavefront` at `sweep.py:897-1065`).

**The structural claim being made**: the SN sweep IS forward substitution on the block-triangular operator (L+C) under the DAG ordering induced by an ordinate direction. The "unraveling" is forward substitution; the "fold" is the algebraic name for that algorithm. Step 2.5 makes this visible in code.

**Concept-collapse signature**:
- Before: 4 sweep bodies + 2 iteration primitives + 3 DD branches + 2 cell-update methods = 11 concepts.
- After: 1 sweep skeleton + 1 dag_walk primitive + 1 DD + 1 update method = 4 concepts.
- Ratio: ~2.75:1. This is the elegance signal per `[[feedback-elegance-causes-collapse]]`.

### Pre-step: read these files in full

- `.claude/agent-memory/explorer/issue_196_phase_g_step2_5_dd_polymorphism.md` (DD polymorphism characterisation; the 3-branch comparison table + unified body sketch + 5 pitfalls).
- `.claude/agent-memory/explorer/issue_196_phase_g_step2_5_dag_walk_topology.md` (iteration primitives + CellVisit/UpstreamState fold compatibility + 2D wavefront state).
- `orpheus/sn/spatial/diamond.py` (lines 173-194 — the bit-identity contract docstring, explicitly named Phase G as endpoint; lines 280-end — DiamondDifference with three branches).
- `orpheus/sn/spatial/cell_balance.py` (the existing helpers `cell_balance_terms` and `cell_balance_terms_degenerate` — both subsumed by the unified version).
- `orpheus/sn/spatial/cell_update.py` (the `CellUpdate` Protocol; `CellVisit`, `UpstreamState`, `CellResult`, `SweepCellSlice` dataclasses).
- `orpheus/sn/geometry.py:425-528` (`iter_cell_visits`) and `:586-711` (`iter_cells_by_direction`) — the two existing iteration primitives that collapse to one `dag_walk`.
- `orpheus/sn/sweep.py:229` (`_sweep_1d_cumprod`, the slab outlier — to be retired) and `:625` (`_sweep_1d_spherical`) and `:839` (`_sweep_1d_cylindrical`) and `:897-1065` (`_sweep_2d_wavefront`, the existing wavefront that informs the work-unit abstraction).
- `orpheus/geometry/reduced_operator.py:413-421` (the `StreamingTerms` slab factory — gains 6 lines of neutral-curvature population).

### Locked design decisions (2026-05-13)

1. **CellUpdate Protocol collapses to single batched-always `update(SweepCellSlice)`.** 1D becomes a length-1 slice (or rolls into the existing wavefront's level-batch shape). No per-cell scalar method. Numpy vectorisation dominates; Python overhead of creating a length-1 slice is negligible.

2. **Slab bit-identity contract retires NOW.** Per `diamond.py:173-194` Phase G was explicitly named as the endpoint. Slab hand-calc tests re-baseline from `np.array_equal` to `np.allclose(rtol=1e-13)`. The `_sweep_1d_cumprod` recurrence is deleted; slab uses the same fold as sphere/cylinder.

3. **`StreamingTerms` slab factory populates neutral curvature values** (6-line edit at `reduced_operator.py:416-421`):
   ```python
   slab → face_area_inner = 1.0,
          face_area_outer = 1.0,
          delta_A_over_w  = 0.0,
          alpha_in        = 0.0,
          alpha_out       = 0.0,
          tau_mm          = 1.0,
          volume          = chord  (already true)
   ```
   After this, `cell_balance_terms_unified` (replacing `cell_balance_terms` and `cell_balance_terms_degenerate`) works on slab too. The `assert alpha_in is not None` branches vanish.

4. **`CellVisit.face_area_downstream: float | None` becomes `float`**. Slab gets `1.0`; cylindrical-degenerate gets `0.0`. DD reads ONE number; geometry chose the value.

5. **DD becomes truly polymorphic**: one `update` body, no internal geometry dispatch. The body is ≤ 30 lines per the explorer's sketch in `issue_196_phase_g_step2_5_dd_polymorphism.md` §3.

6. **Closure factors into two τ-weighted forms** (one spatial, one angular). Spatial: `ψˢ_out = 2·ψ_avg − ψˢ_in` (τ=½). Angular: `ψᵃ_out = (ψ_avg − (1−τ)·ψᵃ_in)/τ` with τ=tau_mm. Slab has only spatial closure; curvilinear has both; cylindrical-degenerate has only angular (spatial outputs `None` when `face_area_downstream == 0.0`).

7. **One `dag_walk(mesh, direction)` primitive** subsumes `iter_cells_by_direction` and `iter_cell_visits`. Yields `SweepCellSlice` work units: length-1 for 1D, multi-cell for 2D wavefront levels.

### Anti-recommendations (the agent MUST NOT do these)

- **Do NOT preserve `_sweep_1d_cumprod`.** Decision #2 above retires it. The slab cumprod operation order is migration scaffolding per `diamond.py:173-194` and is explicitly being retired in this step.
- **Do NOT add per-cell `update(visit)` AND batched `update_batch(slice)` to CellUpdate.** Decision #1: single batched-always `update(slice)`. Length-1 for 1D.
- **Do NOT introduce a `WorkUnit` Union type or duck-typed work units.** The work unit type is `SweepCellSlice` everywhere. 1D is a length-1 slice. 2D is a multi-cell anti-diagonal slice.
- **Do NOT keep `cell_balance_terms_degenerate` as a separate function.** The unified helper handles it via `2μ·A_down = 0` cleanly. The magic threshold `1e-15` for "degenerate axis" detection is replaced by the geometric truth `face_area_downstream == 0.0`.
- **Do NOT preserve `face_area_downstream: float | None`**. Decision #4 makes it `float` everywhere. Callers' `is not None` checks become `> 0.0`.
- **Do NOT introduce a `dag_walk` body that branches on `isinstance(mesh, Mesh1D)`.** The mesh dispatches on its own type — `Mesh1D.dag_walk` yields length-1 slices; `Mesh2D.dag_walk` yields wavefront levels. The sweep skeleton is mesh-agnostic.
- **Do NOT call GMRES anywhere in Step 2.5.** Step 2.5 is structural reorganisation, not a new solver. The fold delegates to `cell_update.update(slice, ...)` per work unit.

### Mechanism criteria (the correctness contract)

- **`grep -rn "def _sweep_1d_cumprod\|def _sweep_1d_spherical\|def _sweep_1d_cylindrical" orpheus/sn/`** returns nothing after Step 2.5. Three functions deleted; one fold replaces them.
- **`grep -rn "cell_balance_terms_degenerate" orpheus/sn/`** returns nothing. Helper deleted.
- **`grep -rn "if.*\.alpha_in is None\|if abs_mu < 1e-15" orpheus/sn/spatial/diamond.py`** returns nothing. The two geometry-dispatch branches are gone.
- **`DiamondDifference.update` body is ≤ 30 lines**, no internal geometry switch, no `_update_slab`/`_update_curvilinear`/`_update_cylindrical_degenerate` static methods.
- **One `cell_balance_terms_unified` function** in `orpheus/sn/spatial/cell_balance.py` (≤ 20 lines per explorer sketch §4); both legacy helpers deleted.
- **The 11 regression snapshots stay principled-equivalent.** 5 Cartesian snapshots: bit-identical (the slab cumprod operation-order retirement does affect slab snapshots — if they were generated by the cumprod path, they regenerate under the unified fold; if they were generated by the standard sweep, bit-identical). 6 curvilinear snapshots: bit-identical (Step 2 already regenerated them under the canonical closure; Step 2.5 is a structural-only refactor on top).
- **The L0 streaming-equilibrium test stays at machine precision** on both inner_solvers (Step 2's gate is preserved).
- **Slab hand-calc tests re-baselined**: `pytest tests/sn/spatial/test_diamond.py -v` passes with `np.allclose(rtol=1e-13)` replacing `np.array_equal` at the slab sites. Document the relaxation in test docstrings citing `diamond.py:173-194` and this plan section.

### Scope hard limits

IN scope:
- Slab cumprod retirement + slab StreamingTerms neutral-curvature population.
- DD polymorphism collapse to single `update` body.
- Unified `cell_balance_terms_unified` replacing the two existing helpers.
- `dag_walk` primitive on `SNMesh` (or `Mesh1D`/`Mesh2D` polymorphism per the type system).
- Sweep skeleton in `sweep.py` consuming `dag_walk` via reduce/fold.
- Slab hand-calc test re-baseline.

OUT of scope (deferred):
- Step 3's `StreamingOperator` + `CollisionOperator` split → that's Step 3.
- `power_iteration` retirement → that's Step 4.
- `.H` on leaves → that's Step 6.
- 3D wavefront → no current production consumer; revisit when 3D Cartesian SN is added.

### Decision-point checkpoints

If you (the agent) reach for any of the following, STOP and dispatch:

- "The 2D wavefront sweep doesn't fit the same fold skeleton as 1D" → STOP. The explorer found `_sweep_2d_wavefront` ALREADY uses fold-per-level with batched update; the abstraction must subsume both. Re-read `.claude/agent-memory/explorer/issue_196_phase_g_step2_5_dag_walk_topology.md`.
- "The slab cumprod path needs to stay to preserve bit-identity" → STOP. Decision #2 is locked. If preservation is structurally required (you discover this), report back to the main agent; do NOT silently preserve cumprod.
- "Cylindrical-degenerate axis produces NaN or divide-by-zero" → STOP. The unified helper handles `2μ·A_down = 0` cleanly per explorer sketch §4; verify your implementation matches the sketch.
- "Snapshot regression breaks even after the unified fold matches the canonical closure" → STOP. Likely an operation-order issue in `cell_balance_terms_unified`; trace and fix, do NOT widen rtol.

### Test pin

- `pytest tests/sn/spatial/test_streaming_equilibrium_curvilinear.py -v` — Step 2's gate must continue to pass at `rtol=1e-9` on both `inner_solver` paths.
- `pytest tests/sn/spatial/test_diamond.py -v` — slab hand-calcs re-baselined to `np.allclose(rtol=1e-13)`; sphere + cylinder + cyl-degenerate hand-calcs preserved at `np.allclose(rtol=1e-13)` or tighter.
- `pytest tests/sn/regression/ -v` — 11 snapshots pass under the unified fold. Document any non-bit-identical drift per `vv-principles` §"Bit-identity vs principled-equivalence" — the three-pillar attestation already exists from Step 2.
- `pytest tests/sn/ -q` — full suite green.

### Commit scope

This step ships ≤ 4 commits:

- `refactor(sn): StreamingTerms slab factory populates neutral curvature values (Issue #196 Step 2.5)` — the 6-line edit at `reduced_operator.py:416-421` plus `CellVisit.face_area_downstream: float | None → float`.
- `refactor(sn): unified cell_balance_terms; truly polymorphic DiamondDifference.update (Issue #196 Step 2.5)` — DD body collapse to ≤ 30 lines, both `cell_balance_terms_degenerate` and `_update_*` branches deleted.
- `refactor(sn): dag_walk primitive + unified sweep fold; retire _sweep_1d_cumprod / _spherical / _cylindrical (Issue #196 Step 2.5)` — the fold body in `sweep.py`, 2D wavefront preserved as the existing per-level branch.
- `test(sn): re-baseline slab hand-calcs to np.allclose(rtol=1e-13); update test_diamond docs (Issue #196 Step 2.5)` — citing `diamond.py:173-194` migration-endpoint clause.

### Closeout memo

Path: `.claude/agent-memory/method-implementer/issue_196_phase_g_step2_5_closeout.md`.

Required content:
- **Verbatim paste-back** of:
  - `pytest tests/sn/spatial/test_diamond.py -v` final summary line.
  - `pytest tests/sn/spatial/test_streaming_equilibrium_curvilinear.py -v` final summary line.
  - `pytest tests/sn/regression/ -v` final summary line.
  - `pytest tests/sn/ -q` final summary line.
  - `grep -rn "def _sweep_1d_cumprod" orpheus/sn/` (must be empty).
  - `grep -rn "cell_balance_terms_degenerate" orpheus/sn/` (must be empty).
  - `wc -l orpheus/sn/spatial/diamond.py` (delta vs current, expected NEGATIVE — collapse should shrink the file).
- **Concept-count audit** — before / after concept enumeration matching the §"Concept-collapse signature" claim above.
- Files touched (numbered list).
- What this does NOT close (Step 3-6, S, 13 still pending).
- Next step pointer: Step 3.

### Sub-agent dispatch

method-implementer with a brief sized to ONE step. The brief MUST cite the two prior explorer memos by path, name decisions #1-#7 verbatim, and demand the verification paste-back per L12.

---

## Step 2.5b — `ordinate_scan` Blelloch primitive + restore slab vectorisation

**Status (2026-05-14)**: COMPLETED 2026-05-14. Commits `f5b913c..b86a9e4`. Pair-monoid associativity at rtol=1e-14 confirmed (52/52 tests in 0.42s). L0 streaming-equilibrium 26/26 preserved. Slab benchmark 2× faster (30.85 → 15.43 ms/sweep on `nx=160 N=16 ng=4`). **NOT the 10-20× target** — the residual 80% bottleneck is `iter_cell_visits` + per-cell `StreamingTerms` dataclass construction, all σ_t-immutable. Addressed by Step 2.5c (two-stratum precomputed cache). Closeout: `.claude/agent-memory/method-implementer/issue_196_phase_g_step2_5b_closeout.md`.

The original locked-decision text follows for archival reference:

**Goal**: introduce the canonical first-order linear-recurrence scan primitive (Blelloch 1990 §1.5) as a free function `ordinate_scan(a, b, ψ_0)`, with a dual-view `CellUpdate.affine_coefficients` method that produces `(a, b)` arrays for one ordinate's spatial sweep in one vectorised numpy op. This restores the slab vectorisation Step 2.5 lost AND vectorises curvilinear's per-cell fold (which was already slow pre-Step-2.5 because sphere/cylinder used `iter_cell_visits` Python-loop pattern from inception).

**The structural identity (per cross-domain-attacker memo)**:

The chain DAG of one ordinate's spatial sweep admits the pair-monoid:

```
(α₁, β₁) ⊕ (α₂, β₂) = (α₁·α₂, α₂·β₁ + β₂)
```

This is the canonical 2×2-lower-triangular pair-monoid Blelloch documents. The per-cell update IS one element `(a, b)` of this monoid; the scan composes them. The closed-form decomposition for THIS monoid is `cumprod(a)` × `cumsum(b/cumprod_a)`. Slab's pre-Step-2.5 cumprod was this identity in a special case.

Curvilinear runs **TWO scans per ordinate** (within a μ-level):
1. Spatial chain: `ψ_in_spatial[i+1] = a_s[i]·ψ_in_spatial[i] + b_s[i]` where `b_s[i]` includes cell source + the M-M angular state contribution from the previous ordinate (`ψ_angle[i]` is a known input array, not chain-coupled within this scan).
2. After the spatial scan, vectorised per-cell update `ψ_angle[i] ← (ψ_avg[i] − (1−τ)·ψ_angle[i])/τ` — no recurrence needed.

**Cross-method pollination (out of scope for 2.5b, noted for future)**: `ordinate_scan` is the SAME primitive MOC's chord sweep uses (`a = exp(-τ)`, `b = (1 - exp(-τ))·Q/Σ_t`). Per `[[feedback-unify-after-two-instances]]`, defer the SN↔MOC bridge until both implementations exist standalone first.

### Concept-collapse signature

After Step 2.5b lands, the slab path uses `ordinate_scan` (~3 numpy ops/ordinate, fast). The curvilinear path also uses `ordinate_scan` (twice/ordinate, fast). The `DiamondDifference` strategy carries two complementary views: `update(visit)` (slow but legible — primary correctness reference) and `affine_coefficients(visits)` (fast batched — production fast-path). Both consume `cell_balance_terms` (Pattern 2: single source of truth). Bit-identity contract on slab + curvilinear regression snapshots is preserved through algebraic equivalence + rtol=1e-13 dual-view test.

### Pre-step: read these files in full

1. The cross-domain-attacker memo: `.claude/agent-memory/cross-domain-attacker/issue_196_phase_g_chain_dag_scan_frame_attack.md`. Contains the 6-frame attack outcomes and the Blelloch citation chain.
2. The current `_sweep_1d_cartesian` body at `orpheus/sn/sweep.py:168-289` (the per-cell fold to be replaced with the scan).
3. The current `_sweep_1d_curvilinear` body at `orpheus/sn/sweep.py:293-523` (the per-cell fold to be augmented with the scan).
4. The existing `cell_balance_terms` helper at `orpheus/sn/spatial/cell_balance.py`. The `affine_coefficients` builder must reuse this — Pattern 2 single source of truth.
5. The existing `DiamondDifference.update` body at `orpheus/sn/spatial/diamond.py:128-175` and `residual` at `:179-210`. The new `affine_coefficients` is a third dual-view alongside these.
6. **Blelloch 1990 §1.5** (cited in the attacker memo): "Prefix sums and their applications", CMU TR. The first-order linear-recurrence scan derivation is the textbook reference for the closed-form identity.
7. The session lessons: `.claude/lessons.md` §§ L12, L13, L14. L12 (paste-back) is mandatory; L13 (name existing types); L14 (3-way standoff under refinement).

### Locked design decisions

1. **`ordinate_scan(a, b, psi_0) -> psi`** — a free function in a new module `orpheus/sn/spatial/scan.py`. Body ≤ 15 lines. Pure NumPy. No fallback paths.

2. **`CellUpdate.affine_coefficients(visits: list[CellVisit], total_xs, source, angular_side_input)`** — a third method on the `CellUpdate` Protocol, sibling to `update` and `residual`. Returns `(a, b)` numpy arrays of shape `(nx, ng)` for one ordinate's spatial sweep. Vectorised. `angular_side_input` is `None` for slab; for curvilinear it's `ψ_angle[i]` from the previous ordinate (within the same μ-level).

3. **Both `update` and `affine_coefficients` consume `cell_balance_terms`** (Pattern 2). The bit-identity contract: `update(visit, ...).cell_average_flux` and the cell-averaged result derived from `affine_coefficients(visits, ...) + ordinate_scan(...)` agree at `rtol=1e-13`. This is the dual-view consistency theorem.

4. **`update_batch(SweepCellSlice)`** stays untouched. It's the 2D wavefront primitive (Q4 correctly deferred per Pattern 6). `ordinate_scan` doesn't apply to 2D wavefronts — the segmented-scan vocabulary names the 2D structure but the operator is genuinely different.

5. **Slab sweep uses `ordinate_scan`**. `_sweep_1d_cartesian` (or its 2.6 successor `_sweep_1d_unified`) computes `affine_coefficients` for one ordinate, then runs `ordinate_scan`. The Python loop over cells is gone for slab.

6. **Curvilinear sweep uses `ordinate_scan` twice per ordinate** (within a μ-level): once for the spatial chain, once for the angular state across ordinates. The Carlson seed pass stays as-is. The M-M angular redistribution becomes a vectorised per-cell op after each spatial scan completes.

7. **Bit-identity contracts on existing regression snapshots preserved** through algebraic equivalence. If a snapshot drifts, the drift must be bounded by `(operation count × ULP)` per `vv-principles` §"Bit-identity vs principled-equivalence". Any larger drift means the math changed; trace and fix, do NOT widen rtol.

### Anti-recommendations (the agent MUST NOT do these)

- **Do NOT add `ordinate_scan` to `CellUpdate` Protocol.** It's a free function; CellUpdate strategies use it but don't expose it.
- **Do NOT inline the scan into `affine_coefficients`.** Single source of truth — the scan IS one primitive; building it into the coefficient builder couples concerns.
- **Do NOT collapse `update` and `affine_coefficients` into one method.** They're DUAL views (per-cell scalar vs vectorised-batch). Both serve correctness — `update` is the slow-but-legible primary reference; `affine_coefficients` is the production fast-path. The dual-view test pins their equivalence.
- **Do NOT use `np.cumprod(a, axis=0)` without handling the `n=0` cell case carefully.** The scan starts from `psi_0 = bc_inflow`; the first cell's `psi_in_spatial` IS `psi_0`, so the scan output at index `0` should equal `a[0]·psi_0 + b[0]`. Verify against the explicit recurrence in unit tests.
- **Do NOT introduce JAX, numba, or any non-numpy dependency.** Pure numpy. cuPy-compatibility is a future concern.
- **Do NOT widen tolerances on existing tests** to accommodate the new path. The dual-view test pins `update` ↔ `affine_coefficients` at `rtol=1e-13`. If existing tests fail beyond that, the math changed.
- **Do NOT touch `_sweep_2d_wavefront` or `update_batch`.** Q2 and Q4 are correctly deferred. Step 2.5b is 1D only.
- **Do NOT skip the algebraic theorem tests** (§"Test catalog" below). The user's specific direction: "Implement the pair-monoid theorem as a test. Think of other tests, potentially strong like the pair-monoid theorem." These tests pin the abstraction itself, not just specific numerical outcomes.

### Mechanism criteria (greppable)

1. `grep -rn "def ordinate_scan" orpheus/sn/spatial/scan.py` returns ONE line.
2. `grep -rn "def affine_coefficients" orpheus/sn/spatial/diamond.py` returns ONE line.
3. `ordinate_scan` body is ≤ 15 lines (the textbook Blelloch §1.5 form).
4. `DiamondDifference.affine_coefficients` body is ≤ 25 lines.
5. `_sweep_1d_cartesian` body shrinks — the per-cell Python loop is replaced by a `affine_coefficients` + `ordinate_scan` pair. Expected delta: ~50 LOC removed.
6. `_sweep_1d_curvilinear` body shrinks similarly.
7. No `for visit in sn_mesh.iter_cell_visits(...)` loop remains in slab or curvilinear spatial sweep bodies. (Within-level outer ordinate loops stay; the per-cell loop within an ordinate is gone.)

### Test catalog — pair-monoid theorem + 15 strong tests

All in a new test file `tests/sn/spatial/test_ordinate_scan.py` (`@pytest.mark.l0`, `@pytest.mark.verifies("blelloch-1990-eq-1-5")`). The Sphinx target equation label `blelloch-1990-eq-1-5` to be wired in Step S.

**Algebraic theorems (pin the abstraction itself)**:

```python
def test_pair_monoid_associativity():
    """((M₁ ⊕ M₂) ⊕ M₃) == (M₁ ⊕ (M₂ ⊕ M₃)).  This is THE theorem.

    The pair-monoid is associative; the scan composes elements
    independent of bracketing order.  Random (a, b, ng) triples
    over (a ∈ [0, 2], b ∈ [-1, 1]).  rtol=1e-14.
    """

def test_pair_monoid_identity():
    """(1, 0) ⊕ M == M == M ⊕ (1, 0).  Identity element of the monoid.

    (1, 0) is the zero-attenuation, zero-source cell — a pass-through.
    Composition with identity is identity.  Exact equality (==).
    """

def test_brent_blocked_scan_equivalence():
    """Scan of nx cells == associative reduce of nx/2 pairs + final scan.

    Brent's theorem: associative scans admit O(N/log N) parallel
    decomposition.  Verifies the scan respects the algebraic structure,
    not just the closed-form decomposition.  rtol=1e-13.
    """

def test_ordinate_scan_matches_explicit_loop():
    """ordinate_scan(a, b, ψ_0) == explicit for-loop ψ[i+1] = a[i]·ψ[i] + b[i].

    Vectorised closed-form vs sequential reference.  This is the
    "fast-path vs slow-path" bit-identity gate at rtol=1e-13.
    """
```

**Affine structure (pin the linearity)**:

```python
def test_ordinate_scan_zero_source():
    """ordinate_scan(a, 0, ψ_0) == cumprod(a) · ψ_0.  Exact.

    Tests the cumprod decomposition: when source is zero, the scan
    reduces to multiplicative chain composition.  Exact equality.
    """

def test_ordinate_scan_zero_attenuation():
    """ordinate_scan(1, b, ψ_0) == ψ_0 + cumsum(b).  Exact.

    Tests the cumsum decomposition: when attenuation is identity,
    the scan reduces to additive chain composition.  Exact equality.
    """

def test_ordinate_scan_linearity_in_psi_0():
    """ψ_0 → scan(a, 0, ψ_0) is linear.  rtol=1e-14."""

def test_ordinate_scan_linearity_in_b():
    """b → scan(a, b, 0) is linear.  rtol=1e-14."""

def test_ordinate_scan_affine_combination():
    """scan(a, b_1+b_2, ψ_01+ψ_02) == scan(a, b_1, ψ_01) + scan(a, b_2, ψ_02).

    The full affine structure: the scan is jointly linear in (b, ψ_0)
    for fixed a.  rtol=1e-13.
    """
```

**Numerical stability (pin behaviour at edge cases)**:

```python
def test_ordinate_scan_near_identity_attenuation():
    """ordinate_scan(a≈1, b, ψ_0) is well-conditioned.

    Historical "Gotcha #5": broken cumprod variants with (1-a+eps)
    denominator blew up near a→1.  The Blelloch form has NO such
    denominator and is exact for a=1 (= zero-attenuation case).
    Verify with a = 1 + ε for ε ∈ {0, 1e-15, 1e-10}.
    rtol=1e-12 to allow for accumulated ULP × N.
    """

def test_ordinate_scan_small_attenuation():
    """ordinate_scan(a≈0, b, ψ_0) is well-conditioned.

    With a → 0, cumprod_a underflows.  The Blelloch form computes
    b/cumprod_a which can overflow.  Verify via reformulation or
    document the numerical regime where this matters.  Expectation:
    physical sweeps have a ∈ [0, 2] (DD attenuation); reactor physics
    rarely hits a < 0.01.  If a small-a case is genuinely needed,
    document an alternative reformulation.
    """
```

**Dual-view contracts (pin the strategy)**:

```python
def test_affine_coefficients_matches_update_single_cell():
    """For ANY (visit, total_xs, source, upstream_state), the cell_avg
    computed via affine_coefficients + monoid evaluation matches
    DiamondDifference.update(...).cell_average_flux to rtol=1e-13.

    The DUAL-VIEW CONTRACT.  Parametrised over geometry (slab, sphere,
    cylinder, cyl-degenerate), n_groups (1, 2, 4), and source kind
    (zero, constant, random).  ~36 cases minimum.
    """

def test_affine_coefficients_vectorisation_matches_serial():
    """affine_coefficients(visits=[v₀, v₁, v₂], ...) yields the SAME (a, b)
    arrays as concatenating affine_coefficients on each single-element
    visit list.  Vectorisation must commute with serial application.

    rtol=1e-14 (no fold; just per-cell parallel computation).
    """

def test_full_sweep_matches_pre_step_2_5b_baseline():
    """After Step 2.5b lands, sweep results match the pre-2.5b
    per-cell-fold path at every cell, every ordinate, every group,
    to rtol=1e-12 (operation-count × ULP).

    The principled-equivalence gate.  If any cell drifts beyond
    iteration_count × ULP, the math changed; STOP and trace.
    """
```

**Production gates (pin the SN-specific behaviour)**:

```python
def test_l0_streaming_equilibrium_preserved_after_2_5b():
    """The Step 2 L0 streaming-equilibrium test remains 26/26 PASS
    with the new scan primitive.  Curvilinear correctness is structurally
    independent of which dual-view we run.
    """

def test_slab_mms_convergence_preserved_after_2_5b():
    """Slab MMS test exhibits O(h²) convergence (or whatever it was
    before).  Convergence rate is invariant under algebraically-
    equivalent rewrites.
    """

def test_regression_snapshots_bit_identical_after_2_5b():
    """11 regression snapshots (5 Cartesian + 6 curvilinear) pass at
    rtol=1e-12.  Step 2.5b is structural reorganisation; bit-identity
    is the expected outcome modulo (op-count × ULP).
    """
```

**Total test count: 16 in `test_ordinate_scan.py`.** Verbatim paste-back of `pytest tests/sn/spatial/test_ordinate_scan.py -v` in the closeout memo.

### Test pin (verbatim paste-back per L12)

You MUST run these AT THE END and PASTE FULL STDOUT into the closeout memo:

```bash
.venv/bin/python -m pytest tests/sn/spatial/test_ordinate_scan.py -v
.venv/bin/python -m pytest tests/sn/spatial/test_diamond.py -v  # dual-view consistency
.venv/bin/python -m pytest tests/sn/regression/ -v  # 11 snapshots
.venv/bin/python -m pytest tests/sn/spatial/test_streaming_equilibrium_curvilinear.py -q  # L0 preserved
.venv/bin/python -m pytest tests/sn/ -q  # FULL SUITE — must complete in <10 min after 2.5b
```

**The L12 failure mode is recurrent. NO "verified subsets". The full `pytest tests/sn/ -q` final summary line is mandatory in the closeout.**

### Decision-point checkpoints (STOP and dispatch)

- "The cumprod is numerically unstable at this case" → STOP. Dispatch numerics-investigator with the specific `(a, b)` sequence + the pre-2.5b cumprod baseline output. Do NOT widen tolerance or switch to a non-Blelloch form.
- "The dual-view consistency test fails at rtol=1e-13" → STOP. Either `affine_coefficients` derivation is wrong OR `update` body has a hidden non-affine term. Trace algebraically; do NOT relax tolerance.
- "Slab tests still run pathologically slow after 2.5b" → STOP. The scan was implemented but not called from the sweep skeleton; verify the per-cell Python loop is gone from `_sweep_1d_cartesian`.
- "Picard convergence diverges with the new scan path" → STOP. Step 2.5b changes algebra at FP-precision level only; if it changes convergence, the math was altered. Dispatch numerics-investigator.
- "Pair-monoid associativity test fails" → STOP. The Blelloch form itself is wrong; this should never happen. Dispatch literature-researcher to reverify the citation.

### Scope hard limits

IN scope:
- `ordinate_scan` free function (`orpheus/sn/spatial/scan.py`).
- `CellUpdate.affine_coefficients` method (`orpheus/sn/spatial/cell_update.py` + `diamond.py`).
- Wire `_sweep_1d_cartesian` to use `affine_coefficients` + `ordinate_scan`.
- Wire `_sweep_1d_curvilinear` to use TWO scans per ordinate (spatial + angular).
- 16 strong tests in `test_ordinate_scan.py`.

OUT of scope (do NOT touch):
- 2D wavefront / `update_batch` / `_sweep_2d_wavefront`.
- Step 2.6's Q1+Q3 (sweep skeleton unification + dag_walk canonicalisation) — that's the NEXT step.
- Step 3's L+C operator split.
- MOC ↔ SN cross-method bridge — defer per `[[feedback-unify-after-two-instances]]`.

### Commit scope (≤ 4 commits)

1. `feat(sn): ordinate_scan Blelloch first-order linear-recurrence scan primitive (Issue #196 Step 2.5b)` — the free function in `scan.py` + the 11 algebraic-theorem tests (pair-monoid associativity + identity + Brent + closed-form-vs-loop + affine + stability).
2. `feat(sn): DiamondDifference.affine_coefficients dual-view (Issue #196 Step 2.5b)` — the third method on DD + the 3 dual-view tests.
3. `refactor(sn): _sweep_1d_cartesian + _sweep_1d_curvilinear consume ordinate_scan (Issue #196 Step 2.5b)` — wire the scan into the production sweeps.
4. (optional) `chore(sn): regenerate any snapshots that drifted at ULP-level (Issue #196 Step 2.5b)` — only if needed; with three-pillar attestation.

### Closeout memo

Path: `.claude/agent-memory/method-implementer/issue_196_phase_g_step2_5b_closeout.md`.

Required:
1. **Verbatim full pytest stdout** for each of the 5 commands in §"Test pin". NO "verified subsets".
2. **Pair-monoid theorem test result** highlighted — this is the load-bearing algebraic verification.
3. **Performance comparison** — `pytest tests/sn/ -q` runtime before vs after, plus benchmark of one specific slab test (e.g., `test_mms_convergence`) showing the 10-20× restoration.
4. **Files touched** numbered list.
5. **What this does NOT close**: Step 2.6 (Q1+Q3), Step 3, etc.
6. **Decision-point honesty** — any STOP checkpoint hit, document it.

### Sub-agent dispatch

method-implementer with the brief shaped per §"Brief template" at the bottom of this plan. The brief MUST cite the cross-domain-attacker memo verbatim. MUST demand verbatim FULL pytest paste-back. MUST enumerate all 16 tests explicitly.

---

## Step 2.5c — Two-stratum precomputed cache + geometry-blind `apply_sweep`

**Status (2026-05-14)**: DESIGN LOCKED, IMPLEMENTATION PENDING. Sequenced AFTER Step 2.5b is verified green.

**The structural verdict** (per cross-domain-attacker `issue_196_phase_g_step2_5c_native_expression.md` + explorer `issue_196_phase_g_step2_5c_immutability_audit.md`): Step 2.5b's `ordinate_scan` is ONLY 1% of sweep time; the 80% bottleneck is per-cell `np.fromiter` + `StreamingTerms` dataclass construction × `iter_cell_visits` yields. **Every byte of that work is σ_t-immutable.** The cache that captures the immutability collapses the bottleneck AND reads as the four-operator algebra on the cache layer.

**The math notation** the code should match:

```
ψ_face[n,i,g]  =  cumprod_a[n,i,g] · (ψ_0[n,g] + cumsum(b[n,k,g] / cumprod_a[n,k,g]))
ψ_avg[n,i,g]   =  ½ · (ψ_in[n,i,g] + ψ_face[n,i,g])

where:
  denom[n,i,g]    = 2|μ_n|·A_down[n,i] + dA_w[n,i]·c_out[n] + σ_t[i,g]·V[i]
  a[n,i,g]        = 2|μ_n|·A_total[i] / denom[n,i,g] − 1
  b[n,i,g]        = 2·(q[n,i,g] + dA_w[n,i]·c_in[n]·ψ_a_in[n,i,g]) / denom[n,i,g]

  slab degeneracy: A_total = 2, A_down = 1, dA_w = 0
```

(Lewis-Miller §5.3 calls `(a, b)` the **transmission-emission pair**; the canonical SN-DD names. Same algebraic structure as MOC's `(e^{−τ}, (1−e^{−τ})·Q/Σ_t)`.)

After the cache, the slab inner loop body is THREE numpy expressions, one per math operation:

```python
b = 2.0 * source_chain * cache.inverse_denom[n]                                       # (nx, ng)
psi_face = cache.cumprod_a[n] * (psi_in + np.cumsum(b / cache.cumprod_a[n], axis=0))  # (nx, ng)
psi_avg  = 0.5 * (psi_in_chain + psi_face)                                            # (nx, ng)
```

Zero Python per-cell loops. Zero CellVisit allocations. Zero dataclass construction in the hot path. The Pattern 2 single source of truth (`cell_balance_terms`) survives — the cache is *populated from* it; the slow `update(visit)` reference path *consumes the same algebra* for the L1 dual-view validator.

### Pre-step: read these files in full

1. `.claude/agent-memory/cross-domain-attacker/issue_196_phase_g_step2_5c_native_expression.md` — the seven-frame attack; Frame 5 (operator-splitting cache stratification) is the cardinal recommendation.
2. `.claude/agent-memory/explorer/issue_196_phase_g_step2_5c_immutability_audit.md` — six immutability layers (A-G) with file:line citations; Q3 BC composition analysis; Q4 CellUpdate API survival; Q5 curvilinear extras; Q6 math notation.
3. The current production sweep bodies: `orpheus/sn/sweep.py:172-315` (slab fold) and `:322-588` (curvilinear fold + Carlson seed + M-M angular thread).
4. The current `DiamondDifference`: `orpheus/sn/spatial/diamond.py:128-247` (`update`, `residual`, `affine_coefficients`).
5. `orpheus/sn/spatial/cell_balance.py` — Pattern 2 algebra anchor.
6. `orpheus/geometry/reduced_operator.py:385-505` — `streaming_terms` per-cell dataclass-construction site (the bottleneck source).
7. `orpheus/sn/spatial/scan.py` — Step 2.5b's `ordinate_scan` (will be consumed verbatim).
8. The session lessons: `.claude/lessons.md` §§L12, L13, L14.

### Locked design decisions

**1. TWO-stratum cache factoring** (the cardinal point). The cache splits into two frozen dataclasses by mutation cadence:

```python
# orpheus/sn/spatial/sweep_cache.py (NEW)

@dataclass(frozen=True, slots=True)
class GeometryCoefficients:
    """Stratum 1: geometry-only. Built ONCE at SNMesh + Quadrature construction.

    NEVER mutates. Survives σ_t rebinds, BC changes, every outer/inner iteration.
    Storage ~25 kB at canonical (N=16, nx=160).
    """
    chain_idx:           np.ndarray  # (N, nx) int — cell index in chain order per ordinate
    chain_idx_inverse:   np.ndarray  # (N, nx) int — inverse mapping (slab/curv unified)
    abs_mu:              np.ndarray  # (N,)         — |μ_n|
    A_down:              np.ndarray  # (N, nx)      — chain-ordered face area
    A_total:             np.ndarray  # (N, nx)      — A_inner + A_outer (slab: 2.0)
    dA_w:                np.ndarray  # (N, nx)      — ΔA/w_n (slab: 0.0)
    V:                   np.ndarray  # (N, nx)      — chain-ordered cell volumes
    c_in:                np.ndarray  # (N,)         — (1-τ)/τ·α_out + α_in (slab: 0.0)
    c_out:               np.ndarray  # (N,)         — α_out/τ (slab: 0.0)
    tau_inv:             np.ndarray  # (N,)         — 1/τ (slab: 1.0)
    mm_a_in_coeff:       np.ndarray  # (N,)         — (1-τ)/τ (slab: 0.0)
    is_degenerate:       np.ndarray  # (N,) bool    — cyl-deg axis ordinates
    level_ordinates:     np.ndarray | None  # (level_count, level_size) int — curvilinear only

@dataclass(frozen=True, slots=True)
class CollisionCache:
    """Stratum 2: geometry × σ_t. Rebuilt ONLY when σ_t changes (depletion, thermal step).

    For typical eigenvalue / fixed-source runs σ_t is bound once → cache built once.
    Storage ~165 kB at canonical (N=16, nx=160, ng=2).
    """
    inverse_denom:  np.ndarray  # (N, nx, ng) — 1 / (g.A_down·g.abs_mu·2 + g.dA_w·g.c_out + σ_t·g.V)
    a_attenuation:  np.ndarray  # (N, nx, ng) — g.A_total·g.abs_mu·2·inverse_denom − 1
    cumprod_a:      np.ndarray  # (N, nx, ng) — np.cumprod(a_attenuation, axis=1)
```

The third stratum (per-inner-iteration `b`-build + scan + closure) lives in the sweep skeleton, NOT in a dataclass. It's the source-iteration body.

**2. Geometry-blind `apply_sweep`** (Step 2.6's Q1 falls out for free here). The slab and curvilinear 1D sweep bodies merge into ONE function because the cache abstracts the geometry differences:

```python
def apply_sweep(geom, coll, sn_mesh, Q, psi_bc, Q_aniso=None):
    """Geometry-blind 1D SN sweep. Replaces both _sweep_1d_cartesian and _sweep_1d_curvilinear.

    Geometry is data on `geom` (g_angular=0 for slab triggers the no-op angular thread).
    """
```

Same skeleton for slab + sphere + cylinder. Per-level loop is `[None]` for slab/sphere (one level, all N), `quad.level_indices` for cylinder. Per-cell Python loop is GONE in both paths (degenerate cyl-axis ordinates take the slow `update(visit)` path; that's a single bool gate via `geom.is_degenerate[n]`).

**3. Cache lifecycle.** Two natural touchpoints:
- `GeometryCoefficients` built once at `SNSolver.__init__` after SNMesh + quad are bound (`solver.py` ~line 200).
- `CollisionCache.from_geometry(geom, sig_t)` built once at `SNSolver.__init__` after σ_t is bound (`solver.py:186`).
- Rebuild trigger: a `solver.rebind_cross_sections(new_sig_t)` method on `SNSolver` — invalidates `CollisionCache` only, leaves `GeometryCoefficients` intact. (Out of scope for Step 2.5c; the API surface is added but the only caller is the test for invalidation.)

**4. Builders delegate to `cell_balance_terms`** (Pattern 2 anchor). The cache populators consume the same algebra the per-cell `update(visit, ...)` does — same `denom` formula, same `a` formula. The L1 dual-view test asserts agreement at `rtol=1e-13`.

**5. Naming.** Per cross-domain Frame 3 (JAX `lax.scan` API): the sweep is geometry-blind. The current names `_sweep_1d_cartesian` and `_sweep_1d_curvilinear` become ONE `_sweep_1d_unified` (or `apply_sweep_1d`). The 2D path `_sweep_2d_wavefront` stays separate (deferred Q2 per Step 2.6; different work-unit shape).

**6. `CellUpdate.affine_coefficients(visits, ...)` retires.** Step 2.5b's vectorised builder collapses to the cache populator. The per-cell `update(visit, ...)` and `residual(...)` survive as the human-legible L1 reference (and the degenerate-axis fallback).

**7. Math notation in code matches Lewis-Miller §5.3 transmission-emission pair.** Sphinx narrative (Step S) names the pair canonically and cites L-M §5.3 + Blelloch §1.5 (the numerical-recurrence side).

### Concept-collapse signature

Step 2.5c (with Q1 of 2.6 absorbed):

- Sweep paths: 3 (`_sweep_1d_cartesian`, `_curvilinear`, `_2d_wavefront`) → 2 (`_sweep_1d_unified`, `_2d_wavefront`).
- Iteration primitives: 3 (`dag_walk` alias, `iter_cell_visits`, `iter_cells_by_direction`) → stays 3 for now (Step 2.6 Q3 canonicalises after this); slow path uses them, fast path doesn't.
- DD per-cell API: `update`, `residual`, `affine_coefficients` → `update`, `residual` (Pattern 2 dual-view) + cache populator function.
- Cache concept: 1 new dataclass (`SweepCoefficientCache`) → 2 (`GeometryCoefficients`, `CollisionCache`) — the operator-algebra layering.

Net: ~14 concepts → 10. Modest collapse on count BUT massive collapse on **per-iteration work** (the 80% bottleneck moves to one-time SNSolver construction).

### Anti-recommendations (the agent MUST NOT do these)

- **Do NOT pack both immutability strata into ONE flat dataclass.** Smell #16 (per cross-domain-attacker memo): cache shape mixing immutability strata. Two dataclasses per decisions §1.
- **Do NOT precompute `b` or `bc_inflow` in the cache.** They are mutable per-sweep / per-inner-iteration by physical necessity. Cache only σ_t-immutable quantities.
- **Do NOT touch BC composition.** Per explorer §Q3 BC apply is already a single-tensor-op `LinearOperator` per realised type (`PermutationOperator.apply = np.take`). BC is per-sweep by reflective-BC dependency; no cache opportunity, no friction.
- **Do NOT delete `CellUpdate.update(visit, ...)`.** It is the L1 dual-view validator. The cache populator derives from `cell_balance_terms`; `update` consumes `cell_balance_terms`; both must agree at `rtol=1e-13` (the Pattern 2 contract).
- **Do NOT touch `_sweep_2d_wavefront`.** 2D wavefront has a genuinely different work-unit shape per Step 2.6's Q2 deferral; the cache pattern does not naturally extend to 2D anti-diagonal scheduling.
- **Do NOT introduce JAX, numba, or any non-numpy dependency.** Pure numpy. The `lax.scan` API match is vocabulary-only per Frame 3 ADOPT-VOCABULARY.
- **Do NOT skip the cache-invariance test.** A test must verify that `CollisionCache` is NOT rebuilt across inner iterations (instrument via a counter). Per cross-domain §"First test" for Frame 5.
- **Do NOT skip the dual-view consistency test.** The cache-driven `apply_sweep(...)` result must equal the per-cell `update(visit, ...)` reference iteration at `rtol=1e-13` over the parametrised geometry × ng × source-kind grid.

### Mechanism criteria (greppable)

After Step 2.5c lands:

1. `grep -rn "def affine_coefficients" orpheus/sn/spatial/diamond.py` returns NOTHING (retired).
2. `grep -rn "np.fromiter\|iter_cell_visits" orpheus/sn/sweep.py` returns NOTHING (per-cell loops gone from hot path).
3. `grep -rn "class GeometryCoefficients\|class CollisionCache" orpheus/sn/spatial/sweep_cache.py` returns TWO definitions.
4. `grep -rn "def _sweep_1d_cartesian\|def _sweep_1d_curvilinear" orpheus/sn/sweep.py` returns NOTHING (unified).
5. `grep -rn "def _sweep_1d_unified\|def apply_sweep_1d" orpheus/sn/sweep.py` returns ONE definition.
6. `_sweep_1d_unified` body is ≤ 80 lines.
7. `GeometryCoefficients.from_mesh_and_quad(sn_mesh)` body is ≤ 60 lines.
8. `CollisionCache.from_geometry(geom, sig_t)` body is ≤ 25 lines (it's mostly broadcasting).
9. Slab benchmark `nx=160 N=16 ng=4`: ≥ 10× faster than Step 2.5b microbench (target: ≤ 1.5 ms/sweep vs 15.43 ms Step 2.5b vs 30.85 ms Step 2.5).
10. Full `pytest tests/sn/ -q` runtime: < 5 min on the same machine (vs Step 2.5's 3 hours @ 28% killed, vs Step 2.5b's ~33 min projected).

### Test catalog — verify the cache invariants

In `tests/sn/spatial/test_sweep_cache.py` (NEW). All `@pytest.mark.l0`.

**Cache structure tests**:

1. `test_geometry_coefficients_built_at_construction` — `GeometryCoefficients.from_mesh_and_quad(sn_mesh)` returns frozen dataclass with all expected fields populated; shapes match.
2. `test_collision_cache_built_at_sigma_t_bind` — `CollisionCache.from_geometry(geom, sig_t)` produces shapes consistent with `geom`; values match the analytic formula `1/(g_streaming + σ_t·V)`.
3. `test_two_strata_independence` — `GeometryCoefficients` shape has NO `ng` axis (proves the strata are separate). `CollisionCache` shape has `(N, nx, ng)`.

**Cache-invariance test (cardinal)**:

4. `test_collision_cache_invariance_under_source_iteration` — run a Picard fixed-point on a homogeneous 1G slab. Instrument `CollisionCache.from_geometry` with a call counter. After convergence (5+ iterations), counter == 1. (Per cross-domain Frame 5 "first test".)
5. `test_geometry_coefficients_invariance_under_sigma_t_change` — after `solver.rebind_cross_sections(new_sig_t)`, `geom is unchanged_geom` (id check). Only `coll` rebuilt.

**Dual-view consistency tests (Pattern 2)**:

6. `test_cache_driven_sweep_matches_per_cell_update` — parametrised over geometry × ng × source kind (36 cases minimum). `apply_sweep_1d(geom, coll, ...)` result and `cell_update.update(visit, ...)` iteration agree at `rtol=1e-13`.
7. `test_cache_populator_matches_cell_balance_terms` — for any cell, `cell_balance_terms(st, ...)` and the cache's per-cell extracted (a, denom) values agree at `rtol=1e-14`.

**Performance gates**:

8. `test_slab_sweep_benchmark_under_2ms` — `nx=160 N=16 ng=4` slab sweep completes in ≤ 2 ms (target ≤ 1.5 ms). Marker `@pytest.mark.slow`; skipped by default but runs in CI.
9. `test_full_sn_suite_under_5min` — separate marker; CI-only; runs full `tests/sn/` and asserts total time < 300s.

**Production gates** (re-verify Step 2.5b + earlier work):

10. `test_l0_streaming_equilibrium_preserved_after_2_5c` — 26/26 PASS.
11. `test_regression_snapshots_bit_identical_after_2_5c` — 11/11 PASS at `rtol=1e-12`.
12. `test_pair_monoid_associativity_still_passes` — Step 2.5b's algebraic anchor stays green.

### Test pin (verbatim full paste-back per L12)

You MUST run AT THE END and paste FULL STDOUT into the closeout memo, inside code fences:

```bash
.venv/bin/python -m pytest tests/sn/spatial/test_sweep_cache.py -v
.venv/bin/python -m pytest tests/sn/spatial/test_ordinate_scan.py -v
.venv/bin/python -m pytest tests/sn/spatial/test_diamond.py -v
.venv/bin/python -m pytest tests/sn/spatial/test_streaming_equilibrium_curvilinear.py -q
.venv/bin/python -m pytest tests/sn/regression/ -v
time .venv/bin/python -m pytest tests/sn/ -q
grep -rn "np.fromiter\|iter_cell_visits" orpheus/sn/sweep.py
grep -rn "def affine_coefficients" orpheus/sn/spatial/diamond.py
wc -l orpheus/sn/sweep.py orpheus/sn/spatial/sweep_cache.py orpheus/sn/spatial/diamond.py
```

**The `time pytest tests/sn/ -q` must complete with a summary line.** L12: no "verified subsets". The full-suite runtime is the load-bearing performance gate; if it's still > 10 min, the cache isn't doing its job.

### Decision-point checkpoints (STOP and dispatch)

- "I want to pack `GeometryCoefficients` and `CollisionCache` into one dataclass for simplicity" → STOP. Smell #16. Decision §1 is locked.
- "I want to cache `b` or `bc_inflow`" → STOP. They're per-sweep mutable; decision §3 documents what's invariant at what cadence.
- "The cache-invariance test (#4) shows `from_geometry` is called multiple times" → STOP. Trace the call sites; one of the sweep paths is rebuilding the cache. Should be once at `__init__` only.
- "The dual-view test (#6) fails at `rtol=1e-13`" → STOP. Either the cache populator derivation is wrong OR `cell_balance_terms` has a hidden non-affine term. Trace algebraically.
- "Slab benchmark is < 5× faster than Step 2.5b" → STOP. Profile; the per-cell `np.fromiter` or `iter_cell_visits` is still in the hot path somewhere.
- "Full suite still takes > 10 min" → STOP. The unified sweep is still going through the slow path for some test class; identify and fix.

### Scope hard limits

IN scope:
- `orpheus/sn/spatial/sweep_cache.py` (NEW) — both dataclasses + populator.
- `orpheus/sn/sweep.py` — unify slab + curvilinear into `_sweep_1d_unified` consuming the cache.
- `orpheus/sn/spatial/diamond.py` — retire `affine_coefficients`; the cache populator subsumes it.
- `orpheus/sn/spatial/cell_update.py` — remove `affine_coefficients` from Protocol.
- `orpheus/sn/solver.py` — build the two caches at `__init__`; expose `rebind_cross_sections` API.
- `tests/sn/spatial/test_sweep_cache.py` (NEW) — 12 tests above.

OUT of scope:
- 2D wavefront / `update_batch` / `_sweep_2d_wavefront` — Q2 deferred.
- `dag_walk` canonicalisation — Step 2.6 Q3 still pending after this.
- Step 3's L+C operator-class split (this work prepares the cache; the LinearOperator wrapper comes in Step 3).
- Multi-physics depletion / thermal feedback consumers of `rebind_cross_sections` — only the API surface is added.

### Commit scope (≤ 4 commits)

1. `feat(sn): GeometryCoefficients + CollisionCache two-stratum cache (Issue #196 Step 2.5c)` — `sweep_cache.py` + cache-structure tests (#1-3) + invariance tests (#4-5).
2. `feat(sn): apply_sweep_1d unified body consuming the two-stratum cache (Issue #196 Step 2.5c)` — unify `_sweep_1d_cartesian` + `_sweep_1d_curvilinear`; consume cache; dual-view tests (#6-7).
3. `refactor(sn): retire DiamondDifference.affine_coefficients; cache populator subsumes it (Issue #196 Step 2.5c)` — Protocol + ABC + class cleanup.
4. `test(sn): performance gates + production-gate re-verification (Issue #196 Step 2.5c)` — tests #8-12.

### Closeout memo

Path: `.claude/agent-memory/method-implementer/issue_196_phase_g_step2_5c_closeout.md`.

Required:
1. **Verbatim full pytest stdout** for ALL EIGHT commands in §"Test pin". NO subsets.
2. **`time pytest tests/sn/ -q`** wall-clock total — the load-bearing performance number.
3. **Slab benchmark microbench**: `nx=160 N=16 ng=4` ms/sweep before (Step 2.5b: 15.43 ms) vs after.
4. **Cache invariance verification**: test #4 PASS line showing `from_geometry` called exactly once.
5. **Dual-view consistency**: sample one parametrised case from test #6 with actual rtol achieved.
6. Files touched, numbered list.
7. **What this does NOT close**: Step 2.6 Q3 (`dag_walk` canonicalisation) still pending; Step 3+.
8. Decision-point honesty — any STOP hit, document.

### Sub-agent dispatch

method-implementer with brief shaped per §"Brief template". Brief MUST cite BOTH explorer + cross-domain memos. MUST demand verbatim full `pytest tests/sn/ -q` paste-back. MUST forbid flat single-dataclass cache (Smell #16).

---

## Step 2.6 — `dag_walk` canonicalisation (Q3 only — Q1 subsumed by 2.5c)

**Status (2026-05-14)**: Q1 (sweep-skeleton unification) SUBSUMED by Step 2.5c — the two-stratum cache abstracts geometry differences as data, so `apply_sweep_1d(geom, coll, ...)` is geometry-blind by construction. Q3 (canonicalise `dag_walk`; delete `iter_cell_visits` + `iter_cells_by_direction`; migrate 6 call sites) REMAINS PENDING after Step 2.5c. With the cache in place, the fast path doesn't use the iter methods at all (cache populator does the gather ONCE); Q3 becomes mostly cleanup — delete the legacy methods + migrate the slow-path / degenerate-axis fallback callers.

The original Step 2.6 locked-decision text (covering both Q1 and Q3) follows for archival reference:

**Goal**: complete the structural collapse Step 2.5 left on the table. Per the further-collapse explorer memo at `.claude/agent-memory/explorer/issue_196_phase_g_step2_5_further_collapse.md`, two of the four candidate unifications are achievable wins; the other two are correct Pattern 6 deferrals. Step 2.6 executes the two achievable wins:

- **Q1**: collapse `_sweep_1d_cartesian` + `_sweep_1d_curvilinear` into `_sweep_1d_unified`.
- **Q3**: canonicalise `dag_walk` as the iteration primitive; delete `iter_cell_visits` + `iter_cells_by_direction`.

Q2 (2D wavefront join 1D skeleton) and Q4 (collapse `update` + `update_batch`) stay deferred — fold-with-accumulator vs fold-with-parallel-reduce are structurally different patterns; forcing the unification would slow 1D ~3-5x. The shared algebra IS already captured (`cell_balance_terms`).

**Concept-count signature post-2.6**:
- Sweep paths: 3 (`_sweep_1d_cartesian` / `_curvilinear` / `_2d_wavefront`) → 2 (`_sweep_1d_unified` / `_2d_wavefront`).
- Iteration methods: 3 (`dag_walk` alias / `iter_cell_visits` / `iter_cells_by_direction`) → 1 (`dag_walk` canonical).
- Total: 11 → 6 (~1.8:1 collapse). Honest realistic floor.

### Pre-step: read these files in full

1. The further-collapse explorer memo: `.claude/agent-memory/explorer/issue_196_phase_g_step2_5_further_collapse.md`. Reads the Q1 unified sweep sketch (≤ 30 lines) and the Q3 canonicalisation migration plan.
2. The two existing 1D sweep bodies that collapse: `orpheus/sn/sweep.py:168-289` (`_sweep_1d_cartesian`) and `orpheus/sn/sweep.py:293-523` (`_sweep_1d_curvilinear`). Identify the stage-by-stage commonalities listed in the memo's §Q1 anatomy table.
3. The three iteration methods on `SNMesh`: `dag_walk` at `geometry.py:425-466` (alias today), `iter_cell_visits` at `geometry.py:468-628`, `iter_cells_by_direction` at `geometry.py:631-711`. The 6 call sites: `operator.py:787, 820, 1027, 1074` (matvec); `sweep.py:263, 483` (SI).
4. The Carlson seed call site at `sweep.py:431` (the only stage that genuinely differs between slab and curvilinear in the unified body).
5. The lessons learned in `.claude/lessons.md` §§ L12, L13, L14. L12 paste-back discipline is mandatory; L13 names existing types to extend; L14 demands 3-way standoff under refinement.

### Locked design decisions

1. **`_sweep_1d_unified(Q, sig_t, sn_mesh, psi_bc, Q_aniso)` replaces both 1D sweep bodies**. The skeleton:
   - Determine `levels`: `[None]` for slab/sphere, `quad.level_indices` for cylinder (existing per `quad.level_indices` attribute).
   - For each `(p_idx, level)` in `enumerate(levels)`:
     - Compute `psi_angle = _seed_angular_state(sn_mesh, level, psi_bc, sig_t, Q)` returning `None` for slab, the per-level Carlson seed result for sphere/cylinder.
     - Iterate ordinates in this level (slab/sphere: all N; cylinder: just `quad.level_ordinate_indices[level]`).
     - For each ordinate `n`, run the inner fold (already identical between the two paths).
   - Cache writes guarded by `if reduced.requires_upstream_angular_state` for pole_key + phi_prev_key.

2. **Factor `_seed_angular_state(sn_mesh, level, psi_bc, sig_t, Q) -> np.ndarray | None`** as a free function in `sweep.py`. Returns `None` for slab (`reduced.requires_upstream_angular_state == False`). For sphere/cylinder, calls the existing `carlson_inward_sweep_from_source` with the level-specific arguments. This is Pattern 2 (single source of truth) — the per-geometry decision lives in one place, not scattered across two sweep bodies.

3. **`dag_walk` becomes the canonical primitive**. New signature: `dag_walk(self, *, ordinate_idx: int | None = None, direction_sign: int | None = None, mu_level_idx: int | None = None) -> Iterator[CellVisit]`. Exactly one of `ordinate_idx` or `direction_sign` is required (XOR; assert this). The body subsumes the logic of both legacy methods.

4. **Delete `iter_cell_visits` and `iter_cells_by_direction`** from `geometry.py` after migrating the 6 call sites to `dag_walk`. No backward-compat aliases — the API is internal.

5. **Migrate 6 call sites**:
   - `sweep.py:263` (in current `_sweep_1d_cartesian`) and `sweep.py:483` (in `_sweep_1d_curvilinear`) become call sites in `_sweep_1d_unified` using `dag_walk(ordinate_idx=n, mu_level_idx=...)`.
   - `operator.py:787, 820, 1027, 1074` (matvec) call `dag_walk(direction_sign=...)`.

6. **Bit-identity contracts preserved.** Q1 + Q3 are pure structural reorganisation — the math is unchanged. All 11 regression snapshots stay bit-identical at `rtol=1e-12`. If a snapshot drifts beyond ULP × iteration count, the refactor introduced an operation-order change; trace and fix.

### Anti-recommendations (the agent MUST NOT do these)

- **Do NOT keep `iter_cell_visits` or `iter_cells_by_direction` as aliases.** Decision #4 deletes them. The 6 call sites migrate.
- **Do NOT fold the 2D wavefront into `_sweep_1d_unified`.** Q2 is correctly deferred per Pattern 6. `_sweep_2d_wavefront` stays as its own body — its work-unit shape is genuinely different.
- **Do NOT collapse `update` and `update_batch`.** Q4 is correctly deferred per Pattern 6. Same reasoning — different algorithm patterns, not different control flow.
- **Do NOT introduce a `WorkUnit` Union or duck-typed work units.** `dag_walk` yields `CellVisit` (the existing 1D type); the 2D path uses `SweepCellSlice` via its own iteration. Different paths, different work-unit types, both correct.
- **Do NOT inline the Carlson seed call into `_sweep_1d_unified`.** Factor it into `_seed_angular_state` per decision #2. Single source of truth.
- **Do NOT alter the Step 2 canonical-closure math.** Step 2.6 is structural reorganisation only. If the L0 streaming-equilibrium test regresses, you've changed math; trace, don't widen rtol.

### Mechanism criteria (greppable)

After Step 2.6 lands, ALL of the following must be true:

1. `grep -rn "def _sweep_1d_cartesian\|def _sweep_1d_curvilinear" orpheus/sn/sweep.py` returns NOTHING. One unified body in place.
2. `grep -rn "def iter_cell_visits\|def iter_cells_by_direction" orpheus/sn/geometry.py` returns NOTHING. `dag_walk` is the only iteration method.
3. `grep -rn "iter_cell_visits\|iter_cells_by_direction" orpheus/` returns NOTHING (all call sites migrated).
4. `_sweep_1d_unified` body is ≤ 80 lines per the explorer §Q1 sketch.
5. `_seed_angular_state` body is ≤ 30 lines.
6. `dag_walk` body is ≤ 50 lines (subsuming the two legacy methods).
7. 11 regression snapshots bit-identical (`rtol=1e-12`) — Step 2.6 is structural-only.
8. L0 streaming-equilibrium test (`tests/sn/spatial/test_streaming_equilibrium_curvilinear.py`) 26/26 PASS at `rtol=1e-9`.

### Test pin (verbatim paste-back required per L12)

You MUST run these commands at the end and PASTE FULL STDOUT into the closeout memo, inside code fences:

```bash
.venv/bin/python -m pytest tests/sn/ -q  # full suite, NO selective subsets
.venv/bin/python -m pytest tests/sn/spatial/test_streaming_equilibrium_curvilinear.py -q
.venv/bin/python -m pytest tests/sn/regression/ -v
grep -rn "def _sweep_1d_cartesian\|def _sweep_1d_curvilinear" orpheus/sn/sweep.py
grep -rn "def iter_cell_visits\|def iter_cells_by_direction" orpheus/sn/geometry.py
grep -rn "iter_cell_visits\|iter_cells_by_direction" orpheus/
wc -l orpheus/sn/sweep.py orpheus/sn/geometry.py
```

**The Step 2.5 closeout's "verified subsets" pattern is the L12 failure mode. DO NOT do that here. Paste the full `pytest tests/sn/ -q` output, including the final summary line.**

### Decision-point checkpoints (STOP and dispatch)

- "Slab adopting `levels = [None]` breaks the slab regression snapshots" → STOP. The wrapper adds no math; if it changes results, you've changed operation order. Trace.
- "`dag_walk(direction_sign=...)` and `dag_walk(ordinate_idx=...)` produce different visit sequences for matvec vs SI" → STOP. They should produce identical visit shapes; matvec was always direction-keyed because the per-ordinate `streaming_terms` differences are handled inside `DiamondDifference.update`, not in the iteration order.
- "The L0 streaming-equilibrium test regresses after the unified body lands" → STOP. The Step 2 canonical closure must not change; if it did, you altered the math.
- "Picard convergence on `_sweep_1d_unified` slows compared to the two-body version" → STOP. Investigate; do NOT widen tolerances.

### Scope hard limits

IN scope:
- Q1 (unified 1D sweep skeleton + `_seed_angular_state` factoring).
- Q3 (canonical `dag_walk` + delete legacy two methods + migrate 6 call sites).

OUT of scope (do NOT touch):
- Q2 (2D wavefront unification) — Pattern 6 deferred.
- Q4 (`update` + `update_batch` collapse) — Pattern 6 deferred.
- Step 3's `StreamingOperator` + `CollisionOperator` split.
- Step 4's `SourceIteration` / `KEigenvalue` consumers.
- Any Step 2 canonical-closure math changes.

### Commit scope (≤ 3 commits)

1. `refactor(sn): _sweep_1d_unified replaces _sweep_1d_cartesian + _sweep_1d_curvilinear via levels=[None] + _seed_angular_state factor (Issue #196 Step 2.6 Q1)` — the unified body + the seed-factor helper.
2. `refactor(sn): canonicalise dag_walk; delete iter_cell_visits + iter_cells_by_direction; migrate 6 call sites (Issue #196 Step 2.6 Q3)` — the API consolidation.
3. (optional) `chore(sn): regenerate any snapshots that drifted at ULP-level after operation-order change (Issue #196 Step 2.6 follow-up)` — only if any snapshot drifts beyond bit-identity; with three-pillar attestation per `vv-principles`.

### Closeout memo

Path: `.claude/agent-memory/method-implementer/issue_196_phase_g_step2_6_closeout.md`.

Required content:
1. **Full pytest stdout pasted verbatim** for the three commands listed in Test pin. No "verified subsets". The L12 failure mode is explicit and recurrent; do NOT repeat it.
2. **Concept-count audit** — 11 → 6, with file:line citations for what was deleted and what replaced it.
3. Files touched, numbered.
4. What this does NOT close: Step 3 (L pure / C separate), Step 4 (SourceIteration consumers), Step 5 (BC audit), Step 6 (`.H` + adjoint).
5. **Q2 + Q4 deferral justification** with explicit Pattern 6 citation.

### Sub-agent dispatch

method-implementer with brief shaped per the §"Brief template" at the bottom of this plan file. Brief MUST cite the further-collapse explorer memo. Brief MUST forbid Q2 and Q4. Brief MUST demand verbatim full `pytest -q` paste-back, not "verified subsets".

---

## Step 3+4 (Unified) — Within-group operator restructure under `S = S_foldable + S_residual`

**Status (2026-05-14)**: DESIGN LOCKED, IMPLEMENTATION PENDING. Replaces the prior separate Step 3 (L/C split) + Step 4 (SourceIteration/KEigenvalue wiring). Bundled per user direction because the same files are touched and "new code never wears the old names". Full Step 3+4 design is in the plan-mode artefact at `~/.claude/plans/crystalline-wondering-token.md` (lines 369-833); summary follows.

### Locked design decisions (user 2026-05-14)

1. **Split scattering as $S = S_{\text{foldable}} + S_{\text{residual}}$** — only P0 within-group folds into $\sigma_r$; Pℓ≥1 self-scatter is inherently anisotropic via $Y_\ell^m(\Omega_n)$. Not an SN peculiarity — inherent to scattering's angular structure for ALL transport solvers.
2. **NO `WithinGroupOperator` class.** Use pure algebra `A_wg = L + C - S.foldable_part()` (an `OperatorSum`). Local variable name `A_wg`.
3. **Local variable `scatter_cycle = A_wg.inverse() @ S_residual`** (math symbol $\mathcal{C}$). Do NOT use `T_off` (clashes with future time operator).
4. **Bundle Step 3+4 into 5 substeps** sharing files; the cosmetic + architectural changes happen together.
5. **Drop F universally from within-group accelerators** (`SourceIteration` / `KrylovAcceleration`). F belongs at Layer 3 (`KEigenvalue` outer) for all transport solver families (SN, CP, MoC, diffusion). The existing `F=ZeroOperator()` papering-over in `KEigenvalue` (`iteration.py:519-531`) gets retired.
6. **Dispatch by accelerator CLASS, not string flag.** `inner_solver=SourceIteration` or `inner_solver=KrylovAcceleration`. Pattern 4.
7. **Aggressive retirement is a deliverable.** Substep 3+4.e is a codebase-wide retirement audit. Old code = noise (memory: `feedback_aggressive_retirement`).

### CORRECTION (2026-05-14, post-revert of `ad37ca0`) — Resolution A: subtractive L decomposition

**Authoritative memo**: `.claude/agent-memory/numerics-investigator/sn_LC_decomposition_derivation.md` (SymPy + numpy derivation).

The first 3+4.b.i attempt implemented `StreamingOperator.apply(ψ) = matvec(ψ, σ_t=0)`. SymPy-derivation proves this is mathematically wrong for **curvilinear**: the Hébert §3.9.4 Eq. (3.434) Carlson coupled-pole seed `phi_aux_i = (dr·σ_t·φ_0/Σw + 2·phi_face) / (dr·σ_t + 2)` is **rational** (not affine) in σ_t — the discrete L (with M-M angular closure) is σ_t-coupled by construction. Cartesian IS affine; curvilinear is not. Empirically `matvec(ψ, σ_t=0)` differs from pure-L by 3-13% rel for sphere/cylinder on random ψ.

**Resolution A (the canonical answer)**: define L subtractively, both L and C carry σ_t at constructor.

```
L: StreamingOperator(sn_mesh, sigma_t).apply(ψ) := M(ψ; σ_t_full) − σ_t_packed ⊙ ψ
C: CollisionOperator(sn_mesh, sigma_t).apply(ψ) := σ_t_packed ⊙ ψ
```

Then `(L + C).apply(ψ) = M(ψ; σ_t_full)` by construction (bit-exact, rel_residual = 0.0 across slab/sphere/cylinder × 3 random seeds, verified in `derivations/diagnostics/diag_LC_decomposition_resolution.py`).

The continuous L (`Ω·∇ψ + (1-μ²)/r · ∂_μ ψ`) is σ_t-independent. The **discrete L** (with Hébert's M-M angular closure) carries σ_t through the Carlson seed — this is a property of the closure choice, not an implementation defect. Analogous to MoC's `exp(-σ_t·s)` characteristic line and the DD coefficient `α_DD(σ_t·dx)`: discrete streaming operators routinely carry parameters their continuous form does not.

**Step 6 impact**: PRESERVED. `(L+C-S).H.apply(x)` distributes via `OperatorSum.apply_transpose` → `L.H.apply(x) − σ_t·x + σ_t·x − S.H.apply(x) = matvec_transpose(x; σ_t) − S.H.apply(x)`. The reverse-direction sweep already required by Step 6 trivially extends with the same subtractive overlay.

**The pseudo-code below (and all subsequent substep specifications) is corrected accordingly**: both L and C take σ_t at constructor; the factory `build_sn_operators(sn_mesh, materials)` ensures both leaves share the same σ_t.

### Math architecture

**Three layers**, F at exactly one:

| Layer | Class / function | Inputs | Math |
|---|---|---|---|
| 3 — eigenvalue outer | `KEigenvalue(L, S, F, *, accelerator_cls)` | `L=A_wg, S=S_residual, F=F` | $(L-S)\psi = F\psi/k$; F enters via `q_fission = F.apply(ψ)/k`. |
| 2 — within-group accelerator | `SourceIteration(L, S)` or `KrylovAcceleration(L, S)` (siblings) | `L=A_wg, S=S_residual` (no F!) | $(I - L^{-1}S)\psi = L^{-1}q_{\text{ext}}$. Picard or GMRES. |
| 1 — per-step primitive | `step(ψ) = inverter(S.apply(ψ))` | shared across L2 accelerators | $\mathcal{C}\psi = A_{\text{wg}}^{-1}\,S_{\text{residual}}\,\psi$. |

**Pseudo-code**:

```python
# Build the algebra at SNSolver.__init__ — pure algebra, no wrapper classes.
# Both L and C carry σ_t at constructor (Resolution A — see CORRECTION above);
# the factory ensures consistency.
L, C, S, F = build_sn_operators(sn_mesh, materials)   # both L and C receive σ_t
A_wg       = L + C - S.foldable_part()         # OperatorSum (math: A_wg)
S_residual = S.residual_part()

# k-eigenvalue route
keff, k_hist, psi = KEigenvalue(
    L=A_wg, S=S_residual, F=F,                  # F at Layer 3 only
    accelerator_cls=SourceIteration,            # or KrylovAcceleration
    inverter=A_wg.solve,
).solve(initial_guess=...)

# Fixed-source route
psi, _ = SourceIteration(
    L=A_wg, S=S_residual, inverter=A_wg.solve,  # no F
).solve(q_ext=q_external)
```

**Convergence**: $\rho(\mathcal{C}) \ll 1$ in essentially all practical regimes. For 1G, $S_{\text{residual}} = 0$ → SI converges in 1 step. For multi-G downscatter, $S_{\text{residual}}$ is strictly triangular → exact in $N_g$ steps. For multi-G upscatter, $\rho \approx$ upscatter ratio — fast.

### Substeps (5 commits)

- **3+4.a** — `ScatteringOperator.foldable_part() / .residual_part() / .foldable_sigma()`. Tests for the algebraic identity `S = S.foldable_part() + S.residual_part()` at `rtol=1e-14`.
- **3+4.b** — `StreamingOperator(sn_mesh, σ_t)` + `CollisionOperator(sn_mesh, σ_t)` leaves under Resolution A (both carry σ_t; L.apply uses subtractive `M(ψ; σ_t) − σ_t⊙ψ`; bit-exact decomposition); `ScatteringOperator.is_foldable_into_sigma_r()` detection; `OperatorSum.solve` fusion hook for the within-group shape (reads σ_t from L/C, builds lazy σ_r cache, routes to `sweep_within_group_1d`); rename `apply_sweep_1d → sweep_within_group_1d`; new `apply_within_group_1d` mirror in `orpheus/sn/sweep.py`. **Split into 3+4.b.i (leaves + foldability detection) and 3+4.b.ii (fusion hook + sweep rename + apply mirror)** for gate-keeper-friendly blast radius. The reverted commit `ad37ca0` is the WRONG approach (`L.apply = matvec(σ_t=0)` is non-physical for curvilinear); see CORRECTION above.
- **3+4.c** — Refactor `SourceIteration` to drop F (signature becomes `(L, S, *, inverter, ...)`); new `KrylovAcceleration` sibling with same shape; `KEigenvalue.accelerator_cls` parameter; retire 4× `_solve_*` + 3× `transport_operator_matvec_*` + 3× `_build_rhs_*` + `SNStreamingOperator` + `power_iteration` use in SN. `solve_sn(..., inner_solver=SourceIteration)` is the public dispatch.
- **3+4.d** — Iteration-count gates (c=0.95 1G cylinder SI → 1 sweep; Krylov → 1 outer × ≤30 inner); 2G upscatter convergence test; perf microbench; full suite `time pytest tests/sn/ -q` < 30 min.
- **3+4.e** — Codebase-wide retirement audit; deletion of aliases/shims/dead tests per `feedback_aggressive_retirement`.

### Reference artefacts (durable)

- **Plan-mode artefact** (full Step 3+4 design, ~470 lines): `~/.claude/plans/crystalline-wondering-token.md` lines 369-833.
- **Explorer memos** (plan-mode artefacts; promote to `.claude/agent-memory/explorer/` on first use):
  - Scattering split audit: `~/.claude/plans/crystalline-wondering-token-agent-a63c996cb3094b0fc.md`.
  - Apply-path + solver coupling audit: `~/.claude/plans/crystalline-wondering-token-agent-ae365fb943dbedc39.md`.
  - Naming sweep audit: `~/.claude/plans/crystalline-wondering-token-agent-aa70cfd7f0e605015.md`.
- **Numerics-investigator memo** (Krylov 27× SLOWER pathology): `.claude/agent-memory/numerics-investigator/krylov_inner_solver_profile_2026_05_14.md`.
- **Nexus bug report** (`retest --scope=branch` returns 0): `.claude/agent-memory/numerics-investigator/nexus_retest_branch_scope_bug.md` — for separate next-session resolution.
- **Memory** (new entries 2026-05-14):
  - `feedback_aggressive_retirement` — old code = noise; retirement is a deliverable.
  - `feedback_unused_args_check_solver_families` — when deciding "keep unused arg for future" vs "drop as dead weight", check math layering across all solver families.

### Anti-recommendations (the agent MUST NOT do these)

- Do NOT introduce `WithinGroupOperator`, `IterationOperator`, or `ScatterCycleOperator` classes — algebra suffices.
- Do NOT use `T_off` or `T` as iteration-operator variable name. Use `scatter_cycle` / 𝒞.
- Do NOT add an `accelerator="picard"|"gmres"` string flag — dispatch by class.
- Do NOT keep F as a constructor argument on within-group accelerators. F is at Layer 3 only.
- Do NOT pass `F=ZeroOperator()` or `F=None` to within-group accelerators — refactored constructors have no F slot.
- Do NOT keep `SNStreamingOperator` / `transport_operator_matvec_*` as deprecation aliases. Delete outright.
- Do NOT add P1+ folding logic — only ℓ=0 within-group folds into σ_r.
- Do NOT regenerate curvilinear snapshots without three-pillar attestation.
- Do NOT skip the retirement audit (3+4.e).

### Verification gates

1. Algebraic identity `S ≡ S.foldable_part() + S.residual_part()` at `rtol=1e-14`.
2. Sweep/apply symmetry `(L+C-S_foldable).apply((L+C-S_foldable).solve(q)) == q` at `rtol=1e-12`.
3. F at exactly one layer: `grep "F\.apply\|F\.solve" orpheus/numerics/iteration.py` matches only inside `KEigenvalue.solve`.
4. Regression bit-identity 11/11 at `rtol=1e-12`.
5. L0 streaming-equilibrium 26/26 PASS on both accelerators.
6. Pair-monoid associativity stays green.
7. L1 SN bridge test passes WITHOUT adapter shims.
8. Iteration-count gates: 1G c=0.95 cyl SI → 1 sweep; Krylov → 1×≤30; 2G upscatter ≤ 5 outers @ ρ ≤ 0.2.
9. Perf gate: `time pytest tests/sn/ -q` < 30 min (pre-3+4: 2h28m).
10. Sphinx -W clean.

### Sub-agent dispatch

method-implementer with brief shaped per §"Brief template". Brief MUST cite this Step 3+4 section, the plan-mode artefact for the full design, the three explorer memos, the Krylov investigator memo, and lessons L12-L16.

### Closeout memo

Path: `.claude/agent-memory/method-implementer/issue_196_phase_g_step3_4_closeout.md`. Required: verbatim full pytest stdout for each substep + final full suite; iteration-count + wall-clock measurements; concept-count audit; retirement audit findings; decision-point honesty (fire/no-fire for each STOP gate).

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
