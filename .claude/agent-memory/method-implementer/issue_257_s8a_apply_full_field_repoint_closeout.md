---
name: issue-257-s8a-apply-full-field-repoint-closeout
description: #257 S8a — re-point every SN operator matvec LEAF (apply/apply_transpose) + the loss_action/loss_action_transpose producers + the M[σ] op from emitting history-bearing TimedFullField to the TIMELESS FullField (cofree base arrow); driver re-attaches via TimedFullField.__add__ recombine (NO advance needed — advance is dead in production). VALUE-NEUTRAL / byte-identical; 0 net-new pyright; the InvertibleOperator.solve resolvent STAYS TimedFullField.
metadata:
  type: project
---

# #257 S8a — apply(FullField) re-point + helper widening (cofree base-arrow codomain)

Branch `feature/field-typed-operator-algebra` (on S8-HEAD `90ea1cd`). **NOT committed.**
VALUE-NEUTRAL (a typing + codomain change only). S8b (`(L+C)−C` drop → pure-L) +
S8c (S/F fibration) NOT done — explicitly out of scope.

## The thesis (cofree finding, the "type error of altitude")
`TimedFullField = Cofree(FullField, d)`. An operator matvec leaf is a base arrow
`FullField → FullField` and must NOT carry a history tail — only the iteration
DRIVER sees the comonad. S8a re-points the leaf OUTPUTS (apply/apply_transpose) +
the `loss_action`/`loss_action_transpose` producers they delegate to + the S3b
`MultiplicationOperator` from emitting `TimedFullField` (with `history_depth=psi.history_depth`)
to the timeless `FullField` (history-free).

## ⭐ STEP-0 ENTANGLEMENT VERDICT = PROCEED (the load-bearing finding)
The driver's history machinery is **INERT in production**: `grep -rn "\.advance("` /
`at_lag` / `_history` across `orpheus/` shows **ZERO production callers** of `advance`
(steady-state SI/Krylov carry `history_depth=2` but `_history=()` always, never advance).
So the re-point is PURE type-plumbing — NOT behavior entanglement. The minimal correct
thing: re-point the **apply/apply_transpose matvec leaves** (the base arrows); the driver
re-attaches the timed type AUTOMATICALLY via `TimedFullField.__add__`'s recombine when a
timeless gain output is added to the timed `rhs` (`rhs = q_ext + S.apply(ψ) + B.apply(ψ)`,
`q_ext` is timed → `TimedFullField.__add__(timed, timeless)` → `_recombine` → timed).
**NO `advance` call is added** (it would be semantically wrong: a gain output is a source
increment, not a state advance; and there's no time-stepping consumer).

## The leaf/resolvent split (the precise S8a boundary, gate spec item #2 ∩ C5(a))
- **apply / apply_transpose = base arrows → re-pointed to FullField** (`StreamingOperator`,
  `InvertibleOperator`, `MultiplicationOperator`/`CollisionOperator`, `SNBoundaryOperator` B,
  the S/F `TimedFullField`-input `@apply.register` arm OUTPUT).
- **`InvertibleOperator.solve` = the iterate-PRODUCING resolvent → STAYS `TimedFullField`**
  (the driver's comonad-carrying iterate; covariant `TimedFullField <: FullField` so binding
  the class on `["FullField"]` is Liskov-safe). solve does NOT call `loss_action` — it builds
  its own `TimedFullField` at the sweep return (`operator.py:~1154`), untouched.
- **`MultiplicationOperator.solve`** (leaf q/f) → FullField (it's a leaf op, not the resolvent).
- Krylov return stays `TimedFullField` independently (templated on `initial_guess` via
  `_unravel_like`/`from_flat`, NOT via operator output type).

## Production diff (7 files — file:line summary)
1. `orpheus/sn/operator.py` — class bindings `StreamingOperator(LinearOperatorMixin["FullField"])`
   + `InvertibleOperator(OperatorSum["FullField"])`; `_require_typed_composite` guard WIDENED
   (`TimedFullField`→`FullField` isinstance, msg "expected FullField" — a `TimedFullField` IS a
   `FullField` so drivers still pass); `StreamingOperator.apply`/`apply_transpose` +
   `InvertibleOperator.apply`/`apply_transpose` sigs `(FullField)->FullField`, bodies build
   `FullField(...)` (drop `_history`/`history_depth`). **`(L+C)−C` arithmetic at :431/473
   UNTOUCHED** (`lpc.bulk.values - self.sigma_t[None]*psi.bulk.values`). Added `FullField` to
   TYPE_CHECKING.
2. `orpheus/sn/loss_representation.py` — 15 `psi/phi:"TimedFullField"` params + 14
   `-> "TimedFullField"` returns (all loss_action-family + `_apply_walk`) → `FullField`; 3
   construction sites (`_OctantWalk.loss_action` / `_OneDimScanWalk._apply_walk`'s `loss_action`
   / curvilinear `loss_action_transpose`) `TimedFullField(...)`→`FullField(...)` (drop history) +
   their 3 local imports. The `initial_guess:"AngularFlux|TimedFullField|None"` annotations + the
   TYPE_CHECKING `TimedFullField` import STAY (the sweep/solve path, not loss_action).
3. `orpheus/transport/multiplication_operator.py` — class `LinearOperatorMixin["FullField"]`;
   apply/solve/apply_transpose `(FullField)->FullField`, bodies build `FullField` (drop history);
   TYPE_CHECKING `TimedFullField`→`FullField`.
4. `orpheus/sn/fission.py` — the `@apply.register def _(self, psi: TimedFullField)` arm RETURN
   `-> "FullField"` (INPUT dispatch annotation `TimedFullField` UNCHANGED = fibration preserved);
   body `FullField(...)`; hoisted `from ...full_field import FullField` to module level.
5. `orpheus/sn/scattering.py` — same as fission (the `TimedFullField`-input arm output → FullField);
   module-level `FullField` import.
6. `orpheus/sn/boundary_operator.py` (B) — `_apply_faces`/`apply`/`apply_transpose`
   `(FullField)->FullField`, `_apply_faces` builds `FullField` (drop history); removed the
   now-unused TYPE_CHECKING `TimedFullField` import, added `FullField`.
7. `orpheus/numerics/spaces/full_field_space.py` — **THE ONE FILE BEYOND THE BRIEF'S LEAF
   SET** (necessary, see SURPRISE below): `FullFieldSpace._rebuild` (the G-adjoint metric helper)
   was hardcoded `replace(x, ..., _history=())` → broke when `x` is now a timeless `FullField`
   operator output (no `_history` field). FIX: route through the polymorphic `x._recombine(...)`
   hook (the algebra-lives-once hook FullField+TimedFullField both expose) — resolves it for
   BOTH carriers, principled (it's exactly what `_recombine` is for). Resolved 7 `test_g_adjoint_reciprocity` failures.

## NEW catcher — C5 (the deliverable I wrote)
`tests/sn/operators/test_apply_full_field_codomain.py` (`@pytest.mark.foundation`, 14 pass + 1 xfail,
all `-O`-firing via `pytest.fail`/`np.testing`/`pytest.raises` — Mode 8 clean):
- **C5a** — every leaf `L/C/S/F/B.apply` (+ `L/C/B.apply_transpose`) returns `type(out) is FullField`
  + `not isinstance(out, TimedFullField)` + no `history_*` attr, per geometry, for every input
  `history_depth`. **Mode-11**: calls `op.apply` DIRECTLY (the matvec leaf has ZERO graph callers).
  **SENTINEL-CONFIRMED** (in-process method-wrap, then removed): C5a fires the re-typed
  `StreamingOperator.apply` leaf, emitting `FullField` (count=1).
- **C5b** — driver re-attach byte-identity: `solve_sn` eigenvalue vs closed-form `k_inf=νΣf/Σa`
  (struct-indep, per geometry × {SI,krylov}; sphere/krylov xfail = #200) + a direct
  `SourceIteration(LC,S,B)` check proving the iterate stays `TimedFullField` (flux warm-start,
  mirroring `SNSolver._solve_source_iteration`).
- **C5c** — `TimedFullField.advance` type-guard still fires (`pytest.raises(TypeError, "new_bulk type")`).

## B legacy-pin re-points (the LARGE scope finding — 41 tests across 9 files broke by design)
Every test asserting the OLD codomain broke: `isinstance(out, TimedFullField)` /
`out.history_depth` / `out._history` on an operator OUTPUT, OR `dataclasses.replace(out,_history=())`.
ALL re-pointed minimally (values byte-identical) — `→ isinstance(out, FullField)` /
`not isinstance(out, TimedFullField)`; the `history_depth`-preservation tests RENAMED to
`test_output_is_timeless_full_field` (the new contract). Files: `test_operators_apply_typed`(10),
`test_g_adjoint_reciprocity`(7, fixed via the prod `_rebuild` change — NO test edit needed),
`test_collision_operator`(2+the earlier-rewritten history test), `test_streaming_operator`(several
+ guard-msg regex `expected TimedFullField`→`expected FullField`), `test_fission_operator`,
`test_scattering_operator`, `test_native_matvec`, `test_bc_extraction_matvec` (`_LpC_apply` return
annot + helper), `test_invertible_operator` (round-trip `LC.solve(LC.apply(ψ))` now wraps the
FullField source into a TimedFullField rhs — production does this), `test_removal_form_matvec_sweep`
(same wrap), `test_krylov_curvilinear_precond_safety` (`psi_typed_warm` annot `TimedFullField|None`→`FullField|None`).

## GATES (all green)
- **pyright `2297 / 19w` = EXACT baseline, 0 net-new, 0 net-new `# type: ignore`.** Measured
  baseline via a `git worktree add -d /tmp/s8a_baseline HEAD` + `.venv` SYMLINK (L28-safe; never
  touched the working tree). The invariance-seam errors at the brief's `operator.py:641/794` had
  ALREADY drifted/resolved pre-S8a (the live op.py reds are 2 pre-existing `AngularFlux|TimedFullField`
  union args at `:960`, NOT the seam); `multiplication_operator.py` went 1→0 (absorbed by relabeling,
  net exactly baseline). The bulk of the apparent diff is RELABELING (`OperatorSum[TimedFullField]`→
  `OperatorSum[FullField]` etc, count-neutral).
- **Regression subset `-O`: EXACTLY the 5 baseline reds, 2032 passed (+14 = C5), 0 non-baseline.**
  ⚠ The baseline-red set is **5, not 2**: the gate spec's `-k "not (sphere_1g_apply or sphere_2g_apply)"`
  filter MISSES 3 `test_bc_extraction_matvec.py::TestVacuumMatvecBitIdentity::*[*-SPH]` (a separate
  #250 stale-SPHERE-snapshot family, `[0/1/2-SPH]`, ~1e15 ULP — CONFIRMED identical on the baseline
  worktree) on top of the 2 #232 mu_y (`test_boundary_conditions::test_2d_mesh_resolution` +
  `test_native_matvec::TestTwoDCartesianRaises::test_two_d_cartesian_loss_action_returns_result`).
- **B snapshot gates**: slab/cylinder `TestT4b/c` PreT4 snapshots GREEN (8 pass — `(L+C)−C` VALUES
  byte-identical); the 2 `sphere_*_apply` snapshots = #250 baseline (confirmed on worktree; they DO
  execute `StreamingOperator.apply` = the matvec leaf, gate spec §E).
- **L1/MMS backstop**: `test_curvilinear_aniso_scattering_p1` + `test_kinf_homogeneous` +
  `test_si_convergence_rate` = 40 pass, 2 xfail (#200 sphere-4g-krylov). Converged limit UNMOVED.
- **Sphinx -W**: `build succeeded`.

## (f) UNTOUCHED (confirmed)
StreamingOperator's `(L+C)−C` arithmetic (`operator.py:431/473`) + the S/F `@singledispatchmethod`
INPUT-dispatch fibration (still dispatches on `TimedFullField`/`ScalarFlux`/`AngularFlux`; only the
TimedFullField-arm OUTPUT codomain changed).

## (g) SURPRISES (classified)
1. **The `full_field_space._rebuild` break (production, beyond the leaf set).** Classified:
   NECESSARY consequence of the codomain change — the G-adjoint metric helper hardcoded
   `_history=()` on the rebuilt composite, which is illegal once operator outputs are timeless
   `FullField`. Fixed via the polymorphic `_recombine` hook (the principled fix; NOT scope creep —
   it's the metric-application analogue of the dunder recombine). Resolved 7 g_adjoint test reds.
2. **41 B-tests broke by design (not 2-3).** Classified: expected test-surface tracking of the
   codomain change (every test asserting `isinstance(out, TimedFullField)` / reading `out.history_depth`).
   ALL are pure type/attr assertions, ZERO value regressions (confirmed — byte-identical converged
   state + slab/cyl snapshots reproduce frozen values). Re-pointed mechanically.
3. **Baseline reds = 5 not 2** (the gate-spec `-k` filter under-covers the bc_extraction `[*-SPH]`
   #250 family). Reconciled via the baseline worktree.

## NO deliverables owed beyond C5
- NO Sphinx stub / `:label:` (S8a mints NO verifiable claim — the timeless-codomain contract is a
  software invariant pinned by the C5 foundation test, not a theory equation; SAME posture as S4.5).
- NO `algebra-of-record` Branch-1 SymPy / L1-derivation manifest (pure type-surface carve, no new
  reference solver). NO new ERR (no bug found — the `_rebuild` break is a typing-tension fix forced
  by the codomain change, not a latent wrong-answer bug).
- NO archivist DISPATCH (the campaign defers the §5.x archivist narrative to a consolidated
  post-S6/S9 pass per the plan; S8a is plumbing/typing, mints no taxonomy).

## LESSON (propose to the algebra-of-record / vv skill ecosystem)
A codomain re-point of operator matvec leaves from a history-bearing comonad to its timeless base
is VALUE-NEUTRAL by construction ONLY because the iteration driver re-attaches the timed type via
the carrier's `__add__` recombine (NO `advance` is needed — verify `advance` is dead in production
first; it is, for steady-state SI/Krylov). The blast radius is THREE layers wider than the named
leaf set: (a) the shared `loss_action` producers (StreamingOp + InvertibleOp both delegate there),
(b) any metric/space helper that hardcodes `_history=()` on a rebuilt composite (route it through
the polymorphic `_recombine` hook instead — it handles both carriers), and (c) EVERY test asserting
`isinstance(out, TimedFullField)` / `out.history_depth` (dozens — they break by design, re-point to
the timeless contract). The resolvent `solve` is NOT a base arrow — it produces the driver's iterate
and STAYS timed (covariant, Liskov-safe under a `["FullField"]` class binding). The `@singledispatchmethod`
fibration is INPUT-dispatch (carrier TYPE) and is ORTHOGONAL to the OUTPUT codomain — re-pointing the
arm's return is NOT touching the fibration.
