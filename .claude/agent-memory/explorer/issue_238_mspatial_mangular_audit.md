---
name: issue-238-mspatial-mangular-audit
description: L20 retirement audit of StreamingOperator.M_spatial / M_angular_redist (the separately-applicable operator-leaf split) on refactor/sn-foundation-cleanup — verdict RETIRE-eligible (no live independent consumer; tests-only), but a thin documented KEEP is defensible. Full dependency surface + rewire map.
metadata:
  type: project
---

# #238 — do M_spatial / M_angular_redist (the operator-leaf split) earn their keep post-#206?

Audited on `refactor/sn-foundation-cleanup` (graph commit f6475e5, rebuilt 2026-06-18).
Line numbers current at that HEAD; re-derive via Nexus if they drift.

**Why:** #206 Phase C moved the 1-D matvec onto `_OneDimScanWalk` and unified
matvec≡sweep on the FUSED `loss_action`. The concern (L20): the
`(M_spatial, M_angular_redist)` separately-applicable-leaf split is now a
parallel path no production code walks.

**How to apply:** this is the documented dependency surface so #238 closes as a
decided audit (RETIRE-with-rewire OR documented-KEEP) without re-deriving.

## The shape (orpheus/sn/operator.py)

- `StreamingOperator.M_spatial` (cached_property, op.py:913) → `_MSpatialOperatorSum`
  (op.py:305), an `OperatorSum` of two `_SpatialSweepDirection` summands (op.py:213,
  signs +1/-1).
- `StreamingOperator.M_angular_redist` (cached_property, op.py:931) → `ZeroOperator`
  for Cartesian; `AngularRedistributionOperator` (op.py:469) for sphere/cylinder.
- The split is produced by `_MSpatialOperatorSum._compute_decomposition` (op.py:344),
  a ONE-LINE delegation to `_OneDimScanWalk.loss_action_decomposed`
  (loss_representation.py:2096), which calls `_apply_walk(..., emit_angular=True)`.
  NO cache (the long-claimed ψ-keyed cache was never implemented — corrected in the
  Phase C doc-fix; each consumer re-walks).

## Production hot path does NOT touch the split (the L20 crux)

`StreamingOperator.apply` (op.py:775, line 834) routes through
`self.loss_representation.loss_action(...)` = `_apply_walk(emit_angular=False)` — the
FUSED single-emission `(L+C)ψ` (skips the ~15% angular-residual store). Same for
`apply_transpose` (loss_action_transpose). The Krylov matvec and the SI rhs both ride
`loss_action`/the InvertibleOperator sweep, never `_compute_decomposition`.

Nexus confirms: `loss_action_decomposed` has exactly ONE caller —
`_compute_decomposition`. `_compute_decomposition`'s only callers are the THREE
standalone leaf applies (`_MSpatialOperatorSum.apply` op.py:444, `_SpatialSweepDirection.apply`
op.py:244 via transient orchestrator, `AngularRedistributionOperator.apply` op.py:551).
`callers()` on all three leaf `.apply` methods = ZERO (they are only reached by reading
the `.M_spatial`/`.M_angular_redist` property and calling `.apply`).

`grep M_spatial|M_angular_redist` in `orpheus/` outside operator.py = docstrings only
(loss_representation.py prose). NO production consumer reads either property.

## The separability capability has NO live consumer

The split exists to let a consumer apply `M_spatial·ψ` OR `M_angular_redist·ψ`
independently. Searched: the only would-be real consumers are OPEN issues, not code:
- #200 (block-inverse face preconditioner — the most likely real consumer of the
  angular/spatial split) — OPEN. Current GMRES preconditioner = `lambda q: q` explicit
  identity (solver.py:313); `_make_sweep_preconditioner` already RETIRED (solver.py:1453).
- #2 (DSA) — OPEN. All DSA/precond mentions in solver.py/scattering.py are
  forward-pointing docstrings.

So today nothing applies the leaves separately — everything uses fused `loss_action`.

## Test surface (all that pins the split)

Module `tests/sn/operators/test_streaming_operator.py` is `pytestmark =
pytest.mark.foundation`. NO `@pytest.mark.verifies(...)` / `@pytest.mark.catches(...)`
on ANY M_spatial/M_angular_redist test → deleting/rewiring them dangles NO equation
`tests` edge and NO ERR-NNN coverage.

Classes that exercise the leaves (the rewire/delete list):
- `TestT4bMSpatialStructure` (op.py:606) — type-level: M_spatial is `_MSpatialOperatorSum`
  with two `_SpatialSweepDirection` summands; slab M_angular_redist is ZeroOperator;
  curvilinear is AngularRedistributionOperator leaf; cached_property identity. PURE
  structure tests of the split — delete if the split goes.
- `TestT4bAlgebraDecompositionInvariantSlab` (op.py:715) — `L = M_spatial − σ_t·ψ`
  bit-exact (slab). Rewire onto `(L+C).apply` (=`_LC_matvec`) − σ_t·ψ.
- `TestT4cAlgebraDecompositionInvariantCurvilinear` (op.py:951) — `(L+C) == M_spat + M_ang`
  ULP-clean (sphere/cylinder) + M_ang writes no boundary + M_ang linearity. The
  decomposition-INVARIANT gate the issue names (TestT4c). If the split retires, the
  `m_ang ≠ 0` curvilinear redistribution still needs a pin — rewire to assert
  `loss_action_decomposed` (the surviving helper) sums to `loss_action`, OR keep just
  the helper-level test.
- `TestT5MaterializeInverseCache` (op.py:1181) — `M_spatial.materialize_inverse_cache()`.
  This method (op.py:395) has NO production caller either; it's an exposed-angle API on
  the orphan. Test compares to `CollisionCache.from_geometry` (the real factory) — rewire
  to call `from_geometry` directly, drop the M_spatial wrapper.
- `TestT4bMSpatialStandaloneApply` (op.py:1236) — per-direction summand sum == unified.
  Pure orphan-internal (tests `_SpatialSweepDirection.apply`) — delete with the split.

`tests/sn/_test_helpers.py:_LC_matvec` (line 322) does NOT depend on the split — it calls
`(L+C).apply` (fused). The graph "reference" edge from `_LC_matvec` → `_compute_decomposition`
is a DOCSTRING mention (line 330), not a call. Survives untouched.

`tests/sn/_fixtures/wave_t_t4/_capture_pre_t4_snapshots.py` — private test-support
(mesh builders `_sphere_mesh`/`_cylinder_mesh`/`_slab_for_snapshot_arm` + a pre-T.4
snapshot capture); consumed ONLY by test_streaming_operator.py. Builders survive; the
snapshot-decomposition docstring (line 6/254) is stale-if-split-goes.

## Internal-to-orphan (delete if RETIRE)

`_SpatialSweepDirection` (op.py:213) + its `.apply`; `_MSpatialOperatorSum` (op.py:305)
incl. `.apply`/`_compute_decomposition`/`materialize_inverse_cache`; `AngularRedistributionOperator`
(op.py:469) + `.apply`; the `M_spatial`/`M_angular_redist` cached_properties (op.py:913/931).
`loss_action_decomposed` (loss_representation.py:2096) becomes the surviving split-emit
helper IF a test needs it (the `emit_angular=True` arm of `_apply_walk` stays — it's the
single source). If NO test needs the decomposed split, `loss_action_decomposed` + the
`emit_angular=True` branch of `_apply_walk` ALSO retire.

Doc nodes pointing at the split (would dangle): `theory/operator_algebra` section
`wave-t-curvilinear-deep-dive` ("Curvilinear M_angular_redist — bespoke leaf") +
`theory/discrete_ordinates` `sn-mms-curvilinear-aniso-verification` (the redistribution
probe — that MMS is independent of the leaf API, stays; only the cross-ref to the leaf
needs an archivist touch).

## Verdict

RETIRE-eligible: no production caller, no live independent consumer, tests-only, no
dangling vv-markers. The discriminator (feedback_aggressive_retirement): same math is
available via the surviving `loss_action_decomposed` helper / the fused `loss_action` —
this is "same-math-available-via-the-helper", not "genuine independent-consumer need".

The honest counter-weight (Cardinal Rule 2 cuts both ways): the split is the cleanest
typed surface a future #200 block-inverse preconditioner (wants `M_angular_redist` alone)
or #2 DSA would consume. A thin documented-KEEP (collapse the standalone `.apply`
plumbing but keep the `M_spatial`/`M_angular_redist` PROPERTIES as named leaves over
`loss_action_decomposed`, pinned by TestT4c) is defensible. The call is the user's:
retire-now-rebuild-when-#200-lands vs keep-the-named-leaf-as-a-typed-anchor. Either way
the standalone `_SpatialSweepDirection.apply` slow path + `materialize_inverse_cache`
(both zero-consumer) should go.
