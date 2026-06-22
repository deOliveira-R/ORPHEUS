---
name: issue-238-mspatial-mangular-retirement-closeout
description: #238 L20 aggressive-retirement of the orphaned StreamingOperator.M_spatial / M_angular_redist separately-applicable operator-leaf split (zero production consumers post-#206); FULL retirement incl. loss_action_decomposed + the emit_angular arm; MMS-covered angular-redistribution V&V judgment; deletions + gate output + grep audit.
metadata:
  type: project
---

# #238 — retire M_spatial / M_angular_redist (the separately-applicable operator-leaf split)

Branch `refactor/sn-foundation-cleanup` 2026-06-18. NOT committed (main agent commits).
Pure deletion of UNUSED code — production byte-identical (proven, see Gate output).
This is an L20 aggressive-retirement of a solver-internal API: no Branch-1/L1/stub algebra-of-record
manifest (nothing new BUILT). The deliverable is the clean deletion + a clean Sphinx build + the
archivist dispatch for the deep-dive narrative rewrite.

## THE V&V JUDGMENT (the load-bearing decision) — FULL RETIREMENT

`loss_action_decomposed` + the `emit_angular=True` arm of `_apply_walk` had exactly ONE caller
(`_MSpatialOperatorSum._compute_decomposition`), which was deleted → they became dead. **Default
to FULL retirement** held: I retired `loss_action_decomposed` + the `emit_angular` arm too.

The discriminator: does the curvilinear MMS cover the angular-redistribution correctness?

- `tests/sn/verification/mms/test_mms_curvilinear.py` does NOT (it uses the ISOTROPIC ansatz
  `ψ_n(r)=A(r)/W` which NULLS the angular redistribution — Mode-7, declared in its docstring).
- `tests/sn/verification/mms/test_curvilinear_aniso_convergence.py` DOES: the ANISOTROPIC ansatz
  `ψ_n(r)=(A(r)+B(r)μ_n)/W` explicitly ACTIVATES `(1-μ²)B/r` (sphere) / `ξ²B/r` (cylinder) — its
  docstring states the pairing rationale verbatim: "If the isotropic companion passes and these
  fail, the bug is in the angular-redistribution path." These run the FULL curvilinear `(L+C)`
  end-to-end via `solve_sn_fixed_source` (NOT the leaf decomposition API), carry `@catches("ERR-026")`,
  and assert O(h²) + magnitude bands ("a re-floored wrong fixed point (ERR-026) would exceed
  5e-3 sphere / 5e-2 cylinder"). The angular term is O(1) in the curvilinear cell balance — a
  wrong/missing term re-floors the ladder or blows the band. **9 passed.**

The deleted decomposition tests (`TestT4c…`) asserted only an internal SPLIT structure
(`(L+C) == M_spat + M_ang`, an algebraic identity; `m_ang ≠ 0`, presence not correctness) that no
longer exists. No specific angular-redistribution bug class is BLIND to the MMS that the
decomposition tests would catch. The MMS is the surviving structurally-independent ground.
Keeping the machinery alive solely to feed a structural test = the same orphan smell one level down.

## DELETIONS

### orpheus/sn/operator.py (the orphan classes + properties)
- `_SpatialSweepDirection` (dataclass + `.apply`) — deleted.
- `_MSpatialOperatorSum` (`.apply`, `_compute_decomposition`, `materialize_inverse_cache`) — deleted.
- `AngularRedistributionOperator` (dataclass + `.apply`) — deleted.
- `StreamingOperator.M_spatial` / `StreamingOperator.M_angular_redist` cached_properties — deleted.
- `ZeroOperator` import (only used by `M_angular_redist`) — removed; `AngularRedistributionOperator`
  removed from `__all__`.
- The Wave-T explanatory comment block (the per-direction-split rationale) — replaced with a #238
  retirement note. Module docstring "geometry-agnostic kernel is `_MSpatialOperatorSum._compute_decomposition`"
  → repointed to `_OneDimScanWalk._apply_walk` + curvilinear-redistribution-is-in-sweep note.
- `cached_property` import KEPT (still used by `loss_representation`).

### orpheus/sn/loss_representation.py (FULL retirement)
- `loss_action_decomposed` (the `(M_spatial, M_angular_redist)` split) — deleted (was the sole
  caller of the `emit_angular=True` arm).
- `_apply_walk` signature: dropped `*, emit_angular: bool`; return type `(np.ndarray, BoundarySourceSink)`
  (was 3-tuple with `m_ang_cell`). Removed `out_ang_g_first` allocation, the 2 `if out_ang_g_first
  is not None:` blocks (curvilinear + degenerate-ordinate arms), `m_ang_cell` computation, docstring
  rewrite (single-emission). `loss_action`'s 3-tuple unpack → 2-tuple.
- Docstring repoints: `loss_action` (ScanMarch) "1-D → operator.M_spatial degeneration" →
  "1-D → `_OneDimScanWalk._apply_walk`"; `out_ang` comment trimmed.

### tests/sn/operators/test_streaming_operator.py (deleted 5 classes, kept 2)
- DELETED: `TestT4bMSpatialStructure`, `TestT4bAlgebraDecompositionInvariantSlab`,
  `TestT4cAlgebraDecompositionInvariantCurvilinear`, `TestT5MaterializeInverseCache`,
  `TestT4bMSpatialStandaloneApply`. Removed the unused `_LC_matvec` import. Rewrote the section
  comment to introduce the surviving snapshot gates + the #238 retirement rationale.
- KEPT (surviving production-path regression gates — use only `StreamingOperator.apply`):
  `TestT4bPreT4RegressionSnapshot` (slab) + `TestT4cPreT4RegressionSnapshotCurvilinear` (curvilinear).
  Their builders `_slab_for_snapshot_arm` / `_sigma_t_from_mat_map` survive (feed the kept gates).
- `TestT5MaterializeInverseCache` → DELETE (NOT rewire): after dropping the `M_spatial` wrapper the
  test would be `from_geometry == from_geometry` (vacuous). The audit's "rewire OR delete if
  redundant" → delete.

### A SURPRISE NOT IN THE BRIEF — a SECOND live consumer file
`tests/sn/operators/test_streaming_operator_decomposition.py` (the brief named only
`test_streaming_operator.py`). Its `TestResolutionADifferentFromPriorWrong::test_subtractive_L_differs_from_matvec_at_zero_sigma_t`
reached into `_MSpatialOperatorSum._compute_decomposition(state)` at σ_t=0 to build a "matvec at
σ_t=0" (bypassing `InvertibleOperator`'s σ>0 validation). REWIRED to
`default_for(sn_mesh).loss_action(sigma_zero, state)` — at σ_t=0, `(L+C)ψ = Lψ` (C ≡ σ·ψ = 0), so
the fused `loss_action` IS the σ_t=0 matvec `M(ψ;0)`; routing through the representation still
bypasses the σ>0 validation. The other two classes in that file (`TestResolutionADecomposition`,
`TestSubtractiveDefinition`) use only `L.apply`/`C.apply`/`_LC_matvec` — untouched. File passes.

### tests/sn/_test_helpers.py + _capture_pre_t4_snapshots.py (stale docstring repoints)
Docstring/comment mentions of `_MSpatialOperatorSum._compute_decomposition` as "the canonical
dual-emission body" were now WRONG (deleted) → repointed to `_OneDimScanWalk._apply_walk`. The
`_capture_LpC_apply` fixture function still captures `(L+C).apply` (valid); only its
decomposition-invariant docstring was updated. Builders survive.

## DOCS (Cardinal Rule 3 — clean build, zero dangling cross-refs)

Repointed ALL live `:class:`/`:attr:`/`:meth:` cross-refs to the deleted symbols (verified ZERO
remain via `grep -rnE ':(class|attr|meth|func|obj):` over deleted names in live .rst). Inline
double-backtick literals in retirement-narrative prose KEPT (deliberate historical mentions).

- `docs/theory/operator_algebra.rst`:
  - The `wave-t-deep-dive` bullets + the per-operator shape-catalogue table streaming rows →
    rewritten to "in-sweep fused matvec; #238 retired the typed-leaf split".
  - The three deep-dive subsections (`wave-t-streaming-deep-dive` / `wave-t-orchestrated-apply` /
    `wave-t-curvilinear-deep-dive`, ~273 lines) → rewritten cohesively. **PRESERVED** all four
    equation labels (`wdd-forward-recurrence`, `mm-half-grid-recurrence`,
    `wave-t-cell-balance-three-terms`, `wave-t-mspat-curvilinear-subtraction` — they carry
    `.. vv-status:` and feed the auto-generated `docs/verification/matrix.rst`) and all three
    section anchors (`wave-t-orchestrated-apply` is referenced by `docs/api/numerics.rst`). The
    physics (WDD + M-M coupling, the MA-Q1 condition) survives; only the typed-leaf framing went.
  - Added a `.. todo:: Archivist` flag in the `wave-t-orchestrated-apply` section for the deeper
    narrative rewrite (the five-future-consumer rationale + Design-B mechanism + curvilinear-
    subtraction smell, all retired).
  - Verification-approach gate list + implementation-map list + the B.5.2 retyped-sites table →
    repointed/condensed.
- `docs/api/numerics.rst`: OperatorSum primitive bullet + consumer-matrix row (2→1 load-bearing) +
  the bespoke-leaves list (3 deleted-class `:class:` entries → 1 retirement `.. note::`) +
  `M_spatial`/`M_angular_redist` properties mention → all repointed.
- `docs/theory/discrete_ordinates.rst`: 4 `:meth:`_MSpatialOperatorSum._compute_LpC[_transpose]``
  cross-refs (already-dead historical names that became dangling once the class was deleted) →
  repointed to `_OneDimScanWalk._apply_walk` / `loss_action_transpose`.

NOTE: I did NOT run a full Sphinx build (the main agent batches it per CLAUDE.md). I verified the
edits removed every dangling symbol cross-ref via grep and that all referenced labels/anchors
(`sn-mms-curvilinear-aniso-verification`, the 4 eq labels, the 3 section anchors) exist.

## GATE OUTPUT (host `.venv/bin/python -O -m pytest`)

- `tests/sn/operators tests/sn/spatial tests/sn/sweep/core tests/sn/solve` + the 2 curvilinear MMS
  files: **7 failed, 1096 passed, 6 skipped, 5 xfailed**. The 7 fails are EXACTLY the documented
  baseline reds (ZERO new): 5 stale sphere snapshots (#250 — `test_vacuum_bulk_bit_identical_1d[{0,1,2}-SPH]`
  + `test_sphere_{1g,2g}_apply_bit_identical`) + 2 mu_y (#232 — `Face 'ymin' requires genuine mu_y`).
- DD strict regression `tests/sn/regression -W "error::…DriftWarning"`: **13 passed, 13 within-tol
  DriftWarnings IDENTICAL to baseline ULP** — incl. the bit-identity proof
  `2d_2g_p1_aniso_dd_8x4_het_si: 6920 ULP / 9.81e-13 rel (within tol)`. This is THE production-
  byte-identical proof: deleting unused code moved NO snapshot, the strict gate passes at the SAME
  ULP as before. (Cross-checks #241 closeout's identical 6920-ULP baseline.)
- MMS angular-redistribution coverage (the V&V judgment ground): **9 passed**.

Why the 2 surviving sphere snapshot gates (which DO exercise my modified curvilinear `_apply_walk`)
are #250 not mine: the CYLINDER snapshots in the SAME class PASS — they also walk the modified
curvilinear arm, so if my `out_ang` removal had changed the fused output they'd break too. They
don't. #250 explicitly names these 2 sphere tests as stale-snapshot reds (the live apply is the
MORE-correct value, frozen snapshot stale since the W1 τ-unclamp `b2d8a6d`).

## FINAL GREP AUDIT

ZERO code references to the deleted symbols in `orpheus/` + `tests/`. Remaining hits are all
comments/docstrings: (a) my deliberate #238 retirement notes, (b) already-dead `_compute_LpC` /
`_compute_LpC_transpose` historical mentions (NOT in the delete list — dead names from the #206 /
O.4a.2 carves, deliberate historical comments). No module-level imports of the deleted symbols
anywhere (test collection safe).

## OWED
- archivist DISPATCH emitted (followup:false) for the deep-dive narrative rewrite (the
  `.. todo::` flag in `operator_algebra.rst` marks the spot).
- Full Sphinx `-W` build — main agent batches.
- Commit — main agent.

## LESSON (candidate for vv-principles / coding-elegance)
An orphan operator-leaf split retirement must audit for SECOND consumer files the brief's dependency
list missed: `_compute_decomposition` was reached not only by the named decomposition tests but by a
σ_t=0-bypass trick in a DIFFERENT test file (`test_streaming_operator_decomposition.py`) that used
the orphan purely as a validation-bypass vehicle. The fix is the same single-source rewire (route
through the surviving fused `loss_action` at σ_t=0), but the lesson is: `mcp__nexus__callers` on the
RETIRED internal method (`_compute_decomposition`) caught only the operator-side caller; a repo-wide
grep for the CLASS NAME (`_MSpatialOperatorSum`) is what surfaced the bypass-trick consumer. Run BOTH
(graph callers + text grep for the class/method name) before declaring a retirement's blast radius.
