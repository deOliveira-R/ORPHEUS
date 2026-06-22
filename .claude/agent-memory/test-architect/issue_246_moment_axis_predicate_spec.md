---
name: issue-246-moment-axis-predicate-spec
description: Verification spec (test-architect, PRE-IMPL) for #246 — the S4-hazard fix in _reframe (thread is_moment_valued from typed origin) + named Kind-A predicates (is_multi_moment / has_spatial_moment_axis). Centerpiece = the d=2 trailing-axis-length-4 collision negative control.
metadata:
  type: project
---

# #246 — moment-axis predicate verification spec (branch `refactor/sn-foundation-cleanup`)

PRE-IMPL test spec. The implementation is NOT written; this defines what
"correct" means for the carve. Authoritative dependency map =
`.claude/agent-memory/explorer/issue_246_moment_axis_predicate.md` (read it for
the site-by-site line numbers — verified live there). This note is the
verification side.

## What #246 actually changes (reframed scope — authoritative)

The literal "replace 4 value-probes with ONE typed predicate" is half-misframed.
Two KINDS:

- **Kind A (boolean "is the moment axis present?")** — routable through a typed
  predicate. Sites: `_reframe`'s array-length half (`sweep_graph.py:100`),
  `_moment_broadcast_sigma`'s rank compare (`loss_representation.py:515`), the
  five `if n_face_moments > 1:` / `if n == 1:` capture-reductions
  (`loss_representation.py:369,1068,1140,1319,1386`).
- **Kind B (the integer width `2^d` / `2^{d-1}`)** — `_n_face_moments`
  (`loss_representation.py:311`) + `_spatial_moment_tail`
  (`loss_representation.py:340`). These SIZE `np.zeros((...,*tail))`
  allocations. They are honest COUNTS, NOT predicates. They **STAY**. Do NOT
  test them as predicates; do NOT collapse them to a boolean.

The three deliverables of the implementation:

1. **The one genuine correctness win (Site 1 — the S4 hazard).** `_reframe`'s
   guard `if frame_signs is None or arr.shape[-1] != frame_signs.shape[0]:` keys
   the moment-axis question on a COINCIDENTAL trailing-axis length. The fix KEEPS
   the already-typed `frame_signs is None ⟺ DD/Step` half (verified live:
   `octant_moment_frame_signs` returns `None` iff `per_axis==1`,
   `_ubld.py:146-147`) and REPLACES the fragile `arr.shape[-1] != 2^d` half with
   a threaded `is_moment_valued: bool` the caller sets from the array's typed
   ORIGIN (is the walk producing a moment-valued buffer), not from a coincidental
   axis length.
2. **Named predicates for the Kind-A consumers.** `DiscretizationSchemeBase.is_multi_moment`
   (= `spatial_basis_per_axis > 1`, the scheme-level UNCONDITIONAL truth — the
   inner-walk sites have `self.mesh.scheme` in scope) and/or
   `BulkField.has_spatial_moment_axis` (= `spatial_moments_per_axis > 1`,
   field-level for boundary readability).
3. **KEEP `_n_face_moments` / `_spatial_moment_tail`** (the counts).

## The field-vs-scheme predicate distinction (LOAD-BEARING — gate design must respect it)

`SpatialMomentSpace` is carried on the typed Field ONLY because the
moment-aware OUTPUT wraps pass `spatial_moments=per_axis`; the default is `1`
everywhere and is NOT auto-read from `mesh.scheme.spatial_basis_per_axis`
(deliberate construct-general/select-narrow, `_bases.py:183-194`). Consequence:

- A hand-built `AngularFlux.from_mesh(values, mesh)` (no `spatial_moments=`) on
  an LD mesh **LACKS** the factor → `field.has_spatial_moment_axis` is `False`
  even though the mesh is LD.
- The SCHEME-level `spatial_basis_per_axis > 1` is the UNCONDITIONAL truth at
  the inner sites.

So the predicate gates MUST encode:
- **Scheme predicate `is_multi_moment` = UNCONDITIONAL** (LD scheme → True
  regardless of any field provenance).
- **Field predicate `has_spatial_moment_axis` = only trustworthy on
  moment-aware-PRODUCER output** (a moment-carrying field → True; a hand-built
  bare field on an LD mesh → False, BY DESIGN). A negative test that builds a
  bare LD field and asserts `has_spatial_moment_axis is False` is a REQUIRED
  documentation-of-intent gate (it pins the construct-general discipline — see
  Gate P3 below).

## Live baseline (probed @ HEAD on `refactor/sn-foundation-cleanup`, Host `.venv/bin/python -O`)

- `octant_moment_frame_signs(signs, per_axis)` → `None` iff `per_axis==1`
  (`_ubld.py:146-147`) — the `frame_signs is None` half is ALREADY a clean
  typed gate. Confirmed.
- `_CellSolve` / `_CellResidual` constructed at FOUR sites
  (`loss_representation.py:1046, 1123, 1302, 1372`), each passing
  `moment_frame_signs=self._moment_frame_signs(signs_eff)`. This is WHERE an
  `is_moment_valued` thread is set (alongside the existing frame-signs pass). The
  1-D scan `_reframe` sites (`loss_representation.py:2341,2353,2999,3042`)
  inherit the signature change; their callers have `frame_signs` (and the scheme)
  in scope.
- Gate files GREEN: `tests/sn/spatial/test_ld_slope_frame.py` (3) +
  `tests/sn/spatial/test_linear_discontinuous.py` (20) → 23 passed -O.
- `tests/sn/regression/` (DD bit-id, 13 cases) PASSES under
  `-W error::...DriftWarning` (44.9 s). NOTE: it emits 13 within-tolerance
  `DriftWarning`s (e.g. `2d_2g_LS4_dd_8x4_het_si: scalar_flux drifted 36293 ULP /
  6.28e-12 rel within tol 1.0e-09`). These are PRE-EXISTING and NOT escalated —
  the `-W error` flag fires only on ABOVE-tolerance drift. See "the bit-id gate"
  below for what #246 must hold.
- `tests/sn/verification/mms/test_mms_ld_2d.py -m "not slow"` GREEN (5p/2desel).

## Expected baseline reds to NOT depend on (per the dispatch brief)

- 5 stale sphere snapshots (#250) + 2 `mu_y` (#232) in `tests/sn/operators`.
- #212 hang at `tests/sn/solve/test_keff_slab::test_heterogeneous_absolute_keff`
  — NEVER run all of `tests/sn`. The gates below NEVER touch `test_keff_slab`.

---

## ⭐ GATE 1 — THE S4-HAZARD NEGATIVE CONTROL (the centerpiece)

**File:** `tests/sn/sweep/core/test_reframe_moment_intent.py` (NEW, unit-level).
**Level:** L0 (term-level — `_reframe` is a single involution applied to one
array; this is the per-term sign/intent verification). `@pytest.mark.l0`.
**`-O`-safe:** `np.testing.assert_array_equal` / `pytest.fail` (function calls
that fire under `-O`). NO bare `assert` (Mode 8).

### Why this gate exists (the failure it catches)

`_reframe`'s OLD probe `arr.shape[-1] == frame_signs.shape[0] == 2^d` collides
when a genuine NON-moment array happens to have a trailing axis of length `2^d`.
On a d=2 mesh `2^d == 4`; a flat source / DD array whose trailing axis is
coincidentally 4 (an anti-diagonal slot count of 4, a flat external with 4 in
the last position) would mis-fire the moment-frame sign-flip involution
(`octant_moment_frame_signs`) — corrupting a non-moment array with a spurious
slope re-sign. It is LATENT today (no flat-source d≥2 LD path reaches `_reframe`
with such an array), so the fix is prophylactic. Per vv-principles §1/§6: the
fix makes the predicate key on the array's INTENT (typed origin), not on a
coincidental axis length — this gate is the negative control that proves it.

### Test rows (all UNIT-level — call `_reframe` directly; production cannot reach
the collision today, so a production-path test is impossible — this is the
correct level per the dispatch brief)

| Test | What it asserts | Regime that exposes the failure |
| ---- | --------------- | ------------------------------- |
| `test_reframe_non_moment_4trailing_passthrough_when_not_moment_valued` | Construct `frame_signs` for a per_axis=2 d=2 scheme (a `(4,)` involution with at least one `-1` entry — use `octant_moment_frame_signs((1,-1), 2)` so it is NOT all-ones), and a NON-moment array `arr` whose trailing axis is length 4 but is NOT a moment vector (e.g. `np.arange(...).reshape(ng, 4)` standing in for a flat/DD buffer). With the NEW signature `_reframe(arr, frame_signs, is_moment_valued=False)` → returns `arr` UNCHANGED (`assert_array_equal(out, arr)`). | The d=2 `2^d==4` collision: trailing axis length 4 matches `frame_signs.shape[0]==4`. The OLD probe `arr.shape[-1] != frame_signs.shape[0]` is FALSE (they're equal) → OLD would multiply by `frame_signs` and flip signs. The NEW `is_moment_valued=False` thread keeps it untouched. |
| `test_reframe_real_moment_applies_involution_when_moment_valued` | Same `frame_signs` (with a `-1`), a REAL moment array `arr` (trailing axis = the `2^d==4` moment vector). `_reframe(arr, frame_signs, is_moment_valued=True)` → returns `arr * frame_signs` (`assert_array_equal(out, arr * frame_signs)`). | The positive half: a genuine moment array at a moment closure DOES get the involution. Proves the predicate is not a blanket "never reframe". |
| `test_reframe_dd_step_passthrough_frame_signs_none` | `frame_signs is None` (DD/Step) with `is_moment_valued` either True or False → always returns `arr` unchanged. | The backward-compat invariant: DD/Step short-circuits on `frame_signs is None` REGARDLESS of the new thread. Both values of `is_moment_valued` must pass through. |
| `test_reframe_old_probe_would_misfire_on_collision` (DOCUMENTING the hazard) | Replicate the OLD guard logic inline (`out_old = arr if (frame_signs is None or arr.shape[-1] != frame_signs.shape[0]) else arr * frame_signs`) on the Gate-1 row-1 inputs and assert `not np.array_equal(out_old, arr)` (the OLD code WOULD have flipped) while the NEW `_reframe(arr, frame_signs, is_moment_valued=False)` returns `arr`. | Explicit demonstration that the OLD probe mis-fires and the NEW thread fixes it — the negative-control PAIR (vv-principles anti-pattern #11: a fix gate needs the broken-instance demonstration AND the fixed-instance pass). |

### Design notes / traps

- **`frame_signs` MUST carry a `-1`.** If you pick `octant_moment_frame_signs`
  for an all-forward octant (e.g. `(1,1)`) it returns all-ones → `arr * signs ==
  arr` and the collision is INVISIBLE (Mode-10: term activated but unconstrained).
  Use a BACKWARD octant like `(1,-1)` or `(-1,-1)` so the involution actually
  flips at least one slot. Verify in the test that `(-1 in frame_signs)` before
  asserting (a premise guard, `pytest.fail` if not — mirrors the pure-z premise
  guard idiom in `test_mms_ld_2d.py:275`).
- **Implementation-shape dependency.** This gate assumes the implementation
  threads a boolean `is_moment_valued` param into `_reframe` (the explorer's
  recommended option-1/elegant-shortcut). IF the implementer instead keeps
  `frame_signs` as the SOLE gate and threads `is_moment_valued` by passing
  `frame_signs=None` for non-moment arrays (the redundant-with-`frame_signs is
  not None` route the explorer flags at `_bases`/`sweep_graph:838-926`), then
  the signature is `_reframe(arr, frame_signs)` UNCHANGED and the "non-moment"
  case is expressed as `frame_signs=None`. In THAT case rows 1/2/4 collapse into
  the existing `frame_signs is None` semantics and the S4 hazard is fixed at the
  CALLER (the caller passes `None` for `Q_cells` in the matvec-zero/flat case).
  **The gate must follow the chosen signature** — coordinate with the
  method-implementer at impl time. The SKELETON below is written for the
  explicit-`is_moment_valued`-param shape (the recommended one); a one-line edit
  retargets it if the signature differs. EITHER WAY the assertions (pass-through
  on non-moment, involution on moment, None-short-circuit on DD) are identical;
  only the call signature changes.

---

## GATE 2 — THE BIT-IDENTITY GATE (existing tests that MUST stay byte-identical)

The change is bit-identical to ALL existing production paths:
- **DD** has `frame_signs is None` → the guard short-circuits BEFORE the
  `arr.shape[-1]` half is ever evaluated → DD never reaches the changed code.
- **LD moment arrays** set `is_moment_valued=True`, matching the old `== 2^d`
  true-branch exactly (the involution still applies to real moment arrays).
- **The Kind-A predicate routing** (`is_multi_moment` / `> 1`) computes the
  IDENTICAL boolean (`spatial_basis_per_axis > 1`) the inline `n_face_moments >
  1` already computed — naming-only, no numerical change (Pattern-6).

So NO snapshot should regenerate and NO converged value should move. The gate is
to RUN these existing suites unchanged after the carve and confirm green +
NO NEW above-tolerance `DriftWarning`.

### Bit-identity suites to run (the existing pins — do NOT modify, just RE-RUN)

| Suite / file | Covers | Invocation | Expected |
| ------------ | ------ | ---------- | -------- |
| `tests/sn/regression/` (`test_dd_regression.py`, 13 cases incl. `2d_1g_LS4_dd_15x15`, `2d_2g_LS4_dd_8x4_het_si`, `slab_fixed_source_dd_n20`, both `*_p1_aniso_*`) | DD matvec/sweep bit-identity in d=1 AND d=2 (the d=2 cases route DD through `_reframe` with `frame_signs=None` → the short-circuit half) + the flat-source slab case | `-O -W "error::tests.sn.regression._regression_assert.DriftWarning"` | 13 passed. The 13 pre-existing WITHIN-tolerance DriftWarnings stay (they are NOT escalated). **NO snapshot regenerates** (the change cannot touch DD). If any case escalates to a FAILURE the carve broke DD bit-id. |
| `tests/sn/spatial/test_ld_slope_frame.py` | ERR-061 — `_reframe` exercised through `solve_sn_fixed_source` (LD slab, both sweep directions); the `catches("ERR-061")` regression for the moment-frame involution | `-O` | 3 passed (no change — LD moment arrays still get the involution via `is_moment_valued=True`). |
| `tests/sn/spatial/test_linear_discontinuous.py` (`TestLDKernel::test_residual_zero_at_solved_cell_avg_2d`, `TestLDLinearExactness`, `TestLDRoundTrip`) | The d=2 LD `cell_kernel_batch`/`residual_kernel_batch` round-trip — the kernel-level matvec twin that consumes the reframed moment arrays in BOTH sweep directions | `-O` | 20 passed (LD moment path unchanged). |
| `tests/sn/verification/mms/test_mms_ld_2d.py` (foundation rows: `test_ld_2d_two_paths_ffw_equals_mfw`, `test_dd_and_ld_2d_converge_to_different_values`, `test_ld_2d_krylov_equals_si_pure_z_quadrature` [ERR-062]) | The 2-D wavefront `_CellSolve`/`_CellResidual` walk (FFW≡MFW), DD≠LD discrimination, the `_moment_broadcast_sigma` pure-z arm (Site 3) | `-O -m "not slow"` | 5 passed / 2 deselected (the predicate routing on Site 3's rank-compare-replacement does not change the broadcast result). |
| `tests/sn/verification/mms/test_mms_ld_2d.py` (slow: `test_ld_2d_stress_krylov_equals_si` [D5b.4], `test_ld_2d_stress_converges_second_order`) | The 2-D LD Krylov-vs-SI matvec twin on the slope-active stress habitat (the `loss_action`/`residual_kernel_batch` path that threads `is_moment_valued` through `_CellResidual`) | `-O` (slow; ~2 min) | green (the matvec twin still reaches the SAME fixed point — `is_moment_valued=True` for the moment probe). |
| `tests/sn/sweep/core/test_phase_c_gates.py` | The 1-D scan `_reframe` sites (`loss_representation.py:2341,2353,2999,3042`) inherit the signature change; this is the 1-D Phase-C matvec/sweep gate | `-O` | green (the 1-D scan moment arrays are moment-valued → `is_moment_valued=True`; DD 1-D scan → `frame_signs=None`). |

### Why DD d=2 is the load-bearing bit-id case (not just d=1)

The S4 collision is a d=2 phenomenon (`2^d==4`). The d=2 DD regression cases
(`2d_1g_LS4_dd_15x15`, `2d_2g_LS4_dd_8x4_het_si`) are the ONLY existing pins that
exercise `_reframe` on a d=2 mesh with `frame_signs=None`. They prove the
short-circuit half still fires for DD in d=2 (the regime where the OLD probe's
`arr.shape[-1] != 4` half would START to matter if the short-circuit ever broke).
KEEP both in the bit-id run.

---

## GATE 3 — THE NAMED-PREDICATE POSITIVE + NEGATIVE (vv anti-pattern #11)

A contract-style predicate needs a positive (real multi-moment → True) AND a
negative (real single-moment → False) test. TWO predicates → two pairs.

**File:** `tests/sn/spatial/test_moment_axis_predicates.py` (NEW).
**Level:** `@pytest.mark.foundation` (software invariant — the predicate is a
typed query over a data structure, NOT an equation `:label:`; no `verifies`).
**`-O`-safe:** `np.testing.assert_*` / `pytest.fail` (the predicates return
`bool`; assert via `if not p: pytest.fail(...)`, never bare `assert`).

### P1/P2 — `DiscretizationSchemeBase.is_multi_moment` (scheme-level, UNCONDITIONAL)

| Test | Asserts |
| ---- | ------- |
| `test_is_multi_moment_true_for_linear_discontinuous` (P1, positive) | `LinearDiscontinuous().is_multi_moment is True` (and `== (spatial_basis_per_axis > 1)`). |
| `test_is_multi_moment_false_for_diamond_difference` (P2, negative) | `DiamondDifference().is_multi_moment is False`. |
| `test_is_multi_moment_false_for_step` (P2', negative) | The Step scheme (if instantiable today — `spatial_basis_per_axis == 1`) → False. IF Step is still a docstring stub (per the D5 note it does not exist as a live class), SKIP with `reason=` pointing at the unlocking class. Verify at impl time. |

**Mutation check (REQUIRED, run once at impl, documented in the test docstring):**
re-define `is_multi_moment` to return the constant `True` and confirm P2 reddens;
return constant `False` and confirm P1 reddens. If a mutation does NOT redden,
the predicate body is not the `spatial_basis_per_axis > 1` it claims (Mode-10).

### P3/P4 — `BulkField.has_spatial_moment_axis` (field-level, PROVENANCE-dependent)

The field predicate is the construct-general/select-narrow one — the negative
test is NOT merely "DD field → False"; it is the LOAD-BEARING documentation that
a moment-LACKING field on an LD mesh is also False (the default-`1` discipline).

| Test | Asserts | Why |
| ---- | ------- | --- |
| `test_has_spatial_moment_axis_true_on_moment_aware_producer_output` (P3, positive) | Take a field produced THROUGH a moment-aware wrap (the explorer names `operator.py:1622`, `scattering.py:786/1163`, the LD sweep OUTPUT — the cleanest reachable producer is the LD `result.angular_flux` from a `solve_sn_fixed_source(..., scheme=LinearDiscontinuous())` call, which carries `spatial_moments=per_axis>1` via the output wrap). `result.angular_flux.bulk.has_spatial_moment_axis is True` (and `== spatial_moments_per_axis > 1`). | The field-keyed predicate is correct on producer output (its only trustworthy regime). |
| `test_has_spatial_moment_axis_false_on_dd_field` (P4, negative) | A DD `solve_sn_fixed_source` output field → `has_spatial_moment_axis is False`. | The plain negative. |
| `test_has_spatial_moment_axis_false_on_hand_built_ld_field` (P4', the CONSTRUCT-GENERAL pin) | Build `AngularFlux.from_mesh(values, sn_mesh)` (NO `spatial_moments=`) on an LD mesh → `has_spatial_moment_axis is False`, DESPITE the mesh being LD. Docstring states this is INTENTIONAL (the factor is producer-set, not mesh-derived; `_bases.py:183-194`). | This pins the deliberate construct-general discipline. WITHOUT it, a future "helpful" change that auto-reads the scheme into the field default would pass silently and break LD byte-id. The test documents WHY the field predicate is provenance-dependent — and is the reason the INNER walk uses the SCHEME predicate, not the field one. |

**Note for the implementer:** because the field predicate is provenance-dependent,
the production INNER-walk routing of the Kind-A `> 1` probes MUST use the
SCHEME-level `is_multi_moment` (or the existing `self._n_face_moments > 1` /
`self._spatial_moment_tail != ()`), NOT a field query. The field predicate is for
method-boundary readability only. Gate P4' is what makes this constraint visible.

---

## GATE 4 — COVERAGE MAP (which sites the gates exercise; blind spots)

### Kind-A consumers and their catcher

| Site | File:line | Kind | Catching gate | Notes |
| ---- | --------- | ---- | ------------- | ----- |
| Site 1 `_reframe` array-length half | `sweep_graph.py:100` | A (the S4 hazard) | **Gate 1** (unit) + Gate 2 DD-d2/LD-2D/slope-frame (production exercise) | The centerpiece. Gate 1 is the ONLY direct unit test of the collision; production cannot reach it today (latent). |
| Site 1 1-D scan `_reframe` ×4 | `loss_representation.py:2341,2353,2999,3042` | A | Gate 2 `test_phase_c_gates.py` (1-D scan) + `test_ld_slope_frame.py` (LD slab through the scan) | Inherits the signature change. The 1-D scan moment arrays ARE moment-valued (`is_moment_valued=True`); DD 1-D → `frame_signs=None`. |
| Site 3 `_moment_broadcast_sigma` | `loss_representation.py:515` | A (rank compare → predicate) | Gate 2 `test_ld_2d_krylov_equals_si_pure_z_quadrature` (ERR-062 — the pure-z arm, the EXACT habitat of this helper, callers `:698`/`:789`) | The pure-z arm is the only path that calls `_moment_broadcast_sigma`. ERR-062's gate already pins sweep≡matvec on it. If #246 threads `is_moment_valued` here, this gate re-confirms the broadcast is unchanged. |
| Site 2 `if n_face_moments > 1:` ×4 + `if n == 1:` | `loss_representation.py:1068,1140,1319,1386,369` | A (naming swap) | Gate 2 `test_mms_ld_2d.py` FFW≡MFW + Krylov≡SI (these walks hit all four `> 1` capture-reductions in d=2 LD) | Pattern-6 naming; `spatial_basis_per_axis > 1` ≡ the old `n_face_moments > 1` by construction. The FFW≡MFW two-paths gate exercises every capture-reduction. |

### Kind-B (counts — STAY; NOT predicate-tested)

| Site | File:line | Catcher (unchanged) |
| ---- | --------- | ------------------- |
| `_n_face_moments` | `loss_representation.py:311` | Gate 2 (the `np.zeros((...,*tail))` allocations it feeds are pinned by the FFW≡MFW shape agreement + the LD kernel round-trip). |
| `_spatial_moment_tail` | `loss_representation.py:340` | Same. |

### Blind spots (documented, not closed by #246)

1. **The production S4 collision is unreachable today** (no flat-source d≥2 LD
   path), so Gate 1 is necessarily UNIT-level. There is NO production-path
   regression for the hazard until #240 D5b-S4 (the non-vanishing boundary
   trace) or a flat-source d≥2 LD path lands. When that path lands, ADD a
   production-path row to Gate 1 (a real `solve_sn_fixed_source` flat-source d=2
   LD whose intermediate array has trailing axis 4). Flag this in the impl's
   honest-scope note.
2. **The 3-D `2^d==8` collision** is even more latent (no 3-D LD production).
   Gate 1 should parametrize the d-dimension (`per_axis=2`, `d in {2,3}` →
   `frame_signs` of length 4 and 8) so the predicate-keys-on-intent claim is
   verified for both — cheap, and the d=3 row is the only coverage 3-D LD will
   get for this until a 3-D consumer arrives.
3. **`_moment_broadcast_sigma`'s rank compare was ALREADY S4-safe** (rank, not
   coincidental length — explorer Site 3). If #246 threads `is_moment_valued`
   into it for consistency, the change is naming-only; the ERR-062 gate is
   sufficient. Do NOT manufacture a new hazard test for Site 3 — there is no
   collision there.

---

## Sequence + invocations (copy-paste, Host `-O`)

```
# Baseline (run BEFORE the carve to confirm green floor):
.venv/bin/python -O -m pytest tests/sn/spatial/test_ld_slope_frame.py \
  tests/sn/spatial/test_linear_discontinuous.py tests/sn/sweep/core/test_phase_c_gates.py -q

# After the carve — Gate 1 (new):
.venv/bin/python -O -m pytest tests/sn/sweep/core/test_reframe_moment_intent.py -q

# After the carve — Gate 3 (new predicates):
.venv/bin/python -O -m pytest tests/sn/spatial/test_moment_axis_predicates.py -q

# After the carve — Gate 2 bit-id (DD strict + LD paths):
.venv/bin/python -O -m pytest tests/sn/regression/ \
  -W "error::tests.sn.regression._regression_assert.DriftWarning" -q   # expect 13 passed, NO new above-tol escalation
.venv/bin/python -O -m pytest tests/sn/spatial/test_ld_slope_frame.py \
  tests/sn/spatial/test_linear_discontinuous.py \
  tests/sn/verification/mms/test_mms_ld_2d.py -m "not slow" \
  tests/sn/sweep/core/test_phase_c_gates.py -q
# (slow LD-2D Krylov≡SI matvec twin — run once before merge:)
.venv/bin/python -O -m pytest tests/sn/verification/mms/test_mms_ld_2d.py -q
```

NEVER run all of `tests/sn` (#212 hang at `test_keff_slab`).

---

## Self-improvement note (fire BEFORE delivery)

No NEW failure mode is introduced by #246's verification plan. The S4 hazard is
a textbook **Mode-6 (convention drift)** instance at the test-DESIGN level — the
OLD probe keyed a moment-axis predicate on a coincidental axis length (definition
site `octant_moment_frame_signs` vs usage site `arr.shape[-1]` disagree on what
"this axis is the moment axis" means). It is ALSO a clean **vv-principles
anti-pattern #11** positive/negative-pair instance (Gate 1's misfire-pair + Gate
3's predicate pairs). Both are already in the skill's tables; no append needed.
The implementation should log the prophylactic S4 fix under an ERR-NNN ONLY IF
a real collision is ever observed in production (it is latent today — per the
"Log every caught bug" directive, an ERR entry is for a CAUGHT bug, not a
prophylactic). Until then it is a documented latent-hazard fix, no ERR.

Links: [[issue_206_phase_c_verification]] (the 1-D scan `_reframe` / Phase-C
gate + DriftWarning escalation path), [[d5b_s3_inc_c_moment_iterate_verification]]
(ERR-061 the moment-frame involution origin; `_reframe` was minted in D5b-S3),
[[d5_nd_polymorphism_verification]] (the `2^d` moment-tensor + `is_multi_moment`
scheme-trait family), [[feedback_vv_tagging]] (foundation NO `verifies()`;
predicate gates are software invariants), [[feedback_regression_tolerance_design]]
(the within-tol DriftWarning vs above-tol escalation distinction).
