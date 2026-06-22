---
name: issue-251-legB-boundary-trace-closeout
description: "#251 Leg B — carry the 2^{d-1} transverse face-slope through the boundary trace for 2-D Cartesian LD (the BOUNDARY twin of #247 Leg A's bulk slope source). Widened geometry.boundary_face_layout (the ONE-site storage lever: append face_moment_tail to each face slot), loss_representation._inflow_to_moments + _octant_face_cochain (rank-discriminated moment pass-through), dropped 4 outflow capture-collapse sites, boundary_source_sink.prescribed_inflow (scalar-or-moment slot assign). 2 xfail-strict gates flipped to live + consumption gate re-targeted onto the public API (Mode-11 closed). DD/Step + scalar-inflow byte-identical. Branch refactor/sn-foundation-cleanup. NOT committed."
metadata:
  type: project
---

# #251 Leg B — boundary transverse-face-slope (the LM-1989 trap BOUNDARY half)

Issue #251 (split from #247), L1 verification. Branch
`refactor/sn-foundation-cleanup` (main checkout, NOT a worktree). Host
`.venv/bin/python -O`. **NOT committed** (main agent commits after elegance + qa
review). The BOUNDARY twin of [[issue-247-legA-moment-source-closeout]] (the BULK
slope source, COMMITTED `d9396a2`).

Implements [[issue-251-legB-boundary-gate-spec]] (test-architect) +
[[issue-251-trace-space-widening]] (explorer). Closes the BOUNDARY half of the
LM-1989 slope-row trap: the manufactured prescribed inflow varies transversely
ALONG each face, but `_inflow_to_moments` zeroed the `2^{d-1}` transverse
face-slope. The transverse face-slope is now CARRIED + sign-constrained.

## THE PRODUCTION CHANGE (4 files, the explorer's minimal site list — held EXACTLY)

The layout-widening path (the explorer's §0 choice): a persistent moment trace
serves BOTH the inflow read AND the outflow store with ONE lever, and
`trace_space.py` needs ZERO changes (its metric + `omega_dot_n` are per-ordinate
and broadcast over the trailing moment axis — "moment-ready by accident").

1. **STORAGE LEVER — `orpheus/sn/geometry.py:boundary_face_layout`** (the ONE
   site that makes the trace moment-valued). Append the scheme's face moment tail
   to each slot shape: `(N, ng, *face_shape(label), *face_moment_tail(per_axis**(ndim-1)))`.
   `self.scheme.spatial_basis_per_axis` is in scope (scheme set at geometry.py:271
   BEFORE the trace at :517). DD/Step → `per_axis^(d-1) == 1` → `face_moment_tail(1)
   == ()` → byte-identical. 1-D slab face is a point (`face_shape == ()` →
   `per_axis^0 == 1` → no tail). VERIFIED live: LD face slots `(24,2,6,2)`/`(24,2,8,2)`
   (the trailing 2 = the `2^{d-1}` transverse-moment axis); DD slots `(24,2,6)`/`(24,2,8)`
   (no axis).

2. **INFLOW lift — `orpheus/sn/loss_representation.py:_inflow_to_moments`** (357).
   The boundary twin of Leg A's `_lift_external_source_to_moments`. Rank-DISCRIMINATES
   against the flat face rank `ndim + 1` (a scalar face `(N_oct,ng,*transverse)` has
   rank `2 + (d-1) = d+1`; transverse carries `d-1` axes) via the Leg-A single-source
   primitive `is_moment_valued_by_flat_rank` (NO third rank spelling — Leg-A NIT-1
   honored). Arms: `n==1` → identity; scalar → widen (slot 0 ← scalar, slopes ZERO —
   the scalar default, the Leg-B asymmetry); moment-resolved → PASS THROUGH (validate
   trailing width == `2^{d-1}`, ValueError naming `2^(d-1)` otherwise — Pattern 4).
   Added `is_moment_valued_by_flat_rank` to the existing `moment_layout` import (131).

3. **FULL-FIELD ORACLE seed — `loss_representation.py:_octant_face_cochain`** (the
   explorer's `_allocate_full_cochain`, 1214). The FFW oracle does NOT call
   `_inflow_to_moments` (it widens inline at the IN-edge seed). Same rank discriminator
   added: `n==1` → original; moment-resolved inflow → seed ALL `2^{d-1}` moments
   (`buf[tuple(in_edge)] = inflow[a]`); scalar → seed slot 0 only.

4. **OUTFLOW store — 4 capture-collapse sites DROPPED** (`loss_representation.py`):
   MFW solve `_sweep_interior` (was 1068-69), MFW apply `_loss_action_interior`
   (1140-41), FFW solve (1319-20), FFW apply (1386-87). Each `if n_face_moments > 1:
   capture = tuple(c[..., AVERAGE_MOMENT] for c in capture)` removed → the capture
   retains the `2^{d-1}` axis. The sheds (708 SOLVE → `boundary_flux.face_view`; 801
   MATVEC → `streamed` = `zeros_like(face_view)`, auto-widens) + the B-residual emit
   (819-823, ordinate-indexed) land the moments AUTOMATICALLY into the now-moment-shaped
   slot. `n_face_moments` is still used (passed to `walk_*`, `_octant_face_cochain`) →
   no unused-var.

5. **PRODUCER — `orpheus/transport/source_sinks/boundary_source_sink.py:prescribed_inflow`**.
   TWO edits, BOTH needed (the architecture bit harder than the explorer's map — see
   below):
   - `view[inflow, :] = arr[inflow, :]` → discriminated assignment: `arr.shape ==
     view.shape` → `view[inflow] = arr[inflow]` (the `[inflow]` spans ALL trailing
     axes; the `, :` assumed exactly one — DD's full slot AND a moment-resolved full
     slot both take this byte-identical arm); `arr.shape == view.shape[:-1]` (a
     SCALAR onto a moment slot) → `view[inflow, ..., AVERAGE_MOMENT] = arr[inflow]`
     (seed the average, slopes zero — the scalar default); else ValueError.
   - Imported `AVERAGE_MOMENT` from `numerics.moment_layout` (downward dep; numerics <
     transport). Docstring `(N, ng)` prose → moment-shape at the 3 prose sites.

## ⚠ WHERE THE ARCHITECTURE BIT HARDER THAN THE EXPLORER'S MAP STATED

**The explorer's `prescribed_inflow` §5 was INCOMPLETE.** The map said the line-274
edit (`, :` → `[inflow]`) was the only producer change and the shape-check would
"auto-validate the wider shape (it compares against the slot shape, which grew)".
That assumed the CALLER supplies a moment-resolved array. But the EXISTING scalar
callers — the 2-D LD MMS `case.prescribed_inflow` (mms/sn.py:1542) and the 1-D
prescribed-inflow MMS — supply a SCALAR `(N, ng, *transverse)` array. After the
trace widened, the LD slot is `(N, ng, *transverse, 2)`, so the unconditional
`arr.shape != view.shape` reject FIRED on every existing scalar LD caller (8 LD-2D
gates went red: `test_ld_2d_stress_krylov_equals_si`, both two-paths, all 4
external-slope #247 gates, scattering M4, the boundary consumption gate). The fix
is the SCALAR-onto-a-moment-slot arm (seed slot 0) — a SECOND producer edit the
explorer's "1 real edit" estimate missed. This is the SAME class of under-scope as
the Leg-A `_inflow_to_moments` field-space layer: a producer that lives ABOVE a
rigid scalar contract needs a typed-union relaxation, not just an indexing fix.
**Lesson: when a carve widens a trace slot, audit the EXISTING SCALAR producers
(the MMS cases), not only the new moment-resolved one — the scalar callers feed the
SAME widened slot and the producer must accept both ranks.**

## THE SURROGATE → PUBLIC-API RE-TARGET (vv Mode-11 closed)

The test-architect's de-risk consumption gate (`test_ld_2d_boundary_slope_sign_mutation_reddens`)
was surrogate-driven: a monkeypatched `_inflow_to_moments` that seeded slot-1 from
a test-supplied slope dict (pinning the CONSUMER but Mode-11-blind to the PRODUCER
trace stamp). Once production landed the moment-resolved trace end-to-end, the
surrogate CONFLICTED (the trace now hands `_inflow_to_moments` a moment-valued
inflow, so the monkeypatch appended a spurious 2nd `(...,2)` axis → frontier-seed
crash). Per the spec §4/§10 mandate, I RE-TARGETED the gate onto the public API:
`_solve_with_boundary_slope` now builds a moment-resolved
`BoundarySourceSink.prescribed_inflow` (slot 0 = centre, slot 1 = ±projected
slope), bundles it with the manufactured bulk source, and solves via the public
`solve_sn_fixed_source` — exercising the producer stamp AND the cochain consumer
END-TO-END (Mode-11 blindness closed; the dedicated stamp catcher is the threading
foundation gate). The monkeypatch is fully removed.

## THE TESTS (`tests/sn/verification/mms/test_mms_ld_2d.py`, the #251 block)

The 2 xfail-strict gates FLIPPED to live (`xfail` markers removed; both are now
plain `@pytest.mark.foundation`):
- `test_ld_2d_boundary_slope_threaded_through_inflow_to_moments` (foundation) — the
  STRUCTURAL teeth (the stamp/production-change proof): `_inflow_to_moments` threads
  the projected moment-resolved face's slot-1 through at machine precision
  (`assert_array_equal`), and recognises the moment-resolved input (no spurious axis).
- `test_ld_2d_boundary_trace_rejects_wrong_transverse_width` (foundation) — D6 the
  negative pin: `_inflow_to_moments` rejects a `(...,3)` trailing width with a
  ValueError naming `2^(d-1)` (Pattern 4 — the moment relaxation does not swallow
  shape bugs).

GREEN-now gates (unchanged or re-targeted):
- `test_face_transverse_legendre_projection_matches_hand_polynomial` (foundation) +
  `test_face_projection_slot0_is_transverse_cell_average` (foundation) — D3, the L11
  face projector (leggauss only). Untouched.
- `test_ld_2d_boundary_scalar_inflow_no_op_negative_control` (foundation) — the
  scalar widening keeps slot-1 EXACTLY zero (the Leg-B asymmetry). Untouched.
- `test_ld_2d_boundary_slope_sign_mutation_reddens` (l1, verifies ld-cartesian-2d) —
  D2 the PRIMARY sign-catcher, RE-TARGETED onto the public API (monkeypatch dropped).

"improves-on-flat" is DROPPED (probed FALSE — the sharper Mode-10; the boundary
slope is sub-floor). The positive verification is the structural lift-threading +
the consumption mutation. No `@catches` (no caught bug; next free ERR-063).

## VERIFICATION (verbatim `-O` pytest summary lines)

```
# Full LD-2-D + symbolic INCL slow, AFTER the change:
tests/sn/verification/mms/test_mms_ld_2d.py tests/derivations/test_sn_mms_ld_2d_stress_symbolic.py
  31 passed, 1 warning in 458.32s   (Leg A baseline was 25 passed, 1 skipped)
# net: the 1 #251 skip stub → 6 #251 gates (2 flipped xfail→live + 4 already green).

# Non-slow LD-2-D + symbolic:
tests/sn/verification/mms/test_mms_ld_2d.py tests/derivations/test_sn_mms_ld_2d_stress_symbolic.py -m "not slow"
  28 passed, 3 deselected, 1 warning in 77.28s

# The 6 #251 gates by name:
  6 passed, 17 deselected, 1 warning in 7.98s

# Non-vacuum prescribed-inflow MMS (the scalar-inflow path — slab + sphere):
tests/sn/verification/analytical/test_mms_prescribed_inflow.py
  4 passed, 1 warning in 0.48s   (baseline 4 passed — byte-identical, n_face_moments==1 → identity)

# Bit-identity DriftWarning strict gate (DD/Step + scalar-inflow byte-identical):
tests/sn/sweep/core tests/sn/solve -W "error::tests.sn.regression._regression_assert.DriftWarning"
  520 passed, 1 skipped, 4 xfailed, 2 warnings in 91.77s
# EXACTLY the spec baseline (520/1/4). NO DriftWarning fired, NO golden moved — the
# moment-tail widening leaves the DD/Step + scalar-inflow trace untouched (the
# negative control: face_moment_tail(1) == ()).

# tests/sn/operators (the 7 KNOWN pre-existing baseline reds):
tests/sn/operators                                          → 7 failed, 492 passed, 4 skipped, 1 xfailed
tests/sn/operators -W "error::...DriftWarning"              → 9 failed, 490 passed, 4 skipped, 1 xfailed
# NON-strict = EXACTLY the 7 pre-existing (5 SPHERE #250 + 2 mu_y #232). The +2 strict
# reds are vacuum_bulk_SLB_seed0/1 = WITHIN-TOL 1-ULP DriftWarnings (drifted 1 ULP /
# 1.50e-16 rel, within tol nulp=5), DOCUMENTED #240-era SLB matvec re-association
# (test_bc_extraction_matvec.py committed at f0d68c3 for #240, PRE-#251). My change
# does NOT touch the 1-D slab path (CumprodScan, no MFW/FFW; 1-D trace has no moment
# tail). ZERO new reds.
```

## MUTATION EVIDENCE (the teeth, mutation-verified per anti-pattern #11 / Mode-10)

Regressing the production `_inflow_to_moments` moment-resolved arm to re-zero the
slope rows (`_mut = face.copy(); _mut[..., 1:] = 0.0`) — the EXACT #251 bug it
closes — reddens the new gates while the scalar no-op stays green, then REVERTED:

| mutation                           | threading gate | consumption gate | scalar no-op control |
|------------------------------------|----------------|------------------|----------------------|
| (baseline, no mutation)            | GREEN          | GREEN            | GREEN                |
| re-zero slope rows (#251 bug)      | RED (Assert)   | RED (Failed)     | **GREEN** (asymmetry)|
| (revert)                           | GREEN          | GREEN            | GREEN                |

The Mode-10 asymmetry: the scalar no-op control feeds a scalar inflow (slot-1 ≡ 0
→ re-zeroing zero is a no-op) → correctly blind. The two O(1) teeth: (1) the
producer threads slot-1 at machine precision (stamp proof — RED on re-zero), (2) a
CONSUMED slope flip moves the converged near-bdy flux ≫ tol (consumption proof —
RED on re-zero because there is no slope left to flip).

## CONSUMPTION-MUTATION NUMBERS (via the REAL public production path, NOT the surrogate)

```
nc=16: +slope vs -slope near-bdy |dphi|/|phi| = 4.1012e-03  global = 3.2724e-03
nc=32: +slope vs -slope near-bdy |dphi|/|phi| = 8.3821e-04  global = 8.9896e-04
_CONSUMPTION_TOL = 1e-8
```
Matches the test-architect's probed surrogate values EXACTLY (4.106e-3/3.273e-3 at
nc=16; 8.383e-4/8.990e-4 at nc=32) → production matches the surrogate. The flip
clears `_CONSUMPTION_TOL=1e-8` by ~5.6 orders and HALVES under refinement (O(h),
boundary-localized) — the genuine-consumption signature. The zero-slope
moment-resolved solve is byte-identical to today's scalar prescribed_inflow solve
(maxdiff 0.0 — the Leg-B asymmetry, the bit-identity guard).

## THE REFLECTIVE-PATH VERDICT (storage-correct; SIGN is a follow-up — do NOT claim verified)

Read `boundary_operator.py:_reflect_trace` (219-227) + `boundary_realizer.py`
(184-207). Confirmed BOTH by reading AND empirically (built a reflective-xmin LD-2D
mesh, seeded a nonzero slot-1, ran `_reflect_trace('apply')` on the moment-shaped
trace `(24,2,4,2)`):
- ✅ The reflective `B` realized law is `PermutationOperator(perm, axis=0) &
  IdentityOperator()` (albedo=1.0 folds to `np.take(x, perm, axis=0)`) — permutes
  ordinates on axis 0, broadcasts UNCHANGED over the new trailing moment axis. The
  write `out_boundary.face_view(face)[target] = full[target]` indexes ordinate axis
  0 only → spans all trailing axes. The moment axis passes through for STORAGE
  WITHOUT a code change and WITHOUT a hard-coded trailing-axis assumption (no crash,
  correct `(24,2,4,2)` output). **NO latent bug introduced by the widening** (unlike
  the `prescribed_inflow:274` `, :` assumption, which I fixed).
- ⚠ The transverse-slope SIGN under reflection across a face is NOT verified — the
  Leg-B MMS is vacuum-BC (H2: vacuum nulls the reflective coupling). PHYSICS: a
  normal-flip reflection preserves the tangent-plane (transverse) coordinate, so the
  transverse slope should reflect WITHOUT a sign flip (pass-through PROBABLY correct)
  — but UNVERIFIED. **FOLLOW-UP for the main agent to file: a reflective-LD MMS leg +
  an `op.H` adjoint check on the transverse-slope reflection** (a genuine Mode-1 sign
  trap the vacuum gates cannot see).

## ELEGANCE (coding-elegance)

- Pattern 7 / single source: `face_moment_tail` is the ONE "append iff > 1" policy
  (the trace lever reuses the same primitive the cell cochain `_n_face_moments` /
  `_spatial_moment_tail` key on); `is_moment_valued_by_flat_rank` is the ONE rank
  discriminator (the Leg-A primitive — NO third rank spelling, NIT-1 honored).
- Pattern 4 (illegal states): `_inflow_to_moments` + `prescribed_inflow` reject a
  wrong trailing moment width / a malformed slot — the moment relaxation does not
  swallow real shape bugs.
- Pattern 6 (defer abstraction): NO production projector — the test owns the
  transverse-face projection (leggauss only, L11); production accepts the
  moment-resolved face, does not compute it (exactly Leg A).
- Negative control byte-identity: the DriftWarning strict gate (520/1/4 unchanged) +
  the 1-D prescribed-inflow MMS (4 green) prove DD/Step + scalar-inflow untouched.

## DELIVERABLE MANIFEST (files changed — NOT committed)

- `orpheus/sn/geometry.py` — `boundary_face_layout` moment-tail widening + docstring.
- `orpheus/sn/loss_representation.py` — `_inflow_to_moments` rank-discriminated
  pass-through + width validation; `_octant_face_cochain` moment-resolved IN-edge
  seed; 4 outflow capture-collapse sites dropped; `is_moment_valued_by_flat_rank`
  import added.
- `orpheus/transport/source_sinks/boundary_source_sink.py` — `prescribed_inflow`
  scalar-or-moment slot assignment + `AVERAGE_MOMENT` import + docstring prose.
- `tests/sn/verification/mms/test_mms_ld_2d.py` — 2 xfail-strict gates flipped to
  live (markers removed); `_solve_with_boundary_slope` re-targeted onto the public
  API (monkeypatch removed); docstrings updated.
- `docs/theory/discrete_ordinates.rst` — the `ld-cartesian-2d-slope-source` note
  Leg B "deferred" → "VERIFIED" + the reflective-sign follow-up + a `.. todo::`
  archivist stub for the rich Leg B narrative.

NO algebra-of-record Branch-1/L1/foundation-module manifest — Leg B is a
production-CONSUMPTION widening of the boundary trace (the moment-ready trace was
"moment-ready by accident"), not a new reference solver. The Branch-1 MMS source +
the symbolic gate already exist from D5b-S4; the L11 face projector is test-side
(leggauss). Targeted Sphinx build (to /tmp) confirms MY section is clean (only the
pre-existing `mesh.py` paramref ERROR; Leg B VERIFIED note + todo rendered in HTML).
The main agent batches the full build + Nexus rebuild.

DISPATCH_REQUEST to archivist emitted (followup:false — the rich narrative goes to
the user).

## LESSON (algebra-of-record / vv Mode-10 reinforcement + the producer-rank carve)

1. A SHARPER Mode-10 (the boundary slope is sub-floor for ANY value claim, not just
   the sign — "improves-on-flat" is UNACHIEVABLE, unlike Leg A's bulk slope) is
   closed by STRUCTURAL teeth ONLY: producer-threading at machine precision +
   consumed-flip ≫ tol, paired with the scalar no-op asymmetry. There is NO
   value-improvement leg (no regime makes a boundary-trace slope the dominant
   forcing — the Mode-10 "O(1)-isolating companion" half is UNAVAILABLE). This
   confirms the test-architect's §13 recommendation: a one-line addition to the vv
   Mode-10 row — "if no regime makes the term the dominant forcing (e.g. a
   boundary-trace slope), the producer-threading + consumed-flip structural pair is
   the COMPLETE resolution — no value-improvement leg to add." (NO new mode; NO ERR.)
2. The producer-rank carve: when a carve WIDENS a trace/field slot, audit the
   EXISTING SCALAR producers (the MMS cases), not only the new moment-resolved one —
   the scalar callers feed the SAME widened slot, so the producer must accept BOTH
   ranks (scalar → seed slot 0; moment-resolved → write the full slot). The
   explorer's "1 real producer edit" under-scoped this by 1 (the SCALAR-onto-a-moment
   arm). Reinforces the Leg-A field-space lesson: a rigid scalar contract sitting
   ABOVE a widened slot needs a typed-union relaxation, not just an indexing fix.
