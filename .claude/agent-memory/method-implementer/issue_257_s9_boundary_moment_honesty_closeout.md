---
name: issue-257-s9-boundary-moment-honesty-closeout
description: "#257 S9 — production MMS boundary-source honesty + the LD coherent-promise gate. SN2DCartesianLDStressMMSCase.prescribed_inflow now EMITS the moment-resolved face slot (slot-0 transverse cell AVERAGE, slot-1 bare transverse P1 slope) via a case-owned leggauss-only _project_inflow_to_face_moments, gated on face_moment_count>1 (DD/Step byte-identical). Diagnostic promoted to tests/sn/verification/mms/test_ld_2d_boundary_promise.py (coherent-promise gate + sub-floor verdict pins + Mode-11 sentinel). #251 GATE-B/C legs re-targeted onto the production producer. NO new field type (#263), NO value gate (sub-floor), NO DD/Step change. Branch feature/field-typed-operator-algebra, HEAD 8e0a2cf. NOT committed."
metadata:
  type: project
---

# #257 S9 — LD boundary moment honesty + the coherent-promise gate

Issue #257 S9 (the final behavioral stage), L1 verification. Branch
`feature/field-typed-operator-algebra`, HEAD `8e0a2cf` (NOT committed — the main
agent commits after elegance + qa). Host `.venv/bin/python -O`. The natural
COMPLETION of [[issue-251-legB-boundary-trace-closeout]]: #251 Leg B landed the
production-CONSUMPTION path (the moment-resolved trace + `_inflow_to_moments`
pass-through end to end); the MMS case `prescribed_inflow` still DROPPED the
slope (built a scalar, hit the producer's scalar branch). S9 closes that
producer-blindness — the case now EMITS the moment slot.

Settled context (NOT re-derived): the boundary transverse-slope is SUB-FLOOR for
the converged flux in every regime (verdict double-confirmed, test-architect +
numerics-investigator, across optical depth + frequency + amplitude). The
coherent promise ("LD is 2nd-order everywhere incl. boundary") is ALREADY TRUE —
delivered by the AVERAGE moment. So verification is STRUCTURAL: NO value/order
gate keyed on the slope; the boundary moment is a PROPERTY not a type (#263).

## DELIVERABLE 1 — production MMS source honesty (`orpheus/derivations/continuous/mms/sn.py`)

`SN2DCartesianLDStressMMSCase.prescribed_inflow` rewritten to gate on
`face_moment_count(sn_mesh.scheme.spatial_basis_per_axis, sn_mesh.ndim)` (the
SAME single-source primitive `SNMesh.boundary_face_layout:1072` keys on):

- **`n_face_moments == 1` (DD/Step):** builds the SCALAR per-face trace
  `(N, ng, n_t)` (cell-CENTRE eval of `(A + μ_x B + μ_y C)/W`), hits the
  producer's scalar branch (`boundary_source_sink.py:286-287`) → **BYTE-IDENTICAL**
  to pre-S9. PROVEN: `np.array_equal(bss_dd.values, <pre-S9 face_coords build>)
  == True` (the DD/Step path produces `(24,2,6)`/`(24,2,8)` slots, no moment axis,
  bit-for-bit the old build).
- **`n_face_moments > 1` (LD):** builds the FULL moment slot
  `(N, ng, n_t, n_face_moments)` via a NEW case-owned
  `_project_inflow_to_face_moments(const_axis, const_val, t_edges, n_face_moments)`,
  hits the producer's full-slot branch (`boundary_source_sink.py:284-285`,
  `arr.shape == view.shape` → `view[inflow] = arr[inflow]`).

`_project_inflow_to_face_moments` projects `ψ_{n,g}(face, t)/W` onto the BARE
per-transverse-cell Legendre moments on `[t_edges[j], t_edges[j+1]]` mapped to
`ξ∈[-1,1]`: slot-0 (`AVERAGE_MOMENT`) `= ⟨ψ,P₀⟩/⟨P₀,P₀⟩` (transverse cell
AVERAGE), slot-1 `= ⟨ψ,P₁⟩/⟨P₁,P₁⟩` (BARE slope, NO θ/h_t — the cochain's
transverse mass `diag(h_t, θ·h_t)` adds them downstream; a θ/h_t-weighted slope
would double-apply the mass = a TRUE bug GATE C catches). **L11:** descends ONLY
from `self._drivers` + `numpy.polynomial.legendre.leggauss(6)` — NEVER
`_inflow_to_moments`/`_ubld`/any LD op/the test-side projectors.

`build_nonvacuum_fixed_source` (`:3552`) passes `case.prescribed_inflow` through
unchanged — VERIFIED no change needed (it just bundles bulk⊕boundary).

**PROBED slot agreement (the producer stamp):** on the correct per-face inflow
ordinate sets (`sn.trace.inflow_indices_for_face`), production slot-1 == the
test-side `_face_transverse_buffers` leggauss reference EXACTLY (0.0 maxdiff) on
ALL 4 faces; outflow slots exactly zero; slot-0 differs only O(h²) (cell-AVERAGE
vs cell-CENTRE, expected). LD slots `(24,2,6,2)`/`(24,2,8,2)`.

## DELIVERABLE 2 — promote the diagnostic to a permanent gate

`derivations/diagnostics/diag_s9_ld_boundary_slope_optical_sweep.py` (untracked,
numerics-investigator) → `tests/sn/verification/mms/test_ld_2d_boundary_promise.py`
(REMOVED the diagnostic; it was untracked, plain `rm`, NOT a git op). Carried
over the 4 tests, Mode-8 fix (bare-`assert` order-band → `pytest.fail`; the
`np.testing` sentinel kept), rich module docstring kept (the durable verdict +
mechanism record). **Import-shim fix:** dropped the `sys.path.insert`; imports
`_solve_with_boundary_slope` as a **package-relative SIBLING**
`from .test_mms_ld_2d import ...` (NOT bare `from test_mms_ld_2d` — the mms dir
has `__init__.py` + pytest `--import-mode=importlib`, so a bare sibling import
`ModuleNotFoundError`s; the relative form resolves). pyright `reportMissingImports`
clears (the new file is 0/0).

- **Coherent-promise gate** `test_first_cell_row_already_second_order` —
  `@l1` + `@verifies("ld-cartesian-2d")` (REUSED label, NOT minted). PROBED:
  first-cell-row flat orders 1.993/2.004/2.001, mom 1.998/2.005/2.003 (both
  O(h²) AT the boundary → the no-asterisk promise, delivered by the average).
- **Sub-floor verdict pins** `test_optical_sweep_slope_never_beats_floor`
  (6-param) + `test_amplified_boundary_slope_still_subfloor` (3-param) —
  `@slow` + `@l1`. Guard the no-value-gate conclusion (red if the LD boundary
  closure ever changes so the slope matters above-floor). Docstrings point to
  #257 S9 + #263.
- **Mode-11 sentinel** `test_slope_toggle_reaches_inflow_to_moments` —
  `@foundation`. PROBED: flat/zero see slot1_max=0; mom/flip see 1.9e-2;
  `phi_sum` differs (consumed); zero==flat byte-identical; `_inflow_to_moments`
  reached 688×/solve.

## DELIVERABLE 3 — re-target #251 Leg-B gates onto the PRODUCTION producer (GATE B + C)

In `tests/sn/verification/mms/test_mms_ld_2d.py`:
- **GATE B** — added a leg to `test_ld_2d_boundary_slope_threaded_through_inflow_to_moments`
  (foundation): the PRODUCTION `case.prescribed_inflow(sn)` emits a slot whose
  SLOT-1 (on each face's inflow ordinates) == the leggauss reference at machine
  precision (`assert_array_equal`). Closes the producer-blindness the #251
  surrogate had (it pinned only the CONSUMER `_inflow_to_moments`; this pins the
  PRODUCER).
- **GATE C** — NEW `test_case_projector_agrees_with_test_face_projector`
  (foundation): the case `_project_inflow_to_face_moments` SLOT-1 == the
  test-side `_face_transverse_legendre` slope on the manufactured trace, per
  cell/ordinate/group at machine precision (`assert_array_equal`). Single-source
  agreement of the two structurally-independent projectors (Cardinal Rule 2 /
  L11) — pins the bare-coeff normalization at the projector level.
- `test_ld_2d_boundary_slope_sign_mutation_reddens` (l1) + the scalar no-op
  control — KEPT GREEN. ⚠ The re-baseline interaction (below) required updating
  `_solve_with_boundary_slope`.

## ⚠ THE RE-BASELINE INTERACTION (the load-bearing test-side fix)

The diagnostic's `_solve_with_boundary_slope(case, nc, slope_sign=None)` USED to
mean "average-only" by calling `build_nonvacuum_fixed_source` → the
AVERAGE-ONLY `case.prescribed_inflow`. After S9 the production `prescribed_inflow`
EMITS the real slope on an LD mesh → `None` would silently become
"average + real slope", collapsing the diagnostic's flat/mom/flip distinction AND
breaking the consumption gate's `phi_zero == phi_today` byte-identity assertion.

FIX (in `test_mms_ld_2d.py`, the test file): `_solve_with_boundary_slope` now
builds ALL branches TEST-SIDE from `_face_transverse_buffers` —
`slope_sign=None` → a SCALAR `(N,ng,n_t)` inflow (producer seeds slot-0, slope
stays zero, the controlled average-only baseline); `±1.0/0.0` → the moment slot
as before. This decouples the diagnostic's average-only baseline from the
production slope honesty (the helper is the CONTROLLED toggle; production's
honesty is the SEPARATE thing under test). `build_nonvacuum_fixed_source` import
STAYS (still used by 3 other consumers at `:332/406/837/940`). Consumption gate
re-run GREEN (`phi_zero == phi_today` byte-identical holds, consumed flip ≫ tol).

## RE-BASELINE CANDIDATES FOR QA (the LD tests that shift sub-floor)

Tests whose BOUNDARY now carries the real projected slope (via
`build_nonvacuum_fixed_source` or `_solve_moment_resolved` → `case.prescribed_inflow`).
ALL confirmed: the shift is SUB-FLOOR (O(h²) preserved); ORDER assertions stay
green; **NO tight value/snapshot pin moved** (the value bands have ~7 decades of
headroom; no `np.array_equal`/snapshot/DriftWarning pins consume the LD-stress
case). qa to VERIFY:

1. `test_ld_2d_stress_converges_second_order` (`:349`, l1/slow,
   verifies ld-cartesian-2d) — order `[>1.95, all>1.85]` + value band
   `1e-9 < err < 1e-2`. Boundary now slope-bearing; shift sub-floor → STAYS green.
2. `test_ld_2d_stress_krylov_equals_si` (`:390`) — both legs share the SAME rhs
   boundary → Krylov≡SI relative comparison UNAFFECTED.
3. `test_ld_2d_external_slope_source_converges_second_order` (`:792`, l1/slow) —
   `_solve_moment_resolved` boundary; order `[>1.9, all>1.8]`. Sub-floor → green.
4. `test_ld_2d_external_slope_source_improves_on_flat` (`:820`, l1) — `phi_mom`
   AND `phi_flat` BOTH go through `case.prescribed_inflow` boundary (now
   slope-bearing on both) → the bulk-moment-vs-flat comparison UNAFFECTED.
5. `test_ld_2d_external_slope_source_sign_mutation_reddens` (`:873`, l1, ×3) —
   both legs share the boundary → relative flip comparison UNAFFECTED.
6. The FFW≡MFW oracle (`_solve_moment_resolved` ancestor, `bss` copy `:469`) —
   both legs share `bss` → relative comparison UNAFFECTED.

ALL 35 gates in `test_ld_2d_boundary_promise.py` + `test_mms_ld_2d.py` PASS
(639.71s). NO test required a value-pin re-baseline — the sub-floor shift is
absorbed within existing bands. qa: confirm the order bands stay green on a fresh
run; flag if any DriftWarning/snapshot pin reds (none expected — no LD-stress
consumers in sweep/core, solve, or regression).

## DELIVERABLE 4 — minimal accessor: SKIPPED (does not earn its keep)

The only slope-index in S9's production code is `slot[..., 1]` inside
`_project_inflow_to_face_moments` (an immediate construction context where the
local LD `[bar, slope]` order is self-documenting + matches the existing
test-side projectors `_face_transverse_legendre`/`_face_transverse_buffers`,
which use literal `[..., 1]` without a named constant). slot-0 IS named
(`AVERAGE_MOMENT`, single-sourced via the local import). There is NO scattered
`slot[..., 1]` read across consumers (the consumer `_inflow_to_moments` threads
the WHOLE slot, not slot-1 by index). A `FACE_SLOPE_MOMENT` constant or a
`BoundaryField` accessor would be gold-plating for one local write. The broad
typed `SpatialMomentSpace` moment-axis predicate is **#246**, NOT S9.
⚠ ELEGANCE DECISION POINT — elegance-enforcer to confirm the skip.

## VERIFICATION (verbatim outputs)

```
# (1) npx pyright — net-new vs baseline 2282:
2282 errors, 19 warnings, 0 informations
#   = EXACTLY baseline 2282, 0 net-new (the diag import-shim +1 cleared on
#   promotion; the new test file is 0/0; 0 new `# type: ignore`).

# (2) .venv/bin/python -O -m pytest -q
#     tests/sn/verification/mms/test_ld_2d_boundary_promise.py
#     tests/sn/verification/mms/test_mms_ld_2d.py
35 passed, 1 warning in 639.71s (0:10:39)
#   (#251 baseline was 31; +4 from the promoted module's 4 gate FUNCTIONS:
#   coherent-promise + 6-param + 3-param verdict pins + Mode-11 sentinel;
#   minus the deselected — net 35.)

# (3) .venv/bin/python -O -m pytest -q
#     tests/sn/verification/mms tests/sn/sweep/core tests/sn/solve
#     --deselect tests/sn/solve/test_keff_slab.py::test_heterogeneous_absolute_keff
590 passed, 1 skipped, 4 xfailed, 3 warnings in 802.26s (0:13:22)
#   ZERO failures, ZERO non-baseline reds (route around #250/#232/#212).

# GATE D — DriftWarning strict (DD/Step + scalar-inflow byte-identity):
# .venv/bin/python -O -m pytest -q tests/sn/sweep/core tests/sn/solve
#   -W "error::tests.sn.regression._regression_assert.DriftWarning"
#   --deselect tests/sn/solve/test_keff_slab.py::test_heterogeneous_absolute_keff
520 passed, 1 skipped, 4 xfailed, 2 warnings in 88.04s (0:01:28)
#   = EXACTLY the #251/#247 baseline (520/1/4); NO DriftWarning fired, NO golden
#   moved — the moment-tail emission affects ONLY LD; DD/Step byte-identical.

# DD/Step byte-identity (focused probe):
#   np.array_equal(case.prescribed_inflow(DD).values, <pre-S9 face_coords build>)
#   == True  (DD slots (24,2,6)/(24,2,8), no moment axis, bit-for-bit).

# 1-D prescribed-inflow MMS (slab/sphere, own prescribed_inflow, untouched):
#   tests/sn/verification/analytical/test_mms_prescribed_inflow.py → 4 passed.
```

## ELEGANCE / DOCS / DECISION POINTS (for elegance + qa)

- **NO new Sphinx label** (`ld-cartesian-2d` REUSED, exists `discrete_ordinates.rst:3084`).
  The rich coherent-promise narrative + the #263 criterion seam is OWED to the
  archivist's consolidated pass (S3b/S5/S6/S8b/S8c + S9 fold) — NOT written here
  (brief: "DO NOT write Sphinx narrative"). DISPATCH_REQUEST emitted (followup:false).
- **Stale-prose updates** (S9 closed the "deferred" state): the case module
  HONEST SCOPE + class docstring + `external_source` docstring + the symbolic
  test `test_ld2d_stress_prescribed_inflow_is_nonvanishing` docstring — all
  "deferred S4 / Leg B deferred / lifts only the AVERAGE moment" prose updated to
  the now-CARRIED state. (In-module honesty, not Sphinx narrative.)
- **NO new field type** (#263), **NO value/order gate** for the slope (sub-floor),
  **NO `BoundaryFlux` retype**, **NO projection duplication** (the production
  projector is structurally independent of the test projectors; GATE C pins their
  agreement, NOT their identity), **NO DD/Step change** (byte-identical, proven).
- **NO algebra-of-record Branch-1/L1 manifest** — S9 is a production-PRODUCER
  completion of the #251 consumption path (the Branch-1 MMS + symbolic gate exist
  from D5b-S4; the L11 projector is case-owned leggauss). Same posture as
  #247/#251/S3a/S3b/S5/S6.
- **NO new ERR** — S9 closes a producer-blindness (the slope was UNVERIFIED at
  the PRODUCER, not WRONG — the #251 consumer already threaded it correctly). No
  caught production bug; ERR-063 stays reserved.

## DELIVERABLE MANIFEST (files changed — NOT committed)

- `orpheus/derivations/continuous/mms/sn.py` — `prescribed_inflow` gated
  scalar-vs-moment build + NEW `_project_inflow_to_face_moments` (leggauss-only,
  L11) + 4 stale-prose docstring updates.
- `tests/sn/verification/mms/test_ld_2d_boundary_promise.py` — NEW (promoted
  diagnostic: coherent-promise gate + 2 verdict pins + Mode-11 sentinel; relative
  sibling import).
- `tests/sn/verification/mms/test_mms_ld_2d.py` — GATE B leg (production-producer
  stamp) + NEW GATE C gate (projector single-source) + `_solve_with_boundary_slope`
  re-baseline fix (None→test-side scalar).
- `tests/derivations/test_sn_mms_ld_2d_stress_symbolic.py` — stale docstring update.
- REMOVED `derivations/diagnostics/diag_s9_ld_boundary_slope_optical_sweep.py`
  (untracked; promoted).

## LESSON (vv Mode-10 / producer-honesty completion)

When a CONSUMPTION-widening (#251 Leg B) lands the moment-resolved trace end to
end but the PRODUCER (the MMS case `prescribed_inflow`) still emits a scalar, the
producer-blindness is closed by making the case EMIT the moment slot — gated on
`face_moment_count > 1` so DD/Step stays byte-identical. The load-bearing
SECOND-ORDER effect: a diagnostic helper whose `slope_sign=None` branch routed
through the (now slope-honest) production `prescribed_inflow` SILENTLY inherits
the new slope, collapsing the controlled flat/mom/flip toggle and breaking the
byte-identity no-op control — the helper's average-only baseline MUST be built
TEST-SIDE (the controlled toggle is the test's, the honesty is production's; do
not let the helper inherit the production change it is meant to be orthogonal to).
This is the THIRD Mode-10 sub-floor-companion-unavailable instance (#240 D5b-S4 →
#247 Leg A → #251 Leg B → S9); the vv Mode-10 row already carries the
"no-dominant-forcing-regime → structural pair is complete" branch — S9 REINFORCES
it, no skill edit. Extends [[issue-251-legB-boundary-trace-closeout]].
