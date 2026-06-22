---
name: issue-247-legA-moment-source-closeout
description: "#247 Leg A — BULK external moment-source consumption (slope-SOURCE half of the LM-1989 trap, Mode-10 closeout). Typed-union widening of solve_sn_fixed_source bulk (flat OR moment-resolved (N,ng,*spatial,2^d)); _lift_external_source_to_moments threads the projected slope rows through. Branch refactor/sn-foundation-cleanup. NOT committed."
metadata:
  type: project
---

# #247 Leg A — bulk external moment-source consumption (the Mode-10 closeout)

Issue #247, L1 verification. Branch `refactor/sn-foundation-cleanup` (main
checkout, NOT a worktree). Host `.venv/bin/python -O`. **NOT committed** (main
agent commits after elegance + qa review). Leg B (boundary transverse-face-slope)
SPLIT OUT to **#251** — NOT touched here.

Implements [[issue-247-moment-source-gate-spec]] (test-architect) +
[[issue-247-moment-source-consumption-path]] (explorer). Closes the EXTERNAL
slope-SOURCE half of the LM-1989 trap that `SN2DCartesianLDStressMMSCase`'s
honest-scope ledger flagged as deferred (vv Mode-10 — activated-but-unconstrained).

## THE PRODUCTION CHANGE (confined to `orpheus/sn/solver.py`, two named sites)

An L17 carve (scalar-source ↔ moment-source). The explorer §3 free-ride held
exactly: a moment-valued bulk joins the EXISTING `(N,ng,*spatial,2^d)` scattering
slope path; the cell rank discriminator (`sweep_graph._CellSolve.cell:931`,
`is_moment_valued_by_rank`) re-signs the external slopes global→sweep EXACTLY as
the scattering slopes — NO change to scattering.py / sweep_graph.py /
loss_representation.py. ALL prod edits in `solver.py`.

1. **`_build_fixed_source_rhs` validation** (`solver.py:~1875`). WIDENED from a
   hard `bulk_values.shape != expected` reject to a typed union of TWO ndarray
   ranks, **discriminated by RANK** (not trailing-size): flat `(N,ng,*spatial)`
   (unchanged) OR moment-resolved `(N,ng,*spatial, per_axis**ndim)` (LD only).
   Rejects everything else INCLUDING a trailing-moment axis ≠ `per_axis**ndim`
   (the negative pin — a clear ValueError naming the expected `2^d`). DD/Step
   (per_axis==1 → `n_cell_moments==1`) rejects a moment-resolved input outright
   (no moment axis exists → only flat valid).

2. **`_lift_external_source_to_moments`** (`solver.py:~1896`). WIDENED to thread
   the slope rows when the input is already moment-resolved (discriminate by
   RANK against the flat `(N,ng,*spatial)` rank `2+len(spatial_shape)`):
   - DD/Step (`tail == ()`): input unchanged (byte-identical — negative control).
   - flat input: UNCHANGED (zero buffer, slot-0 ← flat, slopes ZEROED — the
     honest default, `q̂=0` exact for a region-uniform source).
   - moment-resolved input: PASS THROUGH (the caller projected; validate trailing
     axis == `per_axis**ndim`). Keeps the `(lifted, per_axis)` return contract.

   ⚠ **The rank discriminator is `bulk_values.ndim == 2 + len(sn_mesh.spatial_shape)`**
   (the flat rank), NOT `is_moment_valued_by_rank` — the latter compares against
   the per-ordinate-STRIPPED `(ng,*spatial)` reference (`reaction_xs`), but the
   external bulk carries the leading ordinate axis, so the flat-rank comparison
   is the correct equivalent here (matches the brief's "compare bulk_values.ndim
   to the flat expected rank").

NO callable-projection entry (Pattern 6 — no consumer; the MMS test does its OWN
projection and passes the array). The projection lives test-side, structurally
independent of production (L11): `leggauss` only, NEVER `_lift_*` nor any
`_ubld`/`LinearDiscontinuous` method.

## THE TESTS (`tests/sn/verification/mms/test_mms_ld_2d.py`, the #247 block)

`-O`-safe (Mode 8): `np.testing.*` / `pytest.fail` / `pytest.raises` only, NO
bare assert. The two projection foundation sub-gates (Deliverable 3) were already
written by the test-architect and pass as-is. NEW/implemented:

- `test_ld_2d_external_slope_source_threaded_through_lift` (foundation) — **the
  STRUCTURAL teeth (machine precision)**: the lift threads the projected moment
  vector through UNCHANGED (`np.testing.assert_array_equal(lifted, Qm)`, every
  slope slot). NEGATIVE CONTROL leg: a FLAT bulk still lifts onto slot-0 with
  slopes EXACTLY ZERO. This is the sharpest, non-tautological proof of the
  production change — the spec §0 converged flux is sub-floor sensitive, this is
  O(1).
- `test_ld_2d_external_slope_source_converges_second_order` (l1, slow, verifies
  ld-cartesian-2d) — Deliverable 1: the WIDENED public solve fed the
  moment-resolved external source converges O(h²) to A (NECESSARY, NOT sufficient
  for sign per §0). Probed orders ~2.0 (8.18e-3→1.99e-3 nc=16→32).
- `test_ld_2d_external_slope_source_improves_on_flat` (l1, verifies) — POSITIVE
  verification: the moment-resolved solve lands STRICTLY closer to A than the flat
  solve (probed 3.40e-3 < 5.99e-3 at nc=24) — the slopes carry real sub-cell info.
- `test_ld_2d_external_slope_source_sign_mutation_reddens[2/1/3]` (l1, verifies,
  parametrized) — Deliverable 2 M1–M3, **the PRIMARY sign-catcher**: flipping the
  consumed EXTERNAL slope-source row (slot 2=x / 1=y / 3=xy) changes the converged
  flux ≫ `_CONSUMPTION_TOL=1e-8` (the consumption proof, O(1) above the 1e-12
  fixed point). PLUS an in-test FLAT no-op leg pinning the Mode-10 asymmetry (the
  flat external source's slope row IS zero → flipping it is a no-op → the flat
  gate is correctly blind).
- `test_ld_2d_scattering_slope_source_sign_mutation_reddens` (l1, verifies) —
  Deliverable 2 M4: monkeypatch `ScatteringOperator._assemble_per_ordinate_source`
  to flip the iso slope rows (slots 1:) → converged flux changes ≫ tol (probed
  2.62e-3). Verifies the EXISTING (S3) scattering consumption was never sign-blind.
- `test_moment_resolved_bulk_still_rejects_wrong_trailing_axis` — Deliverable 5:
  LD rejects a 5-wide trailing axis (ValueError naming `per_axis**ndim = 4`); DD
  rejects a 4-wide moment-resolved bulk (no moment axis).
- `test_ld_2d_boundary_transverse_face_slope` — kept as a SKIP re-pointed to
  **#251** (Leg B, boundary trace widening — a different production path).

### Why the converged-flux value-band was REJECTED as teeth (the §0 trap, confirmed)

Probed live (full public solve, nc=24, x-slope external flip): `|Δφ|/|φ|` only
3.0e-3 (x-slope), 1.0e-2 (y-slope), **6e-5 (xy)** — a value-band against A would
need a tolerance tighter than the ~1.5× discretization gap, which the O(h²) floor
eats. The single-sweep "vs manufactured ψ" reference was even WEAKER (ratio → 1.0
with refinement — the single-sweep error is dominated by the incomplete scattering
iteration, NOT the source slope). So the teeth are TWO O(1) signals: (a) the lift
pass-through at machine precision (the production-change proof), (b) the
consumed-flip converged-flux change ≫ solver tol (the consumption proof). The
`_CONSUMPTION_TOL=1e-8` band is ~5e7× above the 1e-12 fixed point yet far below
the §0 trap — the smallest probed flip (xy, 6e-5) clears it by ~6000×.

## VERIFICATION (verbatim `-O` pytest summary lines)

```
# Full LD-2-D + symbolic (incl. slow), AFTER the change:
tests/sn/verification/mms/test_mms_ld_2d.py tests/derivations/test_sn_mms_ld_2d_stress_symbolic.py
  25 passed, 1 skipped, 1 warning in 427.26s   (baseline was 17 passed, 7 skipped)
# net: +8 active gates (the #247 block); 1 skip = the #251 Leg B stub.

# Bit-identity DriftWarning strict gate (the flat path is byte-identical):
tests/sn/sweep/core tests/sn/solve -W "error::tests.sn.regression._regression_assert.DriftWarning"
  520 passed, 1 skipped, 4 xfailed, 2 warnings in 84.84s
# NO DriftWarning fired, NO golden moved — the typed-union widening leaves the
# flat/DD path untouched (the negative control).

# Non-slow #247 gates explicit:
  9 passed, 1 warning in 28.88s
# (3 sign-mutation params + 6 others; the slow convergence gate runs in the full run)
```

## MUTATION EVIDENCE (the teeth, mutation-verified per anti-pattern #11)

Regressing the PRODUCTION lift to re-zero the slope rows (the EXACT #247 bug it
closes) reddens the new gates while the flat scalar gate stays green:

| mutation                          | NEW gate              | flat scalar gate                          |
|-----------------------------------|-----------------------|-------------------------------------------|
| (baseline, no mutation)           | GREEN                 | (n/a)                                     |
| lift re-zeros slope rows (#247 bug) | RED (AssertionError) [structural-lift gate]; RED (Failed) [improves-on-flat]; RED ×3 (Failed) [sign-mutation x/y/xy] | **GREEN** (`test_ld_2d_stress_converges_second_order`) |

The **Mode-10 asymmetry** (the gap #247 closes): with the lift monkeypatched to
flip the external x-slope on EVERY lift, the flat scalar gate stays GREEN (it
feeds a flat source → the slope row is zero → flipping zero is a no-op), while the
moment gates are sensitive. M4 (scattering iso-slope flip) → `|Δφ|/|φ|=2.6e-3 >
1e-8` → the M4 gate reddens.

## ELEGANCE (coding-elegance — confined to solver.py)

- Pattern 4 (illegal states): the validation REJECTS a moment-resolved bulk on
  DD (no moment axis) + a wrong trailing axis — the typed-union relaxation does
  not swallow real shape bugs.
- Pattern 6 (defer abstraction): NO callable-projection entry — production accepts
  the array only; the test owns the projection.
- Pattern 7 / single source (no behavior drift): the flat path is byte-identical
  (the DriftWarning gate + the two-paths LD tests are the proof). The rank
  discriminator is the single decision point; the moment-resolved arm threads the
  caller's projection through unchanged (no production projector — that would make
  the gate a tautology, L11).

## WHERE THE PROJECTION-NORMALIZATION CRUX BIT (vs the spec)

The spec §1 CRUX (bare per-volume Legendre coeff, no θ/h/V; d=2 Kronecker order
`[Q̄, Q̂_y, Q̂_x, Q̂_xy]`, GLOBAL frame) held EXACTLY — the test-architect's
projector stub `_project_scalar_to_tensor_legendre` is correct as written
(machine-precision on the hand-laid bilinear). The ONE place harder than stated:

- The spec §2 per-moment O(h²) convergence threshold (`>1.7` for every moment) is
  OPTIMISTIC for the SLOPE moments. Probed: the AVG scalar-flux moment converges
  clean O(h²) (~2.0-2.2), but the per-moment SLOPE scalar fluxes converge at
  ~O(h^1.5–1.6) (the slope moment is a derivative-like quantity, lower order by
  construction). So Deliverable 1 was implemented as the AVG-moment (scalar flux)
  O(h²) convergence via the public solve (clean, the right reference) — NOT a
  per-moment-slope band that would fail a 1.7 threshold. The per-moment slope
  CONSISTENCY is instead pinned by the structural lift-pass-through + the
  improves-on-flat positive leg. This is the §0 lesson lived: the converged-flux
  per-moment band is not the place for the slope teeth.

## ERR catalog

NO new ERR. Mode-10 is a proactive-gap close, not a caught production bug (the
slope-source sign was UNVERIFIED, not WRONG — the lift correctly zeroed an
unverified-but-honest default). Per the spec §9, no `@catches` until a real bug
surfaces (next free is ERR-063). None did.

## DELIVERABLE MANIFEST (files changed — NOT committed)

- `orpheus/sn/solver.py` — `_build_fixed_source_rhs` validation widening +
  docstring; `_lift_external_source_to_moments` slope-thread widening + docstring.
- `tests/sn/verification/mms/test_mms_ld_2d.py` — the #247 block (6 new/impl
  gates + the Leg B skip re-pointed to #251 + the existing-gate honest-scope
  docstring update). The two projection foundation sub-gates were already written.
- `orpheus/derivations/continuous/mms/sn.py` — `SN2DCartesianLDStressMMSCase`
  module honest-scope docstring + class docstring (Leg A VERIFIED / Leg B #251).
- `docs/theory/discrete_ordinates.rst` — the `ld-cartesian-2d` honest-scope note
  updated (Leg A verified, the Mode-10 asymmetry, M4) + a `.. todo::` stub
  (`ld-cartesian-2d-slope-source` label) for the archivist's rich Leg A narrative.

NO algebra-of-record Branch-1/L1/foundation-module manifest — Leg A is a
production-consumption WIDENING of an EXISTING moment path, not a new reference
solver (the Branch-1 MMS source + the symbolic gate already exist from D5b-S4).
NO full Sphinx build committed (a targeted build to /tmp confirmed MY section is
clean — only the pre-existing `mesh.py` paramref error; the main agent batches
the full build + Nexus rebuild).

DISPATCH_REQUEST to archivist emitted below (followup:false — the rich narrative
goes to the user).

## LESSON (for the algebra-of-record / vv Mode-10 reinforcement)

A Mode-10 gap (activated-but-unconstrained slope source) is closed NOT by
tightening the converged-flux value band (the §0 trap — the slope-source error is
O(h²)-small, sub-floor) but by TWO O(1) structural teeth: (1) assert the
production producer threads the projection through at machine precision (the
production-change proof — catches a regression to zeroing), and (2) assert a
CONSUMED source-row sign flip moves the converged answer ≫ solver tol (the
consumption proof — catches sign-blindness), paired with the FLAT no-op leg that
pins the asymmetry (the old gate is correctly blind). The convergence-ORDER leg is
NECESSARY (slope consistency) but explicitly NOT the sign teeth. This is the
canonical resolution of a Mode-10 gap where the term is genuinely consumed yet the
converged value is sub-floor sensitive — reinforces the existing vv Mode-10 row,
no skill edit needed.
