---
name: issue-257-s9-ld-boundary-slope-optical-verdict
description: "#257 S9 — does the 2-D Cartesian LD boundary TRANSVERSE-SLOPE moment beat the average-only drop above the O(h²) floor in ANY regime (incl. the optically-thin/streaming limit the TA never tested)? VERDICT: NO. flat≡mom≡flip O(h²) across Σ_t·L∈{0.1..2.0} × c∈{0,0.5} × bc_scale∈{1,5,20} × LS-S4/S8 incl nc=128. Slope makes converged flux MONOTONICALLY WORSE as it grows (bc_scale=20: mom 17% worse, still O(h²)). Mechanism: first-cell-row ALREADY O(h²) with avg-only → no O(h) deficiency for slope to repair. Mode-10 companion-unavailable confirmed FUNDAMENTAL; TA right, its codim-1/O(h^1.5) mechanism wording WRONG."
metadata:
  type: project
---

# #257 S9 — LD boundary transverse-slope: optical-depth verdict

Branch `feature/field-typed-operator-algebra`, HEAD `8e0a2cf`. Host `.venv/bin/python`.

## The question & the verdict

Does INCLUDING the boundary transverse-slope moment (vs the average-only drop)
improve the CONVERGED-FLUX error magnitude or spatial ORDER above the O(h²)
floor — in ANY physically-meaningful regime? **VERDICT: NO. No value/order MMS
gate is achievable. The TA's conclusion is CONFIRMED and now ROBUST against the
untested optical-depth axis; its STATED MECHANISM was wrong and is corrected.**

## What I tested (the axis the TA missed)

The TA varied only transverse FREQUENCY (K=3,6). I varied **optical depth** with
a custom uniform-Σ scalable `SN2DCartesianLDStressMMSCase` (constant
`sigma_t_fn`/`sigma_s_fn`; MMS source = PDE residual stays exact for any Σ):

- Σ_t·L ∈ {0.1, 0.5, 1.0, 2.0} (thin/streaming → thick), c ∈ {0, 0.5}
- pure-absorber streaming (c=0) — the regime where the inflow propagates
  ballistically and an inflow-rep error is LEAST boundary-confined (the user's
  intuition's best hope)
- grazing-heavy LS-S8 vs S4; nc = 16/32/64/128 (asymptotic)
- bc_scale ∈ {1,5,20} — amplified the μ-dependent (B,C, slope-carrying) drivers
  so the boundary slope/avg ratio is dominant (the user's exact "non-trivial
  slope at the boundary" hypothesis, stressed to the limit)

**Every regime: flat ≡ mom ≡ flip at O(h²) globally (~2.00) and ~O(h^2.4) at the
boundary ring.** The slope is a sub-floor perturbation that makes the converged
flux SLIGHTLY WORSE, monotonically worse as it grows (bc_scale=20: mom global
1.47e-2 vs flat 1.26e-2 = +17%, STILL O(h²)). It NEVER beats flat. Streaming
limit (Σ_t·L=0.1) is no different from thick.

## The MECHANISM (corrects the TA)

**Decisive scaling probe:** the FIRST-CELL-ROW (i=0, cells directly consuming
the xmin inflow) is ALREADY O(h²) with average-only inflow (orders
1.993/2.004/2.001 across nc 16→128). There is NO O(h) deficiency in the
converged flux for the slope to repair. The slope DOES lift the *inflow
representation* O(h)→O(h²) (the TA's Probe D is real), but that is a 2nd-order
correction to an ALREADY-2nd-order LD cell BALANCE (the cell integrates the
inflow against its own linear basis; the cell-AVERAGE moment is what that
integral needs to O(h²)). So the slope cannot move the converged flux above the
bulk O(h²) floor — including in streaming, because the LD balance integrates the
inflow regardless of how far it propagates.

**TA was RIGHT (no gate) but its mechanism wording was internally inconsistent:**
it claimed "O(h) boundary forcing → O(h^1.5) global, sub-dominant to O(h²)" —
but O(h^1.5) DOMINATES O(h²). The correct statement: the average inflow is
**O(h²)-ADEQUATE for the converged flux in every regime**, so the slope is a
sub-floor refinement, never a deficiency-repair. Canonical Mode-10
companion-unavailable: the boundary transverse-slope has NO regime where it is
the dominant forcing in the converged flux — FUNDAMENTAL to a boundary-trace
moment, not regime-specific.

## Methodology that cracked it

- **Mode-11 sentinel FIRST** (before trusting any null): monkeypatched
  production `_LossRepresentation._inflow_to_moments` to record slot-1 — proved
  flat/zero see slot-1=0, mom/flip see 1.9e-2, and phi_sum DIFFERS (slope
  genuinely consumed: flat 60.8225557, mom 60.8616668, flip 60.7834446, zero==flat
  byte-identical). The toggle is non-vacuous → the null result is REAL.
- **AMPLIFY the suspect term** (bc_scale) — the sharpest disproof of "make the
  slope non-trivial": growing the slope made mom strictly WORSE, the cleanest
  possible refutation of the improves-on-flat hypothesis.
- **First-cell-row order is the mechanism oracle** — if avg-only were O(h)-
  deficient at the boundary cell, the slope COULD repair it; it is O(h²), so it
  can't. This is the decisive scaling argument, not the global L2 (which dilutes).

## Deliverable

`derivations/diagnostics/diag_s9_ld_boundary_slope_optical_sweep.py` (11 tests,
~163s under -O). PROMOTE to `tests/sn/verification/mms/` as the S9 verdict-pin
(if the LD boundary closure ever changes so mom beats flat above-floor, the
verdict must be revisited). ⚠ Promotion MUST convert the bare-`assert` gate
logic to `np.testing`/`pytest.fail` (Mode-8: -O strips bare asserts; verified
the gates fire & pass without -O, but a promoted gate runs under -O).

## S9 consequence

S9 stays STRUCTURAL (the #251 Leg B teeth: machine-precision threading +
consumed-flip mutation + scalar no-op control, ALREADY LIVE). NO value/order MMS
gate. Reframe the motivation from "recover 2nd order at the boundary" (already
true from the average moment) to "type + thread the inflow-representation
refinement the LD closure consumes". This is the 4th Mode-10 sub-floor-
companion-unavailable instance (#240 D5b-S4 → #247 Leg A → #251 Leg B → S9); no
new vv mode, no ERR (no caught production bug).
