---
name: issue-239-scanmarch-coeff-lift
description: #239 2-D ScanMarch coefficient-model lift is ALREADY LANDED on refactor/sn-foundation-cleanup (commit 66dbd9a, #240 D5a); the inline-diamond premise is STALE. LD-2-D is a structural exclusion (bilinear ≠ facewise), not underived research.
metadata:
  type: project
---

# #239 — 2-D ScanMarch coefficient-model lift: ALREADY DONE on this branch

**Fact.** On `refactor/sn-foundation-cleanup` (graph commit `f6475e5`; working tree CLEAN
for all SN spatial/loss files), the #239 lift the issue describes was ALREADY performed —
landed under the #240 D5a campaign in commit **`66dbd9a` "refactor(sn): fold the 2-D
ScanMarch onto the scheme coefficient model (#240 D5a, #239)"**. The GitHub issue #239 is
still OPEN (module:sn, type:improvement), but its premise ("`_scanmarch_row` inlines the
diamond `0.5`/`2.0`; `_sweep_interior` computes DD coefficients inline; `_scanmarch_matvec`
inlines the diamond reconstruction") is STALE — that was the pre-D5a state.

**Why:** the issue text was written before #240 D5a; the natural trigger it names (a 2-D LD
`cell_kernel_batch`) drove the lift early. Issue body still says "DD-2-D-only, inline."
**How to apply:** treat #239 as a CLOSE-VERIFY (regression-pin + close the issue), NOT an
implement-from-scratch. If asked to "do #239," the work is verification + issue hygiene.

## Current wiring (all scheme-generic, no inline diamond)

- `orpheus/sn/spatial/scan.py:312 _scanmarch_row` — scheme-generic: uses blend `w` +
  `DiscretizationSchemeBase.cell_average(in,out,w)` / `.outgoing_face_from_average(ψ̄,ψ_y,w)`.
  No `0.5`/`2.0`. Docstring `.. note::` says "Since #239 the row-march is scheme-generic."
- `loss_representation.py:1562 ScanMarch._sweep_interior` — calls
  `scheme.cartesian_scan_coefficients(s_scan, s_transverse, sig_t) -> (a, inverse_denom, w, (c_y,))`
  then `scheme.source_emission(QV + c_y·ψ_y, inverse_denom, w)`. No inline `2`-factor.
- `loss_representation.py:1629 ScanMarch._loss_action_interior` (the renamed `_scanmarch_matvec`)
  — calls `scheme.reflect_scan_coefficients(ψ̄) -> (α=−1, β=2ψ̄)` + `scheme.residual_kernel_batch(...)`.
  No inline diamond reconstruction.
- `diamond.py:634 DiamondDifference.cartesian_scan_coefficients` — the DD producer (diamond
  `2=1/w_DD` owned HERE; delegates diagonal to `_cartesian_streaming_diagonal`).
- `diamond.py:692 DiamondDifference.reflect_scan_coefficients` — `(α=−1, β=2ψ̄)` via `_reflection_coeffs(ψ̄, _DD_W)`.
- `scheme.py:1031 / :1088` base `cartesian_scan_coefficients` / `reflect_scan_coefficients` —
  ABSTRACT, raise `NotImplementedError`, gated on `transverse_coupling_is_facewise`. So the
  method IS on the polymorphic contract (the 1-D `affine_scan_coefficients` sibling lives at scheme.py:981).

## The LD-2-D question: STRUCTURAL exclusion, NOT underived research

LD sets `is_affine_scannable=True` (linear_discontinuous.py:286 — so LD rides the **1-D**
scan/`CumprodScan`) but leaves `transverse_coupling_is_facewise` at its `False` default
(scheme.py:753; LD docstring line 294). LD does NOT override `cartesian_scan_coefficients`
or `reflect_scan_coefficients` — and never will.

**Reason (linear_discontinuous.py:294-298, the load-bearing sentence):** LD's multi-D
closure is BILINEAR — an independent slope moment per axis (the cross moment ψ̂_xy is
diffusion-limit-load-bearing). The d-D transverse coupling is therefore a 1st-order SLOPE
moment, NOT a 0th-order face value. The row-march `scan(x)∘march(y)` is exact ONLY when the
transverse coupling folds into the scan source as a 0th-order face value (the
`transverse_coupling_is_facewise` precondition). LD's bilinear closure does NOT separate into
per-axis 1-D scans, so LD rides the DAG wavefront in d≥2 (#38 / #240 D5b), by design.

This is NOT a "the LD-2-D scan math is unpublished" exclusion (unlike, say, curvilinear-LD
which raises in `affine_scan_coefficients`). It is the honest #236-style selection: the SEAM
is polymorphic, DD is the only occupant, LD-2-D-scan is a permanent architectural exclusion
because the math structurally can't ride a separable row-march. The `supports` gate
(loss_representation.py:1443) already enforces this: `is_affine_scannable` for the 1-D arm,
the DISTINCT `transverse_coupling_is_facewise` for the 2-D arm.

## The bit-id / equivalence gate (already present)

`tests/sn/sweep/cartesian_2d/test_scan_march_equivalence.py` — `ScanMarch` ≡
`FullFieldWavefront` oracle (principled-equiv, rtol 1e-11 / atol 1e-12; NON-SQUARE x↔y moat;
vacuum+reflective; foundation; verifies `loss-rep-scanmarch{,-solve,-apply}`). This is the
DD principled-equivalence pin through the lift. The issue's named DD-byte-identity DriftWarning
gate + the `w≠1/2` 2-D two-paths gate would only be NEEDED if a facewise non-DD scheme (Step)
were added — none exists today.
