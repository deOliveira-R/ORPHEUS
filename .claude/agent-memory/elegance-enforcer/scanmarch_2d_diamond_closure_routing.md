---
name: scanmarch-2d-diamond-closure-routing
description: #206 ScanMarch-2D diamond-closure routing — the A1 follow-up; routes the surviving 2-D row-march FACE-closure twin through the cell_update seam; PASS-with-nits
metadata:
  type: project
---

# #206 ScanMarch-2D — the 2-D row-march FACE closure single-sourced into the `cell_update` seam

PASS-WITH-NITS (no commit-blockers). The exact follow-up [[a1_diamond_face_closure_seam]]'s CONCERN
section named: route the surviving 2-D ScanMarch row-march twin (`loss_representation.py:1422` `out_y`)
into the `outgoing_face_from_average` seam, leaving the `β=2ψ̄` scan COEFFICIENT at `:1410` as-is.
Implementation is faithful to that recommendation. Working-tree diff; comprehensive 2-D bit-id gate
ran separately (per dispatch constraint, I did NOT run pytest).

**The 3 routed sites (verified exhaustive for the clean ScanMarch-2D closure):**
1. `scan.py:_scanmarch_row` (the SWEEP row): `psi_avg = 0.5*(in_x+out_x)` → `cell_update.cell_average_from_faces(in_x,out_x)`; `out_y = 2.0*psi_avg − psi_y_in` → `cell_update.outgoing_face_from_average(psi_avg, psi_y_in)`. Gains a `cell_update: CellUpdateBase` param.
2. `loss_representation.py:1423` `_loss_action_interior` (the MATVEC row): `out_y = 2.0*psi_bar_row − psi_y_in` → `self.mesh.cell_update.outgoing_face_from_average(psi_bar_row, psi_y_in)`.
3. `loss_representation.py:1315-1318` (the SWEEP caller `_loss_scan_interior`): passes `self.mesh.cell_update` to `_scanmarch_row`.

**EXHAUSTIVENESS PROOF (grep `2\.0 \*` / `0\.5 \*` both files):** `scan.py` now has ZERO inline closure
arithmetic. `loss_representation.py`'s surviving `2.0 *` are ALL scan COEFFICIENTS, NOT closures: `:1310`
`α=2s_x/D−1`, `:1312` `β=2(Q+s_yψ)/D` (sweep x-scan solve coeffs), `:1411` `β=2ψ̄` (matvec apply coeff),
`:1997`/`:2159` `b=2·QV·inv_denom` (CumprodScan chain β). NONE is `2ψ̄−ψ_in` or `½(in+out)`. The diamond-
outgoing-face/cell-average CONCEPT now lives in ONE place (the DD seam `diamond.py:467-497`). The A1
THREE-sites-of-the-concept count is back to ONE. Institutional-pattern-#2 "third site appeared → stop"
is RESOLVED.

**AXIS-1 verdict — leaving `β=2ψ̄` at `:1410/:1411` is CORRECT (PASS).** It is NOT a `ψ_out = 2ψ̄ − ψ_in`
closure; it is the scan-recurrence COEFFICIENT `β` of the apply-direction face recurrence
`out_x = α·in_x + β` with `α=−1, β=2ψ̄` (`_x_scan_faces` docstring `scan.py:264-269` names it explicitly:
"since ψ̄ is KNOWN, `out_x=2ψ̄−in_x` IS a first-order recurrence in the faces"). The `2ψ̄` here is the
SAME diamond relation algebraically, but in its SCAN-COEFFICIENT encoding — not a clean closure call. The
seam takes `(cell_avg, face_in)` and returns one face; `_x_scan_faces` takes scalar β broadcast over a
whole x-chain and ordinate-scans. Routing `:1411` through the seam would be a TYPE/SHAPE mismatch (the seam
is a single-cell relation, the coeff feeds a chain solver). Forcing it would be procedural-transcription-in-
reverse. The honest single source for the `2ψ̄` coefficient identity is the `_x_scan_faces` docstring's
explicit cross-reference, which IS present. ACCEPT.

**AXIS-2 verdict — `scan.py` → `cell_update` TYPE_CHECKING import + method calls is ACCEPTABLE (PASS).**
Both modules live in the `orpheus.sn.spatial` package; the diamond closure is genuinely a scheme/
discretization concern, and `_scanmarch_row` is a scan-substrate primitive that must close the cell it
scans. The dependency direction is scan-primitive → closure-protocol, which is the correct layering
(the scan consumes the closure, not vice-versa). `TYPE_CHECKING`-only import keeps it annotation-only (no
runtime cycle); `cell_update` is passed as a parameter (dependency injection, not a module-level import of
a concrete class). This MIRRORS the established `loss_representation` pattern (`self.mesh.cell_update.*`).
No layering smell.

**AXIS-3 verdict — `_scanmarch_row` should NOT dissolve now; defer to Phase B is CORRECT.** The dispatch
framing "sweep and matvec are nearly symmetric" is true at the ENDPOINTS (`_x_scan_faces` +
`outgoing_face_from_average`) but FALSE in the middle, and the middle is the point:
  - SWEEP (`_scanmarch_row`): ψ̄ is UNKNOWN. Recovers it via `cell_average_from_faces(in_x, out_x)` from
    the SOLVED faces (needs BOTH faces), then emits `out_y`. Uses BOTH seam methods.
  - MATVEC (`_loss_action_interior`): ψ̄ is KNOWN (`probe_oct`). Calls `_x_scan_faces` only for `in_x_row`
    (DISCARDS `_out_x_row`), then builds the residual `LpC_oct = D·ψ̄ − s_x·in_x − s_y·ψ_y_in` (NOT a cell
    average), then emits `out_y`. Uses ONE seam method.
The arms share `_x_scan_faces` (already single-sourced, Pattern-2 satisfied `scan.py:259-269`) and
`outgoing_face_from_average` (now single-sourced). What differs — recover-ψ̄-then-emit-angular-flux vs
residual-from-known-ψ̄ — is the apply/solve duality itself, the SAME asymmetry the C3.x/S6.4(e)
`_CellSolve`/`_CellResidual` split institutionalized at the level-op layer. Dissolving `_scanmarch_row`
into the sweep interior now would abstract over the DIFFERENCE (the unify-after-two trap, anti-pattern #10
+ institutional tell "abstracting over the difference"). The shared-frame extraction is precisely the
declared Phase B; this carve correctly stops at the FACE-closure seam (the genuine shared surface revealed)
and does not over-reach. DEFER endorsed.

**AXIS-4 — no elegance regression from the mechanical edit (PASS).**
  - `_scanmarch_row` docstring REWRITTEN to name the seam + the #206-A single-source rationale + both
    `:meth:` cross-refs (`scan.py:328-333`) — docstring accuracy IMPROVED, not just preserved.
  - `_x_scan_faces` docstring (`scan.py:262-269`) already carried the apply/solve coefficient cross-ref —
    untouched, still correct.
  - Line-wraps clean; the param added trailing on the caller (`loss_representation.py:1317`) and in the
    signature (`scan.py:321`). No lost comments. The `:1423` `out_y` 2-line wrap is idiomatic.
  - Seam docstrings (`cell_update.py:730-750`, `diamond.py:482-497`) already enumerate ScanMarch's
    transverse-y march as a consumer — they predicted this consumer; now it is real.

**NITS (non-blocking, carried from A1 / pre-existing):**
1. (A1-carried, NOT this carve) the `psi_face_in = outgoing_face_from_average(psi_cell, psi_face_in)` rebind
   in `operator.py` 1-D matvec — pre-existing, A2 owns the `psi_face_out` rename. Not in this diff.
2. `_loss_action_interior:1410` discards `_out_x_row` (underscore-named, honest) — the matvec needs only
   the upstream face. Correct, not a smell; noted only so a future Phase-B reader knows the discard is
   deliberate, not a dropped value.

**The β=2ψ̄ removal-trigger A1 asked for:** A1 recommended "decide if β=2ψ̄ wants a named cell_update-owned
coeff helper or a reciprocal cross-ref comment with a removal trigger." This carve chose the cross-ref
route (the `_x_scan_faces` docstring names the coefficient identity) and did NOT mint a coeff helper —
CORRECT (a single-cell `outgoing_face_from_average` cannot type the chain-broadcast coefficient; a separate
`affine_scan_coefficients`-style helper would be unify-after-one). No dangling anti-pattern-#11 TODO.

Files: `scan.py` (251-312 `_x_scan_faces` unchanged single-source; 315-342 `_scanmarch_row` +cell_update
param +routed +docstring; TYPE_CHECKING import 77-82), `loss_representation.py` (1315-1318 sweep caller
passes cell_update; 1423-1425 matvec out_y routed; 1410-1411 β=2ψ̄ coeff LEFT correct), `cell_update.py`
(710-750 Protocol base), `diamond.py` (467-497 DD seam).

Related: [[a1_diamond_face_closure_seam]] (the 1-D precursor + the CONCERN this carve closes),
[[s6_4_e_walk_levelop_collapse]] (the _CellSolve/_CellResidual apply/solve split this row-march mirrors),
[[sweep_strategy_carve]] (CumprodScan/ScanMarch/Compatibility selection layer).
