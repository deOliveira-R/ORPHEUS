---
name: a1-diamond-face-closure-seam
description: #206 Phase A1 — route 1-D diamond FACE closure through cell_update seam; PASS-with-nits; the surviving 2-D ScanMarch row-march twin
metadata:
  type: project
---

# #206 Phase A1 — 1-D diamond face-closure single-sourced into the `cell_update` seam

PASS-WITH-NITS (no commit-blockers). The 1-D analogue of the landed 2-D `_OctantWalk`
single-sourcing (S6.3b `3a79ab3`). Working-tree carve, bit-identity-faithful (recipe 1/2/3 +
A-NEW gates ran separately).

**What A1 did (4 hunks, verified exhaustive in scope):** the diamond FACE closure now lives in
ONE place — `DiamondDifference.cell_average_from_faces` (½(in+out)) + `outgoing_face_from_average`
(2ψ̄−in), declared on the `CellUpdate` Protocol / `CellUpdateBase` (`orpheus/sn/spatial/diamond.py:467-497`,
`cell_update.py:710-750`). Routed: sweep `½(in+out)` at `loss_representation.py:2008,2173`; matvec
`2ψ̄−in` at `operator.py:469,962` (the `_compute_LpC`/`_compute_decomposition` byte-twins, edited
IDENTICALLY — twin not deepened, Phase-C collapse left alone, correct Pattern-6 defer).
`grep "0\.5 \*"` = ZERO across both files → the two sweep sites were the only ones, both routed.

**The `supports` gate (Pattern-4 win, correct invariant + correct layer):** `CumprodScan.supports`/
`ScanMarch.supports` now gate on `mesh.cell_update.is_affine_scannable` (scan schedules consume the
closed-form affine recurrence `ψ_out=a·ψ_in+b` — only exists for affine closures; LD couples two
face moments → `False`). Boolean composition CORRECT: `(is_1d or (is_cartesian and ndim==2)) and
is_affine_scannable` — scannability `and`-conjoined at the OUTER level (property of the closure,
orthogonal to geometry; a non-affine closure disqualifies every geometry). Defense-in-depth NOT
redundancy: gate = frontend-checkable `Compatibility.reason`; the scan-family-triple
`NotImplementedError` defaults (`cell_update.py:703,725,747`) = runtime backstop; `IncompatibleRepresentation`
construction guard catches `default_for` bypass. No helper-extraction warranted (unify-after-one
for two 2-line predicates).

**THE FINDING (CONCERN, out of A1 scope, follow-up):** a surviving FACE-closure twin in the
ScanMarch 2-D row-march — `loss_representation.py:1410` (`2.0*psi_bar_row` as scan coeff `β`) +
`:1422` (`out_y = 2.0*psi_bar_row - psi_y_in`), inside `_loss_action_interior`. `_x_scan_faces`'s
docstring (`scan.py:262-264`) LITERALLY names the apply form `out_x = 2ψ̄ − in_x` → identical
closure, confirmed not a different algebraic form. So the diamond-outgoing-face CONCEPT now lives in
THREE places (the seam + these two un-routed 2-D sites) = institutional pattern #2's "third site
appears → stops being acceptable." NOT VIOLATION/NOT blocking because: (a) explicitly outside A1's
1-D scope; (b) A1 touched none of these lines (only ScanMarch hunk = the supports gate); (c)
structural wrinkle — `:1410`'s `2ψ̄` is the scan COEFFICIENT `β` (half the recurrence), not a clean
`ψ_out`, so only `:1422`'s `out_y` is a mechanical `outgoing_face_from_average` call. Follow-up:
route `:1422` when the row-march is next under the knife; decide if `β=2ψ̄` wants a named
cell_update-owned coeff helper or a reciprocal cross-ref comment with an A-series removal trigger
(anti-pattern #11). Did NOT file an issue (A-series owns it; no-issues-for-inline-fixes).

**NIT (non-blocking, pre-existing):** `psi_face_in = outgoing_face_from_average(psi_cell, psi_face_in)`
(`operator.py:469,962`) rebinds "face_in" to its own OUTPUT (loop-carried inflow=last outflow,
legitimate idiom) — but the named method makes the semantic mismatch louder than the bare arithmetic
did. Pre-A1 identical; mechanical edit didn't cause it. Optional `psi_face_out` rename deferred to A2
(that pass is in this loop anyway; A1 must not risk the bit-id gate).

**A2 is the next step (denom/coefficient single-sourcing) — NOT in A1; the matvec still calling
`cell_balance_for_streaming` for the denom is DELIBERATELY deferred, do not flag.**

Files: `operator.py` (469,962), `loss_representation.py` (2008,2173 routed; 701-704/1216-1220 gates;
1410,1422 surviving twin), `diamond.py` (467-497 seam), `cell_update.py` (351-376/553-559/710-750
Protocol+base+is_affine_scannable), `scan.py` (246-307 `_x_scan_faces` closure-identity proof).

Related: [[s6_4_e_walk_levelop_collapse]] (the 2-D `_CellSolve`/`_CellResidual` level-op single-sourcing
this is the 1-D mirror of), [[sweep_strategy_carve]] (CumprodScan/ScanMarch/Compatibility/default_for
selection layer).
