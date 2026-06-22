---
name: issue-240-phase2-step-d1-outgoing-face
description: #240 Phase 2 D1 review — outgoing_face_from_average single-sources the per-scheme outflow reconstruction; PASS-with-nits; LD w/k algebra verified
metadata:
  type: project
---

#240 Phase 2 Step D1 (branch `feature/sn-space-angle-tier2`, UNCOMMITTED at review time): adds ONE generic op
`outgoing_face_from_average(psi_bar, face_in, w) = (ψ̄−(1−w)·face_in)/w` to `orpheus/sn/spatial/affine_closure.py`
(inverse of the existing `cell_average`), routes 7 previously-inlined per-scheme outflow reconstructions through it.

**VERDICT: PASS-with-nits.** The collapse is GENUINE single-sourcing (Pattern-2 win); grep-clean.

**Why:** D1 follows the established `affine_closure` free-func family pattern (`source_emission`/`cell_average` siblings);
the op is the literal math inverse, named in domain language, correctly LEFT as a free function for D1 (D2 homes all
three onto the DiscretizationScheme Base — premature method-izing would be Pattern-6 violation).

**How to apply (durable D-series facts I verified — trust on re-review):**
- **Collapse is complete.** Grep `2(\.0)?\s*\*\s*psi_(bar|avg)` across `orpheus/sn/` → ONLY survivor is
  `loss_representation.py:1435` `2.0*psi_bar_row` = the β SCAN-source feeding `_x_scan_faces` (CORRECTLY deferred to
  #239; the diff comment names it). The `denom*psi_bar − numer` (diamond.py:393) / `eff_denom*psi_bar − rhs`
  (lin_disc.py:488) are RESIDUAL balance terms, NOT face reconstructions — correctly excluded.
- **LD w/k algebra VERIFIED EXACT.** Reconstruction path (`_kernel_terms`, lin_disc.py:431-433): `g_over_theta=g/θ`,
  `d2=g_over_theta+sigt`, new `w=1/(1+g_over_theta/d2)`. Scan-coefficient path (`affine_scan_coefficients`:557-561):
  `k=p/d2`, `w=1/(1+k)` with `p≡g_over_theta`. SAME `k`, SAME `w` domain quantity (the cell-average blend weight) →
  the reconstruction `w` and the scan `w` are ONE concept, not two spellings. Docstring `k=(g/θ)/D₂` is exact.
- **DD w=½ byte-id is REAL** (÷0.5 = exact ×2, commutes with round-to-nearest). LD w=1/(1+k) = principled ~1-ULP
  (different reduction tree). Both correctly stated in op docstring + test file. Runtime green (strict 505, full 1083,
  3 new foundation unit tests).

**NITs (do-not-block):**
- **NIT-1 (Pattern-2, ACCEPTABLE-FOR-NOW)**: LD computes `w = 1.0/(1.0 + g_over_theta/d2)` VERBATIM at TWO sites —
  `cell_kernel_batch:463` (SOLVE) and `residual_kernel_batch:492` (APPLY). They share `_kernel_terms` for inputs but
  the `w`-line is a 2-line twin. This is the established `_CellSolve`/`_CellResidual` SOLVE/APPLY asymmetry pair (one
  solves for psi_avg, one evals residual at probe psi_bar) — legitimately separate bodies, identical w-tail. Live
  hazard: a future `w(Σ)` Péclet/κ-scheme blend (named in the D6 todo) lands on one-not-the-other. Collapse
  destination: `_kernel_terms` should RETURN `w` (it already returns `d2, g_over_theta`). STOPS being acceptable the
  moment a 3rd LD reconstruction site appears or D2/D6 touches one.
- **NIT-2 (Sphinx, ACCEPTABLE)**: the `.. todo::` stub in `discrete_ordinates.rst` is the anti-pattern-#11 EXCEPTION
  (tracked artifact: names #240 Phase 2 Step D6 + the closeout memo + the unit-gate module). Correct interim.

Carry-forward: D2 homes the three free funcs onto DiscretizationScheme Base — re-check NIT-1 collapses then.
Pin the LD reconstruction at non-flat slope-Q̂ (all D-series LD pins ride Q̂=0; the slope-source path is the recurring
LM-1989 sign-trap blind spot).
