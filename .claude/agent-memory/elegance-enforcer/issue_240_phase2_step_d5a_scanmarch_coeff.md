---
name: issue-240-phase2-step-d5a-scanmarch-coeff
description: #240 Ph2 D5a — fold 2-D ScanMarch inline DD onto the scheme coefficient model; PASS-WITH-NITS (1 BLOCKING). cartesian_scan_coefficients is a legit distinct decomposition (NOT a twin) but the DD cell-balance diagonal + w=½ are now the THIRD/FOURTH transcription across 4 producers → collapse. reflect_scan_coefficients ≡ outgoing_face_from_average in (α,β) form, hard-codes w=½ instead of w-generic base.
metadata:
  type: project
---

# #240 Ph2 D5a — ScanMarch onto coefficient model — PASS-WITH-NITS (1 BLOCKING)

Branch `feature/sn-space-angle-tier2`, reviewed 2026-06-16, pre-commit.
D5a routes the 2-D Cartesian ScanMarch row-march THROUGH two new scheme
coefficient producers so the diamond `2` + blend `w` live in the scheme, not
inline. Principled ~1-ULP re-baseline (qa-verified separately). Closeout:
`.claude/agent-memory/method-implementer/issue_240_phase2_step_d5a_closeout.md`.

## The two new methods (base raising-default + DD override)
- `cartesian_scan_coefficients(*, s_scan, s_transverse, sig_t) → (a, inverse_denom, w, transverse_couplings)` — scheme.py:949 base raise / diamond.py:522 DD. The multi-D scan-march analogue of `affine_scan_coefficients`. DD: `scan_diag=2g_scan`, `c_⊥=2g_⊥`, `S=Σt+scan_diag+Σc_⊥`, `inverse_denom=1/S`, `a=2·scan_diag·inverse_denom−1`, `w=½`.
- `reflect_scan_coefficients(psi_bar) → (α, β)` — scheme.py:1013 base raise / diamond.py:579 DD (`α=−1, β=ψ̄/0.5`).

## RULING concern 1 — cartesian_scan_coefficients is NOT a twin (PASS)
The closeout's "same single source, specialised" claim is HALF right and the right half is load-bearing. The scan-march and the batch kernels are GENUINELY DIFFERENT DECOMPOSITIONS of the shared balance `S·ψ̄=QV+Σ2g_a·ψ_a`:
- `cell_kernel_batch`/`residual_kernel_batch` treat ALL axes KNOWN (DAG upstream-first) → `ψ̄=numer/denom` directly.
- `cartesian_scan_coefficients` has ONE COUPLED unknown axis (x: `in_x(i)=out_x(i−1)`) → must substitute DD closure + solve for `out_x=a·in_x+b`. The transmission `a=2(2g_x)·inv−1` is a SCAN COEFFICIENT no batch kernel computes.
This is the `_CellSolve`/`_CellResidual` asymmetry (apply has known probe ψ̄; solve has coupled inflow). `a` is inherent to the schedule. Reusing cell_kernel_batch to kill `a` is IMPOSSIBLE (batch has no coupled axis). Separate-method-not-generalize-affine_scan is also CORRECT (affine_scan is 1-D-chain-shaped, curvilinear args A_down/A_total/dA_w/c_out, no transverse slot). No common surface yet → Pattern-6 defer governs. PASS.

## NIT BLOCKING — the DD cell-balance diagonal + w=½ are now the 3rd/4th transcription
The COEFFICIENTS differ (concern 1 PASS) but the SUB-EXPRESSION `S=Σt+Σ2g_a` / `inverse_denom=1/S` / `w=½` is now spelled in FOUR producer bodies with NO shared sub-primitive:
- `2.0*s_a` left-fold byte-id between cell_kernel_batch:365 + residual_kernel_batch:407 (closeout admits "apply twin already established"); cartesian_scan_coefficients:568-574 writes a THIRD spelling.
- `w=½` hard-stamped as literal `0.5` in FIVE places (diamond.py:371,413,519,576 + `α=−1`⇔`w=½` at :592). The ONE constant that DEFINES "this scheme is DD" has no name, no home.
The `2g` fold across cell_kernel_batch/residual_kernel_batch PREDATES D5a (was Pattern-2 ACCEPTABLE-FOR-NOW). **D5a adds the THIRD spelling → crosses the institutional "third fold appears" tripwire** → fix belongs to THIS commit.
Bug habitat: the live `w(Σ)` Péclet-blend hazard (flagged in my own D2 memory) lands in 3-of-4 producers, passes every config where the 4th isn't on the scan-march path; bit-id gate is config-bounded.
Collapse dest (the codebase already telegraphs it — D2 predicted "_kernel_terms returns w"): base `_cartesian_streaming_diagonal(sig_t, s_axes)→(denom, inverse_denom)` (preserve explicit-left-fold order = bit-id pin) + a named `_DD_W=0.5` all five sites reference. This is the DD-side re-open of the LD w-fold NIT closed at D2.

## NIT CONCERN (do-now) — reflect_scan_coefficients duplicates outgoing_face_from_average
`out=β+α·in` with `β=ψ̄/w, α=−(1−w)/w` expands EXACTLY to `(ψ̄−(1−w)in)/w` = `outgoing_face_from_average(ψ̄,in,w)` char-for-char. The matvec genuinely needs the SPLIT `(α,β)` form (`_x_scan_faces` consumes coeffs not closure — real distinct need, justified op IN PRINCIPLE). The duplication is the override RE-DERIVING the arithmetic instead of computing from the one op.
(a) COULD be a CONCRETE w-generic BASE method `reflect_scan_coefficients(psi_bar, w): return -(1-w)/w, psi_bar/w` → DD/Step/LD inherit free. BLOCKER: there is NO scheme-level `w` accessor anywhere (`w` is always a PRODUCED coefficient, never a stored attr) — so base must take `w` param + caller must hold it; the apply path (`_loss_action_interior`) currently routes via residual_kernel_batch NOT cartesian_scan_coefficients, so it lacks `w` → that's the ONLY reason DD hard-codes ½.
(c) `w=½→α=−1` "no underflow" is a CONSEQUENCE of w=½, not independent — w-generic base produces α=−1 automatically. Doesn't justify hard-code.
Fix: make it a concrete w-generic base method (retires DD override, pre-empts Step twin) OR if apply-path w-plumbing waits for Step, route the override's ½-arithmetic through a shared `_reflection_coeffs(w)` + add tracked removal trigger (anti-pattern #11). CONCERN-not-BLOCKING only because the w-plumbing may legitimately belong to D5b/Step.

## PASSES (do not touch)
- Protocol placement (check 3): base-only raising-default guarded by `transverse_coupling_is_facewise`, NOT on `@runtime_checkable DiscretizationScheme` Protocol (declares only update/residual). EXACTLY matches affine_scan_coefficients/cell_kernel_batch/residual_kernel_batch precedent + the member-presence-trap NOTE. Consistent with my D5-0 ruling (capability data off runtime_checkable Protocol = pure cost). PASS.
- Consumer fold (check 4): ZERO inline DD survives. `_sweep_interior`:1328 / `_loss_action_interior`:1433 / `_scanmarch_row` scan.py:312 (gained `w` param) — no `2.0*`, no `D_row` division, no hard `0.5` in CODE. Reads as generic coefficient composition. The carve's HEADLINE goal MET. PASS.
- Transverse-source fold (check 5): `source_emission(Q + c_y·ψ_y, inv, w)` reads like the math — transverse-known faces ARE an effective source on the swept row (upstream-known like QV). Doesn't obscure cell balance, expresses it. PASS. (minor: a named `QV_eff` intermediate would read marginally better.)

## Arch opportunity (issue, module:sn type:improvement)
"Single-source the DD Cartesian cell-balance diagonal + blend-weight across the 4 coefficient producers." The 3 CARTESIAN producers share identical `2g` left-fold + `w=½`/`α=−1`; affine_scan_coefficients is CURVILINEAR (different fold order, own snapshots) → cannot fully merge, defer the curvilinear merge to the issue, land the 3-way Cartesian fold + w-naming in THIS commit.
