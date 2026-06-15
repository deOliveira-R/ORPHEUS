---
name: issue-158-increment-b-coefficient-model
description: #158 Increment B coefficient-model carve (affine_scan_coefficients→(a,inverse_denom,w); affine_closure.py; closure-method hooks retired) — PASS-WITH-NITS, the doc/reality conflicts + the apply-formula overstatement to re-check at commit/Inc-C
metadata:
  type: project
---

# #158 Increment B — the coefficient model (affine spatial closure)

Surgical main-agent-direct carve on `feature/sn-space-angle-tier2`. Reframed group-3
of the `CellUpdate` contract: a consistent affine scheme is characterized per cell by
THREE σ_t-epoch coefficients `(a, inverse_denom, w)`; representations (scan / matvec /
future explicit-matrix) apply them with their own schedule. Built on Increment A (the
DAG-family group-2 kernel). Verdict PASS-WITH-NITS; architecture sound.

## What landed (the elegant core — these are RIGHT, do not re-litigate)
- `affine_scan_coefficients` 2-tuple→3-tuple `(a, inverse_denom, w)`; the per-scheme
  closure-method HOOKS `cell_average_from_faces`/`outgoing_face_from_average` DELETED
  from `cell_update.py` base + `diamond.py`. Group 3 = pure coefficients, NO cell math.
- NEW `affine_closure.py`: two generic free funcs `source_emission(QV,inv,w)=QV·inv/w`
  + `cell_average(in,out,w)=(1−w)in+w·out`. Free-functions-not-value-object = CORRECT
  (stateless pointwise ops, no invariant to protect, both consumers import directly;
  NOT exported from spatial/__init__ — internal helpers, fine).
- LD flipped `is_affine_scannable=True` (Schur eliminates slope → single-upstream
  recurrence). LD `affine_scan_coefficients` = ×V Schur `(a,inv,w)` + slab-only guard
  `if any(dA_w!=0) or any(c_out!=0): raise NotImplementedError`. Guard = Pattern-4
  fail-fast at CollisionCache build (geometry-as-data); GOOD.
- `CollisionCache` gained `cell_average_weight` field. DD returns `w=½` broadcast.

## RULINGS (settled — do not re-litigate)
- **Coefficient-model coherence = HONEST.** Scheme owns ONLY `(a,inv,w)`; no cell math
  leaked into any sweep body (grep-verified: scan.py + _run + _apply_walk carry only
  inline DD-w=½ constants in the explicitly-DD-only paths). Pattern-2 SSOT realized.
- **The "matvec=group-2, scan-solve=group-3" SPLIT = CORRECT, not a smell.** Genuine
  apply/solve asymmetry: solve has NO concrete inflow until prefix-scan runs (needs
  factored coeffs); apply has a concrete probe ψ̄ (rides residual_kernel_batch). This
  is the _CellSolve/_CellResidual asymmetry again (cf. scanmarch-2d memory). Folding =
  unify-over-the-difference. The unit gate `test_group3_equals_group2_scan_flat` pins
  the SOLVE consistency at nULP. KEEP.
- **`A_total` unused on LD = contract uniformity, NOT a signature smell.** `affine_scan_
  coefficients` is a geometry-BLIND shared contract (Pattern 7 producer-side, neutral
  slab fields via slab_streaming). DD USES all params (A_down:75, A_total:85, dA_w/c_out:76
  — DD is all-geometry); LD uses A_down + the guard fields, skips A_total (slab: A_total=
  2·A_down redundant). Same geometry-blind-by-data ruling as Increment A. KEEP.
- **2-D ScanMarch inline (sweep scan.py:_scanmarch_row + matvec :1422) = JUSTIFIED
  BOUNDED SEAM.** 2-D is DD-only; the α/β already bake w=½ (the `2` factors), so the
  closure inline is CONSISTENT with the in-row recurrence. The cell_update kwarg was
  dropped from _scanmarch_row (TYPE_CHECKING import gone — clean). Deferred 2-D
  coefficient-model lift documented in the scan.py .note. SHOULD be a filed follow-up
  (user said they plan to file — verify the issue exists at commit).
- **Curvilinear scan SOLVE arm (_run :2766) uses coefficient model and is DD-ONLY in
  practice** (LD curvilinear raises at cache build). DD provides w=½ curvilinearly, so
  source_emission/cell_average apply correctly. The asymmetry "LD guards slab-only, DD
  supports all geometries" is PRINCIPLED (LD's curvilinear closure unpublished).

## NITS / FINDINGS (approval conditions — concrete)
- **NIT-1 (CONCERN, do-now — anti-pattern #11 + Smell-#7, overstated identity):** the
  `affine_closure.py` MODULE docstring (lines 41-45) asserts a precise APPLY formula
  `residual = (ψ̄ − (1−w+w·a)ψ_in)/inverse_denom` and claims it is "verified consistent
  ... by the group-2 ≡ group-3 cross-check," then builds the future-ExplicitMatrix
  guidance (diagonal `1/inverse_denom`, upstream coupling `(1−w+w·a)/inverse_denom`) on
  it. VERIFIED numerically: this formula does NOT reproduce DD's `residual_kernel_batch`
  (DD group-2 disagrees; the `(1−w+w·a)` ψ̄-coupling IS right for both DD+LD, but the
  `/inverse_denom` conflates the ×V SCAN inverse_denom with the ÷V KERNEL diagonal — they
  differ by V and by the source-fold). The ACTUAL test pins the SOLVE direction
  (psi_avg/psi_out), NOT this apply-residual formula. Bug habitat: a future ExplicitMatrix
  implementer transcribes the wrong DD matrix coupling. FIX = correct the formula to the
  ÷V kernel convention (or drop the explicit formula and point at residual_kernel_batch as
  the SSOT) + scope the "verified" claim to what the test actually pins (solve direction).
- **NIT-2 (CONCERN, do-now — Cardinal Rule 3, doc/reality conflict):** `affine_closure.py`
  docstring (lines 47-53 .note) says "The DD regression snapshots RE-BASELINE accordingly."
  FALSE — the gate `tests/sn/{sweep/core,solve}` with DriftWarning-as-error is GREEN at
  505/1skip/4xfail with NO re-baseline (DD byte-identical: w=½ makes generic ops exact
  power-of-2 scalings; x/0.5==x*2.0 bit-exact; residual reassoc bit-id on pinned power-of-2
  meshes). The principled-equiv caveat is RIGHT for LD (two-paths gate explicitly nULP),
  WRONG for DD. FIX = scope the re-baseline note to LD / future non-w=½ schemes; state DD
  stayed byte-identical.
- **NIT-3 (CONCERN, do-now — anti-pattern #11, stale `:meth:` to RETIRED methods, will
  break Sphinx xref + states a now-FALSE fact):**
  - `diamond.py:133` — `is_affine_scannable` docstring cites `:meth:cell_average_from_faces`
    + `:meth:outgoing_face_from_average` (both DELETED) AND says "A Linear-Discontinuous
    closure ... is therefore NOT affine-scannable" — FALSE as of THIS increment. Fix both.
  - `diamond.py:485` — class structural comment calls group 3 "the scan-family TRIPLE
    (`affine_scan_coefficients`/`cell_average_from_faces`/`outgoing_face_from_average`)" —
    now a SINGLE coefficient method. Update to the coefficient-model narrative (headline of
    the increment, yet the class's own overview still describes the retired triple).
  - `loss_representation.py:1886` — `_apply_walk` docstring references
    `cell_update.outgoing_face_from_average` + "reconstruct via the already-seamed diamond
    closure" — now the Cartesian arm rides residual_kernel_batch, curvilinear inlines. Stale.
  - DISCRIMINATOR applied: these are LIVE-role refs on current code = FIX. The 4th hit
    `tests/sn/sweep/core/test_affine_carve_baseline.py:19` names the retired methods in a
    HISTORICAL #206-Phase-A narrative but the gate is still live → softer; update since the
    gate is current, but lowest priority.
- **NIT-4 (record only — algebra-of-record gap):** LD docstring + class claim
  `(a,inv,w)` "SymPy-verified" / "proven numerically to machine ε" but NO checked-in SymPy
  artifact in `derivations/` (searched). Plan `mellow-swinging-breeze.md` may hold it;
  Cardinal Rule 3 wants the canonical derivation in the tree/theory page at close-out.
  Dispatch archivist at Inc-B close (feature not DONE until Sphinx thorough).

## Carry-forward to Increment C (Q̂≠0 slope source)
- The group1≡group2 / group3≡group2 pins are ALL at the Q̂=0 flat slice (the slope-SOURCE
  path — LM-1989 §1.4/§6 sign trap — is UNTESTED). My Increment-A NIT-1 (comment overstates
  the pin at linear_discontinuous.py:360) is STILL LIVE and now ALSO applies to the new
  `affine_closure` cross-check claim. Inc-C must add a Q̂≠0 equivalence pin.

## Gates run this review (all GREEN)
- DD regression: `tests/sn/{sweep/core,solve}` -W DriftWarning-error → 505 passed/1skip/4xfail.
- LD: `test_linear_discontinuous.py` + `test_mms_ld_slab.py` → 24 passed/1xfail.
- Curvilinear DD matvec (inline arm): cylinder+sphere matvec → 6 passed/27 xfailed (the 27
  = pre-existing #206 hand-ref divergence, untouched).
