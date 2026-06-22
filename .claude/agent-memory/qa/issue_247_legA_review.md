---
name: issue-247-legA-review
description: "#247 Leg A QA review (the Mode-10 slope-SOURCE closeout for 2-D Cartesian LD). VERDICT SUPPORTED, NO ERR, NO blocker. Branch refactor/sn-foundation-cleanup, UNCOMMITTED. 2 non-blocking doc-nits (lift docstring overstates eigenvalue sharing; LD:304 slope-label transposition)."
metadata:
  type: project
---

# #247 Leg A — slope-SOURCE half of the LM-1989 trap (Mode-10 closeout) — QA review

Branch `refactor/sn-foundation-cleanup`, UNCOMMITTED. Host `.venv/bin/python -O`.
Review of the typed-union widening of `_build_fixed_source_rhs` +
`_lift_external_source_to_moments` (`orpheus/sn/solver.py`) that threads a
moment-resolved EXTERNAL source's slope rows Q̂ through instead of zeroing
them, + the #247 test block (`tests/sn/verification/mms/test_mms_ld_2d.py`).

## OVERALL VERDICT: SUPPORTED. NO ERR-063. NO blocker. Commit OK.

All 4 claims SUPPORTED with mutation evidence. Two NON-BLOCKING doc-nits
(follow-up, not commit-gating).

## Claim-by-claim

### Claim 1 (Mode-10 gap closed) — SUPPORTED
The slope-SOURCE sign (external Q̂ AND scattering Σ_s·φ̂) is now CONSTRAINED.
Proven by the in-process mutation: with the lift regressed to re-zero slope
rows (the EXACT #247 bug), the 3 new gates RED while the EXISTING flat scalar
gate `test_ld_2d_stress_converges_second_order` stays GREEN (1 passed, 18s) —
the canonical Mode-10 asymmetry. The sign-mutation gate's red message under
the buggy lift was `|Δφ|/|φ|=0.000e+00 ≤ 1e-08 — the slope row is NOT
consumed`: flipping an already-zeroed row is a no-op, which is exactly the gap.

### Claim 2 (flat path byte-identical) — SUPPORTED
Lift probe: LD flat → slot0==flat + slopes EXACTLY zero; DD/Step (per_axis==1,
tail==()) → input returned unchanged AND `is` the same object (true byte-id).
Two-paths bit-id gates (`test_ld_2d_two_paths_ffw_equals_mfw` +
`_stress_`) PASS. Main agent's strict DriftWarning gate 520p/1s/4xf (no drift).
The validation discriminates by EXACT shape (`==expected` / `==moment_expected`)
BEFORE the lift; the moment branch is unreachable for flat/DD.

### Claim 3 (teeth structural, not the fragile value-band) — SUPPORTED
The teeth are TWO O(1) signals, NOT the spec-§0 sub-floor converged-flux band:
(a) the lift pass-through `array_equal(lifted, Qm)` at machine precision (the
production-change proof; under the bug, max |slope| dropped = 0.179 = O(1));
(b) the consumed-flip converged-flux change. Measured live (nc=24, deterministic
solve): NOISE floor (identical re-solve) = **0.000e+00**; signal y-slope 1.00e-2,
x-slope 3.01e-3, xy 5.81e-5. The smallest (xy) clears `_CONSUMPTION_TOL=1e-8`
by ~5.8e3×; noise margin is infinite (noise==0). The band is well-calibrated:
NO false-green risk. The O(h²) convergence leg is correctly labeled
NECESSARY-but-not-sufficient for the sign (a sign flip leaves order ~2).

### Claim 4 (no latent production bug papered over) — SUPPORTED, NO ERR
Independent explorer trace of the now-consumed external slope-source path
(reframe / mass / Kronecker order) found NO sign or magnitude bug:
- External + scattering moment vectors are SUMMED into ONE global-frame array
  before the sweep; the per-octant reframe (`octant_moment_frame_signs`,
  `sweep_graph.py:931-936`) is rank-gated (`is_moment_valued_by_rank`) and is a
  self-inverse involution — it re-signs external slopes EXACTLY as scattering
  slopes (no external-vs-scattering branch, no extra/missing flip, no transpose).
- Mass `R_source = M @ S_moments`, `M=diag(h,θh)` (`_ubld.py` mass_1d,
  `linear_discontinuous.py:500`) — the producer (projection) supplies BARE
  per-volume Legendre coeffs, the kernel's M adds h/θ/V (verified by the
  hand-polynomial foundation sub-gate to 1e-13).
- Kronecker order [bar, y, x, xy], x-slope=slot 2 (`moment_layout.py`,
  projector `_D2_MOMENT_SLOTS`). The slope source was UNVERIFIED before #247,
  not WRONG — the lift correctly zeroed an unverified-but-honest default
  (q̂=0 exact for a region-uniform source). NO ERR (Mode-10 proactive close).

## Specific V&V checks (all PASS)
- **L11 structural independence**: `_project_scalar_to_tensor_legendre` +
  `_manufactured_Q_pointwise` use ONLY `leggauss` + numpy primitives + hand-laid
  algebra. NEVER `_lift_external_source_to_moments`, NEVER `_ubld`/
  `LinearDiscontinuous`. The only `from orpheus` imports are `AngularSourceSink`
  + `TimedFullField` (typed source CONTAINERS, not the LD op nor the lift — no
  contamination). Foundation sub-gate reference is hand-laid polynomial algebra.
- **Anti-pattern #11 (positive AND negative)**: every mutation has both legs.
  Positive: converges O(h²) + improves-on-flat. Negative: 3 sign-flips RED +
  scattering-flip RED. Mode-10 asymmetry (flat slope row zero → flat gate blind)
  proven by running the flat gate GREEN under the buggy lift.
- **Mode 8 (-O safety)**: gates use `np.testing.*`/`pytest.fail`/`pytest.raises`
  only — NO bare assert. Ran under `-O` (the "asserts not executed" warning
  fired) and the gates still fire correctly (4 fast + 5 solve gates pass).
- **`_CONSUMPTION_TOL=1e-8`**: defensible (noise floor 0.0, smallest signal
  5.8e-5 → ~5.8e3× margin; far below the §0 trap). No false-green.

## Mutation evidence (verbatim, the teeth)
Under the buggy lift (re-zero slope rows = the #247 bug), via throwaway plugin:
```
FAILED ...::test_ld_2d_external_slope_source_threaded_through_lift - AssertionError
FAILED ...::test_ld_2d_external_slope_source_improves_on_flat - Failed: ... (7.907e-03) did NOT improve
FAILED ...::test_ld_2d_external_slope_source_sign_mutation_reddens[2-x-slope] - Failed: |Δφ|/|φ|=0.000e+00 ≤ 1e-08
3 failed
# AND under the SAME buggy lift:
test_ld_2d_stress_converges_second_order  → 1 passed   (the Mode-10 asymmetry)
```

## Non-blocking follow-ups (NOT commit-gating)
1. **Lift docstring overstates the eigenvalue sharing.** `_lift_external_
   source_to_moments` docstring + the closeout say "shared by the fixed-source
   AND eigenvalue paths". Grep shows the lift's SOLE production caller is
   `_build_fixed_source_rhs:1931`. The eigenvalue path wraps its SWEEP OUTPUT
   with `spatial_moments=per_axis` (`solver.py:450`) but does NOT call the lift.
   The shared-single-source CLAIM (Cardinal Rule 2) is true for the eq's #240
   D5b-S3 origin; the present wording is stale — trim to "the fixed-source path"
   or cite the eq-path's own wrap site. CR3/stale-doc (L-020/L-028 family).
2. **`linear_discontinuous.py:304` slope-label transposition.** Prose writes the
   d=2 cell unknown `{ψ̄, ψ̂_x, ψ̂_y, ψ̂_xy}` (x before y); canonical Kronecker
   is `[ψ̄, ψ̂_y, ψ̂_x, ψ̂_xy]` (axis0=x OUTER, axis1=y INNER → slot2=x-slope).
   PROSE only — actual slots come from `moment_layout`/`_D2_MOMENT_SLOTS`, no
   code path affected. Pre-existing (not introduced by #247), no ERR.

## Self-improvement check
Both pushbacks are documentation-accuracy (stale docstring + label transpose),
covered by existing CR3/stale-doc lessons (L-020/L-028/L-030). NO novel
anti-pattern → no vv-principles SKILL edit warranted.
