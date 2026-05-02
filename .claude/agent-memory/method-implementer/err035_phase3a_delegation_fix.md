---
name: ERR-035 Phase-3A delegation fix closeout
description: 2026-05-02 Phase-3A symmetric slab heuristic closure refactored as thin wrapper over Phase-3B rank-2 path. Bit-equal agreement at intermediate α; vacuum convergence floor improved from 5e-4 to 9e-6; full 116-test gate green; sphere/cylinder bit-equality preserved.
type: project
---

# ERR-035 Phase-3A symmetric slab → Phase-3B rank-2 delegation fix

**Branch**: `feature/peierls-greens-cylinder` (from commit `43625bd`)
**Date**: 2026-05-02 (same day as ERR-035 was discovered + documented)
**Verdict**: Shipped clean. The recommended-approach delegation
worked on the first try; no structural blockers required the
in-place rewrite alternative.

## The user's contract (verbatim)

> Fix ERR-035 in the Phase 3A symmetric slab Variant α path. The
> bug: Phase 3A's `_apply_operator_slab` uses a heuristic closure
> `ψ_surf = α·B_period/(1−α²·e^{−2τ})` that's mathematically wrong
> at intermediate α (0 < α < 1). The first-principles rank-2
> reduction at α_L = α_R = α gives `ψ_surf = α·B_single/(1−α·e^{−τ})`.

> The cleanest fix: delegate Phase 3A symmetric slab to the Phase 3B
> rank-2 path with α_L = α_R = α.

The delegation approach worked cleanly. Phase-3A's outer iteration
scaffolding did NOT need separate rewriting — it was entirely
absorbed by the asymmetric solver's iteration scaffolding.

## Files modified

| File | Status | Net lines | Purpose |
|------|--------|-----------|---------|
| `orpheus/derivations/continuous/peierls_greens_function/greens_function_slab.py` | rewritten | −495 / +245 | Public API (`SlabGreensResult`, `SlabGreensMGResult`, `solve_greens_function_slab{,_mg}`) preserved as re-wrappers; deleted `_apply_operator_slab`, `_first_leg_chord_slab`, `_bounce_period_chord_slab`, all `scipy.interpolate` and quadrature imports |
| `tests/derivations/test_peierls_greens_function_slab_asymmetric_solver.py` | modified | +20 / −37 | Repurposed `test_rank2_vs_rank1_at_intermediate_alpha_documented_discrepancy` → `test_rank1_path_now_agrees_with_rank2_via_delegation_after_ERR035_fix` (≤ 1e-12 rtol bit-equal gate, tagged `@pytest.mark.catches("ERR-035")`); updated `test_rank2_symmetric_BC_agrees_with_phase3a_at_endpoints` docstring; updated test-strategy header |
| `tests/derivations/test_peierls_greens_function_slab_solver.py` | modified | +37 / −22 | `test_alpha_zero_convergence_floor` re-pinned with new ~9e-6 floor and 5e-5 gate (was 5e-4 / 1e-3); docstring updated with ERR-034 + ERR-035 history |
| `docs/theory/peierls_greens.rst` | modified | +35 / −36 | Updated Phase-3A intro to reference delegation; rewrote `peierls-greens-slab-T` equation as first-principles single-transit form with cross-link to rank-2 monodromy; rewrote "Phase-3A rank-1 algebraic discrepancy" section as `_peierls-greens-slab-asym-err035-fix` with FIXED status; updated convergence-floor narrative; cleaned `alpha_per_period` extension narrative |
| `.claude/skills/vv-principles/error_catalog.md` | modified | +33 / −14 | ERR-035 status flipped from "DOCUMENTED, fix deferred" to "FIXED 2026-05-02", with full fix description + side-effects (vacuum-floor improvement, sphere/cylinder bit-equality preservation) |

**Net code change**: −282 production lines (Phase-3A's heuristic
machinery was ~250 lines of dead-code-after-delegation; the new
wrapper is ~75 lines).

## Test results — full 116-test gate

```
pytest tests/derivations/test_peierls_greens_function_solver.py \
       tests/derivations/test_peierls_greens_function_vacuum.py \
       tests/derivations/test_peierls_greens_function_xverif.py \
       tests/derivations/test_peierls_greens_function_xverif_ps1982.py \
       tests/derivations/test_peierls_greens_function_mg.py \
       tests/derivations/test_peierls_greens_function_symbolic.py \
       tests/derivations/test_peierls_greens_function_cylinder_symbolic.py \
       tests/derivations/test_peierls_greens_function_cylinder_solver.py \
       tests/derivations/test_peierls_variant_alpha_core.py \
       tests/derivations/test_peierls_greens_function_slab_symbolic.py \
       tests/derivations/test_peierls_greens_function_slab_solver.py \
       tests/derivations/test_peierls_greens_function_slab_asymmetric_symbolic.py \
       tests/derivations/test_peierls_greens_function_slab_asymmetric_solver.py
======================= 116 passed in 202.66s (0:03:22) ========================
```

Pre-fix baseline: 116 passed. Post-fix: **116 passed**. Test count
unchanged (one test repurposed). Sphinx `-W` build clean (incremental;
matches pre-fix state — the same 3 SN-MMS unmatched-equation
warnings persist, all unrelated to peierls).

## Accuracy results — pre vs post fix

| Test / regime | Pre-fix floor | Post-fix floor | Status |
|---------------|--------------:|---------------:|--------|
| Sphere V_α1 numerical k_eff = k_inf | 8.3e-17 | **8.3e-17** | bit-equal preserved |
| Cylinder V_α1 numerical k_eff = k_inf | 1.0e-15 | **1.0e-15** | bit-equal preserved |
| Slab V_α1 (α=1, const init) k_eff = k_inf | 1 iter, machine | **1 iter, machine** | preserved |
| Slab V_α1 (α=1, nonuniform init) k_eff = k_inf | < 1e-9 | **< 1e-9** | preserved |
| Slab vacuum α=0, k_eff at finest grid (32,56,128) | 1.30e-1 (with ERR-034) → 1.57e-1 (post-ERR-034) | **1.5744e-1** | k_eff/k_inf ≈ 0.756 |
| Slab vacuum α=0, finest-pair self-consistency | 5e-4 (gate 1e-3) | **8.85e-6** (gate 5e-5) | improved 56× |
| Slab α=0.5 rank-1 vs rank-2 (intermediate α gate) | 1.3e-4 (band [5e-5, 5e-4]) | **0.0e+0** (rtol ≤ 1e-12 bit-equal) | **fixed** — discrepancy eliminated |
| Slab MG 2G asymmetric scattering at α=1 | < 1e-9 | **< 1e-9** | preserved |
| Slab MG G=1 reduction at α=0 | bit-equal | **bit-equal** | preserved |
| Method-of-images symmetry, finest grid (64,64,128) | 6.5e-9 | **6.5e-9** | preserved (asymmetric path untouched) |

Key result: **the Phase-3A symmetric path now produces results
bit-equal to the Phase-3B rank-2 path at all α** because there is
only one underlying computation. The intermediate-α discrepancy
that ERR-035 documented (~1.3e-4 at α=0.5) is eliminated.

## Vacuum convergence floor — re-measured

```
(16, 24, 48): k_eff = 1.5743600374e-01
(20, 32, 64): k_eff = 1.5743967793e-01
(24, 40, 96): k_eff = 1.5744130513e-01
(32, 56, 128): k_eff = 1.5744269900e-01

Refinement diffs:
  (16, 24, 48) -> (20, 32, 64): rel diff = 2.334e-05
  (20, 32, 64) -> (24, 40, 96): rel diff = 1.034e-05
  (24, 40, 96) -> (32, 56, 128): rel diff = 8.853e-06
```

The pre-fix Phase-3A closeout misattributed the slow ~5e-4 floor
to a Gauss-Legendre-on-µ quadrature limitation, predicting that
"a future improvement could use a µ-weighted Gauss-Jacobi
quadrature concentrating nodes near µ = 0". With ERR-034 + ERR-035
fixed, the actual floor is ~9e-6 between the two finest grids, ~56×
tighter. The slow convergence was not a quadrature limitation; it
was the dominant ERR-034 trajectory bug masquerading as quadrature
noise. The GJ-on-µ refinement is no longer load-bearing.

## Repurposed agreement test — tolerance achieved

The original `test_rank2_vs_rank1_at_intermediate_alpha_documented_discrepancy`
asserted a `[5e-5, 5e-4]` disagreement gate at α=0.5 (fail if
agreement gets tighter — silent fix — or looser — regression).
Post-delegation, the achieved agreement is **0.0e+00 relative**
(bit-equal) because both code paths now perform the SAME
underlying computation.

The repurposed test
`test_rank1_path_now_agrees_with_rank2_via_delegation_after_ERR035_fix`
uses `np.testing.assert_allclose(rtol=1e-12, atol=0.0)` — anything
tighter than rtol=1e-12 would be testing IEEE-754 bit-equality,
which is not portably guaranteed even for the same computation
under different reduction orders (although in practice it holds
here because the asymmetric solver is called identically).

The test is now tagged `@pytest.mark.catches("ERR-035")` per the
V&V error-catalog discipline.

## Honest verdict on the delegation

**Did the delegation succeed cleanly, or did Phase-3A's outer
iteration scaffolding require modification?** Cleanly. The
Phase-3A scaffolding (power iteration, fission-rate normalisation,
k_eff update, convergence test) was bit-by-bit identical to the
asymmetric solver's scaffolding — no surprise, since the
asymmetric solver was forked from Phase-3A in the Phase-3B
prototype. The delegation reduces the outer scaffolding to ~25
lines (build kwargs dict; call asymmetric solver; re-wrap result
into the back-compat dataclass).

The only structural concern was the dataclass shape: Phase-3A
exposes `SlabGreensResult` with field `psi`, while the asymmetric
solver returns `SlabAsymmetricGreensResult` with the same field
names. Re-wrapping is a one-liner per dataclass.

**No quadrature node mismatch surfaced**. Phase-3A and Phase-3B
both use the same `(n_x, n_mu, n_traj_quad)` triple with identical
`np.polynomial.legendre.leggauss` ordering, so the delegation is
literally one-to-one.

## Side-effect: code consolidation

Beyond the bug fix, the delegation removes a parallel-path
maintenance burden:

- **Pre-fix**: 695 lines in `greens_function_slab.py` with its
  own `_apply_operator_slab`, `_first_leg_chord_slab`,
  `_bounce_period_chord_slab`, `_scalar_flux_from_psi`, and 1G+MG
  power-iteration scaffolding. Risk: any future fix to the slab
  trajectory machinery had to be applied twice (in Phase-3A and
  Phase-3B). ERR-034 was such a case (fixed in both modules
  symmetrically; if the symmetry had drifted, the bug would
  resurface in one of them).
- **Post-fix**: 250 lines in `greens_function_slab.py` of which
  ~75 are wrapper code and ~175 are docstring. Single underlying
  trajectory implementation in `greens_function_slab_asymmetric.py`.

This collapses the "two-modules, one-API" Phase-3B scope guard
("DO NOT modify Phase 3A symmetric slab") that was load-bearing
during Phase 3B integration but counter-productive afterwards.

## Lessons for the agent's memory

1. **Heuristic closure formulas are an anti-pattern.** ERR-035 was
   the third instance in this project where an algebraic formula
   was generalised across geometries by analogy without
   first-principles re-derivation, and the fourth instance where
   "two closures agree at the corners" was mistaken for "they
   agree generally". The systematic remedy is what `algebra-of-record`
   already prescribes: every Variant α geometry extension MUST
   ship with a Branch-1 SymPy first-principles derivation on a
   non-uniform source, with cross-check against the rank-2
   structurally-correct closure.
2. **When the fix is "delegate to a sibling that's already
   correct", the closeout is always "delete the dead code
   completely".** The brief explicitly called this out:
   > Don't leave them as "back-compat shims" or "alternative
   > paths" — they're known wrong.
   Following that discipline collapsed 250 lines of risk into
   ~75 lines of trivial wrapping.
3. **Convergence-floor narratives are diagnostic-friendly when
   they're tight, diagnostic-hostile when they're loose.** The
   pre-fix Phase-3A vacuum floor was 5e-4 with a 1e-3 gate; the
   closeout memo predicted "GJ-on-µ would help". This was the
   classic ERR-034 sign-pattern (slow apparent convergence from a
   sub-leading linear bias). The post-fix floor is 9e-6 with a
   5e-5 gate — much sharper, and now genuinely diagnostic of the
   underlying numerics rather than a hidden bug. **A quadrature
   convergence rate that doesn't match theoretical expectations is
   a bug fingerprint, not a quadrature limitation, until proven
   otherwise.**
4. **Per the user's contract, the auto-mode bias is to commit and
   move forward, not to pause for review.** The ERR-035 fix was
   the prerequisite for Phase 3C (hollow sphere + annulus) per
   the user's explicit prioritisation. Shipping cleanly without
   user-confirmation cycles is the load-bearing behaviour.

## Status

Branch state: clean. 116/116 tests pass. Sphinx -W clean
(incremental; matches pre-fix). Sphere/cylinder/symmetric-slab
all preserved at machine floor. ERR-035 entry flipped to FIXED
status with full forensic detail. Repurposed regression-prevention
test tagged `@pytest.mark.catches("ERR-035")`. Ready to commit
+ proceed to Phase 3C.

## Recommended next step

Per the Phase-3B closeout memo's recommendation (which has now been
realised by this fix): **proceed to Phase 3C (hollow sphere +
annulus)** on the clean foundation. The rank-2 framework is the
predicted generalisation per the cross-domain frame memo, and the
clean foundation is now structurally ready to receive it.
