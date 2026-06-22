---
name: issue-236-phase3-st5-separability-gate
description: #236 Phase 3 ST5 separability gate — PASS-WITH-NITS; the recurring MMS _l2_1d directory-wide duplication tell
metadata:
  type: project
---

`tests/sn/verification/mms/test_space_angle_separability.py` — #236 Phase 3 (ST5),
the LAST campaign piece. NEW test-only module, 6 tests (2 Cartesian separable / 3
curvilinear-scalar gating / 1 curvilinear per-ordinate L27). VERDICT: **PASS-WITH-NITS**,
1 MUST-FIX. Ran green `6 passed in 2.93s` (`-O`, brief-mandated invocation).

**Why:** characterization gate pinning the (spatial ⊗ angular) error STRUCTURE —
Cartesian SEPARATES additively (`E≈E_space+E_angle`, mixed-2nd-diff≈0), curvilinear
GATES (`E≈max(E_space,E_angle)`, spatial rate quadrature-gated). Mirrors
`test_curvilinear_pole_cell_characterization.py` (#233). L1.

**How to apply (the durable lessons):**

⭐ **MUST-FIX = the recurring `_l2_1d` directory-wide twin.** The volume-weighted L2
body `float(np.sqrt(np.sum(volumes * diff * diff)))` is now byte-identical in FOUR
files in `tests/sn/verification/mms/`: `test_mms_curvilinear.py:49`,
`test_curvilinear_aniso_convergence.py:65`, the new file `:134`, and
`test_curvilinear_pole_cell_characterization.py:155` (as `_l2`, "production-gate
norm" — the naming drift already starting). + same body local-named in
`test_mms_2d.py:33` / `test_mms_ld_2d.py:62`. SIX copies. Rule-of-three passed long
ago. STANDING TELL: every new MMS gate in this dir re-mints `_l2_1d` rather than
importing it — flag on sight. Hoist dest = shared helper (`tests/sn/_test_helpers.py`
already imported by the dir conftest, OR a new `_helpers.py`); conftest.py itself is
hook-oriented so prefer the helper module. The convention is load-bearing (which V
power, sqrt-of-sum vs RMS-normalized, cell vs shell measure) → drift = two gates in
ONE dir disagree on "the error", real regression hides behind a norm mismatch.

**The ladder-recipe is a SECOND (acceptable-for-now) twin:** `_scalar_l2_ladder:140`
inner loop (`build_mesh→external_source→solve_sn_fixed_source(max_inner=500,
inner_tol=1e-13)→scalar_flux[0]→_l2_1d`) == mirror's `_aniso_sphere_l2_ladder:71`.
This file's version is the BETTER form (case-polymorphic `(case, n_cells)` param vs
mirror's hard-coded sphere builder). ACCEPTABLE-FOR-NOW (Pattern 6: don't hoist
before the mirror is refactored to consume it); load-bearing convention = the
`(max_inner=500, inner_tol=1e-13)` pair — pin as named constants when hoisted.
STOPS being acceptable on a 3rd ladder helper or a one-sided tol edit.

**Scalar vs per-ordinate ladder = NOT a collapsible twin (ruled OK, do NOT merge):**
3 genuine structural diffs — (1) `scalar_flux.values[0,:]` (weight-summed) vs
`angular_flux.bulk.values` (N,ng,nx) ordinate-loop; (2) the `1/W` asymmetry is
LOAD-BEARING — `phi_exact` is weight-summed (no /W) but `psi_exact(r,mu_n)` returns
`A+Bμ` WITHOUT /W by documented contract, so per-ordinate MUST `/sum_w` (`:180`);
correct to keep that convention AT the consumer (the only consumer of un-/W psi_exact,
documented as a measured trap `:64-71`+`:157-164`) — the crosswalk "load-bearing
reason required" exception. (3) scalar = one L2/rung, per-ord = max-over-ordinates.
Merging would be a strategy-flag boolean/stringly-dispatch anti-pattern hiding two
measurands. Author kept two named primitives = the elegant call.

**What the file nails (reinforce):** band constants `_CARTESIAN_SEPARABLE_MAX=0.05`/
`_CURVILINEAR_GATING_MIN=0.15` are named thresholds WITH written rationale (3-orders
headroom between measured 0.005 and 0.41 → discriminates, not pins a brittle exact
number) — Pattern-3 textbook, anti-thesis of magic threshold. `_mixed_second_difference`
returns NAMED 4-tuple `(M,dEh,dEN,rel)`, each documented. `@catches("ERR-026")`
placement RIGOROUS — author DECLINED the marker on the cross-term test (no
mutation-proven direct reddening; the proven catcher is the O(h²)-recovery assert)
and explained catches=coverage-claim-not-topic-tag. Measured-evidence-in-docstring
discipline exemplary.

**Micro-nit:** `_mixed_second_difference:191` reduces a general `(rungs,cols)` table
but hard-assumes 2 cols (`L2[i1,1]`); docstring says "2-column" but nothing enforces
→ a 3-col caller silently drops col 2. Add `assert L2.shape[1]==2` (Pattern-4 cheap
guard). Low sev (all callers `column_stack` exactly 2).
