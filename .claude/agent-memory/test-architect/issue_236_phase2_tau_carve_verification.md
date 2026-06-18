---
name: issue-236-phase2-tau-carve-verification
description: #236 Phase 2 τ-carve crosswalk + 4-leg spec — relocate Morel-Montry τ + c_in/c_out from geometry factory to angular closure (bit-identical)
metadata:
  type: project
---

# #236 Phase 2 — τ carve verification (SPEC, branch feature/sn-spatial-angular-product)

⭐⭐ **THE KEY FINDING (live-grepped, supersedes the 2-site brief framing):** τ
(`tau_mm`, BMC-2010 Eq.43 angular weight) is consumed via a **FOUR-site
`c_in`/`c_out` duplication**, not two. `c_out=α_out/τ`, `c_in=(1−τ)/τ·α_out+α_in`
rebuilt byte-for-byte at: **P0** `pole_angular_closure.py:692-693` (matvec owner,
reads `reduced.tau_mm`), **P1** `cell_balance.py:313-314` (scalar solve, reads
`st.tau_mm`), **P2** `diamond.py:306-307` (residual @n_mask=1), **P3**
`sweep_cache.py:309-310` (CumprodScan precompute). `diamond.py:295-303` ALREADY
flags P2 as a "Phase 6 cleanup target" → the carve discharges that debt.

**Producers** (`reduced_operator.py`): SPHERE `:681-688` weight-sum edges,
UNCLAMPED (W1 declamp); CYLINDER `:798-815` per-level η-midpoint edges, CLAMP
`max(0.5,min(1,τ_raw))` RETAINED (structural τ_raw=0 ÷0 at most-inward ord, #229);
SLAB `:495` synthetic neutral `1.0`. Closure binds `quad` (`:648`) → ALREADY HAS
`(mu_x,weights,level_indices)` to produce τ; today just READS `reduced.tau_mm`
(`:660`) / `reduced.tau_mm_per_level` (`:669`).

⭐ **STRUCTURALLY-INDEPENDENT τ REFERENCE** = `contamination.py:156-201`
`morel_montry_weights(quad, geom)` (diff code path, same BMC weight, via
`_cell_edge_cosines:43-92`). **UNCLAMPED both geoms** → SPHERE bit-identical to
producer; CYLINDER equals producer's PRE-CLAMP τ_raw ONLY (`quad.eta==quad.mu_x`
col-0 same data). **No test consumes it as a τ oracle today** (only
`contamination_beta` is tested `test_quadrature.py:375-388`) — Leg-1 producer-eq
is the FIRST consumer. Cylinder gate MUST compare vs `clamp(reference)` + a
NEGATIVE control that the raw reference DISAGREES where clamp bites (anti-pat #11).

**CROSSWALK SCOPE RULING:** present NARROW (move τ-production only, leave
c_in/c_out 4-site dup as documented follow-up — cheap, low-risk, but the
SSoT-violation persists) vs WIDE (closure owns c_in/c_out, all 4 sites read it —
discharges the Phase-6 debt, but P1/P2/P3 re-association may cost bit-identity →
the curvilinear DriftWarning gate becomes the arbiter). Do NOT unilaterally pick.

⭐ **4-LEG SPEC (bit-identical carve → bit-id + producer-eq):**
- **Leg 1 producer-eq (NEW)**: closure-τ == geometry-τ 0-ULP (`np.array_equal`)
  sphere+cyl, ≥2 orders; structurally-indep leg = `morel_montry_weights`
  (sphere exact, cyl vs clamp). NEW file
  `tests/sn/sweep/curvilinear/test_tau_producer_equivalence.py` (skeleton in plan).
- **Leg 2 L14 4-leg standoff (EXISTING, keep green)**: sweep≡matvec≡struct-indep.
  `test_unified_matvec_sphere.py`+`test_unified_matvec_cylinder.py` (2p+27xf,
  route `-k "not (sphere_1g_apply_bit_identical or sphere_2g_apply_bit_identical)"`
  — stale post-ERR-058), `test_streaming_operator_decomposition.py`,
  `test_loss_action_convention.py`, `test_sweep_vs_apply_consistency.py` (THE TWIN;
  ⚠ **SLOW — full set 29p/27xf in ~11 min**), MMS `test_mms_curvilinear.py`+
  `test_curvilinear_aniso_convergence.py` (struct-indep converged value).
- **Leg 3 DriftWarning-strict (EXISTING)**: curvilinear snapshots in
  `test_dd_regression.py` = `sphere_2g_{homogeneous,3reg,p1_aniso}` +
  `cyl_{1g_homogeneous_LS4,1g_homogeneous_product,2g_3reg_LS4}`. ⚠⚠ **CRITICAL
  LIVE FINDING: `-W error::DriftWarning` does NOT escalate in the regression dir**
  — its conftest registers `always::DriftWarning` (`regression/conftest.py:30`)
  which OVERRIDES the CLI `error::`. AND `cyl_1g_homogeneous_product_dd_n20`
  ALREADY drifts 297893 ULP/3.9e-11 at HEAD (not bit-id today) → it CANNOT anchor
  a bit-id claim. The escalating strict gate lives in `tests/sn/sweep/core` +
  `tests/sn/solve` (FINDING-1 from seam note). The sha256 gate
  `test_affine_carve_bit_identity.py` covers slab `si_slab_2g_het`+2D, **NOT
  curvilinear** → no sha256 anchor for sphere/cyl exists; the curvilinear
  bit-id arbiter is Leg-2 matvec TWIN + Leg-3 snapshots run from the qualified
  path.
- **Leg 4 per-ordinate residual (L27, EXISTING)**: τ feeds the angular
  redistribution → residual MUST be per-ordinate (telescopes under Σw else).
  Reference = `test_curvilinear_operator_admits_mms.py` (vol-wt per-ordinate,
  sphere ord≈1.5/cyl≈2.0, catches ERR-058). KEEP green. Any NEW admission check
  per-ordinate.

**Mode-8 LIVE:** `-O` strips bare assert (PytestConfigWarning fires). FOUND:
`contamination_beta` tests `test_quadrature.py:380,387` are bare-assert → INERT
under -O (flag, recommend np.testing migration at touch). `admits_mms` + cyl
unified-matvec have 2 bare asserts each (run those WITHOUT -O or migrate). All
NEW gates MUST be np.testing/pytest.fail/array_equal.

**Adjoint leg:** `test_g_adjoint_reciprocity.py` (sphere/cyl/sphere_2g, -O-safe
`pytest.fail`, +L11 neg control) — c_in/c_out roles unchanged by carve (c_out in
denom diagonal, c_in in recurrence); keep green (14p w/ admits_mms).

**Live baselines @ HEAD (`.venv/bin/python -O`):** closure+quad+W1clamp 66p/5w;
curvilinear regression strict-flag 6p (DriftWarnings fire but DON'T escalate);
unified-matvec sphere+cyl 2p/27xf (routed); g_adjoint+admits_mms 14p; matvec TWIN
full 29p/27xf ~11min. Route-arounds: stale sphere apply bit-id (`-k not`), #212
keff_slab hang (NEVER all tests/sn). Plan file
`.claude/plans/issue_236_phase2_tau_carve_crosswalk.md` (skeletons).

Extends [[feedback_cross_method_protocol]] (contamination as struct-indep ref) +
[[issue_206_phase_c_verification]] (density-vs-scan nULP, DriftWarning escalation
path FINDING-1) + [[feedback_regression_tolerance_design]] (SAFETY×conv_tol).
