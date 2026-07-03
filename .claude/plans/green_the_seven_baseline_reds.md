# Green the 7-and-only-7 baseline reds (before P5)

> **✅ DONE 2026-06-26 — committed `574cff8`, pushed to origin/main, closed #250 + #232.**
> All 7 greened (Class A: 3 SPH `.npy` re-captured; Class B: 2-D→`level_symmetric` quad;
> Class C: sphere-only `.npz` re-capture, 42 keys frozen). numerics-investigator VERDICT:
> CORRECT (L0/L1 + mutation-proven probe). Root cause = STALE snapshots across `b2d8a6d`
> (spherical M-M τ-unclamp, Refs #229) — NOT ERR-058 (the investigator's first attribution;
> corrected vs git+#250, lesson L31). Promoted the non-flat per-ordinate gate
> (`test_curvilinear_operator_admits_anisotropic_mms.py`). Full tree `-m "not slow"`:
> **5475 passed / 0 failed**; Sphinx -W clean; matrix regenerated. This plan is archaeology —
> safe to archive. **NEXT = P5 (#50 energy condensation).**

**Created:** 2026-06-26 (pre-compaction checkpoint). **Branch:** `refactor/operator-inverse-algebra`
@ `932e769` (== `main` — the P4.5+#261+dyad+#265 campaign was ff-merged this session). Tree clean
(only policy-uncommitted `.claude/*` + untracked). **Goal:** retire the long-accepted "7-and-only-7"
red baseline so the full tree is GREEN before starting P5 (#50 energy condensation).

## Why now / framing (READ FIRST)

This session merged the campaign to main and, in doing so, ran the **FULL tree** `-m "not slow"`
(not the habitual `tests/sn tests/numerics` subset — see the gate-scope gap in
`[[reference-test-execution-env]]`). That surfaced **8 reds the subset masked**, and **every one was a
stale-test migration the campaign's refactors left behind** — NOT a production bug:

- 5× `test_bc_universal_invariants` stale `SNMethodSpace` import (`80489f5` follow-up) → fixed `addc59e`
- 2× MMS `qext` indexed 4-D after fields went principled 3-D `(N,ng,nx)` → fixed `5dd2d48`
- 1× `make_mixture` stale fixture (missing `nu`/`chi`) + drifted message → fixed `932e769`

**Strong prior:** the 7 baseline reds are the SAME phenomenon — bit-identity SNAPSHOTS / quadrature
CHOICES that the operator-algebra + field campaign changed structurally but never re-captured. The
test docstrings SAY SO (Class A: "SPH baselines were NOT re-captured"; Class C: "frozen pre-T.4
snapshots"). So greening them ≈ **re-capture the snapshots / fix the stale quad choices — AFTER
verifying the new production value is CORRECT against a structurally-independent reference.**

## ⚠️ THE NON-NEGOTIABLE DISCIPLINE (vv-principles L11 + anti-pattern #6/#10)

**NEVER regenerate a bit-identity snapshot just to make it green.** A snapshot re-captured from a
buggy production path masks the bug forever. For EACH bit-identity red, the gate to re-capture is:
1. **Verify the CURRENT production output is correct** against a STRUCTURALLY-INDEPENDENT reference —
   the campaign's own oracles: the **L1 trajectory-resolvent** (`orpheus/derivations/continuous/
   trajectory_resolvent`, Variant α) and the **L0 streaming-equilibrium analytical** (`φ=Q/(Σ_t(1-c))`,
   `ψ_n=φ/Σw`). These are structurally independent of the production sweep (lessons L14/L27).
2. Only if (1) passes, **re-capture the snapshot** from current production and commit it as the new
   reference, with a commit note citing the verification.
3. If (1) FAILS, it is a **real regression** — STOP, do NOT re-capture, investigate as a bug
   (dispatch `numerics-investigator`).

Per-ordinate, not weight-summed, for any angular-redistribution residual (L27). Multi-group (≥2g)
heterogeneous where possible (L2). Run gates under `python -O` SERIAL (xdist unstable).

## The 7 reds — inventory, classification, approach

### CLASS A — SPH vacuum-matvec bit-identity (3) — STALE SNAPSHOT (self-documented)
- `tests/sn/operators/test_bc_extraction_matvec.py::TestVacuumMatvecBitIdentity::test_vacuum_bulk_bit_identical_1d[0-SPH]`
- `…[1-SPH]`, `…[2-SPH]` (seeds 0/1/2; **only the `SPH` geometry param fails — `SLB`/`CYL` pass**)
- **Error:** `AssertionError: Arrays are not equal to 5 ULP (max is 8.77e+15 / 1.06e+15 / 4.45e+15)`
  (huge ULP ⇒ comparing genuinely-different arrays, i.e. a stale reference, not FP drift).
- **Self-documented:** the test docstring (`test_bc_extraction_matvec.py:273-331`) states SPH carries a
  "pre-existing STRUCTURAL" difference and "SPH baselines were NOT re-captured" (a fixed-seed random-ψ
  snapshot of the `(N,ng,nx,ny)` curvilinear M-M Carlson apply; the slab/cyl baselines WERE captured).
- **Approach:** the snapshot is the byte-exact reference for the *current* curvilinear vacuum-bulk
  matvec. Verify the current SPH bulk residual `(L+C).apply(ψ) − σ_t·ψ` is correct (Class A is a
  *bulk* residual, bit-identical UNCONDITIONALLY per the module header) via the L1/L0 references, then
  re-capture the 3 SPH baselines. Check whether the snapshot is an inline array, a committed `.npy`,
  or a fixture — see the docstring ("committed alongside this test", ~line 106).

### CLASS B — 2-D Cartesian + 1-D quadrature "genuine mu_y" guard (2) — STALE QUAD CHOICE
- `tests/sn/operators/test_boundary_conditions.py::TestSNBCResolution::test_2d_mesh_resolution`
  (fixture `quad = Quadrature.gauss_legendre(4)` at `:26` — a **1-D** quad: all ordinates have mu_y=0)
- `tests/sn/operators/test_native_matvec.py::TestTwoDCartesianRaises::test_two_d_cartesian_loss_action_returns_result`
- **Error:** `ValueError: Face 'ymin' requires genuine mu_y cosines, but every ordinate o[has mu_y=0]`
  (a production GUARD; note `test_boundary_conditions.py:139` carries `@pytest.mark.catches("ERR-052")`).
- **Hypothesis:** these build a **2-D Cartesian** mesh but pass a **1-D** `gauss_legendre` quad that has
  no genuine mu_y — the guard the campaign added correctly rejects it. The fix pattern is established
  elsewhere: `test_removal_form_matvec_sweep.py:170` and `test_mms_ld_2d.py:272` use
  `Quadrature.level_symmetric(...)` *explicitly because it has "genuine mu_y"*. So: **migrate the 2-D
  cases to a genuine-2-D quad (`level_symmetric` / `product_mu_phi` / `lebedev`)** — OR, if the test's
  intent is to assert the guard RAISES (the class is named `…Raises`), update the expected message /
  re-classify the assertion. **DECIDE intent first** by reading the two test bodies (is it a
  should-work case fed a wrong quad, or a should-raise case with a drifted message?).

### CLASS C — curvilinear sphere pre-T.4 StreamingOperator.apply snapshot (2) — STALE SNAPSHOT
- `tests/sn/operators/test_streaming_operator.py::TestT4cPreT4RegressionSnapshotCurvilinear::test_sphere_1g_apply_bit_identical`
- `…::test_sphere_2g_apply_bit_identical` (**sphere fails; slab principled-equivalent surviving gates pass**)
- **Error:** `AssertionError` (snapshot mismatch). Reference file: `PRE_T4_SNAPSHOTS_PATH`
  (`test_streaming_operator.py:565`); these freeze `StreamingOperator.apply` (`= (L+C)ψ − σ_t·ψ`)
  against **pre-T.4 snapshots** (`:547-565`).
- **Approach:** the campaign's curvilinear matvec legitimately evolved past the frozen pre-T.4 sphere
  snapshot (same flavor as Class A). Verify current sphere `StreamingOperator.apply` is correct via the
  L1 trajectory-resolvent + L0 streaming-equilibrium (the file already references
  "structurally-independent references (L1 trajectory-resolvent, analytical)" at `:244`), then
  re-capture the sphere snapshot in `PRE_T4_SNAPSHOTS_PATH`. If the verification FAILS → real
  curvilinear regression, investigate (do not re-capture).

## Suggested execution (post-compaction)

1. **Dispatch `numerics-investigator`** on the curvilinear matvec correctness (Classes A + C share the
   curvilinear sphere apply): "is the CURRENT sphere vacuum-bulk residual + StreamingOperator.apply
   correct vs the L1 trajectory-resolvent and L0 streaming-equilibrium?" — its verdict gates re-capture.
   (Also consider `test-architect` for the snapshot-regeneration discipline / re-capture mechanics.)
2. **Class B** is lighter (test-config) — read the 2 bodies, decide should-work-with-2D-quad vs
   should-raise, fix accordingly. Could be done inline by the main agent.
3. Re-capture snapshots ONLY after the correctness verdict. Commit each fix with the verification cited.
4. **Final gate = the FULL tree** (not the subset): `.venv/bin/python -O -m pytest -m "not slow"
   -p no:xdist --timeout=300 -p no:cacheprovider` → target **0 failed** (down from the 7 baseline).
   pyright `orpheus/` ≤ 412, Sphinx `-W` exit 0 if any docs touched.
5. THEN start P5 (#50 energy condensation) on a green tree.

## Pointers
- Bit-identity vs principled-equivalence + the re-capture discipline: `vv-principles` skill + lessons
  L11/L14/L27. Curvilinear references: lesson L14 (`derivations/continuous/trajectory_resolvent`).
- The 8 stale-test fixes this session are the template: verify-then-migrate, never green-by-fiat.
- ERR-052 (Class B guard origin): `error_catalog.md`.
