---
name: issue-249-l2-norm-singlesource
description: #249 mms/ L2-norm single-source — PASS-WITH-NITS; the slow-set 7th copy the non-slow gate is blind to + the migration-completes-the-followup stale docstring
metadata:
  type: project
---

`tests/sn/_test_helpers.py` + `tests/sn/verification/mms/*` — #249, the MUST-FIX I
(elegance-enforcer) raised in #236 Phase 3 ([[issue-236-phase3-st5-separability-gate]]).
Test-only, bit-identical. VERDICT: **PASS-WITH-NITS** (no blocking; 2 do-now nits).
Ran green: non-slow mms/ `25 passed in 58s`; migrated slow spot-check (delegated sphere
ladder + 2-D `volume_weighted_l2`) `8 passed in 386s`.

**What landed (correct + the home is right):** `volume_weighted_l2(values, reference,
volumes)` is now the SOLE norm; 8 named private copies retired+migrated (`_l2_1d`×2 /
`_l2`×1 / `_l2_2d`×2 / `_l2_error`×3). `scalar_flux_l2_ladder(case, n_cells)` +
`MMS_MAX_INNER=500`/`MMS_INNER_TOL=1e-13` hoisted; `_scalar_l2_ladder` DELETED,
`_aniso_sphere_l2_ladder` now a clean 2-line thin-wrapper (builds sphere-aniso case →
delegates). HOME = `_test_helpers.py` is CORRECT (the SN-wide shared module, 74
consumers, ALREADY solver-coupled — `legacy_proxy_matvec`/`_LC_matvec` build
`StreamingOperator`/`CollisionOperator` and call `.apply`; the `solve_sn_fixed_source`
top-level import does NOT muddy a previously-light module). mms/conftest.py is hook-only
(capability marker) → correctly NOT a helper home, exactly as my Phase-3 note specified.
2-D call-site restructure (`_l2_2d(A-B, V)` → `volume_weighted_l2(A, B, V)`) bit-identical
(`diff*diff == (A-B)²`, same op order).

⭐ **DO-NOW NIT-1 = the SLOW-SET 7th copy the non-slow gate is structurally blind to.**
`test_mms_heterogeneous.py:45` `_cell_l2(err_cells, widths)` = `sqrt(sum(widths*err*err))`
— byte-identical to `volume_weighted_l2`, IN THE SAME `verification/mms/` dir, NOT
migrated. It's a `@pytest.mark.slow` file → invisible to the brief's "25 non-slow" run AND
to my Phase-3 enumeration (which listed 6). This is EXACTLY the standing tell I documented
("every new MMS gate in this dir re-mints `_l2_1d`"). Migrates like the 2-D sites:
`volume_weighted_l2(phi_solver[:,0], phi_ref_g0, mesh.widths)` at :95/:96. **STANDING
LESSON: a de-duplication sweep scoped by a NON-SLOW gate has a slow-set blind spot — when
the MUST-FIX is "retire ALL copies of primitive X in dir D", grep the WHOLE dir
(`grep -rn "def _l2\|sqrt(np.sum(" D`) NOT just what the fast gate exercises. The slow files
are where the un-migrated twin hides.** (Sibling out-of-scope, leave: `_l2_error` in
`verification/analytical/test_mms_prescribed_inflow.py` — same primitive but DIFFERENT dir
+ documents the slab-widths/curv-volumes measure switch; defer to a cross-dir sweep, bound
it by pointing its docstring at `volume_weighted_l2`.)

**DO-NOW NIT-2 = the migration-completes-the-followup stale docstring (anti-pattern #11).**
`volume_weighted_l2` docstring (`_test_helpers.py:60-64`) still says present-tense "Several
mms/ modules carry byte-identical private copies ... the legacy copies migrate here (the
`module:tests` single-source follow-up)." #249 IS that follow-up — it just executed the
migration. Post-commit the sentence misstates state (copies are GONE, not "carry"; the
follow-up is DONE, not pending). Recurring shape: a commit that COMPLETES a migration leaves
the docstring that ANNOUNCED it as future-tense stale. Fix: past-tense + name #249 ("#249
retired the private copies; `_cell_l2` in the slow het file is the last straggler" once
NIT-1 lands → then drop the caveat).

**Reinforce (the call was elegant):** scalar-vs-per-ordinate ladder correctly NOT merged
(still 2 named primitives, the `1/W` asymmetry load-bearing — unchanged from Phase 3);
`MMS_*` knobs named-with-rationale (Pattern 3); the thin-wrapper delegation lost nothing.
