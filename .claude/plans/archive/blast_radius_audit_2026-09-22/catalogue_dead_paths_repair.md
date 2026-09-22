# Error-catalogue dead test paths — repair record (2026-09-21)

Input: `scratch/_blast_radius/catalogue_dead_paths.md` — 14 dead paths, 28 citations.
Output: `docs/theory/verification/error_catalog.rst` only. No commit, no Sphinx build.

## Pre-flight: the list file's claims were re-measured, not trusted

- `ls` on all 14 paths: **14 of 14 MISSING** — the list file's existence claim holds.
- Independent census (regex `tests/[A-Za-z0-9_./-]+\.py` over the 8 189-line file,
  `pathlib.exists` per hit): **95 distinct paths cited, 163 citations, 14 dead / 28 dead
  citations** — reproduces the list file's rows and line numbers exactly.

## The convention this pass applies

The catalogue's `tests/…py` spellings are **navigational**: a reader pastes one into
pytest. So a full path is a promise the file is there.

- **A live gate** is spelled as a full path `tests/…py[::Sym[::sym]]`.
- **A retired gate** is named by module BASENAME plus its former directory
  (`the ``test_cartesian`` module, then directly under ``tests/sn/```), with the
  retirement date and commit hash. History is preserved and recoverable by
  `git show <hash>^:<path>`, and the census's predicate stays meaningful.

This is why the final census can be 0 without deleting any history: the 13 citations
that are genuinely historical (the HISTORY rows below) were re-spelled, not removed.

## Instruments, in the order they decided

1. **`catches` markers, by AST** (not grep). A `grep 'catches("ERR-0NN")'` MISSES the
   multi-arg form: it read **0 catchers for ERR-051**, whose catcher carries
   `@pytest.mark.catches("ERR-039", "ERR-051")`. It also misses module-level
   `pytestmark` lists (ERR-050, ERR-055). The AST census over 612 test files found all.
2. `git log --diff-filter=D` / `--follow` / `git show --stat` per dead path.
3. The cited test-function name, located by **exact** AST name match across `tests/`.
   `grep -l "class TestMultiGroupMultiRegion"` prefix-matched
   `TestMultiGroupMultiRegionSpherical` and produced a WRONG re-point, caught by the
   second census below.
4. The successor module's own docstring, where it names the files it consolidated
   (`test_curvilinear_aniso_convergence.py` names both dead aniso-MMS modules).

## Per-citation disposition

| # | entry:line | old citation | new | instrument |
|---|---|---|---|---|
| 1 | ERR-025:1128 | `tests/sn/test_cartesian.py::test_heterogeneous_absolute_keff` | **RE-POINT** → `tests/sn/eigenvalue/test_keff_slab.py::test_heterogeneous_absolute_keff` | AST `catches("ERR-025")`; `390b9e4d` deleted the host 2026-06-01, `bc8d68f9` landed the slab-keff half. Body AST-diffed: fixture/reference/5e-4 tolerance verbatim, only `result.keff` → `result.outcome.keff` |
| 2 | ERR-026:1212 | `tests/sn/test_sweep_operator_inconsistency.py` | **HISTORY** — retired 2026-05-16, `293e097a` | `--diff-filter=D` |
| 3 | ERR-026:1218 | `tests/sn/test_phase_c_mms.py` (present-tense "additional tagged tests") | **RE-POINT** → `tests/sn/verification/mms/test_curvilinear_aniso_convergence.py` | `bc8d68f9` stat: 175 lines deleted, 272 added in the successor, whose docstring names it |
| 4 | ERR-026:1305 | `:file:`…test_sweep_operator_inconsistency.py`` ("confirms this") | **HISTORY** + live successor `…/test_streaming_equilibrium_curvilinear.py` | `293e097a`; successor confirmed by AST `catches("ERR-026")` |
| 5 | ERR-026:1331 | `:file:`…test_mms_curvilinear_aniso_dd_convergence.py`` (in the input list; missed by my first edit pass, caught by the census re-run) | **HISTORY** + successor | `bc8d68f9` |
| 6 | ERR-026:1389 | `tests/sn/test_snstreamingoperator.py::test_apply_spherical_constant_flux…` (already annotated as history) | **HISTORY, de-spelled** — deleted 2026-05-29, `05864bf6` | already in the prose; hash confirmed |
| 7 | ERR-026:1441 | `tests/sn/test_sweep_operator_inconsistency.py` (Phase-B inventory) | **HISTORY** + successor | `293e097a` |
| 8 | ERR-026:1446 | `tests/sn/spatial/test_pole_angular_closure.py` | **RE-POINT** → `tests/sn/sweep/curvilinear/test_angular_closure.py` (ALIVE) | `105ce125` moved it; P4.9b dropped "pole" 2026-08-28 (the catalogue already knew the successor and left the old path spelled) |
| 9 | ERR-026:1456 | `tests/sn/l1_analytical/test_pole_closure_flat_flux_identity.py` | **RETIRED, no successor** — 2026-05-18, `2f24ae4a`, with the LegacyTau/`BoundaryFaceFlux` objects whose parity it pinned | `--diff-filter=D` + commit subject |
| 10 | ERR-026:1460 | `tests/sn/l1_analytical/test_mms_curvilinear_aniso_dd_convergence.py` | **RE-POINT** → `…/test_curvilinear_aniso_convergence.py` | `bc8d68f9` |
| 11 | ERR-026:1486 | `tests/sn/spatial/test_boundary_face_flux.py` ("deleted" — already past tense) | **HISTORY, de-spelled** — `3fd1302f`, 2026-05-12 | `--diff-filter=D` |
| 12 | ERR-026:1523-4 | two `…test_mms_curvilinear_aniso_dd_convergence.py::…` rows under "tripwires STAY xfail" | **RE-POINT** → `…/test_curvilinear_aniso_convergence.py::` (same two function names) + tense fixed to past | AST: both functions present with `catches("ERR-026")`. Added the true outcome, read from the live modules: isotropic xfails removed 2026-06-12 (#195/ERR-058), aniso 2026-06-13 (#229, `b2d8a6d`) |
| 13 | ERR-026:1629 | `tests/sn/test_snstreamingoperator.py` (Phase-D module list) | **HISTORY** — `05864bf6`, 2026-05-29 | `--diff-filter=D` |
| 14 | ERR-033:2322 | `tests/derivations/test_peierls_slab_multiregion.py` | **RE-POINT** → `tests/derivations/test_peierls_slab_legacy_aggregate.py::TestSlabSingleSurfaceVsPerFaceAggregate` | AST `catches("ERR-033")`. ⚠ `git log` knows NO such path — the old citation was **aspirational, never a record** ("extended … will surface"); said so in the entry |
| 15 | ERR-039:3204 | `tests/numerics/test_projection_operators.py::TestApplyTransposeIsWWeightedAdjoint::…` | **HISTORY, de-spelled** — retired 2026-06-23, `501d4431` (Frame campaign) | `--diff-filter=D` |
| 16 | ERR-039:3257 | same module, Round-2 sync note | **HISTORY, de-spelled** + date/hash added | `501d4431` |
| 17 | ERR-049:4310 | `tests/sn/test_invertible_operator.py` inside a run-log code block | **HISTORY** — class name alone kept in the log; the address mapping moved to the prose below | — |
| 18 | ERR-049:4320 | `…::TestInvertibleSolveBridgeRegression` | **RE-POINT** → `tests/sn/operators/test_streaming_collision_operator.py::TestStreamingCollisionSolveBridgeRegression` | `105ce125` (2026-06-01 move) then `8367346f` (2026-07-28 `InvertibleOperator → StreamingCollisionOperator`). Also re-pointed methods 1 and 2 to their `_composite` names — the unsuffixed twins exist nowhere (AST) |
| 19 | ERR-049:4340 | `…::TestSolve::test_solve_consumes_per_ordinate_rhs` | **RE-POINT**, same class/function, new module | AST: present at `TestSolve` line 484 |
| 20 | ERR-050:4523 | `…::TestSolve::test_solve_forwards_explicit_initial_guess_to_sweep` | **RETIRED, successor named** — gone 2026-07-05 `8cf52153` with the channel it pinned (#280 2.5c retired `initial_guess` threading). Surviving pin: `…::TestSolveTimedFullField::test_rhs_boundary_seeds_the_sweep_inflow` | `git log -S` on the function name; successor read in full (patches `CumprodScan.sweep`, asserts inflow from `rhs.boundary`) |
| 21 | ERR-051:4570 | `tests/numerics/test_projection_operators.py:368-381` | **HISTORY, de-spelled** + `501d4431` | `--diff-filter=D` |
| 22 | ERR-051:4592 | `…::TestGalerkinIdempotencyOnLebedev::test_pi_R_is_identity_on_band_limited` | **HISTORY, de-spelled** + date/hash | `501d4431` |
| 23 | ERR-055:4906 | `tests/sn/test_spherical.py`, `tests/sn/test_cylindrical.py` (Module line) | **HISTORY** + forward-pointer to the test reference | `390b9e4d` |
| 24 | ERR-055:4931 | six tests across the two dead modules (Test reference) | **RE-POINT, all six ALIVE** across three modules — see below | exact-AST name location; `390b9e4d` deleted, `0241df70` created the curvilinear successors |
| 25 | ERR-067:5224 | `tests/sn/operators/test_starting_direction_metric.py::test_derive_gsd_and_close_mode12` | **RE-POINT** → `tests/sn/operators/test_radial_characteristic_metric.py::test_derive_gsd_and_close_mode12` | AST `catches("ERR-067")`; `b015e362` (2026-07-07) renamed the whole `StartingDirection → RadialCharacteristic` family, function name unchanged |

ERR-055's six, at today's addresses:

- `tests/sn/sweep/curvilinear/test_sph_sweep_regression.py::TestSphericalSweepRegression::{test_uniform_source_converges_to_Q_over_sigt, test_single_sweep_all_finite}` (`l0`)
- `tests/sn/eigenvalue/test_keff_curvilinear.py::TestMultiGroupMultiRegionSpherical::test_fixed_source_flux_bounded` (`l2`)
- `tests/sn/sweep/curvilinear/test_cyl_sweep_regression.py::TestCylindricalSweepRegression::test_single_sweep_all_finite` (`l0`)
- `tests/sn/sweep/curvilinear/test_cyl_sweep_regression.py::TestAzimuthalRedistribution::{test_redistribution_telescoping_conservation, test_single_cell_uniform_source_equilibrium}` (`l2`) — class was `TestMultiGroupMultiRegion`, divided by `0241df70`
- plus the explicit catcher added later: `tests/sn/sweep/core/test_sweep_ng2_layout_guard.py`
  (`pytestmark = [pytest.mark.foundation, pytest.mark.catches("ERR-055")]`)

## Tally

**25 citation sites carrying the 28 input citations** (three sites hold two each: ERR-026:1523-4,
ERR-055:4906, ERR-055:4931). Counted from the table above:

- **10 re-pointed** to a live gate (ERR-033's counted here; it is additionally marked as an
  aspirational pointer that never existed);
- **2 retired** — ERR-050's with a named successor pin, ERR-026's flat-flux identity with
  none (its subjects retired with it);
- **13 rewritten as dated history** (5 of them also name the live successor inline).

## Adjacent corrections made in the same sites (verified against the live tree)

- **ERR-033** — heading was `**L1 test that catches it:**`; the catcher is
  `@pytest.mark.foundation`. Per `vv-principles`, `foundation` is ORTHOGONAL to the
  L0–L3 ladder, so the level claim was wrong. Now `**Test that catches it:**` with the
  marker named and the reason stated (an algebraic identity between two primitives, not
  a solver claim).
- **ERR-051** — the entry said the catcher carries `catches("ERR-039")` only, "if the
  marker is added in a follow-up". It has carried `catches("ERR-039", "ERR-051")` since;
  corrected, and noted that it is this entry's sole catcher.
- **ERR-026:1519** — "tripwires STAY xfail" was present-tense-false; all four xfails
  were removed 2026-06-12/13. Past-tensed with the two dates and causes.
- **ERR-026:1529** — dead glob `tests/sn/test_phase_c_*.py` re-spelled as a module family.
- **ERR-067 Module line** — added that `starting_direction_space.py` is today
  `orpheus/numerics/spaces/radial_characteristic_space.py`, `for_levels` on
  `_RadialCharacteristicSubSpace`.

## Verification

Two independently-vocabularied censuses, each with a positive control (X1/X2).

1. **Dead paths** (the task's predicate), controls: the regex must strip a `::selector`;
   a known-dead path must read missing AND a known-live path must read present; the hit
   count must exceed 50 or the harness is broken.
   → **88 distinct paths, 156 citations, DEAD = 0.**
   (One control failed first run — I had assumed `tests/sn/conftest.py`; it does not
   exist. Replaced with `tests/conftest.py`, verified.)
2. **Dead selectors** — every `tests/…py::Sym[::sym]`, brace-aware, resolved against the
   target module's AST. **97 selectors checked.** This census caught a real error in my
   own edit (#24: `grep -l "class TestMultiGroupMultiRegion"` prefix-matched
   `…Spherical`), which was then corrected. 9 residual unresolved: 4 are filter
   artifacts (a module-level constant `_TIGHT_KW`; truncated brace prefixes); **5 are
   genuine pre-existing stale selectors on LIVE files** — see NEEDS.
3. **Structural check on all 36 hunks / 164 changed lines**: indentation ≥ 3 inside every
   `error-entry` body; inline-literal balance (file total even). One flag, a pre-existing
   multi-line `` `` `` span preserved verbatim.
4. **Scope**: `git diff --stat -- docs/theory/` shows `error_catalog.rst` alone. The
   `.claude/` and `docs/development/skills/vv-principles.md` changes in the worktree are
   the concurrent archivist's; untouched here.

No Sphinx build was run (per the brief). Note the gate that build would NOT provide:
a code-xref is plain text at default severity, and nothing under `tests/` is rendered at
all — so this repair was correctness-driven via grep/AST, which is the only instrument
that reaches the `tests/` surface.

---

# Continuation: dead test NAMES (function-level), 2026-09-21

Input: `scratch/_blast_radius/catalogue_dead_functions.md` lists **12 of 77** `tests/…py::A[::B]`
citations whose file exists but whose node does not resolve. Same rules: only
`error_catalog.rst` edited, no commit, no Sphinx build.

## Why my first-pass selector census reported 5, not 12

My earlier census was **too lax**, and the orchestrator's stricter one is correct. It accepted a
name if (a) any def ANYWHERE in the file matched, including a method nested in a class, and (b)
any symbol merely STARTED WITH the cited token. A class-less citation such as
`test_ordinate_scan.py::test_pair_monoid_associativity` therefore "resolved" against
`TestPairMonoidTheorems::test_pair_monoid_associativity`, but pytest cannot run it as spelled.
The census below resolves each citation as a pytest node id: `A` top-level in the module, `B` a
member of `A`.

## Re-measured before building on it

All 12 rows were re-resolved by exact AST name over 612 test files: **12 of 12 fail**.
Five fail only because they omit a class segment, one names a sibling module, and the rest
name tests that were renamed, retired, or never existed. Owning entries were read off the
`.. error-entry::` above each line: ERR-020 (×3), ERR-021, ERR-054 (×2), ERR-060 (×3), ERR-065,
ERR-080 (×2).
Catchers came from an AST pass reading decorators, class decorators and module `pytestmark`
lists, multi-id aware. Positive control: ERR-087's catcher is found.

## Per-row disposition

| line | entry | old | shape | new / instrument |
|---|---|---|---|---|
| 752 | ERR-020 | `test_geometry.py::TestZoneSubdivision::test_equal_volume_single_zone[SPHERICAL]` | **HISTORY** | kept as the symptom report, de-spelled to `TestZoneSubdivision::…`; the class was deleted 2026-05-05 in `81b083be` (`git log -S "class TestZoneSubdivision"`) |
| 753 | ERR-020 | `…::test_equal_volume_multi_zone[SPHERICAL]` | **HISTORY** | same |
| 790 | ERR-020 | `…::TestZoneSubdivision::test_equal_volume_{single,multi}_zone` ("L0 test that catches it") | **RETIRED → successor** | `tests/geometry/test_structured_geometry.py::TestMesh1DFromGeometry::test_equal_volume_{cylindrical,spherical}_invariant` (AST `catches("ERR-020")`, module `foundation`; introduced `b5e85c2d`). Label corrected from L0 to foundation. **The scope is narrower:** single-region cylinder and sphere only. |
| 839 | ERR-021 | `test_ray_tracing.py::test_degenerate_corner_ray` | **RE-POINT** | two tests: `::test_degenerate_corner_ray_box_returns_none` and `::test_degenerate_corner_ray_trace_short_circuits` (AST `catches("ERR-021")`, module `l0`). The cited name **never existed**: fix commit `c0c16d32` shipped the two tests. Prose corrected to what they assert: the box test accepts `None` OR a valid pair, so its contract is "does not crash", not "returns `None`". |
| 4938 | ERR-054 | `test_ordinate_scan.py::test_ordinate_scan_small_attenuation` | **RE-POINT** | `::TestNumericalStability::test_ordinate_scan_small_attenuation` (class segment was missing); its sibling in the same sentence became `TestAffineStructure::test_ordinate_scan_zero_attenuation` |
| 4964 | ERR-054 | `test_ordinate_scan.py::test_pair_monoid_associativity` | **RE-POINT** | `::TestPairMonoidTheorems::test_pair_monoid_associativity` |
| 5182 | ERR-060 | `test_ld_ubld_symbolic.py::test_d2_exact_on_bilinear` | **RE-POINT** | `::TestOracleIIBilinearExactness::test_d2_exact_on_bilinear` (AST `catches("ERR-060")`) |
| 5182 | ERR-060 | `test_ld_ubld_primitive.py::test_d2_exact_on_bilinear` | **RE-POINT** | `::TestPrimitiveMatchesSymbolic::test_d2_exact_on_bilinear` (AST `catches("ERR-060")`) |
| 5182 | ERR-060 | `test_linear_discontinuous.py::test_d2_assembled_matrices_match_symbolic` (the "blind marker" note) | **RE-POINT + HISTORY** | the pin is `test_ld_ubld_primitive.py::TestPrimitiveMatchesSymbolic::test_d2_assembled_matrices_match_symbolic`; it was **never** in `test_linear_discontinuous.py`. The note's future-tense "a coverage-claim error to drop" is DONE: qa's review dropped the marker before `495af604` landed (2026-06-16, commit message), the pin's docstring now explains why it carries none, and the AST confirms 0 ERR-060 markers on it. |
| 5280 | ERR-065 | `tests/moc/test_verification.py::test_n2n_1g_analytical_keff` | **RE-POINT** | `::TestL0N2nReaction::test_n2n_1g_analytical_keff` (states it carries no `catches`); the bare `test_keff_estimator_gate.py::…` in the same sentence gained its directory `tests/sn/eigenvalue/` |
| 6850 | ERR-080 | `test_basis_domain.py::test_e1` | **HISTORY + live name** | the witness THEN was `test_e1_…_the_slabs_pairing_is_refusable` (`9b4a4d9c`, 2026-09-01); the repair `5436184e` (2026-09-02) inverted its slab leg and renamed it `…_the_slab_now_does_too` (`git log -S` on both names). The site sits inside the entry's "everything below is HISTORY" block. |
| 6866 | ERR-080 | `test_quadrature_directional.py::test_q8_4_the_1d_lift_is_still_a_FICTION_and_says_so` | **RETIRED → successor** | the fiction pin went RED as designed; `5436184e` removed it (`-def` in the diff) and added `::test_q8_4_the_1d_rule_ROUTES_its_own_measure_and_binds_the_legendre_basis` (AST `catches("ERR-080")`) |

**Counts per shape (12 rows):** 7 re-pointed (839, 4938, 4964, 5182 ×3, 5280) · 2 retired with a
named successor (790, 6866) · 3 rewritten as dated history (752, 753, 6850). Row 6850 also names
the live test.

## Adjacent corrections (verified, in the same entries)

- **ERR-054**: two lines said the `catches("ERR-054")` marker would land on the two
  `test_si_cyl_20cell_nan_regression` tests "when the fix lands". The fix landed on 2026-06-08
  (`5c373517`), but the marker went to the new `test_ordinate_scan_reset` suite, which was
  deduplicated the same day (`da8f8c9f`). The 20-cell module carries only `regression` and
  `foundation`. Both lines now name the sole catcher,
  `tests/sn/sweep/core/test_ordinate_scan_reset.py::TestOrdinateScanReset::test_ordinate_scan_multiple_and_consecutive_resets`.
- **ERR-020 / ERR-021**: level labels and assertion descriptions now match the tests' markers
  and bodies (see rows 790 and 839).

## Function-level census, my own, over the whole catalogue

Predicate: every `tests/…py::A[::B[::C]]` citation. The file must exist; `A` must be a top-level
def, class or module-level name; each later segment must be a member of the previous class.
`{a,b}` groups are expanded into one node id per alternative, `[…]` parametrize ids are
stripped, and a trailing `*` is read as a prefix.

Controls, 5 of 5 behaving as intended:
- 3 resolve: ERR-087's catcher, a brace pair where both alternatives exist, and the module-level
  `_TIGHT_KW`.
- 2 refuse: the class-less `test_n2n_1g_analytical_keff` form, and a brace pair with one bogus
  alternative, which yields exactly one unresolved id.

- **Post-pass: 79 selector citations, 0 missing files, UNRESOLVED 0.**
- **Reconciliation against the orchestrator:** the same census run on the pre-pass text, rebuilt
  in memory by reversing the replacements, reads **77 citations and 13 unresolved node ids**.
  The 13 ids are exactly the orchestrator's 12 citations, because line 790's
  `{single,multi}` brace is one citation expanding to two ids, where the orchestrator's census
  read it as one prefix. The denominator went 77 → 79 because of +1 (ERR-021 split into two),
  +1 (ERR-054 new marker line), +1 (ERR-054 reference line), +1 (ERR-065 gate gained its
  directory), and −2 (ERR-020 symptom rows are now bare names).
- Regression: **the path census still reads DEAD 0** (87 distinct, 158 citations). Across 48
  hunks, 0 body lines are indented less than 3 spaces, and the `` count is even.
