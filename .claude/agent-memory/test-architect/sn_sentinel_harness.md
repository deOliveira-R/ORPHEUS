---
name: sn-sentinel-harness
description: SN @pytest.mark.sentinel canary regression set — design, the -O assert-stripping hazard, cosmic-ray mutation-validation recipe, per-capability census, and the slab-DD-only-checks-cell_avg gap.
metadata:
  type: project
---

# SN sentinel (canary) regression harness

Plan: `.claude/plans/sn_sentinel_harness.md`. Marker `@pytest.mark.sentinel`
(registered pyproject). Built ON the capability taxonomy
(`sn_test_taxonomy.md`): one cheap sharp sentinel per capability NODE;
a flip localizes via the killed test's `cap()` + the DAG. Tripwire +
localizer, NOT proof. `pytest -m sentinel` = 15 node-IDs, ~4.3 s, green.

## CRITICAL hazard — `-O` strips bare `assert`
The canonical ORPHEUS invocation `python -O -m pytest`
([[feedback_default_test_mode_is_optimize]]) makes bare `assert`
statements NO-OPs. Sentinels using bare `assert` (alpha-dome,
keff homogeneous_exact) become tripwires that CANNOT trip under `-O`.
**The sentinel gate MUST run WITHOUT `-O`.** `np.testing.assert_*` is a
function call and DOES fire under `-O`; bare `assert` does NOT. The set
mixes both, so `-O` is unsafe for the sentinel gate specifically.
**Why:** a canary that can't die is worse than no canary (false green).
**How to apply:** any always-on assert-based gate → drop `-O`; OR
require `np.testing.assert_*` only. This is now **vv-principles Mode 8**
(compiled-out assertion / runtime-mode strip) — the breadcrumb landed in
the skill's failure-mode table; the skill is the canonical home.

## Mutation tool (S0 verdict)
cosmic-ray 8.4.6 over mutmut on Py3.14.3. `local` distributor (NO xdist
— deadlocks per [[sn-taxonomy-reorg-mapping]]). Scope ONE module via
`module-path`; ~1 s/mutant; `dump` carries `definition_name` →
capability node. Does NOT mutate strings. diamond.py FULL per-tier
score = 99.7 % (373/374; lone survivor `-`→`//` in residual() = minor
apply-direction gap). **GOTCHA: `local` leaves the last mutant on disk
if killed → ALWAYS `git checkout -- <module>` after a run.** NEVER run
sentinel-validation tests while a cosmic-ray exec mutates the same
module (it reverts/re-applies concurrently → false green; cost me one
contaminated run).

## Per-capability census (one sentinel/node)
primitives: per-ordinate flat-flux (ERR-006/7) + scatter-magnitude
(ERR-002) + alpha-dome (Sig5) + lebedev-4π (ERR-004). operators:
scattering add_iso_source ref. sweep_core: carlson-seed-linear (ERR-026).
sweep_slab: dd-recurrence-vs-symbolic (ERR-025). sweep_curvilinear:
NEW thin wrapper test_sentinel_sphere_streaming_equilibrium (ERR-026/48).
sweep_cartesian_2d: 2d closed_form_anchor. solve: Q/Σt cylinder (Sig4).
eigenvalue: keff_slab homog 2G (non-degenerate). verification_analytical:
NEW thin wrapper test_sentinel_kinf_slab_2g_krylov (closed-form pillar).
verification_mms: source-vanishes + p1-aniso-degrade.

## Mutation score (cosmic-ray, diamond.py) — S2 RESULT
36 sentinel node-IDs, ~4.1s gate. Sentinel-vs-diamond.py:
**96.8% (362/374; 98.1% excl. 5 equivalent ClassVar/decorator mutants)**
— update 96.1%, update_batch 100%, residual 97.9%. Journey 42.5→81.3→96.8.
GAP-CLOSER PATTERN (gave +54 pts at ~0 wall-clock): the cheap per-NODE
sentinel covered the 2-D update_batch 100% but left update()/residual()
weak — because the node sentinels (slab-DD checks only cell_avg;
curvilinear uses flat flux; NOTHING called residual()). Fix = promote
ALREADY-CHEAP DiamondDifference UNIT tests to sentinels: TestResidual
round-trip(2G/4G het)+linearity, TestBitIdenticalCurvilinear, + assert
outgoing_spatial_flux in test_dd_recurrence. LESSON: a per-capability-NODE
sentinel can leave a module's INTERIOR arithmetic uncovered; mutation-
validate per module and backfill from the module's own cheap unit tests.
12 residual survivors = equivalent (ClassVar flags, @dataclass,
Gt→GtE on structural-presence guards at geometric 0.0, constant tweaks
the round-trip tolerates).

## Known gaps (honest)
- Full heterogeneous keff (test_heterogeneous_absolute_keff, ERR-025)
  is >180 s + non-converging on field-role-typing branch → EXCLUDED;
  noted as gap. eigenvalue node covered by cheaper 2G homog + closed-form.
- S2 mutation score measured for diamond.py ONLY (the spike module).
  Other tier modules (collision/streaming/scattering/sweep_cache/pole_
  angular_closure) NOT yet mutation-validated — per-capability recipe in
  tests/_mutation/README.md; future work to score each tier's module.

## Thin-wrapper pattern for parametrized-matrix sentinels
To pin ONE cheap config without marking a whole stacked-parametrize
matrix (3 stacked @parametrize = cross-product; pytest.param marks a
SLICE not a tuple), add a 6-line `@sentinel` wrapper that delegates to
the verified function with explicit kwargs. Used for curvilinear sphere
+ kinf slab-2g-krylov. For SINGLE-axis parametrize, `pytest.param(...,
marks=pytest.mark.sentinel)` isolates one value cleanly (used for
solve-cylinder + keff-slab-2eg).
