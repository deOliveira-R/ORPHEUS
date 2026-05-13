---
name: issue-196-phase-g-step2-cylinder-fix-closeout
description: Issue #196 Phase G Step 2 cylinder fix closeout. Five `0.5` magic-number sites replaced with geometry-general `/ weights.sum()` (Σw normalisation) per the numerics-investigator's root-cause memo. 12/12 cylinder L0 streaming-equilibrium PASS at rtol=1e-9 (pre-fix 12/12 FAIL with up to 580% rel error). 12/12 sphere streaming-equilibrium PASS (no regression — sphere Σw=2 so the new form is bit-identical). 24/24 new V&V tests in `tests/sn/spatial/test_apply_matvec_cylinder_invariants.py` PASS (apply-matvec flat-flux invariant + 4-leg 3-way standoff). 11/11 regression snapshots PASS (1 cylinder snapshot `cyl_1g_homogeneous_product_dd_n20` regenerated under post-fix solver with three-pillar attestation; 10 stay bit-identical). Phase E flux-shape sentinel 2/2 PASS. Convergence-under-refinement: SI err ~ 1.5e-11, Krylov err ~ 1.5e-11, |SI−K| ~ 3e-12 at n_cells ∈ {20, 40, 80, 160} × n_mu ∈ {2, 4} × n_phi ∈ {2, 4} — fully mesh-independent at FP-noise. ERR-048 manifestation #3 entry written. The fix is a 5-line literal replacement (`0.5` → `1.0 / weights.sum()`) at 3 file sites — no new types, no new modules, no scope creep. Lesson: hardcoded numerical literals that look like `0.5` may be encoding quadrature-dependent normalisations.
metadata:
  type: project
  branch: refactor/sn-operator-algebra
  issue: 196
  phase: G Step 2 cylinder fix
  date: 2026-05-13
  commits_planned:
    - "fix(sn): replace hardcoded sphere Σw=2 normalisation with geometry-general 1/Σw in Carlson seed (Issue #196 Phase G Step 2 cylinder)"
    - "test(sn): promote L0 cylinder apply-matvec flat-flux invariant + 3-way standoff gauntlet (Issue #196 Phase G Step 2 cylinder)"
    - "chore(sn): regenerate cyl_1g_homogeneous_product_dd_n20 snapshot under post-fix solver (Issue #196 Phase G Step 2 cylinder)"
    - "docs(vv): ERR-048 manifestation #3 entry (cylinder Carlson seed Σw normalisation) (Issue #196 Phase G Step 2 cylinder)"
---

# Issue #196 Phase G Step 2 cylinder fix — closeout

## Mission

Apply the surgical 5-line literal replacement identified by the
numerics-investigator's root-cause memo
(`.claude/agent-memory/numerics-investigator/issue_196_phase_g_step2_cylinder_minimal_reproducer.md`):
replace hardcoded `0.5 = 1/Σw_sphere` with geometry-general
`1.0 / weights.sum()` at three file sites covering both the
apply-matvec path (`psi_half_angle_seed.py:569`) and the SI sweep
path (`sweep.py:543/547` sphere; `:754/756` cylinder).

Land verification: pre-fix 12/12 cylinder L0 streaming-equilibrium
FAIL → post-fix 12/12 PASS at `rtol=1e-9`; sphere streaming-
equilibrium 12/12 PASS (no regression); regression suite 11/11 PASS
(1 cylinder snapshot regenerated with three-pillar attestation;
10 stay bit-identical); Phase E flux-shape sentinel 2/2 PASS;
convergence-under-refinement at FP-noise on both inner solvers.

## Files modified

| File | Lines | Change |
|---|---|---|
| `orpheus/sn/spatial/psi_half_angle_seed.py` | 569 | `0.5 * sigma_t * phi_0` → `sigma_t * phi_0 / weights.sum()` |
| `orpheus/sn/sweep.py` | 543–547 | sphere SI: `0.5 * sigma_t_gx * phi_0_prev.T` / `0.5 * Q_1d.T` → `… / Sw` where `Sw = weights.sum()` |
| `orpheus/sn/sweep.py` | 752–756 | cylinder SI: `0.5 * sigma_t_gx * phi_0_prev.T` / `0.5 * Q_1d.T` → `… / Sw_full` where `Sw_full = weights.sum()` |

All three sites are accompanied by Hébert §3.9.4 + Pomraning
inline rationale + explicit reference to Phase G Step 2 (Issue #196)
+ explicit note that the form is bit-identical for sphere
(`Σw_sphere = 2` ⇒ `1/2 = 0.5`).

Mechanism criterion: `grep -n "0\\.5 \\* sigma_t \\* phi_0\\|0\\.5 \\*
sigma_t_gx \\* phi_0\\|0\\.5 \\* Q_1d" …` returns **empty** —
verified.

## Test promotion

New file: `tests/sn/spatial/test_apply_matvec_cylinder_invariants.py`

- **`test_cylinder_apply_matvec_preserves_flat_psi`** (12 cases):
  L0 flat-flux invariant `L·ψ_flat = Σ_t·ψ_flat` per ordinate per
  cell on the apply-matvec path.  Promoted from the numerics-
  investigator's diagnostic
  `derivations/diagnostics/diag_phase_g_step2_cyl_residual_pytest.py`.
- **`test_cylinder_three_way_standoff`** (12 cases):
  4-leg correctness standard — Krylov vs analytical, SI vs
  analytical, Krylov ≡ SI at machine precision on homogeneous
  reflective mixture B 1G cylinder.

Decorators: `@pytest.mark.l0`, `@pytest.mark.verifies("hebert-3-432")`,
`@pytest.mark.catches("ERR-048")`.

Total: 24/24 PASS in 459 s.

## Snapshot regeneration

1 snapshot regenerated: `cyl_1g_homogeneous_product_dd_n20`.

10 snapshots stay **bit-identical** under post-fix:

- slab_2g_homogeneous_dd_n20
- slab_2g_3reg_dd_n40
- sphere_2g_homogeneous_dd_n20
- sphere_2g_3reg_dd_n40
- cyl_1g_homogeneous_LS4_dd_n20  ← cylinder LS4, also stays bit-identical
- cyl_2g_3reg_LS4_dd_n40  ← cylinder MR, also stays bit-identical
- slab_2g_p1_aniso_dd_n20
- sphere_2g_p1_aniso_dd_n20
- 2d_1g_LS4_dd_15x15
- slab_fixed_source_dd_n20

The fact that only ONE cylinder snapshot needed regeneration confirms
the defect was specifically active in the ProductQuadrature path's
per-level Carlson seed; the LS4 cylinder path converged to the
correct fixed point pre-fix despite consuming the same kernel (the
LS4 per-level weight structure made the bug's residual sub-dominant
to the within-group fixed-point convergence).

### Three-pillar attestation for `cyl_1g_homogeneous_product_dd_n20`

Per `vv-principles` §"Bit-identity vs principled-equivalence":

| Pillar | Pre-fix | Post-fix | Pass criterion |
|---|---|---|---|
| **L0 streaming-equilibrium** | `keff=1.5` correct but scalar_flux = `[5.02, 2.46, 1.92, …]` (cv=0.58, max-min=3.5) — non-flat, NOT the k∞ eigenmode | `keff=1.5` exact, scalar_flux = `1.5` flat at every cell (std=1.6e-11, max-min=7.4e-11) — Pomraning isotropy + reflective symmetry | flat scalar flux on homogeneous reflective infinite-medium ✓ |
| **Pomraning pole isotropy** | cv=0.58 (DESTROYS isotropy) | cv=1.04e-11 (passes the 0.01 gate by ~9 orders of magnitude) | cv < 0.01 ✓ |
| **Variant α cross-check** | (was never tested at this snapshot) | `test_phase_e_trajectory_resolvent_flux_shape_crosscheck` 2/2 PASS — including `cyl_2g_3reg_LS4_dd_n40` heterogeneous MR cylinder | shape agreement to Phase E tol_per_cell (sphere 5e-2, cylinder 1.2e-1) ✓ |

The principled refactor passes all three criteria; the bit-identity
contract is narrowed in scope (this one snapshot regenerated) with
documented justification grounded in the structurally-independent
Variant α reference (per ``vv-principles`` §"Bit-identity vs
principled-equivalence" three-criteria gate).

## Verification paste-back

### 1. Pre-fix state — 12/12 cylinder L0 FAILED (interrupted but failures unambiguous)

The pre-fix test was run as a background process and showed 7/12
explicit FAILED markers before being killed for time. The
numerics-investigator's memo
(`.claude/agent-memory/numerics-investigator/issue_196_phase_g_step2_cylinder_minimal_reproducer.md`
§"Empirical 12/12 cylinder FAILURE — verified" lines 426-455) had
already captured the verbatim 12/12 FAIL output from a complete pre-fix
run with `max rel = 5.79` (580%) on the smallest case.

```
$ .venv/bin/python -m pytest tests/sn/spatial/test_streaming_equilibrium_curvilinear.py::test_homogeneous_streaming_equilibrium_cylinder -v
collecting ... collected 12 items

tests/sn/spatial/test_streaming_equilibrium_curvilinear.py::test_homogeneous_streaming_equilibrium_cylinder[source_iteration-4-20] FAILED [  8%]
tests/sn/spatial/test_streaming_equilibrium_curvilinear.py::test_homogeneous_streaming_equilibrium_cylinder[source_iteration-4-40] FAILED [ 16%]
tests/sn/spatial/test_streaming_equilibrium_curvilinear.py::test_homogeneous_streaming_equilibrium_cylinder[source_iteration-4-80] FAILED [ 25%]
tests/sn/spatial/test_streaming_equilibrium_curvilinear.py::test_homogeneous_streaming_equilibrium_cylinder[source_iteration-8-20] FAILED [ 33%]
tests/sn/spatial/test_streaming_equilibrium_curvilinear.py::test_homogeneous_streaming_equilibrium_cylinder[source_iteration-8-40] FAILED [ 41%]
tests/sn/spatial/test_streaming_equilibrium_curvilinear.py::test_homogeneous_streaming_equilibrium_cylinder[source_iteration-8-80] FAILED [ 50%]
tests/sn/spatial/test_streaming_equilibrium_curvilinear.py::test_homogeneous_streaming_equilibrium_cylinder[krylov-4-20] FAILED [ 58%]
tests/sn/spatial/test_streaming_equilibrium_curvilinear.py::test_homogeneous_streaming_equilibrium_cylinder[krylov-4-40] ...
```

Earlier numerics-investigator run captured the full ASSERT trace:

```
Mismatched elements: 320 / 320 (100%)
First 5 mismatches are at indices:
 [0, 0, 0]: 2.4271592681581455 (ACTUAL), 0.7957747154594768 (DESIRED)
 [0, 1, 0]: 2.0387285630992498 (ACTUAL), 0.7957747154594768 (DESIRED)
 [0, 2, 0]: 1.485717987510897  (ACTUAL), 0.7957747154594768 (DESIRED)
Max absolute difference among violations: 4.60885634
Max relative difference among violations: 5.79165969
```

### 2. Post-fix state — 12/12 cylinder L0 PASS at rtol=1e-9

```
$ .venv/bin/python -m pytest tests/sn/spatial/test_streaming_equilibrium_curvilinear.py::test_homogeneous_streaming_equilibrium_cylinder -v
collecting ... collected 12 items

tests/sn/spatial/test_streaming_equilibrium_curvilinear.py::test_homogeneous_streaming_equilibrium_cylinder[source_iteration-4-20] PASSED [  8%]
tests/sn/spatial/test_streaming_equilibrium_curvilinear.py::test_homogeneous_streaming_equilibrium_cylinder[source_iteration-4-40] PASSED [ 16%]
tests/sn/spatial/test_streaming_equilibrium_curvilinear.py::test_homogeneous_streaming_equilibrium_cylinder[source_iteration-4-80] PASSED [ 25%]
tests/sn/spatial/test_streaming_equilibrium_curvilinear.py::test_homogeneous_streaming_equilibrium_cylinder[source_iteration-8-20] PASSED [ 33%]
tests/sn/spatial/test_streaming_equilibrium_curvilinear.py::test_homogeneous_streaming_equilibrium_cylinder[source_iteration-8-40] PASSED [ 41%]
tests/sn/spatial/test_streaming_equilibrium_curvilinear.py::test_homogeneous_streaming_equilibrium_cylinder[source_iteration-8-80] PASSED [ 50%]
tests/sn/spatial/test_streaming_equilibrium_curvilinear.py::test_homogeneous_streaming_equilibrium_cylinder[krylov-4-20] PASSED [ 58%]
tests/sn/spatial/test_streaming_equilibrium_curvilinear.py::test_homogeneous_streaming_equilibrium_cylinder[krylov-4-40] PASSED [ 66%]
tests/sn/spatial/test_streaming_equilibrium_curvilinear.py::test_homogeneous_streaming_equilibrium_cylinder[krylov-4-80] PASSED [ 75%]
tests/sn/spatial/test_streaming_equilibrium_curvilinear.py::test_homogeneous_streaming_equilibrium_cylinder[krylov-8-20] PASSED [ 83%]
tests/sn/spatial/test_streaming_equilibrium_curvilinear.py::test_homogeneous_streaming_equilibrium_cylinder[krylov-8-40] PASSED [ 91%]
tests/sn/spatial/test_streaming_equilibrium_curvilinear.py::test_homogeneous_streaming_equilibrium_cylinder[krylov-8-80] PASSED [100%]

================== 12 passed, 1 warning in 1589.75s (0:26:29) ==================
```

### 3. Sphere no-regression — 12/12 sphere L0 PASS

```
$ .venv/bin/python -m pytest tests/sn/spatial/test_streaming_equilibrium_curvilinear.py::test_homogeneous_streaming_equilibrium_sphere -v
collecting ... collected 12 items

tests/sn/spatial/test_streaming_equilibrium_curvilinear.py::test_homogeneous_streaming_equilibrium_sphere[source_iteration-8-20] PASSED [  8%]
tests/sn/spatial/test_streaming_equilibrium_curvilinear.py::test_homogeneous_streaming_equilibrium_sphere[source_iteration-8-40] PASSED [ 16%]
tests/sn/spatial/test_streaming_equilibrium_curvilinear.py::test_homogeneous_streaming_equilibrium_sphere[source_iteration-8-80] PASSED [ 25%]
tests/sn/spatial/test_streaming_equilibrium_curvilinear.py::test_homogeneous_streaming_equilibrium_sphere[source_iteration-16-20] PASSED [ 33%]
tests/sn/spatial/test_streaming_equilibrium_curvilinear.py::test_homogeneous_streaming_equilibrium_sphere[source_iteration-16-40] PASSED [ 41%]
tests/sn/spatial/test_streaming_equilibrium_curvilinear.py::test_homogeneous_streaming_equilibrium_sphere[source_iteration-16-80] PASSED [ 50%]
tests/sn/spatial/test_streaming_equilibrium_curvilinear.py::test_homogeneous_streaming_equilibrium_sphere[krylov-8-20] PASSED [ 58%]
tests/sn/spatial/test_streaming_equilibrium_curvilinear.py::test_homogeneous_streaming_equilibrium_sphere[krylov-8-40] PASSED [ 66%]
tests/sn/spatial/test_streaming_equilibrium_curvilinear.py::test_homogeneous_streaming_equilibrium_sphere[krylov-8-80] PASSED [ 75%]
tests/sn/spatial/test_streaming_equilibrium_curvilinear.py::test_homogeneous_streaming_equilibrium_sphere[krylov-16-20] PASSED [ 83%]
tests/sn/spatial/test_streaming_equilibrium_curvilinear.py::test_homogeneous_streaming_equilibrium_sphere[krylov-16-40] PASSED [ 91%]
tests/sn/spatial/test_streaming_equilibrium_curvilinear.py::test_homogeneous_streaming_equilibrium_sphere[krylov-16-80] PASSED [100%]

================== 12 passed, 1 warning in 898.44s (0:14:58) ===================
```

### 4. Regression suite — 11/11 PASS (1 regen + 10 bit-identical)

```
$ .venv/bin/python -m pytest tests/sn/regression/ -v
collecting ... collected 11 items

tests/sn/regression/test_dd_regression.py::test_dd_regression[slab_2g_homogeneous_dd_n20] PASSED [  9%]
tests/sn/regression/test_dd_regression.py::test_dd_regression[slab_2g_3reg_dd_n40] PASSED [ 18%]
tests/sn/regression/test_dd_regression.py::test_dd_regression[sphere_2g_homogeneous_dd_n20] PASSED [ 27%]
tests/sn/regression/test_dd_regression.py::test_dd_regression[sphere_2g_3reg_dd_n40] PASSED [ 36%]
tests/sn/regression/test_dd_regression.py::test_dd_regression[cyl_1g_homogeneous_LS4_dd_n20] PASSED [ 45%]
tests/sn/regression/test_dd_regression.py::test_dd_regression[cyl_1g_homogeneous_product_dd_n20] PASSED [ 54%]
tests/sn/regression/test_dd_regression.py::test_dd_regression[cyl_2g_3reg_LS4_dd_n40] PASSED [ 63%]
tests/sn/regression/test_dd_regression.py::test_dd_regression[slab_2g_p1_aniso_dd_n20] PASSED [ 72%]
tests/sn/regression/test_dd_regression.py::test_dd_regression[sphere_2g_p1_aniso_dd_n20] PASSED [ 81%]
tests/sn/regression/test_dd_regression.py::test_dd_regression[2d_1g_LS4_dd_15x15] PASSED [ 90%]
tests/sn/regression/test_dd_regression.py::test_dd_regression[slab_fixed_source_dd_n20] PASSED [100%]

================== 11 passed, 3 warnings in 640.78s (0:10:40) ==================
```

### 5. Phase E flux-shape sentinel — 2/2 PASS

```
$ .venv/bin/python -m pytest tests/sn/test_phase_c_crosscheck.py::test_phase_e_trajectory_resolvent_flux_shape_crosscheck -v
collecting ... collected 2 items

tests/sn/test_phase_c_crosscheck.py::test_phase_e_trajectory_resolvent_flux_shape_crosscheck[sphere_2g_3reg_dd_n40] PASSED [ 50%]
tests/sn/test_phase_c_crosscheck.py::test_phase_e_trajectory_resolvent_flux_shape_crosscheck[cyl_2g_3reg_LS4_dd_n40] PASSED [100%]

=================== 2 passed, 1 warning in 790.62s (0:13:10) ===================
```

Note: the sentinel was already PASS-required (no xfail-strict marker
in the current source); Phase G Step 2 Path C sphere fix had already
removed the xfail-strict gate.  The cylinder case `cyl_2g_3reg_LS4`
PASSES because cylinder LS4 was already correct pre-fix at the
heterogeneous MR snapshot's discretisation level — but the new fix
makes it correct **on every quadrature** including ProductQuadrature.

### 6. Convergence-under-refinement — Krylov ≡ SI ≡ analytical at FP-noise

```
$ .venv/bin/python /tmp/cyl_refinement_postfix_production.py

=== Cylinder post-fix convergence under refinement ===
Production solve_sn_fixed_source, n_mu=2, n_phi=2
n_cells    SI err_ψ         Krylov err_ψ     SI vs K
------------------------------------------------------------
20         1.844e-11        1.464e-11        3.800e-12
40         1.867e-11        1.536e-11        3.308e-12
80         1.887e-11        1.531e-11        3.559e-12
160        1.904e-11        1.455e-11        4.485e-12

Production solve_sn_fixed_source, n_mu=4, n_phi=4
n_cells    SI err_ψ         Krylov err_ψ     SI vs K
------------------------------------------------------------
20         1.689e-11        1.473e-11        2.151e-12
40         1.694e-11        1.540e-11        1.546e-12
80         1.702e-11        1.457e-11        2.449e-12
160        1.714e-11        1.461e-11        2.537e-12
```

**Both inner solvers mesh-independent at FP-noise across the full
refinement sweep.** Compare to pre-fix: SI err_ψ DIVERGED from 7.9
(n=20) to 14.5 (n=80) and Krylov err_ψ from 2.2 to 3.5 (cited in
numerics-investigator memo §"Convergence under refinement — the smoking
gun" lines 105-119).

### 7. 3-way standoff at minimal config

```
$ .venv/bin/python -m pytest derivations/diagnostics/diag_phase_g_step2_cyl_minimal_2x4.py -v -s

=== Cylinder Phase G Step 2 minimal 2x4 (n_cells=2, N=4) ===
weights.sum() = 12.566371 (= 4π = 12.566371)
n_levels = 2
Analytical ψ_n per ordinate per cell = 0.795775
Analytical φ per cell = 10.000000

SI ψ[:, :, 0] (N x nx):
[[0.79577472 0.79577472]
 [0.79577472 0.79577472]
 [0.79577472 0.79577472]
 [0.79577472 0.79577472]]
SI φ[:, 0] (nx,): [10. 10.]

Krylov ψ[:, :, 0]:
[[0.79577472 0.79577472]
 [0.79577472 0.79577472]
 [0.79577472 0.79577472]
 [0.79577472 0.79577472]]
Krylov φ[:, 0]: [10. 10.]

max|ψ_SI - ψ_ref|     = 1.845701e-11   (rel: 0.000%)
max|ψ_K - ψ_ref|      = 2.433342e-11   (rel: 0.000%)
max|ψ_SI - ψ_K|       = 5.876410e-12

SI converges to analytical? True
Krylov converges to analytical? True
SI agrees with Krylov?         True
PASSED
```

### 8. NEW V&V test — `test_apply_matvec_cylinder_invariants.py` 24/24 PASS

```
$ .venv/bin/python -m pytest tests/sn/spatial/test_apply_matvec_cylinder_invariants.py -v
collecting ... collected 24 items

tests/sn/spatial/test_apply_matvec_cylinder_invariants.py::test_cylinder_apply_matvec_preserves_flat_psi[2-2-10] PASSED [  4%]
…(11 more flat-psi cases, all PASS)…
tests/sn/spatial/test_apply_matvec_cylinder_invariants.py::test_cylinder_three_way_standoff[2-2-10] PASSED [ 54%]
…(11 more 3-way-standoff cases, all PASS)…

================== 24 passed, 1 warning in 459.40s (0:07:39) ===================
```

## Mechanism criteria — all GREEN

- ✅ `grep -n "0\\.5 \\* sigma_t \\* phi_0\\|0\\.5 \\* sigma_t_gx \\* phi_0\\|0\\.5 \\* Q_1d" orpheus/sn/spatial/psi_half_angle_seed.py orpheus/sn/sweep.py` → empty
- ✅ `pytest tests/sn/spatial/test_streaming_equilibrium_curvilinear.py -v` → 24/24 PASS (sphere 12 + cylinder 12)
- ✅ `pytest tests/sn/regression/ -v` → 11/11 PASS
- ✅ `pytest tests/sn/test_phase_c_crosscheck.py::test_phase_e_trajectory_resolvent_flux_shape_crosscheck -v` → 2/2 PASS

## Scope discipline

**IN-scope, executed:**
- `psi_half_angle_seed.py:569` literal replacement ✓
- `sweep.py:543, 547, 754, 756` literal replacements ✓
- L0 cylinder gauntlet test promotion ✓
- 1 snapshot regenerated (ProductQuadrature cylinder homogeneous) ✓
- ERR-048 manifestation #3 entry ✓

**OUT-of-scope, deferred:**
- Anything involving L, C, S, F as `LinearOperator` instances → Phase G Step 3
- `transport_operator_matvec_*` BC trace at `operator.py:713`
  (Phase E sentinel doesn't need it — already xpasses; left for the
  Phase G Step 3 unification)
- New types, new strategy classes, new Protocols

## Architectural lessons

Per `coding-elegance` Pattern 7 (normalise at the definition site,
not at every consumer) and Pattern 3 (named intermediates with
domain semantics — physical units):

**Anti-pattern caught (the bug):** The literal `0.5` is the unnamed
arithmetic encoding of the per-level quadrature normalisation
`1/Σw_level`.  The literal works for sphere because `Σw_sphere = 2`
makes `0.5 = 1/2 = 1/Σw_level` coincide.  For 3D quadratures the
identity breaks (`Σw_level = 2π`, not 2) and the literal becomes
arithmetically wrong.

**Pattern caught (the fix):** `1.0 / weights.sum()` consumes the
typed quadrature object directly.  The intermediate IS the
per-level normalising factor — a named, units-aware quantity
(`1/Σw_level` has units of inverse-steradian, the natural
normalisation for `4π`-bounded angular integration).  The form
self-adapts to any quadrature family.

**Lesson for future Phase G Steps:** Audit every hardcoded numerical
literal in the SN operator codebase for "this looks like
`Σw_sphere/2` or `1/Σw_sphere`" patterns.  The cylinder ProductQuadrature
case will detect any remaining instances.

## Discipline lesson

The prior Phase G Step 2 Path C closeout claimed "12/12 cylinder PASS"
without running the cylinder test.  The empirical pre-fix run
(this manifestation's catching test) shows 12/12 FAIL with 580 % rel
error on the smallest case.  **Hallucinated test results are a
session failure** per `vv-principles` § "Log every caught bug" and
Cardinal Rule 1 (correctness is CRITICAL).

The mandate from this brief — VERBATIM paste-back of real test runs
in every closeout — is the structural defense: a future audit can
trace each numeric claim to its real test execution.

## Pointers

- Branch: `refactor/sn-operator-algebra`
- Pre-fix branch tip: `12c9ac9` (numerics-investigator memo committed)
- Numerics-investigator root-cause memo:
  `.claude/agent-memory/numerics-investigator/issue_196_phase_g_step2_cylinder_minimal_reproducer.md`
- Prior closeout (sphere-only fix, cylinder claim unverified):
  `.claude/agent-memory/method-implementer/issue_196_phase_g_step2_path_c_closeout.md`
- Sphere precedent (Phase G Step 2 Path C two-defect fix):
  `.claude/agent-memory/numerics-investigator/issue_196_phase_g_step2_minimal_reproducer.md`
- ERR-048 (manifestations #1+#2 sphere, #3 cylinder):
  `.claude/skills/vv-principles/error_catalog.md`
- Promoted V&V test:
  `tests/sn/spatial/test_apply_matvec_cylinder_invariants.py`
- Promoted V&V (existing, post-fix gates):
  `tests/sn/spatial/test_streaming_equilibrium_curvilinear.py`

## Linked memories

- `[[issue-196-phase-g-step2-path-c-closeout]]` — prior closeout had
  unverified cylinder claim.  Sphere portion CORRECT; cylinder portion
  was never empirically gated.  This closeout supersedes the cylinder
  side.
- `[[issue-196-phase-g-step1-closeout]]` — Phase G Step 1 wrapper-only
  `SNCellOperator` + `AngularRedistribution` promotion completed; the
  Step 2 cylinder fix lands beneath Step 1's architecture without
  touching the LinearOperator surface.
- `[[issue-196-phase-g-step2-replan-blocker]]` — the M-M Picard
  divergence investigation from earlier in this session.  The
  cylinder fix is INDEPENDENT of that line of investigation
  (Path C is the Q_bar normalisation; the Picard-divergence work
  was about angular operator inversion strategies for the broader
  Step 3 architecture).
