# Issue #340, bullet 2 — the systematic sweep for SN gates riding an UNCONVERGED solve

Measured 2026-08-09 at HEAD `4bcce0bd` (clean tree except two docstring-only
edits — see §1.4). Host `.venv` (Py 3.14), SERIAL, `python -O -m pytest`.

Method: escalate `orpheus.numerics.convergence.ConvergenceWarning` (landed
`d9b027d7`) to a hard error and diff the failure set against the known 9-red
baseline. **Measurement, not inspection** — the set of newly-red tests IS the
set of tests whose solve did not converge.

---

## 1. Instrument verification

### 1.1 ⛔ FINDING F0 — the DOCUMENTED escalation flag does not work

The brief, `orpheus/numerics/convergence.py:70`, the `ConvergenceWarning`
class docstring (`:107`), `_warn_if_unconverged`'s own emitted message
(`orpheus/sn/solver.py:454`), and `tests/sn/solve/test_convergence_contract.py:26`
all spell the CI contract as:

```
python -O -m pytest -W error::ConvergenceWarning
```

`[M]` That spelling is **rejected by Python's own warning-option parser**, and
therefore by pytest's:

```
$ .venv/bin/python -c "import warnings; warnings._setoption('error::ConvergenceWarning')"
_OptionError: unknown warning category: 'ConvergenceWarning'

$ .venv/bin/python -O -m pytest tests/sn/solve/test_convergence_contract.py \
      -W error::ConvergenceWarning -q
ERROR: while parsing the following warning configuration:
  error::ConvergenceWarning
AttributeError: module 'builtins' has no attribute 'ConvergenceWarning'
```

`warnings._getcategory` / pytest's `_resolve_warning_category` resolve an
**undotted** category name against `builtins`. `ConvergenceWarning` is not a
builtin, so the option never parses.

Severity: the failure is **loud** (hard `ERROR`, exit non-zero, zero tests
collected), so it cannot silently produce a false "0 new failures". But the
documented CI recipe is **present-tense FALSE at four sites**, one of which is
the runtime warning message the operator is meant to act on. A CI job wired
from any of them fails to start and looks like a config break rather than the
gate it was meant to be.

The working spelling — used for every measurement below:

```
-W error::orpheus.numerics.convergence.ConvergenceWarning
```

### 1.2 Positive control — the escalation actually bites

`tests/sn/solve/test_convergence_contract.py` cannot serve as the control:
`[M]` **all 9 of its tests pass under the flag**, because *every* starved call
in that file is wrapped in its own `warnings.catch_warnings()` /
`pytest.warns()` / `pytest.raises()`. A file that is green under the flag by
construction discriminates nothing.

So the control is an **unprotected** replica of its own starvation fixture —
`_fixed_source(max_inner=50)` on the all-reflective 3-D absorber box, copied
verbatim into a throwaway module outside the repo tree
(`/Users/rodrigo/.claude/jobs/c30e4f25/tmp/i340/test_probe_instrument.py`,
never added to `tests/`):

| leg | flag | result |
|---|---|---|
| starved `max_inner=50` | none | **PASS**, warning emitted (`last residual 1.952e-03`) |
| starved `max_inner=50` | escalated | **FAIL** — `orpheus.numerics.convergence.ConvergenceWarning ... hit max_inner=50 without reaching tol=1.000e-13` raised at `orpheus/sn/solver.py:448` |
| converged `max_inner=4000` | escalated | **PASS** (anti-dud negative control — the flag does not simply redden everything) |

⟹ the instrument is live, category-correct, and discriminating. Any negative
it reports below is trustworthy.

### 1.3 Scope caveat carried into the triage (vv #23 / lessons B6)

The escalation only sees warnings **emitted at a public SN entry** through
`_warn_if_unconverged` (`orpheus/sn/solver.py:412`, four call sites). It is
structurally blind to:

* a solve wrapped in the caller's own `catch_warnings`/`pytest.warns` (the
  contract file — deliberate, and correct);
* a truncated iteration that never reaches a public entry (a bare
  `SourceIteration`, a `GreenOperator`, a Krylov leaf), and
* `@pytest.mark.slow` tests, deselected by the canonical `-m "not slow"`
  (lessons B4). **A slow gate riding a truncation is invisible to this sweep.**

These are named so their absence from the delta is not miscounted as coverage.

### 1.4 A SECOND instrument, because the first has three structural holes

The `-W` sweep is the brief's instrument and it is sound, but its reachable
coverage is exactly *"a `ConvergenceWarning` that escapes the test body"*.
Three classes escape IT:

1. **A suppressed warning** — `warnings.catch_warnings()` / `pytest.warns()`
   in the test body. Correct for the contract file; indistinguishable from
   "did not truncate" in the delta.
2. **An `xfail`ed test** — an xfail absorbs the escalated error and the row
   reports `x` (green). `[M]` 21 xfails collected in `tests/sn -m "not slow"`.
3. ⭐ **INNER truncation at the two EIGENVALUE entries.** `solve_sn` /
   `solve_sn_adjoint` call `_warn_if_unconverged` with
   `budget_name="max_outer"`, and `IterationHistory.converged` there is
   `outcome.converged` from `power_iteration` — the **OUTER** fact only
   (`orpheus/sn/solver.py:2344-2360`). A within-group solve that hits
   `max_inner` inside a power iteration is invisible to the warning *by
   construction*. This is the same defect class #340 is about, one level in.

So a second, independent recorder was built — an in-process pytest plugin
(`/Users/rodrigo/.claude/jobs/c30e4f25/tmp/i340/inner_audit_plugin.py`) that
wraps two production functions and records regardless of pytest's warning
machinery or the test's outcome:

* **A. entry census** — wraps `_warn_if_unconverged`; sees suppressed and
  xfail-absorbed truncations.
* **B. inner census** — wraps `_certify_within_group_exit`, which every
  full-angular inner arm calls with its own `residual_history` and `tol`;
  sees inner truncation under an eigenvalue outer.

`[M]` Plugin verified on the same control pair before use:

```
ENTRY | ...::test_unprotected_starved_solve | solve_sn_fixed_source |
        max_inner=50 | tol=1.000e-13 | last residual 1.952e-03
INNER | ...::test_unprotected_starved_solve | n=1 | worst 1.952e-03 |
        tol=1.000e-13 | [solve_sn_fixed_source[source_iteration]]
```
and **silent** on `test_unprotected_converged_solve` (the anti-dud leg).

Remaining blind spot of BOTH instruments, stated so it is not miscounted:
the **windowed SI arm** skips the certificate (`solver.py:1741 —
`if not windowed:``), and `@pytest.mark.slow` rows are deselected by the
canonical `-m "not slow"` (lessons B4).

### 1.5 Tree state — checked at the END too, because it MOVED

⚠ This is a **shared working tree** and another agent edited it mid-audit
(lessons H7). Reconciled explicitly rather than assumed:

* At audit start, `HEAD = 4bcce0bd` and the only dirty tracked file under
  `orpheus/` or `tests/` was `tests/sn/solve/test_coupled_solve_certificate.py`
  — **docstring-only** (a ⚠ note pointing at the #340 contract file).
* At audit end, `HEAD = 49442156` (*"fix(derivations): the reference solvers
  that KNOW they stopped short now say so"* — the `derivations/` half of #340,
  comment 1). `[M]` `git diff --stat 4bcce0bd..HEAD` touches **nothing** under
  `orpheus/` and only that same docstring under `tests/sn/`.
* `[M]` `orpheus/sn/solver.py` acquired an **uncommitted** `+47/−27` at
  `02:51:32`, between the escalated run (finished `02:18`) and the consequence
  probes (ran `03:05–03:12`). Programmatically classified: all 74 changed lines
  are inside docstrings/comments (the 16 that match a code regex are RST table
  rows and prose about the #341 G-S rate story). **Behaviour-neutral.**

Timing per measurement, so each can be re-derived:

| measurement | ran | production code as of |
|---|---|---|
| `-W` escalated sweep | 01:47–02:18 | `4bcce0bd`, clean |
| D6 budget probe | ~02:35–02:50 | `4bcce0bd`, clean |
| census sweep | 02:26–03:03 | `4bcce0bd` (module imported at start, before the 02:51 docstring edit) |
| P1/P2 consequence probes | 03:05–03:12 | `4bcce0bd` + a docstring-only working-tree edit |

Corroboration that nothing moved underneath: the census run independently
reproduced `9 failed, 2892 passed` — the brief's baseline, exactly.

---

## 2. Raw delta vs the 9-red baseline

```
.venv/bin/python -O -m pytest tests/sn -m "not slow" \
    -W error::orpheus.numerics.convergence.ConvergenceWarning \
    -q --tb=line -rf -p no:cacheprovider

16 failed, 2885 passed, 1 skipped, 114 deselected, 61 xfailed,
6 warnings in 1838.26s (0:30:38)
```

Arithmetic checks against the brief's baseline (2892 passed / 9 failed at the
same HEAD): `2885 + 16 = 2901 = 2892 + 9`. Same collection, no drift.

**The 9 baseline reds are all present and all non-`ConvergenceWarning`** —
so the split is unambiguous, no test had to be re-run to classify it:

| baseline red | exception |
|---|---|
| `test_affine_carve_bit_identity.py::test_converged_flux_bit_identical_after_affine_carve[si_2d_p1_aniso_het]` | `AssertionError` sha mismatch (#333) |
| … `[krylov_2d_p1_aniso_het]` | same |
| … `[si_slab_2g_het]` | same |
| `test_2d_octant_sweep_equivalence.py::…[01_smoke_vacuum_1g_homog_uniformQ_LS4]` | `not equal to 64 ULP (max 1.04e+15)` |
| … `[02_reflective_1g_homog_uniformQ_LS4]` | `5.97e+14` |
| … `[03_l7_trap_mixedBC_2g_het_LS4]` | `1.07e+15` |
| … `[04_vacuum_2g_het_gradientQ_LS6]` | `5.28e+14` |
| … `[05_qaniso_mixedBC_2g_het_LS4]` | `1.13e+15` |
| `test_diamond.py::TestBitIdenticalCurvilinear::test_spherical_inward_bit_identical` | `assert False` (`np.array_equal`) |

### ⟹ THE DELTA — 7 tests, every one a `ConvergenceWarning`

Sizing data for the "derive `max_inner`" work is in the last three columns.

| # | test | entry | budget hit | tol | distance at exit |
|---|---|---|---|---|---|
| D1 | `tests/sn/acceleration/test_dsa_acceleration.py::TestTeeth::test_sign_flipped_correction_breaks_convergence` | `solve_sn_fixed_source` | `max_inner=200` | `1.000e-11` | last residual **8.378e+56** |
| D2 | `tests/sn/acceleration/test_dsa_acceleration.py::TestTeeth::test_zeroed_trace_arm_breaks_the_reflective_case` | `solve_sn_fixed_source` | `max_inner=120` | `1.000e-11` | last residual **3.614e+20** |
| D3 | `tests/sn/solve/test_krylov_curvilinear_precond_safety.py::test_g_d3_3_production_sites_size_restart_from_the_coupled_ravel[eigenvalue]` | `solve_sn` | `max_outer=2` | `1.000e-03` | last \|dk\| **4.663e-15** |
| D4 | `tests/sn/solve/test_si_single_primitive_contract.py::test_fixed_source_si_and_eigenvalue_inner_share_one_primitive[slab]` | `solve_sn_fixed_source` | `max_inner=4` | `1.000e-12` | last residual **5.275e-01** |
| D5 | `…::test_fixed_source_si_and_eigenvalue_inner_share_one_primitive[sphere]` | `solve_sn_fixed_source` | `max_inner=4` | `1.000e-12` | last residual **4.847e-01** |
| D6 | `tests/sn/solve/test_sn_adjoint_certification.py::TestP13Mutations::test_streaming_no_reversal_shifts_k_heterogeneous` | `solve_sn_adjoint` | `max_outer=500` | `1.000e-09` | last \|dk\| **3.331e-16** |
| D7 | `tests/sn/sweep/curvilinear/test_si_cyl_20cell_nan_regression.py::test_si_returns_finite_keff` | `solve_sn` | `max_outer=3` | `1.000e-10` | last \|dk\| **8.628e-02** |

### The two required confirmations

* ✅ **`tests/sn/solve/test_d3_admission.py` does NOT appear** — neither
  `test_d3_pure_absorber_per_ordinate_psi_exact` nor
  `test_d3_scattering_infinite_medium_matches_multigroup_balance`. The
  `11e78430` fix is complete for both.
* ✅ **The `test_sn_adjoint_certification` sphere-500-outer case appears**, at
  exactly the issue's number (`max_outer=500`, `|dk| = 3.331e-16`) — but the
  node that owns it is **D6, a mutation test**, not the certification gate
  itself. See §3 D6: that attribution is an artefact of `functools.lru_cache`,
  and it is itself a finding.

---

## 3. Per-test triage of the delta

**Headline: all 7 are class (a) DELIBERATE. The `-W` delta contains ZERO
class-(b) latent false greens.** Not one of the seven asserts a physical
VALUE against the truncated answer — every one asserts either a structural
fact (types, captured operands), a divergence witness (the residual must be
LARGE), or finiteness. That is a real and reassuring result for the gates the
instrument can see; it is *not* a result about the gates it cannot (§1.4 —
which is why §3.8 exists).

| id | class | why | consequence measurement |
|---|---|---|---|
| D1 | (a) | divergence witness | truncation IS the assertion |
| D2 | (a) | divergence witness | truncation IS the assertion |
| D3 | (a) | `max_outer=2` is structurally unconvergeable | flag can never be True |
| D4 | (a) | constructor spy, "mechanism, not numerics" | no value asserted |
| D5 | (a) | idem | no value asserted |
| D6 | (a)′ | mutation-induced, un-fixable by budget | `[M]` bit-identical at 4× budget |
| D7 | (a) | finiteness-only assertion | truncation cannot make finite → infinite |

### D1 / D2 — `test_dsa_acceleration.py::TestTeeth` (DELIBERATE, documented)

`test_sign_flipped_correction_breaks_convergence` (`max_inner=200`,
last residual **8.378e+56**) and `test_zeroed_trace_arm_breaks_the_reflective_case`
(`max_inner=120`, last residual **3.614e+20**) monkeypatch `DSACorrection.apply`
to break the accelerator and then **require the residual to be LARGE**
(`if last_res < 1e-6: fail` / `if last_res < 1.0: fail`). A converged exit
would FAIL them. Intent is stated in the test, twice each — docstring
("*must break*", "*diverges (measured ρ > 1)*") and an inline comment that
names the mechanism explicitly:

> `# The residual VALUE is the witness (an iteration-count bar is`
> `# off-by-one-prone: a diverging run reports max_inner−1).`

✅ documented. ⚠ but **not suppressed** — see §4/R1.

### D3 — `test_krylov_curvilinear_precond_safety.py::test_g_d3_3_…[eigenvalue]` (DELIBERATE, documented)

`solve_sn(max_outer=2, keff_tol=1e-3, …)`, exit at `|dk| = 4.663e-15`.
Docstring: *"Loose outer tolerances — the claim is the CONSTRUCTION plumbing,
not convergence"*. The test asserts a captured GMRES `restart` integer.

⭐ Sharper than "loose": **`max_outer=2` can NEVER report converged.**
`SNSolver.converged` (`orpheus/sn/solver.py:1559`) opens with
`if iteration <= 2: return False`, so the flag is pinned False by
construction regardless of the physics — `|dk| = 4.66e-15` is 12 orders
*inside* `keff_tol=1e-3`. The warning here is reporting the guard, not a
numerical shortfall. (Same mechanism the contract file uses deliberately for
its starved leg, `_eigenvalue(max_outer=1)`.) The sibling
`[fixed_source]` row does not warn: it takes `max_inner=60` and converges.

### D4 / D5 — `test_si_single_primitive_contract.py` (DELIBERATE, under-documented)

`max_inner=4` at both the eigenvalue inner (`SNSolver(sn_mesh, max_inner=4)`)
and the fixed-source entry; exits at residual **5.275e-01** (slab) /
**4.847e-01** (sphere) against `inner_tol=1e-12`. The whole body asserts
**types and captured operands** — `isinstance(L_eig, SweepOperator)`,
`type(L_eig.inner) is type(L_fs.inner)`, `len(gains) == expected_n`,
`B_eig.block_role is BlockRole.BOUNDARY`. No number from the solve is ever
read. A 4-sweep budget is exactly right for a constructor spy.

⚠ **Finding (documentation).** The docstring says *"the structural identity
of the decomposition (mechanism, not numerics)"*, which makes the intent
*inferable* — but neither `max_inner=4` carries an annotation saying the
truncation is deliberate. Per the brief's rule, an inferable-but-unstated
deliberate truncation is a finding. One comment each fixes it.

### D6 — `test_sn_adjoint_certification.py::TestP13Mutations::test_streaming_no_reversal_shifts_k_heterogeneous`

Class **(a)′** — the truncation is a *consequence of the deliberate mutation*,
undocumented, and **cannot be removed by any budget**.

`[M]` direct measurement (`probe_d6.py`, het slab, GL8, `keff_tol=1e-9`,
`flux_tol=1e-8`, `inner_tol=1e-10`, `max_inner=500`):

| run | k | converged | n_outer | last \|dk\| | \|k − k_fwd\| |
|---|---|---|---|---|---|
| forward `solve_sn` | `1.2161741441974736` | **True** | 5 | — | — |
| adjoint, UNMUTATED | `1.2161741442084186` | **True** | 5 | — | 1.1e-11 |
| adjoint, MUTATED, `max_outer=500` | `-0.6519302852190432` | False | 500 | 3.331e-16 | **1.868104e+00** |
| adjoint, MUTATED, `max_outer=2000` (4×) | `-0.6519302852190432` | False | 2000 | 3.331e-16 | **1.868104e+00** |

**The asserted quantity does not move — bit-identical at 4× budget.** The
gate needs `> 1e-6` and gets `1.87`; there is ~6 orders of margin and no
budget sensitivity whatsoever. Benign.

⭐ *Why* it cannot converge, which is the useful part: the mutated
"adjoint" has a **negative dominant eigenvalue** (`k_mut = −0.652`). A power
iteration on a sign-alternating dominant mode has a flux increment that
oscillates and never decays, so `dphi < flux_tol` is **un-satisfiable** while
`dk` sits at the FP floor (`3.331e-16 < keff_tol = 1e-9`). `SNSolver.converged`
requires BOTH (`solver.py:1564`), so `converged=False` forever. This is not a
budget shortfall at all — raising `max_outer` is futile, and the honest
remedy is `pytest.warns` / `catch_warnings`, never a bigger number.

⛔ **REFUTES the issue-comment attribution.** The #340 comment records
*"`test_sn_adjoint_certification`'s **sphere** run exhausts all 500 outers at
`|dk| = 3.3e-16`"*. At HEAD `4bcce0bd` that is **false**:
`TestP13KEquality::test_sphere_k_equality` is collected under
`-m "not slow"` (verified with `--collect-only`) and **passes** under
escalation, and `[M]` the unmutated adjoint on this module's fixtures
converges in **5** outers. The `|dk| = 3.331e-16` signature belongs to the
MUTATED heterogeneous-slab adjoint, not to the sphere. Same number, different
owner — worth correcting on the issue, because "the sphere solve is
under-budgeted" would send the `max_inner`-derivation work after a
non-existent problem.

### D7 — `test_si_cyl_20cell_nan_regression.py::test_si_returns_finite_keff` (DELIBERATE, documented)

`solve_sn(max_outer=3, max_inner=3)`, exit at `|dk| = 8.628e-02` against
`keff_tol=1e-10` — a genuine, large shortfall, and deliberately so:

> `# NaN appears in the FIRST inner iteration — small caps suffice.`

The only assertion is `np.isfinite(res.keff)`. Truncation cannot turn a
finite value non-finite, so the claim is structurally immune to the budget.
✅ documented, benign.

⚠ **Separate finding in the same file, and it is the worse one.**
`tests/sn/sweep/curvilinear/test_si_cyl_20cell_nan_regression.py:66-67`
executes, at **module import time**:

```python
import warnings
warnings.filterwarnings('ignore')
```

An unqualified, category-less, module-scope mutation of the process-global
warning filters. It happens not to defeat this audit — pytest re-establishes
filters per item inside `catch_warnings` and re-applies `-W` — but that is
pytest saving it, not the code being safe. Any non-pytest consumer that
imports this module inherits a globally silenced interpreter, and this is
precisely the failure the `ConvergenceWarning` design argues against in its
own docstring (*"a warning that always fires is noise, and noise gets
filtered, which is how the next truncation goes unnoticed"* —
`test_convergence_contract.py:196`). It should be a scoped
`pytest.ini`-style filter or a `catch_warnings` block, never a global.

---

### 3.8 What the SECOND instrument found — the census (§1.4)

Same slice, same HEAD, plugin instead of escalation:
`9 failed, 2892 passed, 1 skipped, 114 deselected, 61 xfailed in 2280 s`.
The 9 reds are exactly the baseline ledger ⟹ **the plugin is
behaviour-neutral** (it reproduces the un-escalated baseline the brief
quoted, which is itself an independent confirmation of that baseline).

### A. Public-entry census — 12 rows, and the delta is complete

The 7 delta rows **plus 5 from `test_convergence_contract.py`** — its four
starved legs (`max_inner=50`, residual `1.952e-03`) and its starved
eigenvalue leg (`solve_sn`, `max_outer=1`, `no residual recorded`). All five
are deliberate, documented, and correctly suppressed.

⟹ **Nothing is hidden behind warning suppression or an `xfail` anywhere in
`tests/sn -m "not slow"`.** The `-W` delta is the complete entry-level
picture. (This is the negative that most needed an independent instrument;
now it has one.)

### B. Inner-solve census — 90 raw rows, and ⚠ 44 of them are MY OWN false positives

⚠ **Instrument correction, made before reading anything into the list.** My
wrapper recorded "truncated" whenever `_claims_convergence(history, tol)` was
false, and that predicate is ALSO false for an **empty** `residual_history`.
`KrylovAcceleration.solve`'s docstring (`orpheus/numerics/iteration.py:937`)
states the history is *"Empty if GMRES returned in zero iterations"* — a
warm-started inner whose initial guess already satisfies `rtol`, i.e. the
**opposite** of a truncation. `[M]` all 44 such rows print
`worst last-residual inf` (my `float("inf")` sentinel for an empty list),
and no `KrylovAcceleration ... info=` `RuntimeWarning` appears anywhere in
the run, so scipy reported `info == 0` on every one. **44 rows dropped; 46
genuine.** (vv #17, applied to my own harness: the instrument's NEGATIVES
need verifying too, and this one failed in the flattering direction —
"look how much I found".)

⚠ **Second calibration, before ranking**: an inner solve that truncates on
an EARLY outer under a power iteration that later converges is an
**inexact-Newton posture**, not a defect — the outer re-solves. Census B is a
SCREENING list, not a finding list. What would be a finding is a truncated
inner on the FINAL outer, which biases the converged `k` (the documented
ERR-052 / #200 / `test_kinf_homogeneous_tolerance.py` phenomenon). The
consequence measurements below decide it per case.

The 46 genuine rows, by shortfall factor `last_residual / inner_tol`
(excluding the 12 already triaged in §3 and the direct-call `[noop-test]`
row):

| shortfall | n | worst | inner_tol | test |
|---|---|---|---|---|
| 4.96e+08× | 11 | 4.962e-03 | 1e-11 | `test_d3_admission.py::test_kinf_3d_equals_2d_equals_1d_homogeneous_reflective[2g]` |
| 2.01e+08× | 9 | 2.012e-03 | 1e-11 | … `[4g]` |
| 1.05e+07× | 6 | 1.053e-02 | 1e-09 | `test_kinf_homogeneous_tolerance.py::test_inner_tol_bias_collapses_at_1e_12[2eg-sphere]` |
| 8.92e+06× | 7 | 8.917e-02 | 1e-08 | `test_cache.py::test_collision_cache_invariance_under_source_iteration` |
| 2.06e+06× | 6 | 2.058e-03 | 1e-09 | `…test_inner_tol_bias_collapses_at_1e_12[4eg-sphere]` |
| 1.94e+06× | 6 | 1.941e-03 | 1e-09 | `…[2eg-cylinder]` |
| 1.86e+06× | 4 | 1.860e-05 | 1e-11 | `test_w1_clamp_silent_on_flat.py::test_unclamped_sphere_flux_strictly_positive[2g_R2_S64]` |
| 2.00e+05× | 3 | 2.002e-04 | 1e-09 | `…[4eg-cylinder]` |
| 7.04e+04× | 2 | 7.037e-06 | 1e-10 | `test_dd_regression.py::test_dd_regression[sphere_2g_3reg_dd_n40]` |
| 2.47e+04× | 2–3 | 2.472e-06 | 1e-10 | `test_keff_curvilinear.py::TestMultiGroupMultiRegionSpherical::{test_2g_heterogeneous_converges, test_multigroup_eigenvector_not_flat, test_particle_balance_heterogeneous}` |
| 1.89e+04× | 3 | 1.886e-06 | 1e-10 | `test_dd_regression[cyl_2g_3reg_folded_4x8_dd_n40]` |
| 1.77e+04× | 2 | 1.770e-06 | 1e-10 | `test_dd_regression[slab_2g_3reg_dd_n40]` |
| 2.79e+03× | 1–4 | 2.785e-07 | 1e-10 | `test_keff_curvilinear.py::TestCylinderMultiGroupMultiRegion::{…4 rows}` |
| ~6.1e+02× | 2 | ~6.1e-06 | 1e-08 | `test_krylov_restart_signature.py::test_si_kinf_independent_of_mesh_refinement[5,8,10,16,20,30]` |
| 2.73e+02× | 1 | 2.732e-10 | 1e-12 | `test_streaming_collision_operator.py::…test_si_carve_recovers_analytical_kinf[2eg-sphere]` |
| 45× / 21× / 16× | 1 | 4.5e-11 / 2.1e-11 / 1.6e-11 | 1e-12 | `…[2eg-cylinder]`, `…[2eg-slab]`, `…[4eg-sphere]` |
| 35× | 2 | 3.516e-09 | 1e-10 | `test_dd_regression[sphere_2g_homogeneous_dd_n20]`, `test_w1_clamp_silent_on_flat.py::test_homogeneous_reflective_sphere_iso_unchanged` |
| 9.3× | 2 | 9.254e-08 | 1e-08 | `test_boundary_conditions.py::TestSNBCSweepBehavior::{3 rows}` |
| 8.4× | 2 | 8.375e-08 | 1e-08 | `test_solver_components.py::TestComputeGroupRates::test_homogeneous_keff_matches_analytical_kinf` |
| 3.6× | 1 | 3.573e-08 | 1e-08 | `test_wavefront_cumprod_equivalence.py::test_cumprod_path_hits_analytical_kinf` |
| 1.2× | 1 | 1.247e-10 | 1e-10 | `test_sn_adjoint_certification.py::TestP13KEquality::test_heterogeneous_slab_k_equality` (and the D6 sibling) |

### 3.9 Consequence measurements on the two sharpest screening candidates

`[M]` `probe_consequence.py`, HEAD `4bcce0bd`.

**P1 — `test_d3_admission::test_kinf_3d_equals_2d_equals_1d_homogeneous_reflective[2g]`**
(worst shortfall in the whole census; asserts `|keff − case.k_inf| ≤ atol=1e-8`
against a closed-form matrix eigenvalue, `keff_tol=1e-10`, `inner_tol=1e-11`,
`max_inner` **defaulted to 200**):

| max_inner | keff | \|k − k_inf\| | \|k − k(200)\| | n_outer |
|---|---|---|---|---|
| 200 (shipped) | `1.8750000000029183` | 2.917e-12 | — | 8 |
| 800 (4×) | `1.8750000000026916` | 2.691e-12 | 2.267e-13 | 4 |
| 3000 (15×) | `1.874999999999999` | **1.998e-15** | 2.919e-12 | 3 |

⟹ **the value MOVES but stays 3400× inside the gate's own tolerance**
(2.9e-12 vs `atol=1e-8`). Verdict: benign, and *structurally* so — the
fixture is HOMOGENEOUS all-reflective, where `k_inf` is a material-property
ratio independent of the flux shape, so an inexact flux cannot bias it. The
inner truncation shows up as *outer work* instead (8 → 3 outers), which is
the honest inexact-Newton trade. Not a false green.

**P2 — `test_dd_regression[sphere_2g_3reg_dd_n40]`** (a FROZEN snapshot pin;
`max_inner=300`, `inner_tol=1e-10`, `keff_tol=1e-12`, `flux_tol=1e-10`, pinned
at `SAFETY=10 × conv_tol`):

| max_inner | keff | \|dk\| vs snapshot | max rel \|dφ\| vs snapshot |
|---|---|---|---|
| 300 (shipped) | `1.3816447394326927` | **0.000e+00** | **0.000e+00** |
| 1200 (4×) | `1.3816447394325482` | 1.446e-13 | 1.640e-12 |

⟹ the shipped budget reproduces the snapshot **bit-identically** (it is a
determinism pin and it is doing its job), and the 4×-budget answer stays
inside the pin — `1.4e-13` against `10 × 1e-12 = 1e-11` (69× margin) on
`k_eff`, `1.6e-12` against `10 × 1e-10 = 1e-9` (610× margin) on the flux.
Verdict: benign; the baselines are **not** brittle to the inner truncation,
so a future `max_inner` increase will not force a re-baseline of these cases.

---

## 4. Ranked remediation list

**There are no class-(b) active false greens to report.** Every truncation
reachable by either instrument is deliberate, and the two riskiest
inner-truncation candidates measure benign with 3+ orders of margin. The
list below is therefore ordered by *contract integrity*, not by numerical
danger.

### R1 — ⛔ Fix the escalation recipe: `-W error::ConvergenceWarning` does not parse

Four sites state a CI contract that cannot run (§1.1):
`orpheus/numerics/convergence.py:70`, `:107`; `orpheus/sn/solver.py:454`
(the *emitted message*, i.e. the one a user acts on); and
`tests/sn/solve/test_convergence_contract.py:26`. All must become
`-W error::orpheus.numerics.convergence.ConvergenceWarning`.

⭐ And the gate for it belongs in the contract file, which currently
*asserts escalatability in-process* (`simplefilter("error", …)` +
`pytest.raises`) — a leg that passes with the bare, unusable spelling.
The honest gate is on the **spelling**. `[M]` both APIs discriminate
correctly, so either serves:

```python
from _pytest.config import parse_warning_filter
parse_warning_filter("error::ConvergenceWarning", escape=False)
#   -> UsageError
parse_warning_filter("error::orpheus.numerics.convergence.ConvergenceWarning",
                     escape=False)
#   -> ('error', '', <class '...ConvergenceWarning'>, '', 0)
```

Highest priority because it is the difference between "we have a CI gate"
and "we have a CI job that errors at startup".

### R2 — Make the 7 deliberate rows survive the (now-usable) contract

Once R1 lands and a CI job runs the flag, D1–D7 all go red for entirely
legitimate reasons. Each needs its truncation **declared in the test**, the
way the contract file already does it. Preferred spelling, because it
strengthens the gate instead of merely silencing it:

* **D1, D2** (`test_dsa_acceleration.py::TestTeeth`) — wrap in
  `with pytest.warns(ConvergenceWarning):`. The mutation is *supposed* to
  produce a truncated exit; asserting that it does is a free extra tooth.
* **D6** — same, plus a comment recording the measured mechanism (negative
  dominant eigenvalue ⟹ the flux criterion is un-satisfiable ⟹ do NOT
  "fix" this by raising `max_outer`; `[M]` bit-identical at 4×).
* **D3, D4, D5, D7** — `warnings.catch_warnings()` +
  `simplefilter("ignore", ConvergenceWarning)`, matching
  `test_convergence_contract.py:147`.

### R3 — Unify the warning category: GMRES truncation is NOT escalatable

`KrylovAcceleration.solve` (`orpheus/numerics/iteration.py:1042`) surfaces
scipy's `info != 0` as a **bare `RuntimeWarning`**, so
`-W error::…ConvergenceWarning` does not escalate the Krylov half of exactly
the same defect class — and ERR-053 (the ERR the `convergence.py` module
docstring cites as its own design precedent) lives on that path. Re-point it
at `ConvergenceWarning` (a `RuntimeWarning` subclass, so no caller breaks)
and one flag then covers both primitives. `[M]` no such warning fired in this
run, so the change is inert today and purely a contract repair.

### R4 — `test_si_cyl_20cell_nan_regression.py` sets a GLOBAL warning filter

`tests/sn/sweep/curvilinear/test_si_cyl_20cell_nan_regression.py:66-67`
runs `warnings.filterwarnings('ignore')` — unqualified, category-less, at
module import. pytest's per-item filter reset is what keeps it from
poisoning the run; nothing in the file does. Replace with a scoped filter.

### R5 — Documentation nits found on the way (each one line)

* `test_si_single_primitive_contract.py:130,149` — annotate `max_inner=4`
  as a deliberate spy budget (§3 D4/D5).
* `test_kinf_homogeneous.py:163-167` and `:225-229` — the sphere-4g-krylov
  `pytest.xfail` reason and its twin comment both cite *"exceeds the
  `max_inner=300` budget"*, while `_TIGHT_KW` has said `max_inner=1000`
  since 2026-05-26 (`:129`). Present-tense-false prose in the exact string a
  future reader uses to decide whether #200 is still open.
  ⭐ And it is worse than a nit: both are **imperative** `pytest.xfail(...)`
  calls, which raise before the body runs, so the row can NEVER XPASS. The
  usual "the strict-xfail set is the todo list" self-retirement does not
  apply — if the 1000-inner budget already fixed sphere-4g-krylov, nothing
  in the suite would ever say so. Converting to
  `@pytest.mark.xfail(strict=True, reason=…)` would let the fix announce
  itself, and would immediately answer the stale-budget question.
* Issue #340's own comment attributes the 500-outer truncation to the
  adjoint **sphere** run; `[M]` it is the mutated **heterogeneous-slab**
  adjoint (§3 D6). Worth correcting in the issue so the `max_inner`
  derivation work is not aimed at a non-problem.

### R6 — Sizing data for "derive `max_inner`" (issue #340 open item 1)

From the census, the *unforced* (non-deliberate) inner-budget pressure in
`tests/sn` is concentrated and modest — the worst genuine shortfall under a
converged outer is `4.96e-03` vs `1e-11` on a homogeneous d=3 box where the
truncation provably cannot bias `k`. The load-bearing shape remains the one
already in the issue (`Σ_t·n_inner` invariant; ~an order per dimension;
all-reflective + weak absorption + d=3 is the trigger). This audit adds one
datum: `[M]` on that d=3 homogeneous reflective box, `max_inner` 200 → 3000
buys `|k − k_inf|` `2.9e-12 → 2.0e-15` and cuts the outer count 8 → 3, so a
derived budget would also *pay for itself in outer iterations*, not just in
accuracy.

---

## 5. Refuted, and what I could not determine

Each with its structural reason (per `process-discipline`).

### Refuted

* **"`tests/sn/solve/test_d3_admission.py` is still riding a truncation."**
  REFUTED — neither `test_d3_pure_absorber_per_ordinate_psi_exact` nor
  `test_d3_scattering_infinite_medium_matches_multigroup_balance` appears in
  the escalated delta OR in the entry census. `11e78430` closed both.
* **"`test_sn_adjoint_certification`'s SPHERE run exhausts 500 outers"**
  (issue #340, comment 2). REFUTED at HEAD `4bcce0bd`:
  `TestP13KEquality::test_sphere_k_equality` is collected under
  `-m "not slow"` (`--collect-only`) and passes under escalation, and `[M]`
  the unmutated adjoint on this module's fixtures converges in **5** outers.
  The `|dk| = 3.331e-16` signature belongs to the MUTATED heterogeneous-slab
  adjoint. Structural reason the misattribution is easy: the number is just
  the FP floor of `|dk|` for `k ≈ O(1)`, so it is not distinctive.
* **"Raising the budget will fix D6."** REFUTED by measurement, not by
  argument: `k_mut` is **bit-identical** at `max_outer` 500 and 2000.
  Structural reason: `k_mut = −0.652` is a *negative* dominant eigenvalue,
  so the power iteration's flux increment sign-alternates and
  `dphi < flux_tol` is un-satisfiable while `dk` sits at the FP floor;
  `SNSolver.converged` (`solver.py:1564`) needs both.
* **"The 90-row inner census is 90 findings."** REFUTED — 44 rows are my own
  instrument's false positives (empty GMRES residual history = a
  zero-iteration *converged* exit, not a truncation), and of the remaining
  46, an early-outer truncation under a converged outer is a legitimate
  inexact-Newton posture. Structural reason: `_claims_convergence` returns
  False for an empty history, which is the same predicate the certificate
  uses to mean "no claim was made".
* **`tests/sn/acceleration/test_dsa_rate.py` as a delta member** (the brief
  named it as a known deliberate case). REFUTED — it does not appear in
  either instrument. Its `max_inner=50` rows are S2/P1 **exactness** gates
  that assert the solve converges in 2–3 inners (`res[-1] < 1e-13`); 50 is a
  generous cap, not a starvation. Its one genuinely truncating row
  (`test_accelerated_decade_count_bounded_on_grid`, `max_inner=200`) reads
  the residual *history* and stops at a 10-decade reduction — `[M]` it never
  warned, so those solves do reach `inner_tol`.

### Not determined — with what would settle each

* **`@pytest.mark.slow` rows (114 deselected).** Structurally outside both
  instruments under the canonical `-m "not slow"` (lessons B4), and they are
  the *most* likely to carry big budgets. Settled by one census run with
  `-m slow`; not attempted here because the brief scoped the slice and the
  run cost is unbounded.
* **The WINDOWED SI arm.** `solver.py:1741` skips
  `_certify_within_group_exit` when `windowed`, so census B cannot see a
  truncated 2-D moment-windowed inner. Settled by wrapping
  `SourceIteration.solve` (or `si.solve`'s return) instead of the
  certificate — a one-line plugin change plus a re-run.
* **Whether any of the 46 genuine inner truncations lands on the FINAL
  outer.** Census B records per-call, not per-position, so "biases the
  converged k" cannot be read off it directly. Two rows were settled by
  direct consequence measurement (§3.9, both benign); the rest are
  unmeasured. Settled cheaply by extending the plugin to record the LAST
  certificate call per `solve_*` invocation rather than every call.
* **Non-`tests/sn` trees.** `tests/cp`, `tests/moc`, `tests/diffusion`,
  `tests/derivations` were out of scope by instruction. The `derivations/`
  family is separately known to carry the same defect class (issue #340,
  comment 1: `CriticalSolution.converged: bool = True`), and no SN
  instrument reaches it.

