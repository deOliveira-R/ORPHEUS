# SN Regression-Test Tolerance Redesign + Drift Tripwire (✅ DONE 2026-06-01)

**Landed:** `94d3aa4` (keff_history fixture) · `887f2de` (principled gate) ·
`7851eb2` (regenerate snapshots). Zero production changes. Regression tier =
**11 passed** under normal / `-O` / strict `-W error::DriftWarning`.

**As-built:** `tests/sn/regression/_regression_assert.py::assert_regression`
(single source of truth, `-O`-safe — np.testing + explicit raise). Iterative
gate = `SAFETY × conv_tol` with conv_tol READ off the run config (k_eff→keff_tol
1e-12, eigen flux→flux_tol 1e-10, fixed-source→inner_tol 1e-12); **SAFETY=10**
derived from the iteration-map gap bound `ρ/(1−ρ) ≲ 10` for the roster's `ρ≲0.9`
(principled, not a fudge). Direct branch = `assert_array_almost_equal_nulp`
(present; no current case is direct). `DriftWarning` tripwire: informational by
default (summarized), `-W error::tests.sn.regression._regression_assert.DriftWarning`
escalates to a strict bit-identity gate. Every regenerated snapshot was
corroborated against a structurally-independent reference BEFORE trust
(closed-form k_inf=λ_max(A⁻¹F), particle balance, flat inf-medium flux) — all
reldiffs ≤ ~5e-12.

**Two registry defects found + fixed (test-only):** (1) `slab/sphere_2g_p1_aniso`
were MALFORMED eigen problems (non-multiplying mixture → k=0/0 → snapshots stored
all-NaN, never valid coverage) → reformulated to well-posed **fixed_source** P1
cases (preserves the P1-Galerkin coverage; corroborated vs flat inf-medium flux).
(2) `2d_1g_LS4` raised NotImplementedError (2-D SI deferred) → switched to
`inner_solver="krylov"` (reaches k_inf=1.5 exactly). `test_sweep_cache` stale
`keff_history>=5` (degenerate slab converged in 3) → replaced with a
heterogeneous fuel|mod|fuel 2G slab needing ~7 outer steps (strictly stronger).

**OPEN (user decision, low-stakes):** the P1-aniso eigen→fixed-source
reformulation is a scope expansion forced by the NaN defect. Accepted as the
better option (the eigen form was broken-from-creation; reformulation preserves
P1 coverage). Drop-instead = a one-line registry change if preferred.

---

## (original plan below — design now realized as above)


**Why:** `tests/sn/regression/test_dd_regression.py` uses MAGIC FLOORS
(`rtol=1e-12`/`atol=1e-13` slab; `rtol=5e-6` curvilinear) and a docstring that
claims "bit-for-bit". Diagnosis (2026-06-01):

- **Slab cases drift ~1e-11** (k_eff rel diff 3.3e-12 / 1.4e-11) — past the
  1e-12 floor. Pure FP-non-associativity from the Wave-T reduction-tree change.
  The floor is too tight for a multigroup eigenvalue iteration. NOT a regression.
- **Curvilinear + 2-D cases drift 40–93%** — the snapshots are STALE, frozen
  from the broken curvilinear-k_inf state the docstring itself flagged as "under
  investigation (`test_kinf_homogeneous` failure cluster)". That cluster is NOW
  FIXED + PASSING: `test_kinf_homogeneous` (k_inf=νΣ_f/Σ_a=1.875 exact) +
  `test_keff_curvilinear` = **55 passed** (structurally-independent confirmation).
  So the current solver is verified-correct; the snapshots are stale-from-broken.
  Regenerating bakes the CORRECT output (good news: solver was fixed, not regressed).

**User decision (2026-06-01):** drop the magic floors; regenerate ALL snapshots
fresh from the verified-correct solver; put the whole suite on ONE defensible,
principled tolerance; add a drift tripwire.

## Principled tolerance (no magic numbers)

Tie tolerance to what the computation *promises*:
- **Iterative results** (k_eff, flux from power-iteration / BiCGSTAB):
  `tol = SAFETY × solver_convergence_tol`. The solver only converged to its own
  stopping criterion — tighter is unphysical. Read `conv_tol` off the actual
  solver config (SourceIteration / KEigenvalue / Krylov), don't hardcode.
- **Direct results** (single sweep, no outer iteration):
  `np.testing.assert_array_almost_equal_nulp(actual, expected, nulp=reduction_depth)`
  — the FP-non-associativity bound (reduction_depth ≈ #cells or #ordinates summed).

## Drift tripwire — pytest WARNING (signal) + logging (detail)

`DriftWarning` (a custom `UserWarning` subclass) is the tripwire, NOT logging:
- pytest **summarizes warnings** at run end (visible, deduped) — log lines get buried.
- **Escalatable**: `-W error::DriftWarning` makes ANY sub-tolerance bit-drift a HARD
  fail — so a "pure refactor, zero numerical change" PR runs with it and the
  tripwire becomes a strict bit-identity gate. Default run: informational.
- `logging.debug` is the FORENSIC layer (where the project under-uses logging):
  per-cell ULP breakdown + which array + the reduction-depth bound, surfaced via
  `--log-cli-level=DEBUG` when a DriftWarning fires.

## Reusable helper (single source of truth — coding-elegance Pattern 5)

```python
def assert_regression(actual, expected, *, conv_tol, case_name, kind):
    # 1) hard-fail past the PRINCIPLED correctness tolerance:
    #    iterative -> SAFETY*conv_tol ; direct -> nulp(reduction_depth)
    # 2) else if moved beyond bit-identity (ULP>0): warnings.warn(DriftWarning(
    #       f"{case_name}: drifted {n_ulp} ULP / {rel:.2e} rel (< tol {tol:.2e})"))
    #    + logging.debug(per-cell breakdown)
    # 3) else: silent (bit-identical)
```
All regression tests call it; no per-test magic numbers. Register `DriftWarning`
in `conftest`/`_harness`; document `-W error::DriftWarning` as the strict-refactor gate.

## Execution gating
1. WAIT for the sentinel-harness agent to land (it's editing `tests/sn` +
   `pyproject` markers now — the regression rewrite + snapshot regen would race it).
2. Run the **xdist pin-test** in that same free window (pip-install would race the
   agent's cosmic-ray/mutmut install): test pytest-xdist <3.8 on Py3.14.3; if a
   version avoids the 3.8.0 scheduler deadlock, restore `-n` as the parallel runner
   + update the pyproject guidance; else keep per-tier default.
3. Build `assert_regression` + `DriftWarning`; rewrite `test_dd_regression` onto it;
   `python -m tests.sn.regression._generate_snapshots` to regenerate ALL snapshots
   from the verified-correct solver; confirm green; commit.
4. The 1 `test_sweep_cache::test_collision_cache_invariance` stale `keff_history>=5`
   fixture (trivial slab converges in 3) — fix in the same pass (assert >=3, or
   pick a case that needs ≥5).
