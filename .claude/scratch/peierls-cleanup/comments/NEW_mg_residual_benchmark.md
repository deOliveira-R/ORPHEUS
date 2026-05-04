<!-- New issue body — to be filed 2026-04-30 as part of the post-#138 Peierls cleanup -->
<!-- Source: docs/theory/peierls_nystrom.rst lines 1313-1326 (~13 LoC of existing prose) -->
<!-- Suggested labels: module:peierls, level:L2, type:improvement -->
<!-- Suggested title: "MG: benchmark Peierls multi-group residual against discrete CP modules" -->

## MG: benchmark Peierls multi-group residual against discrete CP modules

**Summary (2026-04-30):**

1. Issue #104 commit 2 (2026-04-24) shipped six 2-group hollow Peierls reference cases (cyl/sph at $r_0/R \in \{0.1, 0.2, 0.3\}$) registered through `_class_a_cases`, but **the L1 parity gate against the discrete-CP MG solvers (`cp_cylinder`, `cp_sphere`) was deferred and never executed**.
2. The cases are **buildable and reproducible** today (gated by `TestMGRegistration` smoke tests), but their MG residuals against an independent reference are **not benchmarked**.
3. This issue tracks running that parity gate (1 % target per Issue #104 AC, matching the existing 1G hollow residual class) so the MG references graduate from "buildable" to "L1-verified".

---

## Context

As of 2026-04-24 (Issue #104 commit 2), `orpheus.derivations.peierls_cases._class_a_cases` registers the following 2-group hollow cells alongside the legacy 1G entries:

| Reference name                          | $r_0 / R$ | Closure                                           |
|-----------------------------------------|----------:|---------------------------------------------------|
| `peierls_cyl1D_hollow_2eg_1rg_r0_10`    |       0.1 | F.4 scalar rank-2 per-face (Ki₃ fold)             |
| `peierls_cyl1D_hollow_2eg_1rg_r0_20`    |       0.2 | same                                              |
| `peierls_cyl1D_hollow_2eg_1rg_r0_30`    |       0.3 | same                                              |
| `peierls_sph1D_hollow_2eg_1rg_r0_10`    |       0.1 | F.4 scalar rank-2 per-face (bare $e^{-\tau}$ kernel) |
| `peierls_sph1D_hollow_2eg_1rg_r0_20`    |       0.2 | same                                              |
| `peierls_sph1D_hollow_2eg_1rg_r0_30`    |       0.3 | same                                              |

All six use the single-region fuel composition (`LAYOUTS[1] = ["A"]`, `get_xs("A", "2g")`) with region-A's 2G XS. Registry is lazy (built on first access to `reference_values.continuous_all` / `continuous_get`); each reference builds on demand at ~1–2 min wall time at default quadrature, so the first access after import is expensive.

The L1 residuals vs an independent MG reference solver have **not** been benchmarked as part of Issue #104 commit 2 — the parity gate against `cp_cylinder` / `cp_sphere` 2G native solvers is the planned follow-up that this issue tracks. The registered cases are buildable and reproducible as established by the smoke tests in `tests.derivations.test_peierls_multigroup.TestMG2GHollowRegistration`.

## What needs benchmarking

For each of the six 2G hollow references above, run an L1 parity gate of the form:

1. Build the Peierls 2G reference $k_{\rm eff}$ and group-flux profiles via the lazy registry.
2. Solve the **same 2G hollow geometry** with the discrete-CP module (`orpheus.cp.cp_cylinder` for cyl, `orpheus.cp.cp_sphere` for sph) at converged spatial discretization.
3. Compute the relative-difference vector of $k_{\rm eff}$ and per-region group-flux residuals.
4. Gate at **1 % rel_diff** (matching the Issue #104 AC and the existing 1G hollow residual class).

Tests should live in a new module `tests/derivations/test_peierls_mg_vs_cp_residual.py` (or extend the existing `test_peierls_multigroup.py` with a `TestMG2GHollowVsDiscreteCP` class). Mark each test with `@pytest.mark.l2` (`level:L2` integration: multi-group + heterogeneous self-convergence) and `@pytest.mark.verifies("peierls-equation")` so the parity gate appears in the V&V matrix.

## Why this is the right gate

Per Cardinal Rule 5, **multi-group (≥2G) is mandatory** for non-degenerate verification: 1-group tests are degenerate ($k = \nu\Sigma_f / \Sigma_a$ regardless of flux shape). The Peierls 1G hollow cases already pass an L1 parity vs the discrete CP modules at 1 % rel_diff (the "1G hollow residual class"); the 2G cases need the same gate before they can claim L1 status. Until then they are at most **L0** (smoke-built but not gated).

Discrete-CP-vs-Peierls is the right reference because:

- Both solvers compute volume-weighted scalar fluxes at the same set of region indices.
- Both are spatially-resolved; agreement at 1 % rules out implementation bugs in either solver's MG kernel.
- The 1 % target is a known-achievable class (1G already passes); a tighter target risks gating on the residual quadrature error of the cheaper of the two solvers, not on physics.

## Acceptance criteria

- [ ] New test module (or class) with one parametrized test per (geometry, $r_0/R$) — six tests total.
- [ ] Each test asserts:
  - $|k_{\rm eff}^{\rm Peierls} - k_{\rm eff}^{\rm CP}| / k_{\rm eff}^{\rm CP} < 0.01$
  - $\max_{r,g}\,|\phi_{r,g}^{\rm Peierls} - \phi_{r,g}^{\rm CP}| / \phi_{r,g}^{\rm CP} < 0.01$
- [ ] All six tests are tagged `@pytest.mark.l2` and `@pytest.mark.verifies("peierls-equation")`.
- [ ] If a case fails the 1 % gate, the failure is investigated to root cause (likely either Peierls quadrature undersampling or CP spatial undersampling); raise discretization until the gate passes or file a child bug if the discrepancy is structural.
- [ ] `python -m tests._harness.audit` shows the six new entries in the V&V matrix at level L2.
- [ ] Sphinx page `docs/theory/peierls_nystrom.rst` "Shipped multi-group references" subsection updates the placeholder text ("L1 residuals not yet benchmarked") to "L1 residuals gated at 1 % vs discrete-CP MG; see `tests/derivations/test_peierls_mg_vs_cp_residual.py`".

## Cross-links

- Predecessor: Issue #104 (Peierls multi-group extension, commit 2 shipped six 2G hollow cases without parity gate)
- Cleanup commit (origin of this issue): <!-- COMMIT-HASH --> (Peierls docs cleanup, post-#138, 2026-04-30)
- Surviving Sphinx anchor: `theory-peierls-multigroup` "Shipped multi-group references" subsection now reads "Multi-group residual benchmarks vs the discrete CP modules are tracked separately under Issue #NNN" (this issue).
- Existing smoke tests: `tests.derivations.test_peierls_multigroup.TestMG2GHollowRegistration`.
- Reference solvers to compare against: `orpheus.cp.cp_cylinder`, `orpheus.cp.cp_sphere` (discrete CP MG drivers).
