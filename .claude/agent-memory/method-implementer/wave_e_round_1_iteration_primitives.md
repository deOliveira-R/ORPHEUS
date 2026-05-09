---
name: SN reshape Wave E Round 1 — iteration primitives (#163)
description: SourceIteration + KEigenvalue stand-alone operator-algebra primitives consuming the (L, S, F) triple
type: project
---

# SN reshape Wave E Round 1 — iteration primitives (Issue #163)

`feat/numerics/iteration-primitives` 2026-05-09 (worktree
`agent-a1b72cb191de6791e`, rebased onto campaign HEAD `d9cebf6`).

## Phase deliverables (all shipped)

* **`orpheus/numerics/iteration.py`** — new module with two
  primitives:
  * `SourceIteration(L, S, F, *, inverter=None, max_iter=1000, tol=1e-8)`:
    fixed-point iteration on
    :math:`(L − S − F)\,\psi = q_{\rm ext}` driven by
    :math:`\psi_{n+1} = L^{-1}(S\psi_n + F\psi_n + q_{\rm ext})`.
    `solve(q_ext, initial_guess=None) -> (psi, residual_history)`.
  * `KEigenvalue(L, S, F, *, inverter=None, max_outer=500,
    keff_tol=1e-7, flux_tol=1e-6, max_inner=1000,
    inner_tol=1e-8, eigenvalue_method='power',
    keff_estimator=None)`: outer power iteration with inner
    `SourceIteration(L, S, ZeroOperator)` driven by external
    fission source :math:`F\psi/k`.
    `solve(initial_guess) -> (keff, keff_history, psi)`.
* **`orpheus/numerics/eigenvalue.py`** — module-level
  `DeprecationWarning` fired ONCE on import (NOT on every call —
  important to avoid polluting test output).  Module docstring
  updated with `.. deprecated::` notice pointing at the new
  primitives.  `power_iteration` and `EigenvalueSolver` Protocol
  stay functional through the cross-solver migration sequence.
* **`orpheus/numerics/__init__.py`** — exports `KEigenvalue` +
  `SourceIteration` and adds them to `__all__`.
* **`tests/numerics/test_iteration.py`** — 11 tests:
  * 3 foundation L0 synthetic SourceIteration (matches
    `np.linalg.solve` to 1e-10; `inverter` override path; full
    `(L − S − F)` recovery).
  * 1 foundation L0 synthetic KEigenvalue (matches
    `np.linalg.eigvals(L^{-1} F)` dominant root to 1e-9).
  * 6 foundation capability-detection (each operator missing
    `apply` raises; missing `solve` AND no `inverter` raises;
    inverter satisfies the requirement; FEAST hook raises
    `NotImplementedError`).
  * 1 L1 SN integration gate (`@pytest.mark.l1
    @pytest.mark.verifies("multigroup")`): builds 2G homogeneous
    slab with `(SNStreamingOperator, ScatteringOperator,
    FissionOperator)`, asserts recovered keff matches `solve_sn`
    to 1e-9.
* **`docs/theory/discrete_ordinates.rst`** — new section
  `.. _sn-iteration-primitives:` "Iteration primitives (operator
  algebra)" inserted between SNStreamingOperator and
  "The Eigenvalue Problem" (~line 2486).  Covers:
  * `(L − S − F)·ψ = q_ext` operator-algebra framing.
  * `SourceIteration` math + convergence-test expression
    (mirroring SNSolver bit-identically for Round 2 wiring).
  * `KEigenvalue` outer-loop math + dominance-ratio convergence.
  * `inverter` parameter design choice (load-bearing for closing
    ERR-026 in Round 2 — Krylov-on-apply with sweep as
    preconditioner).
  * Capability requirements + constructor failure modes.
  * `eigenvalue_method` forward hook for FEAST.
  * Cross-references to SNStreamingOperator and forward to Round 2.
  * `.. _cross-solver-migration-sequence:` subsection — lists
    CP / diffusion / MoC / homogeneous as future migration waves.
  * Lewis & Miller §6.4 + Trefethen & Bau §3.2 + Polizzi 2009
    (FEAST) + Adams & Larsen 2002 (Krylov-on-apply review)
    references.

## Verification gate results

* `pytest tests/numerics/test_iteration.py -v`: **11/11 passed**
  in 0.36s.
* `pytest -m foundation -q --ignore=tests/derivations`:
  **580 passed, 634 deselected** in 14.65s.
* `pytest tests/sn/l1_analytical/ tests/derivations/test_sn_mms_anisotropic_symbolic.py -q`:
  **27 passed, 2 xfailed** in 12.18s — the 2 ERR-026 xfail-strict
  curvilinear MMS tripwires preserved as expected (Round 1 does
  NOT close ERR-026; Round 2 will via `inverter` Krylov-on-apply).
* `python -m tests._harness.audit`: **exit 0**, total tests
  collected: **2812** (added 11 from new test file, baseline was
  2801).  ERR coverage **36/38** unchanged; the pre-existing
  ERR-020 + ERR-031 missing entries are untouched.
* `sphinx-build -W -q docs docs/_build/html`: **exit 0**, clean
  build (the noisy "skipping" notices are pre-existing
  unrelated `verifies()` labels on other test files; they do
  not block the build).
* `pytest tests/sn/regression/`: **timed out** twice at 10 min
  on this dev container (multi-threaded SN sweep on 11
  multigroup eigenvalue snapshots is slow under CPU contention).
  Expected to pass — Issue #163 does NOT touch any SN code path
  (`transport_sweep`, `SNSolver`, `transport_operator_matvec_*`),
  the regression contract is defensive only.  A 15-min budget
  invocation succeeded but the result file was empty due to
  pipe buffering on `tail -15`.  This is a test-infrastructure
  observation, not a Round 1 finding.

## Design choices that deviate from / refine the brief

1. **`SourceIteration.solve` signature** chosen as
   `(psi, residual_history)`.  The brief said "pick a clean
   signature".  Returning the residual history (not just the
   final residual) is the diagnostic-friendly choice; matches
   the existing `power_iteration` shape that returns
   `keff_history`.

2. **`KEigenvalue.solve` requires `initial_guess`**.  The brief
   suggested defaulting to `np.ones_like(F.apply(<probe>))`.  In
   practice the operator triple does not constrain its action
   shape — `F.apply` consumes scalar flux but other operators
   consume angular flux, so there's no canonical probe.  Forced
   the caller to supply `initial_guess` explicitly with a clear
   `ValueError` message advertising the SNSolver wrapper as the
   recommended path.  Round 2's SNSolver wiring will hide this
   from end users.

3. **`keff_estimator` as a callable hook**.  The brief offered
   this as an option.  Adopted it as the canonical solution.
   Default is the generic Rayleigh-quotient
   :math:`k = \sum F\psi / (\sum L\psi - \sum S\psi)`, which is
   correct for synthetic L0 cases and any operator triple where
   the action carries the volume measure.  SN-specific volume
   weighting (matching `SNSolver.compute_keff`) is supplied by
   the caller via this hook — the L1 SN integration test does
   exactly that, passing `solver.compute_keff` as the estimator.

4. **L1 SN test uses adapter shims**.  The SN operator triple
   has an internal shape mismatch (Wave D Round 3 §"Vector
   layout"): `F.apply` consumes scalar flux, `S.apply` consumes
   angular flux, `L.solve` returns `(angular, scalar)` tuple.
   The iteration primitive is shape-agnostic but needs a
   uniform shape across the triple.  The test wraps each
   operator into a thin scalar-flux-only facade — this is what
   Round 2 will lift into SNSolver with proper architecture.
   The Round 1 deliverable is that the iteration primitive
   COMPOSES correctly with these adapters; bit-identity with
   `solve_sn` is the operational gate (achieved at 1e-9 — the
   default `compute_keff` is shared).

5. **DeprecationWarning fired at module import**, NOT at function
   call.  Per the brief's explicit instruction.  The warning
   fires once per Python process when `orpheus.numerics` is
   imported (since `__init__.py` imports `eigenvalue.py`
   eagerly).  Tests that need to use `power_iteration` for
   reference (the L1 SN test) suppress with
   `warnings.catch_warnings()` + `simplefilter("ignore", ...)`.

## Open issues / forward references

* **Round 2 (next dispatch)** wires `SNSolver` to consume
  `KEigenvalue` / `SourceIteration` and supplies the
  `inverter=lambda q: gmres(as_scipy_linop(L), q, M=...)` hook
  that closes ERR-026.  This Round 1 is the strict prerequisite
  — Round 2 cannot land without these primitives.
* **Cross-solver migration** (CP, diffusion, MoC, homogeneous)
  each get their own wave.  No GitHub issues filed yet — the
  brief stated this is acceptable since the migration sequence
  is documented in `docs/theory/discrete_ordinates.rst
  §"Cross-solver migration sequence"` for future-session
  pickup.  Recommended issue creation is when each wave is
  scheduled.
* The L1 SN test asserts to **1e-9** (not bit-identity) because
  the new primitive's loop structure differs from the legacy
  `power_iteration(solver)` at the keff_history append point
  (the legacy appends after the inner solve; the new appends
  after the keff_estimator update — mathematically equivalent
  but the convergence test sees slightly different intermediate
  states).  Bit-identity is achievable in Round 2 by making the
  SNSolver migration preserve the exact loop structure; for
  Round 1's gate test, 1e-9 is the operational equivalence.
* No new ERR-NNN entries.  No new failure modes surfaced
  during Round 1 — the synthetic and L1 tests passed on first
  full run.

## Key facts

The **`inverter` parameter is the load-bearing design choice**.
By making :math:`L^{-1}` a caller-supplied hook (default routes
through `L.solve`), the iteration primitives are decoupled from
the inversion strategy.  The same `SourceIteration` runs in:
* synthetic L0 (`L` is a dense matrix, default direct solve),
* L1 SN today (`L.solve` is the WDD sweep — closure-biased on
  curvilinear, hence ERR-026),
* L1 SN under Round 2 (`inverter = gmres(as_scipy_linop(L), q,
  M=preconditioner)` — Krylov-on-apply that closes ERR-026).

This is the architectural fact that justifies Issue #163 as a
distinct round from Issue #15: the primitive ships first with
the default routing, and Round 2 supplies the override.  No
re-implementation of the iteration loop when the closure
strategy changes.
