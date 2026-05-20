# R-1 Step D — sphere `_solve_krylov` oscillation traced to sweep-as-preconditioner

**Date**: 2026-05-19
**Branch**: `refactor/sn-operator-algebra`
**Commit context**: between Step C (`2c634ab`) and Step D (in progress)

## Symptom

`SNSolver._solve_krylov` carved onto `KrylovAcceleration` (R-1 Step D)
worked on slab (5/5 L1 tests in 0.77s; `keff = 1.875` exact) but
oscillated wildly on sphere:

```
n=0: keff=2.0291238898, dk=1.03e+00
n=1: keff=1.4281719749, dk=6.01e-01
n=2: keff=2.0549067555, dk=6.27e-01
n=3: keff=1.4439832385, dk=6.11e-01
n=4: keff=2.2456821886, dk=8.02e-01
```

Reference `k_inf = 1.875`. The inner GMRES was converging to a
*flux* the outer eigenvalue iteration could not stabilise around.
Per-inner-call time was 12.6 s vs 0.033 s for slab (470× slowdown).

## Root cause — `LC.solve` (the WDD sweep) is not a consistent inverse of `LC.apply` when invoked as a Krylov preconditioner on curvilinear geometry

The carve passes `preconditioner=None` to `KrylovAcceleration`. The
default branch in `orpheus/numerics/iteration.py:596–601` then uses
`L.solve` when `L` advertises `CAP_SOLVE`. Step C's
`InvertibleOperator` advertises `CAP_SOLVE` (it IS sweep-invertible by
construction). So the carve unintentionally runs GMRES with the WDD
sweep as the preconditioner.

On curvilinear geometry the discrete `(L + C)` is **rational in
σ_t** through the Carlson coupled-pole seed (Hébert §3.9.4
Eqs. 3.432–3.435; also `StreamingOperator` docstring at
`orpheus/sn/operator.py:1801–1815`). The Carlson seed is derived from
the **previous-iterate scalar flux** `φ₀(prev)` — threaded into
`InvertibleOperator.solve` via `rhs(1)` (the lag-1 frame of an
`AngularFlux` with attached history). `R-1 Step A` built the typed
history container precisely for this purpose; `R-1 Step B`'s
`SourceIteration.solve` calls `_attach_previous(rhs_raw, psi_prev)` to
attach the previous iterate to the rhs before calling `L.solve`.

But **GMRES feeds the preconditioner a residual vector**, not the
iterate. The residual has no `(1)` history. `InvertibleOperator.solve`
reads `rhs(1) → None` (cold-start branch) and `transport_sweep` falls
back to its in-iteration-source default Carlson seed — which is
**not** the algebraic inverse of `LC.apply` on the discrete
curvilinear operator. The preconditioned operator `A · M⁻¹` is
therefore far from identity; GMRES diverges.

Slab is unaffected because slab's discrete `(L + C)` is **affine** in
σ_t (no Carlson seed needed — Cartesian residual is zero per the same
docstring). The sweep IS an exact inverse of the matvec for slab,
GMRES converges in 47 iters, and the slab L1 passes are genuine
evidence of end-to-end correctness, not coincidence.

## Confirming evidence

`derivations/diagnostics/diag_r1_step_d_probe_b_identity_precond.py`
toggles the preconditioner between the default (sweep) and identity
`M = I`, and prints `(coord, n, keff, dk, inner_iters)`:

| coord  | preconditioner    | n=0 keff | n=1 keff | n=2 keff |
|--------|-------------------|---------:|---------:|---------:|
| slab   | default (sweep)   | 1.875    | 1.875    | 1.875    |
| slab   | identity          | 1.875    | 1.875    | 1.875    |
| sphere | default (sweep)   | 2.029    | 1.428    | 2.055    |
| sphere | identity          | 1.868    | 1.875    | 1.875    |

Sphere with identity preconditioner converges to `keff = 1.875` at
`rtol = 1e-10` within 6 outer iterations (≈ 60 s wall-clock).

## Fix

`orpheus/sn/solver.py` (R-1 Step D `_solve_krylov`): pass an explicit
identity preconditioner instead of `None`:

```python
krylov = KrylovAcceleration(
    LC, S_normalised, ZeroOperator(),
    preconditioner=lambda q: q,  # explicit identity — #200 tracks re-enable
    tol=self.inner_tol, max_iter=self.max_inner,
    restart=min(50, N * ng * nx * ny),
)
```

The user direction for R-1 was "no preconditioner" — the
`preconditioner=None` argument was meant to express that intent, but
`KrylovAcceleration.__init__`'s capability fallback silently routed
through `LC.solve`. The explicit lambda makes the intent
load-bearing and prevents the silent fallback.

Long-term: `KrylovAcceleration` may want a sentinel like
`preconditioner="identity"` so the call site reads explicitly.

## Issue #200 connection

The block-inverse face preconditioner (`σ_r = σ_t - Σ_{s,0}^{g→g}`
folded into the diagonal) is tracked by
[issue #200](https://github.com/deOliveira-R/ORPHEUS/issues/200).
That preconditioner does NOT require the Carlson seed (the
foldable-scattering preconditioner is built directly on
`InvertibleOperator(L, CollisionOperator(σ_r))` — a different
algebraic identity from the full sweep). It will re-enable Krylov
acceleration on curvilinear once #200 lands.

## Regression catch

Promote `diag_r1_step_d_probe_b_identity_precond.py` to a permanent
L1 test (e.g. `tests/sn/test_krylov_curvilinear_preconditioner.py`)
with an assertion that sphere reaches `|keff - 1.875| < 1e-6` after 3
outer iterations under identity preconditioner. If a future change
re-enables `L.solve` as default precond on curvilinear, the test
catches the oscillation immediately.

## Files referenced

- `orpheus/sn/solver.py:~608` (the bug site — `preconditioner=None`)
- `orpheus/numerics/iteration.py:596–601` (the silent
  `None → L.solve` fallback)
- `orpheus/sn/operator.py:2570` (`carlson_seed = rhs(1)` —
  returns `None` for residual vectors)
- `orpheus/sn/operator.py:1801–1822` (Hébert §3.9.4 docstring —
  the structural reason curvilinear `(L + C)` is rational in σ_t)
- `derivations/diagnostics/diag_r1_step_d_probe_b_identity_precond.py`
  (the diagnostic probe)
