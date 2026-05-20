---
name: r1-step-e-invertible-solve-w-bridge
description: R-1 Step E SI carve diagnosis — InvertibleOperator.solve had no ×W convention bridge to transport_sweep, double-divided rhs by W, broke every ng≥2 kinf_homogeneous SI test
metadata:
  type: project
---

# R-1 Step E — `InvertibleOperator.solve` ↔ `transport_sweep` `W` convention bridge

**Date**: 2026-05-19
**Branch**: `refactor/sn-operator-algebra`
**Symptom commit context**: between Step D (`8ff0881`) and Step E (in-progress)

## Symptom

The Step E carve of `_solve_source_iteration` onto
`SourceIteration.solve(q_ext_AngularFlux)` converged to the wrong
fixed point for every ng≥2 case across all coords:

* `slab-2eg-source_iteration`: keff = 1.4843750000 vs ref 1.875000 (drift 0.39)
* `sphere-2eg-source_iteration`: same drift
* `cylinder-2eg-source_iteration`: same drift
* same 4eg signature on every coord
* 1eg passes (degenerate — k = νΣ_f/Σ_a is shape-independent and
  the bug only affects the per-ordinate ψ magnitude, not the
  inferred k from the Rayleigh quotient when scattering is
  symmetric)

The Krylov path (Step D) worked perfectly on the SAME mesh/material
combinations.  Diagnostic-cascade clue: ratio `1.4844/1.875 ≈ 0.7917`
is not a clean `1/W = 0.5`, but a function of `(W² Σ_t - W Σ_s)`
vs `(Σ_t - Σ_s)`.

## Root cause — convention mismatch between operator-algebra and `transport_sweep`

`InvertibleOperator.solve(rhs)` calls `transport_sweep(iso_source=0,
aniso_source=rhs.values, ...)`.  But `transport_sweep` expects the
`aniso_source` in **iso-source magnitude** — it applies the internal
`weight_norm = 1/W` to convert iso → per-ordinate (see
`orpheus/sn/sweep.py:432`: `Q_aniso_p = Q_aniso * weight_norm`).

The operator-algebra layer (`SourceIteration`, `KrylovAcceleration`)
operates in **per-ordinate units** everywhere:

* `S_normalised = self.scattering_op / sum_w` → `S.apply(ψ)` returns
  per-ordinate scatter source
* `q_ext = fission_source / sum_w` (already per-ordinate)
* `rhs = S.apply(ψ) + q_ext` → per-ordinate

So `rhs.values` is per-ordinate, but `transport_sweep` divides it by
`W` AGAIN as if it were iso-magnitude.  Effect: per-ordinate balance
sees source `(Sψ + q)/W²` instead of `(Sψ + q)/W`.

Krylov was unaffected because it uses GMRES with `preconditioner=
lambda q: q` (identity) — it never calls `LC.solve`, only `LC.apply`
(matvec) which operates in consistent per-ordinate units.

## The fix

`orpheus/sn/operator.py` `InvertibleOperator.solve`:

```python
sum_w = float(sn_mesh.quad.weights.sum())
iso_source = sn_mesh.zeros_isotropic_source()
aniso_source = PerOrdinateSource(rhs.values * sum_w, sn_mesh)  # ×W bridge
```

Multiply `rhs.values` by `sum_w` before constructing the aniso
source.  The sweep's internal `/W` then lands in per-ordinate units,
matching the operator-algebra contract.

## Confirming evidence

`derivations/diagnostics/diag_r1_step_e_invertible_solve_w_bridge.py`
(10 tests, all PASS post-fix):

* `test_invertible_solve_slab_uniform_roundtrip` — slab
  `LC.solve(LC.apply(ψ=1)) == ψ=1` to machine zero (streaming-free
  check; sphere/cyl excluded because M-M angular sweep/matvec are
  not strict inverses even for uniform ψ).
* `test_invertible_solve_fixed_source_homogeneous_reflective[slab|sphere|cylinder]`
  — converged `LC.solve(q_per_ord)` on reflective homogeneous medium
  recovers ψ = q/(W Σ_t) per ordinate to rel-tol 1e-6.
* `test_si_carve_recovers_analytical_kinf` × {2eg, 4eg} × {slab,
  sphere, cylinder} — end-to-end SI carve reaches analytical
  k_inf at rel_err < 1e-9.

## Important nuance — what is NOT true

`LC.solve(LC.apply(ψ))` is NOT machine-zero for general ψ even
after the fix.  The WDD sweep (`LC.solve`) and WDD matvec
(`LC.apply`) are **two distinct discrete operators** that share the
same fixed point under SI/Krylov contracts.  They are NOT a strict
matrix inverse of one another.  Only the SOURCE-INTERPRETATION
convention (`rhs.values` units) needed bridging; the discrete-
operator difference between sweep and matvec is structural to the
WDD scheme itself and is the subject of Phase G's eventual full
operator-algebra unification.

For uniform ψ on slab (where streaming vanishes by symmetry) the
round-trip IS machine zero.  For sphere/cylinder even uniform ψ
gives an O(1) round-trip diff because M-M angular closure has
different cell-level discrete forms in `.solve` (sequential angular
recurrence) vs `.apply` (matvec block).

## Sphere also failed (not "secretly passing" as the user thought)

The user's brief stated "sphere-2eg-source_iteration PASSES".
That's not what the test suite shows — sphere-2eg-source_iteration
FAILED with the same ~21% drift as slab and cylinder.  All ng≥2 SI
tests failed identically.  1-group passes because k = νΣ_f/Σ_a is
flux-shape independent (the canonical 1G degeneracy).

## Bug class

ERR-049 candidate: failure mode #4 (convention drift).  Same family
as ERR-022 (signed lethargy across multiple sites), ERR-025 (1/W
normalization convention with sign cancellation).  Mechanism:
operator-algebra layer (Step C/D/E) and `transport_sweep` (Wave A-D)
were authored under different conventions; the SI carve was the
first place they interacted at the per-ordinate-source level.

## Files referenced

* `orpheus/sn/operator.py:~2555-2587` (InvertibleOperator.solve — fix
  site, ×sum_w bridge applied)
* `orpheus/sn/solver.py:~506-525` (_solve_source_iteration R-1 Step
  E carve — consumer of the operator algebra)
* `orpheus/sn/sweep.py:~427-432` (`weight_norm = 1.0 / weights.sum()`
  — the internal /W the bridge compensates)
* `orpheus/sn/scattering.py:~798-799` (docstring noting "the sweep's
  ``1/W`` factor is not applied here — the caller, ``transport_sweep``,
  applies it")
* `derivations/diagnostics/diag_r1_step_e_invertible_solve_w_bridge.py`
  (pinning test — 10 cases)

## Recommendation for promotion

* `test_si_carve_recovers_analytical_kinf` is subsumed by the
  existing `tests/sn/l1_analytical/test_kinf_homogeneous.py::test_kinf_homogeneous[…-source_iteration]`
  parametrisation — once R-1 Step E lands, those 6 tests catch the
  same regression.  Leave the diagnostic version in place as a
  faster cross-check.
* `test_invertible_solve_fixed_source_homogeneous_reflective` is
  NEW coverage — promote to `tests/sn/test_operator.py` under a new
  `TestInvertibleOperator` class, retags as L0 streaming-equilibrium
  (covers L0-SN-001 family for the operator-algebra layer that
  doesn't have a separate L0 test yet).  This is the canonical L0
  probe for the InvertibleOperator class.
* `test_invertible_solve_slab_uniform_roundtrip` is also new L0
  coverage of the convention bridge — promote to the same test
  file as a streaming-free roundtrip check (renamed appropriately).
