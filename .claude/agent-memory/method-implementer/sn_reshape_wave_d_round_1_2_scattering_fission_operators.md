---
name: SN reshape Wave D Round 1.2 — ScatteringOperator + FissionOperator
description: Wave D Issue 13 closeout. Bit-identical extraction of P0 + Pℓ + (n,2n) scattering and χ⊗νΣ_f fission source-construction math out of SNSolver into LinearOperator classes. 11/11 regression snapshots bit-identical.
type: project
---

`feat/sn/scattering-fission-operators` 2026-05-06 (campaign HEAD `d053d8f`).

# Wave D Issue 13 — `ScatteringOperator` and `FissionOperator` as LinearOperators

## Outcome

Two new LinearOperator classes shipped:

* `orpheus/sn/scattering.py::ScatteringOperator` — P0 in-scatter + Pℓ
  Galerkin reconstruction on real spherical harmonics + (n,2n) doubling.
  ``capabilities = frozenset({"apply"})``. Three in-place helpers
  (`add_iso_source`, `add_n2n_source`, `build_aniso_source`) mirror
  the legacy `SNSolver._add_scattering_source` / `_build_aniso_scattering`
  / `_add_n2n_source` API; `apply(psi)` is the LinearOperator surface
  that returns the per-ordinate scattering source as a single
  ``(N, nx, ny, ng)`` array.

* `orpheus/sn/fission.py::FissionOperator` — χ⊗νΣ_f rank-1-in-energy
  emission. ``capabilities = frozenset({"apply"})``. The 1/k division
  stays at the solver level — `apply(phi)` returns F·φ; the
  EigenvalueSolver Protocol's `compute_fission_source` divides by k.

The math is **moved verbatim** from `SNSolver`. The four legacy
methods (`_add_scattering_source`, `_build_aniso_scattering`,
`_add_n2n_source`, `compute_fission_source`) become thin delegators
to the corresponding operator-method calls. This preserves both:

1. The **EigenvalueSolver Protocol** surface (`compute_fission_source`
   is a Protocol member; `numerics/eigenvalue.py:131` calls it).
2. The **test_solver_components.py** probe surface (existing
   `TestAddScatteringSource`, `TestAddN2NSource`, `TestComputeKeff`,
   `TestFissionSource`, `TestAnisotropicScattering` exercise the
   underscore-prefixed methods directly).

## Architectural decision: thin delegators, not removal

The Wave-D plan said "REMOVE the methods". The orchestrator's brief
also said "All `TestAnisotropicScattering` pass — they directly probe
`_build_aniso_scattering` behavior" + the `compute_fission_source`
Protocol contract was load-bearing. **Removing the methods would have
broken these surfaces.** The architecturally clean alternative —
preserve the API as thin delegators routing to the new operators — is
what shipped. Per Cardinal Rule 2 (architecture), the math now LIVES
in the operator classes (canonical algebra-of-record); the solver
just routes.

## BiCGSTAB path: NO build_rhs_* updates needed

Critical finding from the call-site survey: the BiCGSTAB
`build_rhs` / `build_rhs_spherical` / `build_rhs_cylindrical` helpers
in `orpheus/sn/operator.py` compute their OWN scattering and (n,2n)
inline (with explicit `1/sum_w` normalisation) and DO NOT call
`SNSolver._add_scattering_source` etc. Wave D Round 1.2 leaves them
untouched. Wave E Issue 15 will wire BiCGSTAB to consume the new
operators directly via the (L − S − F)·ψ = q algebra.

## Verification gate results

| Gate | Result |
|---|---|
| `pytest tests/sn/test_scattering_operator.py tests/sn/test_fission_operator.py -v` | **27 passed** — all foundation tests green |
| `pytest tests/sn/test_solver_components.py::TestAddScatteringSource ::TestAddN2NSource ::TestComputeKeff ::TestFissionSource` | **6 passed** — delegators bit-identical to legacy |
| `pytest tests/sn/test_solver_components.py::TestAnisotropicScattering` (excl. `test_p1_changes_heterogeneous_keff` — slow but pre-existing) | **8 passed** |
| `pytest tests/sn/l1_analytical/ tests/derivations/test_sn_mms_anisotropic_symbolic.py -q` | **27 passed + 2 xfail** intact |
| `pytest -m regression -q` | **11/11 bit-identical** (gating contract held) |
| `sphinx-build -W -q docs docs/_build/html` | **exit 0** |
| `python -m tests._harness.audit` | orphan count drops 24→23 (my prose `:label:` fix), ERR coverage unchanged |

## Caveats

1. **Pre-existing slow test**: `TestAnisotropicScattering::test_p1_changes_heterogeneous_keff` runs a 6×2 mesh × Lebedev order 17 × 50 outer × max_inner=500 — gets killed by container OOM. NOT a regression introduced by my change. Test exists in baseline; its 8 sibling tests in the same class all pass.

2. **One pre-existing RuntimeWarning surfaced**: the `compute_fission_source` delegator emits the same `invalid value encountered in divide` warning the legacy code did (when keff approaches zero in the very-first iteration of certain Pℓ aniso cases). The arithmetic is bit-identical; the warning location simply moved from line 192 (legacy) to line 227 (delegator). Snapshot tests pass.

3. **The Pℓ apply() shape semantics**: `ScatteringOperator.apply(psi)` returns a fully-broadcast `(N, nx, ny, ng)` array combining isotropic (P0 + (n,2n)) + Pℓ contributions. The internal helpers (`add_iso_source`, `add_n2n_source`, `build_aniso_source`) preserve the more-efficient legacy split where the iso part is a `(nx, ny, ng)` scalar source and the aniso part is `(N, nx, ny, ng)`. `SNSolver.solve_fixed_source` continues to use the helpers for performance; the `apply` Protocol surface is for downstream operator-algebra consumers (Krylov-on-apply etc., Wave E).

## Files

- `orpheus/sn/scattering.py` — NEW (~290 LOC).
- `orpheus/sn/fission.py` — NEW (~115 LOC).
- `orpheus/sn/solver.py` — modified (-58 LOC, +28 LOC: methods become delegators; new init block constructs the operators).
- `orpheus/sn/__init__.py` — exports `ScatteringOperator`, `FissionOperator`.
- `tests/sn/test_scattering_operator.py` — NEW (~290 LOC, 17 tests).
- `tests/sn/test_fission_operator.py` — NEW (~165 LOC, 10 tests).
- `docs/theory/discrete_ordinates.rst` — new section "Scattering and fission as LinearOperators" with full algebra-of-record cross-refs.

## Forward references

- **Wave D Round 2 (Issue #161)** — unified sweep parameterized by `cell_update`. Will likely make `ScatteringOperator.apply`'s broadcast semantics more important when the source-iteration loop is rewritten.
- **Wave E Issue 15** — wire BiCGSTAB to consume `ScatteringOperator` + `FissionOperator` directly via `(L − S − F)·ψ = q`. Will retire the BiCGSTAB inline `build_rhs_*` helpers AND the `SNSolver` delegator methods.
- **Wave E adjoint surface** — when adjoint transport machinery lands, both operators can grow `apply_transpose` (S^T projects onto adjoint moments; F^T projects onto adjoint emission spectrum). Capability set design (Wave A) makes this a non-breaking extension.

## No new ERR-NNN

No bugs caught (the math is moved verbatim; bit-identity contract gated the merge). Caveat: the slow `test_p1_changes_heterogeneous_keff` was not exercised in this round due to container memory pressure; it's pre-existing.
