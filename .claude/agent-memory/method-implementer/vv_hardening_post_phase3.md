---
name: V&V hardening pass — post Phase-3 family
description: 2026-05-02 V&V gap closures on the Variant α Green's-function family — Tasks 2-4 (slab-sym structurally-independent cross-check; off-diagonal intermediate α convergence for the 3 rank-2 geometries; grazing-ray stability across all 6 geometries).
type: project
---

# V&V hardening pass — post Phase-3 family closeout

**Branch**: `feature/peierls-greens-cylinder` (started at HEAD `d63f654`).
**Date**: 2026-05-02.
**Scope**: V&V gap closure post the QA hindsight review (Tasks 2-4;
Task 1 cylinder external-reference cross-check is delegated to a
parallel literature-researcher dispatch and not part of this work).
**Outcome**: 12 new tests added; pre-existing 178 tests preserved at
their pre-change tolerances; full peierls_greens suite (72 tests across
the 6 geometry files including all `slow` variants) passes in
~7m48s.

## Files modified (test files only — production code untouched)

- `tests/derivations/test_peierls_greens_function_solver.py`
  — added `test_grazing_ray_stability_sphere`.
- `tests/derivations/test_peierls_greens_function_cylinder_solver.py`
  — added `test_grazing_ray_stability_cylinder`.
- `tests/derivations/test_peierls_greens_function_slab_solver.py`
  — added Nyström import; added
    `test_alpha_zero_vacuum_agrees_with_nystrom_slab` (Task 2);
    added `test_grazing_ray_stability_slab`.
- `tests/derivations/test_peierls_greens_function_slab_asymmetric_solver.py`
  — added `test_off_diagonal_intermediate_alpha_slab_asym`
    (parametrized ×2);
    added `test_grazing_ray_stability_slab_asym`.
- `tests/derivations/test_peierls_greens_function_hollow_sphere_solver.py`
  — added `test_off_diagonal_intermediate_alpha_hollow_sphere`
    (parametrized ×2);
    added `test_grazing_ray_stability_hollow_sphere`.
- `tests/derivations/test_peierls_greens_function_annulus_solver.py`
  — added `test_off_diagonal_intermediate_alpha_annulus`;
    added `test_grazing_ray_stability_annulus`.

## Per-task verdicts

### Task 2 — Slab-sym secondary structurally-independent cross-check

**Verdict**: candidate (a) — slab-sym Variant α at α=0 vs `peierls_nystrom`
slab solver at vacuum BC. **Landed.**

The slab Variant α and Nyström slab solvers descend from entirely
different mathematical reductions of the slab Peierls equation:

- **Variant α**: angle-resolved Green's function on `(x, μ)` phase
  space with bouncing characteristics, GL on `(x, μ)`, trajectory and
  bounce-period quadrature, then `2π ∫ψ dμ` projection.
- **Nyström**: scalar-flux integral equation
  `φ(x) = ½ ∫ E_1(τ(x,x')) q(x') dx'` with adaptive `mpmath.quad`
  handling the `E_1` log singularity.

Both share `mpmath.expint` / `scipy.special.exp1` BELOW the
trusted-library line — fine per `algebra-of-record` § "Structural
independence applies above the trusted-library line". The
operator-construction code paths above the line are independent.

**Achieved tolerance** at fuel-A τ_L = 5: 1.4e-5 relative agreement
between Variant α @ (24, 40, 96) and Nyström @ (n_panels=8, p=6,
dps=20). Gate set at 5e-5 (~3.5x margin).

**Caveat (logged as future work)**: the Nyström slab supports
`boundary ∈ {"white", "vacuum"}` only — NOT `α ∈ (0, 1)`. The
cross-check therefore covers α=0 only. Intermediate α coverage
requires an α-aware Nyström-style reference (or an SN-slab solver
upgrade). The pre-existing slab-sym ↔ slab-asym delegation gate
(now documented as same-code-twice) covers intermediate α at
self-consistency level.

The pre-existing
`test_rank1_path_now_agrees_with_rank2_via_delegation_after_ERR035_fix`
remains in place as a regression-prevention guard against future
Phase-3A heuristic re-introduction; its same-code-twice character
is now annotated explicitly through the new structurally-independent
sibling.

### Task 3 — Off-diagonal asymmetric intermediate α convergence

**Verdict**: 5 new convergence-to-self tests landed across the 3
rank-2 geometries.

| Geometry      | Off-diagonal pairs              | Achieved finest-pair self-consistency | Gate          |
| ------------- | ------------------------------- | ------------------------------------- | ------------- |
| slab-asym     | (0.7, 0.4) and (0.3, 0.85)      | 9.8e-6 (research-grade)               | 2e-5 (~2x margin) |
| hollow sphere | (0.4, 0.7) and (0.85, 0.3)      | 3.5e-4                                | 1e-3 (matches existing convergence-floor gate) |
| annulus       | (0.4, 0.7)                      | 1.9e-3 (oscillatory about a mean)     | 5e-3          |

**Convergence-floor band project record** (per geometry, fuel-A
material):

- **slab-asym** (1D phase space `(x, μ)`): 1e-5 achievable on
  off-diagonal at modest cost, ladder up to (40, 48, 128).
- **hollow sphere** (2D phase space `(r, μ)` + impact-parameter
  partition surface `b = R_in`): geometric-class quadrature error
  at the partition surface dominates; ~3e-4 finest-pair achievable
  on (48, 48, 96) → (64, 64, 128). 4D phase space with the partition
  effectively makes this 3D in error scaling.
- **annulus** (3D phase space `(r, μ_axial, φ_az)` + impact-
  parameter partition surface + cylinder-specific 3D chord lift
  `τ / sqrt(1 - μ_axial²)`): k_eff appears to oscillate slightly
  about a mean across grid orders rather than converging
  monotonically. 1-2e-3 self-consistency floor at modest cost; tighter
  would require very fine grids.

The convergence noise on hollow sphere and annulus is well-understood
from first principles: the impact-parameter partition surface
`b = R_in` is a phase-space discontinuity that GL approximates at
geometric-class rate (not spectral). All three geometries display
**physically-plausible** k_eff values within `(0, k_inf)` consistently
across the off-diagonal; the gates are calibrated to allow expected
quadrature noise while catching any genuine regression.

### Task 4 — Explicit grazing-ray stability tests

**Verdict**: 6 new tests landed (one per geometry).

| Geometry      | Grazing-ray locus              | Refined dim        | Achieved | Gate |
| ------------- | ------------------------------ | ------------------ | -------- | ---- |
| sphere        | `|μ_surf| → 0`                | n_mu               | < 1e-3   | 1e-3 |
| cylinder      | `μ_axial → ±1` AND `\|sin φ_az\| → 0` | n_phi_az    | < 1e-3   | 1e-3 |
| slab          | `\|μ\| → 0`                   | n_mu               | < 1e-3   | 1e-3 |
| slab-asym     | `\|μ\| → 0` (at intermediate α) | n_mu             | < 1e-3   | 1e-3 |
| hollow sphere | `b → R_in` (partition disc.)   | n_mu               | ~1.5e-3 (32→64), ~3.7e-4 (64→128) | 2e-3 + monotonicity |
| annulus       | `b → R_in` AND `μ_axial → ±1`  | n_phi_az           | ~4.9e-3 (16→24), ~1.1e-3 (24→32) | 6e-3 + monotonicity |

**Stability findings**:
- **All 6 geometries are stable** (no NaN, no oscillation, no blow-up
  at grazing-ray loci under quadrature refinement).
- **Slab and sphere convergence rate** is well-behaved in the grazing
  region (≤1e-3 between consecutive refinements at modest grids).
- **Hollow sphere and annulus** require the looser 2-6e-3 gate
  because of the impact-parameter partition surface (b = R_in)
  that is genuinely discontinuous in phase space; convergence is
  geometric-class rather than spectral. Both monotonically converge,
  triggering the monotonicity check (sign-change with shrinking
  magnitude is also accepted to allow overshoot sequences).
- **No surprising sensitivity** flagged. The annulus 16→24 step is
  genuinely noisy (4.9e-3) but settles to 1.1e-3 by 24→32 — typical
  for 3D phase-space at modest grid orders.

## Acceptance gates met

- [x] All pre-existing tests still pass at the SAME tolerances.
- [x] All 12 new tests pass at their (documented) tolerances.
- [x] The Task-2 cross-check uses a structurally-independent code path
      (Variant α angle-resolved vs Nyström scalar-flux Peierls);
      shared upstream is below the trusted-library line.
- [x] Sphinx `-W` build clean (exit 0).

## Test count progression

- Pre-hardening: 178 tests across the Variant α family (60 of those
  in the 6 geometry test files, the rest in symbolic, MR, MG, vacuum,
  xverif files).
- Post-hardening: 178 + 12 = 190 tests; the 6 geometry files now
  collect 72 tests (60 + 12).

## ERR-NNN catalog updates

None — no production-code bug was caught during this V&V hardening
pass. The achieved tolerances on Tasks 3-4 reveal *structural
properties* of the geometries' phase spaces (impact-parameter
partition discontinuity, 3D phase space dimensionality) rather than
implementation defects. These are documented in the test docstrings
themselves rather than as ERR entries.

## Open follow-ups (logged for future work)

1. **α-aware Nyström slab reference**: extend
   `peierls_nystrom.slab.solve_peierls_eigenvalue` to support
   `boundary` parameter accepting `α ∈ [0, 1]` (interpolating between
   "vacuum" and a generalized α-rescaled E₂ closure). Would unlock
   intermediate-α structural cross-check for slab-sym and slab-asym.
2. **External reference for hollow sphere and annulus off-diagonal**:
   nothing in the current literature corpus provides a closed-form
   k_eff for hollow geometries at off-diagonal intermediate α;
   would be a literature-researcher dispatch task.
3. **Sphere-vs-PS-1982 cross-check**: the existing
   `test_peierls_greens_function_xverif_ps1982.py` covers the
   solid-sphere vacuum case and passes. No gap.

## Commits

A single atomic commit covering all 12 new tests is appropriate
(all are V&V scaffolding, no production-code modifications):

- `test(peierls): V&V hardening — structurally-independent cross-check, off-diagonal intermediate α, grazing-ray stability`
