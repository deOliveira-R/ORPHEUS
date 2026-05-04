---
name: R1 power_iterate_variant_alpha driver
description: Hindsight refactor R1 — collapsed 11 power-iteration loops in trajectory_resolvent into a single driver. Bit-equality preserved; 12-config baseline + full 205+84 test suite identical pre/post.
type: project
---

# R1 — `power_iterate_variant_alpha` driver closeout

`refactor/r1-power-iterate-driver` (off `main`), 2026-05-03.

## Scope

Collapse the structurally-identical Variant α power-iteration outer
loops in `orpheus/derivations/continuous/trajectory_resolvent/`
(sphere 1G + MG + MR; cylinder 1G + MG; slab_asymmetric 1G + MG;
hollow_sphere 1G + MG; annulus 1G + MG) into a single driver
`power_iterate_variant_alpha`. R1 is the lowest-risk highest-leverage
refactor in `.claude/plans/trajectory_resolvent_hindsight_refactor.md`.

## Inventory clarification — 11 power-iteration loops, NOT 12

The plan and brief said "12 power-iteration loops". The 12th loop is in
`solve_greens_function_sphere_mr_fixed_source` (greens_function.py:1233)
— it is **NOT** a power iteration. It has no Rayleigh-quotient update
of `k_eff` and converges on `phi_g` (scalar flux) instead. Refactoring
it into the power-iteration driver would conflate two distinct
algorithms. **Decision (load-bearing variation per
`algebra-of-record`)**: leave it untouched. R1 collapses 11 loops; a
future fixed-source-iteration driver can collapse it if more
fixed-source solvers ship.

## Driver design

`power_iteration.py` exposes:

- `PowerIterationResult` — frozen dataclass: `psi`, `k_eff`,
  `iterations`, `converged`, `k_history`, `residual_history`. The
  history fields are NEW — they emerge for free from the unification.
- `power_iterate_variant_alpha(step, initial_psi, initial_k, *,
  max_iter, tol)` — accepts a per-step callable that does the
  geometry-specific arithmetic in EXACTLY the order the original loop
  did, and returns `(psi_new, Frate_old, Frate_new)`. The driver does
  the universal book-keeping (Rayleigh quotient, normalisation,
  rel-k convergence test, history collection).

The `step` callable is the load-bearing design choice. It encapsulates
the geometry-specific work without forcing the driver to know about
phi/source-profile/geometry-Jacobian shapes. Each call site captures
its operator and grid in a closure.

## Bit-equality preservation

Each call site's `_step` closure runs the SAME FP operations in the
SAME order as the pre-refactor inlined loop. The driver's universal
arithmetic (`k_new = k_eff * Frate_new / Frate_old`,
`psi_normed = psi_new / Frate_new`, `rel_dk = abs(k_new - k_eff) /
max(abs(k_eff), 1e-30)`) is identical to the inlined originals.
IEEE-754 reproducibility guarantees bit-equality of every numerical
output.

**Verification**: 12 baseline configurations captured at HEAD via
exact-bit hex-float repr (`float.hex(...)`). Re-running each post-
refactor reproduces every `k_eff`, `np.sum(psi)`, `np.linalg.norm(psi)`,
and `iterations` to **bit-exact equality** (every digit identical).
Configurations covered: sphere (1G closed/vacuum + MG); cylinder
(1G + MG); slab asymmetric (1G symmetric + 1G α_L=0.5/α_R=0.8 +
MG); hollow sphere (1G + MG); annulus (1G + MG).

## Sphere 1G subtlety

Sphere 1G uses a different `_apply_operator` API (passes `psi` directly
+ a `source_coeff` scalar; computes phi internally) than the other 10
loops (which compute phi explicitly + pass a `source_profile` ndarray).
The driver tolerates this without special-casing because each `_step`
closure encapsulates whichever calling convention its geometry uses.
The phi computation happens twice for sphere 1G (once inside
`_apply_operator`, once outside for the Frate ratio) — bit-equal because
the function is deterministic.

Sphere 1G also did not have a Frate-vanish check pre-refactor; the
driver adds one universally (the message standardises to
"non-multiplying medium" — minor incidental drift on error message
text across the original 11 loops, see "Variation analysis" below).

## Variation analysis (pre-refactor inventory)

| Variation                                | Class              | Disposition                          |
| ---------------------------------------- | ------------------ | ------------------------------------ |
| Frate-vanish error message text          | Incidental drift   | Standardised to one message          |
| Sphere 1G missing Frate-vanish check     | Incidental drift   | Driver adds universal check          |
| Sphere 1G `_apply_operator` API          | Load-bearing       | Closure encapsulates                 |
| MG per-group `for g in range(G)` inner loop | Load-bearing    | Lives inside the `_step` closure     |
| MR per-node σ_s lookup                   | Load-bearing       | Lives inside the `_step` closure     |
| Volume-integral Jacobian (4πr², 2πr, 1)  | Load-bearing       | Captured in geometry's `fission_rate` closure |
| Sphere MR fixed-source iteration         | **Different algorithm** | Out of scope; left untouched     |

## Foundation tests

`tests/derivations/test_trajectory_resolvent_power_iterate.py` — 6
foundation-tagged tests that pin the driver's universal book-keeping
in isolation from any geometry:

1. `test_driver_converges_on_constant_step` — loop runs `max_iter` when
   the residual stays constant.
2. `test_driver_detects_convergence_on_unit_ratio` — `Frate_old ==
   Frate_new` triggers immediate convergence.
3. `test_driver_raises_on_vanished_fission_rate` — `Frate_old < 1e-30`
   raises with the iteration number in the error message.
4. `test_driver_returns_history_lengths_match_iterations` — both
   `k_history` and `residual_history` have exactly `iterations` entries.
5. `test_driver_preserves_bit_equality_on_synthetic_iteration` — the
   load-bearing probe: bit-for-bit equality between the driver's
   output and an inlined reference loop on a synthetic
   geometry-free iteration.
6. `test_driver_returns_powerresult_dataclass` — return type contract.

## Test results

Pre-refactor (HEAD `main`):
- `tests/derivations/test_peierls_greens_function_*`: 205 passed
- `tests/cross_method/`: 84 passed

Post-refactor (`refactor/r1-power-iterate-driver`):
- `tests/derivations/test_peierls_greens_function_*`: **205 passed**
- `tests/cross_method/`: **84 passed**
- `tests/derivations/test_trajectory_resolvent_power_iterate.py`:
  **6 passed** (new — driver foundation tests)

Sphinx `-W -b html`: clean build (no warnings).

## Files changed

- NEW `orpheus/derivations/continuous/trajectory_resolvent/power_iteration.py`
- MOD `orpheus/derivations/continuous/trajectory_resolvent/greens_function.py`
  (sphere 1G, MG, MR loops collapsed to driver calls)
- MOD `orpheus/derivations/continuous/trajectory_resolvent/greens_function_cylinder.py`
  (cylinder 1G + MG loops)
- MOD `orpheus/derivations/continuous/trajectory_resolvent/greens_function_slab_asymmetric.py`
  (slab asym 1G + MG loops)
- MOD `orpheus/derivations/continuous/trajectory_resolvent/greens_function_hollow_sphere.py`
  (hollow sphere 1G + MG loops)
- MOD `orpheus/derivations/continuous/trajectory_resolvent/greens_function_annulus.py`
  (annulus 1G + MG loops)
- NEW `tests/derivations/test_trajectory_resolvent_power_iterate.py`

Net: 11 power-iteration loops collapsed to 1 driver + 11 small
`_step` closures. Geometry files lost ~32 LOC each (init + loop
preamble + loop tail) for an aggregate ~165-line reduction across 5
files; the driver adds ~210 LOC (mostly docstrings).

## Coordination with parallel agent

R0.5 (`MeshTemplate` promotion + cross-method adapter refactor) was
run in parallel on `refactor/r05-mesh-template-promotion`. Zero file
overlap with R1. The R0.5 agent edited `MEMORY.md` mid-session; that
edit is excluded from R1's commit (R0.5's branch will carry it).
Cross-method regression net (84/84) passed identically pre and post
R1 — confirming both refactors are bit-equal independent of each
other.

## Out of scope (defer to R2/R3)

- R2 (`GreensResult` + `PhaseSpaceMesh`) — return-shape unification.
  R1 keeps the existing 12 result dataclasses; the driver's
  `PowerIterationResult` is unwrapped into the existing return shapes
  at each call site.
- R3 (`ChordOracle` Protocol) — the mathematical-leverage refactor.
  Lands after R1 + R2 bed in.
- Fixed-source loop unification — needs ≥2 fixed-source solvers
  (currently only 1).
