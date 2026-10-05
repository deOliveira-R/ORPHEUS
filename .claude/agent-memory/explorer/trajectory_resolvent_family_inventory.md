---
name: trajectory-resolvent-family-inventory
description: Variant-α package shape at a336bde4 — zero production importers, twin map, test claim mix, stale runtime capture
metadata:
  type: project
---

Survey of `orpheus/derivations/continuous/trajectory_resolvent/` at `a336bde4` (2026-10-05); full table in
`scratch/characteristic_architecture/inventory.md`. Re-verify before repeating.

- No module in `orpheus/` outside the package imports it (AST, 1156 files): a test-time reference only.
- Twins: region index x3, chord-at-b x3, sphere first leg x2 (the `greens_function._*` copies dead);
  cylinder first leg = sphere's under mu=cos(phi) to round-off, NOT bitwise (503/2000 differ).
  `SphereChordOracle` = one-region MR sphere to 2.3e-16. `compute_resolvent_T_rank2` has 0 callers.
- Of 298 family test functions: 87 SymPy identities on `origins/`, 17 facade tautologies, 6 driver
  equivalences, 27 k=k_inf (geometry-blind); ~11 independent value pins.
- The `wide_*` coverage captures bind 0 tests to any `orpheus.derivations` node (predate `tests/gates/`).

**Why:** the characteristic-reference re-architecture ([[characteristic-reference-architecture]]) starts from it.
**How to apply:** use as the starting map; re-run the scripts' predicates if the branch moved.
