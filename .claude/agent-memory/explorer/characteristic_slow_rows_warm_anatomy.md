---
name: characteristic-slow-rows-warm-anatomy
description: Where the characteristic slow tier's warm 126 s goes after concept 1 (2026-10-17), and door-versus-client facts for point readings
metadata:
  type: project
---

Measured 2026-10-17 on main `0da80c8c`, using the warm worktree `scratch/characteristic_architecture/p3/c3/oracle/wt_after`.
This records the tree at that date; re-verify before repeating it.

- **The split of the 126 s.**

  | part | wall |
  |---|---|
  | SUT point readings (`point_flux` → `PointRule.row`) | about 83 s |
  | C6 reading operator (648 `PointRule` rows) | 21 s |
  | P1's in-process mpmath `psi_disk` quads | about 13 s |
  | unstored `specular_cylinder_total` | about 6 s |

  Memo hits and system construction cost under 0.1 s per row.
- **About half a cylinder point reading is a batched dense `eigh`** in `attenuation_rule._stieltjes`.
- **`_given` is the door's quantity.** `evaluate(PointValue)` of a pure-absorber `FixedSource` matches it to 1 ulp.
- **`GalerkinSystem` is not a key.** `encode_exact` raises `ContentlessError` because it has eq=False.
- **The main store was not warm for `phi_cylinder`.**
- **A scratch probe run from inside a worktree needs `PYTHONPATH=$PWD`.** Without it, `tests.gates` fails to import.
- **The reading batteries (b5b, p3/product) rebind in process with no `bypass()`.**

**Why:** concept 2 of #405 P3 chooses between the door, a new client, and the oracle extensions.
**How to apply:** start from `scratch/characteristic_architecture/p3/c3/concept2/scoping.md`. See [[traced-memo-miss-cost-anatomy]] and [[traced-memo-client-census]].
