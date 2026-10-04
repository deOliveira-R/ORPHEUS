---
name: traced-memo-client-census
description: #405 P3 traced-memo clients — pickle skips __post_init__ (untraced construction), MR solver result all-data, memo key needs signature binding
metadata:
  type: project
---

Census of 2026-10-04 at `9def3208`; the full report is `scratch/reference_architecture/p3/census.md`. These are claims about the tree on that date, so re-verify them before repeating them.

- **A pickle of a plain frozen dataclass restores `__dict__` and never re-runs `__post_init__`.** A traced child that unpickles `TrajectoryResolventDerivation` therefore never traces `Billiard.__post_init__`, `_route` or `_layered_xs_payload`. Positive control: constructing the derivation in the child traces all of them. `ContentIdentity.__reduce__` routes through the constructor, so a specification is safe.
  **Why:** this is a stale-hit hazard of architecture B that the plan did not list.
  **How to apply:** for any "parent pickles, child traced" design, check whether the boundary objects are rebuilt in the child.
- **The derivation stops pickling after a flux read.** `scalar_flux` is a `cached_property` returning a `functools.cache` closure, so the derivation is unpicklable after any `PointValue` or `FluxIntegral` read.
- **The MR solver results are plain data, and every field is read.** `CylinderGreensMRResult` has 11 fields and the sphere result 10; every field is read by some consumer. A memo at the function serves Billiard and the direct callers bit for bit, but only after `signature.bind` + `apply_defaults`, because Billiard omits `None` kwargs.
- **The reference's rays need more than the plan's solve child.** Rebuilding them needs `r_nodes`, the angle arrays and `region_at_node`, which the plan's solve child (k, density, fission rate) omits. The rebuild itself costs about 0.004 ms.
- **A trace from the child's first line includes about 100 to 120 first-party module bodies** (the import closure), so "a module with a traced function contributes a skeleton" must decide whether `<module>` counts.
- **Census seams Nexus missed:** the `api.SOLVERS` `getattr` dispatch, a bound-function reference (`bare = f if … else g`), and the string-named `s5.c(...)` helpers in `tests/gates/reference`.
