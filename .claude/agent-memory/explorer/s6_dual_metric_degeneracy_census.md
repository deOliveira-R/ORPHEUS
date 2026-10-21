---
name: s6-dual-metric-degeneracy-census
description: S6 (#617/#625/qa F2) census at 3da7fed3 — DualSpace metric reach, composite-dual drop, V*==V on axis-built spaces, degenerate-metric families, override warrants
metadata:
  type: project
---

Tree s6 worktree @ 3da7fed3 (2026-10-21). Report: `scratch/numerics_untangling/packages/S6_census/report.md`, probes beside it.

Durable shape (re-verify before repeating):
- The metric family already holds both faces (DenseMetric stores pinv(G) AND G; `inverse_of` swaps them), so a G⁻¹ dual metric is a face swap, not a new factorisation.
- `DualSpace.of` copies `inner_product_weights`/`axes`/`metric` only: the dual of a `FullFieldSpace`/`CoupledSpace` is EUCLIDEAN (their metric lives in member fields). #617's "dual carries G" is true for leaf spaces only.
- Axis-built space identity = axes, so `V.dual() == V` (same hash) and `R @ R.dual()` composes. Covariance is not in identity.
- `_DualOperator.domain/codomain` mint a fresh DualSpace on every read (~27 600 builds in a numerics+sn/operators run).
- 0 production reads of a DualSpace metric verb; only `.H` reaches the inverse metric in production (legs built only in `AdjointOperator.__init__`).
- Degenerate metrics: trace |Ω·n|w zeros, moment-head layout zeros, rank-deficient dense Parseval heads; handled by silent Moore–Penrose. `A.H.H is A`, the retraction→pullback and pointwise→self overrides are theorems only on non-degenerate ends (measured 0.348 / 0.5 off).
- `FunctionSpace.riesz_lower/raise` were retired in S1; #617's done-when still names them.
- Production PG frames DO reach `basis_space` (`mixture.py` `project`), with 0 metric reads.

Spy recipe that worked: a pytest plugin patching `__post_init__` across the `FunctionSpace.__subclasses__` tree at `pytest_collection_finish`, with an outer-wrapper check `type(self).__post_init__ is wrapper`. xdist breaks on this tree (collection differs per worker), so run serial in 25-file chunks. The cwd-dependent failure `test_axis_generator.py::...G5_6b` shows up when run from /tmp; it is not the spy.

Related: [[adjointability-flag-semantics]], [[moment-metric-slots-and-cross-axis-forms]], [[operator-composite-algebra-seam]].
