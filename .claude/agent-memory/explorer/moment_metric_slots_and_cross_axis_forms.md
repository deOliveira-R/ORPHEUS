---
name: moment-metric-slots-and-cross-axis-forms
description: Where a space metric can live (2026-10-18, main 8e109132) — axis-built spaces refuse any form spanning two axes; the LD moment mass is unit-width diagonal; R pairing is the only per-cell moment-indexed array
metadata:
  type: project
---

Durable shape, measured 2026-10-18 at `8e109132` (re-verify before repeating).

- An axis-built `FunctionSpace` admits a metric object only as a `FactoredMetric` with ONE entry per axis (`_require_positioned_over`, space.py). A block spanning `(spatial, spatial_moment)` is REFUSED there, yet the same object is admitted on an axes-less space of equal shape. So the docs' "machinery available" for the non-separable curvilinear LD mass holds only off the axis-built path. Probe: `scratch/numerics_untangling/frame_bridge/facts/probe_metric_slots.py`.
- `SpatialScheme` does not exist; the owner is `DiscretizationSchemeBase`.
- LD's mass has two spellings: space side `moment_mass_diagonal` (unit-width `[1, θ]`, cell-independent) and kernel side `mass_1d(h, θ)` per-cell `(…,2,2)`; residual kernels divide by its diagonal.
- `ReducedStreamingOperator.redistribution_pairing` `(nx, n_mom, n_thread)` is the only per-cell moment-indexed array; refused at n_mom>1 by the closure. `scheme.mint(chart)` was never built (0 `def mint`); no live work item names it after the 2026-08-30 close-outs (#612, #158, #423 open).
- Census instrument: `space_census_plugin.py` (wraps every FunctionSpace subclass `__post_init__`); xdist collection is non-deterministic here, run per tree with `-p no:xdist`.

See [[harmonic-frame-folded-quadrature-plumbing]], [[affine-operator-split-convention]].
