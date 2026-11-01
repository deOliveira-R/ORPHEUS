---
name: space-membership-and-covector-census
description: Where element-in-space membership is (and is not) checked, the payload-tail convention, the co-vector family and gauges — measured 2026-10-30 at e42fd498
metadata:
  type: project
---

Durable shape, measured 2026-10-30 at `e42fd498` (re-verify before repeating). Detail:
`scratch/layer_dag/ipf_fix/explorer_map.md`.

- Membership is stated on the ELEMENT (`Field.__post_init__` exact shape) and at 10 of 67 operator
  `apply`s (SN `require_member` pair + AxisMarginal), never on the space or metric: `DiagonalMetric.validate_for`
  returns None, so no verb below `_resolved_metric` knows the element shape.
- Space verbs legitimately receive PAYLOAD-TAILED elements (head-only space + trailing (ng, spatial)):
  7 gate functions in 4 files (frame-face `.H`, moment-head metric) rely on it; production 0 in the
  sampled runs. Exact-shape refusal at the verb breaks them (field_algebra half-2 rule).
- IntegratedReactionRate == space.inner_product on every captured space (20 k calls, curvilinear and 2-D/3-D
  included); EnergyAxis refuses weights, but a DenseMetric overlay on the energy block is ADMITTED
  (probe: 83 vs 46) — the one constructible divergence.
- No bound scalar co-vector and no element on a DualSpace exist; dual elements are role leaves
  (SourceSink) on the primal space. Four spellings of the integration co-vector <1,.>_V.
- KernelGauge: one production conformer (LossKernelGauge); projector machinery is generic numerics
  (frame_faces.conjugate + CrossGramInverse), basis is SN-only; diffusion's LU is SILENT on an exactly
  singular loss (1.7e17), homogeneous refuses.

See [[moment-metric-slots-and-cross-axis-forms]], [[adjointability-flag-semantics]].
