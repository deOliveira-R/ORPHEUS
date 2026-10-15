---
name: c2-cylinder-point-row-census
description: #405 P3 C2 measured — 60 cylinder rows cost 239 s/pass (median 2.65 s, not 10-20 s); 40 of 60 under test bypass; row cost is AttenuationRule.of in the exit integral
metadata:
  type: project
---

Tree da74b1ce (2026-10-15); census `scratch/characteristic_architecture/p3/c2/census.md`. Re-verify before repeating.

- 20 of 1739 tests make cylinder rows (60 rows); 40 rows sit under a test `bypass()` (albedos, lift, reference), so
  neither C1 nor a row memo reaches them: bypass (C3) is the lever (630 s/pass incl. 453 s in-process blocks).
- A row's cost is `AttenuationRule.of` inside the exit integral (86 %); `PointRule.of` is ms.
- `point_flux` computes every group's row; the door keeps one (11 of 34 door rows wasted).
- Instrument lesson: a residual (test time minus block time) as a per-stage estimate absorbed in-process blocks and
  mpmath — it overstated rows 4-8x. Time the stage itself with a recording wrapper.
- A dry census's test-level bypass flag needs a COUNTED bypass wrapper installed before test import, separate from the
  session's own bypass.
