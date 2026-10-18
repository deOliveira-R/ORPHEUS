---
name: invariance-kernel-measure-kinds
description: Which DiscreteMeasure kinds the invariance verbs are asked of (production all sphere-family); kernel pads non-quotient supports into R^3 and answers.
metadata:
  type: project
---

Measured at HEAD c9529bb8 (2026-10-18), F6 of the numerics untangling.

- Production asks invariance only of S² and Quotient(Sphere) measures: 20,168 of 20,168 production-entered calls over the numerics/geometry/transport/sn gates. Three entries: `Quadrature.ordinate_permutation`, `Quadrature.quotient`, and `AngularSymmetry.admits_symmetry`.
- The kernel never refuses a non-quotient support. It zero-pads into R³ and answers, so an EnergyGroups measure reads O2_x-invariant. Only Sphere has non-trivial catalogue entries.
- The name `AngularMeasure` is already a Protocol in `sn/angular/redistribution.py`. The prior ruling of 2026-08-29 (space_and_kernel_binding plan) was "no measure subtypes; Quadrature IS the angular measure".
- A second invariance door exists: `rules_product._derived_product_group` calls `RigidMotion.preserves` on the circle and μ factors.
- Probe method: a `-p` plugin wrapping the module attributes catches both the measure.py delegation and the kernel's internal calls. The sn gates take about 40 min under the spy.

**Why:** a later session revisiting an AngularMeasure split would otherwise re-run the 1-hour census.
**How to apply:** re-verify with `scratch/numerics_untangling/f6/inv_spy.py` before quoting the counts.
