---
name: ray-characteristic-machinery-map
description: Where "a straight line through a concentric partition, attenuated integral, billiard closure" is spelled (six bodies, L0 and MoC), and where the Branch-1/2 line falls.
metadata:
  type: reference
---

Surveyed 2026-10-05 at `a336bde4` (characteristic-reference plan, Q2); full table `scratch/characteristic_architecture/existing_machinery.md`, probe `_probe_ray_spellings.py` beside it. Re-verify before repeating counts.

- The line-through-shells optical depth has six bodies that agree to ~1e-14: Variant-α segment helpers, `kernels.chord_half_lengths`, Peierls `CurvilinearGeometry.optical_depth_along_ray` and `_chord_tau_mu_sphere`, MoC `_trace_single_ray`. No L1/L2/input home: the manifold/quotient machinery is angular only (orbit catalogue keys all `(Sphere, ·)`).
- `chord_half_lengths` is imported by production CP AND three references; its docstring's L0 test does not exist. That is the standing X4 case.
- Peierls' walker is the only ray primitive with its own L0 tests; Variant-α imports no other ray body (independent today; the Variant-α vs Peierls xverif gate depends on that).
- The lift is spelled three ways (λ pointwise, MoC `/sin_p`, CP/Peierls Ki_n). The billiard closure five ways (scalar per ray, continuous-μ abandoned, mode-space `BoundaryClosureOperator`, CP Sherman–Morrison, MoC track links).
- Ruling on record: `posing_sequence.md` "The free flight is the MoC Volterra operator", designed once in the MoC campaign (apply + sample), so MC inherits it (#534). #436 = face pairings as one datum.

**How to apply:** for any ray, chord, MoC or CP-geometry question, start here. Separate the pure geometry, which can be shared only if it is independently pinned, from the integral and closure, which each branch must own. Related: [[spatial-transform-category-durable]], [[coordinate-system-group-seam]].
