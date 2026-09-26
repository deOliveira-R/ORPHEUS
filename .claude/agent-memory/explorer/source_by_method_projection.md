---
name: source-by-method-projection
description: Volumetric source per method (SN/diffusion/CP/MoC/0-D/MC) — what each can represent, who projects, the CP/MoC normalising "fixed-source" trap, the curvilinear start-direction fold
metadata:
  type: project
---

Measured 2026-09-25 at `99ac3e66` (#405 P1). Full report: `scratch/reference_architecture/p1_source_by_method.md`.

- The source q(r,Ω,g) is method-INDEPENDENT; what varies is the representable subspace, the projection rule (SN collocates at ordinates; CP/MoC/diffusion/0-D need ℓ=0 region averages; MC needs a sampler) and the measure-owned density (GL W=2 per unit μ, sphere rules W=4π, MoC 1/4π).
- Trap: CP's and MoC's `solve_fixed_source` are eigen inners that RENORMALISE (CP `1/max(phi)` on both paths; MoC `1/total_prod` unless production is 0). Only SN and diffusion's LU resolvent answer an absolute FixedSource.
- Curvilinear SN synthesises the starting-direction (μ=−1) source from ordinate values (`source_from_angular`): exact for polynomials of degree < N, Gibbs for a beam (it goes negative).
- MoC MMS never runs production `MOCSolver`: it runs the reference-side `mms_sweep` clone.
- On concentric-annulus FSRs, a Cartesian linear source cannot carry radial variation (a zero first moment); non-flat MoC there means a radial basis.
- `wide_moc` runtime capture binds 0 exercisers (unattributable); do not use it.

**Why:** the reference-spec `Source` design (Discussion 4, reopened by the user 2026-09-25) hinges on these.
**How to apply:** re-verify the line anchors before citing them; see [[question-and-source-vocabulary]].
