# Issue #132 — Agreement with the reference is not agreement with the physics

**The bug.** The method-of-images approach to white-BC sphere
closure constructed an image-series solution to the Peierls
integral equation on the homogeneous sphere. The series was
well-defined, converged fast (n_max = 5 saturates within 1e-6),
and the iteration ran cleanly under adaptive quadrature. But the
image construction is the mathematical solution to the *specular*
boundary problem (image with angular flip), not the white BC
(Lambertian re-emission with no angular memory). The numerics
solved the wrong problem perfectly.

**Evidence that existed before it was caught.** The vacuum-BC
sanity gate passed: λ_max(M) → 1/σ_t = 1.0 as R → ∞, the correct
asymptotic limit. The image series converged numerically to high
precision. Davison's u = rφ substitution made the spherical
Peierls equation *look* like a 1-D slab, which is a known and
respected reference manipulation. The candidate had analytical
pedigree, numerical convergence, and a passing sanity check —
the conventional triad that licenses "verified."

**Why that evidence didn't catch it.** The analytical target was
itself wrong. White (Mark) BC re-emits with the average angular
distribution, `ψ⁻(r_b, Ω) = J⁺(r_b)/π`, independent of per-ray
Ω. This re-emission CANNOT be reproduced by mirror images of the
source — no method-of-images formulation exists for white BC,
even on a homogeneous sphere. The image series was a perfectly
correct solution to a different physics problem (specular
reflection at the surface), and that different problem happens
to have a clean analytical structure that converged smoothly.
The vacuum-BC sanity gate was structurally unable to distinguish
"white" from "specular" because the BC is irrelevant in the
R → ∞ limit. For sphere 1G/1R (σ_t=1, σ_s=0.5, νσ_f=0.75, R=1),
the image series gave k_eff = 0.704 against the cp_sphere
white-BC reference k_inf = 1.500 — a permanent 53 % error, not a
quadrature artefact.

**What evidence class would have caught it.** Every reference you
trust **MUST** be traced from its boundary-condition mathematics,
not just from numerical convergence and asymptotic limits.
**NEVER** treat "the candidate converged + the candidate matches
a known limit" as verification of the candidate's *physics* —
**instead**, audit the boundary-condition translation step: the
point at which a manipulation ("u = rφ makes the sphere look
like a slab") embeds an implicit BC that may not match the
target. The decisive falsification was a non-asymptotic
comparison against an independently constructed white-BC
reference (cp_sphere) at finite geometry, where the white vs
specular distinction is load-bearing. Sanity gates that wash out
the relevant degree of freedom (here, the BC) are not
verification — they are necessary conditions that admit the
wrong problem. Reference contamination is the structural-failure
case par excellence: numerical convergence and reference
agreement both held; the bug was *which equation was being
solved*. See vv-principles reference.md §6 (reference
contamination as a named failure mode).

**References.** Issue #132 (image-series viability probe
FALSIFIED); numerics-investigator agent memory
`issue_132_davison_image_series.md`; diagnostic scripts
`scratch/derivations/diagnostics/diag_sphere_davison_image_{01..04}_*.py`.
