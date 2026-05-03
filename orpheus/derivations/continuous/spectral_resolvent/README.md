# `spectral_resolvent/` — Meaning-2 reference solver (reserved)

**Status:** empty / placeholder. Implementation deferred.

## Intent

Implements the **closed-form scalar Green's kernel** for slab + sphere via
spectral μ-integration. The construction integrates Sanchez 1986 Eq. A1
(equivalently Pomraning–Siewert 1982 Eq. 21) directly in spectral form,
rather than reconstructing it numerically by ray tracing + multi-bounce
closure.

This folder is the **structurally-independent sister** to
[`trajectory_resolvent/`](../trajectory_resolvent/), which builds the
same scalar Green's kernel by the trajectory route (characteristics +
resolvent T = (I − S)⁻¹). Both folders produce the same physical Green's
function under specialisation to the homogeneous medium, and the
two-path agreement is treated as an L1 cross-check anchor.

The taxonomic split — Meaning 1 (angular Green's via Case spectrum,
[`singular_eigenfunction/`](../singular_eigenfunction/)) vs Meaning 2
(scalar Green's kernel via spectrum vs trajectories) — follows the
Sanchez–Chandrasekhar three-meanings analysis documented in
`.claude/agent-memory/literature-researcher/sanchez_chandrasekhar_three_meanings.md`.

## Canonical references

* **Sanchez, R. (1986).** "Integral form of the equation of transfer for
  a homogeneous sphere with linearly anisotropic scattering."
  *Transport Theory and Statistical Physics* **15**(3), 333–343.
  DOI: 10.1080/00411458608210456. **Eq. A1** is the spectral closed-
  form scalar kernel this folder targets.
  Local PDF: `scratch/literature/Sanchez1986.pdf`.

* **Pomraning, G.C., Siewert, C.E. (1982).** "On the integral form of
  the equation of transfer for a homogeneous sphere." *J. Quant.
  Spectrosc. Radiat. Transfer* **28**(6), 503–506. **Eq. 21** is the
  same kernel for the sphere with diffuse + specular reflection.
  Local PDF: `scratch/literature/PomraningSiewert1982.pdf`.
