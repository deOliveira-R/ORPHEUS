---
status: candidate
project: ORPHEUS
modules: [sn]
primary_criteria: [structure-exposing]
secondary_criteria: [structurally-simpler, expressive]
requires_validation: true
---

# Candidate: Cylindrical α-redistribution as a connection coefficient

## Status

This is a candidate reformulation. It has NOT been validated
in ORPHEUS. High prior for payoff based on standard
differential geometry theory applied to curvilinear
transport, but the work to rewrite cylindrical SN in this
frame has not been done.

Validation path: derive α coefficients from the cylindrical
metric's Christoffel symbols and compare to Bailey / Morel-
Montry weights. Then implement one version in the
connection-based frame, compare against the current code for
(a) the cylindrical test matrix, (b) new per-ordinate
consistency checks. Assess whether signs-from-metric claim
holds algorithmically.

## Problem (current formulation)

In cylindrical SN, the ordinate direction is not constant
along a ray because the coordinate basis rotates with
azimuthal angle. The balance equation acquires an
α-redistribution term that couples ordinates at level m to
ordinates at level m ± 1/2:

streaming + redistribution + collision = source
redistribution ∝ (ΔA / w)(α*{m+1/2} − α*{m-1/2}) ψ_m

The α coefficients are derived by requiring flat-flux
consistency: for ψ = constant, streaming and redistribution
must cancel per-ordinate (Bailey et al. 2009, M-M weights).

Symptoms of the wrong-frame formulation:

- Smell 3 fires: α feels ad-hoc; signs must be memorized
- Smell 11 (variant): the ΔA/w "geometry factor" has a
  mysterious interpretation
- Different sources use different sign conventions
  (absorbing the minus sign into α or not); notation
  conflicts are unresolved by derivation
- A diagnostic check (per-ordinate flat-flux residual) is
  needed to catch implementation bugs — signals that the
  native formulation does not self-enforce the constraint
  structurally

## Structural trigger

From reference.md Part A:

- A.1 Differential geometry: "curvilinear coordinates"

From reference.md Part C:

- Smell 3: magic correction term with hand-waved derivation
- Smell 11 (variant): the geometry factor is effectively a
  stabilization of the flat-flux constraint rather than a
  structural consequence

## Frame

Differential geometry, specifically: the covariant
derivative of a vector field on a curved coordinate system.
The α-redistribution is the connection coefficient
(Christoffel symbol) contribution to the covariant
derivative of the ordinate direction along the streaming
path.

## Reformulation (sketch)

The streaming operator in cylindrical coordinates is the
directional derivative of ψ along the ordinate vector Ω. In
flat (Cartesian) space, Ω is constant and ∂/∂s is a partial
derivative. In cylindrical space, Ω rotates with azimuthal
angle because the basis {ê*r, ê*θ} rotates. The correct
derivative is the covariant derivative:

∇_Ω ψ = Ω^i ∂i ψ + Ω^i Γ^k{ij} Ω^j ∂ψ/∂Ω^k

The second term — the Christoffel symbols Γ — is the
connection contribution. When discretized in the (μ, η, ξ)
basis for cylindrical quadrature, this term becomes the
α-redistribution.

Consequences:

- The signs of α are determined by the metric, not chosen
- The ΔA/w factor is the discrete analog of √g (metric
  determinant) weighting
- The flat-flux consistency is automatic: if ψ is constant,
  ∇_Ω ψ = 0 identically (covariant derivative of a scalar ×
  constant = 0), which implies streaming + redistribution
  cancel per-ordinate structurally
- Different sign conventions in the literature are different
  choices of signature / orientation on the connection; the
  reformulation makes these choices explicit parameters

## Expected elegance payoff

- **Structure-exposing**: α coefficients become a consequence
  of the metric, not a memorized formula. The per-ordinate
  flat-flux consistency becomes structural, not diagnostic.
- **Structurally-simpler**: one formulation (covariant
  derivative) replaces two (streaming + redistribution as
  separate objects). Extends to spherical (two
  redistributions) and general curvilinear without
  rederivation.
- **Expressive**: sign conventions become named choices
  (signature, orientation) rather than conflicting
  literature notations.

## Concrete first test

Derive the α coefficients from the cylindrical metric
Christoffel symbols. Compare to the Bailey M-M weights
entry-by-entry. They must agree.

Pass condition: entry-by-entry agreement (up to documented
sign convention choice). If they do not agree, the
reformulation is wrong (not the classical derivation). If
they agree, proceed to a full rewrite of the cylindrical SN
kernel in the connection frame and measure code-line count
reduction, extensibility to spherical, and whether the
per-ordinate diagnostic test becomes trivially passed
structurally.

Fail condition: disagreement with Bailey weights that cannot
be resolved by sign / orientation convention. In that case,
document the specific disagreement in this file and reject
the frame.

## Literature path

- Cartan, "Leçons sur la théorie des espaces à connexion
  projective" (1937) — original connection theory
- Wald, "General Relativity" (1984), chapter 3 — covariant
  derivative introduction with clear geometric motivation
- Do Carmo, "Riemannian Geometry" (1992) — standard
  reference for Christoffel symbols in arbitrary coordinate
  systems
- Bailey, Morel, Chang, "An analysis of cell-centered
  finite-difference neutron transport schemes on
  orthogonal curvilinear meshes" NSE 162 (2009) — current
  derivation of α
- Morel & Montry, "Analysis and elimination of the
  discrete-ordinates flux dip" TTSP 13 (1984) — flux dip
  analysis
- Modern transport geometry: Adams, Larsen, Morel (various)
  sometimes allude to this connection but do not develop it
  formally

## Transferable pattern

Look for this reformulation when you see:

1. A "correction term" in a curvilinear-coordinate PDE that
   is derived ad-hoc per coordinate system
2. A diagnostic check needed to verify a property that
   should be structural (smell 11 variant)
3. A "geometry factor" whose interpretation is unclear

The connection-coefficient reformulation exposes the
coordinate-independence of the underlying transport and
derives corrections from the metric.

Apply to (candidate targets):

- Spherical SN redistribution (both μ and η redistributions
  should both be Christoffel contributions)
- General-curvilinear SN (mixed coordinate systems for
  complex geometry)
- Transport on manifolds (e.g., accelerator target
  geometries)
