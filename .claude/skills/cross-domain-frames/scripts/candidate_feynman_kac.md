---
status: candidate
project: ORPHEUS
modules: [cp, moc, mc]
primary_criteria: [structure-exposing, expressive]
secondary_criteria: [algorithmic-advantage]
requires_validation: true
---

# Candidate: CP, MOC, MC as discretizations of one Feynman-Kac path integral

## Status

Candidate reformulation. The Feynman-Kac representation of
transport is well-established theory (Spanier & Gelbard,
Lux & Koblinger), but the ORPHEUS solvers are built as
independent implementations rather than as discretizations
of a common object. The reformulation would expose this
commonality and enable structural cross-method borrowing.

Validation path: pick one cross-method borrowing (candidate:
CP-estimated importance as variance reduction for MC
tallies), implement it via the unified path-measure
formulation, and compare variance per unit compute against
both (a) uniform MC and (b) the equivalent naive
implementation without the unifying frame.

## Problem (current formulation)

ORPHEUS has three transport solvers (CP, MOC, MC) built as
independent codebases with distinct internal
representations:

- CP: matrix of first-collision probabilities, Neumann-series
  inner iteration
- MOC: ray tracks, characteristic-form integration along
  rays, angular quadrature of track contributions
- MC: particle histories, stochastic sampling of collision
  sites and outgoing angles

Cross-method borrowings (QMC for MOC rays, CP-estimated
importance for MC variance reduction, etc.) are done ad-hoc
where done at all.

Symptoms of the wrong-frame formulation:

- Smell 1 (variant): three methods, three different "what
  is a contribution" definitions
- Smell 8 (variant): nested loops with algebraic structure
  (angle × space × energy × history) differ across methods
  with no unification of the underlying contribution
- Cross-method code-reuse is minimal
- Variance reduction ideas from MC don't naturally transfer
  to CP / MOC
- Convergence theory for each method developed independently

## Structural trigger

From reference.md Part A:

- A.4 Feynman-Kac / stochastic representation — MASTER BRIDGE
  for any transport problem
- A.4 Measure theory / Radon-Nikodym — for the specific
  borrowing (importance sampling as measure change applies
  across all three methods)

## Frame

Feynman-Kac representation: the transport equation's solution
is an expected value of a functional of a Markov jump process
(the particle history). Each method is a discretization /
quadrature choice on this path space.

- MC: stochastic sampling of the path measure P(path)
- MOC: deterministic quadrature on paths of a specific type
  (straight rays between collisions)
- CP: deterministic quadrature on the first-collision
  kernel — first moment of the path measure
- Woodcock delta tracking: rejection sampling on the path
  measure (change of measure with majorant Σ)

## Reformulation (sketch)

Core object: path measure P on the space of piecewise-linear
trajectories (free flights + collision events). The solution
ψ(x, Ω) is:

ψ(x, Ω) = E_P [ ∫ q(γ(s)) ds | γ(0) = (x, Ω) ]

where γ is a sampled path, q is the source along the path,
and the expectation is over P. All quantities of interest
(flux, reaction rates, etc.) are functionals of γ.

Each method is a choice of how to evaluate the expectation:

- **Direct MC sampling**: draw γ from P, average over N
  samples
- **MOC**: choose a deterministic subset of paths (rays),
  integrate q along each, quadrature in angle
- **CP**: pre-integrate q over all first-collision paths
  from each region, tabulate, solve a linear system for the
  source self-consistency
- **Delta tracking**: sample γ from a majorized measure P'
  with Σ_majorant ≥ Σ, accept with probability
  Σ / Σ_majorant (Radon-Nikodym derivative)

## Expected elegance payoff

- **Structure-exposing**: the three methods become three
  quadratures on one object, not three algorithms. The
  relationships between them become precise (not analogies).
- **Expressive**: variance reduction (importance sampling
  in MC) is a measure change — which translates to a
  preconditioning choice in CP and a weight rule in MOC.
  "Borrowing" becomes a derivation.
- **Algorithmic-advantage**:
  - QMC for MOC rays becomes natural (it was always
    natural; the frame makes it obvious)
  - CP-adjoint-weighted MC tallies (CADIS-like) generalize
    across methods
  - New methods (hybrid deterministic / stochastic) drop
    out as intermediate quadrature choices

## Concrete first test

Pick one cross-method borrowing: CP-estimated importance
applied as variance reduction for MC tallies.

Step 1: Derive the borrowing from the unified path-measure
formulation. The derivation should produce explicit weight
factors from a Radon-Nikodym argument.

Step 2: Implement the borrowing in a minimal test case
(homogeneous slab, 2G, fixed source). Compare variance per
unit compute against (a) uniform MC, (b) independently
derived CADIS on the same problem.

Pass condition: variance reduction factor matches the
CADIS-derived factor to within implementation noise. The
unified frame produces the same answer as the method-
specific derivation, AND the derivation path is shorter /
more general.

Fail condition: the unified-frame derivation produces
different weights, or the implementation shows no variance
reduction. Document the discrepancy and reject the frame
or refine.

## Literature path

- Kac, "On distributions of certain Wiener functionals"
  (1949) — original Feynman-Kac theorem
- Spanier & Gelbard, "Monte Carlo Principles and Neutron
  Transport Problems" (1969) — path-measure treatment of MC
- Lux & Koblinger, "Monte Carlo Particle Transport Methods"
  (1991) — full path-measure formulation with variance
  reduction theory
- Cercignani, "The Boltzmann Equation and Its Applications"
  (1988) — path-integral derivation of the Boltzmann
  equation
- CADIS / FW-CADIS literature (Wagner, Haghighat) — the
  deterministic-estimated-importance approach in MC

## Transferable pattern

Look for this reformulation when you see:

1. Multiple methods for computing the same quantity,
   implemented independently
2. Each method "feels" like it should borrow from the
   others but the code doesn't
3. Convergence or variance analyses developed separately
   per method
4. The underlying physics is stochastic (or can be cast
   stochastically)

The Feynman-Kac reformulation reveals that the methods are
discretizations of one object. Borrowings become natural.

Apply to (candidate targets):

- Kinetics: deterministic (diffusion / transport with time
  derivative) vs stochastic (neutron noise, point kinetics
  with fluctuation) methods
- Uncertainty quantification: deterministic sensitivity
  (adjoint) vs stochastic (Monte Carlo sampling) as two
  quadratures on the same perturbation path space
- Burnup / depletion: deterministic matrix exponentiation
  vs Monte Carlo nuclide sampling as two discretizations of
  the same nuclide trajectory measure
