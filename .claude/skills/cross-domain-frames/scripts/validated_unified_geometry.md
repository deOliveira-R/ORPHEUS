---
status: validated
project: ORPHEUS
modules: [geometry, sn, cp, moc]
date_validated: pre-chat
primary_criteria: [structurally-simpler, structure-exposing]
secondary_criteria: [expressive]
requires_validation: false
---

# Precedent: Unified 1D geometry via topology + metric

## Problem (current / minimal formulation)

Slab, annulus (cylindrical shell), and hollow sphere are
three distinct geometries in 1D transport:

- Slab: flat metric, volume element dV = dx
- Annulus: dV = 2πr dr, annular topology (inner + outer
  radius)
- Hollow sphere: dV = 4πr² dr, spherical topology (inner +
  outer radius)

Historical implementation pattern (pre-reformulation): three
code paths per solver (SN, CP, MOC), each with its own sweep
kernel, its own boundary handling, its own verification tests.

Symptoms of the wrong-frame formulation:

- Smell 1 fires strongly: three near-identical sweep kernels
- Smell 9 fires: boundary handling as three distinct special
  cases
- Additionally: limit transitions (slab → annulus as
  R_inner → ∞ with thickness fixed; solid sphere ← hollow
  sphere as R_inner → 0) are physically meaningful but not
  checkable in the current code because the geometries live
  in separate code paths

## Structural trigger

From reference.md Part A:

- A.1 Topology: "multiple geometry variants with shared
  boundary behavior"
- A.1 Differential geometry: "curvilinear coordinates with
  metric-dependent volume elements"

From reference.md Part C:

- Smell 1 (strongly)
- Smell 9

## Frame

Topology (manifold-with-boundary structure) + differential
geometry (metric / volume element parameterization). Core
objects: a parameterized 1D radial manifold whose
dimensionality and boundary structure are explicit
parameters of the formulation rather than implicit type
distinctions.

## Reformulation (sketch)

A single 1D radial transport kernel parameterized by
geometric descriptors. The specific parameterization used in
ORPHEUS should be documented here; a canonical candidate is:

- `dim ∈ {0, 1, 2}` controlling the volume element
  `dV = r^dim dr` (with 2π, 4π prefactors absorbed into
  normalization)
  - dim = 0: slab (flat)
  - dim = 1: cylindrical shell
  - dim = 2: spherical shell
- `R_inner, R_outer`: inner and outer radii, with
  `R_inner = 0` collapsing annulus → solid cylinder and
  hollow sphere → solid sphere
- Axial-symmetry constraint on the angular quadrature, which
  becomes a parameter of the quadrature rule rather than a
  per-geometry hardcoding

One sweep kernel. Geometry enters via:

- Volume element `r^dim` in the spatial integration
- Surface-area factor `r^dim` at interfaces
- Axial-symmetry constraint on the angular quadrature

NOTE: The specific implementation in ORPHEUS may use a
different parameterization. Edit this section to match the
actual code. The high-level pattern (geometry as parameter,
not type) is the invariant.

## Elegance payoff

- **Structure-exposing**: the common structure (1D radial
  transport under axial symmetry) becomes the primary
  object. Geometry is a parameter of the formulation, not a
  type distinction. The shared physics is made explicit.
- **Structurally-simpler**: one kernel replaces three; test
  suites collapse from three duplicated structures to one
  parameterized structure (with geometry as a fixture).
- **Expressive**: physical limits (solid from hollow sphere;
  curvature limit of cylindrical shell) become expressible
  in the formulation and checkable as regression tests.

## Concrete first test

Retrospective first test that would have predicted success:

Implement the unified kernel. Solve the same homogeneous 1D
transport problem on (a) slab with thickness L, (b)
cylindrical shell with R_outer − R_inner = L and R_inner
small enough that curvature effects are negligible. The
solutions should agree to discretization error. Measure
agreement; measure code-line count vs three-kernel
implementation.

Pass condition: solutions agree; code-line count is less
than sum of three-kernel implementation; all three original
verification suites pass under the unified kernel.

## Literature path

- Lewis & Miller, "Computational Methods of Neutron
  Transport" (1984) — per-geometry derivations in separate
  chapters (typical of the wrong-frame approach)
- Modern transport texts (Hébert; Duderstadt-Hamilton) also
  treat geometries separately; the unified treatment appears
  to be an ORPHEUS novelty as of pre-chat
- Internal ORPHEUS documentation — actual parameterization
  and implementation

## Transferable pattern

Look for this reformulation when you see:

1. N code paths for N variants of a geometry / coordinate
   system
2. The variants share the same underlying physics
3. A small number of parameters (dimension, genus, metric
   signature, radii) distinguishes them
4. Limit transitions between variants are physically
   meaningful but not checkable in the current code

The unification exposes the shared physics as the primary
object and relegates geometry to a parameter.

Apply to (candidate targets):

- 2D: Cartesian vs cylindrical (r-z) vs spherical (r-θ) via
  analogous multi-dimensional parameterization
- Steady-state vs time-dependent as a limit of the kinetic
  problem with time-derivative parameter → 0
- Multi-group → 1G as a limit with scattering kernel
  becoming diagonal
- Fuel pin vs assembly vs core via hierarchical topology
  with explicit scale parameter
