---
status: validated
project: ORPHEUS
modules: [geometry, sn, cp]
date_validated: pre-chat
primary_criteria: [structurally-simpler, algorithmic-advantage]
secondary_criteria: [structure-exposing, expressive]
requires_validation: false
---

# Precedent: Boundary conditions as a tensor network

## Problem (current / minimal formulation)

Reactor physics boundary conditions — vacuum, reflective,
albedo, periodic — are historically implemented as matrix
operators that couple boundary degrees of freedom to interior
degrees of freedom. For an N-cell problem, the BC matrix can
be N × N (albedo with full angular coupling), O(N) storage,
O(N²) or worse for assembly.

Symptoms of the wrong-frame formulation (elegance smells
firing):

- Smell 2: dense matrix explicitly formed even when the BC
  acts locally per boundary segment
- Smell 9: boundary handling as a special case added to
  interior operator logic
- Smell 1 (partial): distinct code paths per BC type, each
  regenerating its own matrix

## Structural trigger

From reference.md Part A:

- A.2 Tensor networks / tensor-train: "compositional boundary
  structure"
- A.1 Topology: "multiple geometry variants with shared
  boundary behavior"

From reference.md Part C:

- Smell 2 fires
- Smell 9 fires
- Smell 1 fires partially

## Frame

Tensor networks. Core objects: small tensors attached to
compositional elements (one per surface segment, per angular
bin, per energy group), connected by contraction edges whose
widths correspond to the information flowing between elements.
For geometries with branching symmetry, a tree tensor network;
for simple 1D with periodic or linear structure, a matrix
product operator (MPO).

## Reformulation (sketch)

The BC is a network of small tensors, one per compositional element:

[angle] ─ [energy] ─ [segment] ─ ... ─ [domain link]

Each tensor captures one local contribution. The full BC
operator is the contraction of this network. Key properties:

- No dense matrix ever formed; matrix-vector products are
  sequential tensor contractions in the contraction ordering
  dictated by the network topology
- BC type change is a local edit to one tensor node, not
  wholesale matrix regeneration
- Composition across geometries: slab / annulus / sphere
  change the network topology (specifically the "segment" leg
  structure) but not the other tensor legs
- Interior operator and BC operator compose as one network;
  the BC is not a separate object grafted onto the interior

## Elegance payoff

- **Structure-exposing**: the compositional structure of the
  BC — previously implicit in the matrix form — becomes the
  primary object. The tensor network diagram is the
  definition, not a visualization aid.
- **Expressive**: a new BC type is written as a small tensor
  plugged into the network, not a new matrix class with its
  own assembly routine.
- **Structurally-simpler**: one framework, four-plus BC
  types, three geometries — previously 12+ code paths,
  collapsed to one framework with local variations.
- **Algorithmic-advantage**: avoids O(N²) storage and O(N³)
  factorization; memory and compute scale with the network's
  maximum bond dimension, which is small (O(1) or O(log N))
  for physical BCs.

## Concrete first test

The validation was done pre-chat. Retrospective first-test
that would have predicted success:

Take a known-working dense-matrix BC implementation for a
specific case (e.g., albedo BC on cylindrical SN). Construct
the tensor network decomposition. Verify that contracting the
network reproduces the dense matrix entry-by-entry. Measure
storage and construction time for both.

Pass condition: entries match to machine precision; tensor
network storage and construction is asymptotically cheaper.

## Literature path

- Orús, "A practical introduction to tensor networks" (2014)
  — MPS / MPO / PEPS introduction
- Evenbly & Vidal — hierarchical tensor networks for
  branching systems
- Internal ORPHEUS documentation — specific parameterization
  and algorithmic details

## Transferable pattern

Look for this reformulation when you see:

1. Multiple variants of the same operator (BC type A, B, C,
   ...; scattering kernel type anisotropic / isotropic;
   quadrature product structure)
2. Each variant implemented as a separate matrix or function
3. Each variant shares structural elements with the others
4. Dense matrix formation when the operator acts locally

The tensor network exposes the shared structure as a common
network, with variants differing only by local tensor
content.

Apply to (candidate targets):

- Scattering kernels (Σₛ decomposition by anisotropy order)
- Quadrature rules (product structure; Lebedev is not a
  product quadrature and would NOT decompose this way)
- Coupled multiphysics interface operators (field coupling
  tensors)
- Response matrices for assembly-homogenization
