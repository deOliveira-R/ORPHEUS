.. _theory-index:

==========================================
Theory and Derivations
==========================================

ORPHEUS theory pages document the mathematics behind every solver and
every reference. They are organised into three branches that
correspond to three distinct **kinds of mathematical artefact**, each
with a different audience, lifecycle, and V&V claim:

.. list-table::
   :header-rows: 1
   :widths: 18 22 30 30

   * - Branch
     - What it documents
     - Production status
     - Audience
   * - :ref:`theory-infrastructure`
     - Cross-cutting foundations: cross-section pipeline, homogeneous
       infinite-medium baseline, transport-method overview.
     - Stable invariants
     - All sessions touching any solver.
   * - :ref:`theory-transport-methods`
     - **Discrete (production) solvers** the user runs for analysis:
       collision probability, discrete ordinates, method of
       characteristics, Monte Carlo (and the support solvers:
       diffusion, fuel behaviour, thermal hydraulics, kinetics).
     - Production-grade. End-user API.
     - Sessions modifying or extending a production solver.
   * - :ref:`theory-reference-solvers`
     - **Continuous (Pillar-2 reference) solvers** that supply the
       analytical / semi-analytical truth values the verification
       suite consumes: Peierls Nyström, trajectory resolvent (Variant
       α), F_N method, singular-eigenfunction expansion, Galerkin
       spectral, Sood-family case registry, plus reserved slots for
       methods queued for implementation.
     - Verification infrastructure. NOT user-facing solvers.
     - Sessions writing a new reference, debugging a verification
       gap, or extending a Pillar-2 truth set.

The split is **load-bearing**, not cosmetic. Discrete solvers and
reference solvers serve different masters:

- Discrete solvers are tuned for **scale** — they must be fast on
  realistic problems with realistic mesh refinement, multi-region
  pin cells, full multi-group cross-section data. Their error budget
  is a discretisation error that decreases under refinement.
- Reference solvers are tuned for **accuracy** — they must be
  arbitrary-precision evaluable on the **specialised problems** that
  yield analytical or semi-analytical solutions. Their error budget
  is a numerical-quadrature floor, typically 10⁻⁸ or better.

Mixing the two confuses readers and corrupts the V&V vocabulary.
This page is the canonical entry point that prevents the confusion.


.. _theory-infrastructure:

Infrastructure
==============

Cross-cutting theory shared across every solver — cross-section data,
the homogeneous infinite-medium reference, and the high-level taxonomy
that places each transport method in the V&V hierarchy.

.. toctree::
   :maxdepth: 1

   cross_section_data
   homogeneous
   verification


.. toctree::
   :maxdepth: 1
   :caption: Discrete (production) solver theory

   transport_methods
   diffusion_1d
   fuel_behaviour
   thermal_hydraulics
   reactor_kinetics


.. toctree::
   :maxdepth: 1
   :caption: Continuous (reference) solver theory

   reference_solvers
