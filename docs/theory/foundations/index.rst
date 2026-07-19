.. _theory-foundations:

===========
Foundations
===========

The mathematics **every** method shares. A concept lives here once; each
method's chapter carries only *its realization* plus a link back.

This part exists because the shared content is genuinely shared **in code,
not by analogy**: :class:`~orpheus.transport.operators.MultiplicationOperator`,
:class:`~orpheus.transport.operators.IsotropicScattering` and the fission
operator are the *same Python objects* imported by S\ :sub:`N`, diffusion and
the infinite-medium solver from :mod:`orpheus.transport.operators`. What
varies between methods is how streaming is represented — not what collision,
scattering, and fission *are*.

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Page
     - What it settles
   * - :doc:`/theory/foundations/path_integral`
     - **The root of the corpus** — the one object all methods discretize
       (the sum over neutron histories), what is invariant
       (:math:`C, S, F`), what varies, the three axes on which methods
       differ, and where each lands. The parent of
       :doc:`/theory/methods/index`. *(Scaffold; authored at Phase H.)*
   * - :doc:`/theory/foundations/operator_algebra`
     - The operator algebra itself: :math:`A = L + C - S - B`, posed
       :math:`A\psi = \tfrac{1}{k}F\psi` (eigenvalue) or :math:`A\psi = q`
       (fixed source). :math:`B` is a first-class **sibling**, not folded
       into :math:`L`; :math:`(L+C)` is the sub-composite whose inverse
       **is** the transport sweep.
   * - :doc:`/theory/foundations/discretization`
     - How a continuous conservation law becomes a finite algebraic system:
       the **cell-balance** invariant (sinks = sources) and its **closures** —
       Step (upwind), Diamond Difference (central), Linear Discontinuous
       (linear upwind), derived once and dimension-agnostic (the same closure
       in space **and** angle).
   * - :doc:`/theory/foundations/frame`
     - Frames, and why projection is **Petrov-Galerkin**: the trial/test
       split, the adjoint, and the realizations (spherical-harmonics
       Galerkin; homogenization and energy condensation as Petrov-Galerkin).
   * - :doc:`/theory/foundations/boundary_conditions`
     - The boundary law :math:`B` as an operator: trace realization,
       reflective / vacuum / white, and the extraction criterion.
   * - :doc:`/theory/foundations/cross_section_data`
     - The cross-section pipeline: mixtures, multigroup data, condensation.
   * - :doc:`/theory/foundations/discrete_measures`
     - Quadrature and measure: axes, weights, and integration.
   * - :doc:`/theory/foundations/spherical_harmonics`
     - The angular basis and the addition theorem.
   * - :doc:`/theory/foundations/structured_geometry`
     - Meshes and structured geometry.
   * - :doc:`/theory/foundations/infinite_medium`
     - The 0-D infinite-medium (:math:`k_\infty`) baseline — the analytical
       anchor every method must reproduce. **Not** spatial homogenization;
       for that see :doc:`/theory/foundations/frame`.

.. toctree::
   :maxdepth: 2

   path_integral
   operator_algebra
   discretization
   frame
   boundary_conditions
   cross_section_data
   discrete_measures
   spherical_harmonics
   structured_geometry
   infinite_medium
