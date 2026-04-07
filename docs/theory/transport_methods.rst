.. _theory-transport-methods:

=================
Transport Methods
=================

Deterministic transport solvers for heterogeneous lattice cells.
Each method discretises the angular variable differently and uses a
geometry-specific kernel or sweep algorithm:

- :ref:`theory-collision-probability` — integral transport via the
  :math:`P_{ij}` matrix (slab, cylindrical, spherical kernels).
- :ref:`theory-discrete-ordinates` — differential transport via angular
  quadrature and spatial sweeps (Cartesian 1-D / 2-D, curvilinear 1-D).
- :ref:`theory-method-of-characteristics` — characteristic ray tracing
  with flat-source approximation (2-D pin cell).
- :ref:`theory-monte-carlo` — stochastic transport via Woodcock delta-tracking
  with analog absorption and weight-based population control.

.. toctree::
   :maxdepth: 2

   collision_probability
   discrete_ordinates
   method_of_characteristics
   monte_carlo
