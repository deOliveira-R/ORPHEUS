Fuel Behaviour — 1-D Radial Thermo-Mechanics
===============================================

Reference for the :mod:`orpheus.fuel` package — a 1-D radial
thermo-mechanical solver for a single fuel pin (pellet + gap + clad).
The solver couples steady-state heat conduction with linear
thermo-elasticity and an iterative gap-conductance closure.

See :ref:`theory-fuel-behaviour` for the governing equations (Fourier
conduction with a volumetric heat source, displacement form of the
equilibrium equations, Ross–Stoute gap model) and the rationale for
the coupling scheme.

Solver
------

.. automodule:: orpheus.fuel.solver
   :members:
   :undoc-members:
   :show-inheritance:
   :noindex:
