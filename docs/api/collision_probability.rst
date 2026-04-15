Collision Probability Solvers
=============================

Reference for the :mod:`orpheus.cp` package — a first-flight collision
probability (CP) solver for 1-D slab and cylindrical geometries with
multi-group scattering and fission. See :ref:`theory-collision-probability`
for the derivation of the :math:`P_{ij}` reciprocity matrix, the
Bickley–Naylor ``Ki_3`` / ``Ki_4`` kernels, and the Gauss–Seidel versus
Jacobi inner-iteration scheme.

Mesh input follows the standard project convention
(:class:`orpheus.geometry.mesh.Mesh1D`) — the CP solver does not own
its own geometry type. Cross-sections come through
:class:`orpheus.data.macro_xs.cell_xs.CellXS`, built from a material-ID
to :class:`orpheus.data.macro_xs.mixture.Mixture` mapping.

Solver
------

.. automodule:: orpheus.cp.solver
   :members:
   :undoc-members:
   :show-inheritance:
   :noindex:
