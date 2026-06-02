Discrete Ordinates Solvers
==========================

Reference for the :mod:`orpheus.sn` package — the discrete-ordinates
(S\ :sub:`N`) transport solvers. Two execution paths share the same
quadrature and geometry layer:

* **Source iteration** via diamond-difference sweeps
  (:mod:`~orpheus.sn.sweep`) — the default path used by
  :class:`~orpheus.sn.solver.SNSolver`.
* **Krylov** via an explicit transport operator
  (:mod:`~orpheus.sn.operator`) — forms ``T: ψ → T·ψ`` as a
  :class:`scipy.sparse.linalg.LinearOperator` so scipy's BiCGSTAB / GMRES
  can drive the inner solve directly.

The theory pages cover the diamond-difference discretisation, the
angular-redistribution term for curvilinear geometry (Bailey et al.
2009), and the source-iteration / Krylov trade-offs.

Solver
------

.. automodule:: orpheus.sn.solver
   :members:
   :undoc-members:
   :show-inheritance:
   :noindex:

Geometry
--------

.. automodule:: orpheus.sn.geometry
   :members:
   :undoc-members:
   :show-inheritance:
   :noindex:

Angular Quadrature
------------------

The angular quadrature now lives in the method-agnostic
:mod:`orpheus.numerics.quadrature` package (re-exported as
:class:`orpheus.sn.Quadrature` for the SN solver's convenience). The
five legacy SN-only subclasses (``AngularQuadrature``,
``GaussLegendre1D``, ``LebedevSphere``, ``LevelSymmetricSN``,
``ProductQuadrature``) collapsed into the single
:class:`~orpheus.numerics.quadrature.Quadrature` value type with
``classmethod`` factories:

* :meth:`~orpheus.numerics.quadrature.Quadrature.gauss_legendre`
  — 1-D Gauss-Legendre on :math:`\mu \in [-1, 1]` (slab / curvilinear
  radial).
* :meth:`~orpheus.numerics.quadrature.Quadrature.level_symmetric`
  — :math:`O_h`-invariant level-symmetric :math:`S_N` (Carlson &
  Lathrop 1968).
* :meth:`~orpheus.numerics.quadrature.Quadrature.lebedev`
  — :math:`O_h`-invariant Lebedev sphere quadrature.
* :meth:`~orpheus.numerics.quadrature.Quadrature.product`
  — Gauss-Legendre :math:`(\mu)` :math:`\times` equispaced
  :math:`(\phi)` product rule.

The per-ordinate angular data is exposed through the cached
:attr:`~orpheus.numerics.quadrature.Quadrature.octants` partition and
the :meth:`~orpheus.numerics.quadrature.Quadrature.spherical_harmonics`
/ :meth:`~orpheus.numerics.quadrature.Quadrature.reflection_index`
methods. The selection driver is
:func:`~orpheus.numerics.quadrature.select_quadrature`, backed by the
:data:`~orpheus.numerics.quadrature.quadrature_registry`. The full
mathematical narrative — the level-symmetric construction, the
selection criterion, and the product-rule cosine layout — lives at
:ref:`discrete-measures` and in the per-module docstrings, accessible
via the standard ``orpheus.numerics.quadrature`` import path. (The
package carries rich ``.. math:: :label:`` docstrings, so it is
cross-referenced here rather than ``automodule``-rendered, to avoid
duplicate-label collisions with the theory pages.)

Transport Sweep
---------------

.. automodule:: orpheus.sn.sweep
   :members:
   :undoc-members:
   :show-inheritance:
   :noindex:

Direct Transport Operator
-------------------------

The :mod:`~orpheus.sn.operator` module forms the transport operator
``T`` explicitly for the Krylov inner solver. The
:class:`~orpheus.sn.operator.EquationMap` dataclass indexes which
``(ordinate, cell)`` pairs are unknowns — ordinates below the equator
and incoming faces at reflective boundaries are excluded and restored
by reflection.

.. automodule:: orpheus.sn.operator
   :members:
   :undoc-members:
   :show-inheritance:
   :noindex:
