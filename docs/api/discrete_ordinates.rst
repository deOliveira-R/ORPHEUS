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

.. automodule:: orpheus.sn.quadrature
   :members:
   :undoc-members:
   :show-inheritance:
   :noindex:

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
