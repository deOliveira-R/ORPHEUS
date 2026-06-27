r"""SN boundary-condition realization — descriptor → operator.

The L3 layer that turns a geometry-declared
:class:`~orpheus.geometry.boundary.BoundaryTraceLaw` into a concrete SN
:class:`~orpheus.numerics.operator.LinearOperator`:

* :mod:`~orpheus.sn.boundary.realizer` — :class:`SNBoundaryRealizer`
  (isinstance-dispatch realizing each ``BoundaryTraceLaw`` concrete) + the
  ``realize_recursively`` walker for composed (rank-N) laws.
* :mod:`~orpheus.sn.boundary.angular` — the angular primitives the realizer
  consumes: :class:`AngularAverageOperator` (white BC) and
  :class:`IncomingSourceOperator` (prescribed inflow).

This is one of the two witnesses to the not-yet-minted ``TransportMethod``
Protocol (the other is the homogenization method-layer in
:mod:`orpheus.transport.mesh.material_mesh`); when a second method's realizer
ships, the walker moves to ``geometry/boundary/`` and the realizer registry
gains its first production consumer (deferred — #219, defer-until-2). Mirrors
the role-keyed layout of :mod:`orpheus.transport`.
"""
