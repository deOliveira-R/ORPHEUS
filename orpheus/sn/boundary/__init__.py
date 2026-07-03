r"""SN boundary-condition realization — descriptor → operator.

The L3 layer that turns a geometry-declared
:class:`~orpheus.geometry.boundary.BoundaryTraceLaw` into a concrete SN
:class:`~orpheus.numerics.operator.LinearOperator`:

* :mod:`~orpheus.sn.boundary.realizer` — :class:`SNBoundaryRealizer`
  (isinstance-dispatch realizing each ``BoundaryTraceLaw`` concrete).
* :mod:`~orpheus.sn.boundary.angular` — the angular primitives the realizer
  consumes: :class:`AngularAverageOperator` (white BC) and
  :class:`IncomingSourceOperator` (prescribed inflow).

This seam was one of the two recorded witnesses of the
:class:`~orpheus.transport.method.TransportMethod` Protocol (the other:
the homogenization method-layer in
:mod:`orpheus.transport.mesh.material_mesh`) — minted at #290 P7b when
the second method's realizer shipped. The rank-N composition walker
moved to :func:`orpheus.geometry.boundary.realize_recursively`
(method-blind; pass ``realizer=SNBoundaryRealizer()``), and the Wave-5
realizer registry dissolved (you hold the method-mesh → you have its
realizer). Mirrors the role-keyed layout of :mod:`orpheus.transport`.
"""
