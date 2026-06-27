r"""SN-specific transport operators — the WDD/Morel–Montry sweep leaves.

The :class:`~orpheus.numerics.operator.LinearOperator`\ s that are specific to
the discrete-ordinates spatial sweep:

* the streaming leaf :math:`L`
  (:class:`~orpheus.sn.operators.streaming.StreamingOperator`),
* the within-group resolvent :math:`L+C` whose ``.solve`` *is* the sweep
  (:class:`~orpheus.sn.operators.streaming.InvertibleOperator`),
* the boundary leaf :math:`B`
  (:class:`~orpheus.sn.operators.boundary.SNBoundaryOperator`).

These are the SN residue #261 deliberately kept out of ``transport/operators/``:
the cross-method *reaction* operators (collision :math:`M[\sigma_t]`, fission,
scattering) live at L2 ``transport`` because every method collides / fissions /
scatters, but the *inversion of streaming* (the sweep stencil + the per-ordinate
curvilinear closures) is method-specific and stays here at L3 ``sn``. Mirrors the
role-keyed layout of :mod:`orpheus.transport.operators`.
"""
