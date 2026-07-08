r"""SN-specific transport operators — the WDD/Morel–Montry sweep leaves.

The :class:`~orpheus.numerics.operator.LinearOperator`\ s that are specific to
the discrete-ordinates spatial sweep:

* the streaming leaf :math:`L`
  (:class:`~orpheus.sn.operators.streaming.StreamingOperator`),
* the within-group resolvent :math:`L+C` whose ``.solve`` *is* the sweep
  (:class:`~orpheus.sn.operators.streaming.InvertibleOperator`),
* its inverse :math:`(L+C)^{-1}` in operator form — what ``(L+C).inverse()``
  returns (:class:`~orpheus.sn.operators.sweep_operator.SweepOperator`); the
  inverse-as-operator carve (#226) so ``K = A_loss.inverse() @ F`` reads like the
  math,
* the boundary leaf :math:`B`
  (:class:`~orpheus.sn.operators.boundary.SNBoundaryOperator`),
* the **ψ½ coupled-block operators** — System B of the augmented within-group
  system (a 2×2 block operator over the transport bulk⊕trace and the
  radial-characteristic ray): the radial straight-characteristic self-block
  :math:`A_{BB}`
  (:class:`~orpheus.sn.operators.radial_characteristic.RadialCharacteristicOperator`),
  the cell-local ray→bulk seed injection :math:`A_{AB}`
  (:class:`~orpheus.sn.operators.radial_characteristic.RadialCharacteristicSeeding`),
  and the ray boundary block :math:`B_b`
  (:class:`~orpheus.sn.operators.boundary.RadialCharacteristicBoundaryOperator`).

These are the SN residue #261 deliberately kept out of ``transport/operators/``:
the cross-method *reaction* operators (collision :math:`M[\sigma_t]`, fission,
scattering) live at L2 ``transport`` because every method collides / fissions /
scatters, but the *inversion of streaming* (the sweep stencil + the per-ordinate
curvilinear closures) is method-specific and stays here at L3 ``sn``. Mirrors the
role-keyed layout of :mod:`orpheus.transport.operators`.
"""
