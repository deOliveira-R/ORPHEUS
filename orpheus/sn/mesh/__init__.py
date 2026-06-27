r"""SN mesh layer — the augmented mesh + the method-space over it.

* :mod:`~orpheus.sn.mesh.augmented_mesh` — :class:`SNMesh`, the axis-primary
  *augmented mesh*: a :class:`~orpheus.transport.mesh.material_mesh.MaterialMesh`
  augmented with the SN angular + spatial discretization machinery (the
  precomputed streaming stencil, the α / geometry-factor / Morel–Montry weights).
  Every transport method carries an ``augmented_mesh`` — that is the cross-method
  name for this file.
* :mod:`~orpheus.sn.mesh.method_space` — :class:`SNMethodSpace`, the realizer's
  argument (mesh + quadrature + trace + face); the precursor to the not-yet-minted
  ``TransportMethod`` Protocol (#219, defer-until-2).

Mirrors :mod:`orpheus.transport.mesh` (the shared ``MaterialMesh`` middle type
this layer extends).
"""
