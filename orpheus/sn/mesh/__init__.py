r"""SN mesh layer — the method-space and the reduced streaming operators the
Problem consumes.

* The Problem itself — :class:`~orpheus.sn.problem.SNProblem`, the axis-primary
  data hub: a :class:`~orpheus.transport.mesh.material_mesh.MaterialMesh`
  augmented with the SN angular + spatial discretization machinery (the
  precomputed streaming stencil, the α / geometry-factor / Morel–Montry weights),
  the posed operators and the pencil — lives one level UP at
  :mod:`orpheus.sn.problem`.  ⛔ Until #412 (2026-09-18) it was ``SNMesh`` in
  this package's ``augmented_mesh.py`` — "augmented mesh" is the cross-method
  name for a mesh carrying a method's machinery (the diffusion family still
  spells its hub :mod:`~orpheus.diffusion.augmented_mesh`), and the SN hub had
  outgrown it: it is the Problem, not a mesh.
* :mod:`~orpheus.sn.mesh.reduced_operator` — the reduced streaming operators
  (the per-chart streaming factories the Problem's stencil is built from).
* :mod:`~orpheus.sn.mesh.method_space` — :class:`SNMethodSpace`, the realizer's
  argument (mesh + quadrature + trace + face); the precursor to the not-yet-minted
  ``TransportMethod`` Protocol (#219, defer-until-2).

Mirrors :mod:`orpheus.transport.mesh` (the shared ``MaterialMesh`` middle type
this layer extends).
"""
