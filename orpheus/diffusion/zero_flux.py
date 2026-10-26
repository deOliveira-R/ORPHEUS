r"""The zero-flux condition, the one exterior condition only diffusion admits (``zero_flux.py`` of ``orpheus/diffusion``).

**The mathematics.** The textbook Dirichlet idealization of a bare surface pins the scalar flux itself to zero at the surface,
:math:`\phi_\Gamma = 0`. In the partial-current basis of the scalar trace, :math:`(J^+, J^-)`, the :math:`P_1` closure
dictionary :math:`\phi_\Gamma = 2(J^+ + J^-)` turns this into the partial-current statement

.. math::

    \phi_\Gamma = 0 \;\Longleftrightarrow\; J^- = -J^+ \;\Longleftrightarrow\; \mathcal{A} = -1 ,

a member of the albedo family :math:`J^- = \mathcal{A}\,J^+` with a NEGATIVE amplitude. It is an excellent test and reference
condition (every analytic sine anchor in ``tests/gates/diffusion`` is derived under it: the bare-slab mode vanishes exactly at
the surface, with no extrapolation length) but it is not a physical transport boundary: no material or void outside the domain
produces it. It must not be confused with vacuum, which is zero INCOMING current (:math:`J^- = 0`, a point of the shared
response cone), where zero flux requires a negative incoming current that no non-negative angular flux can produce.

**The type.** :class:`ZeroFluxCondition` is a :class:`~orpheus.transport.boundary_conditions.BoundaryCondition` that is NOT a
:class:`~orpheus.transport.boundary_conditions.WallCondition`: a method-only condition, with no fields (so every instance is
equal to every other and has one digest), ``return_fraction == -1.0`` (the partial-current statement above, outside the
sub-Markov cone by definition, which is why the cone is ``WallCondition``'s parse and not an invariant of the base class),
``face_local == True`` (it is stated in the face's own frame, so every isometry carries it from face to face and one tag may
hold two faces) and an empty source. It defines no arithmetic: ``ZeroFluxCondition() + x`` is Python's own ``TypeError``
(law L-F), because a sum of a Dirichlet condition and a response has no meaning.

**Who admits it.** Diffusion's table registers it (``DiffusionZeroFlux`` in
:mod:`orpheus.diffusion.boundary_realizer`, tier 1, :math:`\mathcal{A} = -1`) and nothing else does. A table of the
non-negative cone refuses it through the unregistered-unit door on the state-cone axis (a return fraction below zero), which
keeps that refusal without SN, CP, MoC or MC importing this module: an L3 package may not import another L3 package.

**Why the name is ``ZeroFluxCondition`` and the module is here.** The tree once held a class ``ZeroFluxBoundary`` in
``orpheus/geometry/boundary/zero_flux.py`` (the package retired in WF2.4), which production files and tests of the CP, MoC, MC,
derivations and specification slices named until they were re-spelled after this one. A content key names an orpheus type by its
qualname, so a second class under the old name would have made two nameable classes share one, and the qualname gate
(``tests/gates/test_content_identity_census.py::test_no_two_nameable_orpheus_classes_share_a_qualname``) takes no exemption (an exemption would be a dual spelling under
another name). The new type was therefore minted under a fresh name, and keeps it. The module sits in ``diffusion`` and not in a
module ``boundary_conditions.py`` beside the package ``orpheus/transport/boundary_conditions/``: that would be the homonym the
layering rulings refused for the package name.

**Gotchas.**

* ``ZeroFluxCondition() == ZeroFluxCondition()``: no fields, so content equality is trivial, and it is hashable.
* It realizes as ``-IdentityOperator()`` in diffusion and is refused everywhere else; a ``BoundaryRefused`` from another
  method's table is therefore the expected answer, not a defect.
* It is not a response network: ``WallCondition.zero_flux`` does not exist and a mixed wall cannot carry it.

Theory: :ref:`bc-architecture-of-record`.

Layer: ``diffusion`` (L3). Imports ``transport`` (the condition base and the empty source).
"""

from __future__ import annotations

from dataclasses import dataclass

from orpheus.transport.boundary_conditions import BoundaryCondition, InflowTerms

__all__ = ["ZeroFluxCondition"]


@dataclass(frozen=True, eq=False)
class ZeroFluxCondition(BoundaryCondition):
    r"""Zero scalar flux at the face, :math:`\phi_\Gamma = 0`, i.e. :math:`J^- = -J^+` (:math:`\mathcal{A} = -1`).

    No fields. The condition lies outside the sub-Markov cone (its return fraction is negative), which only a method whose
    trace state is a signed partial current can hold; see the module docstring.
    """

    @property
    def return_fraction(self) -> float:
        r""":math:`\rho = J^-/J^+ = -1` (the partial-current statement of :math:`\phi_\Gamma = 0`)."""
        return -1.0

    @property
    def source(self) -> InflowTerms:
        """No affine part."""
        return InflowTerms.empty()

    @property
    def face_local(self) -> bool:
        """Stated in the face's own frame, since every isometry carries it from face to face."""
        return True
