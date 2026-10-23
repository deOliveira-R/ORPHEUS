r"""SN boundary-condition realization — declaration to operator.

The L3 layer that turns a declared exterior (the conditions of
:mod:`orpheus.transport.boundary_conditions`, joined to the faces of a problem by a ``Binding``) into the concrete SN face
operators :math:`\Gamma_+ \to \Gamma_-` consumable by the sweep and the Krylov path:

* :mod:`~orpheus.sn.boundary.realizer` — :data:`SN_REALIZERS`, the sealed
  :class:`~orpheus.transport.boundary_conditions.RealizerTable` of SN: the ONE door that answers "is this declaration admitted
  on this face" and "realize it", filled by five realizer classes (``SNZero``, ``SNSpecularChain``, ``SNDiffuseChain``,
  ``SNSource``, ``SNDecks``).
* :mod:`~orpheus.sn.boundary.realized_operators` — the typed images the table returns: :class:`ZeroMap`,
  :class:`ScaledOrdinatePermutation`, :class:`RankOneReemission` and :class:`RealizedSum`, all
  :class:`RealizedBoundary`. The facts the consumers read (the return fraction, whether the face relabels ordinates, whose
  outflow it consumes, whether it is the zero map) are properties of the TYPE, derived from its parameters.
* :mod:`~orpheus.sn.boundary.angular` — the angular primitives the Lambertian image consumes:
  :class:`~orpheus.sn.boundary.angular.PartialCurrentOperator` and :class:`~orpheus.sn.boundary.angular.IsotropicEmissionOperator`,
  the two links of the factored Lambertian.
* :mod:`~orpheus.sn.boundary.errors` — the two certifications of the realizer's own index sets (ERR-041, ERR-047).
* :mod:`~orpheus.sn.boundary.bound_compat` — ``_BoundBoundaryOperator``, the retired shim that paired a realized operator with
  the law it came from. No consumer imports it; WF2.4 deletes the file.

This seam was one of the two recorded witnesses of the :class:`~orpheus.transport.method.TransportMethod` Protocol (the other: the
homogenization method-layer in :mod:`orpheus.transport.mesh.material_mesh`). Mirrors the role-keyed layout of
:mod:`orpheus.transport`.
"""

from .realized_operators import (
    RankOneReemission,
    RealizedBoundary,
    RealizedSum,
    ScaledOrdinatePermutation,
    ZeroMap,
)
from .realizer import SN_REALIZERS

__all__ = [
    "RankOneReemission",
    "RealizedBoundary",
    "RealizedSum",
    "SN_REALIZERS",
    "ScaledOrdinatePermutation",
    "ZeroMap",
]
