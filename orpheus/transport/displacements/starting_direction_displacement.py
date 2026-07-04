r"""Starting-direction flux displacement :math:`\Delta\psi_{1/2}`.

The difference-space (tangent) sibling of
:class:`~orpheus.transport.fields.starting_direction_flux.StartingDirectionFlux`
— the increment between two ψ½ iterates (the ``starting_direction``
member of an augmented iterate's
:class:`~orpheus.transport.timed_full_field.TimedFullField`). Stored on
the shared R12a-keyed
:class:`~orpheus.numerics.spaces.starting_direction_space.StartingDirectionSpace`
(``mesh.starting_direction_space``); a DISTINCT class. Born ONLY from
``StartingDirectionFlux ⊖ StartingDirectionFlux`` (the
:class:`~orpheus.transport.fields._flux_role.FluxRole` mint, resolved
through the Rep-keyed
:meth:`~orpheus.transport.displacements._displacement.Displacement.sibling_of`
registry — this leaf's existence is FORCED by the composite torsor
algebra: an augmented flux composite's ``ψ₂ − ψ₁`` mints a displacement
per member, the seed block included). See
:class:`~orpheus.transport.displacements._displacement.Displacement`.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import ClassVar

from orpheus.numerics.units import ANGULAR_FLUX_UNITS, Unit
from orpheus.transport.displacements._displacement import Displacement
from orpheus.transport.fields._bases import StartingDirectionField

__all__ = ["StartingDirectionDisplacement"]


@dataclass(frozen=True, eq=False, kw_only=True, repr=False)
class StartingDirectionDisplacement(Displacement, StartingDirectionField):
    r"""L2 starting-direction flux displacement on the shared seed space.

    Like every displacement leaf it shares its flux sibling's space
    (the tangent vector lives in the flux's own function space — here
    the all-zero ghost metric rides along); the CLASS identity is the
    role gate.
    """

    #: Same units as :class:`StartingDirectionFlux` (``1/(cm²·s·sr)`` —
    #: a ψ½ flux-difference).
    UNITS: ClassVar[Unit] = ANGULAR_FLUX_UNITS
