r"""L2 typed flux displacements — the difference vector space :math:`V` of the
affine flux space :math:`A`.

A **displacement** :math:`\Delta\psi = \psi^{(i)} \ominus \psi^{(i-1)}` is the
iterate INCREMENT — the tangent vector between two flux states. It shares its
flux sibling's storage family, :class:`~orpheus.numerics.space.FunctionSpace`,
metric, and units, but is a DISTINCT physical kind: a *step*, not a *state*.
This is the affine completion of the field role grid — the dual of the residual
column (a residual is ``source ⊖ source`` crossing into rate-density; a
displacement is ``flux ⊖ flux`` staying in flux units but changing role from
state to tangent).

Role grid (field vocabulary, issues #205 / #201 / #208; axes per the
#290 P2.5 coherence ruling)
=======================================================================

Three orthogonal axes — "Boundary" is a LOCUS qualifier, never a fourth
family (bulk leaf ``<Family><Role>``, boundary leaf
``<Family>Boundary<Role>``)::

    locus {Bulk, Boundary} × family {Angular, Scalar, Moment}
                           × role {Flux, SourceSink, Residual, Displacement}

       flux         →  orpheus.transport.fields         (…Flux)
       source       →  orpheus.transport.source_sinks   (…SourceSink)
       residual     →  orpheus.transport.residuals      (…Residual)   [rate-density defect]
       displacement →  orpheus.transport.displacements  (…Displacement) [flux tangent]

Per the Field ABC's Layer-1 class-identity gate, ``flux + flux`` is forbidden
(an affine space has no origin); the legitimate operations are ``flux ⊖ flux →
displacement`` (minted here, via :class:`FluxRole`) and ``flux ⊕ displacement →
flux`` (the torsor action). The displacements themselves ARE a vector space
(``+`` / scalar ``·`` / unary ``-`` closed) — inherited from
:class:`~orpheus.numerics.field.Field`.

See :class:`~orpheus.transport.displacements._displacement.Displacement` for the
shared convergence diagnostics (``contraction_ratio`` / ``true_error_estimate``
/ ``where_largest``) and
:class:`~orpheus.transport.fields._flux_role.FluxRole` for the affine algebra.
"""
from __future__ import annotations

from orpheus.transport.displacements._displacement import Displacement
from orpheus.transport.displacements.angular_displacement import AngularDisplacement
from orpheus.transport.displacements.scalar_displacement import ScalarDisplacement
from orpheus.transport.displacements.moment_displacement import MomentDisplacement
from orpheus.transport.displacements.angular_boundary_displacement import AngularBoundaryDisplacement
from orpheus.transport.displacements.scalar_boundary_displacement import (
    ScalarBoundaryDisplacement,
)


__all__ = [
    "Displacement",
    "AngularDisplacement",
    "ScalarDisplacement",
    "MomentDisplacement",
    "AngularBoundaryDisplacement",
    "ScalarBoundaryDisplacement",
]
