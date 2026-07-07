r"""L2 typed transport residuals — every class inherits from
:class:`~orpheus.numerics.field.Field`.

Residuals are the *defect* of a transport balance:
:math:`r = (L + C - S - F)\,\psi - q`. They are Field-shaped objects
that share a dimensional signature with the corresponding **sources**
(both are rate densities, :math:`1/(\mathrm{cm^3 \cdot s \cdot
[sr \cdot] eV})`) yet are a DISTINCT physical kind — a residual is the
imbalance of an equation, not an external drive. Per the Field ABC's
Layer-1 class-identity gate, cross-class arithmetic between a residual
and a source (or a flux) is forbidden *even when units match*; the
balance that combines them goes through a NAMED composition (each
residual leaf's ``from_balance`` factory, B.5 / Issue #201), never a bare
cross-class dunder.

Role grid (field vocabulary, issues #205 / #201)
================================================

This package supplies the **Residual** column of the role grid — three
orthogonal axes ("Boundary" is a LOCUS qualifier, never a fourth
family; #290 P2.5)::

    locus {Bulk, Boundary} × family {Angular, Scalar, Moment}
                           × role {Flux, SourceSink, Residual, Displacement}

       flux     →  orpheus.transport.fields         (…Flux)
       source   →  orpheus.transport.source_sinks   (…SourceSink)
       residual →  orpheus.transport.residuals      (…Residual)

    bulk leaves here: AngularResidual, ScalarResidual;
    boundary leaf: AngularBoundaryResidual;
    starting-direction leaf: RadialCharacteristicResidual (#282 route (a) —
    the ψ½ block of the augmented residual, the third locus's fourth
    role member).

The boundary residual :class:`AngularBoundaryResidual` IS the already-computed
matvec face-defect (the residual of ``γ₋ψ = R·G·γ₊ψ + q``), today
mistyped as :class:`AngularBoundaryFlux`; this package mints the type (B.3) so
the **B.5** operator-output carve can retype the emission sites. See
:mod:`orpheus.transport.residuals.angular_boundary_residual`.

Migration status (field-role-typing plan, step B.3): houses
:class:`AngularResidual`, :class:`ScalarResidual`, and
:class:`AngularBoundaryResidual`.
"""

from __future__ import annotations

from orpheus.transport.residuals.angular_residual import AngularResidual
from orpheus.transport.residuals.scalar_residual import ScalarResidual
from orpheus.transport.residuals.angular_boundary_residual import AngularBoundaryResidual
from orpheus.transport.residuals.radial_characteristic_residual import (
    RadialCharacteristicResidual,
)

__all__ = [
    "AngularResidual",
    "ScalarResidual",
    "AngularBoundaryResidual",
    "RadialCharacteristicResidual",
]
