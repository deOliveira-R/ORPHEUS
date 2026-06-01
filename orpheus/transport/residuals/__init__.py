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
balance that combines them goes through a NAMED composition
(``IterationResidual.from_balance``, B.5 / Issue #201), never a bare
cross-class dunder.

Role grid (field vocabulary, issues #205 / #201)
================================================

This package supplies the **residual** column of the bulk role grid::

    {Angular, Scalar} × {Flux, Source, Residual}
       flux    →  orpheus.transport.fields    (AngularFlux, ScalarFlux)
       source  →  orpheus.transport.sources   (AngularSource, ScalarSource)
       residual →  orpheus.transport.residuals (AngularResidual, ScalarResidual)

The boundary residual column (``BoundaryResidual``) is DEFERRED to the
operator-typing / BC-extraction work (Issue #208), where it gains a
real consumer; see ``field_role_typing_view_g.md`` step B.3 notes.

Migration status (field-role-typing plan, step B.3): houses
:class:`AngularResidual` and :class:`ScalarResidual`.
"""

from __future__ import annotations

from orpheus.transport.residuals.angular_residual import AngularResidual
from orpheus.transport.residuals.scalar_residual import ScalarResidual

__all__ = ["AngularResidual", "ScalarResidual"]
