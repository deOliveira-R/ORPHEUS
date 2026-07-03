r"""L2 typed transport source/sink fields — every class inherits from
:class:`~orpheus.numerics.field.Field`.

Why ``SourceSink`` (not ``Source``)
===================================

A **source** (production: external :math:`q_{\rm ext}`, in-scattering
:math:`S\psi`, fission :math:`F\psi`) and a **sink** (loss: an operator
output such as the streaming+collision rate :math:`L\psi`) are *the same
physical quantity with opposite sign* — a signed rate-density term in the
balance :math:`L\psi + C\psi - S\psi - F\psi = q`. They share one
dimensional signature (:data:`~orpheus.numerics.units.ANGULAR_RATE_UNITS`
/ :data:`~orpheus.numerics.units.SCALAR_RATE_UNITS` in the bulk; flux
units on the trace) and one storage shape, and they add freely. So the
role is one leaf, and the name says so: ``SourceSink``. (The original
``Source`` name read as production-only and hid that operator ``.apply``
outputs — sinks — land in the *same* type; the rename, B.5, removes that
confusion at the signature level.)

These are Field-shaped objects with DIFFERENT physical units than
fluxes. Per Issue #207's architectural principle and Depth B plan §3.6,
cross-class arithmetic between a source/sink and a flux is forbidden by
the Field ABC's Layer 1 class-identity gate — even when shape and units
might align. Cross-class same-units operations (e.g. combining an
isotropic source/sink with a per-ordinate one) go through the canonical
subspace-containment injection or a NAMED composition factory, not an
ad-hoc cross-class dunder.

Role grid (field vocabulary, issues #205 / #201): this package supplies
the **source/sink** column —

* bulk: :class:`ScalarSourceSink`, :class:`AngularSourceSink`,
  :class:`HarmonicMomentSourceSink` (the moment-space source/sink — the
  output of :math:`\Lambda` in the anisotropic in-scatter
  :math:`R\,\Lambda\,M`, completing the moment row of the
  ``(angular ⊗ moment) × (flux ⊗ source)`` carrier grid);
* boundary: :class:`AngularBoundarySourceSink`
  (the prescribed-inflow :math:`q` of the affine boundary law,
  materialised as a stored field on the trace — distinct from the
  geometry-layer :class:`InflowSourceSpec` *generator* it was renamed
  away from; see :mod:`orpheus.transport.source_sinks.angular_boundary_source_sink`).
"""

from __future__ import annotations

from orpheus.transport.source_sinks.scalar_source_sink import ScalarSourceSink
from orpheus.transport.source_sinks.angular_source_sink import AngularSourceSink
from orpheus.transport.source_sinks.angular_boundary_source_sink import AngularBoundarySourceSink
from orpheus.transport.source_sinks.harmonic_moment_source_sink import (
    HarmonicMomentSourceSink,
)

__all__ = [
    "ScalarSourceSink",
    "AngularSourceSink",
    "AngularBoundarySourceSink",
    "HarmonicMomentSourceSink",
]
