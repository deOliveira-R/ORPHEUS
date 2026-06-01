r"""L2 typed transport sources — every class inherits from
:class:`~orpheus.numerics.field.Field`.

Sources are Field-shaped objects with DIFFERENT physical units than
fluxes. Per Issue #207's architectural principle and Depth B plan §3.6,
cross-class arithmetic between sources and fluxes is forbidden by the
Field ABC's Layer 1 class-identity gate — even when shape and units
might align. Cross-class same-units operations (e.g., combining an
isotropic source with a per-ordinate source) go through NAMED
composition factories, not cross-class dunders.

Role grid (field vocabulary, issues #205 / #201): this package supplies
the **source** column —

* bulk: :class:`ScalarSource`, :class:`AngularSource`
  (``orpheus.transport.sources``);
* boundary: :class:`BoundarySource`
  (the prescribed-inflow :math:`q` of the affine boundary law,
  materialised as a stored field on the trace — distinct from the
  geometry-layer :class:`InflowSourceSpec` *generator* it was renamed
  away from; see :mod:`orpheus.transport.sources.boundary_source`).
"""

from __future__ import annotations

from orpheus.transport.sources.scalar_source import ScalarSource
from orpheus.transport.sources.angular_source import AngularSource
from orpheus.transport.sources.boundary_source import BoundarySource

__all__ = ["ScalarSource", "AngularSource", "BoundarySource"]
