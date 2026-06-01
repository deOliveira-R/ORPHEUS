r"""L2 typed transport sources — every class inherits from
:class:`~orpheus.numerics.field.Field`.

Sources are Field-shaped objects with DIFFERENT physical units than
fluxes. Per Issue #207's architectural principle and Depth B plan §3.6,
cross-class arithmetic between sources and fluxes is forbidden by the
Field ABC's Layer 1 class-identity gate — even when shape and units
might align. Cross-class same-units operations (e.g., combining an
isotropic source with a per-ordinate source) go through NAMED
composition factories, not cross-class dunders.

Migration status (Depth B step D-F): currently houses
:class:`ScalarSource` and :class:`AngularSource`. Future steps
will add :class:`AngularSource` (Issue #201).
"""

from __future__ import annotations

from orpheus.transport.sources.scalar_source import ScalarSource
from orpheus.transport.sources.angular_source import AngularSource

__all__ = ["ScalarSource", "AngularSource"]
