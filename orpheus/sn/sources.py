r"""DEPRECATED re-export shim — Sources moved to
:mod:`orpheus.transport.sources`.

This shim exists for ONE merge cycle to give in-flight callers a
window to update their imports (per ``feedback_aggressive_retirement``).
It will be deleted in Depth B step D-K.

New code MUST import from the canonical location:

.. code-block:: python

    from orpheus.transport.sources import IsotropicSource, PerOrdinateSource

Migration rationale: source fields are method-agnostic transport
concepts, so the canonical home is L2 transport, not L3 SN. See plan
§3.1 (layer assignments) and §6 step D-F. Per Issue #207's named-
composition principle, the cross-class ``IsotropicSource +
PerOrdinateSource`` dunder is RETIRED in the L2 forms; callers
combining iso + per-ordinate sources go through
:meth:`PerOrdinateSource.from_isotropic` followed by within-class
addition.
"""

from __future__ import annotations

import warnings as _warnings

from orpheus.transport.sources import (  # noqa: F401
    IsotropicSource,
    PerOrdinateSource,
)

_warnings.warn(
    "orpheus.sn.sources is deprecated; import IsotropicSource and "
    "PerOrdinateSource from orpheus.transport.sources instead "
    "(shim retires in Depth B step D-K). The cross-class "
    "IsotropicSource + PerOrdinateSource dunder is RETIRED in the L2 "
    "forms; use PerOrdinateSource.from_isotropic + within-class add.",
    DeprecationWarning,
    stacklevel=2,
)

__all__ = ["IsotropicSource", "PerOrdinateSource"]
