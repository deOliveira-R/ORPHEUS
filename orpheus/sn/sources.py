r"""DEPRECATED re-export shim — Sources moved to
:mod:`orpheus.transport.source_sinks`.

This shim exists for ONE merge cycle to give in-flight callers a
window to update their imports (per ``feedback_aggressive_retirement``).
It will be deleted in Depth B step D-K.

New code MUST import from the canonical location:

.. code-block:: python

    from orpheus.transport.source_sinks import ScalarSourceSink, AngularSourceSink

Migration rationale: source fields are method-agnostic transport
concepts, so the canonical home is L2 transport, not L3 SN. See plan
§3.1 (layer assignments) and §6 step D-F. Per Issue #207's named-
composition principle, the cross-class ``ScalarSourceSink +
AngularSourceSink`` dunder is RETIRED in the L2 forms; callers
combining iso + per-ordinate sources go through
:meth:`AngularSourceSink.from_isotropic` followed by within-class
addition.
"""

from __future__ import annotations

import warnings as _warnings

from orpheus.transport.source_sinks import (  # noqa: F401
    ScalarSourceSink,
    AngularSourceSink,
)

_warnings.warn(
    "orpheus.sn.sources is deprecated; import ScalarSourceSink and "
    "AngularSourceSink from orpheus.transport.source_sinks instead "
    "(shim retires in Depth B step D-K). The cross-class "
    "ScalarSourceSink + AngularSourceSink dunder is RETIRED in the L2 "
    "forms; use AngularSourceSink.from_isotropic + within-class add.",
    DeprecationWarning,
    stacklevel=2,
)

__all__ = ["ScalarSourceSink", "AngularSourceSink"]
