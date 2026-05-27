r"""DEPRECATED re-export shim — :class:`ScalarFlux` moved to
:mod:`orpheus.transport.fields.scalar_flux`.

This shim exists for ONE merge cycle to give in-flight callers a
window to update their imports (per ``feedback_aggressive_retirement``).
It will be deleted in Depth B step D-K.

New code MUST import from the canonical location:

.. code-block:: python

    from orpheus.transport.fields.scalar_flux import ScalarFlux

Migration rationale: :class:`ScalarFlux` is a method-agnostic
transport concept (a function on the discretized ``(group ×
spatial)`` phase space), so the canonical home is L2 transport,
not L3 SN. See plan §3.1 (layer assignments) and §6 step D-D
(migration step). The dunder algebra is now inherited from the
L1 :class:`~orpheus.numerics.field.Field` ABC.
"""

from __future__ import annotations

import warnings as _warnings

from orpheus.transport.fields.scalar_flux import ScalarFlux  # noqa: F401

_warnings.warn(
    "orpheus.sn.scalar_flux is deprecated; import "
    "ScalarFlux from orpheus.transport.fields.scalar_flux instead "
    "(shim retires in Depth B step D-K).",
    DeprecationWarning,
    stacklevel=2,
)

__all__ = ["ScalarFlux"]
