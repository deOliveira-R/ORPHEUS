r"""DEPRECATED re-export shim — :class:`HarmonicMomentField` moved to
:mod:`orpheus.transport.fields.harmonic_moment_field`.

This shim exists for ONE merge cycle to give in-flight callers a
window to update their imports (per ``feedback_aggressive_retirement``).
It will be deleted in Depth B step D-K.

New code MUST import from the canonical location:

.. code-block:: python

    from orpheus.transport.fields.harmonic_moment_field import HarmonicMomentField

Migration rationale: :class:`HarmonicMomentField` is a method-agnostic
transport concept (a moment field on the spherical-harmonic-tensor
spatial-group phase space), so the canonical home is L2 transport, not
L3 SN. See plan §3.1 (layer assignments) and §6 step D-E (migration
step). The dunder algebra is now inherited from the L1
:class:`~orpheus.numerics.field.Field` ABC; the space is a
:class:`~orpheus.numerics.space.TensorProductSpace` of the form
:math:`\mathrm{SphericalHarmonicSpace}(L) \otimes
\mathrm{CellGroupSpace}` — making this the FIRST production consumer
of D-B's L1 :class:`TensorProductSpace` primitive in a typed field.
"""

from __future__ import annotations

import warnings as _warnings

from orpheus.transport.fields.harmonic_moment_field import (  # noqa: F401
    HarmonicMomentField,
)

_warnings.warn(
    "orpheus.sn.harmonic_moment_field is deprecated; import "
    "HarmonicMomentField from "
    "orpheus.transport.fields.harmonic_moment_field instead "
    "(shim retires in Depth B step D-K).",
    DeprecationWarning,
    stacklevel=2,
)

__all__ = ["HarmonicMomentField"]
