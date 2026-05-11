r"""Prescribed-inflow source for the affine boundary law.

The §16A.1 affine form

.. math::

    \gamma_- \psi \;=\; R\,G\,\gamma_+ \psi \;+\; q

has three terms; this module ships the :math:`q` term -- a
vector-valued quantity on the incoming trace :math:`\Gamma_-`. For
the homogeneous case (vacuum / albedo / reflective / white /
periodic without external injection), :class:`NoSource` is the
sentinel. For a constant scalar inflow, use
:class:`ConstantInflowSource`.

The :class:`BoundarySource` Protocol is :class:`runtime_checkable`
so consumers can ``isinstance(s, BoundarySource)`` to validate the
source contract without subclassing. Custom sources (e.g. a
spatially-varying beam injection) only need to implement
:meth:`evaluate` with the documented signature.

References
----------

* Grand Report v3 §16A.1 (affine boundary form).
* ``.claude/plans/transient-giggling-cake.md`` -- Wave 3 brief.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, Protocol, runtime_checkable

import numpy as np

if TYPE_CHECKING:
    from orpheus.numerics.trace_space import InflowTraceSpace


__all__ = [
    "BoundarySource",
    "NoSource",
    "ConstantInflowSource",
]


@runtime_checkable
class BoundarySource(Protocol):
    r"""Prescribed inflow source :math:`q \in \Gamma_-`.

    The implementor must evaluate the source on a given inflow
    trace and return an array shaped like
    :attr:`InflowTraceSpace.shape` (``(n_ordinates, ng)`` for the
    standard SN trace space, but the contract does not pin
    higher-dimensional shapes -- the only invariant is that the
    returned shape matches the trace's ``shape`` attribute).

    Implementations MUST return ``np.ndarray`` with ``dtype=float``
    so downstream :meth:`apply` paths can add the source without
    casting.
    """

    def evaluate(self, inflow_trace: "InflowTraceSpace") -> np.ndarray:
        """Return source values on the inflow trace."""
        ...


@dataclass(frozen=True)
class NoSource:
    r"""Sentinel for the homogeneous case :math:`q = 0`.

    Default for all rank-1 BCs that don't inject external flux
    (vacuum, albedo, reflective, white, periodic). The evaluate
    method returns a fresh zero array of the trace's shape;
    callers may either skip the source addition entirely (an
    optimization) or add the zeros (uniform code path).
    """

    def evaluate(self, inflow_trace: "InflowTraceSpace") -> np.ndarray:
        return np.zeros(inflow_trace.shape, dtype=float)


@dataclass(frozen=True)
class ConstantInflowSource:
    r"""Constant scalar inflow :math:`q(\Omega) = \text{value}` on every inflow ordinate.

    Useful for fixed-source problems where a uniform incident flux
    is imposed at one face. The output array has the inflow
    trace's full shape; downstream consumers are responsible for
    masking out the outflow side (the source LIVES on
    :math:`\Gamma_-` only -- see
    :class:`~orpheus.geometry.boundary._errors.BoundarySourceNotOnIncomingTraceError`).

    Parameters
    ----------
    value : float
        The scalar inflow level applied to every entry of the
        trace.
    """

    value: float

    def evaluate(self, inflow_trace: "InflowTraceSpace") -> np.ndarray:
        return np.full(inflow_trace.shape, float(self.value), dtype=float)
