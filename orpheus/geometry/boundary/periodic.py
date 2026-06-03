r"""Periodic boundary condition.

See :class:`PeriodicBoundary` for the algebraic definition. This class
was previously named ``PeriodicBoundaryOperator``; the legacy alias
was retired in Wave O step O.4a.1. ``PeriodicBoundary`` (the Grand
Report v3 vocabulary) is the sole live name.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import ClassVar

from ._base import BoundaryTraceLaw


__all__ = ["PeriodicBoundary"]


@dataclass(frozen=True)
class PeriodicBoundary(BoundaryTraceLaw, key="periodic"):
    r"""Periodic boundary: spatial pushforward to the partner face.

    Tensor decomposition :math:`(G_{\text{wrap}}, 1)` where
    :math:`G_{\text{wrap}}` is the spatial pushforward to the
    partner face of the domain (e.g. left <-> right): the incoming flux
    at this face equals the outgoing flux at the partner face for
    every ordinate, with no angular permutation:

    .. math::

        \psi_{\text{in}}(\Omega, x_{\text{this}})
        = \psi_{\text{out}}(\Omega, x_{\text{partner}}).

    Realising periodicity at the *sweep* level requires coupling the
    two faces' boundary-flux buffers -- which is a sweep-orchestration
    concern not modelled by ``apply`` alone. The SN realisation is
    :class:`~orpheus.numerics.operator.PeriodicWrapOperator` (currently
    an angular identity; the spatial-pushforward extension is tracked
    as a follow-up under ``module:sn``).

    This is a **pure descriptor** (Issue #186 / B3 + β2) — it carries
    no ``apply`` method. The two-face plumbing is handled by whoever
    instantiates :class:`PeriodicBoundary` and orchestrates the sweep.
    Realise via :class:`~orpheus.sn.boundary_realizer.SNBoundaryRealizer`.

    Rename history
    --------------
    Previously named ``PeriodicBoundaryOperator``. The legacy alias
    was retired in Wave O step O.4a.1 — ``PeriodicBoundary`` is the
    sole importable name.
    """

    #: Wave-7 sweep-cycle signal (§15A.2). A periodic face couples
    #: opposite ends of the spatial domain, creating a dependency
    #: cycle in the SN sweep DAG that requires Krylov closure rather
    #: than a single sweep.
    creates_sweep_cycle: ClassVar[bool] = True
