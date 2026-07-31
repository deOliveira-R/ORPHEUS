r"""Periodic boundary condition.

See :class:`PeriodicBoundary` for the algebraic definition. This class
was previously named ``PeriodicBoundaryOperator``; the legacy alias
was retired in Wave O step O.4a.1. ``PeriodicBoundary`` (the Grand
Report v3 vocabulary) is the sole live name.
"""

from __future__ import annotations

from dataclasses import dataclass

from ._base import BoundaryTraceLaw
from ._factors import ScalarResponse, SpatialWrap


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
    Realise via :class:`~orpheus.sn.boundary.realizer.SNBoundaryRealizer`.

    Rename history
    --------------
    Previously named ``PeriodicBoundaryOperator``. The legacy alias
    was retired in Wave O step O.4a.1 — ``PeriodicBoundary`` is the
    sole importable name.
    """

    #: The wrap axis (issue #183). Until B1 this law carried **no fields at
    #: all** — so the geometric map it names, a spatial pushforward to the
    #: opposite face, was not expressible on the descriptor. The partner face
    #: is NOT a field: which face is opposite depends on where the law is
    #: installed (configuration), while "wrap along x" is intrinsic. The
    #: realizer derives the partner from the installation face plus this axis.
    axis: str = "x"

    # ── The affine form's two factors (B1) ──────────────────────────────
    @property
    def geometry_map(self) -> "SpatialWrap":
        r""":math:`G` = the spatial pushforward along :attr:`axis`.

        Outflow at one face becomes inflow at the opposite face **at the same
        angle** — which is why a periodic pair closes a sweep cycle from a
        single law, where a lone reflecting face only adds a forward trace
        edge.

        The realized ``PeriodicWrapOperator`` is currently an angular identity
        with the spatial pushforward unbuilt (**#183**); this spec states what
        the law means regardless of that gap.
        """
        return SpatialWrap(axis=self.axis)

    @property
    def response_kernel(self) -> "ScalarResponse":
        r""":math:`R = 1` — a periodic face is loss-free by construction."""
        return ScalarResponse(1.0)

