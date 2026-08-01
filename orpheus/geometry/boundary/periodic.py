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

    The two-face coupling — campaign phase B3.4c
    --------------------------------------------

    This is the only law whose **domain is a different face**. Every other law
    is constitutive: a surface responds to what arrives at its own face, so its
    realized operator consumes :math:`\Gamma_+` of the face it is installed on.
    A wrap consumes the PARTNER's, and that is the whole of what makes it a
    quotient rather than a wall.

    The partner is named by the geometry factor
    (:meth:`~orpheus.geometry.boundary.SpatialWrap.domain_face`) and SUPPLIED
    by the composition
    (:attr:`~orpheus.sn.operators.boundary.SNBoundaryOperator._face_domains`),
    which is why the realized operator is a bare
    :class:`~orpheus.numerics.operator.IdentityOperator`: with the right
    half-trace on the input, ordinate :math:`n` of the partner's outflow IS
    ordinate :math:`n` of this face's inflow, because the two faces' outward
    normals are opposite. The realizer CERTIFIES that identification
    (:math:`\Gamma_+(f') \equiv \Gamma_-(f)`) rather than assuming it — the
    user's B3.4c ruling that the quotient reading becomes a guard, not a mesh
    restructure.

    .. warning::

       Two claims that stood here until B3.4c were **false**, and both said the
       same thing in different words: that "realising periodicity at the sweep
       level is a sweep-orchestration concern not modelled by ``apply``", and
       that "the two-face plumbing is handled by whoever instantiates
       :class:`PeriodicBoundary` and orchestrates the sweep". No such caller
       and no such mechanism ever existed — the composition fed every law its
       OWN face's :math:`\Gamma_+`, so periodic silently returned a face's
       outflow as its own inflow, an O(1) wrong answer (MEASURED 98 % relative
       against the partner-face reference). The responsibility was documented
       as delegated and was in fact simply unimplemented, which is the worst
       shape a doc falsehood can take: it reads as a design decision.

    This is a **pure descriptor** (Issue #186 / B3 + β2) — it carries
    no ``apply`` method. Realise via
    :class:`~orpheus.sn.boundary.realizer.SNBoundaryRealizer`.

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

        Since **B3.4c** the realization matches: the factor names the partner
        face and the composition reads it there. The pushforward that #183
        recorded as unbuilt is this channel, and it is the map's whole content
        — the action ON the trace, once the right half-trace arrives, is the
        identity.
        """
        return SpatialWrap(axis=self.axis)

    @property
    def response_kernel(self) -> "ScalarResponse":
        r""":math:`R = 1` — a periodic face is loss-free by construction."""
        return ScalarResponse(1.0)

