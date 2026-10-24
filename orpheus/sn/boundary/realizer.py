r"""The SN realizer table: ``SN_REALIZERS`` (``realizer.py`` of ``orpheus/sn/boundary``).

**What this module is.** The ONE table through which SN admits and realizes a declared exterior. It replaces
``SNBoundaryRealizer``, a 400-line ``isinstance`` ladder over the string-tagged boundary laws, with a
:class:`~orpheus.transport.boundary_conditions.RealizerTable` keyed by UNIT (a chain of the closed sum
:data:`~orpheus.transport.boundary_conditions.Chain`, the marker ``ZeroNetwork``, the marker ``InflowTerms``) and filled by five
realizer classes that register by class keyword. The table is sealed on the last line of this module, so anyone holding
``SN_REALIZERS`` has executed every registration.

**The functor.** Realization is a FUNCTOR on a wall's response network (:mod:`orpheus.sn.boundary.realized_operators`): the empty
network goes to :class:`~orpheus.sn.boundary.realized_operators.ZeroMap`, ``s * SpecularChain`` to
:class:`~orpheus.sn.boundary.realized_operators.ScaledOrdinatePermutation` (:math:`s\,P_\pi`), ``d * DiffuseChain`` to
:class:`~orpheus.sn.boundary.realized_operators.RankOneReemission` (:math:`d\,(E\circ R_t)`), and a sum of chains to the sum of the
images (:class:`~orpheus.sn.boundary.realized_operators.RealizedSum`, in the canonical chain order). A method never writes a
branch over the shared family's cone point: ``WallCondition`` itself is not a key of the table. The affine part :math:`q` realizes
separately, as an array (:class:`SNSource` has no ``realize``).

* ``SNZero`` (key ``ZeroNetwork``): admitted at every site except the hollow inner face (#511); image :class:`~orpheus.sn.boundary.realized_operators.ZeroMap`.
* ``SNSpecularChain`` (key ``SpecularChain``): ``ordinate_pairing`` certifies the mirror, else ``Refused(ANGULAR_RESOLUTION)``; the hollow inner
  face only at weight 1.0 (the cavity the closure computes, #511); image :class:`~orpheus.sn.boundary.realized_operators.ScaledOrdinatePermutation`.
* ``SNDiffuseChain`` (key ``DiffuseChain``): a radial outer face is refused (#653), the hollow inner face too (#511); image
  :class:`~orpheus.sn.boundary.realized_operators.RankOneReemission`.
* ``SNSource`` (key ``InflowTerms``, a :class:`~orpheus.transport.boundary_conditions.SourceRealizer`): as the diffuse chain, because the radial corner
  action is ruled only for a permuting return or a zero return (#653, #511); the image is an array (``sample``), never an operator.
* ``SNDecks`` (the one :class:`~orpheus.transport.boundary_conditions.DeckRealizer`): a mirror is certified by ``ordinate_pairing``; a translation wrap by
  the equality of its crossing sets; the image is a ``ScaledOrdinatePermutation`` at weight 1.0.

A wall of two chains is the meet of the two chain admissions, and a mixed wall at a radial outer site is refused through its
diffuse chain. A method-only condition of ANOTHER method (diffusion's zero flux) meets the table's unregistered-unit door:
:meth:`~orpheus.transport.boundary_conditions.RealizerTable.admit` returns ``Refused(STATE_CONE)`` because the table's
:class:`~orpheus.transport.boundary_conditions.StateCone` is ``NONNEGATIVE_FLUX`` and the unit's return fraction is negative,
which keeps the refusal's axis without SN importing the diffusion type (SN may not import diffusion: an L3 to L3 edge).

**Site classes.** A :class:`~orpheus.transport.boundary_conditions.FaceSite` reads as one of three geometric classes:
``CART`` (any Cartesian face), ``RAD_OUTER`` (a solid or hollow radial body's outer face) and ``RAD_INNER`` (a hollow radial
body's inner face, a tag with NO trace slot, so it never reaches ``realize``: the shared body meets the table through
``require`` there). The corner action of the radial characteristic is ruled only for a permuting return or a zero return, so the
diffuse chain and every source are refused at ``RAD_OUTER`` AT ADMIT (#653), where they used to construct and refuse at the first
sweep (12 of 12, plan section 4.3 delta 4).

**Two assertions the type system does not replace** (rooted in ``BoundaryError``, :mod:`orpheus.sn.boundary.errors`):
ERR-041 (the claimed inflow set of a hand-built method space against the orientation the face name implies, in ``SNZero.realize``)
and ERR-047 (a source is sampled on the inflow half-trace of its face, in ``SNSource.sample``). Neither is an admission: the
declared orientation they used to catch is unspellable now, and what remains is the realizer's own index sets.

**Retired with the isinstance ladder** (one body per builder, decision 45): ``_narrowed_zero_operator``, ``_deck_kernel``,
``_checked_angular_average`` and ``_specular_crossed_diagnosis`` (their bodies live in
:mod:`~orpheus.sn.boundary.realized_operators`, which owns them), the Gamma-to-Gamma_+ restrictions ``_outflow_restriction`` and
``_partner_outflow_restriction`` (dead: the typed operators read the half-spaces off the trace ladder), and the string
vocabulary their signatures carried (``law_key``, ``crossed_diagnosis``, ``axis: str``, ``outward_sign``).

References
----------
* ``.claude/plans/wave2_boundary_carve.md`` section 2.6 (the table and the typed operators), section 4.2 (the SN columns).
* ``docs/theory/foundations/boundary_conditions.rst`` (``bc-method-realizability``, ``bc-refusal-axes``).
"""

from __future__ import annotations

from typing import assert_never

import numpy as np

from orpheus.geometry.deck import DeckTransformation, FacePairing, MirrorDeck, TranslationDeck, mirror_motion
from orpheus.geometry.faces import AXIS_NAMES, face_normal
from orpheus.geometry.transformation import Permutation
from orpheus.numerics.quadrature.ordinate_pairing import OrdinatePairingError, ordinate_pairing
from orpheus.numerics.spaces.angular_trace_space import TANGENTIAL_EPS, AngularFaceTraceSpace, build_omega_dot_n
from orpheus.sn.boundary.errors import BoundarySourceNotOnIncomingTraceError, VacuumAppliedToOutgoingTraceError
from orpheus.sn.boundary.realized_operators import (
    RankOneReemission,
    RealizedBoundary,
    ScaledOrdinatePermutation,
    ZeroMap,
)
from orpheus.sn.mesh.method_space import SNMethodSpace
from orpheus.transport.boundary_conditions import (
    Admission,
    BoundaryError,
    DeckRealizer,
    DiffuseChain,
    FaceSite,
    InflowTerms,
    Realized,
    RealizerTable,
    RefusalAxis,
    Refused,
    ResponseRealizer,
    SiteClass,
    SourceRealizer,
    SpecularChain,
    StateCone,
    ZeroNetwork,
    sample_terms,
)

__all__ = ["SN_REALIZERS", "SNDecks", "SNDiffuseChain", "SNSource", "SNSpecularChain", "SNZero"]

#: SN's realizer table. ``unregistered_issue=659``: SN registers every shared unit, so its unregistered-unit door meets only
#: another method's own condition (issue #659 tracks the SN and diffusion tables' unregistered units).
SN_REALIZERS: RealizerTable[SNMethodSpace, RealizedBoundary, AngularFaceTraceSpace] = RealizerTable("sn", state_cone=StateCone.NONNEGATIVE_FLUX, unregistered_issue=659)


class _SNSites:
    """What the five realizers share: the readings of a site and of a method space (no module-level helper functions)."""

    #: The hollow radial body's inner face: a tag with no trace and so no ``realize`` (``require`` is its only door).
    _INNER_REFUSAL = "the inner face of a hollow radial body admits only a unit specular return (the cavity its angular closure computes)"
    #: The radial corner action is ruled only for a permuting return or a zero return.
    _CORNER_REFUSAL = "the radial characteristic's corner action is ruled only for a permuting return or a zero return"

    @staticmethod
    def inner_refusal(reason: str | None = None) -> Refused:
        r"""The refusal of every declaration the hollow radial inner face cannot honour, citing #511.

        SCOPE-BOUNDARY[guard] machinery: a radial characteristic that sweeps the inner surface (a trace slot and a corner action on the inner face).
        ruling: wave2_boundary_carve.md section 4.2 (the SN columns), the specification approved by the user 2026-10-23.
        revisit: #511, when the radial axis gains an inner-surface trace.
        """
        return Refused(RefusalAxis.PLUMBING, reason or _SNSites._INNER_REFUSAL, 511)

    @staticmethod
    def corner_refusal() -> Refused:
        return Refused(RefusalAxis.PLUMBING, _SNSites._CORNER_REFUSAL, 653)

    @staticmethod
    def trace_of(space: SNMethodSpace, owner: str):
        """The trace ladder of a method space (the realizers read every half-space off it)."""
        if space is None or space.trace is None:
            raise BoundaryError(
                f"{owner}: realization reads the half-spaces off the trace ladder, and this method space carries none. Construct it via "
                f"SNMethodSpace.for_face(quadrature=..., face=..., trace=...).",
                owner,
            )
        return space.trace

    @staticmethod
    def mirror_pairing(space: SNMethodSpace, face: str, dimension: int) -> Permutation:
        r"""The certified ordinate permutation of the mirror of ``face``'s own plane (``ordinate_pairing`` of ``mirror_motion``)."""
        axis, _ = face_normal(face)
        return ordinate_pairing(space.quadrature, mirror_motion(AXIS_NAMES[axis], dimension))


class SNZero(ResponseRealizer, _SNSites, table=SN_REALIZERS, realizes=ZeroNetwork):
    r"""The image of the empty network: the zero map :math:`\Gamma_+ \to \Gamma_-`.

    Vacuum's whole content is :math:`R = 0`. The same object is the linear image of a pure inflow (the affine :math:`q` travels the
    boundary-source channel) and of any wall of zero weight, so the question "is this face a vacuum" is answered by the REALIZED
    map. Admitted at every site except the hollow radial inner face, where only the unit specular cavity is (#511).
    """

    status = Realized(True)

    def admit(self, unit: ZeroNetwork, weight: float, site: FaceSite, space: SNMethodSpace | None = None) -> Admission:
        return self.inner_refusal() if site.site_class is SiteClass.RAD_INNER else Realized(True)

    def realize(self, unit: ZeroNetwork, weight: float, site: FaceSite, space: SNMethodSpace) -> ZeroMap:
        if space.face is not None and space.inflow_indices is not None:
            # ERR-041: the claimed inflow set and the face NAME are two encodings of one orientation (the annotation-swap bug class
            # desynchronizes them); cross-check against the signed projection the face name alone implies. A method space with no face
            # carries no independent orientation truth, so the guard cannot fire there (the caller owns the claim).
            omega_dot_n = build_omega_dot_n(space.quadrature, (space.face,))[0]
            claimed = np.asarray(space.inflow_indices, dtype=np.intp)
            outgoing = omega_dot_n[claimed] > +TANGENTIAL_EPS
            if np.any(outgoing):
                raise VacuumAppliedToOutgoingTraceError(
                    f"ZeroNetwork at face {space.face!r}: inflow_indices {claimed[outgoing].tolist()} are OUTGOING ordinates (Omega.n > 0) at this "
                    f"face. The zero map fills gamma_-psi on the INCOMING trace only; zeroing outflow slots discards the flux the sweep just "
                    f"computed (ERR-041: swapped inflow/outflow face annotation).",
                    unit,
                )
        return ZeroMap(face=site.face, trace=self.trace_of(space, "SNZero"))


class SNSpecularChain(ResponseRealizer, _SNSites, table=SN_REALIZERS, realizes=SpecularChain):
    r"""The image of ``s * SpecularChain``: :math:`s\,P_\pi`, the scaled ordinate permutation of the face's own mirror plane.

    ``ordinate_pairing`` certifies the pairing (measure-preserving, involutive, inflow to outflow) and RETURNS the permutation the
    image uses, so there is one computation and not a certification beside a second derivation. A rule that is not closed under the
    mirror (an odd-``n_phi`` product rule has no x-mirror closure) is refused on the angular-resolution axis. The hollow radial inner
    face admits the chain only at weight 1.0, the cavity its closure computes (#511).
    """

    status = Realized(True)

    def admit(self, unit: SpecularChain, weight: float, site: FaceSite, space: SNMethodSpace | None = None) -> Admission:
        if site.site_class is SiteClass.RAD_INNER:
            return Realized(True) if weight == 1.0 else self.inner_refusal()
        if space is None:
            return Realized(True)
        try:
            self.mirror_pairing(space, site.face, site.dimension)
        except OrdinatePairingError as error:
            return Refused(RefusalAxis.ANGULAR_RESOLUTION, f"the quadrature does not certify the mirror of {site.face!r}: {error}")
        return Realized(True)

    def realize(self, unit: SpecularChain, weight: float, site: FaceSite, space: SNMethodSpace) -> ScaledOrdinatePermutation:
        trace = self.trace_of(space, "SNSpecularChain")
        return ScaledOrdinatePermutation(
            alpha=weight, permutation=self.mirror_pairing(space, site.face, site.dimension), face=site.face, partner=site.face, trace=trace
        )


class SNDiffuseChain(ResponseRealizer, _SNSites, table=SN_REALIZERS, realizes=DiffuseChain):
    r"""The image of ``d * DiffuseChain``: :math:`d\,(E\circ R_t)`, the factored Lambertian.

    The chain asks no pairing, so every quadrature admits it on a Cartesian face. A radial outer face is refused (#653: the corner
    action is ruled only for a permuting return or a zero return), and so is the hollow inner face (#511).
    """

    status = Realized(True)

    def admit(self, unit: DiffuseChain, weight: float, site: FaceSite, space: SNMethodSpace | None = None) -> Admission:
        match site.site_class:
            case SiteClass.CART:
                return Realized(True)
            case SiteClass.RAD_OUTER:
                return self.corner_refusal()
            case SiteClass.RAD_INNER:
                return self.inner_refusal()
            case _:
                assert_never(site.site_class)

    def realize(self, unit: DiffuseChain, weight: float, site: FaceSite, space: SNMethodSpace) -> RankOneReemission:
        return RankOneReemission(alpha=weight, face=site.face, trace=self.trace_of(space, "SNDiffuseChain"))


class SNSource(SourceRealizer[AngularFaceTraceSpace], _SNSites, table=SN_REALIZERS):
    r"""The source channel: :math:`q` as an ARRAY on :math:`\Gamma_-`, never an operator (no ``realize``).

    Admitted on a Cartesian face. A radial outer face is refused (#653: a source is a free parameter the corner block, a linear
    operator on the ray alone, cannot carry), and so is the hollow inner face (#511). :meth:`sample` runs the ERR-047 placement
    check: the space it samples on must be the inflow half-trace of the face.
    """

    status = Realized(True)

    def admit(self, terms: InflowTerms, site: FaceSite, space: SNMethodSpace | None = None) -> Admission:
        match site.site_class:
            case SiteClass.CART:
                return Realized(True)
            case SiteClass.RAD_OUTER:
                return self.corner_refusal()
            case SiteClass.RAD_INNER:
                return self.inner_refusal()
            case _:
                assert_never(site.site_class)

    def sample(self, terms: InflowTerms, site: FaceSite, space: AngularFaceTraceSpace) -> np.ndarray:
        if getattr(space, "role", None) != "inflow":
            raise BoundarySourceNotOnIncomingTraceError(
                f"the boundary source of face {site.face!r} was sampled on {getattr(space, 'name', space)!r}, which is not an inflow half-trace: the "
                f"affine trace map has q on Gamma_- only, and a source written anywhere else is discarded by the sweep (ERR-047).",
                terms,
            )
        return sample_terms(terms, space)


class SNDecks(DeckRealizer, _SNSites, table=SN_REALIZERS):
    r"""The decks: a mirror plane or a translation wrap, realized as the unit scaled permutation :math:`\alpha = 1`.

    * ``MirrorDeck``: the ordinate permutation of the face's own mirror, certified by ``ordinate_pairing`` (the refusal is on the
      angular-resolution axis). It is the same matrix as a unit specular wall, because a symmetry has no amplitude.
    * ``TranslationDeck``: the identity relabelling between two DISTINCT index sets, :math:`\Gamma_+(f') \to \Gamma_-(f)` with
      ``f'`` the opposite face, EARNED by the crossing check (a translation moves no direction; the outward normals are opposite, so
      :math:`n_f = -n_{f'}` sends outgoing to incoming) rather than assumed. A pair whose crossing sets differ is refused on the
      topological axis: it is not a translation quotient and needs an explicit gluing map (#178).

    The match is exhaustive over :data:`~orpheus.geometry.deck.DeckTransformation`.
    """

    def admit(self, deck: DeckTransformation, pairing: FacePairing, site: FaceSite, space: SNMethodSpace | None = None) -> Admission:
        match deck:
            case MirrorDeck():
                if space is None:
                    return Realized(True)
                try:
                    self.mirror_pairing(space, pairing.face, site.dimension)
                except OrdinatePairingError as error:
                    return Refused(RefusalAxis.ANGULAR_RESOLUTION, f"the quadrature does not certify the mirror of {pairing.face!r}: {error}")
                return Realized(True)
            case TranslationDeck():
                if space is None:
                    return Realized(True)
                trace = self.trace_of(space, "SNDecks")
                try:
                    # The type refuses at construction a pairing that does not cross or that moves measure: admission IS construction.
                    ScaledOrdinatePermutation(
                        alpha=1.0, permutation=Permutation.identity(int(space.quadrature.n_ordinates)),
                        face=pairing.face, partner=pairing.partner, trace=trace,
                    )
                except ValueError as error:
                    return Refused(RefusalAxis.TOPOLOGY, f"the wrap does not carry Gamma_-({pairing.face}) into Gamma_+({pairing.partner}): {error}")
                return Realized(True)
            case _:
                assert_never(deck)

    def realize(self, deck: DeckTransformation, pairing: FacePairing, site: FaceSite, space: SNMethodSpace) -> ScaledOrdinatePermutation:
        trace = self.trace_of(space, "SNDecks")
        match deck:
            case MirrorDeck():
                permutation = self.mirror_pairing(space, pairing.face, site.dimension)
            case TranslationDeck():
                permutation = Permutation.identity(int(space.quadrature.n_ordinates))
            case _:
                assert_never(deck)
        return ScaledOrdinatePermutation(alpha=1.0, permutation=permutation, face=pairing.face, partner=pairing.partner, trace=trace)


SN_REALIZERS.seal()
