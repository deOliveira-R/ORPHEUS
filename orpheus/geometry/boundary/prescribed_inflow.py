r"""Prescribed inflow boundary law.

The §16A.1 affine form

.. math::

    \gamma_- \psi \;=\; R\,G\,\gamma_+ \psi \;+\; q

can carry a nonzero source :math:`q` on the inflow trace.
:class:`PrescribedInflow` implements the rank-0 case where
:math:`R = 0` and only the source matters: :math:`\gamma_- \psi
= q`. Used for fixed-source problems where an external incident
flux is imposed at one or more faces (beam injection, fixed-current
boundary, time-of-flight pulse source).

(Campaign phase **B3.0** corrected the older spelling
":math:`R = G = 0`" — it is the RESPONSE that vanishes. The zero map
is not a bijection and so cannot be a geometry map at all, so
:math:`G` is the identity deck element; writing both as zero spelled
one vanishing twice, once in the wrong tier. See
:attr:`PrescribedInflow.geometry_map`.)

Wave-7 architectural slot
-------------------------

The legacy ``BoundaryOperator`` hierarchy (the pre-refactor BC ABC,
retired in Wave O step O.4a.1) treated source-bearing boundaries by
special-casing them in the solver layer; the affine-form refactor
(Wave 3+) lifted the source into a first-class
:meth:`BoundaryTraceLaw.source` slot defaulting to
:class:`~orpheus.geometry.boundary._source.NoSource`. The
prescribed-inflow law is the smallest concrete BC that exercises
this slot.

The realizer dispatch maps an instance to the **zero morphism**
:math:`\Gamma_+ \to \Gamma_-` — the law is affine,
:math:`\gamma_-\psi = L\,\gamma_+\psi + q`, the realizer tier realizes
:math:`L`, and here :math:`L = 0`. The :attr:`source` reaches a solve
through the boundary-source channel instead
(:meth:`~orpheus.transport.source_sinks.AngularBoundarySourceSink.from_mesh_laws`;
see :ref:`bc-affine-source-channel`).

Until **P3** (2026-08-05) the dispatch produced an
``IncomingSourceOperator`` whose ``apply`` ignored the outgoing flux and
returned :math:`q` — an affine map declared as a linear operator, which
reached the ``B`` block and delivered the source there as well. Campaign
phase **B3.4a** had narrowed that operator's codomain to
:math:`\Gamma_-`, retiring the inflow MASK that used to zero the
non-inflow rows of a full-face emission; that is why the affine form's
:math:`q \in \Gamma_-` holds by TYPING rather than by an erasure
(ERR-047), a property the source channel inherits.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from ._base import BoundaryTraceLaw
from ._factors import SelfPairedDeck, ScalarResponse
from ._source import InflowSourceSpec, NoSource


__all__ = ["PrescribedInflow"]


@dataclass(frozen=True)
class PrescribedInflow(BoundaryTraceLaw, key="prescribed_inflow"):
    r"""Prescribed inflow source on :math:`\Gamma_-`.

    The :attr:`source` is evaluated on the inflow half-trace to
    produce the incoming flux; the outgoing trace is irrelevant. The
    rank-0 case :math:`R = 0` of the affine boundary law (:math:`G` is
    the identity deck element — see :attr:`geometry_map`).

    This is a **pure descriptor** (Issue #186 / B3 + β2) — it carries
    no ``apply`` method. The SN realisation is the zero morphism
    :math:`\Gamma_+ \to \Gamma_-` (since **P3**; narrowed to the two
    half-traces at campaign phase **B3.4a**). Realising still needs a
    **face**, and the demand is structural rather than a leftover: a law
    that cannot name its own domain is not realized, it is guessed, and a
    zero map that cannot name its CODOMAIN degenerates to echoing its
    input. ``SNMethodSpace.minimal`` raises.

    .. code-block:: python

        from orpheus.sn.boundary.realizer import SNBoundaryRealizer
        from orpheus.sn.mesh.method_space import SNMethodSpace
        law = PrescribedInflow(source=ConstantInflowSource(2.5))
        space = SNMethodSpace.for_face(
            quadrature=quad, face="xmax", trace=trace,
        )
        op = SNBoundaryRealizer().realize(law, space)
        psi_in = op.apply(psi_out)   # ignores psi_out; returns q on Γ₋

    Only ``psi_out``'s TRAILING axes are read, to size the group /
    spatial block; its leading axis is ignored entirely, so unlike the
    white arm this operator does not validate that it was handed a
    :math:`\Gamma_+`-shaped argument. Being narrowed and validating
    one's own domain are separate properties.

    Parameters
    ----------
    source : InflowSourceSpec
        Source generator. The default :class:`NoSource` makes this
        BC equivalent to vacuum (:math:`\gamma_- = 0`); use
        :class:`~orpheus.geometry.boundary._source.ConstantInflowSource`
        or a custom
        :class:`~orpheus.geometry.boundary._source.InflowSourceSpec`
        implementation for nonzero incident flux.

    Notes
    -----
    A prescribed inflow does NOT wire outgoing flux back into the
    domain, so it contributes no trace back-edge and cannot take part
    in a sweep cycle. (Whether a *configuration* is acyclic is decided
    by the strongly-connected-component decomposition in
    :mod:`orpheus.derivations.discrete.sn.sweep_acyclicity`, never by a
    per-law flag.)

    Implementation note: field-vs-property collision
    ------------------------------------------------
    :class:`BoundaryTraceLaw` declares ``source`` as a ``@property``
    returning :class:`NoSource()`. Declaring ``source: InflowSourceSpec``
    as a dataclass field on the subclass would collide with the
    inherited descriptor (no setter ⇒ dataclass ``__init__`` raises
    ``AttributeError``). The field is therefore named ``_source``
    and the inherited ``source`` property is overridden to read
    ``self._source``. Callers continue to access ``bc.source``
    uniformly across all concrete BCs.
    """

    _source: InflowSourceSpec = field(default_factory=NoSource)

    # ── The affine form's two factors (B1) ──────────────────────────────
    @property
    def geometry_map(self) -> "SelfPairedDeck":
        r""":math:`G = \mathrm{id}` — the inflow comes entirely from :math:`q`.

        Prescribed inflow is the rank-0 AFFINE case: unlike vacuum it has a
        non-zero source, but like vacuum no outgoing flux is routed back. That
        routing is severed by the **response** (:math:`R = 0`), not by the deck
        transformation — the law identifies no geometry, so :math:`G` is the
        identity element.

        **Corrected in B3.** This property previously returned ``NullMap``
        (:math:`G = 0`), and the then-realized ``IncomingSourceOperator``
        described itself as *"the rank-0 case where R = G = 0"* — the same
        conflation, spelled one vanishing twice: the zero map is not a
        bijection and so cannot be a geometry map, and :math:`R = 0` alone
        already severs the routing. (That operator was itself retired at
        **P3**; the ruling here is unchanged by it.)
        """
        return SelfPairedDeck.identity()

    @property
    def response_kernel(self) -> "ScalarResponse":
        r""":math:`R = 0` — nothing crosses; :attr:`source` carries the law."""
        return ScalarResponse(0.0)

    def __init__(self, source: InflowSourceSpec | None = None) -> None:
        # Accept ``source=`` as the public constructor kwarg; route it
        # to the underlying ``_source`` field (named to avoid the
        # collision with the inherited ``BoundaryTraceLaw.source``
        # property). ``object.__setattr__`` bypasses the frozen
        # guard during construction.
        if source is None:
            source = NoSource()
        object.__setattr__(self, "_source", source)

    @property
    def source(self) -> InflowSourceSpec:
        return self._source

