r"""Prescribed inflow boundary law.

The §16A.1 affine form

.. math::

    \gamma_- \psi \;=\; R\,G\,\gamma_+ \psi \;+\; q

can carry a nonzero source :math:`q` on the inflow trace.
:class:`PrescribedInflow` implements the rank-0 case where
:math:`R = G = 0` and only the source matters: :math:`\gamma_- \psi
= q`. Used for fixed-source problems where an external incident
flux is imposed at one or more faces (beam injection, fixed-current
boundary, time-of-flight pulse source).

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

The realizer dispatch maps an instance to
:class:`~orpheus.sn.boundary.angular.IncomingSourceOperator`, whose
:meth:`apply` ignores the outgoing flux and returns the source
evaluation.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from ._base import BoundaryTraceLaw
from ._source import InflowSourceSpec, NoSource


__all__ = ["PrescribedInflow"]


@dataclass(frozen=True)
class PrescribedInflow(BoundaryTraceLaw, key="prescribed_inflow"):
    r"""Prescribed inflow source on :math:`\Gamma_-`.

    The :attr:`source` is evaluated against an inflow trace to
    produce the incoming flux; the outgoing trace is irrelevant. The
    rank-0 case :math:`R = G = 0` of the affine boundary law.

    This is a **pure descriptor** (Issue #186 / B3 + β2) — it carries
    no ``apply`` method. The SN realisation is
    :class:`~orpheus.sn.boundary.angular.IncomingSourceOperator(source)`.
    Realise via:

    .. code-block:: python

        from orpheus.sn.boundary.realizer import SNBoundaryRealizer
        from orpheus.sn.mesh.method_space import SNMethodSpace
        law = PrescribedInflow(source=ConstantInflowSource(2.5))
        op = SNBoundaryRealizer().realize(law, SNMethodSpace.minimal(quad))
        psi_in = op.apply(psi_out)   # ignores psi_out; returns source

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

