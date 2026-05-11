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

The legacy ``BoundaryOperator`` hierarchy treated source-bearing
boundaries by special-casing them in the solver layer; the new
affine-form refactor (Wave 3+) lifts the source into a first-class
:meth:`BoundaryTraceLaw.source` slot defaulting to
:class:`~orpheus.geometry.boundary._source.NoSource`. The
prescribed-inflow law is the smallest concrete BC that exercises
this slot.

The realizer dispatch (Wave 7 / C7.5) maps an instance to
:class:`~orpheus.sn.angular_operator.IncomingSourceOperator`, whose
:meth:`apply` ignores the outgoing flux and returns the source
evaluation.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING, ClassVar

import numpy as np

from orpheus.numerics.operator import CAP_APPLY

from ._base import BoundaryTraceLaw
from ._source import BoundarySource, NoSource

if TYPE_CHECKING:
    from orpheus.sn.quadrature import AngularQuadrature


__all__ = ["PrescribedInflow"]


@dataclass(frozen=True)
class PrescribedInflow(BoundaryTraceLaw, key="prescribed_inflow"):
    r"""Prescribed inflow source on :math:`\Gamma_-`.

    The :attr:`source` is evaluated against an inflow trace to
    produce the incoming flux; the outgoing trace is irrelevant (the
    :meth:`apply` ignores ``psi_out``). The rank-0 case
    :math:`R = G = 0` of the affine boundary law.

    Parameters
    ----------
    source : BoundarySource
        Source generator. The default :class:`NoSource` makes this
        BC equivalent to vacuum (:math:`\gamma_- = 0`); use
        :class:`~orpheus.geometry.boundary._source.ConstantInflowSource`
        or a custom
        :class:`~orpheus.geometry.boundary._source.BoundarySource`
        implementation for nonzero incident flux.

    Notes
    -----
    :pydata:`creates_sweep_cycle` is ``False``: a prescribed inflow
    does NOT wire outgoing flux back into the domain, so the SN
    sweep DAG remains acyclic for this BC. The §15A.2 cycle detector
    treats prescribed-inflow faces as plain Dirichlet faces.

    Source evaluation
    -----------------
    The :meth:`apply` synthesises a probe
    :class:`~orpheus.numerics.trace_space.InflowTraceSpace` carrying
    the shape of the incoming-flux buffer (``psi_out.shape``) and
    delegates to :meth:`BoundarySource.evaluate`. Sources that need
    richer trace metadata (face-tagged inflow injection, per-ordinate
    masks) require Wave 8's full
    :class:`~orpheus.sn.method_space.SNMethodSpace` wiring; the
    Wave-7 ship-state covers the :class:`NoSource` /
    :class:`~orpheus.geometry.boundary._source.ConstantInflowSource`
    cases that only need ``shape``.

    Implementation note: field-vs-property collision
    ------------------------------------------------
    :class:`BoundaryTraceLaw` declares ``source`` as a ``@property``
    returning :class:`NoSource()`. Declaring ``source: BoundarySource``
    as a dataclass field on the subclass would collide with the
    inherited descriptor (no setter ⇒ dataclass ``__init__`` raises
    ``AttributeError``). The field is therefore named ``_source``
    and the inherited ``source`` property is overridden to read
    ``self._source``. Callers continue to access ``bc.source``
    uniformly across all concrete BCs.
    """

    capabilities: ClassVar[frozenset[str]] = frozenset({CAP_APPLY})
    creates_sweep_cycle: ClassVar[bool] = False

    _source: BoundarySource = field(default_factory=NoSource)

    def __init__(self, source: BoundarySource | None = None) -> None:
        # Accept ``source=`` as the public constructor kwarg; route it
        # to the underlying ``_source`` field (named to avoid the
        # collision with the inherited ``BoundaryTraceLaw.source``
        # property). ``object.__setattr__`` bypasses the frozen
        # guard during construction (same pattern as
        # :class:`MixedBoundaryOperator`).
        if source is None:
            source = NoSource()
        object.__setattr__(self, "_source", source)

    @property
    def source(self) -> BoundarySource:
        return self._source

    @property
    def kind(self) -> str:
        return "prescribed_inflow"

    def apply(
        self,
        psi_out: np.ndarray,
        quadrature: "AngularQuadrature",
    ) -> np.ndarray:
        # Wave 7 (C7.5): delegate to the realizer's IncomingSourceOperator.
        # The realizer dispatches on isinstance(self, PrescribedInflow)
        # and returns IncomingSourceOperator(self.source).
        from orpheus.sn.boundary_realizer import (
            SNBoundaryRealizer,
            SNMethodSpace,
        )
        op = SNBoundaryRealizer().realize(
            self, SNMethodSpace.minimal(quadrature),
        )
        return op.apply(psi_out)
