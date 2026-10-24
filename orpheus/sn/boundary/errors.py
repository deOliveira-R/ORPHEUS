r"""The SN realizer's own certifications (``errors.py`` of ``orpheus/sn/boundary``).

One guard of the SN realizer checks a support that the type system does not replace, so it stays an SN-side assertion, rooted
in transport's :class:`~orpheus.transport.boundary_conditions.errors.BoundaryError` (the refused object travels as
``.declaration``, not as a string ``.law``):

* :class:`BoundarySourceNotOnIncomingTraceError` (**ERR-047**). The affine trace map :math:`\gamma_-\psi = R\,\gamma_+\psi + q`
  has :math:`q \in \Gamma_-`. A sampler is handed :math:`\Gamma_-(f)` itself and returns exactly its shape, so the support
  holds BY CONSTRUCTION; what the realizer must still know is that the space it samples on IS the inflow tier of the face.
  A source sampled on any other tier would write into slots the sweep silently discards (the total inflow SHORT of the user's
  intent), so ``SNSource.sample`` refuses a space that is not an inflow half-trace.

ERR-041 (the zero map fills only the incoming trace) had a realizer guard, ``VacuumAppliedToOutgoingTraceError``, that
cross-checked a hand-supplied index set against the face name; the typed operators read their half-spaces off the trace ladder,
the index sets retired (WF2.4), and the invariant is witnessed on the ladder
(``tests/gates/numerics/test_angular_trace_space.py::TestTheInflowHalfSpaceIsTheIncomingOrdinates``). The strict inflow/outflow
classification (ERR-040) has no production raise site, so no class is minted for it here.

Layer: ``sn`` (L3). Imports transport's error root only.
"""

from __future__ import annotations

from orpheus.transport.boundary_conditions.errors import BoundaryError

__all__ = ["BoundarySourceNotOnIncomingTraceError"]


class BoundarySourceNotOnIncomingTraceError(BoundaryError):
    r"""ERR-047: a boundary source was sampled on a space that is not the inflow half-trace :math:`\Gamma_-` of its face."""
