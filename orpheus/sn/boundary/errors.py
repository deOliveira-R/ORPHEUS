r"""The SN realizer's own certifications (``errors.py`` of ``orpheus/sn/boundary``).

Two guards of the SN realizer check index sets that the type system does not replace, so they stay SN-side assertions, rooted
in transport's :class:`~orpheus.transport.boundary_conditions.errors.BoundaryError` (the refused object travels as
``.declaration``, not as a string ``.law``):

* :class:`VacuumAppliedToOutgoingTraceError` (**ERR-041**). The zero map fills the INFLOW half-trace of a face. The claimed
  ``inflow_indices`` of a hand-built method space and the face NAME are two independently supplied encodings of one orientation
  (the annotation-swap bug class desynchronizes them: indices built from the OUTGOING set of the face), so the realizer
  cross-checks the claimed set against the signed projection :math:`\Omega\cdot n_f` the face name alone implies, through the
  single :func:`~orpheus.numerics.spaces.angular_trace_space.build_omega_dot_n` primitive. Zeroing outflow slots would discard the
  flux the sweep just computed, while every shape-level test stays green (``apply`` returns a same-shaped array regardless of
  orientation). On the canonical ``SNProblem`` path both encodings derive from the trace space, so the guard is green by
  construction; it bites on hand-built method spaces, which is exactly where the annotation swap lives.
* :class:`BoundarySourceNotOnIncomingTraceError` (**ERR-047**). The affine trace map :math:`\gamma_-\psi = R\,\gamma_+\psi + q`
  has :math:`q \in \Gamma_-`. A sampler is handed :math:`\Gamma_-(f)` itself and returns exactly its shape, so the support
  holds BY CONSTRUCTION; what the realizer must still know is that the space it samples on IS the inflow tier of the face.
  A source sampled on any other tier would write into slots the sweep silently discards (the total inflow SHORT of the user's
  intent), so ``SNSource.sample`` refuses a space that is not an inflow half-trace.

The third SN-side certification the plan names, the strict inflow/outflow classification (ERR-040), has no production raise site
today (``[M]`` the old ``assert_inflow_outflow_classification`` was reached only by the retiring law base), so no class is minted
for it here.

Layer: ``sn`` (L3). Imports transport's error root only.
"""

from __future__ import annotations

from orpheus.transport.boundary_conditions.errors import BoundaryError

__all__ = ["BoundarySourceNotOnIncomingTraceError", "VacuumAppliedToOutgoingTraceError"]


class VacuumAppliedToOutgoingTraceError(BoundaryError):
    r"""ERR-041: the zero map was built over claimed inflow ordinates that are OUTGOING at the face (:math:`\Omega\cdot n_f > 0`)."""


class BoundarySourceNotOnIncomingTraceError(BoundaryError):
    r"""ERR-047: a boundary source was sampled on a space that is not the inflow half-trace :math:`\Gamma_-` of its face."""
