r"""Periodic boundary condition.

See :class:`PeriodicBoundary` for the algebraic definition. The legacy
``PeriodicBoundaryOperator`` name is re-exported as a deprecated alias
from the package ``__init__.py`` (Wave 7 rename per Grand Report v3
vocabulary).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, ClassVar

import numpy as np

from orpheus.numerics.operator import CAP_APPLY

from ._base import BoundaryTraceLaw

if TYPE_CHECKING:
    from orpheus.sn.quadrature import AngularQuadrature


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
    concern not modelled by ``apply`` alone (the
    :class:`BoundaryTraceLaw` consumes one face's outgoing flux at
    a time). The primitive here therefore returns ``psi_out``
    unchanged: the contract is "the incoming side equals the outgoing
    flux you pass in", and the *spatial* coupling is handled by whoever
    instantiates :class:`PeriodicBoundary` and orchestrates the two-face
    plumbing. This is why periodic-BC support in :func:`solve_sn` is a
    downstream wave (it requires sweep changes); this primitive ships
    so that downstream code has the algebraic object to depend on.

    Wave-7 rename note
    ------------------
    Previously named ``PeriodicBoundaryOperator``. The legacy name is
    preserved as a deprecated alias.
    """

    capabilities: ClassVar[frozenset[str]] = frozenset({CAP_APPLY})

    def apply(
        self,
        psi_out: np.ndarray,
        quadrature: "AngularQuadrature",
    ) -> np.ndarray:
        # Wave 7 (C7.3): delegate to the Wave-5 SNBoundaryRealizer.
        # The realizer dispatches on isinstance(self, PeriodicBoundary)
        # and returns the Wave-0 ``PeriodicWrapOperator``, whose
        # apply body is ``psi_out.copy()`` (bit-equivalent to the
        # pre-Wave-7 legacy body).
        #
        # Local imports break the cycle ``orpheus.geometry.boundary``
        # → ``orpheus.sn.boundary_realizer`` → ``orpheus.geometry.boundary``.
        from orpheus.sn.boundary_realizer import (
            SNBoundaryRealizer,
            SNMethodSpace,
        )
        op = SNBoundaryRealizer().realize(
            self, SNMethodSpace.minimal(quadrature),
        )
        return op.apply(psi_out)
