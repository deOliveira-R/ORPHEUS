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

    #: Wave-7 sweep-cycle signal (§15A.2). A periodic face couples
    #: opposite ends of the spatial domain, creating a dependency
    #: cycle in the SN sweep DAG that requires Krylov closure rather
    #: than a single sweep.
    creates_sweep_cycle: ClassVar[bool] = True

    def apply(
        self,
        psi_out: np.ndarray,
        quadrature: "AngularQuadrature | None" = None,
    ) -> np.ndarray:
        r"""Identity pass-through at the per-face apply level.

        The angular structure is identity (periodic wrap is a
        *spatial* pushforward; the angular components of
        :math:`\psi` are unchanged). The ``quadrature`` argument is
        accepted for backward-compat with the legacy 2-arg form but
        is unused; the spatial-pushforward semantics are handled by
        the sweep orchestrator that couples opposite faces (see the
        class docstring).

        Issue #176 / C176.3: ``quadrature`` is now optional and
        ignored. Body inlined as ``psi_out.copy()`` (matches the
        realizer's :class:`PeriodicWrapOperator.apply`).
        """
        return psi_out.copy()
