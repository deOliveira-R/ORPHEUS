r"""Pure albedo boundary condition.

See :class:`AlbedoBoundaryOperator` for the algebraic definition.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, ClassVar

import numpy as np

from orpheus.numerics.operator import CAP_APPLY

from ._base import BoundaryOperator

if TYPE_CHECKING:
    from orpheus.sn.quadrature import AngularQuadrature


__all__ = ["AlbedoBoundaryOperator"]


@dataclass(frozen=True)
class AlbedoBoundaryOperator(BoundaryOperator, key="albedo"):
    r"""Pure albedo boundary: scalar multiple of the outgoing flux.

    Tensor decomposition :math:`(I, \alpha)` where :math:`I` is the
    angular identity and :math:`\alpha \in [0, 1]` is the albedo:

    .. math::

        \psi_{\text{in}}(\Omega) = \alpha \, \psi_{\text{out}}(\Omega).

    No angular redistribution. Useful as a *building block* for
    :class:`~orpheus.geometry.boundary.mixed.MixedBoundaryOperator`
    (where albedo and specular shares are independent parameters),
    and as a stand-alone primitive when the boundary is a pure
    attenuator with no angular structure.

    Parameters
    ----------
    albedo : float
        Albedo coefficient. ``0`` is vacuum, ``1`` is perfect
        same-direction return.
    """

    capabilities: ClassVar[frozenset[str]] = frozenset({CAP_APPLY})

    albedo: float = 0.0

    def apply(
        self,
        psi_out: np.ndarray,
        quadrature: "AngularQuadrature",
    ) -> np.ndarray:
        return self.albedo * psi_out
