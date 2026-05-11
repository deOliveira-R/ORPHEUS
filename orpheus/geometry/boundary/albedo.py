r"""Pure albedo boundary condition.

See :class:`AlbedoBoundary` for the algebraic definition. The legacy
``AlbedoBoundaryOperator`` name is re-exported as a deprecated alias
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


__all__ = ["AlbedoBoundary"]


@dataclass(frozen=True)
class AlbedoBoundary(BoundaryTraceLaw, key="albedo"):
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

    Wave-7 rename note
    ------------------
    Previously named ``AlbedoBoundaryOperator``. The legacy name is
    preserved as a deprecated alias.

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
