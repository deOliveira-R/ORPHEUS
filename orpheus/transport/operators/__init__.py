r"""Cross-method transport operators — the §5.6 operator algebra.

The method-agnostic :class:`~orpheus.numerics.operator.LinearOperator`\ s that
map a flux carrier to a source/sink — collision :math:`M[\sigma_t]`
(:class:`~orpheus.transport.operators.MultiplicationOperator`),
fission, scattering — plus the
:class:`~orpheus.transport.operators.IntegralKernelOperator`
Protocol they satisfy.

**This package IS the public surface** — import the algebra's members from
here, not from the submodule that happens to define one
(``from orpheus.transport.operators import ScatteringOperator``). The split
across ``transfer`` / ``scattering`` / ``n2n`` / ``fission`` /
``isotropic_scattering`` / ``multiplication_operator`` is file organisation,
not API: which file a member lives in has moved before (#261 relocated the
whole family out of ``orpheus.sn``; #426 step 2 moved the transfer bindings'
shared body into ``transfer``) and may move again, whereas the algebra's
membership is stable.
Re-exported 2026-08-03 — until then the package exported nothing at all while
the theory corpus referenced these names at package level across 16 sites, so
every one of those cross-references was silently dead. Every transport method (SN, CP, MoC, …) collides,
fissions and scatters, so these live at L2 ``transport`` — NOT in a method
package. #261 relocated them out of ``orpheus.sn`` (they are reaction operators,
not SN-specific; only the WDD sweep — ``StreamingOperator`` / ``StreamingCollisionOperator``
/ ``LossRepresentation`` — stays SN-specific). The §5.6 *Functionals*
:class:`~orpheus.transport.reaction_rate_functional.ReactionRateFunctional`
and :class:`~orpheus.transport.reaction_rate_functional.IntegratedReactionRate`
are a distinct abstraction (flux→scalar) and deliberately stay at the
``transport`` top level.
"""

from .fission import FissionOperator
from .integral_kernel_operator import IntegralKernelOperator
from .isotropic_scattering import (
    IsotropicFission,
    IsotropicN2N,
    IsotropicScattering,
    IsotropicTransfer,
)
from .multiplication_operator import MultiplicationOperator
from .n2n import N2NOperator
from .scattering import ScatteringOperator
from .transfer import LegendreMomentTransfer, TransferOperator

__all__ = [
    "FissionOperator",
    "IntegralKernelOperator",
    "IsotropicFission",
    "IsotropicN2N",
    "IsotropicScattering",
    "IsotropicTransfer",
    "LegendreMomentTransfer",
    "MultiplicationOperator",
    "N2NOperator",
    "ScatteringOperator",
    "TransferOperator",
]
