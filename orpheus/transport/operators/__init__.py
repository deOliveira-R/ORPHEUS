r"""Cross-method transport operators — the §5.6 operator algebra.

The method-agnostic :class:`~orpheus.numerics.operator.LinearOperator`\ s that
map a flux carrier to a source/sink — collision :math:`M[\sigma_t]`
(:class:`~orpheus.transport.operators.multiplication_operator.MultiplicationOperator`),
fission, scattering — plus the
:class:`~orpheus.transport.operators.integral_kernel_operator.IntegralKernelOperator`
Protocol they satisfy. Every transport method (SN, CP, MoC, …) collides,
fissions and scatters, so these live at L2 ``transport`` — NOT in a method
package. #261 relocated them out of ``orpheus.sn`` (they are reaction operators,
not SN-specific; only the WDD sweep — ``StreamingOperator`` / ``InvertibleOperator``
/ ``LossRepresentation`` — stays SN-specific). The §5.6 *Functionals*
:class:`~orpheus.transport.reaction_rate_functional.ReactionRateFunctional`
and :class:`~orpheus.transport.reaction_rate_functional.IntegratedReactionRate`
are a distinct abstraction (flux→scalar) and deliberately stay at the
``transport`` top level.
"""
