"""Test-side oracle of the live projector of a :class:`~orpheus.numerics.space.FunctionSpace` (moved from production, dead-B row 74)."""

from __future__ import annotations

from typing import Any

from orpheus.numerics.space import FunctionSpace


def project_live(function_space: FunctionSpace, x: Any) -> Any:
    r"""The G-orthogonal projector onto the live part, :math:`P = G^{+}G`, of ``function_space`` applied to ``x``.

    :math:`P` is the orthogonal projector onto :math:`\operatorname{range} G = (\ker G)^{\perp}`; it is the identity on a
    non-degenerate space (to rounding: `[M]` at most 2.7e-16 relative, the two metric faces applied in turn) and zeroes the
    kernel coordinates of a degenerate one. Two elements are the same element of the quotient :math:`V / \ker G` iff their
    live projections agree, which is how an adjoint law that holds modulo the kernel is stated and gated (ruling 74; e.g.
    ``A.H.H`` and ``A`` agree after the codomain's projection). One definition for every space: the composites' block-wise
    metric verbs make it block-wise.
    """
    return function_space.apply_inverse_metric(function_space.apply_metric(x))
