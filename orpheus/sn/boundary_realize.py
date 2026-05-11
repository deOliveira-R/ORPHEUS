r"""Tree-walking realisation for Wave-0-composed boundary laws.

The :meth:`~orpheus.sn.boundary_realizer.SNBoundaryRealizer.realize`
method dispatches by ``isinstance(law, BoundaryTraceLaw)`` and returns
a 1-arg :class:`~orpheus.numerics.operator.LinearOperator`. But users
may build a boundary law by composing multiple BCs via Wave-0 algebra:

.. code-block:: python

    law = (
        0.3 * ReflectiveBoundary(axis="x")
        + 0.7 * WhiteBoundary(axis="x", outward_sign=+1)
    )

Here ``law`` is an :class:`~orpheus.numerics.operator.OperatorSum`
containing :class:`~orpheus.numerics.operator.ScaledOperator` wrappers
around the leaf :class:`~orpheus.geometry.boundary.BoundaryTraceLaw`
instances. Passing this directly to
:meth:`SNBoundaryRealizer.realize` raises
:class:`~orpheus.geometry.boundary.BoundaryError` because the realiser
does NOT dispatch on Wave-0 composers — they are not
``BoundaryTraceLaw`` subclasses.

:func:`realize_recursively` walks the expression tree, realising each
leaf ``BoundaryTraceLaw`` against the method space, then re-assembling
the result through the same Wave-0 composers:

.. code-block:: python

    realised = realize_recursively(law, method_space)
    # realised is OperatorSum(ScaledOperator(0.3, realised_spec),
    #                        ScaledOperator(0.7, realised_white))
    # where realised_spec and realised_white are 1-arg Wave-0
    # primitives consumable by the SN sweep / Krylov path.

Wave-11 background
==================

Pre-Wave-11, the rank-N composition was handled by a dedicated
:class:`MixedBoundaryOperator` class whose ``apply`` delegated to the
realiser's ``isinstance(law, MixedBoundaryOperator)`` branch (itself a
loop computing
``sum(coeff * realize(prim, method_space) for coeff, prim in components)``
and reassembling via ``OperatorSum``). Wave 11 deleted that composer:
Wave-0's :class:`~orpheus.numerics.operator.OperatorSum` /
:class:`~orpheus.numerics.operator.ScaledOperator` algebra is sufficient,
and the per-BC class added no value over the algebra dunders inherited
by every :class:`BoundaryTraceLaw` via
:class:`~orpheus.numerics.operator.LinearOperatorMixin`. This module
houses the tree-walker that replaces the deleted recursive realiser
branch.

The walker recognises the Wave-0 composers in
:mod:`orpheus.numerics.operator` (``ScaledOperator``, ``OperatorSum``,
``OperatorProduct``, ``TensorProductOperator``,
``SumOfTensorProductsOperator``); leaves MUST be
:class:`BoundaryTraceLaw` instances. Unknown node types raise
:class:`~orpheus.geometry.boundary.BoundaryError` naming the offending
type.

Placement note
==============

The walker currently lives at :mod:`orpheus.sn.boundary_realize`
because the only shipped functional realiser is
:class:`~orpheus.sn.boundary_realizer.SNBoundaryRealizer`. Per the
"unify after two instances" memory rule, the cross-method
generalisation (a method-agnostic walker that resolves the realiser
via :class:`~orpheus.geometry.boundary.BoundaryRealizerRegistry`) is
deferred until the second functional realiser (MoC, MC, CP, or
diffusion) ships. When that happens, the walker will move to
``orpheus.geometry.boundary.realize`` or similar.

References
----------

* ``.claude/plans/transient-giggling-cake.md`` Wave 11 — C11.2.
* Grand Report v3 §16A.3 lines 2841–2860 (realiser-as-third-layer
  motivation), §16A.10 (boundary-trace decomposition).
* The Wave-0 operator algebra: :mod:`orpheus.numerics.operator`.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from orpheus.geometry.boundary import BoundaryError, BoundaryTraceLaw
from orpheus.numerics.operator import (
    LinearOperator,
    OperatorProduct,
    OperatorSum,
    ScaledOperator,
    SumOfTensorProductsOperator,
    TensorProductOperator,
)

if TYPE_CHECKING:
    from orpheus.sn.method_space import SNMethodSpace


__all__ = ["realize_recursively"]


def realize_recursively(
    op: LinearOperator | BoundaryTraceLaw,
    method_space: "SNMethodSpace",
) -> LinearOperator:
    r"""Realise a tree of :class:`BoundaryTraceLaw` leaves + Wave-0 composers.

    Walks the expression tree, realising every
    :class:`~orpheus.geometry.boundary.BoundaryTraceLaw` leaf via the
    :class:`~orpheus.sn.boundary_realizer.SNBoundaryRealizer` and
    reassembling the result through the same Wave-0 composers
    (preserving the tree shape). Pure tree-walking — no math; the
    structurally-independent reference for tests is the leaf
    realisation (``SNBoundaryRealizer().realize(leaf, method_space)``)
    plus the Wave-0 algebra of the composers, neither of which this
    function touches.

    Parameters
    ----------
    op : :class:`LinearOperator` or :class:`BoundaryTraceLaw`
        The expression tree. Leaves MUST be
        :class:`BoundaryTraceLaw` instances. Internal nodes are
        recognised Wave-0 composers:
        :class:`~orpheus.numerics.operator.OperatorSum`,
        :class:`~orpheus.numerics.operator.OperatorProduct`,
        :class:`~orpheus.numerics.operator.ScaledOperator`,
        :class:`~orpheus.numerics.operator.TensorProductOperator`,
        :class:`~orpheus.numerics.operator.SumOfTensorProductsOperator`.

    method_space : :class:`SNMethodSpace`
        Method space passed verbatim to
        :meth:`SNBoundaryRealizer.realize` at every leaf.

    Returns
    -------
    :class:`LinearOperator`
        The realised 1-arg operator with the composer structure of
        ``op`` preserved (every leaf replaced by its realisation).

    Raises
    ------
    :class:`BoundaryError`
        If a tree node is neither a :class:`BoundaryTraceLaw` nor one
        of the five recognised Wave-0 composers. The error's ``law``
        field carries the offending type's name.

    Examples
    --------
    Realise a Wave-0 ``ScaledOperator``-around-leaf::

        >>> from orpheus.geometry.boundary import ReflectiveBoundary
        >>> from orpheus.sn.boundary_realizer import (
        ...     SNBoundaryRealizer, SNMethodSpace,
        ... )
        >>> from orpheus.sn.quadrature import GaussLegendre1D
        >>> ms = SNMethodSpace.minimal(GaussLegendre1D.create(4))
        >>> law = 0.5 * ReflectiveBoundary(axis="x")
        >>> realised = realize_recursively(law, ms)  # ScaledOperator
        ...                                          # around realised
        ...                                          # PermutationOperator

    Realise the standard Marshak mixed boundary::

        >>> from orpheus.geometry.boundary import WhiteBoundary
        >>> law = (
        ...     0.3 * ReflectiveBoundary(axis="x")
        ...     + 0.7 * WhiteBoundary(axis="x", outward_sign=+1)
        ... )
        >>> realised = realize_recursively(law, ms)  # OperatorSum
        ...                                          # of ScaledOperators
    """
    # Import lazily to avoid the import cycle
    # orpheus.geometry.boundary → orpheus.sn.boundary_realizer →
    # orpheus.geometry.boundary. The SN realiser is the only consumer
    # of ``SNMethodSpace`` and the walker is itself SN-specific today.
    from orpheus.sn.boundary_realizer import SNBoundaryRealizer

    if isinstance(op, BoundaryTraceLaw):
        # Leaf: dispatch through the SN realiser.
        return SNBoundaryRealizer().realize(op, method_space)

    if isinstance(op, ScaledOperator):
        # ScaledOperator exposes ``scalar`` + ``op``. Recurse on the
        # inner operator and wrap the result with the same scalar.
        inner = realize_recursively(op.op, method_space)
        return ScaledOperator(op.scalar, inner)

    if isinstance(op, OperatorSum):
        # OperatorSum is binary: exposes ``a`` and ``b``. Recurse on
        # each operand and reassemble via the same binary composer.
        a = realize_recursively(op.a, method_space)
        b = realize_recursively(op.b, method_space)
        return OperatorSum(a, b)

    if isinstance(op, OperatorProduct):
        # OperatorProduct is binary: exposes ``a`` and ``b``.
        a = realize_recursively(op.a, method_space)
        b = realize_recursively(op.b, method_space)
        return OperatorProduct(a, b)

    if isinstance(op, TensorProductOperator):
        # TensorProductOperator exposes ``ops`` (tuple of factors).
        ops = tuple(realize_recursively(o, method_space) for o in op.ops)
        return TensorProductOperator(ops)

    if isinstance(op, SumOfTensorProductsOperator):
        # SumOfTensorProductsOperator exposes ``summands`` (tuple of
        # TensorProductOperator instances). Each summand IS a
        # TensorProductOperator (constructor enforces this); recursive
        # realisation of a TensorProductOperator returns a
        # TensorProductOperator (per the branch above), so the result
        # remains a valid SumOfTensorProductsOperator.
        summands = tuple(
            realize_recursively(s, method_space) for s in op.summands
        )
        return SumOfTensorProductsOperator(summands)

    raise BoundaryError(
        f"realize_recursively cannot handle node of type "
        f"{type(op).__name__}. Expected BoundaryTraceLaw (leaf) or a "
        f"recognised Wave-0 composer (OperatorSum, OperatorProduct, "
        f"ScaledOperator, TensorProductOperator, "
        f"SumOfTensorProductsOperator).",
        law=type(op).__name__,
    )
