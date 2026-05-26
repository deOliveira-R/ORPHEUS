r"""Tree-walking realisation for descriptor-tree-composed boundary laws.

The :meth:`~orpheus.sn.boundary_realizer.SNBoundaryRealizer.realize`
method dispatches by ``isinstance(law, BoundaryTraceLaw)`` and returns
a 1-arg :class:`~orpheus.numerics.operator.LinearOperator`. But users
may build a boundary law by composing multiple BCs via the
descriptor-tree algebra:

.. code-block:: python

    law = (
        0.3 * ReflectiveBoundary(axis="x")
        + 0.7 * WhiteBoundary(axis="x", outward_sign=+1)
    )

Here ``law`` is a
:class:`~orpheus.geometry.boundary._composition.LawSum` containing
:class:`~orpheus.geometry.boundary._composition.LawScaled` wrappers
around the leaf
:class:`~orpheus.geometry.boundary.BoundaryTraceLaw` instances.
Passing this directly to :meth:`SNBoundaryRealizer.realize` raises
because the realiser does NOT dispatch on
:class:`LawSum` / :class:`LawScaled` — those are descriptor-tree
composers, not laws.

:func:`realize_recursively` is the **type-transformer**: it walks
the descriptor tree (whose nodes are
``BoundaryTraceLaw | LawSum | LawScaled``) and returns an operator
tree whose nodes are Wave-0
:class:`~orpheus.numerics.operator.OperatorSum` /
:class:`~orpheus.numerics.operator.ScaledOperator` composers around
1-arg realised
:class:`~orpheus.numerics.operator.LinearOperator` leaves:

.. code-block:: python

    realised = realize_recursively(law, method_space)
    # realised is OperatorSum(ScaledOperator(0.3, realised_spec),
    #                         ScaledOperator(0.7, realised_white))
    # where realised_spec and realised_white are 1-arg Wave-0
    # primitives consumable by the SN sweep / Krylov path.

This is the **only** place a descriptor becomes an operator. The
§16A.3 three-layer architecture (descriptor / realizer / operator)
is enforced by the type system: the input lives in the descriptor
type family, the output in the operator type family, and the
function name advertises the transformation.

Wave-11 background
==================

Pre-Wave-11, rank-N composition was handled by a dedicated
``MixedBoundaryOperator`` class whose ``apply`` delegated to a
realiser branch. Wave 11 deleted that composer, briefly replacing
it with raw Wave-0 ``OperatorSum`` / ``ScaledOperator`` over BC
descriptors. Issue #186 (B3 + β2) then formalised the composition
tree as its own type family
(:class:`~orpheus.geometry.boundary._composition.LawSum` /
:class:`~orpheus.geometry.boundary._composition.LawScaled`) so that
"this is a law, that is an operator" is checkable by inspection
rather than convention.

The walker dispatches on three node types only:

* :class:`BoundaryTraceLaw` — leaf, dispatches to
  :class:`SNBoundaryRealizer.realize`.
* :class:`LawScaled` — recurse on ``inner``, wrap with
  :class:`ScaledOperator`.
* :class:`LawSum` — recurse on ``a`` and ``b``, reassemble with
  :class:`OperatorSum`.

Any other node raises :class:`TypeError`. In particular, Wave-0
operator-tree composers (:class:`OperatorProduct`,
:class:`TensorProductOperator`,
:class:`SumOfTensorProductsOperator`) are NOT valid inputs — those
are *operator*-tree composers, not *law*-tree composers, and they
cannot appear in the descriptor side. Tests that previously walked
an operator-tree with BC leaves are migrated to the new descriptor
form in C-B3.12.

Placement note
==============

The walker currently lives at :mod:`orpheus.sn.boundary_realize`
because the only shipped functional realiser is
:class:`~orpheus.sn.boundary_realizer.SNBoundaryRealizer`. Per the
"unify after two instances" memory rule, the cross-method
generalisation (a method-agnostic walker that resolves the realiser
via :class:`~orpheus.geometry.boundary.BoundaryRealizerRegistry`) is
deferred until the second functional realiser (MoC, MC, CP, or
diffusion) ships.

References
----------

* ``.claude/plans/bc-trace-law-descriptor-cleanup.md`` Issue #186
  (B3 + β2) — this rewrite.
* ``.claude/plans/transient-giggling-cake.md`` Wave 11 — C11.2
  (the predecessor walker over Wave-0 composers).
* Grand Report v3 §16A.3 lines 2841–2860 (realiser-as-third-layer
  motivation), §16A.10 (boundary-trace decomposition).
* The Wave-0 operator algebra:
  :mod:`orpheus.numerics.operator`.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from orpheus.geometry.boundary import BoundaryTraceLaw, LawScaled, LawSum
from orpheus.numerics.operator import (
    LinearOperator,
    OperatorSum,
    ScaledOperator,
)

if TYPE_CHECKING:
    from orpheus.geometry.boundary import LawNode
    from orpheus.sn.method_space import SNMethodSpace


__all__ = ["realize_recursively"]


def realize_recursively(
    node: "LawNode",
    method_space: "SNMethodSpace",
) -> LinearOperator:
    r"""Walk a descriptor tree and realise it as a Wave-0 operator tree.

    The descriptor tree's leaves are
    :class:`~orpheus.geometry.boundary.BoundaryTraceLaw` instances;
    internal nodes are
    :class:`~orpheus.geometry.boundary._composition.LawSum` /
    :class:`~orpheus.geometry.boundary._composition.LawScaled`. The
    output tree has the same shape with
    :class:`~orpheus.numerics.operator.OperatorSum` /
    :class:`~orpheus.numerics.operator.ScaledOperator` composers and
    1-arg :class:`~orpheus.numerics.operator.LinearOperator` leaves
    (realised by
    :class:`~orpheus.sn.boundary_realizer.SNBoundaryRealizer`).

    This is the ONLY function that transforms a descriptor into an
    operator. The §16A.3 three-layer split is type-checkable through
    this signature: input in the law-type family, output in the
    operator-type family.

    Parameters
    ----------
    node
        Descriptor-tree root. Must be a
        :class:`BoundaryTraceLaw` leaf or a :class:`LawSum` /
        :class:`LawScaled` composer.
    method_space
        Method space passed verbatim to
        :meth:`SNBoundaryRealizer.realize` at every leaf.

    Returns
    -------
    :class:`LinearOperator`
        The realised 1-arg operator with the composer structure of
        ``node`` preserved (every leaf replaced by its realisation).

    Raises
    ------
    TypeError
        If a node is not a :class:`BoundaryTraceLaw`,
        :class:`LawSum`, or :class:`LawScaled`. The error names the
        offending type.

    Examples
    --------
    Realise a scaled-leaf descriptor::

        >>> from orpheus.geometry.boundary import ReflectiveBoundary
        >>> from orpheus.sn.boundary_realizer import (
        ...     SNBoundaryRealizer, SNMethodSpace,
        ... )
        >>> from orpheus.numerics.quadrature import Quadrature
        >>> ms = SNMethodSpace.minimal(Quadrature.gauss_legendre(4))
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

    if isinstance(node, BoundaryTraceLaw):
        # Leaf: dispatch through the SN realiser.
        return SNBoundaryRealizer().realize(node, method_space)

    if isinstance(node, LawScaled):
        # Recurse on the inner descriptor, wrap with ScaledOperator.
        inner_op = realize_recursively(node.inner, method_space)
        return ScaledOperator(node.scalar, inner_op)

    if isinstance(node, LawSum):
        # Recurse on both operands, reassemble via OperatorSum.
        a_op = realize_recursively(node.a, method_space)
        b_op = realize_recursively(node.b, method_space)
        return OperatorSum(a_op, b_op)

    raise TypeError(
        f"realize_recursively expected BoundaryTraceLaw | LawSum | "
        f"LawScaled (the descriptor-tree type family), got "
        f"{type(node).__name__}. Operator-tree composers (e.g. "
        f"OperatorProduct, TensorProductOperator) are not valid "
        f"inputs — realize each descriptor first, then compose the "
        f"results via Wave-0 operator algebra."
    )
