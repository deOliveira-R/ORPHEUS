r"""Descriptor-tree composition for boundary laws (Issue #186 / B3 + β2).

Two frozen dataclasses (:class:`LawSum`, :class:`LawScaled`) form a
closed algebra over the union type
``BoundaryTraceLaw | LawSum | LawScaled`` — collectively the
:data:`LawNode` alias. The tree is a **pure descriptor**: it has no
``apply`` method; calling it before realization is a static type
error (not a runtime ``AttributeError``).

The §15.2 sum-of-tensor-products form

.. math::

    R \;=\; \sum_\alpha c_\alpha \, G_\alpha

is realized at the SN level by walking the descriptor tree via
:func:`orpheus.sn.boundary_realizer.realize_recursively`, which turns
each leaf into a 1-arg Wave-0 :class:`LinearOperator` and re-assembles
them under Wave-0 :class:`~orpheus.numerics.operator.OperatorSum` /
:class:`~orpheus.numerics.operator.ScaledOperator` composers. The
descriptor → operator transformation is the **only** place an
unrealized law becomes callable.

Algebra closure
---------------

* ``LawScaled * scalar`` folds the scalar (``LawScaled(2, LawScaled(3, x))``
  → ``LawScaled(6, x)``); the inner node never re-nests.
* ``LawSum + node`` / ``node + LawSum`` produces a new :class:`LawSum`
  with the new node positioned as a sibling — the tree is **not**
  flattened. ``(a + b) + c`` is ``LawSum(LawSum(a, b), c)``, distinct
  from ``LawSum(a, LawSum(b, c))``. The realizer walks both shapes
  and the algebraic value is identical up to floating-point
  non-associativity; tests should not assert specific tree shapes
  beyond what is algebraically meaningful.
* ``LawSum - node`` rewrites as ``LawSum(self, LawScaled(-1, node))``.
* ``-LawSum`` wraps as ``LawScaled(-1, self)``.

Mixing :class:`LawNode` instances with already-realized
:class:`~orpheus.numerics.operator.LinearOperator` instances is
**not** supported: callers MUST realize the descriptor tree first
(via :func:`realize_recursively`) before composing with Wave-0
operator-tree composers.

References
----------

* Grand Report v3 §15.2 lines 2015–2044 — tensor-decomposition form.
* Grand Report v3 §16A.3 lines 2841–2860 — three-layer architecture
  (descriptor / realizer / operator).
* ``.claude/plans/bc-trace-law-descriptor-cleanup.md`` — this scope.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, Union

if TYPE_CHECKING:
    from ._base import BoundaryTraceLaw


LawNode = Union["BoundaryTraceLaw", "LawSum", "LawScaled"]
r"""Union type for descriptor-tree nodes: a :class:`BoundaryTraceLaw`
leaf or one of the two recognised internal composers (:class:`LawSum`,
:class:`LawScaled`)."""


__all__ = ["LawNode", "LawScaled", "LawSum"]


@dataclass(frozen=True)
class LawScaled:
    r"""Scalar coefficient times an inner law-descriptor node.

    Represents :math:`\alpha \cdot \mathrm{inner}` in the §15.2
    tensor-decomposition form. The :class:`LawScaled` algebra is
    closed and includes a **constant-folding** optimization:
    multiplying a :class:`LawScaled` by another scalar collapses the
    chain, so ``LawScaled(α, LawScaled(β, x))`` never appears at
    rest — it folds to ``LawScaled(α * β, x)`` at construction.

    Adding a :class:`LawScaled` to anything produces a :class:`LawSum`;
    subtracting rewrites via ``LawScaled(-1, other)``.

    Parameters
    ----------
    scalar : float
        Coefficient. Stored as :class:`float`; integer inputs are
        widened.
    inner : :data:`LawNode`
        The descriptor being scaled. Must itself be a
        :class:`~orpheus.geometry.boundary.BoundaryTraceLaw` leaf or
        a :class:`LawSum` node (a :class:`LawScaled` would have folded
        into ``self`` via the ``__mul__`` constructor).
    """

    scalar: float
    inner: "LawNode"

    def __add__(self, other: "LawNode") -> "LawSum":
        return LawSum(self, other)

    def __radd__(self, other: "LawNode") -> "LawSum":
        return LawSum(other, self)

    def __sub__(self, other: "LawNode") -> "LawSum":
        return LawSum(self, LawScaled(-1.0, other))

    def __rsub__(self, other: "LawNode") -> "LawSum":
        return LawSum(other, LawScaled(-1.0, self))

    def __mul__(self, scalar: float) -> "LawScaled":
        if not isinstance(scalar, (int, float)):
            return NotImplemented
        return LawScaled(float(scalar) * self.scalar, self.inner)

    __rmul__ = __mul__

    def __truediv__(self, scalar: float) -> "LawScaled":
        if not isinstance(scalar, (int, float)):
            return NotImplemented
        return LawScaled(self.scalar / float(scalar), self.inner)

    def __neg__(self) -> "LawScaled":
        return LawScaled(-self.scalar, self.inner)


@dataclass(frozen=True)
class LawSum:
    r"""Sum of two law-descriptor nodes.

    Represents :math:`a + b` in the §15.2 tensor-decomposition form.
    The :class:`LawSum` algebra is closed; the tree is **not**
    flattened (``(a + b) + c`` remains ``LawSum(LawSum(a, b), c)``).

    Parameters
    ----------
    a : :data:`LawNode`
        Left operand.
    b : :data:`LawNode`
        Right operand.
    """

    a: "LawNode"
    b: "LawNode"

    def __add__(self, other: "LawNode") -> "LawSum":
        return LawSum(self, other)

    def __radd__(self, other: "LawNode") -> "LawSum":
        return LawSum(other, self)

    def __sub__(self, other: "LawNode") -> "LawSum":
        return LawSum(self, LawScaled(-1.0, other))

    def __rsub__(self, other: "LawNode") -> "LawSum":
        return LawSum(other, LawScaled(-1.0, self))

    def __mul__(self, scalar: float) -> "LawScaled":
        if not isinstance(scalar, (int, float)):
            return NotImplemented
        return LawScaled(float(scalar), self)

    __rmul__ = __mul__

    def __truediv__(self, scalar: float) -> "LawScaled":
        if not isinstance(scalar, (int, float)):
            return NotImplemented
        return LawScaled(1.0 / float(scalar), self)

    def __neg__(self) -> "LawScaled":
        return LawScaled(-1.0, self)
