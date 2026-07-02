r"""The sweep inverse operator :math:`(L+C)^{-1}` — the inverse of an
:class:`~orpheus.sn.operators.streaming.InvertibleOperator`, in operator form.

The inverse-as-operator carve (#226) replaces ``A.solve(b)`` (a gated method
call) with ``A.inverse().apply(b)`` / ``A.inverse() @ b`` (apply the inverse
OPERATOR) — so the forward view ``A.apply`` and the inverse view
``A.inverse().apply`` become the two views of ONE operator, exactly the way
``A`` and ``A.H`` are. ``(L+C).inverse()`` returns this :class:`SweepOperator`.

Per the taxonomy (§12 step 4, landed), a structured / triangular operator
returns a DIRECT-substitution inverse (this — the WDD forward sweep); a
general sum returns the preconditioned-splitting
:class:`~orpheus.numerics.green_operator.GreenOperator`. Sweep-vs-Green is
*which mathematical OBJECT* ``.inverse()`` returns, keyed by the operand's
structure (the realization — direct substitution vs Richardson splitting —
rides the object, never the name); the interface (the canonical seeded
``apply``) is the same.

**Operator vs application context (Grand Report §38).** The semantic object is
"the inverse of :math:`(L+C)`"; the WDD sweep's per-cell coefficient cache is its
APPLICATION context, not its identity — the SAME ``SweepOperator`` inverts any
``rhs``. So this is a thin typed wrapper that DELEGATES to
:meth:`~orpheus.sn.operators.streaming.InvertibleOperator.solve`, NOT a re-home
of the sweep machinery (which stays on the forward operator where the σ lives).

**The invertibility axis (taxonomy §13 I2, step 1).** The involution
``A.inverse().inverse() == A`` now has a consumer — the universal
functoriality gate — so :meth:`~SweepOperator.inverse` returns the wrapped
:class:`~orpheus.sn.operators.streaming.InvertibleOperator` ITSELF (the law
holds by OBJECT IDENTITY), ``is_invertible`` is ``True``, and ``solve`` on the
inverse is the forward matvec ``inner.apply`` (solving
:math:`(L+C)^{-1} y = b` IS applying :math:`(L+C)`), keeping the faithfulness
keystone ``is_invertible ≡ CAP_SOLVE`` honest here too.

**Still deferred (no consumer).** The adjoint-inverse
``A.H.inverse() == A.inverse().H`` (the transpose-solve) is issue #280 —
``.H`` / ``is_adjointable`` inherit the base
:class:`~orpheus.numerics.operator.LinearOperator` defaults until then.

**Wrap-delegate back-half.** The back-half (``capabilities`` /
domain↔codomain swap / ``solve→inner.apply`` / ``is_invertible`` /
``inverse()→inner``) is inherited from
:class:`~orpheus.numerics.operator.InverseWrapMixin` — the extraction this
class and :class:`~orpheus.numerics.operator.InverseOperator` recorded as
their collapse trigger, fired by the 3rd sibling
(:class:`~orpheus.numerics.green_operator.GreenOperator`, taxonomy §12
step 4). This class keeps only :meth:`~SweepOperator.apply` (the sweep
delegation, threading the Carlson seed) and ``__repr__``; its ctor guard
is the ``SweepInvertible`` TYPE itself (schedule-triangularity is what
makes the forward sweep-invertible — no value check needed).
"""
from __future__ import annotations

from typing import TYPE_CHECKING, Union

from orpheus.numerics.operator import InverseWrapMixin, LinearOperator

if TYPE_CHECKING:
    from orpheus.sn.operators.scheduled_invertible import (
        ScheduledInvertibleOperator,
    )
    from orpheus.sn.operators.streaming import InvertibleOperator
    from orpheus.transport.full_field import FullField
    from orpheus.transport.timed_full_field import TimedFullField

    #: The sweep-invertible forwards this class inverts: the plain WDD
    #: composite ``(L+C)`` and the schedule-folded ``M = (L+C−B_lower)``
    #: (#226 step 2) — one wrapper for the whole schedule-triangular family.
    SweepInvertible = Union[InvertibleOperator, ScheduledInvertibleOperator]


class SweepOperator(
    InverseWrapMixin["SweepInvertible"], LinearOperator["FullField"],
):
    r"""The inverse operator :math:`A^{-1}` of a schedule-triangular forward
    ``A`` — :class:`InvertibleOperator` ``(L+C)`` or
    :class:`ScheduledInvertibleOperator` ``M = (L+C−B_lower)`` (#226 step 2).

    :meth:`apply` runs the forward-substitution sweep by delegating to
    ``inner.solve`` — BIT-IDENTICAL to ``inner.solve(rhs, initial_guess=...)``.
    Endomorphic on the composite ``FullField`` carrier (an inverse swaps
    domain/codomain, which are equal here because the forward is endomorphic).

    The wrap-delegate back-half — ``capabilities`` (``apply`` inverts,
    ``solve`` un-inverts; ``apply_transpose`` stays deferred, the
    adjoint-inverse is #280), the domain↔codomain swap, ``solve`` = the
    forward matvec, and the object-identity involution ``inverse() →
    inner`` — is inherited from
    :class:`~orpheus.numerics.operator.InverseWrapMixin`.
    """

    def apply(
        self, rhs: "FullField", *, initial_guess: "FullField | None" = None,
    ) -> "TimedFullField":
        r"""Return :math:`(L+C)^{-1}\,\text{rhs}` via the WDD sweep.

        Delegates to
        :meth:`~orpheus.sn.operators.streaming.InvertibleOperator.solve`
        (bit-identical); ``initial_guess`` threads the previous-iterate seed
        (the curvilinear Morel–Montry Carlson coupled-pole start) through.
        """
        return self.inner.solve(rhs, initial_guess=initial_guess)

    def __repr__(self) -> str:
        return f"SweepOperator({self.inner!r})"
