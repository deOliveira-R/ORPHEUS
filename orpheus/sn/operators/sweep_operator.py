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
contract (``is_invertible=True``, a working ``solve``) honest here too.

**The adjoint-inverse (#280 Phase 2.5c).** :meth:`~SweepOperator.apply_transpose`
is the reverse-scan transpose-solve ``(A⁻¹)ᵀ = (Aᵀ)⁻¹`` (delegating to the
inner ``(L+C)``'s :meth:`~orpheus.sn.operators.streaming.InvertibleOperator.solve_transpose`,
the 2.5b reverse-scan), and ``is_adjointable`` flips ``True`` on the
``InvertibleOperator`` arm. This is what makes the swap law
``A.H.inverse() ≡ A.inverse().H`` (via
:meth:`~orpheus.numerics.operator._AdjointOperator.inverse`) an identity of
the algebra. The schedule-folded
:class:`~orpheus.sn.operators.scheduled_invertible.ScheduledInvertibleOperator`
arm keeps the base ``is_adjointable = False`` (its reverse-scan is the #280
sibling deferral — no consumer).

**Wrap-delegate back-half.** The back-half (domain↔codomain swap /
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

from orpheus.numerics.operator import (
    InverseWrapMixin,
    LinearOperator,
    MissingAdjoint,
)

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

    The wrap-delegate back-half — (``apply`` inverts,
    ``solve`` un-inverts), the domain↔codomain swap, ``solve`` = the
    forward matvec, and the object-identity involution ``inverse() →
    inner`` — is inherited from
    :class:`~orpheus.numerics.operator.InverseWrapMixin`.
    :meth:`apply_transpose` / :attr:`is_adjointable` are the #280 Phase 2.5c
    adjoint-inverse wiring (the reverse-scan ``(A⁻¹)ᵀ``), added here.
    """

    def apply(
        self, rhs: "FullField", *, initial_guess: "FullField | None" = None,
    ) -> "TimedFullField":
        r"""Return :math:`(L+C)^{-1}\,\text{rhs}` via the WDD sweep.

        Delegates to
        :meth:`~orpheus.sn.operators.streaming.InvertibleOperator.solve`
        (bit-identical). ``initial_guess`` is the inverse family's canonical
        :class:`~orpheus.numerics.iteration.SupportsSeededApply` keyword — the
        driver threads the previous iterate UNIFORMLY every step (#285) — but
        the WDD sweep is a DIRECT, exact inverse: there is nothing to seed (the
        curvilinear :math:`\psi_{1/2}` starting direction is computed DIRECTLY
        from the source since #282 route (a), 2.5d). So it is accepted and
        DROPPED, exactly as the other exact inverses do
        (:class:`~orpheus.numerics.matrix_inverse_operator.MatrixInverseOperator`
        / :class:`~orpheus.sn.operators.windowing.WindowedSweep` /
        the numerics ``_SeededExactApply``). The warm start lives at the
        ITERATION layer (:meth:`~orpheus.numerics.iteration.SourceIteration.solve`'s
        ``initial_guess`` :math:`x_0`) and the ITERATIVE
        :class:`~orpheus.numerics.green_operator.GreenOperator`, never on a
        direct sweep.
        """
        del initial_guess  # direct exact inverse — nothing to seed (#282/2.5d)
        return self.inner.solve(rhs)

    @property
    def is_adjointable(self) -> bool:
        r"""Whether the reverse-scan transpose-solve ``(A⁻¹)ᵀ = (Aᵀ)⁻¹`` exists.

        ``True`` iff the inner forward is the ``(L+C)``
        :class:`~orpheus.sn.operators.streaming.InvertibleOperator` — which
        carries :meth:`~orpheus.sn.operators.streaming.InvertibleOperator.solve_transpose`
        (the 2.5b reverse-scan) — AND that forward is itself adjointable, so
        the two-factor geometry gate rides ``inner.is_adjointable`` (DD-1D
        yes; LD / multi-D Cartesian defer). The schedule-folded
        :class:`~orpheus.sn.operators.scheduled_invertible.ScheduledInvertibleOperator`
        arm has no ``solve_transpose`` (its reverse-scan is the #280 sibling
        deferral — no consumer), so a ``SweepOperator`` over it stays
        non-adjointable.
        """
        from orpheus.sn.operators.streaming import InvertibleOperator

        return (
            isinstance(self.inner, InvertibleOperator)
            and self.inner.is_adjointable
        )

    def apply_transpose(self, b: "FullField") -> "FullField":
        r"""Return ``(L+C)⁻ᵀ b`` — the reverse-scan transpose-solve (#280 2.5c).

        ``SweepOperator`` is the inverse operator ``A⁻¹ = (L+C)⁻¹`` (endomorphic
        on the composite), so its Euclidean transpose is
        ``(A⁻¹)ᵀ = (Aᵀ)⁻¹`` — exactly the inner's
        :meth:`~orpheus.sn.operators.streaming.InvertibleOperator.solve_transpose`
        (the 2.5b reverse-scan of the forward WDD sweep-scan). Wiring this is
        what makes the swap law ``A.H.inverse() ≡ A.inverse().H`` (via
        :meth:`~orpheus.numerics.operator._AdjointOperator.inverse`) an
        object identity of the algebra, and the metric adjoint-solve
        ``A.inverse().H.apply(b) = G⁺·apply_transpose(G·b)`` fall out of
        :meth:`~orpheus.numerics.operator._AdjointOperator.apply` for free.

        This is the plain EUCLIDEAN transpose; the metric conjugation of the
        physical G-adjoint ``.H`` is applied AROUND it by
        :class:`~orpheus.numerics.operator._AdjointOperator`.

        The ``isinstance`` narrowing is the direct-call backstop (the eager
        ``.H`` gate already refuses via :attr:`is_adjointable`): a
        ``SweepOperator`` over the schedule-folded
        :class:`~orpheus.sn.operators.scheduled_invertible.ScheduledInvertibleOperator`
        raises :class:`~orpheus.numerics.operator.MissingAdjoint` rather than
        an ``AttributeError`` on the missing ``solve_transpose``.
        """
        from orpheus.sn.operators.streaming import InvertibleOperator

        if not isinstance(self.inner, InvertibleOperator):
            raise MissingAdjoint(
                f"SweepOperator over {type(self.inner).__name__} has no "
                f"reverse-scan transpose-solve — #280 sibling scope wires only "
                f"the (L+C) InvertibleOperator arm; the schedule-folded "
                f"M = (L+C−B_lower) reverse-scan is deferred (no consumer)."
            )
        return self.inner.solve_transpose(b)

    def __repr__(self) -> str:
        return f"SweepOperator({self.inner!r})"
