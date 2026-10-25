r"""The sweep inverse operator :math:`(L+C)^{-1}` — the inverse of an
:class:`~orpheus.sn.operators.streaming.StreamingCollisionOperator`, in operator form.

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
:meth:`~orpheus.sn.operators.streaming.StreamingCollisionOperator.solve`, NOT a re-home
of the sweep machinery (which stays on the forward operator where the σ lives).

**The invertibility axis (taxonomy §13 I2, step 1).** The involution
``A.inverse().inverse() == A`` now has a consumer — the universal
functoriality gate — so :meth:`~SweepOperator.inverse` returns the wrapped
:class:`~orpheus.sn.operators.streaming.StreamingCollisionOperator` ITSELF (the law
holds by OBJECT IDENTITY), ``is_invertible`` is ``True``, and ``solve`` on the
inverse is the forward matvec ``inner.apply`` (solving
:math:`(L+C)^{-1} y = b` IS applying :math:`(L+C)`), keeping the faithfulness
contract (``is_invertible=True``, a working ``solve``) honest here too.

**The adjoint-inverse (#280 Phase 2.5c).** The sweep's ``apply_transpose``
is the inverse family's one law
(:meth:`~orpheus.numerics.operator.InverseWrapMixin.apply_transpose`): the
reverse-scan transpose-solve ``(A⁻¹)ᵀ = (Aᵀ)⁻¹``, the inner ``(L+C)``'s
:meth:`~orpheus.sn.operators.streaming.StreamingCollisionOperator.solve_transpose`
(the 2.5b reverse scan), and ``is_adjointable`` reads ``True`` exactly where
that forward is transpose-invertible (the 1-D scan family; multi-D is #616). This is what makes the swap law
``A.H.inverse() ≡ A.inverse().H`` (via
:meth:`~orpheus.numerics.operator.AdjointOperator.inverse`) an identity of
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
delegation; the canonical ``initial_guess`` keyword is accepted and dropped,
since a direct sweep has nothing to seed) and ``__repr__``; its ctor guard is
the nominal :class:`ScheduleTriangularForward` type (schedule-triangularity is
what makes the forward sweep-invertible, and it is declared, not checked).
"""
from __future__ import annotations

from abc import ABC, abstractmethod
from typing import TYPE_CHECKING, Any, Optional

from orpheus.numerics.operator import InverseWrapMixin, LinearOperator

if TYPE_CHECKING:
    from orpheus.numerics.frame import FrameBase
    from orpheus.numerics.space import FunctionSpace
    from orpheus.transport.full_field import FullField
    from orpheus.transport.timed_full_field import TimedFullField


class ScheduleTriangularForward(ABC):
    r"""The NOMINAL type of a forward :math:`A` whose inverse is ONE forward-substitution scan.

    The schedule-triangular family: the plain WDD composite ``(L+C)``
    (:class:`~orpheus.sn.operators.streaming.StreamingCollisionOperator`) and the
    schedule-folded ``M = (L+C−B_lower)``
    (:class:`~orpheus.sn.operators.scheduled_invertible.ScheduledInvertibleOperator`,
    #226 step 2). Membership is declared by inheriting this class, never inferred
    from the members, because no member witnesses triangularity: a structural
    protocol over ``solve`` would admit any operator with a solve (a diagonal, a
    generic inverse), which the sweep does not invert. The abstract members are
    the ones the family's inverses read (the function-space pair the inverse
    swaps, the forward matvec, the exact ``solve`` :class:`SweepOperator`
    delegates to, and the composite solve body
    :class:`~orpheus.sn.operators.windowing.WindowedSweep` emits its moments
    through); the
    mixin's transpose arms read the forward through its
    :class:`~orpheus.numerics.operator.LinearOperator` face.

    The forward is upstream of its inverse: each member imports this module and
    its ``inverse()`` mints a :class:`SweepOperator` on itself, and this module
    names no member.
    """

    @property
    @abstractmethod
    def domain(self) -> Optional["FunctionSpace"]: ...

    @property
    @abstractmethod
    def codomain(self) -> Optional["FunctionSpace"]: ...

    @abstractmethod
    def apply(self, x: Any, /) -> Any: ...

    @abstractmethod
    def solve(self, b: Any, /) -> Any: ...

    @abstractmethod
    def _solve_timed_full_field(
        self, rhs: "FullField", *, moment_frame: "FrameBase | None" = None,
    ) -> "TimedFullField":
        """The family's composite solve body, which the windowed sweep reads to emit its moment frame."""


class SweepOperator(
    InverseWrapMixin[ScheduleTriangularForward], LinearOperator["FullField"],
):
    r"""The inverse operator :math:`A^{-1}` of a schedule-triangular forward
    ``A`` — :class:`StreamingCollisionOperator` ``(L+C)`` or
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
    The transpose axis is the mixin's law
    :math:`(A^{-1})^{\mathsf T} = (A^{\mathsf T})^{-1}` (#280 Phase 2.5c): the
    sweep transposes iff its forward is transpose-invertible
    (:attr:`~orpheus.numerics.operator.LinearOperator.is_transpose_invertible`),
    and its ``apply_transpose`` is the forward's reverse scan
    :meth:`~orpheus.sn.operators.streaming.StreamingCollisionOperator.solve_transpose`.
    For ``(L+C)`` that needs the matvec transpose's factors AND a loss
    representation that realises the reverse scan (the 1-D scan family;
    multi-D defers). Read off the matvec factors alone, the flag was ``True``
    on 2-D poses whose ``.H.apply`` raised (`[M]` 2026-10-18,
    ``tests/gates/numerics/test_adjointability_is_honest.py``). The
    schedule-folded
    :class:`~orpheus.sn.operators.scheduled_invertible.ScheduledInvertibleOperator`
    arm has no reverse scan (the #280 sibling deferral, no consumer), so a
    sweep over it reads ``False`` and its ``apply_transpose`` raises
    :class:`~orpheus.numerics.operator.MissingAdjoint`. The JOINT adjoint of
    a carrying mesh lives on the M grid
    (:meth:`~orpheus.numerics.coupled_system.CoupledSubstitutionOperator.apply_transpose`,
    the transposed block substitution).
    """

    def apply(
        self, rhs: "FullField", *, initial_guess: "FullField | None" = None,
    ) -> "TimedFullField":
        r"""Return :math:`(L+C)^{-1}\,\text{rhs}` via the WDD sweep.

        Delegates to
        :meth:`~orpheus.sn.operators.streaming.StreamingCollisionOperator.solve`
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

    def __repr__(self) -> str:
        return f"SweepOperator({self.inner!r})"
