r"""The sweep inverse operator :math:`(L+C)^{-1}` — the inverse of an
:class:`~orpheus.sn.operators.streaming.InvertibleOperator`, in operator form.

The inverse-as-operator carve (#226) replaces ``A.solve(b)`` (a gated method
call) with ``A.inverse().apply(b)`` / ``A.inverse() @ b`` (apply the inverse
OPERATOR) — so the forward view ``A.apply`` and the inverse view
``A.inverse().apply`` become the two views of ONE operator, exactly the way
``A`` and ``A.H`` are. ``(L+C).inverse()`` returns this :class:`SweepOperator`.

Per the Grand Report v3 §5.7 operator taxonomy, a structured / triangular
operator returns a DIRECT-substitution inverse (this — the WDD forward sweep);
a general sum returns a Krylov ``GreenOperator`` (deferred — Pattern 6, no
consumer yet). Sweep-vs-Krylov is simply *which kind* of inverse operator
``.inverse()`` returns; the interface (``apply``) is the same.

**Operator vs application context (Grand Report §38).** The semantic object is
"the inverse of :math:`(L+C)`"; the WDD sweep's per-cell coefficient cache is its
APPLICATION context, not its identity — the SAME ``SweepOperator`` inverts any
``rhs``. So this is a thin typed wrapper that DELEGATES to
:meth:`~orpheus.sn.operators.streaming.InvertibleOperator.solve`, NOT a re-home
of the sweep machinery (which stays on the forward operator where the σ lives).

**Deferred (no Phase-2 consumer).** The adjoint-inverse
``A.H.inverse() == A.inverse().H`` (the transpose-solve) is issue #280, and the
round-trip ``A.inverse().inverse() == A`` has no consumer yet, so ``.H`` /
``.inverse()`` / ``is_invertible`` / ``is_adjointable`` are NOT provided here:
they inherit the base :class:`~orpheus.numerics.operator.LinearOperator` defaults
(``is_*`` ``False``; ``.H`` raises until #280). The operator advertises only
:pydata:`~orpheus.numerics.operator.CAP_APPLY`.
"""
from __future__ import annotations

from typing import TYPE_CHECKING, Optional

from orpheus.numerics.operator import CAP_APPLY, LinearOperator

if TYPE_CHECKING:
    from orpheus.numerics.space import FunctionSpace
    from orpheus.sn.operators.streaming import InvertibleOperator
    from orpheus.transport.full_field import FullField
    from orpheus.transport.timed_full_field import TimedFullField


class SweepOperator(LinearOperator["FullField"]):
    r"""The inverse operator :math:`(L+C)^{-1}` of an :class:`InvertibleOperator`.

    :meth:`apply` runs the WDD forward-substitution sweep by delegating to
    ``inner.solve`` — BIT-IDENTICAL to ``(L+C).solve(rhs, initial_guess=...)``.
    Endomorphic on the composite ``FullField`` carrier (an inverse swaps
    domain/codomain, which are equal here because ``L+C`` is endomorphic).
    """

    #: The inverse advertises only ``apply``; ``solve`` / ``apply_transpose`` are
    #: deferred (the adjoint-inverse is #280). Plain class attr (NOT a field).
    capabilities = frozenset({CAP_APPLY})

    def __init__(self, inner: "InvertibleOperator") -> None:
        #: The forward operator :math:`(L+C)` this is the inverse of.
        self.inner = inner

    @property
    def domain(self) -> Optional["FunctionSpace"]:
        # An inverse maps the inner operator's CODOMAIN back to its DOMAIN;
        # InvertibleOperator is endomorphic, so both are the FullField space.
        return self.inner.codomain

    @property
    def codomain(self) -> Optional["FunctionSpace"]:
        return self.inner.domain

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
