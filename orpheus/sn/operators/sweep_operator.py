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

**Wrap-delegate back-half twin (collapse trigger).** The back-half here
(``capabilities`` / domain↔codomain swap / ``solve→inner.apply`` /
``is_invertible`` / ``inverse()→inner``) is deliberately byte-identical to
:class:`~orpheus.numerics.operator.InverseOperator`'s — two witnesses, kept
twinned per defer-until-≥2. TRIGGER: at the 3rd sibling
(``GreenOperator`` / ``MatrixInverseOperator``, taxonomy §12 steps 4–5),
extract a shared mechanism mixin; do not hand-re-derive it again.
"""
from __future__ import annotations

from typing import TYPE_CHECKING, Optional

from orpheus.numerics.operator import CAP_APPLY, CAP_SOLVE, LinearOperator

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

    #: ``apply`` inverts; ``solve`` un-inverts (the forward matvec, see
    #: :meth:`solve`) — so ``CAP_SOLVE`` rides with ``is_invertible`` (the
    #: faithfulness keystone). ``apply_transpose`` stays deferred (the
    #: adjoint-inverse is #280). Plain class attr (NOT a field).
    capabilities = frozenset({CAP_APPLY, CAP_SOLVE})

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

    def solve(self, b: "FullField", /) -> "FullField":
        r"""Solve :math:`(L+C)^{-1}\,y = b`, i.e. return :math:`(L+C)\,b`.

        The inverse of the inverse acts forward: this is the matvec
        ``inner.apply``, delegated — the ``CAP_SOLVE`` face that keeps the
        faithfulness keystone ``is_invertible ≡ CAP_SOLVE`` honest on this
        object (taxonomy §13 I2 / step 1).
        """
        return self.inner.apply(b)

    @property
    def is_invertible(self) -> bool:
        return True  # ((L+C)^{-1})^{-1} = (L+C) — the wrapped operator itself

    def inverse(self) -> "InvertibleOperator":
        r"""Return :math:`((L+C)^{-1})^{-1} = (L+C)` — the wrapped operator, by identity.

        The involution law (taxonomy §13 I2) holds as an OBJECT-IDENTITY
        fact: ``A.inverse().inverse() is A``.
        """
        return self.inner

    def __repr__(self) -> str:
        return f"SweepOperator({self.inner!r})"
