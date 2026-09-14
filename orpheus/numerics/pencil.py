r"""The operator PENCIL — the degree-1 operator family :math:`\mathcal{A}(\sigma) = A - \sigma M`.

Layer 1 of the terminal object (consumers campaign step 2, R-cc2 / shape (B),
``.claude/plans/consumers_step2_design.md`` §3.5): the Problem's LAST step is the
PAIR ``(A, M)`` on ONE space — what an eigenvalue system IS (regular iff
:math:`\det(A - \lambda M) \not\equiv 0`; exactly :math:`\operatorname{rank} M`
finite eigenvalues by Weierstrass–Kronecker; the spectrum = the poles of the
resolvent :math:`R(\sigma) = \mathcal{A}(\sigma)^{-1}`) and what a shifted source
system consumes (``at(σ)`` handed to a source posing: the subcritical multiplying
source is ``SourcePosing(pencil.at(1), q)``; a Laplace/noise family is a PATH of
source posings over one pencil).  It carries NO physics, NO inverse and NO
resolvent method: the resolvent is the Strategy's composition per method (SN: the
splitting of the Problem's factors; diffusion/homogeneous:
``MatrixInverseOperator(pencil.at(σ))`` — never ``A.inverse()``; every ``at(σ)``
invalidates a factorisation the Strategy holds).

**Degree contract.** The family is AFFINE in :math:`\sigma` — ``at`` is linear,
``at(0)`` is the ``lhs`` object itself (``[M]`` ``ScaledOperator(0.0, ·)`` is
refused by the algebra, so the σ = 0 member cannot be spelled as a scaling), and
``at(σ + τ) − at(σ) == −τ·M`` as operators (bit-identically in matrix form; the
APPLIED form differs at the 1e-14 level — ``[M]`` 144/200 draws — because the
two sums associate differently).  A delayed-neutron α problem is RATIONAL in α;
it needs a ``linearize()`` seam (a companion pencil on an enlarged space), never a
silent assumption that this type's ``at`` is enough.

**Which methods can express their problem as a pencil.** SN, diffusion,
homogeneous, CP (``(I − K, F)`` over two opaque operators) — MC cannot: it samples
the Neumann series and holds no operator pair; stated, not papered over.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Generic, TypeVar

import numpy as np

from orpheus.numerics.operator import LinearOperator, Vector

V = TypeVar("V", bound=Vector)  # the operator algebra's own carrier bound


@dataclass(frozen=True)
class OperatorPencil(Generic[V]):
    """The pair ``(lhs, rhs)`` = :math:`(A, M)` on one space; ``at(σ) = A − σ·M``.

    The ENDS LAW at construction: both operators must agree on domain and on
    codomain (the first witness — ``OperatorPencil(rec.loss, hub.fission)`` on
    SN pairs a ``CoupledSpace`` with a ``FullFieldSpace`` and is REFUSED; the
    posed pair ``(loss, production)`` is admitted).  Ends an operator does not
    declare (``None``) are not compared — an OPAQUE pair (CP's) is admitted on
    the caller's word, which the docstring says rather than a flag.
    """

    lhs: LinearOperator[V, V]
    rhs: LinearOperator[V, V]
    degree: int = 1  # the affine contract (a ClassVar in spirit; a field so it is dataclass-visible)

    def __post_init__(self) -> None:
        for end in ("domain", "codomain"):
            a, m = getattr(self.lhs, end), getattr(self.rhs, end)
            if a is not None and m is not None and a != m:
                raise ValueError(
                    f"OperatorPencil ends law: lhs.{end} and rhs.{end} must be ONE "
                    f"space; got {a!r} and {m!r}. A pencil is a pair of operators "
                    f"on one space — pose the rhs on the lhs's ends first."
                )

    # ── the family ──────────────────────────────────────────────────────
    def at(self, sigma: float) -> LinearOperator[V, V]:
        r"""The member :math:`\mathcal{A}(\sigma) = A - \sigma M`; ``at(0.0) is lhs``."""
        sigma = float(sigma)
        if sigma == 0.0:
            return self.lhs
        return self.lhs - sigma * self.rhs

    @property
    def H(self) -> "OperatorPencil[V]":
        r"""The adjoint pencil :math:`(A^\dagger, M^\dagger)` — the same spectrum, NO extra datum."""
        return OperatorPencil(self.lhs.H, self.rhs.H)

    # ── the queries ─────────────────────────────────────────────────────
    @property
    def is_regular(self) -> bool:
        r"""Regular iff :math:`\det(A - \lambda M) \not\equiv 0` — certified here by
        the σ = 0 member being invertible (a sufficient condition; the one the
        transport pencils satisfy because ``at(0)`` is what the sweep inverts)."""
        return bool(self.lhs.is_invertible)

    @property
    def rhs_rank(self) -> int:
        r"""The number of FINITE eigenvalues (Weierstrass–Kronecker) = :math:`\operatorname{rank} M`.

        Dense (``rhs.as_matrix()``): a diagnostic for the fixtures that can
        assemble; the SN composite refuses ``as_matrix`` (``[M]`` D-limit).
        """
        return int(np.linalg.matrix_rank(np.asarray(self.rhs.as_matrix())))

    def apply_at(self, sigma: float, x: Any) -> Any:
        """``at(σ).apply(x)`` spelled without materialising the member (the applied law's spelling)."""
        return self.at(sigma).apply(x)
