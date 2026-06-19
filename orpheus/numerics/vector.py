r"""The vector type of the transport operator algebra.

The operator algebra in :mod:`orpheus.numerics.operator` reduces every
transport solve to compositions of a small set of linear operators
acting on a discrete state :math:`\psi`:

.. math::

    (L + C - S - F - B)\,\psi \;=\; q .

Historically the :class:`~orpheus.numerics.operator.LinearOperator`
Protocol typed that action as ``apply(x: np.ndarray) -> np.ndarray``.
That annotation is a *documented falsehood*: the production SN operators
consume and return a :class:`~orpheus.transport.timed_full_field.TimedFullField`
(``bulk`` :math:`\oplus` ``boundary`` :math:`\oplus` history), and the
flat ``np.ndarray`` appears **only** at the scipy-Krylov serialization
boundary (``to_flat`` / ``from_flat``). The algebra never actually
speaks ``np.ndarray``.

This module names the structural contract the algebra *does* speak. A
:class:`Vector` is anything that supports the vector-space operations
the source-iteration / Krylov / power-iteration loop requires — vector
addition, vector subtraction, left-multiplication by a scalar
(``scalar * ψ`` — the form every composer uses, see
:class:`~orpheus.numerics.operator.ScaledOperator` and
:class:`~orpheus.numerics.operator.ZeroOperator`), and division by a
scalar (the eigenvalue renormalisations ``F ψ / k`` and ``ψ / p``):

.. math::

    \psi_{n+1} \;=\; L^{-1}\!\Bigl(\textstyle\sum_i g_i\,\psi_n + q\Bigr) .

Everything that gets added in that loop (``g_i ψ``, ``q``, ``ψ``) must
be the same vector type, which is why the algebra is an *endomorphism*
algebra on one carrier (see :ref:`operator-algebra`). :class:`Vector`
is that carrier's type.

Why a *structural* (duck-typed) Protocol, not an ABC
----------------------------------------------------

Three different things satisfy :class:`Vector` for free, none of which
could share an inheritance ancestor:

* ``np.ndarray`` — has ``__add__`` / ``__sub__`` / scalar ``__mul__``
  natively. It is the serialization wire format and the element type of
  the flat axis-primitives (``DiagonalOperator``, ``PermutationOperator``,
  …) that act on one tagged axis of a field's ``.values``.
* every :class:`~orpheus.numerics.field.Field` leaf — ``AngularFlux``,
  ``ScalarFlux``, ``HarmonicMomentField``, the boundary leaves — via the
  dunders gated by ``Field._check_partner`` (class-identity-is-units,
  same-space, the #208 affine torsor gate).
* :class:`~orpheus.transport.timed_full_field.TimedFullField` — the
  composite carrier — via its delegate dunders that propagate to
  ``bulk`` + ``boundary``.

A nominal ABC would force ``np.ndarray`` to inherit from an ORPHEUS base
(impossible) and would buy nothing the runtime ``_check_partner`` gates
do not already enforce. The structural Protocol matches the existing
``LinearOperator`` style ("any object providing ``apply``") and promotes
the ad-hoc ``_is_ravellable`` duck-type in :mod:`orpheus.numerics.iteration`
to a named contract.

Layering note
-------------

:mod:`orpheus.numerics` sits *below* :mod:`orpheus.transport`, so this
Protocol cannot name ``TimedFullField`` directly. The structural
:class:`Vector` is precisely how a transport-level carrier conforms to a
numerics-level abstraction without ``numerics`` importing ``transport``.
The :data:`V` type variable lets :meth:`LinearOperator.apply` be typed
``apply(self, x: V) -> V`` — honest about the endomorphism, agnostic to
which concrete carrier (flux, scalar, or moment state) flows through.
"""

from __future__ import annotations

from typing import Protocol, Self, TypeVar, runtime_checkable

__all__ = ["Vector", "V"]


@runtime_checkable
class Vector(Protocol):
    r"""Structural contract for an element of the algebra's vector space.

    A :class:`Vector` supports the four vector-space operations the
    algebra and its iteration drivers compose: addition with a like
    vector, subtraction of a like vector, left-multiplication by a
    scalar, and division by a scalar. Scalar multiplication is contracted
    as ``__rmul__`` (``scalar * vector``) because that is the form the
    algebra actually invokes —
    :class:`~orpheus.numerics.operator.ScaledOperator` returns
    ``scalar * op.apply(x)`` and
    :class:`~orpheus.numerics.operator.ZeroOperator` returns ``0.0 * x``.
    Scalar division (``__truediv__``) is contracted because the
    eigenvalue drivers renormalise the carrier directly: ``F ψ / k`` in
    :class:`~orpheus.numerics.iteration.KEigenvalue` and the
    production-rate renormalisation ``ψ / p`` in
    :func:`~orpheus.numerics.eigenvalue.power_iteration`. The additive
    identity is obtained structurally as ``0.0 * x``, so no abstract
    ``zero`` member is required.

    The contract is the operations of a vector space (plus the
    by-a-scalar division the eigenvalue renormalisations spell as ``/``)
    — NOT carrier utilities such as ``.copy()`` or ``.to_flat()``; those
    belong to the concrete leaf, not to this structural type.

    Satisfied — without inheritance — by ``np.ndarray``, by every
    :class:`~orpheus.numerics.field.Field` leaf, and by
    :class:`~orpheus.transport.timed_full_field.TimedFullField`. The
    algebra's ``apply`` / ``solve`` speak :class:`Vector`, never
    ``np.ndarray`` (which is confined to the scipy serialization
    boundary).

    The contract is intentionally minimal — it pins only the operations
    the operator composers and iteration drivers need, leaving the units
    / space / mesh invariants to the concrete leaf's own
    ``_check_partner`` gate, which raises on any cross-class or
    cross-space combination at runtime.
    """

    def __add__(self, other: Self) -> Self:
        """Return ``self + other`` (vector addition with a like vector)."""
        ...

    def __sub__(self, other: Self) -> Self:
        """Return ``self - other`` (vector subtraction of a like vector)."""
        ...

    def __rmul__(self, scalar: float) -> Self:
        """Return ``scalar * self`` (left-multiplication by a scalar)."""
        ...

    def __truediv__(self, scalar: float) -> Self:
        """Return ``self / scalar`` (division by a scalar).

        Used by the eigenvalue drivers to renormalise the carrier in
        place: ``F ψ / k`` (:class:`~orpheus.numerics.iteration.KEigenvalue`)
        and ``ψ / p`` (:func:`~orpheus.numerics.eigenvalue.power_iteration`).
        """
        ...


#: Type variable bound to :class:`Vector`, used to type the endomorphism
#: signature ``apply(self, x: V) -> V`` on
#: :class:`~orpheus.numerics.operator.LinearOperator`. A single ``V``
#: (not a per-operator generic) is honest because every operator in the
#: ``(L + C - S - F - B)`` algebra maps a carrier to a same-typed
#: carrier; the inner leaf may change role (flux → source) at runtime,
#: but the vector type does not.
V = TypeVar("V", bound=Vector)
