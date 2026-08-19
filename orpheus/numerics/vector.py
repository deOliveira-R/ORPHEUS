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

Everything *summed at a given stage* of that loop shares one carrier
type — the closed source/sink algebra adds sources to sources (the
``g_i ψ`` terms and ``q``), the flux space adds fluxes to fluxes (the
``ψ`` iterates). :class:`Vector` is the structural type of *a* carrier;
it does NOT assert that an operator maps a carrier to the *same*
carrier. The operators are honestly
:class:`~orpheus.numerics.operator.LinearOperator`\ ``[Din, Cout]`` —
``S`` and ``F`` map a flux carrier to a *distinct* source/sink carrier,
and ``L`` maps a flux composite to a source composite (see
:ref:`operator-algebra`). The endomorphic majority (``C``, the loss
solve) is the special case ``Din == Cout``.

Why a *structural* (duck-typed) Protocol, not an ABC
----------------------------------------------------

Three different things satisfy :class:`Vector` for free, none of which
could share an inheritance ancestor:

* ``np.ndarray`` — has ``__add__`` / ``__sub__`` / scalar ``__mul__``
  natively. It is the serialization wire format and the element type of
  the flat axis-primitives (``DiagonalOperator``, ``PermutationOperator``,
  …) that act on one tagged axis of a field's ``.values``.
* every :class:`~orpheus.numerics.field.Field` leaf — ``AngularFlux``,
  ``ScalarFlux``, ``HarmonicMomentFlux``, the boundary leaves — via the
  dunders gated by ``Field._check_partner`` (class-identity-is-units,
  same-space, mesh-bound where the leaf is — the fiber discipline; the
  #208 affine torsor gate retired at campaign 1 CS3).
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
The :data:`V` type variable is the structural carrier placeholder; the
:class:`~orpheus.numerics.operator.LinearOperator` Protocol is
parameterised over a *pair* ``[Din, Cout]`` so :meth:`apply` is typed
``apply(self, x: Din) -> Cout`` — honest that an operator may map one
carrier (flux) to a *distinct* carrier (source/sink), while ``[V] ≡
[V, V]`` (a PEP-696 ``Cout = Din`` default) keeps the endomorphic
majority spelled with a single parameter.
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

    Satisfied — without inheritance — by ``np.ndarray`` **at runtime**, by
    every :class:`~orpheus.numerics.field.Field` leaf, and by
    :class:`~orpheus.transport.timed_full_field.TimedFullField`.

    .. warning::

       **``np.ndarray`` does NOT satisfy this protocol STATICALLY** — measured
       2026-07-31 against pyright + the shipped numpy stubs. The members here
       are ``Self``-typed, and numpy's overloaded ``__add__`` will not bind to
       ``Self``: the checker tries ``ndarray[_AnyShape, dtype[bool]]`` and
       reports ``dtype[float64]`` is not a subtype of it. Consequences worth
       knowing before you fight one of them:

       * Operators that speak bare arrays are declared **unparameterized**
         (``class PermutationOperator(LinearOperator)``), which sidesteps the
         bind entirely — this is why they typecheck and is a deliberate
         spelling, not an omission.
       * A **generic** operator's callback hooks have no such escape. The one
         instance in the tree (``ZeroOperator``'s ``codomain_zero`` /
         ``transpose_zero`` used at the SN boundary trace) therefore carries a
         narrowed ``# type: ignore[reportArgumentType]`` naming this gap.

       Do not "fix" a call site by weakening the protocol; the gap is upstream
       in numpy's stubs, and the runtime conformance the docstring claims is
       genuine.

    The
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


#: Type variable bound to :class:`Vector` — the structural carrier
#: placeholder. It is the *invariant* parameter the operator-algebra
#: composers (:class:`~orpheus.numerics.operator.OperatorSum`,
#: :class:`~orpheus.numerics.operator.OperatorProduct`, …) are generic
#: over, and the input side of the honest two-parameter
#: :class:`~orpheus.numerics.operator.LinearOperator`\ ``[Din, Cout]``
#: Protocol (whose ``Cout`` defaults to the input, so the endomorphic
#: majority — ``C``, the loss solve — is still spelled ``[V]``).
#: The single ``V`` is NOT an endomorphism claim: ``S``/``F`` genuinely
#: map a flux carrier to a *distinct* source/sink carrier
#: (``apply(x: Din) -> Cout`` with ``Cout ≠ Din``). #208 made that role
#: change a typed fact; P4.5 (#65) made it the operator's *type*.
V = TypeVar("V", bound=Vector)

#: The DRIVER-BOUNDARY carrier placeholder — deliberately UNBOUNDED.
#: Where ``V`` constrains the operator algebra's internal composition
#: slots to :class:`Vector`, the Protocol-conformance boundary
#: (:class:`~orpheus.numerics.eigenvalue.EigenvalueSolver` and its
#: consumer :func:`~orpheus.numerics.eigenvalue.power_iteration`) must
#: bind a bare ``np.ndarray`` carrier concretely — five solver families
#: drive it with plain arrays.  (:class:`~orpheus.numerics.iteration
#: .KEigenvalue` itself stays on the bounded ``V``: a class-scoped
#: TypeVar infers per-call and is never forced against the ndarray
#: stubs — only the Protocol solve is.)  ``np.ndarray`` satisfies
#: :class:`Vector` at runtime, but numpy's stubs cannot prove the
#: structural match (the ``__add__`` bool-dtype self-binding overloads
#: reject it), so a ``bound=Vector`` Protocol would statically reject
#: every ndarray implementer the moment the TypeVar is solved.  The
#: missing bound is that stub limitation, NOT a wider contract — the
#: runtime contract is identical to ``V``'s (#276 A4).  Distinct from
#: :mod:`~orpheus.numerics.space`'s ``Carrier`` (a ``default=Any``
#: payload slot on :class:`~orpheus.numerics.space.FunctionSpace`) —
#: same name, different axis: this one is the driver's ITERATE type.
Carrier = TypeVar("Carrier")
