r"""The functional type of the transport operator algebra (the §5.6 suffix law).

The grand-report §5.6 *suffix law* partitions the maps of the transport
algebra by what they produce:

* an **Operator** maps a field to a field, :math:`\psi \mapsto A\,\psi`
  (the :class:`~orpheus.numerics.operator.LinearOperator` surface
  ``apply(x: V) -> V``);
* a **Kernel** integrates a field against a measure on one or more axes
  (the :class:`~orpheus.numerics.operator.SumOfTensorProductsOperator`
  scattering / fission form — still field → field, but nonlocal);
* a **Functional** maps a field to a **scalar** (or, fiberwise over
  space, to a scalar-field) — the surface this module names.

:class:`Functional` is the co-vector companion of
:class:`~orpheus.numerics.vector.Vector`. Where a :class:`Vector` is an
element of the algebra's carrier space, a :class:`Functional` is a
*linear-or-nonlinear form* on that space: it consumes a carrier and
returns a number. The reaction-rate form
:math:`r = \langle \Sigma_a, \psi\rangle = \int \Sigma_a\,\phi\,dV`, the
Rayleigh-quotient :math:`k`-estimator, and the per-cell fission emission
density :math:`p(\vec r) = \sum_{g'} \nu\Sigma_{f,g'}\,\phi_{g'}` are all
Functionals in this sense — the contraction that ends in a scalar (or a
scalar field) is the suffix law's defining move.

Why two type variables (an input carrier *and* a result)
--------------------------------------------------------

The §5.6 suffix law is stated "field → scalar", and for a reaction-rate
or :math:`k`-eigenvalue functional the result is a scalar ``float``. But
the *per-cell* production rate
:math:`p(\vec r) = \sum_{g'} \nu\Sigma_{f,g'}(\vec r)\,\phi_{g'}(\vec r)`
is a Functional applied **fiberwise over space**: it contracts the group
axis only, returning a scalar-*field* (one number per cell), not a lone
scalar. So the result type is genuinely distinct from the input carrier —
the output of the per-cell production functional is a group-collapsed
density, not a member of the input flux carrier space. A single
``float | V`` union would be a *documented falsehood*: the group-collapsed
density is neither a ``float`` nor the input carrier (it has lost the group
axis the carrier carries). The result typevar :data:`R_co` is therefore
left **unbounded** — the concrete functional names what it returns
(``float`` for a reaction-rate / keff functional, an ``np.ndarray``
scalar-field for the per-cell density).

The input typevar :data:`V_contra` is bound to
:class:`~orpheus.numerics.vector.Vector` — the same carrier the operator
algebra acts on — and is **contravariant** (the variance of a consumed
argument), the co-vector mirror of the invariant
:data:`~orpheus.numerics.vector.V` (which the operator endomorphism both
consumes and returns). :data:`R_co` is **covariant** (the variance of a
produced result).

Why this is a SIBLING of LinearOperator, not a subclass
-------------------------------------------------------

A :class:`Functional` has exactly one method, :meth:`evaluate`. It
deliberately carries **none** of the
:class:`~orpheus.numerics.operator.LinearOperator` surface — no
``apply``, no ``capabilities``, no ``solve``, no ``apply_transpose``,
no ``H`` / ``domain`` / ``codomain``. The two Protocols share no member:
a Functional speaks ``evaluate`` (field → scalar/scalar-field); a
LinearOperator speaks ``apply`` (field → field) plus a capability set.
This disjointness IS the category's defining property — it is what makes
"the production rate is a Functional, not an operator" a *structural*
fact the type system enforces, rather than prose. The
intrinsic-property gate (``tests/transport/test_functional_category.py``)
verifies the partition in both directions: a Functional is not a
LinearOperator, and a LinearOperator is not a Functional.

Why a *structural* (duck-typed) Protocol, not an ABC
----------------------------------------------------

The same reasoning that makes :class:`~orpheus.numerics.vector.Vector` a
structural Protocol applies here. A transport-level concrete functional
(``ReactionRateFunctional`` in :mod:`orpheus.transport`) conforms to
this numerics-level abstraction *without inheritance* — it simply
exposes ``evaluate``. A nominal ABC would force every concrete functional
to inherit from a numerics base, which would either pull a transport
import down into ``numerics`` (a layering violation) or force the
transport leaf to subclass across the layer boundary. The structural
Protocol matches the existing :class:`~orpheus.numerics.operator.LinearOperator`
and :class:`~orpheus.numerics.vector.Vector` style ("any object providing
the method participates"), and ``@runtime_checkable`` lets
``isinstance(x, Functional)`` discriminate the category at runtime.

Layering note
-------------

:mod:`orpheus.numerics` sits *below* :mod:`orpheus.transport`, so this
Protocol names no transport type — it speaks only of the
:class:`~orpheus.numerics.vector.Vector` carrier (the input) and an
unbounded result :math:`R`. The structural :class:`Functional` is
precisely how a transport-level concrete functional (the §5.6
``ReactionRateFunctional``, which carries :math:`\Sigma_x` as a
:class:`~orpheus.transport.fields.cross_section_field.CrossSectionField`)
conforms to a numerics-level abstraction without ``numerics`` importing
``transport``.

References
----------

* Grand Report v3 §5.6 — the suffix law (Operator / Kernel / Functional
  / Projection / Reconstruction).
* ``.claude/agent-memory/cross-domain-attacker/coefficient_field_promotion_frames.md``
  — Frame 3 (the production rate as a Functional in
  :math:`F = M_\chi \circ \mathrm{ProductionRate} \circ M_{\nu\Sigma_f}`).
* ``.claude/plans/issue_257_coefficient_field_promotion.md`` — S5.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from typing import Protocol, TypeVar, runtime_checkable

from orpheus.numerics.field import Field
from orpheus.numerics.vector import Vector

__all__ = ["Functional", "InnerProductFunctional", "V_contra", "R_co"]

#: Input carrier type variable for :class:`Functional`, bound to
#: :class:`~orpheus.numerics.vector.Vector` (the same carrier the operator
#: algebra acts on). It is **contravariant**: a functional that accepts a
#: more general carrier is usable wherever a functional accepting a more
#: specific one is required — the standard variance of a consumed argument.
#: (This is the co-vector mirror of :data:`~orpheus.numerics.vector.V`,
#: which is invariant because the operator endomorphism both consumes and
#: returns it; a functional only *consumes* its carrier, so the input
#: typevar is contravariant.)
V_contra = TypeVar("V_contra", bound=Vector, contravariant=True)

#: Result type variable for :class:`Functional`. Deliberately **unbounded**
#: and **covariant** (the standard variance of a produced result) — a
#: Functional may return a scalar ``float`` (a reaction-rate or
#: :math:`k`-eigenvalue functional) or, fiberwise over space, a
#: scalar-*field* (the per-cell production density, an ``np.ndarray`` with
#: the group axis collapsed). The result is NOT in general the input
#: carrier :class:`~orpheus.numerics.vector.Vector` (the per-cell density
#: has lost the group axis), so a ``float | V`` union would mistype it.
R_co = TypeVar("R_co", covariant=True)


@runtime_checkable
class Functional(Protocol[V_contra, R_co]):
    r"""Structural contract for a §5.6 functional: a field → scalar map.

    A :class:`Functional` consumes a :class:`~orpheus.numerics.vector.Vector`
    carrier :math:`x` and returns a result :math:`R` (a scalar ``float``,
    or a scalar-field when the functional acts fiberwise over space). It
    is the co-vector companion of :class:`~orpheus.numerics.vector.Vector`
    and a *sibling* of :class:`~orpheus.numerics.operator.LinearOperator`,
    NOT a subclass: it carries the single :meth:`evaluate` method and
    deliberately none of the operator surface (no ``apply``, no
    ``capabilities``, no ``solve`` / ``apply_transpose`` / ``H`` /
    ``domain`` / ``codomain``). The two Protocols share no member —
    that disjointness is the category's defining property.

    Satisfied — without inheritance — by any object exposing
    ``evaluate``; the canonical transport instance is
    :class:`orpheus.transport.reaction_rate_functional.ReactionRateFunctional`,
    which carries :math:`\Sigma_x` and computes the per-cell reaction-rate
    density :math:`r_x(\vec r) = \sum_{g'} \Sigma_{x,g'}\,\phi_{g'}` (the
    production rate for :math:`\Sigma_x = \nu\Sigma_f`).
    """

    def evaluate(self, x: V_contra, /) -> R_co:
        r"""Return the functional's value on the carrier :math:`x`.

        The §5.6 suffix-law surface: a field → scalar (or, fiberwise over
        space, field → scalar-field) map. Mandatory. Every
        :class:`Functional` implements exactly this — it is the entire
        contract, and the deliberate absence of an ``apply`` /
        ``capabilities`` surface is what distinguishes the functional
        category from the operator category.
        """
        ...


@dataclass(frozen=True, eq=False)
class InnerProductFunctional:
    r"""The generic co-vector :math:`\langle w, \cdot\rangle` — the canonical §5.6 Functional.

    Contracts a carrier against a fixed **weight** :math:`w` over one axis:

    .. math::

        \langle w, x\rangle_j \;=\; \sum_{k} w_k\, x_k
        \qquad\text{(sum over }\texttt{axis}\text{)} .

    This is the concrete, numerics-level realisation of the structural
    :class:`Functional` Protocol — a co-vector (a *row*), the dual companion
    of a :class:`~orpheus.numerics.vector.Vector` (a *column*). Its defining
    role in the operator algebra is as the **row-factor of a rank-1 operator**:
    the dyad :math:`|v\rangle\langle w|` is built by
    :func:`~orpheus.numerics.operator.outer`\ ``(reconstruction, functional)``
    and applies as ``reconstruction * functional.evaluate(x)`` — the functional
    *is* the contraction, not a parallel description of one. The transport-level
    :class:`~orpheus.transport.reaction_rate_functional.ReactionRateFunctional`
    **specialises** this type, carrying the weight as a domain-typed
    :class:`~orpheus.transport.fields.cross_section_field.CrossSectionField`
    (production = ``⟨νΣf,·⟩``, absorption = ``⟨Σa,·⟩``) — exactly as
    :class:`~orpheus.transport.frames.harmonic_frame.HarmonicFrame` specialises
    :class:`~orpheus.numerics.frame.GalerkinFrame`.

    The contracted axis is kept as length-1 (``keepdims=True``) so the output
    ``(1, *complement)`` broadcasts cleanly against the rank-1 reconstruction
    (the dyad's column) — the leading-1 axis is load-bearing for that
    broadcast, not incidental.

    Parameters
    ----------
    weight : NDArray
        The co-vector weight :math:`w`. Shape must align with the carrier on
        the contracted ``axis`` (and broadcast on the complement). A
        :class:`~orpheus.numerics.field.Field` weight is accepted (its
        ``.values`` is read).
    axis : int, default 0
        The contraction axis. For the multigroup reaction rate this is the
        leading **group** axis (0); the spatial axes are the surviving fiber.

    Notes
    -----
    ``eq=False`` (identity equality / hash): the ``weight`` is an
    :class:`numpy.ndarray`, for which value-equality is element-wise (not a
    ``bool``) and hashing is undefined — the functional is identified by
    object identity, like :class:`~orpheus.transport.fields.cross_section_field.CrossSectionField`.
    """

    weight: NDArray
    axis: int = 0

    def evaluate(self, x: "Field | NDArray", /) -> np.ndarray:
        r"""Return :math:`\langle w, x\rangle` contracted over :attr:`axis`.

        .. math::

            \langle w, x\rangle \;=\; \sum_{k} w_k\, x_k
            \qquad\text{(over }\texttt{axis}\text{, }\texttt{keepdims=True}\text{)} .

        Accepts a typed :class:`~orpheus.numerics.field.Field` carrier (its
        ``.values`` is read) or a raw :class:`numpy.ndarray`. The reduction is
        the single numpy primitive ``(w * x).sum(axis, keepdims=True)`` — the
        same one the legacy ``RankOneOperator.apply`` performed inline, so a
        dyad built ``outer(v, InnerProductFunctional(w, axis))`` reproduces the
        old rank-1 action bit-for-bit.
        """
        w_arr = self.weight.values if isinstance(self.weight, Field) else np.asarray(self.weight)
        x_arr = x.values if isinstance(x, Field) else np.asarray(x)
        return (w_arr * x_arr).sum(axis=self.axis, keepdims=True)
