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
(``ProductionRateFunctional`` in :mod:`orpheus.transport`) conforms to
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
``ProductionRateFunctional``, which carries :math:`\nu\Sigma_f` as a
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

from typing import Protocol, TypeVar, runtime_checkable

from orpheus.numerics.vector import Vector

__all__ = ["Functional", "V_contra", "R_co"]

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
    ``orpheus.transport.production_rate_functional.ProductionRateFunctional``,
    which carries :math:`\nu\Sigma_f` and computes the per-cell fission
    emission density :math:`p(\vec r) = \sum_{g'} \nu\Sigma_{f,g'}\,\phi_{g'}`.
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
