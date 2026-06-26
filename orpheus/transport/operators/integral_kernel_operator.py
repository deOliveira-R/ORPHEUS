r"""The integral-kernel type of the transport operator algebra (§5.6 Kernel).

The grand-report §5.6 *suffix law* partitions the maps of the transport
algebra by **what they produce and how they reach across the domain**:

* an **Operator** maps a field to a field **locally** — its action at a
  point depends only on the field at that point. The canonical instance
  is the §5.7 multiplication operator
  :class:`~orpheus.transport.operators.multiplication_operator.MultiplicationOperator`
  (``M[f]`` = ``C = M[σ_t]``): pointwise multiply, **diagonal** in
  ``(cell, group, ordinate)``, no integration against a measure.
* a **Kernel** maps a field to a field **nonlocally** — its action
  integrates the field against a measure on one or more axes. The
  Boltzmann scattering and fission emission operators are the named
  instances: the fission rate at :math:`(\vec r, g)` reads the flux at
  *every* group :math:`g'` at :math:`\vec r` (a contraction over the
  group axis), and the :math:`P_\ell` anisotropic scattering source at
  ordinate :math:`\hat\Omega_n` reads the flux at *every* ordinate
  (a projection-then-reconstruction over the angular axis). This module
  names that category.
* a **Functional** maps a field to a **scalar** (or, fiberwise over
  space, a scalar-field) — the
  :class:`~orpheus.numerics.functional.Functional` surface (S5,
  ``evaluate(x) -> R``), disjoint from the operator surface.

The locality criterion (Frame 3)
=================================

The single discriminator between an **Operator** and a **Kernel** is
**locality**, in the sense of the cross-domain-attacker memo's Frame 3
(``coefficient_field_promotion_frames.md``):

* a :class:`~orpheus.transport.operators.multiplication_operator.MultiplicationOperator`
  is **local / diagonal** — :math:`(M[f]\,\psi)(\vec r,\hat\Omega,g)
  = f(\vec r, g)\,\psi(\vec r,\hat\Omega,g)`; the output at a point is a
  pointwise function of the input at that *same* point, so the operator
  has no integral structure to expose;
* an :class:`IntegralKernelOperator` is **nonlocal** — its action sums
  the field against a measure over at least one axis (groups for
  fission, ordinates for :math:`P_\ell` scattering), so there is a
  genuine *integral kernel* :math:`K` such that
  :math:`(A\,\psi)(x) = \int K(x, x')\,\psi(x')\,d\mu(x')`. The
  :attr:`kernel` member exposes that kernel as a typed, composable
  :class:`~orpheus.numerics.operator.LinearOperator`.

A Kernel REFINES LinearOperator (unlike the disjoint Functional)
================================================================

Crucially the Kernel category is a **refinement of**
:class:`~orpheus.numerics.operator.LinearOperator`, NOT disjoint from
it. A Kernel still maps a field to a field — it still has ``apply`` +
``capabilities`` (the matvec arms of fission / scattering are untouched;
the integral structure is the *reading* of the same action). This is the
asymmetry the type partition encodes:

* a Kernel **is** a LinearOperator (it has the full operator surface)
  AND **adds** the :attr:`kernel` member;
* but only **some** LinearOperators are Kernels — a local operator
  (:class:`~orpheus.transport.operators.multiplication_operator.MultiplicationOperator`,
  :class:`~orpheus.numerics.operator.IdentityOperator`) is a clean
  LinearOperator with **no** ``kernel`` and is therefore NOT an
  :class:`IntegralKernelOperator`.

Contrast the S5 :class:`~orpheus.numerics.functional.Functional`, which
shares **no** member with LinearOperator (it speaks ``evaluate``, not
``apply``): the Functional is a *sibling* of the operator category,
while the Kernel is a *refinement* of it. The intrinsic-property gate
(``tests/transport/test_integral_kernel_category.py``) verifies both
relationships: a Kernel is still a LinearOperator, but a kernel-less
LinearOperator (Identity / Multiplication) is not a Kernel, and a
Functional is neither.

What the ``kernel`` member returns
==================================

The :attr:`kernel` is the operator's **separable / composite integral
structure**, typed as the common
:class:`~orpheus.numerics.operator.LinearOperator` supertype (NOT a
forced single concrete return type, because the two named instances
expose structurally different kernels):

* **Fission** —
  :class:`~orpheus.numerics.operator.TensorProductOperator`
  ``RankOneOperator(χ, νΣf, axis=0) & IdentityOperator()`` — the rank-1
  ``χ ⊗ νΣf`` emission kernel: a group-axis outer product (contract
  ``νΣf·φ`` over groups, broadcast across groups by ``χ``).
* **Scattering** —
  :class:`~orpheus.numerics.operator.OperatorProduct`
  ``R ∘ Λ ∘ M`` — the anisotropic Legendre redistribution
  (``M`` projects ``ψ`` onto harmonic moments, ``Λ`` applies the per-ℓ
  group-transfer cross sections, ``R`` reconstructs the per-ordinate
  source via the addition theorem).

The union of the two kernel shapes is just
:class:`~orpheus.numerics.operator.LinearOperator`; the Protocol does
not force a narrower common type.

Why a structural (duck-typed) Protocol, not an ABC
==================================================

The same reasoning that makes
:class:`~orpheus.numerics.operator.LinearOperator`,
:class:`~orpheus.numerics.vector.Vector`, and the S5
:class:`~orpheus.numerics.functional.Functional` structural Protocols
applies here. The two named Kernels (``FissionOperator`` /
``ScatteringOperator``) are
:class:`~orpheus.numerics.operator.LinearOperator` subclasses built
with :func:`functools.singledispatchmethod` dispatch arms and
``@dataclass`` ; a nominal ABC base would collide with that MRO and
force a redundant inheritance edge. A structural Protocol lets each
operator participate simply by exposing ``kernel`` — and
``@runtime_checkable`` lets ``isinstance(op, IntegralKernelOperator)``
discriminate the category at runtime, checking the full LinearOperator
surface (``apply`` / ``capabilities`` / ``domain`` / ``codomain``) plus
the ``kernel`` refinement.

Layering note
=============

The §5.6 *Kernel* is **transport** vocabulary (the integral kernels of
the Boltzmann equation), shared across SN / CP / MoC, so it lives in
:mod:`orpheus.transport` (L2) alongside
:class:`~orpheus.transport.operators.multiplication_operator.MultiplicationOperator`
(S3b, the local Operator) and
:class:`~orpheus.transport.production_rate_functional.ProductionRateFunctional`
(S5, the Functional). The named Kernels
(:class:`~orpheus.transport.operators.fission.FissionOperator` /
:class:`~orpheus.transport.operators.scattering.ScatteringOperator`) ALSO live in
:mod:`orpheus.transport` (L2) and satisfy this Protocol from the same
layer — #261 relocated them out of the ``sn`` package (they are
cross-method reaction operators — every method scatters / fissions —
not SN-specific), realising the S6 carrier-agnostic-core plan.

References
----------

* Grand Report v3 §5.6 — the suffix law (Operator / Kernel / Functional
  / Projection / Reconstruction); §15 (the rank-1 fission ``χ ⊗ νΣf``
  and the ``Σ_ℓ Pℓ ⊗ Σ_{s,ℓ}`` scattering forms).
* ``.claude/agent-memory/cross-domain-attacker/coefficient_field_promotion_frames.md``
  — Frame 3 (the locality criterion separating the multiplication
  operator from the integral kernels).
* ``.claude/plans/issue_257_coefficient_field_promotion.md`` — S6.
"""

from __future__ import annotations

from typing import Protocol, TypeVar, runtime_checkable

from orpheus.numerics.operator import LinearOperator
from orpheus.numerics.vector import Vector

__all__ = ["IntegralKernelOperator"]

#: Carrier type variable for :class:`IntegralKernelOperator`, bound to
#: :class:`~orpheus.numerics.vector.Vector` and **invariant** — inherited
#: verbatim from :class:`~orpheus.numerics.operator.LinearOperator`
#: (``apply(x: V) -> V`` both consumes and returns the carrier, so ``V``
#: is invariant by dual use, the same as the operator endomorphism).
V = TypeVar("V", bound=Vector)


@runtime_checkable
class IntegralKernelOperator(LinearOperator[V], Protocol[V]):
    r"""Structural contract for a §5.6 *Kernel*: a NONLOCAL linear operator.

    An :class:`IntegralKernelOperator` is a
    :class:`~orpheus.numerics.operator.LinearOperator` whose action is
    **nonlocal** — it integrates the carrier field against a measure on
    one or more axes (groups for fission, ordinates for :math:`P_\ell`
    scattering) — and which exposes that integral structure through a
    single :attr:`kernel` property returning a
    :class:`~orpheus.numerics.operator.LinearOperator`.

    It is a **refinement of** LinearOperator (it inherits the full
    operator surface ``apply`` / ``capabilities`` / ``domain`` /
    ``codomain`` and adds ``kernel``), NOT disjoint from it — unlike the
    S5 :class:`~orpheus.numerics.functional.Functional`, which shares no
    member with LinearOperator. A LOCAL / diagonal operator (the §5.7
    :class:`~orpheus.transport.operators.multiplication_operator.MultiplicationOperator`,
    or :class:`~orpheus.numerics.operator.IdentityOperator`) is a clean
    LinearOperator with no ``kernel``, and so is NOT an
    :class:`IntegralKernelOperator` — the refinement is strict.

    Satisfied — without inheritance — by any object that carries the
    LinearOperator surface AND a ``kernel`` member. The canonical
    transport instances are
    :class:`orpheus.transport.operators.fission.FissionOperator` (kernel =
    ``χ ⊗ νΣf`` rank-1 TensorProductOperator) and
    :class:`orpheus.transport.operators.scattering.ScatteringOperator` (kernel =
    ``R ∘ Λ ∘ M`` OperatorProduct).
    """

    @property
    def kernel(self) -> LinearOperator:
        r"""The operator's integral-kernel structure, as a LinearOperator.

        The §5.6 nonlocal-action surface: the separable / composite
        integral kernel :math:`K` such that
        :math:`(A\,\psi)(x) = \int K(x, x')\,\psi(x')\,d\mu(x')`. It is
        the member that distinguishes a Kernel from a local
        :class:`~orpheus.transport.operators.multiplication_operator.MultiplicationOperator`
        (which has no integral structure to expose). Returned as the
        common :class:`~orpheus.numerics.operator.LinearOperator`
        supertype — a
        :class:`~orpheus.numerics.operator.TensorProductOperator` for
        the fission rank-1 ``χ ⊗ νΣf`` form, an
        :class:`~orpheus.numerics.operator.OperatorProduct` for the
        scattering ``R ∘ Λ ∘ M`` redistribution.
        """
        ...
