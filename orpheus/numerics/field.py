r"""L1 algebraic base of every typed transport field.

A :class:`Field` is the pair ``(values: ndarray, space: FunctionSpace)``
with closed same-CLASS, same-SPACE arithmetic. It is the codomain-side
sibling of :class:`~orpheus.numerics.operator.LinearOperator`: where
operators carry ``(domain, codomain)`` tags to enable algebraic
composition, Fields carry a single ``space`` tag that anchors the
inner product, dimensional label, and shape contract used by the
dunder algebra.

Why this lives at L1
====================

The Grand Report v3 (§5.5, §32.5) prescribes that every typed
transport field — ``AngularFlux``, ``ScalarFlux``, ``HarmonicMomentField``,
``BoundaryFaceFlux``, ``IsotropicSource``, ``PerOrdinateSource``,
``AngularSource``, ``AngularResidual`` — share a single algebraic
base. The pre-Depth-B codebase had six concrete classes each carrying
an identical hand-coded dunder skeleton plus a ``_validate_partner``
helper (see ``orpheus/sn/scalar_flux.py`` pre-D-D). The repetition
IS the architectural smell — Cardinal Rule 2 ("if you see shared
CONCEPTS in 2 places, stop and reconsider"). :class:`Field` is the
consolidation.

It lives at L1 because it knows nothing about transport: it is just
"values + space + algebra". The mesh-bound transport types (AngularFlux,
etc.) lift this ABC at L2 (``orpheus/transport/fields/...``) by adding
their domain-specific fields (e.g. ``mesh: SNMesh``, ``boundary``).

The three-layer dimensional enforcement story
=============================================

Dimensional consistency is enforced at three layers, each with a
different cost/coverage trade-off (see plan §7.5):

* **Layer 1 — class identity.** ``Field._check_partner`` rejects
  ``type(self) is not type(other)`` before any value comparison.
  This is the *primary* gate: even when units match (e.g. an
  ``AngularSource`` and an ``AngularResidual`` may both carry
  ``1/(cm²·s·sr·eV)``), the cross-class arithmetic raises by
  construction. Same units gives PERMISSION to add in linear
  algebra; it does not give MEANING. The "complex-frequency"
  analog: ``i ω − γ = s`` is the named complex-frequency
  construction; naive ``ω + γ`` is forbidden by convention.
* **Layer 2 — dimensional check at solver construction.**
  ``SourceIteration.__init__`` / ``KrylovAcceleration.__init__``
  do a single O(1) ``pint.Unit.dimensionality`` comparison per
  operator to verify the operator algebra is dimensionally sound
  before any iteration runs. Cost: microseconds per solver build.
  ALWAYS runs (both ``pytest`` and ``python -O``). See plan §5.
* **Layer 3 — defensive assert in dunders.** Inside
  :meth:`_check_partner` an ``assert self.space.units ==
  other.space.units`` catches the rare class/units misdesign
  (two instances of the same class whose spaces nonetheless
  carry inconsistent units — e.g. one in ``1/cm²``, one in
  ``1/m²``). Stripped in ``-O`` mode; defense in depth during
  development.

Together these layers make dimensional-mismatch bugs unrepresentable
without paying the cost of full pint-Quantity arithmetic on every
ndarray operation (which would be prohibitive).

Class identity for cross-class same-units operations
====================================================

When two distinct Field subclasses share a dimensional signature
(e.g. ``AngularResidual`` and ``AngularSource`` both have units
``1/(cm²·s·sr·eV)``), arithmetic between them is REQUIRED to go
through an explicit *named composition* — a factory method that
constructs the result with a definite physical interpretation. The
canonical example (planned for Issue #201):

.. code-block:: python

    # FORBIDDEN — cross-class arithmetic raises TypeError:
    iteration_residual = angular_residual - angular_source

    # REQUIRED — named composition with explicit physical meaning:
    iteration_residual = IterationResidual.from_balance(
        lhs=angular_residual, rhs=angular_source,
    )

The named-composition discipline IS what makes the Field algebra
sound under physical interpretation. ``coding-elegance`` Pattern 4
(illegal states unrepresentable) takes its strictest form here.

References
----------

* Grand Report v3 §5.5 (Field hierarchy), §32.5 (Field primitive
  spec). Located at ``.claude/plans/neutron_transport_grand_report_v3.md``.
* Depth-B plan §3.2 (the ABC spec) and §7.5 (three-layer enforcement
  story). Located at ``.claude/plans/depth_b_field_on_function_space.md``.
* ``coding-elegance`` Pattern 1 (read-as-the-math via dunder),
  Pattern 2 (single source of truth), Pattern 4 (illegal states
  unrepresentable).
"""

from __future__ import annotations

from abc import ABC
from dataclasses import dataclass, replace
from typing import TypeVar

import numpy as np
from numpy.typing import NDArray

from orpheus.numerics.space import FunctionSpace

__all__ = ["Field"]


T = TypeVar("T", bound="Field")


@dataclass(frozen=True, eq=False, kw_only=True)
class Field(ABC):
    r"""Algebraic base of every typed transport field.

    A :class:`Field` is the pair ``(values, space)`` with
    same-CLASS, same-SPACE closed arithmetic. Concrete subclasses
    add domain-specific fields (``mesh``, ``boundary``, etc.) on top
    of this ABC; the dunder algebra is inherited unchanged.

    Parameters
    ----------
    values : NDArray
        Field values; ``values.shape`` MUST equal ``space.shape``.
    space : FunctionSpace
        The function space this field lives on. Carries the inner-
        product metric, the dimensional label (``units``), and the
        shape contract.

    Notes
    -----
    The dataclass is ``frozen=True`` (immutability — values must not
    be mutated externally) and ``eq=False`` (auto-generated ``__eq__``
    on ndarray fields is ambiguous and triggers truthiness errors).
    Field equality is intentionally NOT defined — two fields are
    "the same" only by identity (``is``), not by content equality.
    Use :meth:`linf` / :meth:`l2` / :meth:`inner_product` for
    quantitative comparison.

    ``kw_only=True`` lets subclasses declare additional required
    fields without running afoul of the dataclass field-order
    constraint (mandatory fields after defaulted ones).
    """

    values: NDArray
    space: FunctionSpace

    def __post_init__(self) -> None:
        if type(self) is Field:
            raise TypeError(
                "Field is abstract; instantiate a concrete subclass "
                "(e.g., AngularFlux, ScalarFlux, ...) instead."
            )
        if self.values.shape != self.space.shape:
            raise ValueError(
                f"{type(self).__name__}: values.shape {self.values.shape!r} "
                f"does not match space.shape {self.space.shape!r}"
            )

    # ── Algebra ─────────────────────────────────────────────────────────

    def _check_partner(self, other: object) -> None:
        r"""Reject ill-formed binary operations.

        Three-layer dimensional enforcement (see module docstring):

        * **Layer 1** (always): ``type(self) is type(other)``. The
          primary gate — cross-class arithmetic is forbidden even
          when units match.
        * **Layer 2** is at solver-construction (not here).
        * **Layer 3** (stripped in ``-O``): defensive assert that
          ``space.units`` agree exactly. Tautological if the class
          hierarchy is well-designed; catches class/units misdesign
          where two instances of the same class carry different
          units (e.g. one in ``1/cm²``, one in ``1/m²``).
        """
        if type(self) is not type(other):
            raise TypeError(
                f"{type(self).__name__} arithmetic requires a same-class "
                f"partner; got {type(other).__name__}. Cross-class "
                f"arithmetic (even with matching units) requires an "
                f"explicit named composition (e.g. "
                f"IterationResidual.from_balance(lhs, rhs))."
            )
        if self.space != other.space:
            raise ValueError(
                f"{type(self).__name__} arithmetic requires equal space; "
                f"got {self.space!r} vs {other.space!r}"
            )
        # Layer 3 — defensive units assertion, stripped in -O mode.
        assert self.space.units == other.space.units, (
            f"internal: class identity matched but units mismatch "
            f"({self.space.units} vs {other.space.units}) — "
            f"class/units inconsistency in {type(self).__name__}"
        )

    def __add__(self: T, other: T) -> T:
        self._check_partner(other)
        return replace(self, values=self.values + other.values)

    def __sub__(self: T, other: T) -> T:
        self._check_partner(other)
        return replace(self, values=self.values - other.values)

    def __neg__(self: T) -> T:
        return replace(self, values=-self.values)

    def __mul__(self: T, scalar: float) -> T:
        return replace(self, values=self.values * float(scalar))

    def __rmul__(self: T, scalar: float) -> T:
        return self.__mul__(scalar)

    def __truediv__(self: T, scalar: float) -> T:
        return replace(self, values=self.values / float(scalar))

    # ── Diagnostics ─────────────────────────────────────────────────────

    @property
    def linf(self) -> float:
        r"""Return :math:`\lVert\text{values}\rVert_\infty`."""
        return float(np.abs(self.values).max())

    @property
    def l2(self) -> float:
        r"""Return the space-induced L² norm
        :math:`\sqrt{\langle x, x \rangle_{\text{space}}}`."""
        return float(self.space.norm(self.values))

    def inner_product(self: T, other: T) -> float:
        r"""Return the space-induced inner product
        :math:`\langle \text{self}, \text{other} \rangle_{\text{space}}`.

        Requires ``other`` to be same-class and same-space — the
        partner check is identical to addition's.
        """
        self._check_partner(other)
        return self.space.inner_product(self.values, other.values)

    def copy(self: T) -> T:
        r"""Return a deep copy carrying an owned ndarray."""
        return replace(self, values=self.values.copy())
