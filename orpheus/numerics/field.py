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
``BoundaryFlux``, ``ScalarSource``, ``AngularSource``,
``ScalarResidual``, ``AngularResidual`` — share a single algebraic
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

Dimensional enforcement under View-G
====================================

Units do NOT live on the :class:`~orpheus.numerics.space.FunctionSpace`
(the "View-G" decision, issues #205 / #207): a space is pure geometry
and is role-/dimension-agnostic. Units are a property of the *quantity*
(the field) and the *map* (the operator), enforced at two layers:

* **Layer 1 — class identity.** ``Field._check_partner`` rejects
  ``type(self) is not type(other)`` before any value comparison.
  This is the runtime gate: even when units match (e.g. an
  ``AngularSource`` and an ``AngularResidual`` may both carry
  ``1/(cm³·s·sr·eV)``), the cross-class arithmetic raises by
  construction. Same units give PERMISSION to add in linear
  algebra; they do not give MEANING. The "complex-frequency"
  analog: ``i ω − γ = s`` is the named complex-frequency
  construction; naive ``ω + γ`` is forbidden by convention. Because
  each role-leaf carries its units as a class constant ``UNITS``
  (Phase B of the field-vocabulary plan), class identity *is* units
  identity — no per-instance units check is needed, and the old
  space-level ``units`` assert retired together with the field.
* **Layer 2 — dimensional check at operator construction.** The
  dimensional soundness of an operator algebra (the ``L``, ``C``,
  ``S``, ``F`` unit-gains composing to a consistent residual) is
  checked once when the composite operator is assembled, via the
  operator's unit-gain — O(1) per build, runs under ``python -O``.
  This is the operator-side counterpart to the field-side ``UNITS``
  and lands with the operator typing (issue #208).

Together these make dimensional-mismatch bugs unrepresentable without
paying the cost of full pint-Quantity arithmetic on every ndarray
operation (which would be prohibitive).

Class identity for cross-class same-units operations
====================================================

When two distinct Field subclasses share a dimensional signature
(e.g. ``AngularResidual`` and ``AngularSource`` both have units
``1/(cm³·s·sr·eV)``), arithmetic between them is REQUIRED to go
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
from typing import TYPE_CHECKING, ClassVar, TypeVar

import numpy as np
from numpy.typing import NDArray

from orpheus.numerics.space import FunctionSpace

if TYPE_CHECKING:
    from orpheus.numerics.units import Unit

__all__ = ["Field"]


T = TypeVar("T", bound="Field")


def _unit_label(obj: object) -> str:
    r"""Format an object's ``UNITS`` for a diagnostic message.

    Reads the class constant ``UNITS`` (a :class:`pint.Unit`, set on every
    role leaf — see :mod:`orpheus.numerics.units`) and renders it short +
    pretty (``"1/cm³/s/sr"``). Returns ``"<no units>"`` for anything without
    a ``UNITS`` constant (the abstract bases, or a non-:class:`Field`
    operand). Diagnostics-only — never called on the arithmetic hot path.
    """
    units = getattr(type(obj), "UNITS", None)
    if units is None:
        return "<no units>"
    return f"{units:~P}"


@dataclass(frozen=True, eq=False, kw_only=True, repr=False)
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
        product metric and the shape contract. Under View-G it does
        NOT carry units — those are the field's :attr:`UNITS` constant.

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

    #: Dimensional identity of the quantity (View-G). Every CONCRETE role
    #: leaf MUST set this to one of the signatures in
    #: :mod:`orpheus.numerics.units` (e.g. ``ANGULAR_FLUX_UNITS``); the
    #: abstract bases leave it unset. This is the units-side companion of
    #: the per-leaf ``_SPACE_NAME`` contract — declaring it here makes the
    #: obligation typed and visible: a leaf that forgets ``UNITS`` raises on
    #: ``.UNITS`` access (so #208's operator unit-gain check fails loudly
    #: rather than silently inheriting ``None``). Class identity *is* units
    #: identity (one ``UNITS`` per leaf class), so the class-identity
    #: arithmetic gate doubles as a units gate WITHOUT ever reading
    #: ``UNITS`` on the hot path.
    UNITS: ClassVar[Unit]

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

    def __repr__(self) -> str:
        r"""Concise, units-aware repr — class, shape, space name, units.

        Replaces the dataclass auto-repr (which dumps the full ``values``
        ndarray and the ``space`` / ``mesh`` objects). Every concrete leaf
        sets ``repr=False`` on its ``@dataclass`` so it inherits this. The
        ``units=`` field makes the role-leaf's dimensional identity visible
        at a glance (and in tracebacks).
        """
        return (
            f"{type(self).__name__}(shape={self.values.shape}, "
            f"space={self.space.name!r}, units={_unit_label(self)})"
        )

    # ── Algebra ─────────────────────────────────────────────────────────

    def _check_partner(self, other: object) -> None:
        r"""Reject ill-formed binary operations.

        Dimensional enforcement under View-G (see module docstring):

        * **Layer 1** (always, here): ``type(self) is type(other)``.
          The runtime gate — cross-class arithmetic is forbidden even
          when units match. Because each role-leaf's units are a class
          constant ``UNITS``, class identity *is* units identity.
        * **Layer 2** (operator construction, not here): the operator
          unit-gain dimensional check (issue #208).
        """
        if type(self) is not type(other):
            raise TypeError(
                f"{type(self).__name__} [{_unit_label(self)}] arithmetic "
                f"requires a same-class partner; got {type(other).__name__} "
                f"[{_unit_label(other)}]. Cross-class arithmetic is forbidden "
                f"even when the units match — same units grant permission to "
                f"add in linear algebra, not meaning. Use an explicit named "
                f"composition (e.g. IterationResidual.from_balance(lhs, rhs))."
            )
        if self.space != other.space:
            raise ValueError(
                f"{type(self).__name__} arithmetic requires equal space; "
                f"got {self.space!r} vs {other.space!r}"
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
