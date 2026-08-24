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
transport field — ``AngularFlux``, ``ScalarFlux``, ``HarmonicMomentFlux``,
``AngularBoundaryFlux``, ``ScalarSourceSink``, ``AngularSourceSink``,
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
  ``AngularSourceSink`` and an ``AngularResidual`` may both carry
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

Two mechanisms work together when distinct Field subclasses share a
dimensional signature (e.g. ``AngularResidual`` and ``AngularSourceSink``
both carry ``1/(cm³·s·sr)``):

* **The class gate forbids cross-class arithmetic** even when units
  match — ``angular_residual - angular_source`` RAISES. Same units grant
  permission to add in linear algebra; they do not grant meaning.
* **A role transition goes through a named composition.** The transport
  balance :math:`r = A\psi - q` differences two *same-class*
  ``AngularSourceSink`` operands (the operator output :math:`A\psi` and
  the external source :math:`q`), but the *result* is a different role —
  a residual (a defect), not a source. Bare same-class subtraction would
  typecheck yet MIS-TYPE the defect as ``AngularSourceSink``; the named
  composition :meth:`Field._from_balance` (exposed as each residual
  leaf's ``from_balance`` factory) makes the transition explicit and
  lands the result in the correct ``AngularResidual`` class:

  .. code-block:: python

      # same-class subtraction typechecks but MIS-TYPES the defect:
      wrong = operator_output - q_ext           # AngularSourceSink (!)

      # named composition — correctly typed as the residual role:
      residual = AngularResidual.from_balance(lhs=operator_output, rhs=q_ext)

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
from typing import TYPE_CHECKING, ClassVar, Self, TypeIs, TypeVar

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
    #: abstract bases leave it unset. Declaring it here makes the
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

    # ── Construction ─────────────────────────────────────────────────────

    @classmethod
    def zeros(cls: type[T], space: FunctionSpace, **fields: object) -> T:
        r"""Allocate a zero-filled field of this class on ``space``.

        The single shared zero-allocation primitive (B.5.A): ``values`` is
        ``np.zeros(space.shape)``; any subclass-specific dataclass fields
        (``L``, ``spatial_moments``, ...) pass through ``**fields``. Since
        CS4b S5 this IS the allocator call sites spell, on the carrier's
        cached space mints (the mesh-keyed sugar tier retired; the moment
        family's keyed ``zeros_for_mesh_and_L`` still delegates here until
        its S6 re-home) — the zero-construction lives in exactly one place
        (``coding-elegance`` Pattern 2).
        """
        return cls(values=np.zeros(space.shape), space=space, **fields)  # type: ignore[call-arg]

    # ── Algebra ─────────────────────────────────────────────────────────

    def _is_same_class(self, other: "Field") -> TypeIs[Self]:
        r"""The static face of the Layer-1 gate: ``type(other) is type(self)``.

        A ``TypeIs`` so that the one runtime fact the gate establishes —
        the partner is *exactly* this class — is available to the type
        checker wherever the gate has run (parse, don't validate).
        """
        return type(other) is type(self)

    def _check_partner(self, other: "Field") -> Self:
        r"""Reject an ill-formed binary partner; return it proven as ``Self``.

        Dimensional enforcement under View-G (see module docstring):

        * **Layer 1** (always, here): ``type(self) is type(other)``.
          The runtime gate — cross-class arithmetic is forbidden even
          when units match. Because each role-leaf's units are a class
          constant ``UNITS``, class identity *is* units identity.
        * **Layer 2** (operator construction, not here): the operator
          unit-gain dimensional check (issue #208).

        The return is the parse-don't-validate face: subclass overrides
        chain ``partner = super()._check_partner(other)`` and receive the
        partner re-typed as ``Self``, so their own guards (mesh identity,
        ``L`` match, layout) need no narrowing ceremony.
        """
        if not self._is_same_class(other):
            raise TypeError(
                f"{type(self).__name__} [{_unit_label(self)}] arithmetic "
                f"requires a same-class partner; got {type(other).__name__} "
                f"[{_unit_label(other)}]. Cross-class arithmetic is forbidden "
                f"even when the units match — same units grant permission to "
                f"add in linear algebra, not meaning. Use an explicit named "
                f"composition (e.g. AngularResidual.from_balance(lhs, rhs))."
            )
        if self.space != other.space:
            raise ValueError(
                f"{type(self).__name__} arithmetic requires equal space; "
                f"got {self.space!r} vs {other.space!r}"
            )
        return other

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

    # ── Named composition ────────────────────────────────────────────────

    @classmethod
    def _from_balance(cls: type[T], lhs: "Field", rhs: "Field") -> T:
        r"""Shared engine for the residual leaves' ``from_balance`` factory.

        The ONE sanctioned way to cross the Layer-1 class gate: build a
        field of class ``cls`` (always a *residual* role leaf) as the signed
        difference ``lhs − rhs`` of two same-class operands whose units match
        ``cls``. This is the dual of :meth:`__sub__` — ``__sub__`` keeps the
        class fixed; ``_from_balance`` permits a *class transition* (e.g.
        ``AngularSourceSink`` operands → ``AngularResidual`` result) precisely
        because the transition is dimensionally sound and explicitly named.

        Each residual leaf exposes this as its public ``from_balance``
        classmethod — its bespoke construction story (a residual is the
        *defect of a balance*, never built from thin air), parallel to
        ``AngularSourceSink.from_isotropic``. So
        the public API reads ``AngularResidual.from_balance(lhs=Aψ, rhs=q)``;
        the leaves delegate here so the check-and-reconstruct logic lives in
        exactly one place (``coding-elegance`` Pattern 2).

        Three guards, in order:

        1. **Same-class operands.** ``type(lhs) is type(rhs)`` — the two
           sides of one balance are the same kind before differencing.
        2. **Dimensional soundness.** ``lhs.UNITS == rhs.UNITS == cls.UNITS``,
           compared by EXACT unit equality (``sr``-sensitive — the ERR-039
           guard; see :mod:`orpheus.numerics.units`), so a residual cannot be
           formed from operands whose signature differs from its own.
        3. **Same space.** delegated to the operands' :meth:`_check_partner`
           (its class branch is already satisfied, so the space-content
           check fires).

        The result rides the operands' space directly (CS4b S4 — the Q1
        flip): since S2a every role leaf of one family on one carrier
        shares the SAME carrier-cached space (role is CLASS identity,
        never space identity), so ``lhs.space`` IS the residual's own
        space, and the construction is the plain L1 constructor —
        ``cls(values=lhs.values − rhs.values, space=lhs.space)``.
        The pre-S4 spelling routed through ``cls.from_mesh(…, lhs.mesh)``,
        an L2 detour whose two ``type: ignore`` died with the mesh field.

        **The recorded choice this flip creates (Q1, ruled 2026-08-22):**
        a ``MomentResidual`` is now one small override away — the moment
        family's extra fields (``L``, ``spatial_moments``) thread from the
        operand (``L=lhs.L, spatial_moments=lhs.spatial_moments``), no
        3-arg-factory blocker remains. It stays deliberately UNMINTED
        because the residual path's contract is per-ordinate (the
        operators consume reconstructed angular flux, never the moment
        iterate) — mint it WITH its first consumer, not before. The LD
        certificate un-skip this may also enable is #402.
        """
        if type(lhs) is not type(rhs):
            raise TypeError(
                f"{cls.__name__}.from_balance requires two same-class operands "
                f"(the two sides of one balance, differenced before the role "
                f"transition); got {type(lhs).__name__} and "
                f"{type(rhs).__name__}."
            )
        if not (lhs.UNITS == rhs.UNITS == cls.UNITS):
            raise TypeError(
                f"{cls.__name__}.from_balance: operand units "
                f"[{lhs.UNITS:~P}] must match the residual's units "
                f"[{cls.UNITS:~P}] (exact, sr-sensitive). A balance of two "
                f"{lhs.UNITS:~P} quantities yields a {cls.UNITS:~P} residual "
                f"only when the signatures agree."
            )
        lhs._check_partner(rhs)  # same space content (class branch already ok)
        return cls(values=lhs.values - rhs.values, space=lhs.space)

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

    def _index_tuples(self, positions: "np.ndarray") -> list[tuple[int, ...]]:
        r"""Decode flat positions into this leaf's ``values``-layout index
        tuples — the single spelling of the layout decode shared by every
        index-reporting diagnostic (``where_largest``, ``cone_violations``)."""
        shape = np.asarray(self.values).shape
        return [tuple(int(i) for i in np.unravel_index(p, shape)) for p in positions]

    @property
    def principal_bulk_leaf(self) -> "Field":
        r"""The leaf whose space norm carries this iterate's convergence
        diagnostics — for a leaf field, itself.

        The convention every carrier answers for itself (campaign 1 CS3-R,
        relocated from a structural walk in ``numerics.iteration``): a
        composite names its interior; a coupled block vector names its
        first system's; a leaf IS its own bulk. The iteration layer reads
        this property and nothing else about the iterate's anatomy.
        """
        return self

    def where_largest(self, k: int = 1) -> list[tuple[int, ...]]:
        r"""The ``k`` index tuples with the largest :math:`|\text{values}|`.

        The per-entry magnitude map: WHICH entries dominate — for an iterate
        increment, the cells / groups / ordinates that are not converging
        (pole-cell resonance, material-interface slow modes, a lagging
        group); for a residual, where the equation defect concentrates.
        Indices are into this leaf's ``values`` layout (e.g. ``(n, g, ix)``
        for an angular field), largest first.

        Promoted from the retired displacement diagnostics surface
        (campaign 1 CS3, 2026-08-19): the map needs only ``values``, so it
        is a property of ANY field, not of difference-ness.

        Parameters
        ----------
        k : int
            How many of the largest-magnitude entries to return (default 1).
        """
        flat = np.abs(np.asarray(self.values)).ravel()
        k = max(1, min(int(k), flat.size))
        top = np.argpartition(flat, -k)[-k:]
        top = top[np.argsort(flat[top])[::-1]]  # largest first
        return self._index_tuples(top)

    def cone_violations(self, k: int | None = None) -> list[tuple[int, ...]]:
        r"""Index tuples of entries OUTSIDE the positive cone ``K = {v ≥ 0}``.

        The element-level cone-membership predicate of the campaign-1 CS3
        ruling (flux lives in the positive cone :math:`K \subset V` of an
        ordered vector space; membership is a PREDICATE of the element,
        never a constructor invariant — diamond difference does not
        preserve K, so a ψ≥0 type would refuse production output).
        Emptiness IS membership; the returned structure makes the
        predicate's own correctness assertable (``vv`` anti-#14 — a bool
        could not say WHERE).

        An entry violates iff ``not (value >= 0.0)`` — so ``-0.0`` is a
        member (IEEE: ``-0.0 >= 0.0``), and a non-finite ``nan`` is a
        violation (an unordered entry is not in any cone), while ``+inf``
        is admitted. Violations are ordered most-negative first (``nan``
        entries last, at the magnitude ``nan`` sorts to); ``k`` caps how
        many are returned — the emptiness answer is exact under any cap.

        **What a green reading does NOT claim** (say it here so no audit
        reads it otherwise): production does not KEEP fields in K — there
        is no fixup, no clipping, no warning on the SN spatial path; a
        violating solve is not refused; DD is not repaired; and
        cone-PRESERVATION is the realization's own flag
        (``DiscretizationScheme.is_positivity_preserving``), not this
        observer's.

        **The basis-kind consult (campaign 1, CS1 step 4).** The
        per-component sign test is only meaningful on a COORDINATE cone,
        and whether the space has one is the space's own structural
        answer (:attr:`FunctionSpace.has_coordinate_cone`, three-valued):

        * ``False`` (axis-built, any MODAL factor) ⟹ **REFUSE** with a
          typed error — components are expansion coefficients, a positive
          function may have negative coefficients, so answering would
          manufacture violations out of a basis choice;
        * ``True`` (axis-built, all NODAL) ⟹ answer, same arithmetic as
          the legacy path;
        * ``None`` (legacy space, ``axes is None``) ⟹ the pre-CS1
          behavior unchanged, deliberately — the question cannot be
          answered structurally, and refusing would fire on every legacy
          space in the tree.
        """
        if self.space.has_coordinate_cone is False:
            raise ValueError(
                f"cone_violations is meaningless on {self.space.name!r}: "
                f"the space carries a MODAL factor, so its components are "
                f"expansion COEFFICIENTS — a positive function may have "
                f"negative coefficients, and a per-component sign test "
                f"would manufacture violations out of a basis choice. "
                f"Evaluate the field in a nodal realization first."
            )
        flat = np.asarray(self.values).ravel()
        bad = np.flatnonzero(~(flat >= 0.0))
        if bad.size == 0:
            return []
        bad = bad[np.argsort(flat[bad])]  # most-negative first
        if k is not None:
            bad = bad[: max(1, int(k))]
        return self._index_tuples(bad)

    def copy(self: T) -> T:
        r"""Return a deep copy carrying an owned ndarray."""
        return replace(self, values=self.values.copy())
