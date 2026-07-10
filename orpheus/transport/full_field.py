r"""The composite carrier: a generic ``interior ⊕ boundary`` block field.

L2 (transport, method-agnostic). This module holds TWO things:

* :class:`Composite` — the **generic** two-block composite
  ``Composite[Interior, Boundary]``: an ``interior`` field paired with its
  ``boundary`` partner, carrying the vector-space algebra ONCE. It is the
  structural object the transport operator algebra
  :math:`(L + C - S - F - B)\,\psi = q` acts on — every operator leaf maps a
  composite to a composite (the inner role may change, flux → source, but the
  carrier type does not). The name is **structural, not domain-role**: a
  domain meaning arises from the *specialization* — ``SN`` is
  ``Composite[AngularFlux, AngularBoundaryFlux]``, diffusion / CP are
  ``Composite[ScalarFlux, ScalarBoundaryFlux]``, MoC its own leaves.
* :class:`FullField` — the **SN specialization**
  ``Composite[BulkField, BoundaryField]`` plus the OPTIONAL curvilinear
  starting-direction (ψ½) third block. This third block is transport-SN
  -specific; it is a *temporary* extension the coupled-block campaign's Phase B
  evicts into its own System-B composite (after which ``FullField`` collapses
  to a pure ``Composite``).

Why the interior + boundary split is the right L2 abstraction
=============================================================

The pre-D-H :class:`~orpheus.sn.angular_flux.AngularFlux` conflated volumetric
flux values with boundary trace values. Per Cardinal Rule 2 (shared concepts →
shared abstraction), the interior + boundary split is NOT SN-specific — every
transport method has the same pair:

* **SN**: interior =
  :class:`~orpheus.transport.fields.angular_flux.AngularFlux` on
  ``(N, ng, *spatial)``; boundary =
  :class:`~orpheus.transport.fields.angular_boundary_flux.AngularBoundaryFlux`
  on the flat face-trace layout.
* **CP / MoC / diffusion**: their own interior + boundary leaves
  (diffusion / CP already build ``Composite[ScalarFlux, ScalarBoundaryFlux]``).

The optional starting-direction block (#282 route (a), 2.5d)
============================================================

On a mesh whose Morel–Montry thread genuinely consumes independent
starting-direction state (the R12a predicate — the 1-D sphere; see
:mod:`orpheus.numerics.spaces.radial_characteristic_space`), the SN
:class:`FullField` carries a THIRD block: the per-level half-angle flux ψ½ as a
typed :class:`~orpheus.transport.fields._bases.RadialCharacteristicField` leaf.
This dissolves the #282 back edge — the seed stops being a lagged
solver-internal estimate and becomes state the operator reads and the solve
produces. ``radial_characteristic`` is ``None`` exactly when the mesh carries no
seed level (Cartesian; every production cylinder rule), presence is resolved
MESH-side by :meth:`FullField.zeros`, and mixed-presence arithmetic raises —
three spellings of the same illegal-states-unrepresentable discipline. The ψ½
block lives on the SN subclass, NOT on the generic :class:`Composite`, because
it is a curvilinear-SN augmentation, not part of every method's interior ⊕
boundary structure.

The cofree-comonad framing (the #217 split)
===========================================

:class:`~orpheus.transport.timed_full_field.TimedFullField` is the **cofree
comonad** ``Cofree(FullField, depth=d)`` over the base: it pairs the *current*
timeless frame with a rotating history buffer of prior timeless frames. Only the
iteration / time-stepping drivers see the comonad (the history); the operator
algebra is blind to it — it reads a timeless frame in and writes a timeless
frame out.

This is why the split is structurally FORCED rather than aesthetic:

* A **static source** ``q = q_int ⊕ q_∂`` has no time — it never ``advance``\ s.
  Typing it as a history-bearing ``TimedFullField`` hands it verbs (``advance`` /
  ``at_lag``) it must never use ("a type error of altitude").
* An **iterating flux state** ``ψ^n`` DOES advance through iterations / time
  steps — it is the comonad, ``TimedFullField``.

So the operator-algebra carrier is the timeless :class:`FullField`; the
driver-level iterate is the timed :class:`TimedFullField`. This module holds the
algebra ONCE (DRY) on :class:`Composite` — every subclass inherits it.

Algebra contract
================

Same-class arithmetic propagates to every block:

.. code-block:: python

    a + b = Composite(interior=a.interior + b.interior,
                      boundary=a.boundary + b.boundary)

The six vector-space dunders (``+``, ``-``, unary ``-``, scalar ``*``,
``scalar *``, ``/ scalar``) live ONCE on :class:`Composite` and route through two
small per-shape hooks — :meth:`Composite._map_binary` (elementwise over two
operands' blocks) and :meth:`Composite._map_unary` (elementwise over one
operand's blocks) — then rebuild via the polymorphic :meth:`Composite._recombine`
hook. A subclass carrying an extra block (:class:`FullField`'s ψ½) overrides ONLY
those three hooks to thread the extra block; the dunders themselves are never
duplicated. :class:`~orpheus.transport.timed_full_field.TimedFullField` overrides
:meth:`_recombine` alone to rebuild a ``TimedFullField`` with empty history
(#217: "algebra results carry empty history"). One definition of the algebra, the
correct concrete return type for each subclass.

Cross-class arithmetic is rejected at two layers: :meth:`Composite._check_partner`
rejects a partner that is not a :class:`Composite`; the member-level leaf dunders
enforce role / units / mesh / space (the #208 affine torsor gate ``flux + flux →
TypeError``, the cross-mesh guard, the leaf-type match ``AngularFlux + ScalarFlux
→ TypeError``) by delegation, so the composite does NOT pre-check member types
(that would block the legitimate composite torsor ``flux + displacement → flux``).

Grep signal
===========

``Composite`` — the generic *structural* interior ⊕ boundary carrier.
``FullField`` — the SN specialization (the *full* SN domain, *timeless*), plus
the optional ψ½ block. Its history-bearing subclass keeps the strong three-token
grep signal ``TimedFullField`` (Timed + Full + Field).

References
==========

* GH **issue #217** — the timeless-``FullField`` extraction (the composite source
  is the first timeless consumer).
* ``.claude/plans/coupled_block_operator_campaign.md`` — the ``Composite``
  generalization (Phase A2) + the ψ½-eviction (Phase B).
* Grand Report v3 §5.5 (Field hierarchy), §5.3 (``DirectSumSpace`` /
  :class:`~orpheus.numerics.spaces.full_field_space.FullFieldSpace`).
* ``coding-elegance`` Pattern 2 (the algebra lives ONCE in the base via the
  recombine + block-map hooks) + Pattern 5 (``Composite`` is the right primitive;
  ``FullField`` composes it with the ψ½ block, ``TimedFullField`` with history).
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from typing import TYPE_CHECKING, Callable, Generic, Self, TypeVar

import numpy as np

from orpheus.transport.fields._bases import (
    BulkField,
    BoundaryField,
    RadialCharacteristicField,
)

if TYPE_CHECKING:
    from numpy.typing import NDArray

    from orpheus.transport.mesh.material_mesh import MaterialMesh


__all__ = ["Composite", "FullField"]

#: The interior (volumetric) leaf type — any
#: :class:`~orpheus.transport.fields._bases.BulkField` subclass
#: (``AngularFlux`` for SN, ``ScalarFlux`` for diffusion / CP, ``RayField`` for
#: MoC). The generic parameter makes ``Composite[AngularFlux, ...].interior``
#: read as the PRECISE leaf type, not the ``BulkField`` base.
Interior = TypeVar("Interior", bound=BulkField)
#: The boundary (trace) leaf type — any
#: :class:`~orpheus.transport.fields._bases.BoundaryField` subclass.
Boundary = TypeVar("Boundary", bound=BoundaryField)
#: A composite flavor — the template type of :meth:`Composite.from_flat` (its
#: return follows the template, so the reconstruction preserves the concrete
#: subclass: a ``TimedFullField`` template yields a ``TimedFullField``).
CompositeT = TypeVar("CompositeT", bound="Composite")


@dataclass(frozen=True, kw_only=True)
class Composite(Generic[Interior, Boundary]):
    r"""The generic ``interior ⊕ boundary`` composite carrier.

    A structural two-block composite: an ``interior`` volumetric field paired
    with its ``boundary`` trace partner. Generic over the two leaf types, so the
    specialization carries the domain meaning
    (``Composite[AngularFlux, AngularBoundaryFlux]`` = an SN frame,
    ``Composite[ScalarFlux, ScalarBoundaryFlux]`` = a diffusion / CP frame),
    while THIS type is method-agnostic. Holds the vector-space algebra ONCE (the
    six dunders route through the :meth:`_map_binary` / :meth:`_map_unary` /
    :meth:`_recombine` hooks a subclass overrides to add extra blocks).

    Parameters
    ----------
    interior : Interior
        The volumetric field. Any
        :class:`~orpheus.transport.fields._bases.BulkField` subclass.
    boundary : Boundary
        The boundary partner field on the trace of ``interior``'s domain. Any
        :class:`~orpheus.transport.fields._bases.BoundaryField` subclass.

    Notes
    -----
    NOT a :class:`~orpheus.numerics.field.Field` subclass at the typed-class
    level — its natural backing is a
    :class:`~orpheus.numerics.spaces.full_field_space.FullFieldSpace`
    (direct-sum) Field. Ships as a structured composite with delegate-style
    dunders that propagate to every block.
    """

    interior: Interior
    boundary: Boundary

    # ── Construction ─────────────────────────────────────────────────

    @classmethod
    def zeros(
        cls,
        *,
        interior: "type[Interior]",
        boundary: "type[Boundary]",
        mesh: "object",
    ) -> "Self":
        r"""Allocate a zero composite from the interior + boundary leaf TYPES.

        Generic over the method's leaf types: the caller passes the interior and
        boundary :class:`~orpheus.numerics.field.Field` *subclasses* (SN passes
        :class:`~orpheus.transport.fields.angular_flux.AngularFlux` /
        :class:`~orpheus.transport.fields.angular_boundary_flux.AngularBoundaryFlux`;
        diffusion / CP pass their scalar leaves), and each is zero-allocated on
        ``mesh`` via its own :meth:`zeros_on`. This keeps the cross-method-generic
        container free of any hard-wired leaf type.
        """
        return cls(
            interior=interior.zeros_on(mesh),  # type: ignore[attr-defined]
            boundary=boundary.zeros_on(mesh),  # type: ignore[attr-defined]
        )

    # ── Construction validation ──────────────────────────────────────

    def __post_init__(self) -> None:
        if not isinstance(self.interior, BulkField):
            raise TypeError(
                f"{type(self).__name__}: bulk must be a BulkField; got "
                f"{type(self.interior).__name__}"
            )
        if not isinstance(self.boundary, BoundaryField):
            raise TypeError(
                f"{type(self).__name__}: boundary must be a BoundaryField "
                f"(an AngularBoundaryField / ScalarBoundaryField family leaf); "
                f"got {type(self.boundary).__name__}"
            )
        # Mesh-identity check (where both members carry a ``mesh`` attribute —
        # the cross-method generic contract). For SN both AngularFlux and
        # AngularBoundaryFlux carry ``mesh: SNMesh``; other methods follow the
        # same convention.
        bulk_mesh = getattr(self.interior, "mesh", None)
        boundary_mesh = getattr(self.boundary, "mesh", None)
        if bulk_mesh is not None and boundary_mesh is not None:
            if bulk_mesh is not boundary_mesh:
                raise ValueError(
                    f"{type(self).__name__}: bulk and boundary must share "
                    "mesh identity (both bound to the same mesh instance); "
                    f"got bulk.mesh={bulk_mesh!r}, "
                    f"boundary.mesh={boundary_mesh!r}"
                )

    # ── The composite's single mesh ───────────────────────────────────

    @property
    def mesh(self) -> "MaterialMesh":
        r"""The one mesh both leaves are bound to — the ``__post_init__``
        mesh-identity invariant made readable.

        Read off the BOUNDARY leaf (either leaf works — the invariant guarantees
        identity). The static type is the method-agnostic
        :class:`~orpheus.transport.mesh.material_mesh.MaterialMesh` base: each
        trace family narrows its OWN ``mesh`` declaration
        (:class:`~orpheus.transport.fields._bases.AngularBoundaryField` →
        ``SNMesh``), so an SN consumer that needs quadrature / sweep data reads it
        off the SNMesh its operator was constructed with (the #226 F2 pattern) or
        off the family-typed boundary leaf — never through this method-generic
        surface.
        """
        return self.boundary.mesh

    # ── Polymorphic recombine + block-map hooks (Pattern 2) ───────────

    def _recombine(
        self, *, interior: "BulkField", boundary: "BoundaryField",
    ) -> "Self":
        r"""Rebuild a composite of the SAME concrete type from recombined blocks.

        The polymorphic hook the block-map methods route their result through.
        The base spelling is ``replace(self, ...)`` — provably ``Self`` (and
        ``replace`` re-runs ``__post_init__`` so the block invariants re-fire for
        free). A subclass carrying an extra block (:class:`FullField`'s ψ½)
        OVERRIDES this to accept + thread it;
        :class:`~orpheus.transport.timed_full_field.TimedFullField` overrides it to
        rebuild a ``TimedFullField`` with an EMPTY history (#217).
        """
        return replace(self, interior=interior, boundary=boundary)  # type: ignore[arg-type]

    def _map_binary(
        self, other: "Composite", op: "Callable[[object, object], object]",
    ) -> "Self":
        r"""Apply a binary elementwise ``op`` across the blocks of ``self`` and
        ``other``, then recombine.

        The single place a two-operand dunder (``+`` / ``-``) spells its block
        propagation. A subclass with an extra block overrides this to also
        combine that block (:class:`FullField` threads ψ½ via
        :meth:`FullField._combine_radial_characteristic`).
        """
        return self._recombine(
            interior=op(self.interior, other.interior),  # type: ignore[arg-type]
            boundary=op(self.boundary, other.boundary),  # type: ignore[arg-type]
        )

    def _map_unary(
        self, op: "Callable[[object], object]",
    ) -> "Self":
        r"""Apply a unary elementwise ``op`` across the blocks of ``self``, then
        recombine.

        The single place a one-operand transform (unary ``-``, scalar ``*`` / ``/``,
        :meth:`copy`) spells its block propagation. A subclass with an extra block
        overrides this to also transform that block.
        """
        return self._recombine(
            interior=op(self.interior),  # type: ignore[arg-type]
            boundary=op(self.boundary),  # type: ignore[arg-type]
        )

    # ── Algebra (propagates to every block via the map hooks) ─────────

    def _check_partner(self, other: object) -> None:
        r"""Reject a partner that is not a :class:`Composite`.

        Layer 1 at the CONTAINER level only. The member-level leaf algebra (the
        affine gate ``flux + flux → TypeError``, the torsor ``flux + displacement
        → flux``, the displacement mint ``flux − flux → displacement``, the
        leaf-type / cross-mesh rejection) is the SINGLE SOURCE OF TRUTH on the
        leaves — ``__add__`` / ``__sub__`` delegate to ``self.interior ±
        other.interior`` and ``self.boundary ± other.boundary``, where the leaf
        dunders enforce role-correctness (#208). Pre-checking member types here
        would BLOCK the legitimate composite torsor (``flux`` interior +
        ``displacement`` interior), so the redundant member pre-checks are
        intentionally absent.

        The accepted partner is ANY :class:`Composite` flavor — a timeless base or
        a timed subclass. This is load-bearing for the time-derivative stencil
        ``state.at_lag(0) - state.at_lag(1)``: the current frame is a
        :class:`~orpheus.transport.timed_full_field.TimedFullField` while a
        historical frame is a timeless :class:`FullField` snapshot, and the two
        must subtract. The CONCRETE result type is governed by :meth:`_recombine`,
        so accepting a base partner does not weaken the "algebra results carry
        empty history" guarantee. ``state + 42`` (a non-``Composite`` partner) is
        still rejected here; the leaf gate then rejects any leaf/mesh/units
        mismatch (including a cross-method ``AngularFlux + ScalarFlux``) by
        delegation.
        """
        if not isinstance(other, Composite):
            raise TypeError(
                f"{type(self).__name__} arithmetic requires a same-class "
                f"partner; got {type(other).__name__}."
            )

    # ``other`` is deliberately the BASE ``Composite`` (not ``Self``): the partner
    # rule is "any Composite flavor" (see ``_check_partner`` — load-bearing for the
    # timed − timeless time-derivative stencil), and the RESULT flavor is governed
    # by ``self``'s ``_recombine`` hook alone.
    def __add__(self, other: "Composite") -> "Self":
        self._check_partner(other)
        return self._map_binary(other, lambda a, b: a + b)  # type: ignore[operator]

    def __sub__(self, other: "Composite") -> "Self":
        self._check_partner(other)
        return self._map_binary(other, lambda a, b: a - b)  # type: ignore[operator]

    def __neg__(self) -> "Self":
        return self._map_unary(lambda a: -a)  # type: ignore[operator]

    def __mul__(self, scalar: float) -> "Self":
        return self._map_unary(lambda a: a * float(scalar))  # type: ignore[operator]

    def __rmul__(self, scalar: float) -> "Self":
        return self.__mul__(scalar)

    def __truediv__(self, scalar: float) -> "Self":
        return self._map_unary(lambda a: a / float(scalar))  # type: ignore[operator]

    # ── Flat-vector protocol (Krylov / scipy.gmres adapter) ──────────
    #
    # Direct-sum flat representation: ``concat(interior.values.ravel(),
    # boundary.values)``. The boundary values are already a flat 1-D ndarray (per
    # the trace leaf's flat-backing storage); the interior values are reshaped via
    # :meth:`ndarray.ravel`. The Krylov adapter at
    # :mod:`orpheus.numerics.iteration` consumes this flat representation as the
    # GMRES iterate vector; round-trip exactness is the load-bearing invariant.

    def to_flat(self) -> "NDArray":
        r"""Pack the composite blocks into a flat 1-D vector.

        The packed layout is ``[interior.values.ravel(), boundary.values, ...]``
        — the direct-sum representation, its ordered blocks supplied by
        :meth:`_flat_parts` (which a subclass with an extra block overrides). Lives
        ONCE here; :meth:`from_flat` is its inverse.
        """
        return np.concatenate(self._flat_parts())

    def _flat_parts(self) -> "list[NDArray]":
        r"""The ordered 1-D block arrays :meth:`to_flat` concatenates.

        Base: ``[interior.values.ravel(), boundary.values]``. A subclass carrying
        an extra block (:class:`FullField`'s ψ½) overrides to append it, so the
        flat protocol stays defined ONCE on :class:`Composite`.
        """
        return [
            self.interior.values.ravel(),
            self.boundary.values,  # already 1-D (flat trace storage)
        ]

    @classmethod
    def from_flat(cls, flat: "NDArray", template: "CompositeT") -> "CompositeT":
        r"""Reconstruct a composite from a flat 1-D vector + template.

        The ``template`` provides the shapes, types, AND the concrete composite
        class: reconstruction is delegated to the template's :meth:`_from_flat`
        instance hook, so the result is the SAME concrete type as ``template`` (a
        :class:`~orpheus.transport.timed_full_field.TimedFullField` template yields
        a ``TimedFullField`` with preserved ``history_depth`` and an EMPTY history).

        Parameters
        ----------
        flat : np.ndarray
            1-D vector matching ``template.to_flat()`` in size.
        template : Composite
            Source of structural metadata (shapes, spaces, meshes) AND the
            concrete return type.
        """
        return template._from_flat(flat)

    def _from_flat(self, flat: "NDArray") -> "Self":
        r"""Rebuild a composite of ``self``'s type from a flat vector (self = template).

        The instance-hook inverse of :meth:`_flat_parts`: slices
        ``interior | boundary`` from the template layout, rebuilds each leaf with
        the template's ``space`` / ``mesh``, and routes the pair through
        :meth:`_recombine`. A subclass with an extra block overrides to slice +
        thread it. Instance method (not classmethod) so the ``Self`` override is
        Liskov-clean.
        """
        n_interior = self.interior.values.size
        n_boundary = self.boundary.values.size
        expected_total = n_interior + n_boundary
        if flat.size != expected_total:
            raise ValueError(
                f"{type(self).__name__}.from_flat: flat.size = {flat.size} does "
                f"not match template total size {n_interior} + {n_boundary} = "
                f"{expected_total}"
            )
        interior_values = flat[:n_interior].reshape(self.interior.values.shape)
        boundary_values = flat[n_interior : n_interior + n_boundary]
        return self._recombine(
            interior=replace(self.interior, values=interior_values),
            boundary=replace(self.boundary, values=boundary_values),
        )

    # ── Diagnostics ──────────────────────────────────────────────────

    def copy(self) -> "Self":
        r"""Return a deep copy with owned ndarrays.

        Snapshots every carried block with owned ndarrays (via :meth:`_map_unary`,
        so a subclass's extra block is copied too). Used by callers that need a
        stable iterate without aliasing. Routes through :meth:`_recombine`, so a
        :class:`~orpheus.transport.timed_full_field.TimedFullField` copy is a
        ``TimedFullField`` with EMPTY history (the existing ``copy`` drops history
        — bit-identical behaviour).
        """
        return self._map_unary(lambda a: a.copy())  # type: ignore[attr-defined]


@dataclass(frozen=True, kw_only=True)
class FullField(Composite[BulkField, BoundaryField]):
    r"""The SN composite: ``interior ⊕ boundary`` plus the optional ψ½ block.

    The :class:`Composite` specialization the SN operator algebra acts on —
    ``interior`` an :class:`~orpheus.transport.fields.angular_flux.AngularFlux`,
    ``boundary`` an
    :class:`~orpheus.transport.fields.angular_boundary_flux.AngularBoundaryFlux`
    — extended with the OPTIONAL curvilinear starting-direction (ψ½) third block.
    Bound to the ``BulkField`` / ``BoundaryField`` ABCs (not a single leaf) so the
    same class serves every SN geometry.

    See the module docstring for the ψ½ block's role (#282 route (a)) and the
    cofree-comonad split (:class:`~orpheus.transport.timed_full_field.TimedFullField`).
    The ψ½ block lives HERE, not on :class:`Composite`, because it is a
    curvilinear-SN augmentation — the coupled-block campaign's Phase B evicts it
    into its own System-B composite, after which ``FullField`` is a pure
    ``Composite``.

    Parameters
    ----------
    interior : BulkField
        The volumetric / bulk field (typically ``AngularFlux`` for SN).
    boundary : BoundaryField
        The boundary partner field on the trace of ``interior``'s domain.
    radial_characteristic : RadialCharacteristicField or None
        The OPTIONAL third block (#282 route (a), campaign phase 2.5d): the
        per-μ-level starting-direction flux ψ½ of a curvilinear Morel–Montry
        thread, typed as composite state. ``None`` ⟺ the mesh has NO seed-carrying
        level (the R12a predicate — Cartesian always; every production cylinder
        rule; only the sphere carries it today). Presence is MESH-keyed, never
        caller-keyed: :meth:`zeros` consults the mesh predicate, so a seed block on
        a Cartesian composite is unrepresentable through the factory.
        Mixed-presence arithmetic raises (a partner without the block cannot
        silently drop it).
    """

    #: The optional starting-direction block (#282 route (a)). ``None`` ⟺ the mesh
    #: carries no seed-carrying level (R12a) — absence is a structural fact of the
    #: phase space, not a lazily-unfilled slot.
    radial_characteristic: RadialCharacteristicField | None = None

    # ── Construction (override — the mesh-keyed ψ½ third block) ───────

    @classmethod
    def zeros(  # type: ignore[override]
        cls,
        *,
        interior: type[BulkField],
        boundary: type[BoundaryField],
        mesh: "object",
        radial_characteristic: type[RadialCharacteristicField] | None = None,
    ) -> "Self":
        r"""Allocate a zero SN composite from the leaf TYPES (+ the ψ½ leaf).

        Extends :meth:`Composite.zeros` with the mesh-keyed ψ½ block: the
        starting-direction leaf is allocated iff ``radial_characteristic`` is
        given AND ``mesh.radial_characteristic_space`` is non-``None`` (the R12a
        predicate), so SN call sites pass the class UNIFORMLY across geometries and
        a Cartesian / cylinder composite still (correctly) carries ``None``.
        Methods without the concept omit it.

        Parameters
        ----------
        interior, boundary : type[BulkField] / type[BoundaryField]
            The leaf CLASSES to instantiate (must expose ``zeros_on(mesh)``).
        mesh : object
            The phase-space carrier passed through to each leaf's ``zeros_on``.
        radial_characteristic : type[RadialCharacteristicField], optional
            The starting-direction leaf CLASS (SN passes
            :class:`~orpheus.transport.fields.radial_characteristic_flux.RadialCharacteristicFlux`
            for a flux composite, its SourceSink sibling for a source).
        """
        seed_leaf: RadialCharacteristicField | None = None
        if (
            radial_characteristic is not None
            and getattr(mesh, "radial_characteristic_space", None) is not None
        ):
            seed_leaf = radial_characteristic.zeros_on(mesh)  # type: ignore[arg-type]
        return cls(
            interior=interior.zeros_on(mesh),  # type: ignore[attr-defined]
            boundary=boundary.zeros_on(mesh),  # type: ignore[attr-defined]
            radial_characteristic=seed_leaf,
        )

    # ── Construction validation (override — super + ψ½ checks) ────────

    def __post_init__(self) -> None:
        super().__post_init__()
        if self.radial_characteristic is not None and not isinstance(
            self.radial_characteristic, RadialCharacteristicField
        ):
            raise TypeError(
                f"{type(self).__name__}: radial_characteristic must be a "
                f"RadialCharacteristicField leaf or None; got "
                f"{type(self.radial_characteristic).__name__}"
            )
        boundary_mesh = getattr(self.boundary, "mesh", None)
        if self.radial_characteristic is not None and boundary_mesh is not None:
            if self.radial_characteristic.mesh is not boundary_mesh:
                raise ValueError(
                    f"{type(self).__name__}: radial_characteristic must share "
                    "mesh identity with interior/boundary; got "
                    f"radial_characteristic.mesh="
                    f"{self.radial_characteristic.mesh!r}, "
                    f"boundary.mesh={boundary_mesh!r}"
                )

    # ── Hooks (override — thread the ψ½ third block) ──────────────────

    def _recombine(  # type: ignore[override]
        self,
        *,
        interior: BulkField,
        boundary: BoundaryField,
        radial_characteristic: RadialCharacteristicField | None,
    ) -> "Self":
        r"""Rebuild an SN composite from recombined blocks (interior ⊕ boundary ⊕ ψ½).

        ``radial_characteristic`` is a REQUIRED keyword (no default): every caller
        must state the recombined third block explicitly, so a caller that forgets
        it fails loudly instead of silently dropping the seed — the silent-drop bug
        class the 2.5d carrier gates pin (§16.A A3).
        """
        return replace(
            self,
            interior=interior,
            boundary=boundary,
            radial_characteristic=radial_characteristic,
        )

    def _map_binary(  # type: ignore[override]
        self, other: "Composite", op: "Callable[[object, object], object]",
    ) -> "Self":
        return self._recombine(
            interior=op(self.interior, other.interior),  # type: ignore[arg-type]
            boundary=op(self.boundary, other.boundary),  # type: ignore[arg-type]
            radial_characteristic=self._combine_radial_characteristic(other, op),
        )

    def _map_unary(  # type: ignore[override]
        self, op: "Callable[[object], object]",
    ) -> "Self":
        sd = self.radial_characteristic
        return self._recombine(
            interior=op(self.interior),  # type: ignore[arg-type]
            boundary=op(self.boundary),  # type: ignore[arg-type]
            radial_characteristic=None if sd is None else op(sd),  # type: ignore[arg-type]
        )

    def _combine_radial_characteristic(
        self,
        other: "Composite",
        combine: "Callable[[object, object], object]",
    ) -> RadialCharacteristicField | None:
        r"""Combine the two operands' starting-direction blocks, or ``None``.

        Presence must MATCH: a seeded ⊕ unseeded pair raises — silently dropping
        (or fabricating) the ψ½ block is exactly the bug class the carrier gates
        pin (§16.A A3). When both are absent the result is absent; when both are
        present the combined leaf comes from the member-level algebra (which
        enforces role / mesh / space exactly as for interior and boundary — e.g.
        seed ``flux − flux`` mints a ``RadialCharacteristicDisplacement``).
        """
        mine = self.radial_characteristic
        theirs = getattr(other, "radial_characteristic", None)
        if (mine is None) != (theirs is None):
            raise ValueError(
                f"{type(self).__name__} arithmetic with MIXED "
                f"starting-direction presence: one operand carries the ψ½ block "
                f"and the other does not (self: {mine is not None}, other: "
                f"{theirs is not None}). On a seed-carrying mesh (R12a) every "
                f"composite must carry the block — a partner without it cannot "
                f"silently drop it."
            )
        if mine is None:
            return None
        return combine(mine, theirs)  # type: ignore[return-value]

    # ── Flat-vector protocol (override the hooks — the ψ½ tail) ───────

    def _flat_parts(self) -> "list[NDArray]":
        parts = [
            self.interior.values.ravel(),
            self.boundary.values,  # already 1-D (flat trace storage)
        ]
        if self.radial_characteristic is not None:
            parts.append(self.radial_characteristic.values)  # already 1-D
        return parts

    def _from_flat(self, flat: "NDArray") -> "Self":
        r"""Rebuild an SN composite from a flat vector (ψ½-aware; self = template).

        Slices ``interior | boundary | ψ½`` from the template layout and routes the
        rebuilt blocks through :meth:`_recombine`, so a
        :class:`~orpheus.transport.timed_full_field.TimedFullField` template yields
        a ``TimedFullField`` with an EMPTY history. The ψ½ tail is present exactly
        when the template carries the block.
        """
        n_interior = self.interior.values.size
        n_boundary = self.boundary.values.size
        template_seed = self.radial_characteristic
        n_seed = 0 if template_seed is None else template_seed.values.size
        expected_total = n_interior + n_boundary + n_seed
        if flat.size != expected_total:
            raise ValueError(
                f"{type(self).__name__}.from_flat: flat.size = {flat.size} does "
                f"not match template total size {n_interior} + {n_boundary} + "
                f"{n_seed} = {expected_total}"
            )
        interior_values = flat[:n_interior].reshape(self.interior.values.shape)
        boundary_values = flat[n_interior : n_interior + n_boundary]
        new_seed = (
            None
            if template_seed is None
            else replace(template_seed, values=flat[n_interior + n_boundary :])
        )
        return self._recombine(
            interior=replace(self.interior, values=interior_values),
            boundary=replace(self.boundary, values=boundary_values),
            radial_characteristic=new_seed,
        )
