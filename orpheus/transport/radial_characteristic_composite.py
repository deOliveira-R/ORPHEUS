r"""System B as an independent composite: ``Composite[interior ⊕ boundary]`` of ψ½.

The coupled-block campaign's Phase B poses the curvilinear ψ½ ray as **System B** —
its own ``interior ⊕ boundary`` composite, exactly parallel to System A's
:class:`~orpheus.transport.full_field.FullField`
(``Composite[AngularFlux, AngularBoundaryFlux]``). Here that composite is
realized:

.. code-block:: python

    RadialCharacteristicComposite
        = Composite[RadialCharacteristicInteriorField,    # the marched cells (A_BB)
                    RadialCharacteristicBoundaryField]     # the r = R corner (B_b)

It is a **trivial** subclass — no hook overrides — because System B adds no third
block (unlike ``FullField``'s optional ψ½): the generic
:class:`~orpheus.transport.full_field.Composite` base (its slot bounds relaxed to
``Field`` in Phase B) holds the entire 2-block vector-space algebra, and System B
inherits it whole. The subclass adds only the concrete-locus ``__post_init__``
guard, an ergonomic ``from_mesh`` factory (presence-gated by the leaf factories),
and the split-fidelity **bridge** (:meth:`from_unified` / :meth:`to_unified`) to
and from the still-live unified ψ½ leaves.

**Role-erased slots (B.2b DP2, the FullField precedent).** The static parameters
bind the locus FIELD BASES, not the flux leaves — exactly as ``FullField`` binds
``BulkField``/``BoundaryField`` — so ONE composite class carries the flux state
(the iterate), the source emission (an operator ``.apply`` output: ``A_BA`` /
``B_b``), and the displacement (minted per block by ``⊖``). Role identity lives
on the MEMBERS (the Field class-identity gate rejects cross-role sums); a
consumer that needs a specific role parses it off the member (the #289-F2
discipline).

The bridge is the additive-migration seam: while the split coexists with the
unified ψ½ leaf (until Phase C (4e) re-writes the walk slots and retires the
unified), ``to_unified(from_unified(u)) == u`` is the "the split IS the unified,
re-typed" contract — the retirement licence. It is **role-preserving** via the
exact-class leaf table :data:`_UNIFIED_TO_SPLIT` (flux ↔ flux pair, source ↔
source pair, displacement ↔ displacement pair) — one bridge body, no per-role
twins.

References
----------

* ``.claude/plans/coupled_block_operator_campaign.md`` — Phase B (pose System B
  as an independent composite); B.1a (base relaxation), B.1b (split spaces),
  B.1c (split leaves), B.1d (this composite + the bridge).
* ``coding-elegance`` Pattern 2 (the 2-block algebra lives ONCE on ``Composite``;
  System B is a no-op extension) + Pattern 4 (the leaf-type guard makes a wrong
  pairing unconstructable).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

from orpheus.transport.displacements.radial_characteristic_boundary_displacement import (
    RadialCharacteristicBoundaryDisplacement,
)
from orpheus.transport.displacements.radial_characteristic_displacement import (
    RadialCharacteristicDisplacement,
)
from orpheus.transport.displacements.radial_characteristic_interior_displacement import (
    RadialCharacteristicInteriorDisplacement,
)
from orpheus.transport.fields._bases import (
    RadialCharacteristicBoundaryField,
    RadialCharacteristicField,
    RadialCharacteristicInteriorField,
)
from orpheus.transport.fields.radial_characteristic_boundary_flux import (
    RadialCharacteristicBoundaryFlux,
)
from orpheus.transport.fields.radial_characteristic_flux import RadialCharacteristicFlux
from orpheus.transport.fields.radial_characteristic_interior_flux import (
    RadialCharacteristicInteriorFlux,
)
from orpheus.transport.full_field import Composite
from orpheus.transport.source_sinks.radial_characteristic_boundary_source_sink import (
    RadialCharacteristicBoundarySourceSink,
)
from orpheus.transport.source_sinks.radial_characteristic_interior_source_sink import (
    RadialCharacteristicInteriorSourceSink,
)
from orpheus.transport.source_sinks.radial_characteristic_source_sink import (
    RadialCharacteristicSourceSink,
)

if TYPE_CHECKING:
    from orpheus.sn.mesh.augmented_mesh import SNMesh

__all__ = ["RadialCharacteristicComposite"]

#: The ψ½ direction legs, in canonical (DAG) order: inward first, then the
#: pole-continued outward leg (the same order the split spaces lay out).
_SIGNS: tuple[int, int] = (-1, +1)

#: The role-preserving bridge table: unified ψ½ leaf class → its split
#: ``(interior, boundary)`` leaf pair. EXACT-class keyed — a role IS a class
#: identity (the Field arithmetic gate), so the bridge maps flux ↔ flux pair,
#: source ↔ source pair, displacement ↔ displacement pair through ONE body
#: (B.2b DP2). A new unified role joins by adding a row, not a bridge twin.
_UNIFIED_TO_SPLIT: dict[
    type[RadialCharacteristicField],
    tuple[
        type[RadialCharacteristicInteriorField],
        type[RadialCharacteristicBoundaryField],
    ],
] = {
    RadialCharacteristicFlux: (
        RadialCharacteristicInteriorFlux,
        RadialCharacteristicBoundaryFlux,
    ),
    RadialCharacteristicSourceSink: (
        RadialCharacteristicInteriorSourceSink,
        RadialCharacteristicBoundarySourceSink,
    ),
    RadialCharacteristicDisplacement: (
        RadialCharacteristicInteriorDisplacement,
        RadialCharacteristicBoundaryDisplacement,
    ),
}

#: The inverse view, keyed by the split INTERIOR leaf class (the composite's
#: role witness); carries the boundary partner so :meth:`to_unified` can
#: verify the pairing (mixed-role composites are constructable — the slots
#: are role-erased — but never bridgeable).
_SPLIT_TO_UNIFIED: dict[
    type[RadialCharacteristicInteriorField],
    tuple[
        type[RadialCharacteristicField],
        type[RadialCharacteristicBoundaryField],
    ],
] = {
    interior: (unified, boundary)
    for unified, (interior, boundary) in _UNIFIED_TO_SPLIT.items()
}


@dataclass(frozen=True, kw_only=True)
class RadialCharacteristicComposite(
    Composite[RadialCharacteristicInteriorField, RadialCharacteristicBoundaryField],
):
    r"""System B: the ψ½ ray as an ``interior ⊕ boundary`` composite.

    ``interior`` is a
    :class:`~orpheus.transport.fields._bases.RadialCharacteristicInteriorField`
    (the marched cells — flux state, source emission, or, after
    ``composite ⊖ composite``, an interior displacement); ``boundary`` a
    :class:`~orpheus.transport.fields._bases.RadialCharacteristicBoundaryField`
    (the r = R corner). Role-erased slots (the FullField precedent — see the
    module docstring): role identity lives on the members. The 2-block algebra
    (``±``, scalar ``·``, ``to_flat`` / ``from_flat``, ``copy``, the affine
    torsor propagated to both blocks) is inherited whole from
    :class:`~orpheus.transport.full_field.Composite`.
    """

    def __post_init__(self) -> None:
        # System B narrows the generic slots to the ψ½ split loci (the concrete
        # guard belongs with the concrete specialization, as FullField guards
        # BulkField / BoundaryField). Guard the FIELD BASE, not a role leaf, so
        # flux, source (operator emissions), and displacement (from ⊖)
        # composites are all admitted — role-erased slots, B.2b DP2.
        if not isinstance(self.interior, RadialCharacteristicInteriorField):
            raise TypeError(
                f"{type(self).__name__}: interior must be a "
                f"RadialCharacteristicInteriorField (its flux / source / "
                f"displacement leaf); got {type(self.interior).__name__}"
            )
        if not isinstance(self.boundary, RadialCharacteristicBoundaryField):
            raise TypeError(
                f"{type(self).__name__}: boundary must be a "
                f"RadialCharacteristicBoundaryField (its flux / source / "
                f"displacement leaf); got {type(self.boundary).__name__}"
            )
        super().__post_init__()  # base: Field + mesh-identity.

    # ── Construction ─────────────────────────────────────────────────

    @classmethod
    def from_mesh(cls, mesh: "SNMesh") -> "RadialCharacteristicComposite":
        r"""A zero ψ½ flux composite on ``mesh`` (presence-gated).

        Builds the two zero flux leaves from ``mesh``'s split spaces — so on a
        NON-carrying mesh (no R12a seed level) the leaf factories raise, i.e.
        System B is unconstructable exactly where it does not exist (presence =
        block existence).
        """
        return cls(
            interior=RadialCharacteristicInteriorFlux.zeros_on(mesh),
            boundary=RadialCharacteristicBoundaryFlux.zeros_on(mesh),
        )

    # ── The split-fidelity bridge (to / from the unified ψ½ leaf) ─────

    @classmethod
    def from_unified(
        cls, unified: RadialCharacteristicField,
    ) -> "RadialCharacteristicComposite":
        r"""Split a unified ψ½ leaf into System B — role-preserving.

        The unified leaf's ``cells`` legs become the interior member, its
        ``corner`` legs the boundary member — the additive-migration seam. The
        member CLASSES follow the unified leaf's role through the exact-class
        table :data:`_UNIFIED_TO_SPLIT` (flux → flux pair, source → source
        pair, displacement → displacement pair). Paired with :meth:`to_unified`
        it round-trips exactly per role (the "split IS the unified, re-typed"
        retirement licence).
        """
        split = _UNIFIED_TO_SPLIT.get(type(unified))
        if split is None:
            raise TypeError(
                f"RadialCharacteristicComposite.from_unified: no split leaf "
                f"pair for the unified role {type(unified).__name__} — the "
                f"bridge is role-preserving and knows "
                f"{sorted(c.__name__ for c in _UNIFIED_TO_SPLIT)}. A new "
                f"unified role joins by adding a _UNIFIED_TO_SPLIT row."
            )
        interior_cls, boundary_cls = split
        mesh = unified.mesh
        interior = interior_cls.zeros_on(mesh)
        boundary = boundary_cls.zeros_on(mesh)
        for level in unified.levels:
            for sign in _SIGNS:
                interior.cells(level, sign)[...] = unified.cells(level, sign)
                boundary.corner(level, sign)[...] = unified.corner(level, sign)
        return cls(interior=interior, boundary=boundary)

    def to_unified(self) -> RadialCharacteristicField:
        r"""Recompose System B into its unified ψ½ leaf — role-preserving.

        The inverse of :meth:`from_unified`: the interior ``cells`` + boundary
        ``corner`` legs are gathered back onto the unified leaf's
        ``(level, sign, part)`` layout. The unified CLASS follows the
        composite's role, witnessed by the interior member through the inverse
        table (and the boundary member must carry the SAME role — mixed-role
        composites are constructable, the slots being role-erased, but have no
        unified counterpart and refuse loudly here).
        """
        entry = _SPLIT_TO_UNIFIED.get(type(self.interior))
        if entry is None:
            raise TypeError(
                f"RadialCharacteristicComposite.to_unified: no unified role "
                f"for the interior member {type(self.interior).__name__} — "
                f"the bridge is role-preserving and knows "
                f"{sorted(c.__name__ for c in _SPLIT_TO_UNIFIED)}."
            )
        unified_cls, boundary_cls = entry
        if type(self.boundary) is not boundary_cls:
            raise TypeError(
                f"RadialCharacteristicComposite.to_unified: mixed-role "
                f"composite — interior {type(self.interior).__name__} pairs "
                f"with boundary {boundary_cls.__name__}, got "
                f"{type(self.boundary).__name__}. A mixed-role composite has "
                f"no unified counterpart."
            )
        mesh = self.mesh
        unified = unified_cls.zeros_on(mesh)
        for level in self.interior.levels:
            for sign in _SIGNS:
                unified.cells(level, sign)[...] = self.interior.cells(level, sign)
                unified.corner(level, sign)[...] = self.boundary.corner(level, sign)
        return unified
