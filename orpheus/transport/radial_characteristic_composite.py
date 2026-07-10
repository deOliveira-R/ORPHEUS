r"""System B as an independent composite: ``Composite[interior ⊕ boundary]`` of ψ½.

The coupled-block campaign's Phase B poses the curvilinear ψ½ ray as **System B** —
its own ``interior ⊕ boundary`` composite, exactly parallel to System A's
:class:`~orpheus.transport.full_field.FullField`
(``Composite[AngularFlux, AngularBoundaryFlux]``). Here that composite is
realized:

.. code-block:: python

    RadialCharacteristicComposite
        = Composite[RadialCharacteristicInteriorFlux,     # the marched cells (A_BB)
                    RadialCharacteristicBoundaryFlux]      # the r = R corner (B_b)

It is a **trivial** subclass — no hook overrides — because System B adds no third
block (unlike ``FullField``'s optional ψ½): the generic
:class:`~orpheus.transport.full_field.Composite` base (its slot bounds relaxed to
``Field`` in Phase B) holds the entire 2-block vector-space algebra, and System B
inherits it whole. The subclass adds only the concrete-locus ``__post_init__``
guard, an ergonomic ``from_mesh`` factory (presence-gated by the leaf factories),
and the split-fidelity **bridge** (:meth:`from_unified` / :meth:`to_unified`) to
and from the still-live unified
:class:`~orpheus.transport.fields.radial_characteristic_flux.RadialCharacteristicFlux`.

The bridge is the additive-migration seam: while the split coexists with the
unified ψ½ leaf (until later Phase-B steps re-type the operators + retire the
unified), ``to_unified(from_unified(u)) == u`` is the "the split IS the unified,
re-typed" contract — the retirement licence.

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

from orpheus.transport.fields._bases import (
    RadialCharacteristicBoundaryField,
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

if TYPE_CHECKING:
    from orpheus.sn.mesh.augmented_mesh import SNMesh

__all__ = ["RadialCharacteristicComposite"]

#: The ψ½ direction legs, in canonical (DAG) order: inward first, then the
#: pole-continued outward leg (the same order the split spaces lay out).
_SIGNS: tuple[int, int] = (-1, +1)


@dataclass(frozen=True, kw_only=True)
class RadialCharacteristicComposite(
    Composite[RadialCharacteristicInteriorFlux, RadialCharacteristicBoundaryFlux],
):
    r"""System B: the ψ½ ray as an ``interior ⊕ boundary`` composite.

    ``interior`` is a
    :class:`~orpheus.transport.fields._bases.RadialCharacteristicInteriorField`
    (the marched cells — flux state or, after ``composite ⊖ composite``, an
    interior displacement); ``boundary`` a
    :class:`~orpheus.transport.fields._bases.RadialCharacteristicBoundaryField`
    (the r = R corner). The 2-block algebra (``±``, scalar ``·``, ``to_flat`` /
    ``from_flat``, ``copy``, the affine torsor propagated to both blocks) is
    inherited whole from :class:`~orpheus.transport.full_field.Composite`.
    """

    def __post_init__(self) -> None:
        # System B narrows the generic slots to the ψ½ split loci (the concrete
        # guard belongs with the concrete specialization, as FullField guards
        # BulkField / BoundaryField). Guard the FIELD BASE, not the flux role, so
        # a displacement composite (from ⊖) is also admitted.
        if not isinstance(self.interior, RadialCharacteristicInteriorField):
            raise TypeError(
                f"{type(self).__name__}: interior must be a "
                f"RadialCharacteristicInteriorField (its flux / displacement "
                f"leaf); got {type(self.interior).__name__}"
            )
        if not isinstance(self.boundary, RadialCharacteristicBoundaryField):
            raise TypeError(
                f"{type(self).__name__}: boundary must be a "
                f"RadialCharacteristicBoundaryField (its flux / displacement "
                f"leaf); got {type(self.boundary).__name__}"
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
        cls, unified: RadialCharacteristicFlux,
    ) -> "RadialCharacteristicComposite":
        r"""Split a unified :class:`RadialCharacteristicFlux` into System B.

        The unified ψ½ leaf's ``cells`` legs become the interior flux, its
        ``corner`` legs the boundary flux — the additive-migration seam. Paired
        with :meth:`to_unified` it round-trips exactly (the "split IS the unified,
        re-typed" retirement licence).
        """
        mesh = unified.mesh
        interior = RadialCharacteristicInteriorFlux.zeros_on(mesh)
        boundary = RadialCharacteristicBoundaryFlux.zeros_on(mesh)
        for level in unified.levels:
            for sign in _SIGNS:
                interior.cells(level, sign)[...] = unified.cells(level, sign)
                boundary.corner(level, sign)[...] = unified.corner(level, sign)
        return cls(interior=interior, boundary=boundary)

    def to_unified(self) -> RadialCharacteristicFlux:
        r"""Recompose System B into a unified :class:`RadialCharacteristicFlux`.

        The inverse of :meth:`from_unified`: the interior ``cells`` + boundary
        ``corner`` legs are gathered back onto the unified leaf's
        ``(level, sign, part)`` layout.
        """
        mesh = self.mesh
        unified = RadialCharacteristicFlux.zeros_on(mesh)
        for level in self.interior.levels:
            for sign in _SIGNS:
                unified.cells(level, sign)[...] = self.interior.cells(level, sign)
                unified.corner(level, sign)[...] = self.boundary.corner(level, sign)
        return unified
