r"""The emission channels of a material and the cell coefficient (#405 P1 step 8).

A material's cross sections form a grid of cells; one cell is the pair
``(material id, channel)``. A :class:`CellCoefficient` is a DIRECTION in the
system's parameter space: the set of cells whose coefficients are scaled
together. The k-eigenvalue is the direction of every fission-emission cell,
the classical c-eigenvalue (secondaries per collision) that of every emission
cell, and a single cell is a one-element set.

**Only the emission channels are cells today** (the user's ruling of
2026-10-02). Every method's fission operator is the fission emission
:math:`\chi \otimes \nu\Sigma_f` alone, and the scattering and (n,2n)
emissions are the other gains, so scaling an emission cell moves exactly the
operator it names. A removal channel (capture, an absorber search) is not a
cell yet: :math:`\Sigma_t` is stored on the ``Mixture`` beside its parts, so
scaling a part would not move the collision operator. The removal cells
arrive with the reaction grid of posing unit 3 (#526), which derives every
total where it is used.

**A key, not a number.** A :class:`CellCoefficient` is the opaque, hashable
parameter key the question values of :mod:`orpheus.numerics.question` hold
(``Eigen(CellCoefficient.every(Channel.FISSION_EMISSION))`` is the k
question). The reference specification RESOLVES it against its materials
(:meth:`CellCoefficient.resolve`): the written key may say "every material"
(:data:`EVERY_MATERIAL`), and the resolved key lists the explicit non-zero
cells, so a cache key never holds a quantifier.
"""

from __future__ import annotations

import enum
from collections.abc import Iterable
from dataclasses import dataclass
from typing import TYPE_CHECKING, Any

from orpheus.numerics.content import ContentIdentity
from orpheus.numerics.scalars import parse_integer

if TYPE_CHECKING:
    from orpheus.data.macro_xs.mixture import Mixture
    from orpheus.data.materials import Materials

__all__ = ["EVERY_MATERIAL", "CellCoefficient", "Channel", "MaterialQuantifier"]


class Channel(enum.Enum):
    """The emission channels whose scaling today's operators realise (a closed set)."""

    FISSION_EMISSION = "fission emission"
    SCATTERING_EMISSION = "scattering emission"
    N2N_EMISSION = "(n,2n) emission"

    def is_carried_by(self, mixture: "Mixture") -> bool:
        r"""Whether the cell ``(mixture, self)`` is non-zero.

        A material CARRIES a cell when the cell's coefficient is non-zero:
        every channel field exists on every ``Mixture``, so existence would
        refuse nothing (the orchestrator's ruling 5 of step 8). Fission
        emission is carried by a producing mixture (:math:`\nu\Sigma_f > 0`,
        :attr:`~orpheus.data.macro_xs.mixture.Mixture.is_producing`, the
        predicate the emission-spectrum law keys on); scattering and (n,2n)
        emission by a mixture with a non-zero entry in any Legendre block of
        the stack.
        """
        if self is Channel.FISSION_EMISSION:
            return mixture.is_producing
        stack = mixture.SigS if self is Channel.SCATTERING_EMISSION else mixture.Sig2
        return any(block.count_nonzero() > 0 for block in stack)


class MaterialQuantifier(enum.Enum):
    """The one quantifier a cell's material may be: every material of the resolving declaration."""

    EVERY = "every material"


EVERY_MATERIAL = MaterialQuantifier.EVERY
"""Stands for every material of the declaration that resolves the key, and carries the channel."""


def _admit_cell(cell: Any) -> tuple[int | MaterialQuantifier, Channel]:
    if not isinstance(cell, tuple) or len(cell) != 2:
        raise TypeError(f"CellCoefficient: a cell is a (material id, Channel) pair, got {cell!r}")
    material, channel = cell
    if not isinstance(channel, Channel):
        raise TypeError(f"CellCoefficient: the channel of {cell!r} is a Channel, got a {type(channel).__name__}")
    if material is EVERY_MATERIAL:
        return material, channel
    return parse_integer(material, f"CellCoefficient: the cell {cell!r}", "a material id"), channel


@dataclass(frozen=True, eq=False)
class CellCoefficient(ContentIdentity):
    """The direction that scales a set of ``(material id, Channel)`` cells together.

    ``cells`` is any iterable of pairs, frozen at construction; a material id
    is an ``int`` or :data:`EVERY_MATERIAL`. The set is the content: the order
    and repetition of the pairs given are not.
    """

    cells: Iterable[tuple[int | MaterialQuantifier, Channel]]

    def __post_init__(self) -> None:
        cells = frozenset(_admit_cell(cell) for cell in self.cells)
        if not cells:
            raise ValueError("CellCoefficient: a direction names at least one cell")
        object.__setattr__(self, "cells", cells)

    @classmethod
    def every(cls, *channels: Channel) -> "CellCoefficient":
        """The direction of the given channels in every material that carries them."""
        if not channels:
            raise ValueError("CellCoefficient.every: name at least one channel")
        return cls((EVERY_MATERIAL, channel) for channel in channels)

    def resolve(self, materials: "Materials") -> "CellCoefficient":
        """The explicit non-zero cells this direction scales in ``materials``.

        :data:`EVERY_MATERIAL` becomes every material that carries the
        channel; an explicit cell the material does not carry is dropped
        (it scales a zero). A material id ``materials`` does not declare is
        refused, and so is a direction with no non-zero cell left (a zero
        direction has no pole to find). Resolving a resolved key returns it.
        """
        explicit: set[tuple[int, Channel]] = set()
        for material, channel in self.cells:
            if material is EVERY_MATERIAL:
                explicit |= {(i, channel) for i, mixture in materials.items() if channel.is_carried_by(mixture)}
                continue
            assert isinstance(material, int)  # type narrowing: the quantifier was handled above
            if material not in materials:
                raise ValueError(
                    f"CellCoefficient: material {material} is not declared "
                    f"(declared ids: {sorted(materials.ids)})"
                )
            if channel.is_carried_by(materials[material]):
                explicit.add((material, channel))
        if not explicit:
            named = sorted((str(m), c.value) for m, c in self.cells)
            raise ValueError(
                f"CellCoefficient: every cell {named} is zero in the declaration, "
                f"so the key is a zero direction"
            )
        return CellCoefficient(explicit)
