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
question). It holds explicit cells and, separately, the channels it names in
EVERY material (:meth:`CellCoefficient.every`). The reference specification
RESOLVES it against the materials of its problem
(:meth:`CellCoefficient.resolve`): the resolved key lists only explicit
non-zero cells and names no channel in every material, so a cache key never
holds a quantifier.
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

__all__ = ["CellCoefficient", "Channel"]


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
        predicate the emission-spectrum law keys on: a producing mixture's
        spectrum is a probability simplex, so its emission is non-zero);
        scattering and (n,2n) emission by a mixture with a non-zero entry in
        any Legendre block of the stack.
        """
        match self:
            case Channel.FISSION_EMISSION:
                return mixture.is_producing
            case Channel.SCATTERING_EMISSION:
                return _any_nonzero(mixture.SigS)
            case Channel.N2N_EMISSION:
                return _any_nonzero(mixture.Sig2)


def _any_nonzero(stack: Iterable[Any]) -> bool:
    return any(block.count_nonzero() > 0 for block in stack)


def _admit_cell(cell: Any) -> tuple[int, Channel]:
    if not isinstance(cell, tuple) or len(cell) != 2:
        raise TypeError(f"CellCoefficient: a cell is a (material id, Channel) pair, got {cell!r}")
    material, channel = cell
    return parse_integer(material, f"CellCoefficient: the cell {cell!r}", "a material id"), _admit_channel(channel)


def _admit_channel(channel: Any) -> Channel:
    if not isinstance(channel, Channel):
        raise TypeError(f"CellCoefficient: a channel is a Channel, got a {type(channel).__name__} ({channel!r})")
    return channel


@dataclass(frozen=True, eq=False)
class CellCoefficient(ContentIdentity):
    """The direction that scales a set of ``(material id, Channel)`` cells together.

    ``cells`` is an iterable of explicit pairs and
    ``channels_in_every_material`` an iterable of channels named in every
    material that carries them (spelled :meth:`CellCoefficient.every`); both are frozen at construction, and
    their order and repetition are not content. A direction names at least
    one cell or channel.
    """

    cells: Iterable[tuple[int, Channel]] = ()
    channels_in_every_material: Iterable[Channel] = ()

    def __post_init__(self) -> None:
        cells = frozenset(_admit_cell(cell) for cell in self.cells)
        channels = frozenset(_admit_channel(channel) for channel in self.channels_in_every_material)
        if not cells and not channels:
            raise ValueError("CellCoefficient: a direction names at least one cell")
        object.__setattr__(self, "cells", cells)
        object.__setattr__(self, "channels_in_every_material", channels)

    @classmethod
    def every(cls, *channels: Channel) -> "CellCoefficient":
        """The direction of the given channels in every material that carries them."""
        if not channels:
            raise ValueError("CellCoefficient.every: name at least one channel")
        return cls(channels_in_every_material=channels)

    def resolve(self, materials: "Materials") -> "CellCoefficient":
        """The explicit non-zero cells this direction scales in ``materials``.

        ``materials`` are the materials of the problem (a specification's,
        restricted to those its geometry assigns). A channel named in every
        material becomes the cells of every material that carries it; an
        explicit cell the material does not carry is dropped (it scales a
        zero). A cell on a material outside ``materials`` is refused, and so
        is a direction with no non-zero cell left (a zero direction has no
        pole to find). Resolving a resolved key returns it.
        """
        for material, _ in self.cells:
            if material not in materials:
                raise ValueError(
                    f"CellCoefficient: material {material} is not among the problem's materials "
                    f"(ids: {sorted(materials.ids)})"
                )
        explicit = {(i, channel) for i, channel in self.cells if channel.is_carried_by(materials[i])}
        explicit |= {(i, channel) for channel in self.channels_in_every_material for i, mixture in materials.items() if channel.is_carried_by(mixture)}
        if not explicit:
            named = sorted([f"({m}, {c.value})" for m, c in self.cells] + [f"(every material, {c.value})" for c in self.channels_in_every_material])
            raise ValueError(
                f"CellCoefficient: every cell {named} is zero in the problem's materials, "
                f"so the key is a zero direction"
            )
        return CellCoefficient(explicit)
