r"""The walls of a concentric body: each boundary point with what its face role returns.

A wall is a boundary point of the domain :math:`[r_0, r_n]`, named by its
breakpoint index (:math:`0` or :math:`n`, never by its position in a tuple),
with the three numbers the characteristic closure needs from its face:

* the **specular** amplitude, returned along the reflected line (or, across a
  wrap, along the translated line);
* the **diffuse** amplitude, returned isotropically (Lambertian);
* the **partner**, the breakpoint at which the returned path re-enters: the
  wall itself for a mirror and for diffuse re-emission, the opposite wall for
  a periodic wrap.

The first two are one value, :class:`~orpheus.specification.exterior.WallAmplitudes`
(the cone coordinates :math:`(s, d)` of the wall's linear part), which :class:`Wall`
COMPOSES, so the wall's bounds have one owner and this module states no wall invariant of its own.

**One reader, two plain inputs.** A boundary face of the geometry carries a
role, and :meth:`Walls.of` reads it by one exhaustive ``match``:

* a :class:`~orpheus.geometry.deck.MirrorDeck` is a mirror returning
  everything along the reflected line, :math:`(1, 0)`, to the wall itself;
* a :class:`~orpheus.geometry.deck.TranslationDeck` is a wrap returning
  everything along the translated line, :math:`(1, 0)`, to the opposite wall;
* an integer tag is a wall the specification's ``exterior`` states: the
  ``WallAmplitudes`` of that tag, returned to the wall itself.

A deck states its own response (a symmetry of the domain needs no entry), and
a tag's response is the reference's own closed value; this reference imports
no production response and no typed boundary law (the user's ruling of
2026-10-06, P1 step (b), kept by the wave-2 carve: references closed,
production versatile). The old string-tag registry, which parsed a tag into
a typed law of the retired law family, is gone with the family: the
interface between a reference and production is plain data, crossed by the
one test-side loop ``tests/gates/_reference_pairing.py``.

What the vocabulary cannot state (a prescribed inflow, a combined
specular-and-diffuse wall, a zero flux) has no spelling in
``WallAmplitudes``, and the specification refuses it where it is built; a
curved face's wrap is a geometry refusal
(:class:`~orpheus.geometry.deck.DeckError`) before any wall is read.
"""

from __future__ import annotations

import operator
from collections.abc import Mapping
from dataclasses import dataclass, replace
from typing import NoReturn

import numpy as np

from orpheus.geometry.chart import Chart
from orpheus.geometry.chord import ConcentricPartition
from orpheus.geometry.deck import FaceRole, TranslationDeck
from orpheus.geometry.structured_geometry import StructuredGeometry
from orpheus.derivations.common.reference_body import amplitudes_at
from orpheus.specification.exterior import WallAmplitudes, admit_exterior


def end_on(breakpoint: np.ndarray | int, partition: ConcentricPartition) -> np.ndarray:
    r"""A wall's breakpoint re-keyed onto ``partition``, a refinement of the body's with the same ends.

    A wall is an end of the domain (:class:`Walls`), so :math:`0` stays
    :math:`0` and the body's last breakpoint becomes the refinement's,
    :math:`m`; the one law :meth:`Walls.on` and
    :meth:`~.closure.DiffuseWalls.on` re-key by.
    """
    return np.where(np.asarray(breakpoint) == 0, 0, partition.n_regions)


def _refuse(what: str, missing: str) -> NoReturn:
    """Refuse a boundary declaration the characteristic reference does not serve, naming why."""
    raise NotImplementedError(f"the characteristic reference does not serve a {what}: {missing}")


@dataclass(frozen=True)
class Wall:
    r"""One boundary point and what it returns: its amplitudes, and where the returned path re-enters.

    The amplitudes are one :class:`~orpheus.specification.exterior.WallAmplitudes`, which owns the wall's
    invariant (each amplitude in :math:`[0, 1]`, never more returned than received, specular or diffuse and
    never both), so a wall that violates it cannot be built and this class states none of it. ``specular``
    and ``diffuse`` are read-only views of it.

    Attributes
    ----------
    breakpoint:
        The breakpoint index of the boundary point, :math:`0` or :math:`n`.
    amplitudes:
        The cone coordinates :math:`(s, d)` of the wall: the amplitude returned along the reflected (or translated)
        line and the amplitude returned isotropically.
    partner:
        The breakpoint at which the returned path re-enters.
    """

    breakpoint: int
    amplitudes: WallAmplitudes
    partner: int

    def __post_init__(self) -> None:
        # one value, one memo key: the indices as ``int``; the amplitudes arrive parsed (a ``WallAmplitudes`` is canonical)
        object.__setattr__(self, "breakpoint", operator.index(self.breakpoint))
        object.__setattr__(self, "partner", operator.index(self.partner))
        # an input parse, not a second statement of the bounds: a bare pair is not the value that owns them
        if not isinstance(self.amplitudes, WallAmplitudes):
            raise TypeError(
                f"Wall: the amplitudes at breakpoint {self.breakpoint} are a WallAmplitudes, got a {type(self.amplitudes).__name__}"
            )

    @property
    def specular(self) -> float:
        """The amplitude returned along the reflected (or translated) line."""
        return self.amplitudes.specular

    @property
    def diffuse(self) -> float:
        """The amplitude returned isotropically."""
        return self.amplitudes.diffuse


@dataclass(frozen=True)
class Walls:
    r"""The walls of a body, inner first.

    The ``*_at`` lookups take breakpoint indices (as :class:`~orpheus.geometry.chord.Transits`
    reports them) and raise :class:`IndexError` on a breakpoint that is not a
    wall: an interior breakpoint, a solid body's centre, or the absent-transit
    code :math:`n + 1`.

    Attributes
    ----------
    walls:
        One :class:`Wall` per boundary point, inner first.
    n_regions:
        The number :math:`n` of regions of the body.
    chart:
        The body's chart, which decides whether a translation is a symmetry of its level sets.
    """

    walls: tuple[Wall, ...]
    n_regions: int
    chart: Chart

    def __post_init__(self) -> None:
        ends = {0, self.n_regions}
        breakpoints = [w.breakpoint for w in self.walls]
        if not set(breakpoints) <= ends or len(set(breakpoints)) != len(breakpoints):
            raise ValueError(f"walls sit at distinct breakpoints among {sorted(ends)}; got {breakpoints}")
        by_breakpoint = {w.breakpoint: w for w in self.walls}
        for wall in self.walls:
            if wall.partner == wall.breakpoint:
                continue
            if self.chart.acts_on_kept_space:
                _refuse(
                    "periodic wrap on a cylinder or sphere",
                    "a translation maps a radial chart's level sets onto none of them",
                )
            partner = by_breakpoint.get(wall.partner)
            if partner is None or partner.partner != wall.breakpoint:
                _refuse(
                    f"periodic wrap at breakpoint {wall.breakpoint} whose partner does not wrap back",
                    "a wrap identifies two walls, so both carry it",
                )

    @classmethod
    def of(cls, geometry: StructuredGeometry, exterior: Mapping[int, WallAmplitudes]) -> "Walls":
        r"""The walls of ``geometry``: each boundary face's role read, a deck as itself and a tag through ``exterior``.

        ``exterior`` is admitted as an exact cover of ``geometry.tags`` first
        (:func:`~orpheus.specification.exterior.admit_exterior`), so a missing or an unused tag is refused by number.
        """
        exterior = admit_exterior(exterior, geometry, owner="Walls.of")
        n = len(geometry.breakpoints) - 1
        indices = (n,) if len(geometry.boundary_points) == 1 else (0, n)
        walls = tuple(
            _wall_at_face(role, k, opposite=n - k, exterior=exterior)
            for role, k in zip(geometry.face_roles, indices, strict=True)
        )
        return cls(walls, n, Chart(geometry.coord))

    def on(self, partition: ConcentricPartition) -> "Walls":
        r"""The same walls keyed on ``partition``, a refinement of the body's with the same ends.

        A wall is an end of the domain, so breakpoint :math:`n` becomes the
        refinement's last index and :math:`0` stays :math:`0`; a partner is
        re-keyed with it. That the refinement keeps the body's ends is the
        refinement's own invariant (:class:`~.basis.PanelBasis`).
        """
        if partition.chart != self.chart:
            raise ValueError(f"walls on a {self.chart.coord} body are not re-keyed onto a {partition.chart.coord} partition")
        rekeyed = tuple(
            replace(w, breakpoint=int(end_on(w.breakpoint, partition)), partner=int(end_on(w.partner, partition)))
            for w in self.walls
        )
        return Walls(rekeyed, partition.n_regions, self.chart)

    @property
    def _slot(self) -> np.ndarray:
        slot = np.full(self.n_regions + 1, len(self.walls))      # out of range for the per-wall arrays
        slot[[w.breakpoint for w in self.walls]] = np.arange(len(self.walls))
        return slot

    def _at(self, values: list[float] | list[int], breakpoint: np.ndarray, where: np.ndarray | None) -> np.ndarray:
        breakpoint = np.asarray(breakpoint)
        if where is None:
            return np.asarray(values)[self._slot[breakpoint]]
        looked_up = np.asarray(values)[self._slot[np.where(where, breakpoint, self.walls[0].breakpoint)]]
        return np.where(where, looked_up, 0)

    def specular_at(self, breakpoint: np.ndarray, where: np.ndarray | None = None) -> np.ndarray:
        """The specular amplitude of the wall at each breakpoint; 0 where ``where`` is False."""
        return self._at([w.specular for w in self.walls], breakpoint, where)

    def diffuse_at(self, breakpoint: np.ndarray, where: np.ndarray | None = None) -> np.ndarray:
        """The diffuse amplitude of the wall at each breakpoint; 0 where ``where`` is False."""
        return self._at([w.diffuse for w in self.walls], breakpoint, where)

    def partner_at(self, breakpoint: np.ndarray, where: np.ndarray | None = None) -> np.ndarray:
        """The breakpoint at which the path returned at each breakpoint re-enters; 0 where ``where`` is False."""
        return self._at([w.partner for w in self.walls], breakpoint, where)

# ── the one reader: a face role to a wall ────────────────────────────────


def _wall_at_face(role: FaceRole, breakpoint: int, *, opposite: int, exterior: Mapping[int, WallAmplitudes]) -> Wall:
    r"""The wall at ``breakpoint`` that a face role states: its amplitudes, and where the returned path re-enters.

    The amplitudes are :func:`~orpheus.derivations.common.reference_body.amplitudes_at`'s; the partner is this reader's
    own concern: a translation re-enters at the ``opposite`` breakpoint, every other role at its own.
    """
    partner = opposite if isinstance(role, TranslationDeck) else breakpoint
    return Wall(breakpoint, amplitudes_at(role, exterior), partner)


__all__ = ["Wall", "Walls", "end_on"]
