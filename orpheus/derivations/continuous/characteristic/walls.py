r"""The walls of a concentric body: each boundary point with what its law returns.

A wall is a boundary point of the domain :math:`[r_0, r_n]`, named by its
breakpoint index (:math:`0` or :math:`n`, never by its position in a tuple),
with the three numbers the characteristic closure needs from its law:

* the **specular** amplitude, returned along the reflected line (or, across a
  wrap, along the translated line);
* the **diffuse** amplitude, returned isotropically (Lambertian);
* the **partner**, the breakpoint at which the returned path re-enters: the
  wall itself for a mirror and for diffuse re-emission, the opposite wall for
  a periodic wrap.

**One reader.** Every wall is read from its law's two factors, the deck
(``geometry_map``) and the response (``response_kernel``), by one table
(``_wall_of``). A ``BC`` tag is first parsed into the typed law it names
by the reference's own registry (:data:`TAG_REGISTRY`), so a tag and its law
cannot be read differently.

**Why the registry is the reference's own.** A tag is a name each method
resolves through its own registry. Production's parse
(``orpheus.transport.method._law_from_tag``) is versatile and changes with
production, and a closed reference must not move when production does (the
user's ruling of 2026-10-06, P1 step (b)). The two parses are a declared
duplicate across the branch line, recorded in the conceptual view.
"""

from __future__ import annotations

from collections.abc import Callable
from typing import NoReturn
from dataclasses import dataclass, replace

import numpy as np

from orpheus.geometry.boundary import (
    BC,
    AlbedoBoundary,
    BoundaryTraceLaw,
    LambertianReemission,
    NoSource,
    PairedDeck,
    PeriodicBoundary,
    ReflectiveBoundary,
    ScalarResponse,
    SelfPairedDeck,
    SpecularReemission,
    SpecularReturn,
    VacuumInflow,
    WhiteBoundary,
)
from orpheus.geometry.chart import Chart
from orpheus.geometry.chord import ConcentricPartition
from orpheus.geometry.structured_geometry import StructuredGeometry



def _refuse(what: str, missing: str) -> NoReturn:
    """Refuse a boundary declaration the characteristic reference does not serve, naming why."""
    raise NotImplementedError(f"the characteristic reference does not serve a {what}: {missing}")


@dataclass(frozen=True)
class Wall:
    r"""One boundary point and what its law returns.

    Attributes
    ----------
    breakpoint:
        The breakpoint index of the boundary point, :math:`0` or :math:`n`.
    specular:
        The amplitude returned along the reflected (or translated) line, in :math:`[0, 1]`.
    diffuse:
        The amplitude returned isotropically, in :math:`[0, 1]`.
    partner:
        The breakpoint at which the returned path re-enters.
    """

    breakpoint: int
    specular: float
    diffuse: float
    partner: int

    def __post_init__(self) -> None:
        if not (0.0 <= self.specular <= 1.0 and 0.0 <= self.diffuse <= 1.0):
            _refuse(
                f"wall at breakpoint {self.breakpoint} returning (specular {self.specular}, "
                f"diffuse {self.diffuse})",
                "a returned amplitude outside [0, 1] is not a physical wall",
            )


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
    def of(cls, geometry: StructuredGeometry) -> "Walls":
        r"""The walls of ``geometry``: each boundary law read through its factors."""
        n = len(geometry.breakpoints) - 1
        indices = (n,) if len(geometry.boundary_points) == 1 else (0, n)
        walls = tuple(
            _wall_of(_law_of(declared, outward_sign=-1 if k == 0 else +1), k, opposite=n - k)
            for declared, k in zip(geometry.boundaries, indices, strict=True)
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
        n, m = self.n_regions, partition.n_regions
        index = {0: 0, n: m}
        rekeyed = tuple(replace(w, breakpoint=index[w.breakpoint], partner=index[w.partner]) for w in self.walls)
        return Walls(rekeyed, m, self.chart)

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

# ── the reference's tag registry ─────────────────────────────────────────


def _parameters(tag: BC, *names: str) -> tuple[float, ...]:
    r"""The tag's parameters ``names``, exactly: a missing or an undeclared one is refused, never dropped."""
    if set(tag.params) != set(names):
        _refuse(
            f"boundary tag {tag!r}",
            f"the tag kind {tag.kind!r} takes exactly the parameters {list(names)}",
        )
    return tuple(float(tag.params[name]) for name in names)


#: Each tag kind the reference admits, with the wall's outward sign, to the law it names.
#: Partial white has no tag here: it is spelled as the typed ``WhiteBoundary(albedo=...)``.
TAG_REGISTRY: dict[str, Callable[[BC, int], BoundaryTraceLaw]] = {
    "vacuum": lambda tag, sign: (_parameters(tag), VacuumInflow())[1],
    "reflective": lambda tag, sign: (_parameters(tag), ReflectiveBoundary(axis="x"))[1],
    "partial": lambda tag, sign: AlbedoBoundary(*_parameters(tag, "albedo"), SpecularReturn(axis="x")),
    "white": lambda tag, sign: (_parameters(tag), WhiteBoundary(axis="x", outward_sign=sign))[1],
    "periodic": lambda tag, sign: (_parameters(tag), PeriodicBoundary(axis="x"))[1],
}


def _law_of(declared: BC | BoundaryTraceLaw, *, outward_sign: int) -> BoundaryTraceLaw:
    r"""The typed law a boundary declaration names: a law as it is, a tag through the registry."""
    if isinstance(declared, BoundaryTraceLaw):
        return declared
    parse = TAG_REGISTRY.get(declared.kind)
    if parse is None:
        _refuse(
            f"boundary tag {declared!r}",
            f"the reference admits the tag kinds {sorted(TAG_REGISTRY)}",
        )
    return parse(declared, outward_sign)


# ── the one reader: a law's two factors to a wall ────────────────────────


def _wall_of(law: BoundaryTraceLaw, breakpoint: int, *, opposite: int) -> Wall:
    r"""The wall at ``breakpoint`` that ``law`` declares, read from its deck and its response.

    **SCOPE-BOUNDARY[guard]** machinery: a boundary-source question (an inflow the closure adds), and an unstated re-emission shape.
    ruling: the user, 2026-10-06, P1 step (b) first rung (`.claude/plans/characteristic_reference_architecture.md`).
    revisit: when the reference vocabulary gains a boundary-source question, the source arm is served.
    """
    if not isinstance(law.source, NoSource):
        _refuse(
            f"boundary law {law!r} with an inflow source",
            "the reference's questions carry no boundary source",
        )
    deck, response = law.geometry_map, law.response_kernel
    match deck, response:
        case PairedDeck(), ScalarResponse() if response.amplitude == 1.0:
            wall = Wall(breakpoint, 1.0, 0.0, opposite)
        case SelfPairedDeck(), ScalarResponse() if not deck.is_identity and response.amplitude == 1.0:
            wall = Wall(breakpoint, 1.0, 0.0, breakpoint)
        case SelfPairedDeck(), ScalarResponse() if deck.is_identity and response.is_zero:
            wall = Wall(breakpoint, 0.0, 0.0, breakpoint)
        case SelfPairedDeck(), SpecularReemission() if deck.is_identity:
            wall = Wall(breakpoint, response.amplitude, 0.0, breakpoint)
        case SelfPairedDeck(), LambertianReemission() if deck.is_identity:
            wall = Wall(breakpoint, 0.0, response.amplitude, breakpoint)
        case _:
            _refuse(
                f"boundary law {law!r}",
                f"no wall reads from the deck {type(deck).__name__} with the response "
                f"{response!r}: a wall is vacuum (the identity deck, nothing returned), a mirror "
                "or a wrap returning everything, or a specular or Lambertian re-emission on the "
                "identity deck; anything else is an unstated re-emission shape or no transport wall",
            )
    return wall


__all__ = ["TAG_REGISTRY", "Wall", "Walls"]
