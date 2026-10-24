r"""The walls of a reference body, declared as plain values and posed as face roles plus an ``exterior`` (wave 2, S5; plan 6.0 R1 to R3).

Not a test module. Before the carve a reference test stated a body's walls as boundary laws (``BC.vacuum``, ``BC.reflective``,
``BC("partial", {"albedo": a})``, ``AlbedoBoundary(...)``) held by the geometry. After it, a face carries a ROLE and the
reference's walls are plain data of the specification: a role is a deck (a symmetry of the domain: ``MirrorDeck``,
``TranslationDeck``) or an integer tag, and the tags are bound to a :class:`~orpheus.specification.exterior.WallAmplitudes`
by the specification's ``exterior``. This module is the ONE place the tests' declarations become that pair, so that a file
states "a vacuum wall and a mirror" and the rules of the carve (``.claude/plans/wave2_boundary_carve.md`` 6.0) are applied once:

* **R1 tags.** Equal amplitude declarations on one body share one tag, numbered ``0, 1, ...`` in face order (a slab with
  ``vacuum | vacuum`` carries the one tag ``0``; ``vacuum | mirror-wall`` carries two).
* **R2.** A :data:`MIRROR` on a Cartesian face is a ``MirrorDeck()``; on a curved face it is a tag bound to the specular wall
  ``WallAmplitudes(1, 0)`` (the carve's meaning change: a reflecting curved wall is a response, never a deck).
* **R3.** A :data:`WRAP` on a Cartesian face is a ``TranslationDeck()``; on a curved face it stays a ``TranslationDeck`` and the
  geometry refuses it (:class:`~orpheus.geometry.deck.DeckError`), the refusal the curved-wrap gate pins.

The declarations are the values themselves, never strings: :data:`VACUUM`, :data:`WHITE`, :data:`MIRROR`, :data:`WRAP`,
:func:`partial` (a specular return of a fraction) and :func:`white` (a diffuse return of a fraction).
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any, NamedTuple

from orpheus.data.materials import Materials
from orpheus.geometry import CoordSystem, MirrorDeck, StructuredGeometry, TranslationDeck
from orpheus.specification.exterior import WallAmplitudes
from orpheus.specification.specification import GeometrySpecification

__all__ = [
    "MIRROR", "VACUUM", "WHITE", "WRAP", "Posed", "cylinder", "exterior_of", "from_thicknesses", "geometry", "partial", "posed", "roles_of", "slab", "spec", "sphere", "white",
]

Declared = MirrorDeck | TranslationDeck | WallAmplitudes

VACUUM = WallAmplitudes(0.0, 0.0)
"""Nothing returned."""
WHITE = WallAmplitudes(0.0, 1.0)
"""Everything returned isotropically."""
MIRROR = MirrorDeck()
"""Everything returned along the reflected line: a deck on a slab, the specular wall ``(1, 0)`` on a curved face."""
WRAP = TranslationDeck()
"""Everything returned through the opposite face: a deck, Cartesian only."""


def partial(albedo: float) -> WallAmplitudes:
    """The fraction ``albedo`` returned specularly."""
    return WallAmplitudes(float(albedo), 0.0)


def white(albedo: float) -> WallAmplitudes:
    """The fraction ``albedo`` returned isotropically."""
    return WallAmplitudes(0.0, float(albedo))


class Posed(NamedTuple):
    """A geometry with its ``exterior``: the two halves of a body a reference specification takes."""

    geometry: StructuredGeometry
    exterior: dict[int, WallAmplitudes]


def _curved(coord: CoordSystem) -> bool:
    return coord is not CoordSystem.CARTESIAN


def _normalised(coord: CoordSystem, declared: Sequence[Declared]) -> tuple[Declared, ...]:
    """R2: a mirror on a curved face is the specular wall ``(1, 0)``; every other declaration is itself."""
    mirror_wall = WallAmplitudes(1.0, 0.0)
    return tuple(mirror_wall if isinstance(wall, MirrorDeck) and _curved(coord) else wall for wall in declared)


def roles_of(coord: CoordSystem, declared: Sequence[Declared]) -> tuple[Any, ...]:
    """The face roles of ``declared`` (inner first): decks as decks (R2, R3), amplitudes as tags numbered by first appearance (R1)."""
    tags: dict[WallAmplitudes, int] = {}
    return tuple(
        tags.setdefault(wall, len(tags)) if isinstance(wall, WallAmplitudes) else wall
        for wall in _normalised(coord, declared)
    )


def exterior_of(coord: CoordSystem, declared: Sequence[Declared]) -> dict[int, WallAmplitudes]:
    """The ``exterior`` that binds the tags :func:`roles_of` numbers: each tag to the amplitudes declared at its first face."""
    amplitudes = [wall for wall in _normalised(coord, declared) if isinstance(wall, WallAmplitudes)]
    tags: dict[WallAmplitudes, int] = {}
    for wall in amplitudes:
        tags.setdefault(wall, len(tags))
    return {tag: wall for wall, tag in tags.items()}


def _faces(coord: CoordSystem, breakpoints: Sequence[float]) -> int:
    return 2 if not _curved(coord) or breakpoints[0] > 0.0 else 1


def geometry(coord: CoordSystem, breakpoints: Sequence[float], mat_ids: Sequence[int], declared: Sequence[Declared] | None = None) -> StructuredGeometry:
    """The geometry with ``declared`` walls as face roles; ``None`` is vacuum on every boundary face."""
    walls = (VACUUM,) * _faces(coord, breakpoints) if declared is None else declared
    return StructuredGeometry(coord=coord, breakpoints=tuple(breakpoints), mat_ids=tuple(mat_ids), face_roles=roles_of(coord, walls))


def posed(coord: CoordSystem, breakpoints: Sequence[float], mat_ids: Sequence[int], declared: Sequence[Declared] | None = None) -> Posed:
    """The geometry and its ``exterior``, from one list of declared walls."""
    walls = (VACUUM,) * _faces(coord, breakpoints) if declared is None else tuple(declared)
    return Posed(geometry(coord, breakpoints, mat_ids, walls), exterior_of(coord, walls))


def from_thicknesses(coord: CoordSystem, thicknesses: Sequence[float], mat_ids: Sequence[int], declared: Sequence[Declared] | None = None, *, r_0: float = 0.0) -> Posed:
    """The body whose intervals have the given widths (a registry's spelling), posed with its walls."""
    breakpoints = [r_0]
    for t in thicknesses:
        breakpoints.append(breakpoints[-1] + t)
    walls = (VACUUM,) * _faces(coord, breakpoints) if declared is None else tuple(declared)
    built = StructuredGeometry.from_thicknesses(
        coord=coord, thicknesses=thicknesses, mat_ids=mat_ids, face_roles=roles_of(coord, walls), r_0=r_0,
    )
    return Posed(built, exterior_of(coord, walls))


def slab(breakpoints: Sequence[float], mat_ids: Sequence[int], left: Declared, right: Declared) -> Posed:
    """A slab posed with its two walls."""
    return posed(CoordSystem.CARTESIAN, breakpoints, mat_ids, (left, right))


def sphere(breakpoints: Sequence[float], mat_ids: Sequence[int], outer: Declared, inner: Declared | None = None) -> Posed:
    """A sphere posed with its outer wall (and its inner one when hollow)."""
    return posed(CoordSystem.SPHERICAL, breakpoints, mat_ids, (outer,) if inner is None else (inner, outer))


def cylinder(breakpoints: Sequence[float], mat_ids: Sequence[int], outer: Declared, inner: Declared | None = None) -> Posed:
    """A cylinder posed with its outer wall (and its inner one when hollow)."""
    return posed(CoordSystem.CYLINDRICAL, breakpoints, mat_ids, (outer,) if inner is None else (inner, outer))


def spec(body: Posed, mixtures: Mapping[int, Any] | Materials, question: Any) -> GeometrySpecification:
    """The specification of a posed body."""
    materials = mixtures if isinstance(mixtures, Materials) else Materials(dict(mixtures))
    return GeometrySpecification(materials, body.geometry, question, exterior=body.exterior)
