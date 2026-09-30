r"""The boundary faces of a mesh and the law each one carries.

A mesh's boundary is a finite set of faces, and each face carries one boundary
law. Both the set and its names follow from the mesh's coordinate system and
its extent, by one rule, :func:`face_inventory`:

* along the first axis, the faces are the
  :meth:`~orpheus.geometry.coord.CoordSystem.boundary_points` of its edges:
  both ends on a slab and on a hollow cylinder or sphere, only the outer
  surface on a solid one, whose centre :math:`r = 0` is an interior point and
  carries no law;
* along every further axis (the y of (x, y), the z of (r, z)), both ends;
* each face is named by :attr:`~orpheus.mesh.axis.FaceLabel.face_name`, the
  names S\ :sub:`N`'s boundary table uses: ``xmin``, ``xmax``, ``ymin``,
  ``ymax``, a solid radial axis's outer surface being ``xmax``.

:class:`FaceLaws` is the value both :class:`~orpheus.mesh.structured.Mesh1D`
and :class:`~orpheus.mesh.structured.Mesh2D` store: a mapping from face name
to law over exactly the inventory, in inventory order. A law per face, rather
than per side, is the seed for a side whose faces carry different laws.
"""

from __future__ import annotations

from collections.abc import Iterator, Mapping
from typing import TYPE_CHECKING

import numpy as np

from orpheus.geometry.coord import CoordSystem
from orpheus.geometry.structured_geometry import parse_boundary_law
from orpheus.mesh.axis import FaceLabel

if TYPE_CHECKING:
    from orpheus.geometry.boundary import BC, BoundaryTraceLaw


def face_inventory(
    coord: CoordSystem, first_axis_edges: np.ndarray, dimension: int,
) -> tuple[str, ...]:
    """The names of a mesh's boundary faces, axis by axis, inner face first."""
    along_first = coord.boundary_points(
        float(first_axis_edges[0]), float(first_axis_edges[-1]),
    )
    first = ("min", "max") if len(along_first) == 2 else ("outer",)
    endpoints = (first,) + (("min", "max"),) * (dimension - 1)
    return tuple(
        FaceLabel(axis_index, endpoint).face_name
        for axis_index, axis_endpoints in enumerate(endpoints)
        for endpoint in axis_endpoints
    )


def _why(coord: CoordSystem, inventory: tuple[str, ...], given: tuple[str, ...]) -> str:
    """The reason a declaration's faces differ from the inventory, when there is one."""
    if coord is CoordSystem.CARTESIAN:
        return ""
    if "xmin" in given and "xmin" not in inventory:
        return (
            f" (the centre r = 0 of a solid {coord.name.lower()} body is an "
            f"interior point and carries no law)"
        )
    if "xmin" in inventory and "xmin" not in given:
        return (
            f" (a hollow {coord.name.lower()} body has an inner surface, "
            f"which needs its own law)"
        )
    return ""


class FaceLaws(Mapping[str, "BC | BoundaryTraceLaw"]):
    """One boundary law per face of a mesh: an ordered, frozen, picklable mapping.

    Built by :meth:`over`, which checks the declaration against the mesh's
    :func:`face_inventory`. Equality is a mapping's (the same faces carrying
    equal laws); iteration yields the face names in inventory order.
    """

    __slots__ = ("_items",)

    _items: tuple[tuple[str, "BC | BoundaryTraceLaw"], ...]

    def __init__(self, items: "tuple[tuple[str, BC | BoundaryTraceLaw], ...]") -> None:
        object.__setattr__(self, "_items", tuple(items))

    @classmethod
    def over(
        cls,
        inventory: tuple[str, ...],
        laws: object,
        where: str,
        coord: CoordSystem = CoordSystem.CARTESIAN,
    ) -> "FaceLaws":
        """The laws of ``laws`` (a mapping from face name to law) over ``inventory``, or a keyed refusal."""
        if not isinstance(laws, Mapping):
            raise TypeError(
                f"{where} is a mapping from face name to law over the faces "
                f"{inventory}, got {type(laws).__name__}"
            )
        given = tuple(laws)
        if set(given) != set(inventory):
            raise ValueError(
                f"{where}: this {coord.name.lower()} mesh has the boundary faces "
                f"{inventory}, one law each; got the faces {given}"
                + _why(coord, inventory, given)
            )
        return cls(tuple(
            (face, parse_boundary_law(laws[face], f"{where}[{face!r}]"))
            for face in inventory
        ))

    def __getitem__(self, face: str) -> "BC | BoundaryTraceLaw":
        for name, law in self._items:
            if name == face:
                return law
        raise KeyError(face)

    def __iter__(self) -> Iterator[str]:
        return (name for name, _ in self._items)

    def __len__(self) -> int:
        return len(self._items)

    def __setattr__(self, name: str, value: object) -> None:
        raise AttributeError("FaceLaws is immutable")

    def __reduce__(self):
        return (FaceLaws, (self._items,))

    def __repr__(self) -> str:
        body = ", ".join(f"{name!r}: {law!r}" for name, law in self._items)
        return f"FaceLaws({{{body}}})"


__all__ = ["FaceLaws", "face_inventory"]
