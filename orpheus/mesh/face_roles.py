r"""The boundary faces of a mesh and the role each one carries.

A mesh's boundary is a finite set of faces, and each face carries one role: a
deck (the domain continues across it) or an ``int`` tag (the domain ends there;
the problem binds the tag to an exterior condition). Both the set and its names follow from the mesh's coordinate system and
its extent, by one rule, :func:`face_inventory`:

* along the first axis, the faces are the
  :meth:`~orpheus.geometry.coord.CoordSystem.boundary_points` of its edges:
  both ends on a slab and on a hollow cylinder or sphere, only the outer
  surface on a solid one, whose centre :math:`r = 0` is an interior point and
  carries no role;
* along every further axis (the y of (x, y), the z of (r, z)), both ends;
* each face is named by :attr:`~orpheus.mesh.axis.FaceLabel.face_name`, the
  names S\ :sub:`N`'s boundary table uses: ``xmin``, ``xmax``, ``ymin``,
  ``ymax``, a solid radial axis's outer surface being ``xmax``.

:class:`FaceRoles` is the value both :class:`~orpheus.mesh.structured.Mesh1D`
and :class:`~orpheus.mesh.structured.Mesh2D` store: a mapping from face name
to role over exactly the inventory, in inventory order. A role per face, rather
than per side, is the seed for a side whose faces carry different roles.
"""

from __future__ import annotations

from collections.abc import Mapping

import numpy as np

from orpheus.geometry.coord import CoordSystem
from orpheus.geometry.deck import DeckTransformation, FaceRole, parse_face_roles
from orpheus.mesh.axis import FaceLabel
from orpheus.values.content import FrozenMapping


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
            f"interior point and carries no role)"
        )
    if "xmin" in inventory and "xmin" not in given:
        return (
            f" (a hollow {coord.name.lower()} body has an inner surface, "
            f"which needs its own role)"
        )
    return ""


class FaceRoles(FrozenMapping[str, FaceRole]):
    """One boundary role per face of a mesh: an ordered, frozen, picklable mapping.

    Built by :meth:`over`, which checks the declaration against the mesh's
    :func:`face_inventory`. Iteration yields the face names in inventory
    order. A :class:`~orpheus.values.content.FrozenMapping`, so equality
    and hash are content identity: the same faces carrying equal roles,
    whatever the order the caller gave; it is equal only to another
    ``FaceRoles``.
    """

    __slots__ = ()

    @classmethod
    def over(
        cls,
        inventory: tuple[str, ...],
        roles: object,
        where: str,
        *,
        coord: CoordSystem,
    ) -> "FaceRoles":
        """The roles of ``roles`` (a mapping from face name to role) over ``inventory``, or a keyed refusal.

        ``coord`` is required: a default of ``CARTESIAN`` would admit a mirror on a radial face for a caller that
        omits it. Each role is a deck or an ``int`` tag; a deck must be one the face's axis admits on ``coord``
        (:func:`~orpheus.geometry.deck.admitted_decks`) and the pairings close
        (:func:`~orpheus.geometry.deck.parse_face_roles`).
        """
        if not isinstance(roles, Mapping):
            raise TypeError(
                f"{where} is a mapping from face name to role over the faces "
                f"{inventory}, got {type(roles).__name__}"
            )
        given = tuple(roles)
        if set(given) != set(inventory):
            raise ValueError(
                f"{where}: this {coord.name.lower()} mesh has the boundary faces "
                f"{inventory}, one role each; got the faces {given}"
                + _why(coord, inventory, given)
            )
        parsed = parse_face_roles({face: roles[face] for face in inventory}, coord, where)
        return cls(tuple(parsed.items()))

    @property
    def tags(self) -> frozenset[int]:
        """The tags its faces carry (the decks excluded)."""
        return frozenset(role for role in self.values() if not isinstance(role, DeckTransformation))


__all__ = ["FaceRoles", "face_inventory"]
