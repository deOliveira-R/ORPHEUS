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


from orpheus.geometry.coord import CoordSystem
from orpheus.geometry.deck import FaceRole, parse_face_roles, tags_of
from orpheus.geometry.faces import face_inventory, why_faces_differ
from orpheus.values.content import FrozenMapping


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
                + why_faces_differ(coord, inventory, given)
            )
        parsed = parse_face_roles({face: roles[face] for face in inventory}, coord, where)
        return cls(tuple(parsed.items()))

    @property
    def tags(self) -> frozenset[int]:
        """The tags its faces carry (the decks excluded)."""
        return tags_of(self.values())


__all__ = ["FaceRoles", "face_inventory"]  # face_inventory is the geometry's, re-exported
