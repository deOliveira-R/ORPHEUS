r"""Declare the exterior conditions of a test mesh from the tags its faces carry (wave 2, S2b codemod rules R1 to R8).

Not a test module (no ``test_`` prefix). A face's role is a deck (a symmetry of the domain: ``MirrorDeck``, ``TranslationDeck``) or
an ``int`` tag of the exterior, and the condition a tag stands for is declared where the method is known: the keyword-only
``boundary_conditions`` of ``SNProblem`` and of every ``solve_sn*`` entry (``.claude/plans/wave2_boundary_carve.md`` section 2.4,
2.10 and 6.0, rule R1: faces whose old declaration was one value at one site share one tag, and the solve at that site receives
``boundary_conditions={tag: condition}``).

A test file declares ONE table, ``{tag: condition}``, for the tags its meshes carry (the codemod numbers a tag per distinct
declaration of the file: the same condition is the same tag, and every ``PrescribedInflow`` of the old spelling is its own tag, so
the codemod never manufactures a symmetry claim the old code did not make, law L-C). :func:`conditions_for` then selects, for one
carrier, the entries of the table its faces carry. It exists because the exact cover is strict (a condition for a tag no face carries
is refused, ``UnusedConditionError``) while a test file builds meshes of several kinds (deck-only meshes need no condition, tagged
ones need theirs), and because a mesh travels through helper functions between the line that builds it and the line that solves it.
The selection is mechanical and visible: it never fills in a missing tag (a tag the table lacks is a ``KeyError`` naming the tag), and
it never invents a condition for a tag.
"""

from __future__ import annotations

from collections.abc import Mapping
from typing import Any

from orpheus.geometry.deck import DeckTransformation
from orpheus.mesh.axis import roles_of_axes

__all__ = ["conditions_for", "tags_of_carrier"]


def tags_of_carrier(carrier: Any) -> tuple[int, ...]:
    """The integer tags the faces of ``carrier`` carry, ascending: a mesh or problem (``face_roles``) or a tuple of axes."""
    roles = carrier.face_roles if hasattr(carrier, "face_roles") else roles_of_axes(tuple(carrier))
    return tuple(sorted({int(role) for role in roles.values() if not isinstance(role, DeckTransformation)}))


def conditions_for(carrier: Any, table: Mapping[int, Any]) -> dict[int, Any]:
    """The entries of ``table`` for the tags ``carrier``'s faces carry (empty for a deck-only carrier).

    Raises
    ------
    KeyError
        A tag the carrier carries has no entry in ``table`` (the table is the file's declaration of what its tags mean; a missing
        tag is a file that built a face this table does not describe, never a tag to be filled in).
    """
    tags = tags_of_carrier(carrier)
    missing = [tag for tag in tags if tag not in table]
    if missing:
        raise KeyError(f"the table declares {sorted(table)} but the carrier's faces carry the tags {list(tags)}: no condition for {missing}")
    return {tag: table[tag] for tag in tags}
