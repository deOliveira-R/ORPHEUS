r"""Binding the declared conditions to a carrier's faces (``binding.py`` of the conditions package).

**The mathematics.** Each boundary face has exactly one role: glued (a deck: the domain continues) or facing the
exterior (a tag: the domain ends, ``FaceRole = DeckTransformation | int``). A *binding* is a map
:math:`\beta : T \to \text{BoundaryCondition}` from the set :math:`T` of tags the faces carry, supplied where the method
is known. It is lawful when (a) it is an EXACT COVER of the tags (every tag a face carries has a condition, every
condition's tag is carried by a face) and (b) it is INVARIANT: for every non-identity :math:`h` in the problem's
colouring group :math:`H` that moves a face of tag :math:`t`, :math:`h_*\beta(t) = \beta(t)`. Every shared kernel is
face-local, so (b) holds for it by a theorem (in one dimension; see ``condition.py``); it can fail only for a condition
whose content is stated in a frame the face does not supply (a prescribed inflow written in laboratory directions).

**The types.**

* :class:`BoundaryCarrier` (a ``Protocol``): the object that owns the faces AND their colouring symmetry (decision 51;
  the carrier ruling of 2026-10-23). ``face_roles`` is a ``Mapping[str, FaceRole]`` (the contract is what ``bind``
  reads: it iterates face to role; the tags predicate is the free function ``orpheus.geometry.deck.tags_of``) and
  ``symmetry`` is the colouring group of THE SAME faces. ``Mesh1D`` and ``Mesh2D`` conform with no edit (they store a
  ``FrozenMapping`` of roles), ``MaterialMesh`` conforms through ``roles_of_axes`` (slice S2b). Passing the carrier,
  not its two derived values, makes a binding checked against another carrier's group unspellable (a G-set is the
  group together with the faces it acts on).
* :class:`Binding`: the result. ``decks`` and ``tags`` are two typed maps with disjoint keys; ``declared`` the
  conditions; ``symmetry`` the group. The face-to-condition join is DERIVED from ``(tags, declared)`` and never stored,
  so a binding that joins a face to a condition of another tag cannot be built. The parse re-checks the exact cover
  (``set(tags.values()) == declared.tags``): the refusal text of :func:`bind` (face-aware) and the parse of the value
  (structural) are the same fact stated at the door and at the type.
* :func:`bind`: the ONE shared body every method calls (it replaces ``resolve_boundary_conditions``, ``_law_from_tag``
  and the six string tables). It coerces ``conditions`` once through ``BoundaryConditions.of`` and refuses, in this
  order: a tag a face carries with no condition (:class:`~.errors.UnboundTagError`, naming the face and the tag); a
  condition whose tag no face carries (:class:`~.errors.UnusedConditionError`: the exact cover catches a typo that a
  reachable-subset rule would hide); a condition that is not invariant (:class:`~.errors.NonInvariantConditionError`).

**Uncomputed groups.** A carrier whose group is :class:`~orpheus.geometry.colouring.UncomputedSymmetries` (a
``Mesh2D``) has no group machinery yet: ``bind`` refuses, for any tag with two or more faces, a condition that is not
face-local, with a message naming the 2-D colouring group as the missing machinery (a declared scope boundary, tagged at
the guard in :func:`_refuse_non_invariant`). Reducing the group instead of refusing (ruling 2's
alternative) is not built: nothing consumes the problem's group today, so a reduced group has no reader.

**Gotchas.**

* A deck-only carrier needs no condition: ``BoundaryConditions.empty()`` binds it.
* ``Binding.sources`` reads ``BoundaryCondition.source`` of EVERY condition type (no ``isinstance`` on a leaf), so a
  method-only condition with a source is not missed.
* ``Binding`` holds contentless inflow samplers by identity: its digest raises ``ContentlessError`` when a condition
  holds a sampled term (the uniform-only case is ordinary content).

Layer: ``transport`` (L2). Imports ``geometry`` (the deck and colouring types), ``values`` and its siblings.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from typing import Protocol, runtime_checkable

from orpheus.geometry.colouring import Symmetries, UncomputedSymmetries
from orpheus.geometry.deck import DeckTransformation, FaceRole, tags_of
from orpheus.geometry.faces import FACE_NAMES
from orpheus.transport.boundary_conditions.condition import BoundaryCondition
from orpheus.transport.boundary_conditions.declared import BoundaryConditions
from orpheus.transport.boundary_conditions.errors import NonInvariantConditionError, UnboundTagError, UnusedConditionError
from orpheus.transport.boundary_conditions.inflow import InflowTerms
from orpheus.values.content import ContentIdentity, FrozenMapping

__all__ = ["Binding", "BoundaryCarrier", "bind"]


@runtime_checkable
class BoundaryCarrier(Protocol):
    """The object that owns the faces and their colouring symmetry."""

    @property
    def face_roles(self) -> Mapping[str, FaceRole]:
        """The face inventory with each face's role (a deck or a tag)."""
        ...

    @property
    def symmetry(self) -> Symmetries:
        """The colouring group of THE SAME faces."""
        ...


@dataclass(frozen=True, eq=False)
class Binding(ContentIdentity):
    """The declared conditions joined to a carrier's faces: ``decks``, ``tags``, ``declared`` and the ``symmetry``."""

    decks: FrozenMapping[str, DeckTransformation]
    tags: FrozenMapping[str, int]
    declared: BoundaryConditions
    symmetry: Symmetries

    def __post_init__(self) -> None:
        if set(self.decks) & set(self.tags):
            raise ValueError(f"Binding: the faces {sorted(set(self.decks) & set(self.tags))} are both a deck and a tag")
        if set(self.tags.values()) != set(self.declared.tags):
            raise ValueError(f"Binding: the tags {sorted(set(self.tags.values()))} the faces carry are not the declared tags {sorted(self.declared.tags)}")

    @property
    def faces(self) -> tuple[str, ...]:
        """Every face, decks and tagged faces merged, in the inventory order of ``FACE_NAMES``."""
        held = {*self.decks, *self.tags}
        return tuple(face for face in FACE_NAMES if face in held)

    def condition_of(self, face: str) -> BoundaryCondition:
        """The condition of a tagged face: ``declared[tags[face]]``."""
        return self.declared[self.tags[face]]

    def faces_of(self, tag: int) -> tuple[str, ...]:
        """The faces carrying ``tag``, in inventory order."""
        return tuple(face for face in self.faces if self.tags.get(face) == tag)

    @property
    def sources(self) -> FrozenMapping[str, InflowTerms]:
        """The tagged faces whose condition has a nonzero source, read through ``BoundaryCondition.source`` of every condition."""
        return FrozenMapping({face: self.condition_of(face).source for face in self.faces if face in self.tags and not self.condition_of(face).source.is_zero})


def bind(carrier: BoundaryCarrier, conditions: BoundaryConditions | Mapping[int, BoundaryCondition]) -> Binding:
    """Bind ``conditions`` to the faces of ``carrier``: the exact cover, then the invariance (see the module docstring)."""
    declared = BoundaryConditions.of(conditions)
    roles, symmetry = carrier.face_roles, carrier.symmetry
    decks = {face: role for face, role in roles.items() if isinstance(role, DeckTransformation)}
    tags = {face: role for face, role in roles.items() if face not in decks}
    for face, tag in tags.items():
        if tag not in declared:
            raise UnboundTagError(f"face {face!r} carries tag {tag} and no condition is declared for it", tag=tag, face=face)
    carried = tags_of(roles.values())
    for tag in declared:
        if tag not in carried:
            raise UnusedConditionError(f"a condition is declared for tag {tag} and no face carries it (carried tags: {sorted(carried)})", declared[tag], tag=tag)
    _refuse_non_invariant(tags, declared, symmetry)
    return Binding(FrozenMapping(decks), FrozenMapping(tags), declared, symmetry)


def _refuse_non_invariant(tags: Mapping[str, int], declared: BoundaryConditions, symmetry: Symmetries) -> None:
    """Raise :class:`NonInvariantConditionError` for a condition that a non-identity element of the group does not fix."""
    if isinstance(symmetry, UncomputedSymmetries):
        # SCOPE-BOUNDARY[guard] (#529) machinery: the 2-D colouring group of a Mesh2D (today UncomputedSymmetries), consumed by a
        # quotient or a cache. ruling: wave2_boundary_carve.md section 2.4, the orchestrator, 2026-10-22 (reducing the group instead of
        # refusing has no reader). revisit: the first consumer of the problem's symmetry group.
        for tag, condition in declared.items():
            faces = tuple(face for face, t in tags.items() if t == tag)
            if len(faces) >= 2 and not condition.face_local:
                raise NonInvariantConditionError(
                    f"tag {tag} on {faces} holds a condition that is not face-local, and the colouring group of {symmetry.carrier} is not computed "
                    "(missing machinery: the 2-D colouring group, #529)",
                    condition, tag=tag, isometry=None, faces=faces,
                )
        return
    for isometry in symmetry.elements:
        if isometry.is_identity:
            continue
        for tag, condition in declared.items():
            faces = tuple(face for face, t in tags.items() if t == tag)
            if any(isometry.faces[face] != face for face in faces) and not condition.invariant_under(isometry):
                raise NonInvariantConditionError(f"the condition of tag {tag} on {faces} is not invariant under {isometry}", condition, tag=tag, isometry=isometry, faces=faces)
