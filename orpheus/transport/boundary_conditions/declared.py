r"""The conditions a problem declares (``declared.py`` of the conditions package).

:class:`BoundaryConditions` is the exact twin of ``orpheus.data.materials.Materials``: where ``Materials`` maps a
material id (an ``int``) to a material, ``BoundaryConditions`` maps a TAG (an ``int``, the colour a face carries) to
a :class:`~.condition.BoundaryCondition`. The geometry says which faces share an exterior model and never what the
model is (tags are a declaration, not an observation); this value says what the model is, and ``bind`` joins the two.

**Parse at the boundary.** Keys are ``int`` (a ``bool`` is refused with ``TypeError``, as in ``Materials``); every
value must be a ``BoundaryCondition`` (a ``WallDisplacement``, a deck and a ``str`` are refused with a ``TypeError``
naming the type: this is the one narrowing site for displacements). The mapping is stored as a
:class:`~orpheus.values.content.FrozenMapping` in declared order.

**One difference from ``Materials``: empty is admitted.** A carrier of decks only (such as a half-cell slab whose two
faces are mirrors) has no tag, so ``BoundaryConditions.empty()`` is a valid declaration, and the default of the method
heads. The exact cover performed by ``bind`` is the witness of the default: a tagged carrier with no binding is
refused loudly, naming each tag.

**No ``restrict``.** ``Materials.restrict`` has four production callers; the specification names none for a conditions
twin, so it would be dead weight (decision 56).

Layer: ``transport`` (L2). Imports ``values`` and :mod:`.condition`.
"""

from __future__ import annotations

import operator
from collections.abc import Iterator, KeysView, Mapping, ValuesView
from dataclasses import dataclass
from typing import Any

from orpheus.transport.boundary_conditions.condition import BoundaryCondition
from orpheus.values.content import ContentIdentity, FrozenMapping

__all__ = ["BoundaryConditions"]


@dataclass(frozen=True, eq=False)
class BoundaryConditions(ContentIdentity):
    """The conditions declared for a problem: tag -> condition, in declared order (the twin of ``Materials``)."""

    conditions: Mapping[int, BoundaryCondition]

    def __post_init__(self) -> None:
        declared: dict[int, BoundaryCondition] = {}
        for tag, condition in self.conditions.items():
            if isinstance(tag, bool):
                raise TypeError(f"BoundaryConditions: a tag is an int, got the bool {tag!r}")
            if not isinstance(condition, BoundaryCondition):
                raise TypeError(f"BoundaryConditions: tag {tag!r} holds a {type(condition).__name__}, not a BoundaryCondition")
            declared[operator.index(tag)] = condition
        object.__setattr__(self, "conditions", FrozenMapping(declared))

    @classmethod
    def empty(cls) -> BoundaryConditions:
        """No condition: the declaration of a carrier whose faces are all decks."""
        return cls({})

    @classmethod
    def of(cls, conditions: BoundaryConditions | Mapping[int, BoundaryCondition]) -> BoundaryConditions:
        """Parse-at-boundary coercion: a ``BoundaryConditions`` passes through, a mapping is parsed."""
        return conditions if isinstance(conditions, BoundaryConditions) else cls(conditions)

    def __getitem__(self, tag: int) -> BoundaryCondition:
        return self.conditions[tag]

    def __iter__(self) -> Iterator[int]:
        return iter(self.conditions)

    def __len__(self) -> int:
        return len(self.conditions)

    def __contains__(self, tag: object) -> bool:
        return tag in self.conditions

    def keys(self) -> KeysView[int]:
        return self.conditions.keys()

    def values(self) -> ValuesView[BoundaryCondition]:
        return self.conditions.values()

    def items(self) -> Any:
        return self.conditions.items()

    def get(self, tag: int, default: BoundaryCondition | None = None) -> BoundaryCondition | None:
        return self.conditions.get(tag, default)

    @property
    def tags(self) -> frozenset[int]:
        """The declared tags."""
        return frozenset(self.conditions)
