r"""The declared conditions: the exact twin of ``Materials`` (wave 2, slice S2a, gate 10).

Spec: ``.claude/plans/wave2_boundary_carve.md`` section 2.4 ("``declared.py``"), decision 56 (no ``restrict``).

``BoundaryConditions`` is what a problem DECLARES: a map from an integer tag to a ``BoundaryCondition``, exactly as
``Materials`` (``orpheus/data/materials.py:72``) is a map from an integer material id to a ``Mixture``. The claims
(``vv-principles``), THEOREM on the type, each with the input that reddens it:

* the keys are parsed to ``int`` (an ``np.integer`` is the same tag; a ``bool`` is refused, ``TypeError``);
* every value is a ``BoundaryCondition``: a deck, a ``str``, a number and a ``WallDisplacement`` are refused
  (``TypeError`` naming the type), the parse that makes a displacement unbindable and a deck unspellable here;
* EMPTY IS ADMITTED, the one difference from ``Materials`` (a carrier of decks only, such as ``pwr_slab_half_cell``,
  has no tag);
* the order of declaration is behaviour (iteration order) and is not content (equal, one digest);
* the mapping surface the readers use (``[]``, ``in``, ``iter``, ``len``, ``keys``, ``values``, ``items``, ``get``,
  ``tags``), and ``of`` as the parse-at-boundary coercion; NO ``restrict`` (no consumer: dead weight).

First red, at the base tip ``4ec95bc6``: ``ModuleNotFoundError``. Mutations: store a ``dict`` instead of a
``FrozenMapping`` (mutable after the fact, the identity rows redden); drop the ``bool`` refusal; drop the
``BoundaryCondition`` check on values; refuse the empty map; add a ``restrict`` method.
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.geometry.deck import MirrorDeck, TranslationDeck
from orpheus.values.content import FrozenMapping, content_digest
from tests.gates._content_identity_helpers import require
from tests.gates.transport._boundary_conditions_helpers import bc, local, node

pytestmark = pytest.mark.foundation

_PARSE = node("test_wall_condition.py", "test_wall_condition_parse")


@pytest.mark.rests_on(_PARSE)
def test_conditions_map() -> None:
    """The twin of ``Materials`` with the differences the section states.

    Reddens on: the mutations of the module docstring. The population of non-conditions refused is listed (a deck
    of each type, a ``str``, an ``int``, ``None``, a displacement, a tuple), so the check is not one example.
    """
    vac, white = bc.WallCondition.vacuum(), bc.WallCondition.white(0.5)
    declared = bc.BoundaryConditions({3: white, 1: vac})
    require(isinstance(declared.conditions, FrozenMapping), "the map is stored as a FrozenMapping, in declared order")
    require(list(declared) == [3, 1] and list(declared.keys()) == [3, 1], f"the declared order is behaviour: {list(declared)}")
    require(list(declared.values()) == [white, vac] and list(declared.items()) == [(3, white), (1, vac)], "values and items follow the order")
    require(declared[3] == white and declared.get(1) == vac and declared.get(7) is None and declared.get(7, vac) == vac, "[] and get")
    require(1 in declared and 3 in declared and 2 not in declared and len(declared) == 2, "in and len")
    require(declared.tags == frozenset({1, 3}) and type(declared.tags) is frozenset, "tags is a frozenset of the keys")
    with pytest.raises(KeyError):
        declared[2]

    reordered = bc.BoundaryConditions({1: vac, 3: white})
    require(reordered == declared and content_digest(reordered) == content_digest(declared) and hash(reordered) == hash(declared),
            "the order of declaration is content")
    require(bc.BoundaryConditions({1: vac}) != declared and bc.BoundaryConditions({1: white, 3: white}) != declared, "a different map is a different value")

    for key in (np.int64(3), np.int32(3)):
        parsed = bc.BoundaryConditions({key: white})  # type: ignore[dict-item]  # numpy-int keys: the coercion row
        require(list(parsed) == [3] and type(next(iter(parsed))) is int, f"{key!r} is not coerced to the int tag")
        require(parsed == bc.BoundaryConditions({3: white}), "an np.integer tag is the same tag")
    for key in (True, False):
        with pytest.raises(TypeError):
            bc.BoundaryConditions({key: vac})
    with pytest.raises((TypeError, ValueError)):
        bc.BoundaryConditions({"1": vac})  # type: ignore[dict-item]  # a str is not a tag

    refused = [
        MirrorDeck(), TranslationDeck(), "vacuum", 1, None, (vac,), 2.0 * bc.WallCondition.specular(0.6),
        bc.WallDisplacement.of((1.0, vac)),
    ]
    for value in refused:
        with pytest.raises(TypeError) as caught:
            bc.BoundaryConditions({1: value})
        require(type(value).__name__ in str(caught.value), f"the refusal of {type(value).__name__} does not name the type: {caught.value}")
    require(len(refused) >= 8, "the refused population shrank")
    probe = local().Probe()
    require(len(bc.BoundaryConditions({1: probe})) == 1, "a method-only BoundaryCondition is a declarable value")

    empty = bc.BoundaryConditions.empty()
    require(len(empty) == 0 and list(empty) == [] and empty.tags == frozenset(), "empty is admitted")
    require(empty == bc.BoundaryConditions({}) and content_digest(empty) == content_digest(bc.BoundaryConditions({})), "empty() is the empty map")
    require(not hasattr(bc.BoundaryConditions, "restrict"), "BoundaryConditions has no restrict (decision 56: no consumer)")

    require(bc.BoundaryConditions.of(declared) is declared, "of() of a declaration is the declaration itself")
    coerced = bc.BoundaryConditions.of({3: white, 1: vac})
    require(isinstance(coerced, bc.BoundaryConditions) and coerced == declared, "of() coerces a bare mapping")
    with pytest.raises(TypeError):
        bc.BoundaryConditions.of({1: "vacuum"})  # type: ignore[dict-item]  # coercion parses once, at the boundary
