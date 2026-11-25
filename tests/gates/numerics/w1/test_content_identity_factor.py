r"""The content-identity roster of the factor classes (N03, spec 7.4 roster row).

``Factor``, ``Subquotient`` and ``Restriction`` are ``ContentIdentity`` classes, so each owns one roster ``Entry`` and runs the legs of
``tests/gates/_content_identity_helpers.py``: equal content is one value (equal ``==``, ``hash`` and digest), a moved part moves the digest, the pickle
round trip. ``tests/gates/test_content_identity_census.py`` imports ``ROSTER`` from here for its population gate S5.9; that wiring is asserted by
``test_the_census_imports_this_roster``. First red at ``ff305e0f``: the three classes do not exist (import-only), and the census names
``['Factor', 'Restriction', 'Subquotient']`` as having no roster entry. Mutation: drop an entry (the population gate reddens).
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.numerics.axis import Axis, BasisKind, Declared, Restriction
from orpheus.numerics.declaration import Coord
from orpheus.numerics.factor import Factor, Kind, Subquotient, Variance, VarianceChart
from orpheus.numerics.metric import DiagonalMetric
from tests.gates._content_identity_helpers import (
    Entry,
    check_equal_pair,
    check_perturbation,
    check_pickle,
    check_population,
    leg,
    pair_ids,
    param_id,
    perturbation_ids,
)

pytestmark = pytest.mark.foundation


def _axis(w: tuple[float, ...] = (1.0, 2.0, 3.0, 0.0), label: str = "a") -> Axis:
    return Axis(label, (len(w),), origin=Declared(BasisKind.NODAL, DiagonalMetric(np.array(w))))


def _sk(S: tuple[int, ...] = (0, 1, 2), K: tuple[int, ...] = ()) -> Subquotient:
    return Subquotient(Coord(4, S), Coord(4, K))


_MD = Variance(Kind.MEASURE, VarianceChart.DENSITY)

_SUBQUOTIENT = Entry(
    cls=Subquotient, base=lambda: _sk((0, 1, 2), (1,)), parts=("S", "K"),
    perturb={
        "S": (leg("another S", lambda: _sk((0, 1, 2, 3), (1,))),),
        "K": (leg("another K", lambda: _sk((0, 1, 2), (2,))),),
    },
    pairs=(("two builds", lambda: _sk((0, 1, 2), (1,)), lambda: _sk((0, 1, 2), (1,))),),
)

_FACTOR = Entry(
    cls=Factor, base=lambda: Factor(_axis(), _MD, _sk()), parts=("axis", "variance", "sk"),
    perturb={
        "axis": (leg("another axis", lambda: Factor(_axis((1.0, 2.0, 5.0, 0.0)), _MD, _sk())),),
        "variance": (leg("another variance", lambda: Factor(_axis(), Variance(Kind.FUNCTION, VarianceChart.DENSITY), _sk())),),
        "sk": (leg("another sk", lambda: Factor(_axis(), _MD, _sk((0, 1)))),),
    },
    pairs=(("two builds", lambda: Factor(_axis(), _MD, _sk()), lambda: Factor(_axis(), _MD, _sk())),),
)

_RESTRICTION = Entry(
    cls=Restriction, base=lambda: Restriction("parent-digest", (0, 2)), parts=("parent_key", "indices"),
    perturb={
        "parent_key": (leg("another parent", lambda: Restriction("other-digest", (0, 2))),),
        "indices": (leg("another index set", lambda: Restriction("parent-digest", (0, 3))),),
    },
    pairs=(("two builds", lambda: Restriction("parent-digest", (0, 2)), lambda: Restriction("parent-digest", (0, 2))),),
)

ROSTER: tuple[Entry, ...] = (_FACTOR, _SUBQUOTIENT, _RESTRICTION)


@pytest.mark.parametrize("entry", ROSTER, ids=lambda e: e.id)
def test_population(entry: Entry) -> None:
    check_population(entry)


@pytest.mark.parametrize(
    "entry,part,the_leg",
    [pytest.param(e, p, l, id=param_id(e.id, p, l[0])) for e, p, l in perturbation_ids(ROSTER)],
)
def test_a_moved_part_moves_the_digest(entry: Entry, part: str, the_leg) -> None:
    check_perturbation(entry, part, the_leg)


@pytest.mark.parametrize("entry,pair", [pytest.param(e, p, id=param_id(e.id, p[0])) for e, p in pair_ids(ROSTER)])
def test_equal_content_is_one_value(entry: Entry, pair) -> None:
    check_equal_pair(entry, pair)


@pytest.mark.parametrize("entry", ROSTER, ids=lambda e: e.id)
def test_pickle_round_trip(entry: Entry) -> None:
    check_pickle(entry)


def test_the_census_imports_this_roster() -> None:
    """The three entries are in the census' ROSTER (the population gate S5.9 ranges over it); dropping the import reds here."""
    from tests.gates import test_content_identity_census as census

    classes = {e.cls for e in census.ROSTER}
    missing = [c.__name__ for c in (Factor, Subquotient, Restriction) if c not in classes]
    if missing:
        pytest.fail(f"tests/gates/test_content_identity_census.py ROSTER has no entry for {missing}: import ROSTER from this file")
