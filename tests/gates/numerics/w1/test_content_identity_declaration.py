r"""The content-identity roster of the declaration classes (N01, spec 6.2, 7.2 G01f roster legs).

``Coord``, ``Spanned`` and ``Ann`` are ``ContentIdentity`` classes (a declaration is identity content of ``DenseMetric`` now and of
``Factor`` from N03), so each owns one roster ``Entry`` and runs the legs of ``tests/gates/_content_identity_helpers.py``: equal
content is one value (equal ``==``, ``hash`` and digest), a moved part moves the digest, the pickle round trip.
``tests/gates/test_content_identity_census.py`` imports ``ROSTER`` from here for its population gate S5.9; that wiring is
asserted by ``test_the_census_imports_this_roster`` (the census reddens, as a population gate, the moment the classes exist without
an entry; this row reddens if the import is dropped).
"""

from __future__ import annotations

import pytest

from orpheus.numerics.declaration import Ann, Coord, Spanned
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

_COORD = Entry(
    cls=Coord, base=lambda: Coord(5, (0, 2, 3)), parts=("n", "indices"),
    perturb={
        "n": (leg("another n", lambda: Coord(6, (0, 2, 3))),),
        "indices": (leg("another subset", lambda: Coord(5, (0, 2, 4))), leg("a smaller subset", lambda: Coord(5, (0, 2)))),
    },
    pairs=(("two builds", lambda: Coord(5, (0, 2, 3)), lambda: Coord(5, (0, 2, 3))),),
)

_SPANNED = Entry(
    cls=Spanned, base=lambda: Spanned(5, ("frame-gram-kernel", "ab12"), 2), parts=("n", "key", "dim"),
    perturb={
        "n": (leg("another n", lambda: Spanned(6, ("frame-gram-kernel", "ab12"), 2)),),
        "key": (leg("another key", lambda: Spanned(5, ("frame-gram-kernel", "cd34"), 2)),),
        "dim": (leg("another dim", lambda: Spanned(5, ("frame-gram-kernel", "ab12"), 3)),),
    },
    pairs=(("two builds", lambda: Spanned(5, ("frame-gram-kernel", "ab12"), 2), lambda: Spanned(5, ("frame-gram-kernel", "ab12"), 2)),),
)

_ANN = Entry(
    cls=Ann, base=lambda: Ann(5, Spanned(5, ("k", 1), 2)), parts=("n", "of"),
    perturb={
        "n": (leg("another n", lambda: Ann(6, Spanned(6, ("k", 1), 2)), "of"),),
        "of": (leg("another subspace", lambda: Ann(5, Spanned(5, ("k", 2), 2))),),
    },
    pairs=(("two builds", lambda: Ann(5, Spanned(5, ("k", 1), 2)), lambda: Ann(5, Spanned(5, ("k", 1), 2))),),
)

ROSTER: tuple[Entry, ...] = (_COORD, _SPANNED, _ANN)


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


@pytest.mark.rests_on("tests/gates/numerics/w1/test_declaration_normal_form.py::test_annihilator_normal_form")
def test_the_census_imports_this_roster() -> None:
    """The three entries are in the census' ROSTER (the population gate S5.9 ranges over it); dropping the import reds here."""
    from tests.gates import test_content_identity_census as census

    classes = {e.cls for e in census.ROSTER}
    missing = [c.__name__ for c in (Coord, Spanned, Ann) if c not in classes]
    if missing:
        pytest.fail(f"tests/gates/test_content_identity_census.py ROSTER has no entry for {missing}: import ROSTER from this file")
