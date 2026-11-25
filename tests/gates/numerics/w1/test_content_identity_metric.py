r"""The content-identity roster of the metric arms ``DiagonalMetric``, ``DenseMetric`` and ``DualMetric`` (N04, spec 7.5 row G04q, decision ``[ON04-6]``).

Each arm is a ``ContentIdentity`` class, so each owns one roster ``Entry`` and runs the legs of ``tests/gates/_content_identity_helpers.py``: equal content is one value (equal
``==``, ``hash`` and digest), a moved part moves the digest, the pickle round trip. ``tests/gates/test_content_identity_census.py`` imports ``ROSTER`` from here for its population
gate S5.9; that wiring is asserted by ``test_the_census_imports_this_roster``. First red at ``cf924aaa``: the census names ``['DenseMetric', 'DiagonalMetric', 'DualMetric']`` as having no roster
entry once the arms are on ``ContentIdentity``. Mutation: drop an entry (the population gate reddens). Exact on any host (small stated matrices).
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.numerics.declaration import Coord, Spanned
from orpheus.numerics.metric import DenseMetric, DiagonalMetric, DualMetric, Face
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


def _diag(w: tuple[float, ...] = (1.0, 2.0, 3.0)) -> DiagonalMetric:
    return DiagonalMetric(np.array(w))


def _dense(g: tuple[tuple[float, float], tuple[float, float]] = ((1.0, 1.0), (1.0, 1.0)), face: Face = Face.METRIC, key: str = "a") -> DenseMetric:
    return DenseMetric(np.array(g), kernel=Spanned(2, (key,), 1), face=face)


_DIAGONAL = Entry(
    cls=DiagonalMetric, base=_diag, parts=("weights",),
    perturb={"weights": (leg("another weight", lambda: _diag((1.0, 2.0, 3.5))),)},
    pairs=(("two builds", _diag, _diag),),
)

_DENSE = Entry(
    cls=DenseMetric, base=_dense, parts=("stated", "kernel", "face"),
    perturb={
        "stated": (leg("another matrix", lambda: _dense(((2.0, 2.0), (2.0, 2.0)))),),
        "kernel": (leg("another kernel", lambda: _dense(key="b")),),
        "face": (leg("another face", lambda: _dense(face=Face.GRAM)),),
    },
    pairs=(("two builds", _dense, _dense),),
)

_DUAL = Entry(
    cls=DualMetric, base=lambda: DualMetric(_diag()), parts=("primal",),
    perturb={"primal": (leg("another primal", lambda: DualMetric(_diag((1.0, 2.0, 3.5)))),)},
    pairs=(("two builds", lambda: DualMetric(_diag()), lambda: DualMetric(_diag())),),
    fields_are_parts=False,
)

ROSTER: tuple[Entry, ...] = (_DIAGONAL, _DENSE, _DUAL)


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
    missing = [c.__name__ for c in (DiagonalMetric, DenseMetric, DualMetric) if c not in classes]
    if missing:
        pytest.fail(f"tests/gates/test_content_identity_census.py ROSTER has no entry for {missing}: import ROSTER from this file")


def test_spanned_and_coord_kernels_are_content() -> None:
    """A dense form's kernel is a ``Spanned`` or a ``Coord``: ``==`` is a ``bool`` and the hash is stable on each (a contentless part would fall back to identity)."""
    g = np.diag([2.0, 3.0, 0.0])
    for kernel in (Coord(3, (2,)), Spanned(2, ("a",), 1)):
        m = g if isinstance(kernel, Coord) else np.ones((2, 2))
        a = DenseMetric(m, kernel=kernel, face=Face.METRIC)
        b = DenseMetric(m.copy(), kernel=kernel, face=Face.METRIC)
        if not ((a == b) is True and hash(a) == hash(b)):
            pytest.fail(f"a dense form with kernel {kernel!r} is not content")
        if not ((DualMetric(a) == DualMetric(b)) is True and hash(DualMetric(a)) == hash(DualMetric(b))):
            pytest.fail(f"the dual of a dense form with kernel {kernel!r} is not content")
