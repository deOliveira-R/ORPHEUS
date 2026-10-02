r"""Content identity of the mesh-free functions (#405 P1 step 6, S6.9 and S6.16).

:class:`~orpheus.numerics.phase_space_function.RegionwiseConstant` and
:class:`~orpheus.numerics.phase_space_function.Symbolic` join the content
rosters; S5.9's population row (``test_content_identity.py``) reads this
one, so a new mesh-free type added without an entry reds there.

The legs that are the encoder's own rules, re-read on these types: a table
entry moved by one ULP is another value; a ``(2, 3)`` and a ``(3, 2)`` table
with the same bytes are two values (the shape is content); ``-0.0`` is
``0.0`` and an integer table is its float twin. A ``Symbolic`` is its text and
the version of the SymPy that wrote it.
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.numerics.phase_space_function import RegionwiseConstant, Symbolic
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

_TABLE = np.array([[1.5, 0.25, 2.0], [0.0, 3.0, 0.75]])


def _symbolic(*extra: str, version: str | None = None) -> Symbolic:
    import sympy as sp

    expressions = (1 + Symbolic.mu * sp.exp(-Symbolic.r), Symbolic.r**2 + sp.cos(Symbolic.phi))
    return Symbolic.from_srepr(tuple(sp.srepr(e) for e in expressions) + extra, version)


def _symbolic_moved_coefficient() -> Symbolic:
    import sympy as sp

    expressions = (1 + Symbolic.mu * sp.exp(-Symbolic.r), Symbolic.r**2 + 2 * sp.cos(Symbolic.phi))
    return Symbolic.of(*expressions)


_REGIONWISE = Entry(
    cls=RegionwiseConstant, base=lambda: RegionwiseConstant(_TABLE.copy()), parts=("values",),
    perturb={
        "values": (
            leg("one entry by one ULP", lambda: RegionwiseConstant(np.where(
                np.arange(_TABLE.size).reshape(_TABLE.shape) == 4, np.nextafter(_TABLE, np.inf), _TABLE))),
            leg("the same bytes in the transposed shape", lambda: RegionwiseConstant(_TABLE.reshape(3, 2).copy())),
        ),
    },
    pairs=(
        ("two builds", lambda: RegionwiseConstant(_TABLE.copy()), lambda: RegionwiseConstant(_TABLE.copy())),
        ("-0.0 vs 0.0", lambda: RegionwiseConstant(np.where(_TABLE == 0.0, -0.0, _TABLE)),
         lambda: RegionwiseConstant(_TABLE.copy())),
        ("an integer table vs its float twin", lambda: RegionwiseConstant(np.array([[1, 2], [3, 4]])),
         lambda: RegionwiseConstant(np.array([[1.0, 2.0], [3.0, 4.0]]))),
    ),
)

_SYMBOLIC = Entry(
    cls=Symbolic, base=_symbolic, parts=("srepr", "sympy_version"),
    perturb={
        "srepr": (
            leg("one coefficient", _symbolic_moved_coefficient),
            leg("a third group", lambda: _symbolic("Integer(1)")),
        ),
        "sympy_version": (leg("another SymPy", lambda: _symbolic(version="0.0.0-other")),),
    },
    pairs=(
        ("two builds", _symbolic, _symbolic),
        ("text and expressions", _symbolic, lambda: Symbolic.of(*_symbolic().expressions)),
    ),
)

ROSTER: tuple[Entry, ...] = (_REGIONWISE, _SYMBOLIC)


@pytest.mark.parametrize("entry", ROSTER, ids=lambda e: e.id)
def test_s6_9_population(entry: Entry) -> None:
    check_population(entry)


@pytest.mark.parametrize(
    "entry,part,the_leg",
    [pytest.param(e, p, l, id=param_id(e.id, p, l[0])) for e, p, l in perturbation_ids(ROSTER)],
)
def test_s6_9_a_moved_part_moves_the_digest(entry: Entry, part: str, the_leg) -> None:
    check_perturbation(entry, part, the_leg)


@pytest.mark.parametrize(
    "entry,pair", [pytest.param(e, p, id=param_id(e.id, p[0])) for e, p in pair_ids(ROSTER)]
)
def test_s6_9_equal_content_is_one_value(entry: Entry, pair) -> None:
    check_equal_pair(entry, pair)


@pytest.mark.parametrize("entry", ROSTER, ids=lambda e: e.id)
def test_s6_9_pickle_round_trip(entry: Entry) -> None:
    check_pickle(entry)
