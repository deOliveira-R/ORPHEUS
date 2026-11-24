r"""N02, the seventh cutoff: the three readers of a cylindrical problem's degenerate radial cosine agree with the cosines, exactly.

A cylindrical ``SNProblem`` decides which ordinates have a radial cosine of exactly zero (the tangential ordinates, which cross no
radial face) at three readers: ``SNProblem._iter_cylindrical_visits`` (the visit's downstream face area), ``SNProblem._representative_ordinate``
(the candidate of a sign class) and ``StreamingCoefficientCache.is_degenerate`` (the sweep cache's flag). Commit ``26f9815a`` deleted the class
constant ``_DEGENERATE_ABS_ETA_THRESHOLD = 1e-15``, which made the first two read ``|eta| < 1e-15`` while the cache read the exact zero of
``half_range``, so a cosine in ``(0, 1e-15)`` was degenerate at two readers and genuine at the third.

Fixtures, both on ``Quadrature.folded_product(2, 6)`` (the one rule of the tree whose levels have exactly-zero radial cosines and that a
cylindrical problem admits):

* the shipped rule: radial cosines ``mu_x = [0.7071, 0, -0.7071, 0.7071, 0, -0.7071]``, exactly zero at ordinates 1 and 4;
* the nudged rule: the same nodes with each exactly-zero ``mu_x`` replaced by ``9.5e-16`` (outward, and inside the old band).

Both problems are posed with vacuum on the outer wall. A specular wall refuses the nudged rule at construction (the pairing about the axis maps
the ``9.5e-16`` ordinate to itself, ERR-045), which is the admission predicate working and not a weakened gate: the readers under test do not
read the wall, so the two fixtures carry the same wall and differ only in the cosines.

Every expected value is computed here from the cosines by ``np.sign`` against the exact zero, never from ``half_range``.

First red on the tip ``363976fb`` is by mutation (the rows are green there, since the commit that deleted the constant is the fix): the
threshold put back into ``_iter_cylindrical_visits`` reddens the nudged visit row, and put back into ``_representative_ordinate`` reddens the
nudged representative row; see the unit's report for the red lines.
"""

from __future__ import annotations

from dataclasses import replace
from typing import Callable

import numpy as np
import pytest

from orpheus.geometry import StructuredGeometry
from orpheus.mesh import CellsByCount, Mesher
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.problem import SNProblem
from orpheus.sn.sweep.cache import StreamingCoefficientCache
from orpheus.transport.boundary_conditions import WallCondition
from tests.gates._exterior import conditions_for
from tests.gates._half_range_helpers import census_rules
from tests.gates.sn._test_helpers import placeholder_materials

pytestmark = pytest.mark.foundation

NX = 4
NUDGE = 9.5e-16

_OWNER_TESTS = "tests/gates/numerics/w1/test_half_range_one_owner.py"


def _shipped() -> Quadrature:
    return Quadrature.folded_product(n_mu=2, n_phi=6)


def _nudged() -> Quadrature:
    r"""The shipped rule with every exactly-zero ``mu_x`` replaced by ``9.5e-16``, through the measure's nodes."""
    shipped = _shipped()
    nodes = np.array(shipped.measure.nodes, copy=True)
    nodes[nodes[:, 0] == 0.0, 0] = NUDGE
    return replace(shipped, measure=replace(shipped.measure, nodes=nodes))


def _cylinder(quad: Quadrature) -> SNProblem:
    mesh = (
        Mesher(StructuredGeometry.cylinder((0.0, 1.5), (0,), outer=0))
        .partition(CellsByCount.uniform_width(NX))
        .mesh
    )
    return SNProblem(
        mesh, quad, placeholder_materials(), boundary_conditions=conditions_for(mesh, {0: WallCondition.vacuum()})
    )


FIXTURES: dict[str, Callable[[], Quadrature]] = {"shipped": _shipped, "nudged": _nudged}


def _radial_sign(quad: Quadrature) -> np.ndarray:
    r"""Hand classification of every ordinate: the sign of its radial cosine ``mu_x`` against the exact zero, ``(N,)`` int."""
    return np.sign(np.asarray(quad.mu_x, dtype=np.float64)).astype(np.int64)


def _levels(quad: Quadrature) -> list[list[int]]:
    return [[int(g) for g in level] for level in quad.level_indices]


def _visits(problem: SNProblem, within: int, level: int) -> tuple[list[float], list[int]]:
    r"""``(downstream face area per visit, cell index per visit)`` of one ordinate, ``within`` its level ``level``."""
    visits = list(problem._iter_cylindrical_visits(within, level))
    return [float(v.face_area_downstream) for v in visits], [int(v.cell_idx) for v in visits]


def _visit_violations(quad: Quadrature) -> list[str]:
    r"""What ``_iter_cylindrical_visits`` gets wrong against the hand classification of ``quad`` (an empty list is correct)."""
    problem = _cylinder(quad)
    sign = _radial_sign(quad)
    bad: list[str] = []
    forward, backward = list(range(NX)), list(range(NX - 1, -1, -1))
    for lv, level in enumerate(_levels(quad)):
        for within, n in enumerate(level):
            areas, cells = _visits(problem, within, lv)
            if len(areas) != NX:
                bad.append(f"ordinate {n}: {len(areas)} visits, wanted {NX}")
                continue
            all_zero = all(a == 0.0 for a in areas)
            if all_zero != (sign[n] == 0):
                bad.append(f"ordinate {n} (mu_x={float(quad.mu_x[n])!r}): all downstream areas zero is {all_zero}, wanted {sign[n] == 0}")
            if sign[n] != 0 and not any(a != 0.0 for a in areas):
                bad.append(f"ordinate {n}: no nonzero downstream area although the cosine is nonzero")
            if cells != (backward if sign[n] < 0 else forward):
                bad.append(f"ordinate {n}: cell order {cells}, wanted {'backward' if sign[n] < 0 else 'forward'}")
    return bad


def _representative_violations(quad: Quadrature) -> list[str]:
    r"""What ``_representative_ordinate`` gets wrong against the hand classification of ``quad`` (an empty list is correct).

    The return value is the index within the level. The hand answer is the first ordinate of the level whose cosine has the asked sign.
    """
    problem = _cylinder(quad)
    sign = _radial_sign(quad)
    bad: list[str] = []
    for lv, level in enumerate(_levels(quad)):
        for direction in (+1, -1):
            got = int(problem._representative_ordinate(direction, lv))
            hand = [w for w, n in enumerate(level) if sign[n] == direction]
            if not 0 <= got < len(level):
                bad.append(f"level {lv} sign {direction:+d}: returned {got}, outside the level of {len(level)}")
            elif sign[level[got]] != direction:
                bad.append(f"level {lv} sign {direction:+d}: returned ordinate {level[got]} of sign {sign[level[got]]:+d}")
            elif got != hand[0]:
                bad.append(f"level {lv} sign {direction:+d}: returned within-level {got}, the first of the sign class is {hand[0]}")
    return bad


def _cache_violations(quad: Quadrature) -> list[str]:
    r"""What the sweep cache's degenerate flags get wrong against the hand classification of ``quad`` (an empty list is correct)."""
    flags = np.asarray(StreamingCoefficientCache.from_mesh_and_quad(_cylinder(quad)).is_degenerate)
    want = _radial_sign(quad) == 0
    if flags.dtype != np.bool_ or flags.shape != want.shape:
        return [f"flags have dtype {flags.dtype} and shape {flags.shape}, wanted bool {want.shape}"]
    return [] if np.array_equal(flags, want) else [f"flags {flags.tolist()}, wanted {want.tolist()}"]


def _degenerate_sets(quad: Quadrature) -> dict[str, tuple[int, ...]]:
    r"""The ordinates each of the three readers calls degenerate, as sorted global indices."""
    problem = _cylinder(quad)
    visit_set: list[int] = []
    for lv, level in enumerate(_levels(quad)):
        for within, n in enumerate(level):
            if all(a == 0.0 for a in _visits(problem, within, lv)[0]):
                visit_set.append(n)
    flags = StreamingCoefficientCache.from_mesh_and_quad(problem).is_degenerate
    return {
        "visits": tuple(sorted(visit_set)),
        "cache": tuple(int(i) for i in np.flatnonzero(flags)),
        "hand": tuple(int(i) for i in np.flatnonzero(_radial_sign(quad) == 0)),
    }


# ── row a: the shipped rule ──────────────────────────────────────────


def test_the_shipped_rule_has_exactly_zero_radial_cosines_at_ordinates_one_and_four() -> None:
    r"""Control of the fixture: the degenerate ordinates of the shipped rule are 1 and 4, and no other cosine is within ``1e-15`` of zero."""
    mu = np.asarray(_shipped().mu_x)
    assert np.flatnonzero(mu == 0.0).tolist() == [1, 4], mu
    assert np.flatnonzero(np.abs(mu) < 1e-15).tolist() == [1, 4]
    assert float(np.min(np.abs(mu[mu != 0.0]))) > 0.5


@pytest.mark.rests_on(f"{_OWNER_TESTS}::test_the_owner_is_exact_on_a_round_off_cosine")
def test_visits_of_the_shipped_rule_have_zero_downstream_area_exactly_at_the_zero_cosines() -> None:
    r"""On the shipped rule every visit of ordinates 1 and 4 has downstream face area ``0.0``, and no other ordinate has all areas zero.

    The inward ordinates end at the axis, where the inner face has zero area, so a single zero area does not mark a degenerate ordinate:
    the predicate is that every visit is zero. Cell order is forward for a cosine of ``>= 0`` and backward for a negative one.
    """
    assert _visit_violations(_shipped()) == []


@pytest.mark.rests_on(f"{_OWNER_TESTS}::test_the_owner_is_exact_on_a_round_off_cosine")
def test_representative_ordinate_of_the_shipped_rule_is_never_a_zero_cosine() -> None:
    r"""``_representative_ordinate(+1, level)`` and ``(-1, level)`` return an ordinate of cosine ``> 0`` and ``< 0``, the first of its class."""
    assert _representative_violations(_shipped()) == []


@pytest.mark.rests_on(f"{_OWNER_TESTS}::test_the_owner_is_exact_on_a_round_off_cosine")
def test_the_sweep_cache_marks_the_shipped_zero_cosines() -> None:
    r"""``StreamingCoefficientCache.is_degenerate`` is true exactly at ordinates 1 and 4 of the shipped rule."""
    assert _cache_violations(_shipped()) == []


# ── row b: the nudged rule (the row the 1e-15 threshold fails) ───────


def test_the_nudged_rule_has_no_zero_cosine_and_two_cosines_inside_the_old_band() -> None:
    r"""Control of the fixture: ordinates 1 and 4 carry ``9.5e-16``, strictly positive, inside the retired ``1e-15`` band."""
    mu = np.asarray(_nudged().mu_x)
    assert not np.any(mu == 0.0)
    assert np.flatnonzero(np.abs(mu) < 1e-15).tolist() == [1, 4]
    assert [float(mu[i]) for i in (1, 4)] == [NUDGE, NUDGE]


@pytest.mark.rests_on(f"{_OWNER_TESTS}::test_the_owner_is_exact_on_a_round_off_cosine")
def test_visits_of_the_nudged_rule_have_no_degenerate_ordinate() -> None:
    r"""On the nudged rule no ordinate has all downstream areas zero: the ``9.5e-16`` ordinates are outward, iterate forward over the
    cells and carry the nonzero face areas of the outer faces. The retired ``|eta| < 1e-15`` test read them as degenerate."""
    assert _visit_violations(_nudged()) == []


@pytest.mark.rests_on(f"{_OWNER_TESTS}::test_the_owner_is_exact_on_a_round_off_cosine")
def test_representative_ordinate_of_the_nudged_rule_is_the_first_of_its_class() -> None:
    r"""On the nudged rule the first outward ordinate of each level is the ``9.5e-16`` one (within-level index 1), and ``_representative_ordinate``
    returns it. This is the one sign-class row the retired ``> 1e-15`` test fails, and it pins the choice of the first candidate: any
    ordinate of the class (the ``0.7071`` one) has the same cell order, so the sign class alone cannot see the threshold on this rule."""
    assert _representative_violations(_nudged()) == []


@pytest.mark.rests_on(f"{_OWNER_TESTS}::test_the_owner_is_exact_on_a_round_off_cosine")
def test_the_sweep_cache_marks_no_nudged_ordinate() -> None:
    r"""``StreamingCoefficientCache.is_degenerate`` is false at every ordinate of the nudged rule."""
    assert _cache_violations(_nudged()) == []


@pytest.mark.rests_on(
    f"{_OWNER_TESTS}::test_the_owner_is_exact_on_a_round_off_cosine",
)
@pytest.mark.parametrize("name", list(FIXTURES))
def test_the_three_readers_agree_on_the_degenerate_ordinates(name: str) -> None:
    r"""The visit reader, the sweep cache and the hand classification name the same degenerate ordinates: ``(1, 4)`` on the shipped rule
    and none on the nudged rule."""
    sets = _degenerate_sets(FIXTURES[name]())
    assert sets["visits"] == sets["cache"] == sets["hand"], sets
    assert sets["hand"] == ((1, 4) if name == "shipped" else ()), sets


# ── row c: the three cosine spellings are the same bytes ─────────────

#: The 1-D rules of the census (17) beside the folded products of the fixtures; of these, those a cylindrical problem admits are the
#: folded products. The census is the tree's list of shipped constructors (``tests/gates/_half_range_helpers.py::census_rules``).
_CENSUS = census_rules()
_BYTES_RULES = [(n, f) for n, f in _CENSUS] + [("shipped", _shipped), ("nudged", _nudged)]


@pytest.mark.parametrize("name, build", _BYTES_RULES, ids=[n for n, _ in _BYTES_RULES])
def test_mean_axis_cosine_mu_x_and_eta_are_the_same_bytes(name: str, build: Callable[[], Quadrature]) -> None:
    r"""``quad.mean_axis_cosine(0)`` (read by the owner), ``quad.mu_x`` and ``quad.eta`` (read by the walks) are one array's bytes.

    Population: the 17 census rules of ``census_rules()`` plus the two fixtures, 19 rules; a cylindrical ``SNProblem`` admits 3 of the 17
    (the folded products) and both fixtures.
    """
    q = build()
    a = np.asarray(q.mean_axis_cosine(0), dtype=np.float64).tobytes()
    assert a == np.asarray(q.mu_x, dtype=np.float64).tobytes(), f"{name}: mean_axis_cosine(0) and mu_x differ"
    assert a == np.asarray(q.eta, dtype=np.float64).tobytes(), f"{name}: mean_axis_cosine(0) and eta differ"


def test_the_row_c_population_is_17_census_rules_and_the_cylinder_admits_the_folded_products() -> None:
    r"""The population statement of row c is a measurement: 17 census rules, of which exactly the 3 folded products build a cylindrical problem."""
    assert len(_CENSUS) == 17
    admitted: list[str] = []
    for name, build in _CENSUS:
        try:
            _cylinder(build())
        except ValueError:
            continue
        admitted.append(name)
    assert admitted == ["fold2x4", "fold4x8", "fold8x8"], admitted
