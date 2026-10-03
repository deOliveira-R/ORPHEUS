r"""Content identity of the observables (#405 P2 step 3, R3.4-R3.5).

Specified by the test-architect (2026-10-03, ``.claude/plans/reference_p2_spec.md``
§1.3). Its ``ROSTER`` joins the union in
``tests/gates/numerics/test_content_identity.py`` (S5.9 walks
``orpheus.numerics`` and reds while a content class has no roster entry;
S5.6's route gate then covers the four types).

An observable is a cache key's part (a reference's reading is stored per
observable), so: every part moves the digest (the population is the type's
``content_parts``); a ratio is ORDERED (``Ratio(a, b)`` and ``Ratio(b, a)``
are two values); equal content is one value; the digests are the same in every
process; a RECORD row pins their bytes after landing.
"""

from __future__ import annotations

import importlib
import os
import subprocess
import sys
from pathlib import Path
from typing import Any

import numpy as np
import pytest

from orpheus.numerics.content import content_digest
from orpheus.numerics.mesh_free_function import RegionwiseConstant, Symbolic
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
    require,
)

pytestmark = pytest.mark.foundation

_HERE = "tests/gates/numerics/test_content_identity_observable.py"
_ROOT = Path(__file__).resolve().parents[3]
_MODULE = "orpheus.numerics.observable"
_TABLE = np.array([[1.5, 0.25], [0.0, 3.0]])


def _c(name: str) -> Any:
    return getattr(importlib.import_module(_MODULE), name)


def _table(values: np.ndarray = _TABLE) -> RegionwiseConstant:
    return RegionwiseConstant(values.copy())


def _ulp_table() -> RegionwiseConstant:
    moved = _TABLE.copy()
    moved[1, 1] = np.nextafter(moved[1, 1], np.inf)
    return RegionwiseConstant(moved)


def _flux(weight: Any = None) -> Any:
    return _c("FluxIntegral")(_table() if weight is None else weight)


def _point(position: float = 0.75, group: int = 1) -> Any:
    return _c("PointValue")(position, group)


def _ratio(numerator: Any = None, denominator: Any = None) -> Any:
    return _c("Ratio")(_flux() if numerator is None else numerator, _point() if denominator is None else denominator)


ROSTER: tuple[Entry, ...] = (
    Entry(
        cls=_c("FluxIntegral"),
        base=_flux,
        parts=("weight",),
        perturb={"weight": (
            leg("one ULP of one coefficient", lambda: _flux(_ulp_table())),
            leg("a symbolic weight of the same group count", lambda: _flux(Symbolic.of(1 + Symbolic.r, Symbolic.r**2))),
        )},
        pairs=(("two builds", _flux, _flux),
               ("int and float table", lambda: _flux(RegionwiseConstant(np.array([[1, 0], [2, 3]]))),
                lambda: _flux(RegionwiseConstant(np.array([[1.0, 0.0], [2.0, 3.0]]))))),
    ),
    Entry(
        cls=_c("Ratio"),
        base=_ratio,
        parts=("numerator", "denominator"),
        perturb={
            "numerator": (leg("another weight", lambda: _ratio(numerator=_flux(_ulp_table()))),
                          leg("a point value", lambda: _ratio(numerator=_point(0.25, 0)))),
            "denominator": (leg("another group", lambda: _ratio(denominator=_point(group=0))),
                            leg("a flux integral", lambda: _ratio(denominator=_flux()))),
        },
        pairs=(("two builds", _ratio, _ratio),),
    ),
    Entry(cls=_c("Eigenvalue"), base=lambda: _c("Eigenvalue")(), parts=(), perturb={},
          pairs=(("two builds", lambda: _c("Eigenvalue")(), lambda: _c("Eigenvalue")()),)),
    Entry(
        cls=_c("PointValue"),
        base=_point,
        parts=("position", "group"),
        perturb={
            "position": (leg("one ULP", lambda: _point(position=float(np.nextafter(0.75, 1.0)))),
                         leg("the sign", lambda: _point(position=-0.75))),
            "group": (leg("another group", lambda: _point(group=0)),),
        },
        pairs=(("int and float position", lambda: _point(position=1), lambda: _point(position=1.0)),
               ("-0.0 and 0.0 position", lambda: _point(position=-0.0), lambda: _point(position=0.0))),
    ),
)


@pytest.mark.parametrize("entry", ROSTER, ids=[e.id for e in ROSTER])
def test_r3_4_the_population_is_the_parts(entry: Entry) -> None:
    check_population(entry)


@pytest.mark.parametrize(
    "entry, part, the_leg",
    perturbation_ids(ROSTER),
    ids=[param_id(e.id, p, lg[0]) for e, p, lg in perturbation_ids(ROSTER)],
)
def test_r3_4_each_part_moves_the_digest(entry: Entry, part: str, the_leg: Any) -> None:
    check_perturbation(entry, part, the_leg)


@pytest.mark.parametrize("entry, pair", pair_ids(ROSTER), ids=[param_id(e.id, p[0]) for e, p in pair_ids(ROSTER)])
def test_r3_4_equal_content_is_one_value(entry: Entry, pair: Any) -> None:
    check_equal_pair(entry, pair)


@pytest.mark.parametrize("entry", ROSTER, ids=[e.id for e in ROSTER])
def test_r3_4_pickle(entry: Entry) -> None:
    check_pickle(entry)


@pytest.mark.rests_on(f"{_HERE}::test_r3_4_each_part_moves_the_digest")
def test_r3_4_a_ratio_is_ordered_and_kinds_never_merge() -> None:
    """``Ratio(a, b) != Ratio(b, a)`` (two digests, two set members); and the
    kinds never merge: a ``FluxIntegral`` of a table, the ``Eigenvalue``, a
    ``PointValue`` and a ``Ratio`` of the first two are four values, none
    equal to the weight it holds."""
    a, b = _flux(), _point()
    forward, backward = _c("Ratio")(a, b), _c("Ratio")(b, a)
    require(forward != backward and len({forward, backward}) == 2, "a ratio and its inverse merged")
    require(content_digest(forward) != content_digest(backward), "a ratio and its inverse share a digest")
    values = [a, _c("Eigenvalue")(), b, forward]
    require(len(set(values)) == 4, "two kinds of observable merged in a set")
    require(a != a.weight, "a flux integral equals its weight")


# ═════════════════════════════════════════════════════════════════════════════
# R3.5 — the same digests in every process; the RECORD bytes
# ═════════════════════════════════════════════════════════════════════════════

_SEED_SCRIPT = """
import numpy as np
import orpheus.numerics.observable as m
from orpheus.numerics.content import content_digest
from orpheus.numerics.mesh_free_function import RegionwiseConstant, Symbolic
t = RegionwiseConstant(np.array([[1.5, 0.25], [0.0, 3.0]]))
f = m.FluxIntegral(t)
print(m.__file__)
for v in (f, m.FluxIntegral(Symbolic.of(1 + Symbolic.r, Symbolic.r**2)), m.Eigenvalue(),
          m.PointValue(0.75, 1), m.Ratio(f, m.PointValue(0.75, 1)), m.Ratio(m.PointValue(0.25, 0), f)):
    print(content_digest(v).hex(), hash(v))
print('control', hash('salted'))
"""


def _seed_run(seed: str) -> list[str]:
    env = {**os.environ, "PYTHONHASHSEED": seed, "PYTHONPATH": str(_ROOT)}
    out = subprocess.run([sys.executable, "-O", "-c", _SEED_SCRIPT], cwd=_ROOT, env=env, capture_output=True, text=True)
    require(out.returncode == 0, out.stderr)
    return out.stdout.strip().splitlines()


@pytest.mark.rests_on("tests/gates/numerics/test_content_identity.py::test_s5_1_digests_and_hashes_are_seed_stable")
def test_r3_5_digests_and_hashes_are_seed_stable() -> None:
    """``PYTHONHASHSEED`` 1 and 2 print the same six digests and hashes (a
    regionwise and a symbolic weight, the eigenvalue, a point value, a ratio,
    a ratio of a point value to a flux integral); the control line must differ."""
    one, two = _seed_run("1"), _seed_run("2")
    require(one[0].startswith(str(_ROOT)), f"the subprocess imported {one[0]} (L22)")
    require(len(one) == 8, f"activation: {len(one)} lines printed")
    require(one[1:-1] == two[1:-1], f"digests or hashes moved with the seed:\n{one}\n{two}")
    require(one[-1] != two[-1], "control: the str hash did not move, so the seeds did not take effect")


_PIN_AFTER_LANDING = "pin after landing"
_PINNED: dict[str, str] = {
    "flux integral of a table": "24c25dfbe46136af67720f38b891ab450d72c6b3b0bd5619697b36ff182a033a",  # pinned 2026-10-03 at landing
    "eigenvalue": "344afce0a7478f8f83fa0f179987b067df0bbcd83430aace6be94e80df338d1e",  # pinned 2026-10-03 at landing
    "point value": "cff6ea6dfc995c344ca4cc27f892e4284d5217031b132e2975c5781888dee720",  # pinned 2026-10-03 at landing
    "ratio of a flux integral to a point value": "c5e2a64cb1b15d4980fdfafcbb0a20afe661ccb1f9822a80e1ba235db2708691",  # pinned 2026-10-03 at landing
}


def _fingerprint_values() -> dict[str, Any]:
    return {
        "flux integral of a table": _flux(),
        "eigenvalue": _c("Eigenvalue")(),
        "point value": _point(),
        "ratio of a flux integral to a point value": _ratio(),
    }


@pytest.mark.rests_on(f"{_HERE}::test_r3_5_digests_and_hashes_are_seed_stable")
def test_r3_5_the_digest_bytes_are_recorded() -> None:
    """RECORD (designed to red on any edit of the observable schema or the
    encoder): the producer fingerprint the P3 cache keys rest on. Until
    pinned it fails with ``pin after landing`` and prints the digests; re-pin
    only with the reason, in the commit message."""
    got = {name: content_digest(value).hex() for name, value in _fingerprint_values().items()}
    if any(pin == _PIN_AFTER_LANDING for pin in _PINNED.values()):
        pytest.fail(f"pin after landing: replace the placeholders of _PINNED with {got}")
    for name, pin in _PINNED.items():
        require(got[name] == pin, f"{name}: the digest bytes moved ({got[name]} != {pin}): every cache key "
                                  f"holding an observable is invalidated; re-pin only with the reason")
